"""Drafts, for every study that has no row, the kind of its trait: binary, with its cases and controls, or quantitative.

  uv run python -m cases [RELEASE]

RELEASE defaults to $SUMSTATS_RELEASE. A study of a binary trait has cases and controls, and its power goes by its
effective N, 4 / (1 / cases + 1 / controls), not by their sum. They are those of a record of the study's sample:

  1. the GWAS Catalog's record of the study, if it is one (data/catalog.tsv), or of the source that its corrected
     copy was taken from (data/corrections.tsv);
  2. else the record of its paper (data/samples.tsv) whose sample is the study's, by its N: the N of the inventory,
     or the largest N of its file's SNPs (data/metrics.tsv), that of a SNP every cohort has. No SNP of a file has
     more individuals than its sample: a record that the file's largest N exceeds is not the study's;
  3. else, for a study of the Pan-UK Biobank, the phenotype of its N in the Pan-UK Biobank's manifest (downloaded to
     .cache/), if one phenotype only has that N.

A record that counts individuals, not cases, makes the study quantitative; so does a paper of which the Catalog
counts individuals only, and a Pan-UK Biobank N that only measurements have. The Catalog counts individuals for some
binary traits too, as for the yes-or-no traits of UK Biobank in Loh et al. 2018: read the draft, and correct such a
row by hand.

The result is data/cases.tsv: study ID, kind (binary or quantitative), cases, controls, source (the record they are
from) and how they were found. A row that says "check" could not be drafted. Complete a row by hand and write "by
hand" in `how`, with where its numbers are from: a later run keeps the rows that are there.
"""
import argparse
import csv
import gzip
import os
import tomllib
import urllib.request
from pathlib import Path

import pandas as pd

from build import merge, overlay, release as rel
from corrections import fix
from corrections.__main__ import planned
from papers import cite
from qc.checks import ANCESTRIES, recorded

DATA = merge.DATA
CACHE = DATA.parent / ".cache" / "panukb"
MANIFEST = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz"  # every phenotype of the Pan-UK Biobank, with its cases and controls by ancestry
MEASURED = {"continuous", "biomarkers"}  # the manifest's kinds of phenotype that are measurements
HEAD = ["id", "kind", "cases", "controls", "source", "how"]
ABOVE = 0.001  # the largest N of a file may be above its sample by this share, for rounding, and no more


def manifest(cache: Path = CACHE, refresh: bool = False) -> pd.DataFrame:
    """The Pan-UK Biobank's manifest of phenotypes, fetched into `cache` unless it is there."""
    path = cache / "phenotype_manifest.tsv.bgz"
    if refresh or not path.exists():
        cache.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(urllib.request.Request(MANIFEST, headers={"User-Agent": "gwassumstats"}), timeout=600) as response:
            path.write_bytes(response.read())
    with gzip.open(path, "rt", encoding="utf-8") as file:
        return pd.read_csv(file, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE)


def phenotypes(table: pd.DataFrame, population: str) -> pd.DataFrame:
    """Of each phenotype of the manifest that `population` has: its name, whether it is a measurement, and its cases
    and controls there (a measurement has individuals, counted as cases, and no controls)."""
    count = lambda column: pd.to_numeric(table[f"{column}_{population}"], errors="coerce").fillna(0).astype(int)
    named = table["description"].str.strip() + table["coding_description"].map(lambda text: f": {text.strip()}" if text.strip() not in ("", "NA") else "")
    p = pd.DataFrame({"name": named, "measured": table["trait_type"].isin(MEASURED), "cases": count("n_cases"), "controls": count("n_controls")})
    return p[p["cases"] > 0]


def draft(studies: pd.DataFrame, catalog: pd.DataFrame, samples: pd.DataFrame, sources: dict[str, str], phenotypes: dict[str, pd.DataFrame], same: float) -> dict[str, dict]:
    """{study: kind, cases, controls, source, how}. `studies`: anc, n, the largest N of its file's SNPs (`largest`),
    the PMID of its paper (`pmid`) and the N of the release's inventory (`released`), by ID; `catalog`: the Catalog's
    record of a study that is one; `samples`: the Catalog's records of the papers cited; `sources`: {study: the
    accession its corrected copy was taken from}; `phenotypes`: by ancestry, those of the Pan-UK Biobank; `same`: a
    sample is a study's if their N differ by no more than this share."""
    by_accession = samples.drop_duplicates("accession").set_index("accession")
    row = lambda kind="", cases="", controls="", source="", how="": {"kind": kind, "cases": cases, "controls": controls, "source": source, "how": how}
    of = lambda record, source, how: row("binary", record["cases"], record["controls"], source, how) if record["cases"] and record["controls"] else row("quantitative", source=source, how=how)
    drafted = {}
    for study, s in studies.iterrows():
        ancestry = ANCESTRIES.get(s["anc"], s["anc"])
        near = lambda size, ns=(s["n"], s["largest"]): bool(size) and ns[1] <= size * (1 + ABOVE) and any(abs(n / size - 1) <= same for n in ns)
        if study in catalog.index:
            c = catalog.loc[study]
            record = recorded(pd.Series({"accession": study, "n": c["n"], "ancestries": f"{c['ancestry']} {c['n']}", "sample": c["sample"]}), ancestry)
            drafted[study] = of(record, study, "its record in the GWAS Catalog") if bool(record["cases"]) == bool(record["controls"]) else row(source=study, how="check: its record in the GWAS Catalog has cases without controls")
            continue
        if sources.get(study) in by_accession.index:
            record = recorded(pd.Series({"accession": sources[study], **by_accession.loc[sources[study]]}), ancestry)
            whole = near(record["cases"] + record["controls"]) if record["cases"] else near(record["own"]) or near(record["n"])
            drafted[study] = of(record, sources[study], "the record of its source in the GWAS Catalog") if whole else row(source=sources[study], how="check: the record of its source in the GWAS Catalog is not of its N")
            continue
        records = [recorded(r, ancestry) for _, r in samples[samples["pmid"] == s["pmid"]].iterrows()] if s["pmid"] else []
        mine = [r for r in records if (near(r["cases"] + r["controls"]) if r["cases"] and r["controls"] else not r["cases"] and (near(r["own"]) or near(r["n"])))]
        if len({(r["cases"], r["controls"]) for r in mine}) == 1:
            drafted[study] = of(mine[0], mine[0]["accession"], "the record of its paper with its N in the GWAS Catalog")
        elif records and not any(r["cases"] for r in records):
            drafted[study] = row("quantitative", how="the GWAS Catalog counts individuals, not cases, in every record of its paper")
        elif study.startswith("PANUKBB_") and s["anc"] in phenotypes:
            p = phenotypes[s["anc"]]
            mine = p[p["cases"] + p["controls"] == s["released"]]
            if len(mine) and mine["measured"].all():
                drafted[study] = row("quantitative", source="Pan-UK Biobank manifest", how="every phenotype of its N in the Pan-UK Biobank's manifest is a measurement")
            elif len(mine) == 1:
                drafted[study] = row("binary", mine["cases"].iloc[0], mine["controls"].iloc[0], "Pan-UK Biobank manifest", f"the one phenotype of its N in the Pan-UK Biobank's manifest: {mine['name'].iloc[0]}")
            else:
                drafted[study] = row(how=f"check: {len(mine)} phenotypes of its N in the Pan-UK Biobank's manifest")
        else:
            drafted[study] = row(how="check: no record of its sample" if not records else f"check: none of the {len(records)} records of its paper is of its N" if not mine else "check: records of its N disagree")
    return drafted


def main():
    parser = argparse.ArgumentParser(prog="python -m cases", description="Drafts the kind of every study that has no row: binary, with its cases and controls, or quantitative.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--refresh", action="store_true", help="download the Pan-UK Biobank's manifest again")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    d = rel.studies(Path(args.release))
    plan = planned([row for row in cite.table(DATA / "corrections.tsv") if row["id"] in d.index])
    read = lambda name: pd.read_csv(DATA / name, sep="\t", dtype=str, keep_default_na=False)
    measured = pd.concat([read(name)[lambda t: t["file"] == "h2"].set_index("id")["n_max"].astype(float) for name in ("metrics.tsv", "metrics_corrected.tsv") if (DATA / name).exists()])
    studies = d[["anc"]].assign(n=d["n"].astype(int), released=d["n"].astype(int), pmid=[pmids[0] if pmids else "" for pmids in merge.sources(d, DATA)[0]],
                                largest=measured[~measured.index.duplicated(keep="last")].reindex(d.index).fillna(d["n"]))  # a corrected copy's, where it was measured
    for study, n in overlay.counts(plan, DATA).items():
        studies.loc[study, ["n", "largest"]] = n
    sources = {study: fix.sourced(fixes[0]["value"])["source"] for study, fixes in plan.items() if fixes[0]["fix"] == "source"}
    listed = manifest(refresh=args.refresh) if d.index.str.startswith("PANUKBB_").any() else None
    new = draft(studies, read("catalog.tsv").set_index("id"), read("samples.tsv"), sources, {anc: phenotypes(listed, anc) for anc in studies["anc"].unique()} if listed is not None else {},
                tomllib.loads((DATA / "checks.toml").read_text(encoding="utf-8"))["file"]["same_n"])

    path = DATA / "cases.tsv"
    rows = cite.table(path) if path.exists() else []
    kept = {row["id"] for row in rows}
    rows += [{"id": study, **row} for study, row in new.items() if study not in kept]
    order = {study: k for k, study in enumerate(d.index)}
    rows.sort(key=lambda row: order.get(row["id"], len(order)))
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, HEAD, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} studies, {len(rows) - len(kept)} drafted now: {sum(row['kind'] == 'binary' for row in rows)} binary, {sum(row['kind'] == 'quantitative' for row in rows)} quantitative; "
          f"{sum(row['how'].startswith('check') for row in rows)} to check")


if __name__ == "__main__":
    main()
