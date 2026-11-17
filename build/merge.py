"""What the repository adds to a release, from data/: the paper of each study, its short name, the trait it is of,
what the quality checks say of it, and its cases and controls."""
import json
import re
import tomllib
from pathlib import Path

import pandas as pd

from build import recommend
from papers import cite
from qc import checks

DATA = Path(__file__).parents[1] / "data"
NAMED = 3  # the authors of a paper that the site names before "et al."
NAME = re.compile(r"[A-Za-z]+_\d{4}_[A-Za-z0-9]+")  # a short name: FirstAuthor_Year_Trait


def sources(d: pd.DataFrame, data: Path) -> tuple[list[list[str]], dict, dict, dict]:
    """The papers each study cites, as PMIDs, the paper of its statistics first; PubMed's record of each paper; for
    each paper, how the inventory writes it for the studies whose statistics are its own; and the groups' short names.

    A study that the release gives a PMID cites that paper. Another cites what data/references.tsv says of its
    reference, for this study or for all. A row there that still says "check" is left out: its studies cite nothing yet."""
    rows = {(row["reference"], row["study"]): [pmid for pmid in row["pmid"].split(";") if pmid]
            for row in cite.table(data / "references.tsv") if "check" not in row["how"]}
    cited = [[pmid] if pd.notna(pmid) else rows.get((ref, study)) or rows.get((ref, "")) or [] for study, pmid, ref in zip(d.index, d["pmid"], d["ref"])]
    records = json.loads((data / "papers.json").read_text(encoding="utf-8"))
    missing = sorted({pmid for pmids in cited for pmid in pmids} - set(records))
    if missing:
        raise ValueError(f"data/papers.json lacks the paper {missing[0]}: run python -m papers")
    references = {pmid: [] for pmids in cited for pmid in pmids}
    for pmids, ref in zip(cited, d["ref"]):
        if pmids:
            references[pmids[0]].append(ref)
    return cited, records, references, {row["group"]: row["short"] for row in cite.table(data / "consortia.tsv")}


def citations(d: pd.DataFrame, data: Path) -> tuple[list[list[str]], dict[str, dict]]:
    """The papers each study cites (see sources), and what the site shows of every paper cited, {PMID: ...}, by PMID:
    its label, the day it was published, the first authors, the bibliographic record, and its BibTeX entry under a
    key of its own."""
    cited, records, references, groups = sources(d, data)
    keys = cite.keys({pmid: (cite.who(records[pmid], refs, groups)[0], cite.year(records[pmid], refs)) for pmid, refs in references.items()})
    papers = {}
    for pmid in sorted(references, key=int):
        paper = records[pmid]
        papers[pmid] = {"key": keys[pmid], "label": cite.label(paper, references[pmid], groups), "published": paper["published"],
                        "authors": paper["authors"][:NAMED], "more": max(0, len(paper["authors"]) - NAMED),
                        **{field: paper[field] for field in ("title", "journal", "year", "volume", "issue", "pages", "doi")}, "bib": cite.bibtex(keys[pmid], pmid, paper)}
    return cited, papers


def names(d: pd.DataFrame, data: Path) -> pd.Series:
    """The short name of each study, from data/names.tsv. Every study has one, of the form FirstAuthor_Year_Trait,
    and no two studies the same, whatever the case of its letters."""
    rows = cite.table(data / "names.tsv")
    given = pd.Series({row["id"]: row["name"] for row in rows}, dtype=object)
    if len(given) != len(rows):
        raise ValueError("data/names.tsv names a study twice")
    lacking = d.index.difference(given.index)
    if len(lacking):
        raise ValueError(f"data/names.tsv lacks the short name of {lacking[0]}: run python -m names")
    named = given.reindex(d.index)
    odd, twice = named[~named.str.fullmatch(NAME.pattern)], named[named.str.lower().duplicated()]
    if len(odd):
        raise ValueError(f"data/names.tsv: {odd.iloc[0]} is not FirstAuthor_Year_Trait")
    if len(twice):
        raise ValueError(f"data/names.tsv: {twice.iloc[0]} names two studies")
    return named


def traits(d: pd.DataFrame, data: Path) -> tuple[pd.Series, dict[str, list[str]]]:
    """The trait each study is of (data/traits.tsv); and of each trait, by name, the EFO term that fits ("" if none
    does) and its category. Every study has a trait."""
    rows = {row["id"]: row for row in cite.table(data / "traits.tsv")}
    lacking = d.index.difference(list(rows))
    if len(lacking):
        raise ValueError(f"data/traits.tsv lacks the trait of {lacking[0]}: run python -m traits")
    of = pd.Series({study: rows[study]["trait"] for study in d.index}, dtype=object)
    return of, {rows[study]["trait"]: [rows[study]["efo"], rows[study]["category"]] for study in sorted(d.index, key=lambda study: rows[study]["trait"].lower())}


def quality(d: pd.DataFrame, tables: dict, data: Path, corrected: dict[str, list] | None = None) -> pd.DataFrame:
    """What the quality checks say of each study (qc/checks.py): `qc` (pass, warn or fail) and `notes`, its findings;
    and what was measured of its file (data/metrics.tsv): `hm3`, the share it has of the HapMap3 SNPs that LDSC
    regresses on, `chi2`, the mean of Z², and `lambda`, the median of Z² over its expectation. `tables`: the genetic
    correlations of each ancestry. The GWAS Catalog's records (data/catalog.tsv) are compared if they are there. What
    a person found is counted too: a label of the release that is wrong (data/labels.tsv), and any other finding
    (data/findings.tsv).

    `corrected`: {study: its corrections as [fix, what it gives, why]}, when `d` and `tables` show the studies as
    corrected (build/overlay.py); a corrected study is then judged as corrected. Of a study with a corrected copy, the
    copy is judged: what was measured of it (data/metrics_corrected.tsv), and the sign test with the copies in
    (data/signs_corrected.tsv). A study with a corrected label is judged under that label, like one of
    data/labels.tsv: it warns, for its ID still says the other trait. Of what a person found of a corrected study,
    only what holds whatever its file stays (`of` says "study"): the rest was of the release's file, and is why the
    study was corrected."""
    corrected = corrected or {}
    copied = {study for study, fixes in corrected.items() if any(kind != "relabel" for kind, _, _ in fixes)}
    metrics = pd.read_csv(data / "metrics.tsv", sep="\t").set_index(["id", "file"])
    if copied:
        metrics = metrics[~metrics.index.get_level_values("id").isin(copied)]
        if (data / "metrics_corrected.tsv").exists():
            metrics = pd.concat([metrics, pd.read_csv(data / "metrics_corrected.tsv", sep="\t").set_index(["id", "file"])])
    known = data / "catalog.tsv"
    catalog = pd.read_csv(known, sep="\t", dtype=str, keep_default_na=False).set_index("id") if known.exists() else pd.DataFrame(columns=["n", "ancestry", "author", "published"])
    reviewed = {row["id"]: (row["label"], row["why"]) for row in cite.table(data / "labels.tsv")}
    reviewed |= {study: (value, why) for study, fixes in corrected.items() for kind, value, why in fixes if kind == "relabel"}
    by_hand = [(row["id"], row["level"], row["finding"]) for row in cite.table(data / "findings.tsv") if row["id"] not in corrected or row.get("of") == "study"]
    recorded = data / "samples.tsv"  # the samples the GWAS Catalog records for the papers cited, if python -m candidates wrote them
    samples = pd.read_csv(recorded, sep="\t", dtype=str, keep_default_na=False) if recorded.exists() else None
    cited = {study: pmids[0] for study, pmids in zip(d.index, sources(d, data)[0]) if pmids}
    tested = data / ("signs_corrected.tsv" if corrected else "signs.tsv")  # the sign test, if python -m qc.signs wrote it
    compared = pd.read_csv(tested, sep="\t") if tested.exists() else None
    judged = checks.judge(d, tables, metrics, catalog, reviewed, by_hand, tomllib.loads((data / "checks.toml").read_text(encoding="utf-8")), samples, cited,
                          compared, traits(d, data)[0].to_dict(), names(d, data).to_dict())
    file = metrics.xs("h2", level="file").reindex(d.index)
    return judged.assign(hm3=(file["in_hm3"] / file["hm3"]).clip(upper=1).round(3), chi2=file["chi2_mean"].round(3), **{"lambda": file["lambda_gc"].round(3)}, largest=file["n_max"])


def sizes(d: pd.DataFrame, largest: pd.Series, data: Path) -> pd.DataFrame:
    """Of each study, from data/cases.tsv: whether its trait is binary (`kind`: "binary", "quantitative", or "" where
    the table does not say yet); its `cases` and `controls`, missing where they are not known; the record they are
    from (`cases_from`) and what a person noted of them (`cases_how`); and `neff`, its effective N
    (build/recommend.py), missing where it is not known. `largest`: the largest N of the SNPs of each study's file.
    Every study has a row."""
    rows = {row["id"]: row for row in cite.table(data / "cases.tsv")}
    lacking = d.index.difference(list(rows))
    if len(lacking):
        raise ValueError(f"data/cases.tsv lacks {lacking[0]}: run python -m cases")
    if odd := [study for study in d.index if rows[study]["kind"] not in ("binary", "quantitative", "") or bool(rows[study]["cases"]) != bool(rows[study]["controls"])]:
        raise ValueError(f"data/cases.tsv: {odd[0]} is binary or quantitative, with cases and controls or with neither")
    number = lambda column: pd.Series([float(rows[study][column]) if rows[study][column] else float("nan") for study in d.index], index=d.index)
    said = lambda study: rows[study]["how"].removeprefix("by hand: ") if rows[study]["how"].startswith("by hand: ") and rows[study]["kind"] == "binary" else None
    t = pd.DataFrame({"kind": [rows[study]["kind"] for study in d.index], "cases": number("cases"), "controls": number("controls"),
                      "cases_from": [rows[study]["source"] if rows[study]["kind"] == "binary" and rows[study]["source"] else None for study in d.index],
                      "cases_how": [said(study) for study in d.index]}, index=d.index)
    known = [recommend.effective(kind, n, cases, controls, most if pd.notna(most) else None) for kind, n, cases, controls, most in zip(t["kind"], d["n"], t["cases"], t["controls"], largest.reindex(d.index))]
    return t.assign(neff=pd.Series(known, index=d.index, dtype=float).round())
