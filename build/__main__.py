"""Builds the site's data from a release folder.

  SUMSTATS_PASSPHRASE=... uv run python -m build RELEASE

RELEASE defaults to $SUMSTATS_RELEASE; without $SUMSTATS_PASSPHRASE the build asks for the passphrase. Neither is
ever written to the repository. The build writes, in site/data/:

  studies.json         the inventory with h², each study's short name, the paper it cites,
                       what the quality checks say of it, its cases and controls, and
                       whether it is the study to use of its trait; whether its sample
                       includes UK Biobank; the release's counts                         public
  papers.json          every paper cited: authors, title, journal, DOI, BibTeX entry     public
  candidates.json      the candidates for the next release, and the rule's backtest      public
  rg_{ancestry}.json   every genetic correlation of an ancestry                          public
  paths.enc.json       where each study's file is, encrypted (crypt.py)                  readable after login

The papers are those of data/references.tsv and data/papers.json (python -m papers); the short names, those of
data/names.tsv (python -m names). A study without a short name stops the build. The quality checks (qc/checks.py)
judge each study from its file as measured (data/metrics.tsv: python -m qc.scan), its results and its label; a
corrected study, as corrected (its copy as measured: data/metrics_corrected.tsv), with what they say of the release's
file kept beside. A study's cases and controls are those of data/cases.tsv (python -m cases); of the studies of a
trait and an ancestry that pass, the one of the largest effective N is recommended, unless one whose effective N is
not known has a larger h² z (build/recommend.py). The
candidates are those of data/candidates.json (python -m candidates), as they are.

A release that has not changed rewrites the same files, and keeps the encrypted one as it is, so readers stay logged in.
"""
import argparse
import getpass
import gzip
import json
import os
import re
from pathlib import Path

from cryptography.exceptions import InvalidTag

import pandas as pd

from build import crypt, liability, merge, overlay, recommend, release as rel

SITE_DATA = Path(__file__).parents[1] / "site" / "data"
BLOB = "paths.enc.json"


def dumps(value) -> str:
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"))


def holds(path: Path, plain: bytes, passphrase: str, iterations: int) -> bool:
    """Whether the blob at `path` already encrypts `plain` under this passphrase."""
    try:
        blob = json.loads(path.read_text())
        return blob["iterations"] == iterations and crypt.decrypt(blob, passphrase) == plain
    except (FileNotFoundError, InvalidTag, KeyError, ValueError):
        return False


def build(release: Path, out: Path, passphrase: str, iterations: int = crypt.ITERATIONS, data: Path = merge.DATA, copies: Path | None = None) -> list[Path]:
    """Writes the site's data for the release at `release` into `out`, with what `data` adds to it; returns the files.
    The release's own counts are of the release as it is; the studies and their genetic correlations are shown with
    the corrections laid over them (build/overlay.py), and checked so: what the checks say of a corrected study as
    released is kept in its `was`. `copies`: the folder of the corrected copies, for the paths that only a logged-in
    reader sees."""
    theirs = rel.studies(release)
    pairs = {anc: rel.pairs_of(release, anc) for anc in theirs["anc"].unique()}
    as_released = {anc: rel.correlations(table, theirs, anc) for anc, table in pairs.items()}
    cited, papers = merge.citations(theirs, data)
    judged = merge.quality(theirs, as_released, data)
    of, traits = merge.traits(theirs, data)
    d, tables, fixed, was = overlay.laid(theirs, pairs, data, copies)
    rg = {anc: rel.correlations(table, d, anc) for anc, table in tables.items()}
    released = judged  # what the checks say of the release as it is; then of the studies as shown: the corrected ones as corrected
    if was:
        judged = merge.quality(d, rg, data, fixed)
        for study in fixed:
            if (released.at[study, "qc"], released.at[study, "notes"]) != (judged.at[study, "qc"], judged.at[study, "notes"]):
                was[study] = was.get(study, {}) | {"qc": released.at[study, "qc"], "notes": released.at[study, "notes"]}
    for anc, t in rg.items():  # the release's rg and its standard error, for the pairs estimated again: by the pair's position
        place, n = {study: k for k, study in enumerate(t["ids"])}, len(t["ids"])
        before = {min(place[a], place[b]) * (2 * n - min(place[a], place[b]) - 1) // 2 + abs(place[a] - place[b]) - 1: value for (a, b), value in ((pair.split(" "), value) for pair, value in was.get(anc, {}).items())}
        if before:
            t["was"] = {str(position): before[position] for position in sorted(before)}
    named = "|".join(re.escape(study) for study in sorted(d.index, key=len, reverse=True))
    said = lambda text: re.sub(rf"(?<![\w.])({named})(?![\w.])", r"{\1}", text)  # a study that a correction's reason names, as the site links it
    name, sized = merge.names(d, data), merge.sizes(d, judged["largest"], data)
    published = pd.Series([papers[pmids[0]]["published"] if pmids else "" for pmids in cited], index=d.index)  # of the paper of a study's statistics
    recommended, superseded = recommend.chosen(pd.DataFrame({"of": of, "anc": d["anc"], "qc": judged["qc"], "n": d["n"], "neff": sized["neff"], "h2_z": d["h2_z"], "name": name, "published": published}))
    binary = sized["kind"] == "binary"  # a quantitative study's effective N is its N: only a binary study's cases, controls and effective N are written
    counted = sized[["cases", "controls", "neff"]].where(binary).astype("Int64")
    prevalence = liability.prevalences(data)
    scaled = pd.DataFrame([liability.scaled(kind, proxy, cases, controls, n, h2, se, prevalence.get((trait, anc), {}).get("k"))
                           for kind, proxy, cases, controls, n, h2, se, trait, anc in zip(sized["kind"], sized["proxy"], sized["cases"], sized["controls"], d["n"], d["h2"], d["h2_se"], of, d["anc"])],
                          index=d.index, columns=["h2l", "h2l_se"]).astype(float).round(5)
    biobank = merge.biobank(d, data)
    public = (d.assign(name=name, of=of).join(judged[["qc", "hm3", "chi2", "lambda"]]).join(counted).join(scaled).join(biobank)
              .assign(kind=sized["kind"].where(binary), proxy=sized["proxy"].where(sized["proxy"]), cases_from=sized["cases_from"], cases_how=sized["cases_how"], rec=recommended, sup=superseded)
              [["name", "of", *rel.PUBLIC, "h2l", "h2l_se", "qc", "hm3", "chi2", "lambda", "kind", "cases", "controls", "neff", "proxy", "cases_from", "cases_how", "ukb", "ukb_from", "ukb_how", "rec", "sup"]].reset_index())
    records = public.astype(object).where(public.notna(), None).to_dict("records")
    records = [{k: v for k, v in r.items() if v is not None} | ({"cite": pmids} if pmids else {}) | ({"notes": notes} if notes else {})
               | ({"fixed": [[kind, value, said(why)] for kind, value, why in fixed[r["id"]]]} if r["id"] in fixed else {}) | ({"was": was[r["id"]]} if r["id"] in was else {})
               for r, pmids, notes in zip(records, cited, judged["notes"])]
    summary = rel.summary(release, theirs, as_released) | {"qc": {level: int((judged["qc"] == level).sum()) for level in ("pass", "warn", "fail")},
                                                          "traits": {"all": len(traits), "efo": sum(bool(efo) for efo, _ in traits.values())},
                                                          "binary": {"studies": int(binary.sum()), "counted": int((binary & sized["neff"].notna()).sum()),
                                                                     "liability": int(scaled["h2l"].notna().sum()), "proxy": int((binary & sized["proxy"]).sum())},
                                                          "ukb": {"yes": int((biobank["ukb"] == True).sum()), "no": int((biobank["ukb"] == False).sum()), "unknown": int(biobank["ukb"].isna().sum())},
                                                          "recommended": {"of": int(pd.DataFrame({"of": of, "anc": d["anc"]}).drop_duplicates().shape[0]), "studies": int(recommended.notna().sum()),
                                                                          "by_z": int((recommended == "z").sum()), "superseded": int(superseded.notna().sum())}}
    if fixed:
        summary["corrected"] = {"studies": len(fixed), "copies": sum(any(kind != "relabel" for kind, _, _ in fixes) for fixes in fixed.values()),
                                "qc": {level: int((released["qc"] == level).sum()) for level in ("pass", "warn", "fail")},  # of the release's files
                                "signed": int((d["cross"] & ~theirs["cross"]).sum()), "pairs": sum(len(was.get(anc, {})) for anc in rg),
                                "panel": {anc: {"out": list(theirs.index[(theirs["anc"] == anc) & theirs["panel"].notna() & d["panel"].isna()]),
                                                "in": list(d.index[(d["anc"] == anc) & d["panel"].notna() & theirs["panel"].isna()])} for anc in rg}}
    lines = lambda items: ",\n".join(items)  # a study, a paper a line, for readable diffs
    cited_as = {pmid: {"label": p["label"], "published": p["published"]} for pmid, p in papers.items()}  # what the table of studies needs of a paper
    used = {key: row for key, row in prevalence.items() if key in set(zip(of[binary], d["anc"][binary]))}  # the prevalence of each binary trait and ancestry that has one: [K, its paper, its words]
    stated = {trait: {anc: [row["k"], row["source"].removeprefix("PMID "), row["quote"]] for (t, anc), row in used.items() if t == trait} for trait, _ in used}
    if unknown := [row["source"] for row in used.values() if row["source"].removeprefix("PMID ") not in papers]:
        raise ValueError(f"data/prevalence.tsv: {unknown[0]} is not a paper of the collection")
    texts = {"studies.json": f'{{\n"release":{dumps(summary)},\n"papers":{dumps(cited_as)},\n"prevalence":{dumps(stated)},\n"traits":{{\n' + lines(f"{dumps(trait)}:{dumps(term)}" for trait, term in traits.items()) + '\n},\n"studies":[\n'
                             + lines(map(dumps, records)) + "\n]}\n",
             "papers.json": "{\n" + lines(f"{dumps(pmid)}:{dumps(paper)}" for pmid, paper in papers.items()) + "\n}\n",
             "candidates.json": dumps(json.loads((data / "candidates.json").read_text(encoding="utf-8"))) + "\n"}
    for anc, t in rg.items():
        texts[f"rg_{anc}.json"] = "{\n" + ",\n".join(f"{dumps(k)}:{dumps(v)}" for k, v in t.items()) + "\n}\n"
    out.mkdir(parents=True, exist_ok=True)
    for name, text in texts.items():
        (out / name).write_text(text, encoding="utf-8")
    plain = dumps(rel.private(release, d)).encode()
    if not holds(out / BLOB, plain, passphrase, iterations):
        (out / BLOB).write_text(json.dumps(crypt.encrypt(plain, passphrase, iterations), indent=1) + "\n")
    return [out / name for name in [*texts, BLOB]]


def main():
    parser = argparse.ArgumentParser(prog="python -m build", description="Builds the site's data from a release folder.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--out", type=Path, default=SITE_DATA, help="where to write (default: site/data)")
    parser.add_argument("--data", type=Path, default=merge.DATA, help="what the repository adds to the release (default: data)")
    parser.add_argument("--copies", type=Path, default=os.environ.get("SUMSTATS_CORRECTIONS"), help="the folder of the corrected copies, for their paths (default: $SUMSTATS_CORRECTIONS)")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    passphrase = os.environ.get("SUMSTATS_PASSPHRASE") or getpass.getpass("Passphrase: ")
    if not passphrase:
        parser.error("the passphrase is empty")
    for path in build(Path(os.path.abspath(args.release)), args.out, passphrase, data=args.data, copies=args.copies and Path(os.path.abspath(args.copies))):
        size = path.stat().st_size
        print(f"{path.name:<16} {size / 1e3:>8,.0f} kB  ({len(gzip.compress(path.read_bytes(), 6)) / 1e3:,.0f} kB gzipped)")


if __name__ == "__main__":
    main()
