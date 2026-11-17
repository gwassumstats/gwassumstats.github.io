"""Drafts, for every study that has no row, whether its sample includes participants of UK Biobank.

  uv run python -m ukb [RELEASE]

RELEASE defaults to $SUMSTATS_RELEASE. A study's sample includes UK Biobank ("yes"), or not ("no"), as:

  1. its ID says, for a study of UK Biobank alone (UKB_460K, the Pan-UK Biobank);
  2. else the cohorts of the GWAS Catalog's record of its sample say: the record of the study, of the source its
     corrected copy was taken from, or the record that data/cases.tsv takes its cases from. A list that names UK
     Biobank says yes; one that names other cohorts only says no, unless it names "other", "multiple" or "NR";
  3. else its paper's day of publication says no, if it is before May 2015, when UK Biobank first released its
     genotypes (Bycroft et al. 2018, PMID 30305743: "An interim release of genotype data on approximately 150,000 UK
     Biobank participants in May 2015").

The rest says "check": read the paper, and write "yes", "no", or nothing where it does not say, with the source
(PMID) and "by hand: " and why in `how`. A later run keeps the rows that are there. The result is data/ukb.tsv: study
ID, ukb, source and how.
"""
import argparse
import csv
import os
import re
from pathlib import Path

import pandas as pd

from build import merge, release as rel
from corrections.__main__ import planned
from corrections import fix
from papers import cite

DATA = merge.DATA
HEAD = ["id", "ukb", "source", "how"]
OWN = ("UKB_460K.", "PANUKBB_")  # the IDs of studies of UK Biobank alone
NAMED = re.compile(r"\bUKB|UK ?Biobank|\bUKBB", re.I)
VAGUE = re.compile(r"^(other|multiple|NR)$", re.I)  # a cohort list that names one of these may leave UK Biobank out
RELEASED = "2015-05"  # UK Biobank's first release of genotypes
BYCROFT = "PMID 30305743"


def said(cohorts: str) -> str:
    """What a record's list of cohorts says: "yes" if it names UK Biobank, "no" if it names other cohorts only, ""
    if it says nothing that can be read so."""
    names = [name.strip() for name in cohorts.split(",") if name.strip()]
    if any(NAMED.search(name) for name in names):
        return "yes"
    return "no" if names and not any(VAGUE.match(name) for name in names) else ""


def draft(studies: pd.DataFrame, records: dict[str, str]) -> dict[str, dict]:
    """{study: ukb, source, how}. `studies`: by ID, the record of its sample (`record`, an accession or "") and the
    day its paper was published (`published`, "" if none is known); `records`: the cohorts of each record."""
    drafted = {}
    for study, s in studies.iterrows():
        if study.startswith(OWN):
            drafted[study] = {"ukb": "yes", "source": "", "how": "a study of UK Biobank alone, by its ID"}
        elif s["record"] and (answer := said(records.get(s["record"], ""))):
            drafted[study] = {"ukb": answer, "source": s["record"], "how": f"the cohorts of its record in the GWAS Catalog: {records[s['record']]}"}
        elif s["published"] and s["published"] < RELEASED:
            drafted[study] = {"ukb": "no", "source": BYCROFT, "how": "its paper was published before UK Biobank first released genotypes, in May 2015"}
        else:
            drafted[study] = {"ukb": "", "source": "", "how": "check: no record of its cohorts"}
    return drafted


def main():
    parser = argparse.ArgumentParser(prog="python -m ukb", description="Drafts whether the sample of every study that has no row includes UK Biobank.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    d = rel.studies(Path(args.release))
    read = lambda name: pd.read_csv(DATA / name, sep="\t", dtype=str, keep_default_na=False)
    catalog, samples = read("catalog.tsv").set_index("id"), read("samples.tsv")
    records = {**samples.set_index("accession")["cohort"].to_dict(), **catalog["cohort"].to_dict()}
    plan = planned([row for row in cite.table(DATA / "corrections.tsv") if row["id"] in d.index])
    sources = {study: fix.sourced(fixes[0]["value"])["source"] for study, fixes in plan.items() if fixes[0]["fix"] == "source"}
    cases = {row["id"]: row["source"] for row in cite.table(DATA / "cases.tsv")}
    record = lambda study: study if study in catalog.index else sources.get(study) if sources.get(study) in records else cases.get(study, "") if cases.get(study, "") in records else ""
    cited, papers = merge.citations(d, DATA)
    studies = pd.DataFrame({"record": [record(study) for study in d.index], "published": [papers[pmids[0]]["published"] if pmids else "" for pmids in cited]}, index=d.index)
    new = draft(studies, records)

    path = DATA / "ukb.tsv"
    rows = cite.table(path) if path.exists() else []
    kept = {row["id"] for row in rows}
    rows += [{"id": study, **row} for study, row in new.items() if study not in kept]
    order = {study: k for k, study in enumerate(d.index)}
    rows.sort(key=lambda row: order.get(row["id"], len(order)))
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, HEAD, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} studies, {len(rows) - len(kept)} drafted now: {sum(row['ukb'] == 'yes' for row in rows)} with UK Biobank, {sum(row['ukb'] == 'no' for row in rows)} without; "
          f"{sum(row['how'].startswith('check') for row in rows)} to check")


if __name__ == "__main__":
    main()
