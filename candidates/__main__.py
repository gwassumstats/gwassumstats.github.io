"""Finds GWAS with public summary statistics that the collection lacks, in the GWAS Catalog, and backtests the rule.

  uv run python -m candidates [--refresh]

It downloads the Catalog's tables of studies and ancestries (to .cache/, once; --refresh fetches them again) and
writes, from them and the site's data:

  data/candidates.json   the candidates for the next release; the papers the rule sets aside; and the backtest: what
                         the same rule finds in the window of the last release, against what that release added
  data/catalog.tsv       the Catalog's record of each study of the collection that it knows, for the quality checks
                         and the short names: its paper, its sample, its cohorts, and whose statistics it holds
                         (women's, men's or both: from the metadata file beside its summary statistics)

The rule and the window are in data/candidates.toml.
"""
import argparse
import json
import tomllib
from pathlib import Path

import pandas as pd

from candidates import catalog

ROOT = Path(__file__).parents[1]
DATA = ROOT / "data"


def backtest(studies: pd.DataFrame, taken: pd.Index, since: str, until: str, rule: dict) -> dict:
    """What the rule finds among the studies published in the window and in the Catalog by its end, one condition
    at a time, against the studies that the release added (`taken`); what it misses, and why; what it finds beyond."""
    met = catalog.conditions(studies, rule)
    window = (studies["published"] >= since) & (studies["published"] <= until) & (studies["added"] <= until)
    found, steps = window, []
    for k, text in enumerate(catalog.STEPS):
        found = found & met[k]
        got = studies[found]
        steps.append({"rule": text.format(**rule), "studies": len(got), "papers": got["pmid"].nunique(),
                      "found": int(got.index.isin(taken).sum()), "found_papers": studies.loc[got.index.intersection(taken), "pmid"].nunique()})
    why = lambda i: ("not in the Catalog" if i not in studies.index else
                     f"published {studies.at[i, 'published'].replace('-', '/')}, outside the window" if not window[i] else
                     f"not {catalog.STEPS[int(met.loc[i].idxmin())].format(**rule)}")
    beyond = studies[found & ~studies.index.isin(taken)]
    known = set(studies.loc[studies.index.intersection(taken), "pmid"])
    return {"since": since, "until": until, "taken": len(taken), "taken_papers": len(known), "steps": steps,
            "missed": [{"id": i, "why": why(i)} for i in taken if i not in studies.index or not found[i]],
            "beyond": {"studies": len(beyond), "papers": beyond["pmid"].nunique(), "in_taken_papers": int(beyond["pmid"].isin(known).sum()),
                       "new_papers": catalog.papers(beyond[~beyond["pmid"].isin(known)])}}


def main():
    parser = argparse.ArgumentParser(prog="python -m candidates", description="Finds candidates for the next release in the GWAS Catalog.")
    parser.add_argument("--refresh", action="store_true", help="download the Catalog's tables again")
    parser.add_argument("--cache", type=Path, default=ROOT / ".cache" / "gwas-catalog", help="where the Catalog's tables are kept")
    args = parser.parse_args()
    setup = tomllib.loads((DATA / "candidates.toml").read_text())
    rule, since, until = setup["rule"], setup["since"], setup["until"]
    site = json.loads((ROOT / "site" / "data" / "studies.json").read_text())
    have = pd.Index([s["id"] for s in site["studies"]])
    taken = pd.Index([s["id"] for s in site["studies"] if s["release"] == site["release"]["year"]])

    fetched = catalog.download(args.cache, args.refresh)
    studies = catalog.load(args.cache)
    met = catalog.conditions(studies, rule)
    new = studies[(studies["added"] > until) & met[0] & ~studies.index.isin(have)]  # public, and not there for the last release
    passing = met.loc[new.index].all(axis=1)
    aside = [{"why": "more than {max_studies} studies".format(**rule), "papers": catalog.papers(new[~met.loc[new.index, 1]])},
             {"why": "another journal", "papers": catalog.papers(new[met.loc[new.index, 1] & ~met.loc[new.index, 4]])}]
    for group in aside:  # a paper set aside is named with its number of studies, not with each of them
        group["papers"] = [paper | {"studies": len(paper["studies"])} for paper in group["papers"]]
    out = {"catalog": {"fetched": fetched, "studies": len(studies), "public": int(studies["sumstats"].sum()),
                       "published_to": studies["published"].max(), "added_to": studies["added"].max()},
           "rule": rule, "since": since, "until": until,
           "candidates": catalog.papers(new[passing]), "aside": [group for group in aside if group["papers"]],
           "other": int((~passing & met.loc[new.index, 1] & met.loc[new.index, 4]).sum()),  # in a paper that passes, but another ancestry or too small
           "backtest": backtest(studies, taken, since, until, rule)}
    (DATA / "candidates.json").write_text(json.dumps(out, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
    known = studies.loc[have.intersection(studies.index)]
    known = known.assign(sex=catalog.sexes(known, args.cache / "meta"))[["pmid", "author", "published", "journal", "n", "ancestry", "sex", "trait", "sample", "cohort", "title"]]
    known.rename_axis("id").to_csv(DATA / "catalog.tsv", sep="\t")
    test = out["backtest"]
    print(f"Catalog fetched {fetched}: {len(studies):,} studies, published to {out['catalog']['published_to']}, added to {out['catalog']['added_to']}")
    print(f"candidates: {sum(len(p['studies']) for p in out['candidates'])} studies in {len(out['candidates'])} papers; set aside: "
          + "; ".join(f"{len(g['papers'])} papers ({g['why']})" for g in out["aside"]))
    for step in test["steps"]:
        print(f"backtest: {step['rule']}: {step['studies']:,} studies in {step['papers']} papers; found {step['found']}/{test['taken']} studies, {step['found_papers']}/{test['taken_papers']} papers")
    print(f"backtest: missed {len(test['missed'])}: {sorted({m['why'] for m in test['missed']})}; beyond: {test['beyond']['studies']} studies in {test['beyond']['papers']} papers")
    print(f"catalog.tsv: {len(known)} studies of the collection")


if __name__ == "__main__":
    main()
