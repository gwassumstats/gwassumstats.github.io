"""The candidates for the next release: the rule on a made-up catalog, its backtest, and what the site shows of both."""
import json
import tomllib

import pandas as pd
import pytest

from candidates import catalog
from candidates.__main__ import backtest
from conftest import DATA, ROOT, write

RULE = {"journals": ["Nat Genet", "Nature"], "max_studies": 2, "ancestries": ["European", "East Asian"], "min_n": 50_000}


def study(i, pmid, published, journal="Nat Genet", n=100_000, ancestry="European", sumstats=True, added=None):
    return {"id": i, "pmid": pmid, "author": f"Author{pmid} A", "published": published, "added": added or published, "journal": journal, "title": f"Paper {pmid}", "trait": f"Trait {i}",
            "sample": "", "cohort": "", "sumstats": sumstats, "location": "", "n": n, "ancestry": ancestry}


@pytest.fixture
def made_up():
    """A catalog of five papers in the window of a release, one before it and two after."""
    return pd.DataFrame([
        study("G1", "10", "2024-05-01"), study("G2", "10", "2024-05-01", n=40_000),  # one large enough, one too small
        study("G3", "11", "2025-01-01", journal="Some J"),  # another journal
        study("G4", "12", "2025-06-01"), study("G5", "12", "2025-06-01"), study("G6", "12", "2025-06-01"),  # three studies: more than the rule takes from a paper
        study("G7", "13", "2026-01-01", ancestry="African American or Afro-Caribbean, European"), study("G8", "13", "2026-01-01", sumstats=False),
        study("G9", "14", "2023-12-31"),  # before the window
        study("G10", "15", "2026-02-01", added="2026-10-01"), study("G11", "16", "2026-10-02", ancestry="East Asian"),  # in the Catalog after the window's end
        study("G12", "17", "2025-03-01"),  # nothing against it
    ]).set_index("id")


def test_the_catalog_s_tables_are_read_into_a_row_a_study(tmp_path):
    columns = ["STUDY ACCESSION", "PUBMED ID", "FIRST AUTHOR", "DATE", "DATE ADDED TO CATALOG", "JOURNAL", "STUDY", "DISEASE/TRAIT", "INITIAL SAMPLE SIZE", "COHORT",
               "FULL SUMMARY STATISTICS", "SUMMARY STATS LOCATION", "MAPPED_TRAIT", "MAPPED_TRAIT_URI"]
    write(tmp_path / "studies.tsv", columns, [
        ["G1", "10", "Kim A", "2024-05-01", "2024-06-01", "Nat Genet", "A paper", "Trait 1", "1,500 European ancestry cases, 500 East Asian ancestry cases, 40 cases", "UKB|BBJ", "yes", "http://example.org/G1", "trait one", "http://www.ebi.ac.uk/efo/EFO_1"],
        ["G2", "10", "Kim A", "2024-05-01", "2024-06-01", "Nat Genet", "A paper", "Trait 2", "700 individuals", "", "no", "NA", "", ""]])
    write(tmp_path / "ancestry.tsv", ["STUDY ACCESSION", "STAGE", "NUMBER OF INDIVDUALS", "BROAD ANCESTRAL CATEGORY"], [
        ["G1", "initial", 1500, "European"], ["G1", "initial", 500, "East Asian"], ["G1", "replication", 9000, "African unspecified"], ["G1", "initial", 40, "NR"],
        ["G2", "initial", 400, "European, NR"], ["G2", "initial", 300, "European, NR"]])
    table = catalog.load(tmp_path)
    assert table.loc["G1", ["pmid", "n", "ancestry", "ancestries", "sumstats", "cohort", "efo"]].tolist() == ["10", 2040, "East Asian, European", "European 1500; East Asian 500; NR 40", True, "UKB, BBJ", "EFO_1"]
    assert table.loc["G2", ["n", "ancestry", "ancestries", "sumstats"]].tolist() == [700, "European", "European, NR 700", False]  # the first stage only; a group as the Catalog names it


def test_the_rule_is_five_conditions(made_up):
    met = catalog.conditions(made_up, RULE)
    assert list(met.columns) == list(range(len(catalog.STEPS)))
    assert {i: [k for k in met.columns if not met.at[i, k]] for i in made_up.index if not met.loc[i].all()} == {
        "G2": [3], "G3": [4], "G4": [1], "G5": [1], "G6": [1], "G7": [2], "G8": [0]}  # too small; another journal; too many in a paper; two ancestries; no summary statistics
    assert met.at["G7", 1] and met.at["G8", 1]  # a paper is counted by its studies with summary statistics: one here


def test_papers_group_their_studies_the_latest_first(made_up):
    papers = catalog.papers(made_up.loc[["G1", "G2", "G11"]])
    assert [(p["pmid"], p["author"], p["published"], [s["id"] for s in p["studies"]]) for p in papers] == [("16", "Author16 A", "2026-10-02", ["G11"]), ("10", "Author10 A", "2024-05-01", ["G1", "G2"])]
    assert papers[0]["studies"][0] == {"id": "G11", "trait": "Trait G11", "n": 100_000, "ancestry": "East Asian"}


def test_the_backtest_counts_what_the_rule_finds_misses_and_finds_beyond(made_up):
    taken = pd.Index(["G1", "G2", "G9", "G0"])  # what a release added: one the rule finds, one too small, one published before the window, one that is not in the Catalog
    test = backtest(made_up, taken, "2024-01-01", "2026-08-26", RULE)
    assert (test["taken"], test["taken_papers"]) == (4, 2)
    assert [(step["studies"], step["papers"], step["found"], step["found_papers"]) for step in test["steps"]] == [(8, 5, 2, 1), (5, 4, 2, 1), (4, 3, 2, 1), (3, 3, 1, 1), (2, 2, 1, 1)]
    assert test["steps"][1]["rule"] == "in a paper of at most 2 such studies" and test["steps"][3]["rule"] == "of at least 50,000 individuals"
    assert test["missed"] == [{"id": "G2", "why": "not of at least 50,000 individuals"}, {"id": "G9", "why": "published 2023/12/31, outside the window"}, {"id": "G0", "why": "not in the Catalog"}]
    assert (test["beyond"]["studies"], test["beyond"]["papers"], test["beyond"]["in_taken_papers"]) == (1, 1, 0)  # G12, which the release did not take
    assert [(p["pmid"], [s["id"] for s in p["studies"]]) for p in test["beyond"]["new_papers"]] == [("17", ["G12"])]


def test_a_study_s_sex_is_read_from_its_metadata_file(made_up, tmp_path):
    for name, text in (("G1", "gwas_id: G1\nsex: F\n"), ("G2", "sex: combined\n"), ("G3", "sex: NR\n"), ("G4", "trait_description:\n  - sex: M\n"), ("G5", "")):
        (tmp_path / f"{name}.yaml").write_text(text)
    assert list(catalog.sexes(made_up.loc[["G1", "G2", "G3", "G4", "G5"]], tmp_path)) == ["F", "combined", "", "", ""]  # as told; not reported; another key's line; no file


# ---------------------------------------------------------------- the committed result


def test_the_samples_are_those_of_the_papers_cited_and_those_that_a_study_s_cases_are_from():
    samples = pd.read_csv(ROOT / "data" / "samples.tsv", sep="\t", dtype=str, keep_default_na=False)
    named = set(pd.read_csv(ROOT / "data" / "cases.tsv", sep="\t", dtype=str, keep_default_na=False)["source"])  # a record of another paper may hold a study's sample
    cited = samples["pmid"].isin(set(json.loads((ROOT / "data" / "papers.json").read_text())))
    assert (cited | samples["accession"].isin(named)).all() and (~cited).sum() == 12 and not samples["accession"].duplicated().any()
    assert samples.groupby("pmid").size().max() <= tomllib.loads((ROOT / "data" / "candidates.toml").read_text())["rule"]["max_studies"]  # among more, one always fits
    by_ancestry = samples["ancestries"].str.findall(r" (\d+)(?:;|$)").map(lambda counts: sum(map(int, counts)))
    assert (by_ancestry == samples["n"].astype(int)).all() and (samples["n"].astype(int) > 0).all()  # the individuals of each ancestry add up to the sample


def test_the_site_shows_the_candidates_that_the_script_found():
    found, shown = json.loads((ROOT / "data" / "candidates.json").read_text()), json.loads((DATA / "candidates.json").read_text())
    assert shown == found
    test, have = found["backtest"], {d["id"]: d for d in json.loads((DATA / "studies.json").read_text())["studies"]}
    assert test["taken"] == sum(d["release"] == 2026 for d in have.values()) == 205 and test["taken_papers"] == 35
    last = test["steps"][-1]
    assert last["found"] == test["taken"] - len(test["missed"]) and last["found"] >= 0.95 * test["taken"]  # the rule finds nearly all that the release added
    assert all(m["id"] in have and have[m["id"]]["release"] == 2026 for m in test["missed"])
    assert [step["studies"] for step in test["steps"]] == sorted((step["studies"] for step in test["steps"]), reverse=True)  # each condition narrows
    assert not [s["id"] for paper in found["candidates"] for s in paper["studies"] if s["id"] in have]  # a candidate is not in the collection yet
    assert all(paper["added"] > found["until"] for paper in found["candidates"])
