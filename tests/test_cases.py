"""A study's cases and controls: how they are drafted, and what the committed table says."""
import json
import re
import tomllib

import pandas as pd

from build import merge
from build.recommend import effective
from cases.__main__ import HEAD, draft, phenotypes
from conftest import DATA
from corrections import fix
from corrections.__main__ import planned
from papers import cite

STUDIES = json.loads((DATA / "studies.json").read_text())["studies"]
SAME = 0.01


def sample(pmid, accession, n, text, ancestries=None):  # a record of the GWAS Catalog, as data/samples.tsv has it
    return {"pmid": pmid, "accession": accession, "trait": "x", "n": str(n), "ancestries": ancestries or f"European {n}", "sample": text}


def test_cases_and_controls_are_drafted_from_the_record_of_a_study_s_sample():
    study = lambda n, pmid="", largest=None, anc="EUR", released=None: {"anc": anc, "n": n, "largest": largest or n, "pmid": pmid, "released": released or n}
    studies = pd.DataFrame.from_dict({
        "GCST1": study(300), "GCST2": study(500), "GCST3": study(100),
        "taken": study(300, "p0"), "taken_other_n": study(900, "p0"),
        "of_its_n": study(300, "p1"), "of_its_largest_n": study(280, "p1", largest=400), "above_every_record": study(300, "p2", largest=330), "two_of_its_n": study(600, "p3"),
        "measured": study(777, "p4"), "east": study(150, "p5", anc="EAS"),
        "PANUKBB_One": study(1000), "PANUKBB_Measured": study(2000), "PANUKBB_Several": study(3000), "PANUKBB_Relabelled": study(1234, released=2000), "nothing": study(50, "p9")}, orient="index")
    catalog = pd.DataFrame.from_dict({"GCST1": {"n": "300", "ancestry": "European", "sample": "100 European ancestry cases, 200 European ancestry controls"},
                                      "GCST2": {"n": "500", "ancestry": "European", "sample": "500 European ancestry individuals"},
                                      "GCST3": {"n": "100", "ancestry": "European", "sample": "100 European ancestry cases"}}, orient="index")
    samples = pd.DataFrame([
        sample("p0", "GCST9", 300, "70 European ancestry cases, 50 European ancestry proxy cases, 180 European ancestry controls"),
        sample("p1", "GCST10", 300, "90 European ancestry cases, 210 European ancestry controls"), sample("p1", "GCST11", 400, "150 European ancestry cases, 250 European ancestry controls"),
        sample("p2", "GCST12", 301, "100 European ancestry cases, 201 European ancestry controls"),
        sample("p3", "GCST13", 600, "100 European ancestry cases, 500 European ancestry controls"), sample("p3", "GCST14", 600, "200 European ancestry cases, 400 European ancestry controls"),
        sample("p4", "GCST15", 800, "800 European ancestry individuals"), sample("p4", "GCST16", 900, "900 European ancestry individuals"),
        sample("p5", "GCST17", 450, "50 East Asian ancestry cases, 100 East Asian ancestry controls, 100 European ancestry cases, 200 European ancestry controls", "East Asian 150; European 300")])
    listed = pd.DataFrame({"name": ["Asthma", "Height", "Weight", "Gout", "Stroke"], "measured": [False, True, True, False, False], "cases": [100, 2000, 2000, 30, 60], "controls": [900, 0, 0, 2970, 2940]})
    got = draft(studies, catalog, samples, {"taken": "GCST9", "taken_other_n": "GCST9"}, {"EUR": listed}, SAME)
    said = lambda name: (got[name]["kind"], got[name]["cases"], got[name]["controls"], got[name]["source"], got[name]["how"])
    assert [name for name in got if got[name]["proxy"]] == ["taken"] and got["taken"]["proxy"] == "yes"  # a record that counts proxy cases, relatives of a case, says so
    # a study that is a record of the Catalog: its record's cases and controls; individuals make it quantitative; cases alone are to check
    assert said("GCST1") == ("binary", 100, 200, "GCST1", "its record in the GWAS Catalog") and said("GCST2") == ("quantitative", "", "", "GCST2", "its record in the GWAS Catalog")
    assert said("GCST3") == ("", "", "", "GCST3", "check: its record in the GWAS Catalog has cases without controls")
    # a corrected study: the record of the source its copy was taken from, if it is of its N
    assert said("taken") == ("binary", 120, 180, "GCST9", "the record of its source in the GWAS Catalog") and got["taken_other_n"]["how"] == "check: the record of its source in the GWAS Catalog is not of its N"
    # another: the record of its paper with its N, the inventory's or the largest of its file's SNPs; no record that its file's largest N exceeds
    assert said("of_its_n") == ("binary", 90, 210, "GCST10", "the record of its paper with its N in the GWAS Catalog") and said("of_its_largest_n")[:4] == ("binary", 150, 250, "GCST11")
    assert got["above_every_record"]["how"] == "check: none of the 1 records of its paper is of its N" and got["two_of_its_n"]["how"] == "check: records of its N disagree"
    assert said("east")[:4] == ("binary", 50, 100, "GCST17")  # of a sample of several ancestries, the cases and controls of the study's own
    assert said("measured") == ("quantitative", "", "", "", "the GWAS Catalog counts individuals, not cases, in every record of its paper")
    # a study of the Pan-UK Biobank: the phenotype of its N in the manifest, if there is one only, or if all are measurements; by the N of the release, whatever a correction made of it
    assert said("PANUKBB_One") == ("binary", 100, 900, "Pan-UK Biobank manifest", "the one phenotype of its N in the Pan-UK Biobank's manifest: Asthma")
    assert said("PANUKBB_Measured")[:2] == ("quantitative", "") and said("PANUKBB_Relabelled")[0] == "quantitative" and got["PANUKBB_Several"]["how"] == "check: 2 phenotypes of its N in the Pan-UK Biobank's manifest"
    assert got["nothing"] == {"kind": "", "cases": "", "controls": "", "proxy": "", "source": "", "how": "check: no record of its sample"}


def test_the_manifest_s_phenotypes_are_read_for_an_ancestry():
    table = pd.DataFrame({"trait_type": ["phecode", "continuous", "biomarkers", "categorical"], "description": ["Gout", "Height", "Urate", "Operation code"], "coding_description": ["NA", "", "NA", "appendicectomy"],
                          "n_cases_EUR": ["30", "2000", "1500", "NA"], "n_controls_EUR": ["970", "NA", "NA", "NA"], "n_cases_EAS": ["3", "NA", "20", "5"], "n_controls_EAS": ["97", "NA", "NA", "95"]})
    assert phenotypes(table, "EUR").to_dict("records") == [{"name": "Gout", "measured": False, "cases": 30, "controls": 970}, {"name": "Height", "measured": True, "cases": 2000, "controls": 0},
                                                           {"name": "Urate", "measured": True, "cases": 1500, "controls": 0}]  # a measurement has individuals, and no controls
    assert list(phenotypes(table, "EAS")["name"]) == ["Gout", "Urate", "Operation code: appendicectomy"]  # those that the ancestry has; an answer or a code is named with its field


def test_every_study_is_binary_or_quantitative_and_a_binary_study_s_cases_are_from_a_record():
    rows = cite.table(merge.DATA / "cases.tsv")
    by_id, by_study = {d["id"]: d for d in STUDIES}, {row["id"]: row for row in rows}
    assert [row["id"] for row in rows] == [d["id"] for d in STUDIES] and list(rows[0]) == HEAD
    assert all(row["kind"] in ("binary", "quantitative") and row["how"] and "check" not in row["how"] for row in rows)  # nothing left to look at
    counted = [row for row in rows if row["cases"] or row["controls"]]
    assert all(row["kind"] == "binary" and re.fullmatch(r"[1-9]\d*", row["cases"]) and re.fullmatch(r"[1-9]\d*", row["controls"]) and row["source"] for row in counted)
    assert all(re.fullmatch(r"|GCST\d+|PMID \d+|Pan-UK Biobank manifest", row["source"]) for row in rows)
    # a record named is one of the committed tables, and a paper named is the study's own: whoever doubts a number can look it up
    recorded = {row["accession"] for row in cite.table(merge.DATA / "samples.tsv")} | {row["id"] for row in cite.table(merge.DATA / "catalog.tsv")}
    assert not [row["id"] for row in rows if row["source"].startswith("GCST") and row["source"] not in recorded]
    assert not [row["id"] for row in rows if row["source"].startswith("PMID ") and row["source"][5:] not in by_id[row["id"]]["cite"]]
    assert all(row["source"] == "Pan-UK Biobank manifest" for row in rows if row["id"].startswith("PANUKBB_") and row["kind"] == "binary" and row["cases"])
    # a binary study without cases says why; those of UK Biobank's 460K are not public
    without = [row for row in rows if row["kind"] == "binary" and not row["cases"]]
    assert len(without) == 37 and all(row["how"].startswith("by hand: ") and re.search(r"not (public|known)|in no record read", row["how"]) for row in without)
    assert sum("Its cases are not public" in row["how"] for row in without) == 28 == sum(row["id"].startswith("UKB_460K.") for row in without)
    # the site shows what the table says, and the effective N that follows from it
    read = lambda name: pd.read_csv(merge.DATA / name, sep="\t").query("file == 'h2'").set_index("id")["n_max"]
    largest = pd.concat([read("metrics.tsv"), read("metrics_corrected.tsv")]).groupby(level=0).last()
    for d in STUDIES:
        row = by_study[d["id"]]
        assert d.get("kind") == (row["kind"] if row["kind"] == "binary" else None) and (str(d.get("cases", "")), str(d.get("controls", ""))) == (row["cases"], row["controls"])
        assert d.get("cases_from") == (row["source"] or None if row["kind"] == "binary" else None)
        expected = effective(row["kind"], d["n"], float(row["cases"] or "nan"), float(row["controls"] or "nan"), largest[d["id"]])
        assert d.get("neff") == (round(expected) if row["kind"] == "binary" and expected is not None else None)
        assert "neff" not in d or d["neff"] <= max(d["n"], largest[d["id"]]) and d["neff"] <= 4 / (1 / d["cases"] + 1 / d["controls"]) + 1  # no more than its individuals, nor than its cases and controls allow
    assert sum(d.get("kind") == "binary" for d in STUDIES) == 244 and sum("neff" in d for d in STUDIES) == 207


def test_a_row_that_no_person_wrote_is_what_the_draft_gives():
    """The rows drafted from the committed tables are drafted again from them, and are the same: a table that changed
    under a row shows here. The rows read from the Pan-UK Biobank's manifest, which is not committed, are left out."""
    rows = [row for row in cite.table(merge.DATA / "cases.tsv") if not row["how"].startswith("by hand") and "Pan-UK Biobank" not in row["how"]]
    plan = planned(cite.table(merge.DATA / "corrections.tsv"))
    sources = {study: fix.sourced(fixes[0]["value"])["source"] for study, fixes in plan.items() if fixes[0]["fix"] == "source"}
    read = lambda name: pd.read_csv(merge.DATA / name, sep="\t").query("file == 'h2'").set_index("id")["n_max"]
    largest = pd.concat([read("metrics.tsv"), read("metrics_corrected.tsv")]).groupby(level=0).last()
    studies = pd.DataFrame({"anc": [d["anc"] for d in STUDIES], "n": [d["n"] for d in STUDIES], "pmid": [d["cite"][0] for d in STUDIES],
                            "released": [d.get("was", {}).get("n", d["n"]) for d in STUDIES]}, index=[d["id"] for d in STUDIES])
    studies["largest"] = [d["n"] if "n" in d.get("was", {}) or d["id"] in sources and d["id"] not in largest.index else largest[d["id"]] for d in STUDIES]
    table = lambda name: pd.read_csv(merge.DATA / name, sep="\t", dtype=str, keep_default_na=False)
    again = draft(studies, table("catalog.tsv").set_index("id"), table("samples.tsv"), sources, {}, tomllib.loads((merge.DATA / "checks.toml").read_text())["file"]["same_n"])
    assert len(rows) > 400 and not [(row["id"], again[row["id"]]) for row in rows if {key: str(value) for key, value in again[row["id"]].items()} != {key: row[key] for key in HEAD[1:]}]
