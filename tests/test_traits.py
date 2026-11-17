"""The complex trait of each study: how one is drafted, and what the committed table says."""
import re
from collections import defaultdict

import pandas as pd

from build import merge
from conftest import DATA
from papers import cite
from traits.__main__ import category, draft, key, spelled

import json

STUDIES = json.loads((DATA / "studies.json").read_text())["studies"]


def test_a_short_name_without_its_qualifiers_is_the_key_of_its_trait():
    assert [key(name) for name in ("Kanai_2018_SBPEAS", "Venkatesh_2024_BMIFemale", "Pozarickij_2024_SBPUnadjBMIEAS", "Kany_2025_AorticStenosisMTAG", "Park_2024_MetSNoUKB",
                                   "Venkatesh_2025_FemaleInfertility", "Speliotes_2010_BMI", "Gu_2025_HbA1cFemaleEAS", "Karczewski_2025_MAPAdjMed")] \
        == ["SBP", "BMI", "SBP", "AorticStenosis", "MetS", "FemaleInfertility", "BMI", "HbA1c", "MAP"]  # a qualifier at the end only: FemaleInfertility stays


def test_a_trait_is_named_after_its_most_frequent_spelling_cleaned():
    assert spelled(["BMI (Speliotes 2010)", "Body Mass Index", "Body mass index", "BMI"]) == "Body mass index"  # two spell it out: more than write BMI once
    assert spelled(["Cigarattes per Day"]) == "Cigarettes per day" and spelled(["Long Sleep Duation"]) == "Long sleep duration"  # the release's slips
    assert spelled(["Serum urate levels"]) == "Serum urate" and spelled(["Pulse Rate Automated Reading"]) == "Pulse rate"
    assert spelled(["Alzheimer’s disease", "Alzheimers Disease"]) == "Alzheimer's disease" and spelled(["HbA1c"]) == "HbA1c"  # an abbreviation keeps its capitals


def test_traits_are_drafted_with_term_and_category():
    studies = pd.DataFrame({"name": ["Li_2024_BMI", "Li_2024_BMIFemale", "Ito_2020_BMIEAS", "Li_2024_Gout", "Roe_2019_Tanning", "Roe_2019_Oddity"],
                            "trait": ["Body mass index", "Body mass index", "BMI", "Gout", "Tanning", "An oddity"]}, index=list("abcdef"))
    catalog = pd.DataFrame({"efo": ["EFO_1", "MONDO_9", "EFO_5, EFO_6"]}, index=["a", "d", "f"])  # several terms for one study: not taken
    uri = "http://www.ebi.ac.uk/efo/"
    mappings = pd.DataFrame({"Disease trait": ["BMI", "Gout", "Skin tanning", "Tanning"], "EFO term": ["body mass index", "gout", "suntan", "suntan"],
                             "EFO URI": [uri + "EFO_1", uri + "MONDO_9", uri + "EFO_3", uri + "EFO_3"], "Parent term": ["Body measurement", "Other disease", "Other measurement", "Other measurement"]})
    rules = [("gout", "renal"), ("skin|tanning", "skin and hair"), ("Body measurement", "anthropometric")]
    got = draft(studies, catalog, mappings, rules, {"BMI": ("Body mass index", "body mass index")})
    assert got.loc["a"].to_dict() == {"trait": "Body mass index", "efo": "EFO_1", "category": "anthropometric", "how": "the EFO term of this name in the GWAS Catalog"}
    assert got.loc[["a", "b", "c"], "trait"].nunique() == 1  # one trait, whatever the sex and the ancestry
    assert got.loc["d"].to_dict() == {"trait": "Gout", "efo": "MONDO_9", "category": "renal", "how": "the EFO term of this name in the GWAS Catalog"}
    assert got.loc["e"].to_dict() == {"trait": "Tanning", "efo": "EFO_3", "category": "skin and hair", "how": "the EFO term that the GWAS Catalog maps a trait of this name to"}
    assert got.loc["f"].to_dict() == {"trait": "An oddity", "efo": "", "category": "other", "how": "check: no EFO term found; check: no category"}
    assert category("Hematological measurement", "Anything", [("Hematological", "blood")]) == "blood" and category("", "Anything", []) == "other"


def test_every_study_has_one_trait_and_a_trait_one_term_and_one_category():
    rows = cite.table(merge.DATA / "traits.tsv")
    assert [row["id"] for row in rows] == [d["id"] for d in STUDIES] and all(row["trait"] and row["category"] and row["how"] for row in rows)
    said = defaultdict(set)
    for row in rows:
        said[row["trait"]].add((row["efo"], row["category"]))
    assert not [trait for trait, values in said.items() if len(values) > 1]
    assert all(re.fullmatch(r"(EFO|MONDO|HP|OBA|GO|NCIT|Orphanet|CHEBI|UBERON)_\w+", row["efo"]) for row in rows if row["efo"])
    assert all("check" in row["how"] or "no EFO term fits: " in row["how"] for row in rows if not row["efo"])  # a trait without a term is one to look at, or says why none fits
    assert sum(not row["efo"] for row in rows) == 22 and len({row["trait"] for row in rows if not row["efo"]}) == 18
    trait = {row["id"]: row["trait"] for row in rows}
    # studies of one trait whatever their paper, ancestry, sex or adjustment
    assert len({trait[study] for study in ("PASS_BMI1", "UKB_460K.body_BMIz", "PASS_BMI_Yengo2018", "PASS_BMI_Akiyama2017_EAS", "GCST90429765")}) == 1
    assert trait["PANUKBB_LDLDirectAdjustedByMedication"] == trait["PASS_FEV1FVC_Shrine2019"]  # the trait of what the file holds, not of its label
    assert 300 < len(said) < 400 and sum(len(values) for values in said.values()) == len(said)
