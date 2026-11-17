"""Whether a study's sample includes UK Biobank: python -m ukb drafts it."""
import pandas as pd

from ukb.__main__ import draft, said


def test_a_record_s_cohorts_say_yes_no_or_nothing():
    assert said("UKB, FinnGen") == said("UK Biobank") == said("deCODE, UKBB") == "yes"
    assert said("BBJ, KoGES") == "no" and said("") == said("multiple") == said("NR") == said("CHARGE, other") == ""  # a list that may leave some out says nothing


def test_a_study_is_drafted_from_its_id_its_record_or_its_day():
    studies = pd.DataFrame({"record": ["", "GCST1", "GCST2", "GCST3", "", ""], "published": ["2018", "2020", "2020", "2014-01", "2014-01", "2015-06"]},
                           index=["UKB_460K.body_BMI", "with", "without", "vague_old", "old", "new"])
    got = draft(studies, {"GCST1": "UKB, HUNT", "GCST2": "BBJ", "GCST3": "multiple"})
    assert [got[study]["ukb"] for study in studies.index] == ["yes", "yes", "no", "no", "no", ""]
    assert got["old"]["source"] == "PMID 30305743" and got["with"]["source"] == "GCST1" and got["new"]["how"].startswith("check")
