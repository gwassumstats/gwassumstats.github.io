"""Corrected copies: what each fix does to a file's rows, the record kept of it, and the table that asks for them."""
import json
import re

import pandas as pd
import pytest

from conftest import DATA, ROOT, panel
from corrections import fix
from corrections.__main__ import correct, planned, sha256
from papers import cite

STUDIES = {d["id"] for d in json.loads((DATA / "studies.json").read_text())["studies"]}
TEXT = "SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t100.50\t1.2300\nrs2\tC\tT\t100.50\t-0.5\nrs2\tC\tT\t90\t0.7\nrs3\tG\tA\t50.25\tinf\nrs4\tT\tC\t100.50\t0\nrs5\tA\tC\t100.50\tNA\n"


@pytest.fixture
def source(tmp_path):
    path = tmp_path / "release" / "X.sumstats"
    path.parent.mkdir()
    path.write_text(TEXT)
    return path


def test_a_flip_turns_the_sign_and_keeps_every_digit(source):
    rows, changed = fix.flip(fix.read(source).iloc[[0, 1, 4]])
    assert list(rows["Z"]) == ["-1.2300", "0.5", "0"] and changed == {"z_turned": 2} and list(rows["N"]) == ["100.50"] * 3


def test_a_clean_drops_what_is_no_number_then_what_is_listed_twice(source):
    rows, changed = fix.clean(fix.read(source))
    assert list(rows["SNP"]) == ["rs1", "rs2", "rs4"] and list(rows["Z"]) == ["1.2300", "-0.5", "0"]  # of rs2, the first row stays
    assert changed == {"z_not_a_number": 2, "listed_twice": 1}


def test_another_n_scales_each_snp_s_n(source):
    rows, changed = fix.rescale(fix.read(source), 201)  # the largest N, 100.50, becomes 201: a SNP that every cohort has
    assert list(rows["N"]) == ["201.000", "201.000", "180.000", "100.500", "201.000", "201.000"] and changed == {"n_scaled_by": 2.0, "n_largest_before": 100.5}


def test_a_copy_is_written_with_its_record_and_the_source_is_left_alone(source, tmp_path):
    out = tmp_path / "scratch"
    plan = planned([{"id": "X", "fix": "clean", "value": "", "why": "rows listed twice"}, {"id": "X", "fix": "flip", "value": "", "why": "signed for A2"},
                    {"id": "Y", "fix": "relabel", "value": "Another trait", "why": "its N"}])
    assert list(plan) == ["X"] and [row["fix"] for row in plan["X"]] == ["clean", "flip"]  # a relabel changes no file
    lines = correct(plan, {"X": str(source)}, out, "the commit")
    copy, record = out / "sumstats" / "X.sumstats", json.loads((out / "provenance" / "X.json").read_text())
    assert copy.read_text() == "SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t100.50\t-1.2300\nrs2\tC\tT\t100.50\t0.5\nrs4\tT\tC\t100.50\t0\n" and source.read_text() == TEXT
    assert record["source"] == {"path": str(source), "sha256": sha256(source), "rows": 6} and record["copy"] == {"path": str(copy), "sha256": sha256(copy), "rows": 3}
    assert record["changed"] == {"z_not_a_number": 2, "listed_twice": 1, "z_turned": 2} and record["code"] == "the commit" and [f["why"] for f in record["fixes"]] == ["rows listed twice", "signed for A2"]
    assert lines == [{"id": "X", "fixes": "clean; flip", "rows_before": 6, "rows": 3, "sha256": sha256(copy), "z_not_a_number": 2, "listed_twice": 1, "z_turned": 2}]
    with pytest.raises(ValueError, match="mend is no fix"):
        planned([{"id": "X", "fix": "mend", "value": "", "why": ""}])


def test_a_source_says_where_the_statistics_are_and_how_to_read_them():
    assert fix.sourced("GCST006906; a1=Allele1; a2 = Allele2; n=446696") == {"source": "GCST006906", "file": "", "n": 446696.0, "roles": {"a1": "Allele1", "a2": "Allele2"}}
    assert fix.sourced("https://example.org/a.txt.gz;file=x.tsv;n=Weight;") == {"source": "https://example.org/a.txt.gz", "file": "x.tsv", "n": None, "roles": {"n": "Weight"}}
    assert fix.sourced("GCST90292609") == {"source": "GCST90292609", "file": "", "n": None, "roles": {}}
    with pytest.raises(ValueError, match="harmonised of the source GCST1 is not key=value"):
        fix.sourced("GCST1; harmonised")
    with pytest.raises(ValueError, match="the source of X comes after another fix"):
        planned([{"id": "X", "fix": "flip", "value": "", "why": ""}, {"id": "X", "fix": "source", "value": "GCST1", "why": ""}])


def test_a_study_is_taken_from_its_source_then_fixed(source, tmp_path):
    theirs, out = tmp_path / "theirs.out", tmp_path / "scratch"
    theirs.write_text("MarkerName Allele1 Allele2 Effect StdErr\nrs1 a g 0.2 0.1\nrs2 t c 0.3 0.1\nrs0 a g 0.1 0.1\n")  # as its authors gave it; rs0 is not a HapMap3 SNP
    plan = planned([{"id": "X", "fix": "source", "value": f"{theirs.as_uri()}; a1=Allele1; a2=Allele2", "why": "the release's file has other alleles"},
                    {"id": "X", "fix": "flip", "value": "", "why": "to show that a fix follows"}])
    lines = correct(plan, {"X": str(source)}, out, "the commit", *panel(tmp_path), sizes={theirs.as_uri(): 5000.0})
    copy, record = out / "sumstats" / "X.sumstats", json.loads((out / "provenance" / "X.json").read_text())
    assert copy.read_text() == "SNP\tA1\tA2\tN\tZ\nrs1\tA\tG\t5000\t-2\nrs2\tT\tC\t5000\t-3\n" and source.read_text() == TEXT
    assert (out / "sources" / "by-url" / "theirs.out").read_text() == theirs.read_text()  # the download is kept
    assert record["source"]["url"] == theirs.as_uri() and record["source"]["sha256"] == sha256(theirs) and record["in_place_of"] == {"path": str(source), "sha256": sha256(source), "rows": 6}
    assert record["formatting"]["written"] == 2 and record["formatting"]["n_from"] == "the N given for the study" and record["formatting"]["z_from"] == "Effect / StdErr"
    assert record["changed"] == {"rows_of_the_source": 3, "snps_also_in_the_release_file": 2, "z_turned": 2} and record["copy"] == {"path": str(copy), "sha256": sha256(copy), "rows": 2}
    assert lines == [{"id": "X", "fixes": f"source {theirs.as_uri()}; a1=Allele1; a2=Allele2; flip", "rows_before": 6, "rows": 2, "sha256": sha256(copy), **record["changed"]}]


def test_the_table_of_corrections_names_studies_of_the_release_and_says_why():
    rows = cite.table(ROOT / "data" / "corrections.tsv")
    assert {row["id"] for row in rows} <= STUDIES and all(row["why"] and row["fix"] in {*fix.FIXES, "source", "relabel"} for row in rows)
    assert all(bool(row["value"]) == (row["fix"] in ("n", "relabel", "source")) for row in rows)
    assert all(re.fullmatch(r"GCST\d+|https://\S+", fix.sourced(row["value"])["source"]) for row in rows if row["fix"] == "source") and planned(rows)
    assert sum(row["fix"] == "flip" for row in rows) == 6 and sum(row["fix"] == "relabel" for row in rows) == 6
