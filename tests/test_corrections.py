"""Corrected copies: what each fix does to a file's rows, the record kept of it, and the table that asks for them."""
import json

import pandas as pd
import pytest

from conftest import DATA, ROOT
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


def test_the_table_of_corrections_names_studies_of_the_release_and_says_why():
    rows = cite.table(ROOT / "data" / "corrections.tsv")
    assert {row["id"] for row in rows} <= STUDIES and all(row["why"] and row["fix"] in {*fix.FIXES, "relabel"} for row in rows)
    assert all(bool(row["value"]) == (row["fix"] in ("n", "relabel")) for row in rows)
    assert sum(row["fix"] == "flip" for row in rows) == 6 and sum(row["fix"] == "relabel" for row in rows) == 6
