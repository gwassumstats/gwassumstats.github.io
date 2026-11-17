"""Corrected copies: what each fix does to a file's rows, the record kept of it, and the table that asks for them."""
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

from conftest import DATA, ROOT, panel
from corrections import fix, ldsc
from corrections.__main__ import correct, planned, sha256, standing
from corrections.panel import selected
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


def test_a_file_of_one_pair_of_letters_gets_each_snp_s_alleles(source):
    alleles = pd.DataFrame({"A1": ["G", "C", "T"], "A2": ["A", "T", "G"]}, index=pd.Index(["rs1", "rs2", "rs4"], name="SNP"))
    rows, changed = fix.lettered(fix.read(source), alleles)
    assert rows[["SNP", "A1", "A2", "Z"]].values.tolist() == [["rs1", "G", "A", "1.2300"], ["rs2", "T", "C", "-0.5"], ["rs2", "T", "C", "0.7"], ["rs4", "T", "G", "0"]]  # the second by the alphabet first; Z as it was
    assert changed == {"alleles_from_the_list": 4, "not_in_the_list": 2}  # rs3 and rs5 are not of the list


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
    assert standing("X", plan["X"], out) and not standing("X", plan["X"][:1], out) and not standing("Z", plan["X"], out)  # there, as the table asks; other fixes; no copy


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


H2_LOG = """Reading summary statistics from X.sumstats ...
Read summary statistics for 1171029 SNPs.
Total Observed scale h2: 0.1461 (0.009)
Categories: baseL2_0 Coding_UCSCL2_0
Lambda GC: 1.2
Mean Chi^2: 1.3
Intercept: 0.7683 (0.0123)
Ratio < 0 (usually indicates GC correction).
"""
RG_LOG = """Beginning analysis ...
Computing rg for phenotype 2/4
Heritability of phenotype 1
Total Observed scale h2: 0.07 (0.004)
Intercept: 1.01 (0.008)
Genetic Covariance
Total Observed scale gencov: -0.0116 (0.005)
Intercept: 0.004 (0.006)
Genetic Correlation
Genetic Correlation: -0.1243 (0.0566)
Z-score: -2.1969
P: 0.028

Computing rg for phenotype 3/4
Total Observed scale gencov: 0.0204 (0.0057)
Genetic Correlation: nan (nan) (rg out of bounds)
Z-score: nan (nan) (rg out of bounds)
P: nan (nan) (rg out of bounds)
WARNING: One of the h2's was out of bounds.

Computing rg for phenotype 4/4
ERROR computing rg for phenotype 4/4, from file Z.sumstats.

Summary of Genetic Correlation Results
p1 p2 rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se
"""


def test_the_logs_of_ldsc_are_read_as_the_release_reads_them():
    assert ldsc.heritability(H2_LOG) == {"h2": "0.1461", "h2_se": "0.009", "h2_z": "16.23333", "intercept": "0.7683", "m": "1171029"}
    call = "Call:\n./ldsc.py \\\n--ref-ld-chr L. \\\n--out X \\\n--rg x.sumstats,a.sumstats,b.sumstats,c.sumstats \\\n--w-ld-chr W. \n\n"  # as LDSC repeats its call
    assert ldsc.listed(call + RG_LOG) == ["x.sumstats", "a.sumstats", "b.sumstats", "c.sumstats"]
    assert ldsc.correlations(call + RG_LOG) == {"a.sumstats": ["-0.1243", "0.0566", "-2.1969", "0.028", "-0.0116", "0.005"], "b.sumstats": ["nan"] * 6, "c.sumstats": ["nan"] * 6}  # out of bounds; failed
    assert ldsc.correlations(call + RG_LOG.split("Summary")[0]) == {} and ldsc.correlations(call.replace(",c.sumstats", "") + RG_LOG) == {}  # LDSC did not finish; another list of files


def test_ldsc_is_run_on_a_copy_with_the_release_s_options_and_not_again_while_its_logs_stand(tmp_path):
    (tmp_path / "h2.log").write_text(H2_LOG)
    (tmp_path / "rg.log").write_text(RG_LOG)
    stand_in = tmp_path / "ldsc.py"  # writes the log that LDSC would, after the options it was given, and counts its runs
    stand_in.write_text("import sys, pathlib\nargs = sys.argv[1:]\nhere = pathlib.Path(__file__)\nlog = here.with_name('h2.log' if '--h2' in args else 'rg.log').read_text()\n"
                        "pathlib.Path(args[args.index('--out') + 1] + '.log').write_text(' '.join(args) + '\\n' + log)\nhere.with_name('runs').open('a').write('h2\\n' if '--h2' in args else 'rg\\n')\n")
    settings = {"python": sys.executable, "ldsc": str(stand_in), "EUR": {"baseline": "B.", "frequencies": "F.", "ldscores": "L.", "weights": "W."}}
    files = {name: str(tmp_path / f"{name}.sumstats") for name in "XABCY"}
    for k, file in enumerate(files.values()):
        Path(file).write_text("SNP\n")
        os.utime(file, (1_000_000 + k, 1_000_000 + k))  # all older than any log
    scratch, partners = tmp_path / "scratch", {name: files[name] for name in "ABC"}
    runs = lambda: (tmp_path / "runs").read_text().split()
    found = ldsc.rerun("X", "EUR", files["X"], partners, scratch, settings)
    assert found == {"study": "X", "h2": ldsc.heritability(H2_LOG), "rg": {"A": ["-0.1243", "0.0566", "-2.1969", "0.028", "-0.0116", "0.005"], "B": ["nan"] * 6, "C": ["nan"] * 6}} and runs() == ["h2", "rg"]
    assert (scratch / "ldsc" / "h2" / "X.log").read_text().splitlines()[0] == f"--h2 {files['X']} --ref-ld-chr B. --w-ld-chr W. --overlap-annot --print-coefficients --frqfile-chr F. --out {scratch / 'ldsc' / 'h2' / 'X'}"
    assert (scratch / "ldsc" / "rg" / "X.log").read_text().splitlines()[0] == f"--rg {','.join(files[name] for name in 'XABC')} --ref-ld-chr L. --w-ld-chr W. --out {scratch / 'ldsc' / 'rg' / 'X'}"
    assert ldsc.rerun("X", "EUR", files["X"], partners, scratch, settings) == found and runs() == ["h2", "rg"]  # the logs stand: nothing is run again
    os.utime(files["C"])  # C's file changed after the log: its pair no longer stands
    assert ldsc.rerun("X", "EUR", files["X"], partners, scratch, settings, covered=frozenset("C"))["rg"] == found["rg"] | {"C": None} and runs() == ["h2", "rg"]  # C is corrected: its own task has the pair
    os.utime(files["C"], (1_000_003, 1_000_003))
    os.utime(files["X"])  # the study's own file changed: both logs fall
    assert ldsc.rerun("X", "EUR", files["X"], partners, scratch, settings) == found and runs() == ["h2", "rg", "h2", "rg"]
    assert ldsc.rerun("Y", "EUR", files["Y"], {}, scratch, settings)["rg"] == {} and runs()[-1] == "h2" and not (scratch / "ldsc" / "rg" / "Y.log").exists()  # no partners: h² only
    stand_in.write_text("raise SystemExit('LDSC fails')")
    with pytest.raises(ValueError, match="LDSC did not finish for X"):  # another partner, and LDSC fails
        ldsc.rerun("X", "EUR", files["X"], partners | {"Y": files["Y"]}, scratch, settings)
    Path(files["Y"]).touch()
    with pytest.raises(ValueError, match="LDSC did not finish for Y"):
        ldsc.rerun("Y", "EUR", files["Y"], {}, scratch, settings)
    assert ldsc.lacking("Y", files["Y"], {}, scratch) and not ldsc.lacking("X", files["X"], partners, scratch) and ldsc.lacking("X", files["X"], partners | {"Y": files["Y"]}, scratch)
    assert not ldsc.lacking("X", files["X"], partners | {"Y": files["Y"]}, scratch, frozenset("Y")) and ldsc.found("X", files["X"], partners | {"Y": files["Y"]}, scratch)["rg"]["Y"] is None


def test_the_panel_takes_the_strongest_study_and_leaves_what_resembles_it():
    z = {"A": 30.0, "B": 20.0, "C": 10.0, "D": 9.0, "E": 5.0, "F": 50.0}
    rg = {frozenset(pair): value for pair, value in {("A", "B"): 0.9, ("A", "C"): 0.1, ("A", "D"): -0.4, ("A", "E"): 0.0, ("B", "C"): 0.2, ("B", "D"): 0.0, ("B", "E"): 0.0,
                                                    ("C", "D"): 0.31, ("C", "E"): 0.0, ("D", "E"): 0.0, ("F", "A"): None, ("F", "B"): None, ("F", "C"): 0.0, ("F", "D"): 0.0, ("F", "E"): 0.0}.items()}
    assert selected(z, rg) == ["A", "C"]  # F, in two pairs that failed, is left out first; B and D resemble A (rg² above 0.1); E is too weak
    estimated = rg | {frozenset(("F", "A")): 0.5, frozenset(("F", "B")): 0.0}
    assert selected(z, estimated) == ["F", "B", "C", "D"]  # F stands now, and A resembles it; C and D do not resemble each other enough: 0.31² is 0.096
    assert selected(z, estimated | {frozenset(("C", "D")): 0.32}) == ["F", "B", "C"] and selected(z | {"F": 5.0}, estimated) == ["A", "C"]


def test_the_panel_rule_gives_the_release_s_own_panels(release):
    for ancestry, size in (("EUR", 85), ("EAS", 17)):
        read = lambda name: pd.read_csv(release / name, sep="\t", dtype=str, keep_default_na=False)
        h2, pairs, theirs = read(f"sldsc_heritabilities_{ancestry}_2026.txt"), read(f"ldsc_genetic_correlations_trait_pairs_{ancestry}_2026.txt"), read(f"non_redundent_traits_{ancestry}_2026.txt")
        estimated = pairs.iloc[:, 2:].map(ldsc.finite).all(axis=1)
        rg = {frozenset((a, b)): float(value) if good else None for a, b, value, good in zip(pairs["trait1"], pairs["trait2"], pairs["genetic_correlation"], estimated)}
        assert selected({study: float(z) for study, z in zip(h2["trait_identifier"], h2["h2_z"]) if ldsc.finite(z)}, rg) == list(theirs["trait_identifier"]) and len(theirs) == size


def test_the_table_of_corrections_names_studies_of_the_release_and_says_why():
    rows = cite.table(ROOT / "data" / "corrections.tsv")
    assert {row["id"] for row in rows} <= STUDIES and all(row["why"] and row["fix"] in {*fix.FIXES, "source", "relabel"} for row in rows)
    assert all(bool(row["value"]) == (row["fix"] in ("n", "relabel", "source", "alleles")) for row in rows)
    assert all(re.fullmatch(r"GCST\d+|https://\S+", fix.sourced(row["value"])["source"]) for row in rows if row["fix"] == "source") and planned(rows)
    assert sum(row["fix"] == "flip" for row in rows) == 6 and sum(row["fix"] == "relabel" for row in rows) == 6
