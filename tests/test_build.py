"""The build, on a made-up release (conftest.mini): what it writes, what it refuses, and the encryption."""
import json

import pytest
from cryptography.exceptions import InvalidTag

import numpy as np

from build import crypt
from build.__main__ import BLOB, build
from build.release import tree
from conftest import NAN, RG, write

FAST = 1000  # PBKDF2 rounds, for tests that do not need the real cost


def built(release, out, passphrase="correct horse"):
    build(release, out, passphrase, FAST)
    return {path.name: json.loads(path.read_text()) for path in out.iterdir()}


def test_datasets_carry_the_inventory_h2_panel_and_paper(mini, tmp_path):
    data = built(mini, tmp_path / "out")["datasets.json"]
    old, *_, new = data["datasets"][:4]
    assert [d["id"] for d in data["datasets"]] == ["A_unsigned", "B", "C", "GCST0001", "D_EAS", "GCST0002", "GCST0003"]
    assert old == {"id": "A_unsigned", "trait": "Height (old)", "anc": "EUR", "n": 1000, "m": 900, "origin": "Old 2015", "ref": "Smith et al., 2010 Nature",
                   "signed": False, "alleles": True, "cross": False, "release": 2024, "h2": 0.2, "h2_se": 0.02, "h2_z": 10.0, "intercept": 1.01}
    assert new == {"id": "GCST0001", "trait": "Trait new", "anc": "EUR", "n": 4000, "m": 930, "origin": "New 2026", "ref": "Ortín et al. 2025 Nat Genet",
                   "signed": True, "alleles": True, "cross": True, "release": 2026, "h2": 0.1461, "h2_se": 0.009, "h2_z": 16.23333, "intercept": 0.7683,
                   "panel": 1, "m_rg": 926, "pmid": "111", "paper": "https://pubmed.ncbi.nlm.nih.gov/111/",
                   "catalog": "https://www.ebi.ac.uk/gwas/studies/GCST0001", "paper_id": "NG_Ortin_2025_x"}
    assert {d["id"]: d["panel"] for d in data["datasets"] if "panel" in d} == {"GCST0001": 1, "C": 2, "GCST0002": 1}  # the panels' own order
    assert data["datasets"][1]["trait"] == "Alzheimer’s disease" and data["datasets"][1]["h2"] == -0.0058


def test_release_counts(mini, tmp_path):
    release = built(mini, tmp_path / "out")["datasets.json"]["release"]
    assert release == {
        "year": 2026, "packaged": "2026-09-07", "papers": 2, "ancestries": ["EUR", "EAS"],
        "all": {"datasets": 7, "before": 4, "added": 3, "signed": 6, "signed_before": 3, "unsigned": 1, "unsigned_added": 0, "no_alleles": 0},
        "EUR": {"datasets": 4, "before": 3, "added": 1, "signed": 3, "signed_before": 2, "pairs": 3, "failed": 1,
                "worst": {"id": "C", "failed": 1, "pairs": 2}, "panel": 2, "panel_added": 1, "pairs_before": 1, "panel_before": 2, "panel_kept": 1},
        "EAS": {"datasets": 3, "before": 1, "added": 2, "signed": 3, "signed_before": 1, "pairs": 3, "failed": 1,
                "worst": {"id": "GCST0002", "failed": 1, "pairs": 2}, "panel": 1, "panel_added": 1, "pairs_before": 0}}


def test_pairs_are_stored_by_position_whatever_their_order_in_the_release(mini, tmp_path):
    out = built(mini, tmp_path / "out")
    assert out["rg_EUR.json"] == {  # the pairs of ids, row by row: B–C, B–GCST0001, C–GCST0001
        "ids": ["B", "C", "GCST0001"], "rg": [7.6455e-06, -0.048, None], "se": [0.096, 0.0826, None], "z": [7.9634e-05, -0.5809, None],
        "p": [0.0, 0.5613, None], "gcov": [1.3498e-06, -0.0037, -0.0023], "gcov_se": [0.0169, 0.0065, 0.0007], "notes": {"2": "rg out of bounds"},
        # the heatmaps' order: B and GCST0001 are the closest (|rg| 0.048), C joins them; the panel is C and GCST0001, without an estimate
        "heat": {"all": {"order": [0, 2, 1], "merges": [[0, 1, 0.952], [3, 2, 1.0]]}, "panel": {"order": [1, 2], "merges": [[0, 1, 1.0]]}}}
    assert out["rg_EAS.json"] == {
        "ids": ["D_EAS", "GCST0002", "GCST0003"], "rg": [0.5, 1.1027, None], "se": [0.05, 0.0168, None], "z": [10.0, 65.6, None],
        "p": [5e-06, 0.0, None], "gcov": [0.01, 0.02, None], "gcov_se": [0.001, 0.002, None], "notes": {},
        # rg 1.1027 counts as 1: distance 0; then GCST0002 at the mean of 1 - 0.5 and 1 - 0 (no estimate)
        "heat": {"all": {"order": [0, 2, 1], "merges": [[0, 1, 0.0], [3, 2, 0.75]]}, "panel": {"order": [1], "merges": []}}}
    assert "NaN" not in (tmp_path / "out" / "rg_EUR.json").read_text()  # not JSON, and a browser would refuse the file


def test_clustering_is_average_linkage():
    pairs = np.array([[0, .1, .8, .9], [.1, 0, .7, .6], [.8, .7, 0, .2], [.9, .6, .2, 0]])  # two close pairs
    assert tree(pairs) == {"order": [0, 1, 2, 3], "merges": [[0, 1, 0.1], [2, 3, 0.2], [4, 5, 0.75]]}  # then their mean distance
    chain = np.array([[0, .9, .2], [.9, 0, .5], [.2, .5, 0]])
    assert tree(chain, ["x", "y", "z"]) == {"order": ["x", "z", "y"], "merges": [[0, 1, 0.2], [3, 2, 0.7]]}  # 0.7: neither 0.5 (nearest) nor 0.9 (farthest)
    assert tree(np.zeros((1, 1)), ["only"]) == {"order": ["only"], "merges": []} and tree(np.zeros((0, 0))) == {"order": [], "merges": []}


def test_paths_are_only_in_the_encrypted_file(mini, tmp_path):
    out = tmp_path / "out"
    blob = built(mini, out)[BLOB]
    for path in out.iterdir():
        assert "/lab/" not in path.read_text() and str(mini) not in path.read_text()
    private = json.loads(crypt.decrypt(blob, "correct horse"))
    assert private["release"] == str(mini) and "summary_statistics_2026.txt" in private["tables"]
    assert private["folders"] == ["/lab/old/sumstats", "/lab/new/sumstats/EUR", "/lab/new/sumstats/EAS", "/lab/new/sumstats_rg/EUR", "/lab/new/sumstats_rg/EAS"]
    assert private["files"]["A_unsigned"] == [[0, "A_unsigned.sumstats"], None]  # unsigned: no rg file
    assert private["files"]["B"] == [[0, "B.sumstats"], [0, "B.sumstats"]]
    assert private["files"]["GCST0003"] == [[2, "GCST0003.sumstats"], [4, "renamed.txt"]]
    assert len(private["files"]) == 7
    with pytest.raises(InvalidTag):
        crypt.decrypt(blob, "correct horse ")


def test_a_rebuild_changes_nothing_and_a_new_passphrase_only_the_blob(mini, tmp_path):
    out = tmp_path / "out"
    read = lambda: {path.name: path.read_bytes() for path in out.iterdir()}
    build(mini, out, "one", FAST)
    first = read()
    build(mini, out, "one", FAST)
    assert read() == first  # the blob too: its salt is kept, so readers stay logged in
    build(mini, out, "two", FAST)
    assert {name for name, content in read().items() if content != first[name]} == {BLOB}
    assert crypt.decrypt(json.loads(read()[BLOB]), "two") == crypt.decrypt(json.loads(first[BLOB]), "one")


@pytest.mark.parametrize("rows, problem", [
    ([["B", "C", *NAN], ["B", "GCST0001", *NAN]], "a pair is missing"),
    ([["B", "C", *NAN], ["C", "B", *NAN], ["B", "GCST0001", *NAN]], "a pair is there twice"),
    ([["B", "C", *NAN], ["B", "GCST0001", *NAN], ["C", "A_unsigned", *NAN]], "an unsigned dataset"),
    ([["B", "C", *NAN], ["B", "GCST0001", *NAN], ["C", "D_EAS", *NAN]], "a dataset of the other ancestry"),
    ([["B", "C", *NAN], ["B", "GCST0001", *NAN], ["C", "C", *NAN]], "a dataset with itself"),
])
def test_a_pair_table_that_is_not_every_pair_once_stops_the_build(mini, tmp_path, rows, problem):
    write(mini / "ldsc_genetic_correlations_trait_pairs_EUR_2026.txt", RG, rows)
    with pytest.raises(ValueError, match="pair table"):
        build(mini, tmp_path / "out", "x", FAST)


@pytest.mark.parametrize("name, old, new, problem", [
    ("sldsc_heritabilities_EAS_2026.txt", "GCST0003", "GCST0009", "heritability"),  # an h² for a dataset that is not in the inventory
    ("summary_statistics_2026.txt", "\tD_EAS\t", "\tC\t", "occurs twice"),
    ("non_redundent_traits_EAS_2026.txt", "GCST0002", "B", "panel"),  # an EUR dataset in the EAS panel
    ("_audit/sumstats_file_index.tsv", "GCST0002\t2026", "GCST0009\t2026", "file index"),
    ("_audit/new_study_provenance.tsv", "https://www.ebi.ac.uk/gwas/studies/GCST0003", "javascript:alert(1)", "link"),
])
def test_a_release_that_disagrees_with_itself_stops_the_build(mini, tmp_path, name, old, new, problem):
    text = (mini / name).read_text()
    assert old in text
    (mini / name).write_text(text.replace(old, new))
    with pytest.raises(ValueError, match=problem):
        build(mini, tmp_path / "out", "x", FAST)


def test_encryption_round_trip():
    blob = crypt.encrypt(b"where the files are", "a passphrase", FAST)
    assert crypt.decrypt(blob, "a passphrase") == b"where the files are"
    assert b"where" not in json.dumps(blob).encode() and "passphrase" not in json.dumps(blob)
    for wrong in ("A passphrase", "a passphrase ", ""):
        with pytest.raises(InvalidTag):
            crypt.decrypt(blob, wrong)
    with pytest.raises(InvalidTag):  # a changed blob is refused, not decrypted to something else
        crypt.decrypt(blob | {"data": blob["data"][:-4] + ("AAA=" if blob["data"][-4:] != "AAA=" else "BBB=")}, "a passphrase")
    again = crypt.encrypt(b"where the files are", "a passphrase", FAST)
    assert again["salt"] != blob["salt"] and again["iv"] != blob["iv"] and again["data"] != blob["data"]  # a new key each time
