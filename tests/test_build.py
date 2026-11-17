"""The build, on a made-up release (conftest.mini): what it writes, what it refuses, and the encryption. And what
finds the papers, drafts the short names and judges the studies, which the build merges in."""
import json

import pytest
from cryptography.exceptions import InvalidTag

import numpy as np
import pandas as pd

from build import crypt
from build.__main__ import BLOB, build as build_with
from build.release import tree
from conftest import MEASURED, NAMES, NAN, PAPERS, RG, data_of, measured, write
from names.__main__ import draft, trait
from papers import cite, pubmed
from papers.__main__ import choose, parse
from qc import checks, labels
from qc.scan import measure

FAST = 1000  # PBKDF2 rounds, for tests that do not need the real cost


def build(release, out, passphrase, iterations=FAST):
    """The build of a made-up release, with the made-up papers and names beside it."""
    return build_with(release, out, passphrase, iterations, data=data_of(release))


def built(release, out, passphrase="correct horse"):
    build(release, out, passphrase)
    return {path.name: json.loads(path.read_text()) for path in out.iterdir()}


def test_studies_carry_the_inventory_h2_panel_and_paper(mini, tmp_path):
    data = built(mini, tmp_path / "out")["studies.json"]
    old, *_, new = data["studies"][:4]
    assert [d["id"] for d in data["studies"]] == ["A_unsigned", "B", "C", "GCST0001", "D_EAS", "GCST0002", "GCST0003"]
    assert old == {"id": "A_unsigned", "name": "Smith_2010_Height", "trait": "Height (old)", "anc": "EUR", "n": 1000, "m": 900, "origin": "Old 2015",
                   "ref": "Smith et al., 2010 Nature", "signed": False, "alleles": True, "cross": False, "release": 2024, "h2": 0.2, "h2_se": 0.02, "h2_z": 10.0,
                   "intercept": 1.01, "qc": "pass", "hm3": 0.9, "chi2": 1.3, "lambda": 1.2, "cite": ["1"]}
    assert new == {"id": "GCST0001", "name": "OrtinVela_2025_TraitNew", "trait": "Trait new", "anc": "EUR", "n": 4000, "m": 930, "origin": "New 2026",
                   "ref": "Ortín et al. 2025 Nat Genet", "signed": True, "alleles": True, "cross": True, "release": 2026, "h2": 0.1461, "h2_se": 0.009,
                   "h2_z": 16.23333, "intercept": 0.7683, "panel": 1, "m_rg": 926, "pmid": "111", "paper": "https://pubmed.ncbi.nlm.nih.gov/111/",
                   "catalog": "https://www.ebi.ac.uk/gwas/studies/GCST0001", "paper_id": "NG_Ortin_2025_x", "qc": "pass", "hm3": 0.93, "chi2": 1.3, "lambda": 1.2,
                   "cite": ["111"], "notes": [["note", "LDSC intercept 0.77: the statistics look corrected for inflation (genomic control), which lowers h²"]]}
    assert {d["id"]: d["panel"] for d in data["studies"] if "panel" in d} == {"GCST0001": 1, "C": 2, "GCST0002": 1}  # the panels' own order
    assert data["studies"][1]["trait"] == "Alzheimer’s disease" and data["studies"][1]["h2"] == -0.0058


def test_a_study_cites_its_paper(mini, tmp_path):
    out = built(mini, tmp_path / "out")
    data, papers = out["studies.json"], out["papers.json"]
    # the release's PMID; the reference's paper; a paper and the resource it names; nothing while the row says "check"
    assert {d["id"]: d.get("cite") for d in data["studies"]} == {"A_unsigned": ["1"], "B": ["2"], "C": None, "GCST0001": ["111"], "D_EAS": ["3", "4"],
                                                                  "GCST0002": ["222"], "GCST0003": ["222"]}
    # who and which year a paper is cited by: its first author; a group by its short name; the inventory's year (2025), not that of the issue (2026)
    assert {pmid: paper["label"] for pmid, paper in data["papers"].items()} == {
        "1": "Smith et al. 2010 Nature", "2": "TEC 2021 Nat Genet", "3": "Kim et al. 2018 Nat Genet", "4": "Kim 2018 Nature", "111": "Ortín Vela et al. 2025 Nat Genet", "222": "Ito 2024 Nat Commun"}
    # the day each was published, as precisely as PubMed dates it: the table of studies sorts its references by it
    assert {pmid: paper["published"] for pmid, paper in data["papers"].items()} == {"1": "2010-10-10", "2": "2021-03", "3": "2018-02-05", "4": "2018", "111": "2025-12-19", "222": "2024-06-01"}
    assert {pmid: {"label": paper["label"], "published": paper["published"]} for pmid, paper in papers.items()} == data["papers"] and list(papers) == ["1", "2", "3", "4", "111", "222"]
    assert {pmid: paper["key"] for pmid, paper in papers.items()} == {"1": "Smith2010", "2": "TEC2021", "3": "Kim2018a", "4": "Kim2018b", "111": "OrtinVela2025", "222": "Ito2024"}
    assert papers["1"] == {"key": "Smith2010", "label": "Smith et al. 2010 Nature", "published": "2010-10-10", "authors": ["Smith J", "Doe A Jr", "Roe B"], "more": 1, "title": "Height of 1,000 people",
                           "journal": "Nature", "year": "2010", "volume": "1", "issue": "2", "pages": "10-9", "doi": "10.1000/one",
                           "bib": "@article{Smith2010,\n  author = {Smith, J. and Doe, Jr, A. and Roe, B. and Poe, C.},\n  title = {{Height of 1,000 people}},\n  journal = {Nature},\n"
                                  "  year = {2010},\n  volume = {1},\n  number = {2},\n  pages = {10--19},\n  doi = {10.1000/one},\n  pmid = {1}\n}\n"}
    assert papers["2"]["bib"] == ("@article{TEC2021,\n  author = {{The Example Consortium}},\n  title = {{Loci \\& genes for Alzheimer's disease}},\n  journal = {Nat Genet},\n"
                                  "  year = {2021},\n  volume = {53},\n  doi = {10.1000/two},\n  pmid = {2}\n}\n")
    assert papers["111"]["year"] == "2026" and "year = {2026}" in papers["111"]["bib"] and papers["111"]["authors"] == ["Ortín Vela S", "Li X"] and papers["111"]["more"] == 0


def test_a_paper_without_its_record_stops_the_build(mini, tmp_path):
    (data_of(mini) / "papers.json").write_text(json.dumps({pmid: paper for pmid, paper in PAPERS.items() if pmid != "222"}))
    with pytest.raises(ValueError, match="lacks the paper 222"):
        build(mini, tmp_path / "out", "x")


def test_a_row_for_one_study_goes_before_the_row_of_its_reference(mini, tmp_path):
    path = data_of(mini) / "references.tsv"
    path.write_text(path.read_text() + "Kim et al. 2018 Nat Genet\tD_EAS\t1\t\tby hand: one reference, two papers\tSmith J 2010 Nature\n")
    assert built(mini, tmp_path / "out")["studies.json"]["studies"][4]["cite"] == ["1"]


def test_every_study_has_a_short_name(mini, tmp_path):
    assert [d["name"] for d in built(mini, tmp_path / "out")["studies.json"]["studies"]] == [name for _, name, _ in NAMES]


@pytest.mark.parametrize("rows, problem", [
    (NAMES[:-1], "lacks the short name of GCST0003"),
    ([*NAMES[:-1], ["GCST0003", "Ito_2024_Trait-F", "by hand"]], "Ito_2024_Trait-F is not FirstAuthor_Year_Trait"),
    ([*NAMES[:-1], ["GCST0003", "Ito_24_TraitF", "by hand"]], "is not FirstAuthor_Year_Trait"),
    ([*NAMES[:-1], ["GCST0003", "ITO_2024_traiteeas", "by hand"]], "ITO_2024_traiteeas names two studies"),  # the same letters as another, in another case
    ([*NAMES, ["GCST0003", "Ito_2024_Other", "by hand"]], "names a study twice"),
])
def test_short_names_that_are_missing_odd_or_shared_stop_the_build(mini, tmp_path, rows, problem):
    write(data_of(mini) / "names.tsv", ["id", "name", "how"], rows)
    with pytest.raises(ValueError, match=problem):
        build(mini, tmp_path / "out", "x")


def test_the_checks_judge_every_study(mini, tmp_path):
    studies = built(mini, tmp_path / "out")["studies.json"]["studies"]
    assert {d["id"]: d["qc"] for d in studies} == {"A_unsigned": "pass", "B": "fail", "C": "warn", "GCST0001": "pass", "D_EAS": "warn", "GCST0002": "warn", "GCST0003": "warn"}
    assert {d["id"]: d["notes"] for d in studies if "notes" in d} == {
        "B": [["fail", "h² -0.0058: not above 0, so nothing for LDSC to work with"]],
        "C": [["warn", "40% of the HapMap3 SNPs that LDSC regresses on"]],
        "GCST0001": [["note", "LDSC intercept 0.77: the statistics look corrected for inflation (genomic control), which lowers h²"]],  # a note: it passes
        # rg 1.10 between Trait D and Trait F: one of the two labels is wrong, and each names the other as {ID}
        "D_EAS": [["warn", "label in doubt: rg 1.10 with {GCST0003} (Trait F), a study of another trait by its label; one of the two labels is wrong"],
                  ["note", "nearly identical to 1 other study (rg of at least 0.95): {GCST0003}, N 7,000"]],
        "GCST0002": [["warn", "5 SNPs listed twice"]],
        "GCST0003": [["warn", "LDSC intercept 1.20 with a mean χ² of 1.25: most of the inflation is not polygenic"],  # the worst first; then in the order of the checks
                     ["warn", "label in doubt: rg 1.10 with {D_EAS} (Trait D), a study of another trait by its label; one of the two labels is wrong"],
                     ["note", "the median N of its SNPs is 4,000; the inventory says 7,000"],
                     ["note", "nearly identical to 1 other study (rg of at least 0.95): {D_EAS}"]]}
    assert [(d["hm3"], d["chi2"], d["lambda"]) for d in studies[2:4]] == [(0.4, 1.3, 1.2), (0.93, 1.3, 1.2)]


@pytest.mark.parametrize("row, finding", [
    (measured("B", 910, n=2000, columns="SNP N Z"), ["fail", "its columns are SNP N Z, not SNP A1 A2 N Z"]),  # the inventory says it has alleles
    (measured("B", 910, n=2000, z_positive=0.97), ["fail", "signed in the inventory, but 97% of its Z are positive: it holds |Z|"]),
    (measured("B", 910, n=2000, z_positive=0.7), ["warn", "70% of its Z are positive"]),
    (measured("B", 910, n=2000, z_missing=100), ["fail", "100 SNPs with a Z that is not a number"]),  # more than 1% of its rows
    (measured("B", 910, n=2000, z_missing=3), ["warn", "3 SNPs with a Z that is not a number"]),
    (measured("B", 910, n=2000, bad_alleles=2), ["fail", "2 SNPs with alleles that are not two of A, C, G, T"]),
    (measured("B", 910, n=2000, not_rs=7), ["warn", "7 SNPs without an rsID"]),
    (measured("B", 911, n=2000), ["warn", "911 SNPs in the file; the inventory says 910"]),
    (measured("B", 910, n=2000, as_list=0, list_reversed=909), ["warn", "its alleles are those of the HapMap3 list, the other way round, at 100% of its SNPs: the mark of a formatting that set the alleles "
                                                                       "without regard to the effect allele. Where such a file was compared with its source, its Z went with its alleles at half of the SNPs only"]),
    (measured("B", 910, n=2000, in_hm3=299), ["fail", "only 299 of the 1,000 HapMap3 SNPs that LDSC regresses on: too few for LDSC"]),
    (measured("B", 910, n=2000, chi2_mean=1.01), ["warn", "mean χ² 1.01: no signal above chance"]),
    (measured("B", 910, n=2000, chi2_mean=0.8), ["warn", "mean χ² 0.80: no signal above chance; the statistics are deflated"]),
    (measured("B", 910, "rg", n=2000, ambiguous=4), ["warn", "its copy for rg: 4 strand-ambiguous SNPs"]),
    (measured("A_unsigned", 900, z_positive=0.5), ["warn", "unsigned in the inventory, but 50% of its Z are positive: they look signed"]),
    (None, ["warn", "its file was not measured"]),
    (measured("B", 920, n=3000, in_hm3=400), ["warn", "by every number measured, its file is that of {C}: one file under two IDs"]),  # C's numbers
])
def test_what_the_checks_find_in_a_file(mini, tmp_path, row, finding):
    study = row[0] if row else "B"
    rows = [line.split("\t") for line in (data_of(mini) / "metrics.tsv").read_text().splitlines()[1:]]
    write(data_of(mini) / "metrics.tsv", MEASURED, [r for r in rows if (r[0], r[1]) != (study, row[1] if row else "h2")] + ([row] if row else []))
    notes = {d["id"]: d["notes"] for d in built(mini, tmp_path / "out")["studies.json"]["studies"] if "notes" in d}
    assert finding in notes[study]


def test_the_catalog_s_record_is_compared_where_there_is_one(mini, tmp_path):
    record = lambda study, n, ancestry, author, published: [study, "1", author, published, "Nat Genet", n, ancestry, "", "A trait", "", "", "A title"]
    write(data_of(mini) / "catalog.tsv", ["id", "pmid", "author", "published", "journal", "n", "ancestry", "sex", "trait", "sample", "cohort", "title"], [
        record("GCST0001", 4000, "European", "Ortin S", "2025-03-01"),  # as the inventory has it, but for an accent: nothing to say
        record("GCST0002", 6500, "European", "Ito K", "2023-12-30")])
    notes = {d["id"]: d.get("notes", []) for d in built(mini, tmp_path / "out")["studies.json"]["studies"]}
    assert [note for note in notes["GCST0001"] if "Catalog" in note[1]] == []
    assert [note for note in notes["GCST0002"] if "Catalog" in note[1]] == [["fail", "ancestry EAS; the GWAS Catalog says European"], ["warn", "N 6,000; the GWAS Catalog says 6,500"],
                                                                         ["warn", "reference Ito et al. 2024 Nat Commun; the GWAS Catalog says Ito K 2023"]]


SAMPLES = ["pmid", "accession", "trait", "n", "ancestries", "sample"]
EUROPEAN = lambda accession, n, sample="", pmid="2": [pmid, accession, "A trait", n, f"European {n}", sample or f"{n:,} European ancestry individuals"]
BOTH = lambda accession, sample: ["2", accession, "A trait", 3000, "European 2050; East Asian 950", sample]  # a sample of two ancestries


@pytest.mark.parametrize("n, records, finding", [
    # an N that is not whole numbers, and whose largest is just below half the effective N of 1,500 cases and 3,000 controls, 2,000
    ((1500.5, 1900.25, 1900.25), [EUROPEAN("G20", 4500, "1,500 European ancestry cases, 3,000 European ancestry controls")],
     ["warn", "its N looks like half an effective N, as PGC files carry (NEFFDIV2): 1,900.25 at most, where the 1,500 cases and 3,000 controls of G20 give 2,000, and 4,500 in all. h² is 2.37 times what their count gives"]),
    # an effective N that is near the count of cases and controls: a note; of two samples that fit alike, the one of the study's ancestry alone is named
    ((1800.5, 1960.5, 1960.5), [BOTH("G20", "up to 900 European ancestry cases, 1,150 European ancestry controls, 400 East Asian ancestry cases, 550 East Asian ancestry controls"),
                                EUROPEAN("G21", 2050, "500 European ancestry male cases, 400 European ancestry female cases, 1,150 European ancestry controls")],
     ["note", "its N looks like an effective N: 1,960.50 at most, where the 900 cases and 1,150 controls of G21 give 2,020, and 2,050 in all. h² is 1.05 times what their count gives"]),
    ((1400.7, 2050.0, 2050.0), [EUROPEAN("G20", 2050, "900 European ancestry cases, 1,150 European ancestry controls")], None),  # the largest N is the count of a sample: nothing to say
    ((1900.5, 2060.5, 2060.5), [EUROPEAN("G20", 2100)], ["note", "its N is not whole numbers (2,060.50 at most): the weights of a meta-analysis, or an effective N"]),
    ((1400.7, 2000.0, 2060.0), [], ["note", "its N is not whole numbers (1,400.70 at least): the weights of a meta-analysis, or an effective N"]),  # the Catalog has nothing of its paper
    ((2000, 2000, 2000), [["2", "G20", "A trait", 2000, "European 1500; East Asian 500", "1,500 European ancestry individuals, 500 East Asian ancestry individuals"], EUROPEAN("G21", 1500)],
     ["warn", "its N, 2,000 at most, is that of a sample of several ancestries: G20 has 2,000 individuals, at most 1,500 of them European, and no sample recorded for its paper has more of these. "
              "The file looks like the analysis of all ancestries together"]),
    ((2000, 2000, 2000), [["2", "G20", "A trait", 2000, "European 1500; East Asian 500", ""], EUROPEAN("G21", 1950)], None),  # a sample of the paper has about as many Europeans
    ((2000, 2000, 2000), [["2", "G20", "A trait", 2000, "European, East Asian 2000", ""]], None),  # the Catalog does not say how many are European: they all may be
    ((2000, 2000, 2000), [EUROPEAN("G20", 700)], ["note", "the median N of its SNPs, 2,000, is above every sample the GWAS Catalog records for its paper (700): it may count samples that the statistics do not hold, "
                                                           "such as a replication stage"]),
    ((2000, 2000, 2000), [EUROPEAN("G20", 9000), EUROPEAN("G21", 12000)], ["note", "the median N of its SNPs, 2,000, is far from every sample the GWAS Catalog records for its paper (9,000 to 12,000): "
                                                                                  "a part of the study, or an effective N"]),
    ((2000, 2000, 2000), [EUROPEAN("G20", 9000), EUROPEAN("G21", 2500)], None),  # one sample of the paper is near enough
    ((2000, 2000, 2000), [EUROPEAN("G20", 9000, pmid="7")], None),  # of another paper
])
def test_the_n_of_a_file_is_compared_with_the_samples_the_catalog_records(mini, tmp_path, n, records, finding):
    rows = [line.split("\t") for line in (data_of(mini) / "metrics.tsv").read_text().splitlines()[1:]]
    write(data_of(mini) / "metrics.tsv", MEASURED, [r for r in rows if r[0] != "B"] + [measured("B", 910, n_min=n[0], n_median=n[1], n_max=n[2])])  # B cites the paper 2; the inventory says N 2,000
    write(data_of(mini) / "samples.tsv", SAMPLES, records)
    notes = {d["id"]: d.get("notes", []) for d in built(mini, tmp_path / "out")["studies.json"]["studies"]}
    assert [note for note in notes["B"] if note[0] != "fail"] == ([finding] if finding else [])  # its fail is of its h²
    assert notes["C"] == [["warn", "40% of the HapMap3 SNPs that LDSC regresses on"]]  # C cites no paper: nothing to compare


def test_a_study_that_is_a_record_of_the_catalog_is_compared_with_that_record(mini, tmp_path):
    write(data_of(mini) / "samples.tsv", SAMPLES, [["222", "GCST0002", "Trait E", 6000, "East Asian 6000", ""], ["222", "GCST0003", "Trait F", 100, "East Asian 100", ""]])
    notes = {d["id"]: [note[1] for note in d.get("notes", [])] for d in built(mini, tmp_path / "out")["studies.json"]["studies"]}
    assert "the median N of its SNPs, 4,000, is above every sample the GWAS Catalog records for it (100): it may count samples that the statistics do not hold, such as a replication stage" in notes["GCST0003"]
    assert not [note for note in notes["GCST0002"] if "Catalog" in note]  # 6,000, as its file has


def test_a_sample_is_counted_for_a_study_s_ancestry():
    row = {"accession": "G1", "n": "3000", "ancestries": "European 1800; NR 250; East Asian 950", "sample": "up to 900 European ancestry cases, 1,150  European ancestry controls, 400 Han Chinese ancestry cases, 550 East Asian ancestry controls"}
    assert checks.recorded(row, "European") == {"accession": "G1", "n": 3000, "own": 2050, "cases": 900, "controls": 1150}  # those of no reported ancestry may be European
    assert checks.recorded(row, "East Asian") == {"accession": "G1", "n": 3000, "own": 1200, "cases": 0, "controls": 550}  # counted by the words that name the ancestry
    japanese = {"accession": "G2", "n": "36792", "ancestries": "East Asian 36792", "sample": "8,180 Japanese ancestry cases, 28,612 Japanese ancestry controls"}
    assert checks.recorded(japanese, "East Asian") == {"accession": "G2", "n": 36792, "own": 36792, "cases": 8180, "controls": 28612}  # a sample of one ancestry: every case and control
    assert checks.recorded(japanese, "European")["own"] == 0


def test_a_label_that_a_person_found_wrong_warns_and_is_read_as_corrected(mini, tmp_path):
    write(data_of(mini) / "labels.tsv", ["id", "label", "why"], [["GCST0003", "Trait D", "its N is that of Trait D"]])
    studies = {d["id"]: d for d in built(mini, tmp_path / "out")["studies.json"]["studies"]}
    assert ["warn", "the release calls it Trait F; its file holds Trait D: its N is that of Trait D"] in studies["GCST0003"]["notes"]
    assert not [note for d in studies.values() for note in d.get("notes", []) if "label in doubt" in note[1]] and studies["D_EAS"]["qc"] == "pass"  # the two agree now


def test_what_a_person_found_counts_like_a_finding_of_the_checks(mini, tmp_path):
    write(data_of(mini) / "findings.tsv", ["id", "level", "finding"], [["A_unsigned", "fail", "its N is that of {B}"], ["A_unsigned", "note", "looked at by hand"], ["Gone", "fail", "of a study that left"]])
    studies = {d["id"]: d for d in built(mini, tmp_path / "out")["studies.json"]["studies"]}
    assert studies["A_unsigned"]["qc"] == "fail" and studies["A_unsigned"]["notes"] == [["fail", "its N is that of {B}"], ["note", "looked at by hand"]]
    write(data_of(mini) / "findings.tsv", ["id", "level", "finding"], [["A_unsigned", "bad", "of no level"]])
    with pytest.raises(ValueError, match="made by hand is a bad"):
        build(mini, tmp_path / "out", "x")


def test_the_candidates_are_passed_on_as_they_are(mini, tmp_path):
    assert built(mini, tmp_path / "out")["candidates.json"] == json.loads((data_of(mini) / "candidates.json").read_text())
    assert "\n" not in (tmp_path / "out" / "candidates.json").read_text().strip()  # on one line: small


def test_release_counts(mini, tmp_path):
    release = built(mini, tmp_path / "out")["studies.json"]["release"]
    assert release == {
        "year": 2026, "packaged": "2026-09-07", "papers": 2, "ancestries": ["EUR", "EAS"], "qc": {"pass": 2, "warn": 4, "fail": 1},
        "all": {"studies": 7, "before": 4, "added": 3, "signed": 6, "signed_before": 3, "unsigned": 1, "unsigned_added": 0, "no_alleles": 0},
        "EUR": {"studies": 4, "before": 3, "added": 1, "signed": 3, "signed_before": 2, "pairs": 3, "failed": 1,
                "worst": {"id": "C", "failed": 1, "pairs": 2}, "panel": 2, "panel_added": 1, "pairs_before": 1, "panel_before": 2, "panel_kept": 1},
        "EAS": {"studies": 3, "before": 1, "added": 2, "signed": 3, "signed_before": 1, "pairs": 3, "failed": 1,
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
    assert sorted(path.name for path in out.iterdir()) == ["candidates.json", "papers.json", BLOB, "rg_EAS.json", "rg_EUR.json", "studies.json"]
    for path in out.iterdir():
        assert "/lab/" not in path.read_text() and str(mini) not in path.read_text()
    private = json.loads(crypt.decrypt(blob, "correct horse"))
    assert private["release"] == str(mini) and "summary_statistics_2026.txt" in private["tables"]
    assert private["folders"] == ["/lab/old/sumstats", "/lab/new/sumstats/EUR", "/lab/new/sumstats/EAS"]
    assert private["files"]["A_unsigned"] == [0, "A_unsigned.sumstats"] and private["files"]["B"] == [0, "B.sumstats"]
    assert private["files"]["GCST0003"] == [2, "GCST0003.sumstats"]  # its file, not the copy that rg was estimated on
    assert len(private["files"]) == 7 and "sumstats_rg" not in json.dumps(private) and "renamed.txt" not in json.dumps(private)
    with pytest.raises(InvalidTag):
        crypt.decrypt(blob, "correct horse ")


def test_a_rebuild_changes_nothing_and_a_new_passphrase_only_the_blob(mini, tmp_path):
    out = tmp_path / "out"
    read = lambda: {path.name: path.read_bytes() for path in out.iterdir()}
    build(mini, out, "one")
    first = read()
    build(mini, out, "one")
    assert read() == first  # the blob too: its salt is kept, so readers stay logged in
    build(mini, out, "two")
    assert {name for name, content in read().items() if content != first[name]} == {BLOB}
    assert crypt.decrypt(json.loads(read()[BLOB]), "two") == crypt.decrypt(json.loads(first[BLOB]), "one")


@pytest.mark.parametrize("rows, problem", [
    ([["B", "C", *NAN], ["B", "GCST0001", *NAN]], "a pair is missing"),
    ([["B", "C", *NAN], ["C", "B", *NAN], ["B", "GCST0001", *NAN]], "a pair is there twice"),
    ([["B", "C", *NAN], ["B", "GCST0001", *NAN], ["C", "A_unsigned", *NAN]], "an unsigned study"),
    ([["B", "C", *NAN], ["B", "GCST0001", *NAN], ["C", "D_EAS", *NAN]], "a study of the other ancestry"),
    ([["B", "C", *NAN], ["B", "GCST0001", *NAN], ["C", "C", *NAN]], "a study with itself"),
])
def test_a_pair_table_that_is_not_every_pair_once_stops_the_build(mini, tmp_path, rows, problem):
    write(mini / "ldsc_genetic_correlations_trait_pairs_EUR_2026.txt", RG, rows)
    with pytest.raises(ValueError, match="pair table"):
        build(mini, tmp_path / "out", "x")


@pytest.mark.parametrize("name, old, new, problem", [
    ("sldsc_heritabilities_EAS_2026.txt", "GCST0003", "GCST0009", "heritability"),  # an h² for a study that is not in the inventory
    ("summary_statistics_2026.txt", "\tD_EAS\t", "\tC\t", "occurs twice"),
    ("non_redundent_traits_EAS_2026.txt", "GCST0002", "B", "panel"),  # an EUR study in the EAS panel
    ("_audit/sumstats_file_index.tsv", "GCST0002\t2026", "GCST0009\t2026", "file index"),
    ("_audit/new_study_provenance.tsv", "https://www.ebi.ac.uk/gwas/studies/GCST0003", "javascript:alert(1)", "link"),
])
def test_a_release_that_disagrees_with_itself_stops_the_build(mini, tmp_path, name, old, new, problem):
    text = (mini / name).read_text()
    assert old in text
    (mini / name).write_text(text.replace(old, new))
    with pytest.raises(ValueError, match=problem):
        build(mini, tmp_path / "out", "x")


# ---------------------------------------------------------------- papers: how one is found and cited


def test_a_reference_of_the_inventory_names_author_year_and_journal():
    assert parse("Speliotes et al., 2010 Nat Genet") == ("Speliotes", "2010", "Nat Genet")
    assert parse("NG_Kanai_2018_BBJ") == ("Kanai", "2018", "Nat Genet")
    assert parse("Karlsson Linner et al. 2019 Nat Genet") == ("Karlsson Linner", "2019", "Nat Genet")
    assert parse("Manning et. al., 2012 Nat Genet") == ("Manning", "2012", "Nat Genet")
    assert parse("Kessler et al. Nature") == ("Kessler", "", "Nature")
    assert parse("Sanchez-Roige et al. 2018 American Journal of Psychiatry") == ("Sanchez-Roige", "2018", "Am J Psychiatry")
    assert parse("PAN UKBB") == ("", "", "")  # no first author: left for a hand


def test_the_day_a_paper_was_published_is_as_precise_as_pubmed_gives_it():
    assert [pubmed.date(text) for text in ("2010 Oct 10", "2010 Nov", "2019 Jan-Feb", "2018 Dec 15-31", "2010", "2019 Spring", "2014 Feb 6", "")] \
        == ["2010-10-10", "2010-11", "2019-01", "2018-12-15", "2010", "2019", "2014-02-06", ""]
    # the earliest of the day online and the day in print; of two that may be the same day, the more precise
    assert pubmed.earliest("2014-02-20", "2013-12-25") == "2013-12-25" and pubmed.earliest("2026-01", "2025-12-19") == "2025-12-19"
    assert pubmed.earliest("2019-01", "2019-01-14") == "2019-01-14" and pubmed.earliest("2019-01", "2019-02-02") == "2019-01"
    assert pubmed.earliest("2018-07", "") == "2018-07" and pubmed.earliest("2019", "2019-06-01") == "2019-06-01" and pubmed.earliest("", "") == ""


def test_of_several_papers_the_one_whose_title_names_the_trait():
    found = {"1": {"title": "Genome-wide analysis of insomnia in 1,331,010 individuals"}, "2": {"title": "Genome-wide meta-analysis identifies new loci for Alzheimer's disease risk"}}
    assert choose(["Insomnia"], found) == ("1", "of 2 papers, the one whose title names the trait")
    assert choose(["Alzheimers Disease"], found)[0] == "2"
    assert "check" in choose(["Height"], found)[1]  # neither: a hand decides
    assert choose(["Insomnia"], {"1": found["1"]}) == ("1", "the only paper: its title names the trait")
    assert choose(["Height"], {"3": {"title": "A report"}}) == ("3", "the only paper; check: its title does not name the trait")


def test_who_and_which_year_a_paper_is_cited_by():
    groups = {"The Example Consortium": "TEC"}
    paper = lambda *authors, year="2019": {"authors": list(authors), "year": year, "journal": "Nat Genet"}
    assert cite.who(paper("de Lange KM", "Moutsianas L"), [], groups) == ("de Lange", True)
    assert cite.who(paper("Karlsson Linnér R", "Biroli P"), [], groups) == ("Karlsson Linnér", True) and cite.who(paper("Ito K"), [], groups) == ("Ito", False)
    assert cite.who(paper("Taylor HA Jr", "Ito K"), [], groups) == ("Taylor", True)
    assert cite.who(paper("{The Example Consortium}"), ["Ruderfer et al. 2018 Cell"], groups) == ("TEC", False)  # a group that is listed
    assert cite.who(paper("{Another Working Group}"), ["Ruderfer et al. 2018 Cell"], groups) == ("Ruderfer", True)  # else whom the inventory names
    assert cite.who(paper("{Another Working Group}"), ["NC_deRojas_2021_ADGWAS"], groups) == ("deRojas", True)
    assert cite.who(paper("{Another Working Group}"), ["PAN UKBB"], groups) == ("Another Working Group", False)
    # the year the inventory gives the paper; PubMed's if it gives none, or a preprint's, or several with PubMed's among them
    assert cite.year(paper("Kany S", year="2026"), ["Kany et al. 2025 Nat Genet"] * 8) == "2025"
    assert cite.year(paper("Okada Y", year="2014"), ["Okada et al., 2014 Nature", "Nature_Okada_2013_JapaneseRAGWAS"]) == "2014"
    assert cite.year(paper("{A League}", year="2023"), ["Berkovic et al. 2022 medRxiv"]) == "2023" and cite.year(paper("Kessler MD", year="2022"), ["Kessler et al. Nature"]) == "2022"
    assert cite.year(paper("Bycroft C", year="2018"), []) == "2018" and cite.year(paper("Li X", year="2020"), ["Li 2018", "Li 2019", "Li 2019"]) == "2019"
    assert cite.label(paper("Kany S", "Li X", year="2026"), ["Kany et al. 2025 Nat Genet"], groups) == "Kany et al. 2025 Nat Genet"
    assert cite.label(paper("{The Example Consortium}"), [], groups) == "TEC 2019 Nat Genet"


def test_a_bibtex_key_is_who_and_year_and_a_letter_among_several():
    assert cite.keys({"30617256": ("Jansen", "2019"), "30804565": ("Jansen", "2019"), "5": ("de Lange", "2017"), "6": ("Ortín Vela", "2024"), "7": ("PGC", "2014")}) \
        == {"30617256": "Jansen2019a", "30804565": "Jansen2019b", "5": "deLange2017", "6": "OrtinVela2024", "7": "PGC2014"}
    assert cite.keys({"100": ("Li", "2020"), "20": ("Li", "2020"), "3": ("Li", "2020")}) == {"3": "Li2020a", "20": "Li2020b", "100": "Li2020c"}  # by PMID, as numbers


def test_pages_and_authors_as_bibtex_writes_them():
    assert [cite.pages(p) for p in ("937-48", "441-7", "1214-1231.e11", "2380-2391", "e1003356", "1234", "")] == ["937--948", "441--447", "1214--1231.e11", "2380--2391", "e1003356", "1234", ""]
    assert cite.person("Speliotes EK") == ("Speliotes", "EK", "") and cite.person("Taylor HA Jr") == ("Taylor", "HA", "Jr") and cite.person("van der Harst P") == ("van der Harst", "P", "")
    assert cite.person("Dan Rujescu") == ("Dan Rujescu", "", "")  # not surname and initials: left whole
    entry = cite.bibtex("Key2020", "42", {"authors": ["Mägi R", "Longstreth WT Jr", "Dan Rujescu", "{A Group}"], "title": "100% of {loci}_x", "journal": "PLoS Genet", "year": "2020",
                                          "volume": "", "issue": "", "pages": "", "doi": "10.1371/journal.pgen.1"})
    assert entry == ("@article{Key2020,\n  author = {Mägi, R. and Longstreth, Jr, W. T. and Dan Rujescu and {A Group}},\n  title = {{100\\% of loci\\_x}},\n  journal = {PLoS Genet},\n"
                     "  year = {2020},\n  doi = {10.1371/journal.pgen.1},\n  pmid = {42}\n}\n")


# ---------------------------------------------------------------- short names: how one is drafted


PHRASES = {"body mass index": "BMI", "bmi": "BMI", "ldl direct adjusted by medication": "LDLAdjMed", "ldldirect": "LDL", "forced vital capacity": "FVC", "(mtag)": "MTAG", "(uk biobank)": "",
           "(female)": "Female", "type 2 diabetes": "T2D", "systolic blood pressure": "SBP", "(bmi adjusted)": "AdjBMI", "x sex interaction": "xSex", "alzheimer's disease": "Alzheimers",
           "average axial diffusivity (ad)": "AxialDiffusivity"}


@pytest.mark.parametrize("name, short", [
    ("Body Mass Index", "BMI"), ("BMI (UK Biobank)", "BMI"), ("LDLdirect", "LDL"), ("LDL Direct Adjusted By Medication", "LDLAdjMed"),  # a phrase of the list, the longest first
    ("Fasting Glucose", "FastingGlucose"), ("Age of Drinking Initiation", "AgeDrinkingInitiation"), ("Drinks per Week", "DrinksPerWeek"),  # words kept, small ones dropped
    ("Hemoglobin (HGB)", "HGB"), ("Lung Volume (LungVol)", "LungVol"), ("Glycated haemoglobin (HbA1c)", "HbA1c"),  # a name that ends in its code
    ("Aortic stenosis (MTAG)", "AorticStenosisMTAG"), ("Average axial diffusivity (AD)", "AxialDiffusivity"),  # but a qualifier is no code, and a listed name goes first
    ("ADHD", "ADHD"), ("IGF1", "IGF1"), ("Hypothyroidism NOS", "HypothyroidismNOS"), ("CystatinC", "CystatinC"),  # capitals inside a word are kept
    ("Breast Cancer (female)", "BreastCancerFemale"), ("Systolic blood pressure x sex interaction", "SBPxSex"), ("Systolic blood pressure (BMI adjusted)", "SBPAdjBMI"),
    ("Alzheimer’s disease", "Alzheimers"), ("Type 2 Diabetes (UK Biobank)", "T2D"), ("Forced vital capacity (FVC)", "FVC"), ("Raynaud's syndrome", "RaynaudsSyndrome"),
])
def test_a_trait_name_as_one_word(name, short):
    assert trait(name, PHRASES) == short


def test_short_names_are_drafted_from_paper_and_trait_and_told_apart():
    study = lambda i, name, anc="EUR", n=1000, ref="Li et al. 2024 Nat Genet": {"id": i, "trait": name, "anc": anc, "n": n, "ref": ref}
    studies = [study("a", "Body Mass Index"), study("b", "Body mass index"), study("c", "BMI", n=3000), study("d", "Body Mass Index", anc="EAS"),
               study("e", "Type 2 Diabetes"), study("f", "Type 2 Diabetes", n=500), study("g", "Type 2 Diabetes", n=700),
               study("h", "Fasting Glucose", anc="EAS"), study("i", "Height", ref="Smith et al., 2010 Nature"), study("j", "Height", ref="PAN UKBB"),
               study("k", "Breast Cancer (female)"), study("l", "Height"), study("m", "Preeclampsia"), study("n", "Prostate cancer")]
    cited = {s["id"]: ("Li", "2024") for s in studies} | {"i": None, "j": None}
    catalog = {"a": {"sex": "F", "sample": ""}, "b": {"sex": "", "sample": "1,000 European ancestry male cases"}, "c": {"sex": "combined", "sample": "500 female cases, 500 male cases"},
               "h": {"sex": "F", "sample": ""}, "k": {"sex": "F", "sample": ""}, "l": {"sex": "M", "sample": ""}, "m": {"sex": "F", "sample": ""}, "n": {"sex": "M", "sample": ""}}
    women, men = "drafted: of women only, by the GWAS Catalog", "drafted: of men only, by the GWAS Catalog"
    assert draft(studies, cited, PHRASES, catalog) == {
        "a": ("Li_2024_BMIFemale", women), "b": ("Li_2024_BMIMale", men), "c": ("Li_2024_BMI", "drafted"),  # of one sex: by the submitters' word, or else the sample's
        "d": ("Li_2024_BMIEAS", "drafted"), "h": ("Li_2024_FastingGlucoseFemaleEAS", women),  # an East Asian study ends in EAS
        "e": ("Li_2024_T2D", "drafted"), "f": ("Li_2024_T2DN0k", "check: told apart by its sample size only"), "g": ("Li_2024_T2DN1k", "check: told apart by its sample size only"),
        "i": ("Smith_2010_Height", "check: no paper found"), "j": ("Unknown_0000_Height", "check: no paper found"),  # without a paper: whom and when the inventory names
        "k": ("Li_2024_BreastCancerFemale", "drafted"), "l": ("Li_2024_HeightMale", men),  # a name that says so already; a single-sex study alone of its trait says so too
        "m": ("Li_2024_Preeclampsia", "drafted"), "n": ("Li_2024_ProstateCancer", "drafted")}  # but not of a trait that one sex only has


# ---------------------------------------------------------------- the checks: a study's label, a file's numbers


LABEL_RULES = {"twin_rg": 0.95, "twin_se": 0.05, "agree_rg": 0.7, "conflict_rg": 0.6, "conflict_h2_z": 4, "ignore": ["adjusted", "by", "medication", "uk", "biobank"],
               "related": [["hemoglobin", "hematocrit"]]}


def label_check(traits, rg, se=0.02, h2_z=10):
    """The findings for studies a, b, c, ... with these traits, of one ancestry, and these rg by pair ("ab": 0.99)."""
    ids = list("abcdefgh"[:len(traits)])
    pairs = [x + y for i, x in enumerate(ids) for y in ids[i + 1:]]
    table = {"ids": ids, "rg": [rg.get(pair) for pair in pairs], "se": [se if pair in rg else None for pair in pairs]}
    return labels.check({i: {"trait": trait, "h2_z": h2_z} for i, trait in zip(ids, traits)}, {"EUR": table}, LABEL_RULES)


def test_a_study_that_is_another_trait_s_twin_fails_when_that_trait_is_confirmed():
    found = label_check(["LDL Direct Adjusted By Medication", "FEV1-FVC Ratio", "FEV1/FVC ratio", "LDLdirect"], {"ab": 0.99, "ac": 1.0, "bc": 0.98, "ad": 0.08, "bd": 0.01, "cd": 0.02})
    assert found == {"a": [("fail", "label in doubt: rg 0.99 with {b} (FEV1-FVC Ratio), a study of another trait by its label; other studies confirm that label, none this study's"),
                           ("fail", "label in doubt: rg 1.00 with {c} (FEV1/FVC ratio), a study of another trait by its label; other studies confirm that label, none this study's"),
                           ("warn", "label in doubt: rg 0.08 only with {d} (LDLdirect), a study of the same trait by its label")]}  # d is not blamed: a failed one explains it


def test_a_study_known_to_be_at_fault_explains_a_disagreement():
    pair = (["Atrial Fibrillation And Flutter", "Atrial Fibrilation"], {"ab": -0.19})
    assert set(label_check(*pair)) == {"a", "b"}
    ids = ["a", "b"]
    table = {"EUR": {"ids": ids, "rg": [-0.19], "se": [0.02]}}
    found = labels.check({i: {"trait": trait, "h2_z": 10} for i, trait in zip(ids, pair[0])}, table, LABEL_RULES, faulty={"b"})
    assert set(found) == {"b"}  # b's Z does not go with its alleles: a is not blamed for disagreeing with it
    twins = labels.check({"a": {"trait": "Operation Tonsillectomy Adenoids", "h2_z": 10}, "b": {"trait": "Varicose Veins", "h2_z": 10}}, {"EUR": {"ids": ids, "rg": [1.03], "se": [0.01]}}, LABEL_RULES, faulty={"a"})
    assert set(twins) == {"a"}  # a holds b's trait: b's label is not in doubt


def test_two_twins_of_different_labels_both_warn_when_nothing_confirms_either():
    found = label_check(["Operation Tonsillectomy Adenoids", "Varicose Veins"], {"ab": 1.03})
    assert found == {"a": [("warn", "label in doubt: rg 1.03 with {b} (Varicose Veins), a study of another trait by its label; one of the two labels is wrong")],
                     "b": [("warn", "label in doubt: rg 1.03 with {a} (Operation Tonsillectomy Adenoids), a study of another trait by its label; one of the two labels is wrong")]}


def test_two_studies_of_one_label_that_disagree_warn_unless_another_confirms_one():
    two = label_check(["Atrial Fibrillation And Flutter", "Atrial Fibrilation"], {"ab": -0.19})  # one label, but for a slip of the pen and a longer name
    assert {study: [level for level, _ in findings] for study, findings in two.items()} == {"a": ["warn"], "b": ["warn"]}
    three = label_check(["Pulse Pressure Adjusted By Medication", "Pulse pressure", "Pulse pressure (UK Biobank)"], {"ab": 0.11, "ac": 0.09, "bc": 0.97})
    assert set(three) == {"a"} and [level for level, _ in three["a"]] == ["warn", "warn"]  # b and c confirm each other
    assert label_check(["Height", "Height"], {"ab": 0.3}, se=0.2) == {}  # too uncertain to say they disagree
    assert label_check(["Height", "Height"], {"ab": 0.3}, h2_z=3) == {}  # too weak a study to say so
    assert label_check(["Height", "Height"], {}) == {}  # no estimate, no finding
    assert labels.check({i: {"trait": "Height", "h2_z": 10} for i in "ab"}, {"EUR": {"ids": ["a", "b"], "rg": [0.3], "se": [0.02]}}, LABEL_RULES | {"distinct": [["b", "a"]]}) == {}  # a pair known to differ


def test_related_labels_may_be_twins():
    assert label_check(["Hemoglobin", "Hematocrit (HCT)"], {"ab": 0.97}) == {}  # a pair that is listed
    assert label_check(["Knee Osteoarthritis", "Osteoarthritis in Knee (KneeOA)"], {"ab": 0.99}) == {}  # a word in common
    assert label_check(["Thyroid", "Thyroid stimulating hormone levels"], {"ab": 0.2}) == {}  # not the same label: a word alone is a label of its own
    assert label_check(["Hair Balding Pattern", "Hair Balding Pattern 2"], {"ab": 0.14}) == {}  # nor these: a number makes another label
    assert label_check(["LDL", "Height"], {"ab": 0.99}, se=0.2) == {}  # too uncertain to call them twins


def test_a_file_is_measured_snp_by_snp(tmp_path):
    path = tmp_path / "X.sumstats"
    path.write_text("SNP\tA1\tA2\tN\tZ\n" + "".join(f"{row}\n" for row in (
        "rs1\tA\tG\t1000\t2.0", "rs2\tA\tT\t1000\t-1.0", "rs2\tC\tG\t900\t1.0", "rs3\tA\tA\t800\tinf", "chr1:5\tG\tT\t1000\t-3.0", "rs6\tN\tC\t\t0.5")))
    alleles = pd.DataFrame({"A1": ["G", "A", "T"], "A2": ["A", "T", "G"]}, index=pd.Index(["rs1", "rs2", "chr1:5"], name="SNP"))
    m = measure(str(path), pd.Index(["rs1", "rs2", "rs8", "rs9"]), alleles)
    assert (m["as_list"], m["list_reversed"]) == (1, 2) and "as_list" not in measure(str(path), pd.Index(["rs1"]))  # rs2 once as the list has it; rs1 and chr1:5 the other way round
    assert (m["columns"], m["rows"], m["duplicates"], m["not_rs"], m["in_hm3"], m["hm3"]) == ("SNP A1 A2 N Z", 6, 1, 1, 2, 4)  # rs2, listed twice, counts once
    assert (m["bad_alleles"], m["ambiguous"], m["n_missing"], m["n_median"], m["z_missing"]) == (2, 2, 1, 1000, 1)  # A/A and N/C; A/T and C/G; an infinite Z is counted and left out
    assert m["z_positive"] == 0.6 and m["chi2_mean"] == pytest.approx((4 + 1 + 1 + 9 + 0.25) / 5) and m["lambda_gc"] == pytest.approx(1 / 0.4549364)
    plain = tmp_path / "Y.sumstats"
    plain.write_text("SNP\tN\tZ\nrs1\t10\t1.5\n")
    assert measure(str(plain), pd.Index(["rs1"]))["columns"] == "SNP N Z" and "bad_alleles" not in measure(str(plain), pd.Index(["rs1"]))


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
