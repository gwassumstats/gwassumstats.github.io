"""From a source to the release's format: the rules of the formatting, on made-up files, and where a source is found."""
import gzip

import pandas as pd
import pytest

from conftest import panel, write
from sources import fetch
from sources.format import formatted, named, reference

SSF = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value", "rs_id", "n"]  # as the GWAS Catalog names the columns


@pytest.fixture
def hapmap3(tmp_path):
    """The made-up reference of conftest.panel, as the formatting reads it."""
    return reference(*panel(tmp_path))


def test_a_file_is_formatted_by_the_release_s_rules(tmp_path, hapmap3):
    write(tmp_path / "study.tsv", SSF, [
        ["1", "100", "G", "A", "0.2", "0.1", "0.05", "rs1", "1000"],          # its alleles the other way round: kept as the file has them
        ["CHR01", "200.0", "g", "a", "-0.15", "0.07", "0.03", "rs2", "900.5"],  # on the other strand, in small letters
        ["1", "300", "A", "T", "0.3", "0.1", "0.01", "rs3", "1000"],           # not the alleles of rs3
        ["1", "400", "T", "G", "0.1", "0.1", "0.3", "NA", "1000"],             # no rsID: the SNP at its position
        ["1", "500", "A", "G", "0.1", "0.2", "0.6", "rs90000", "1000"],          # an rsID that is not of HapMap3: the SNP at its position
        ["1", "600", "A", "G", "0.1", "0.2", "0.6", "rs0", "1000"],            # not a HapMap3 SNP, nor at the position of one
        ["1", "100", "A", "G", "9", "1", "0.5", "rs1", "1000"],                # rs1 again: the first row stands
        ["2", "500", "C", "A", "0.4", "0.1", "0.001", "2:500", "1000"],        # two SNPs at this position: the one whose alleles fit
        ["2", "500", "A", "T", "0.4", "0.1", "0.001", "2:500", "1000"],        # the alleles of neither
        ["2", "800", "C", "T", "0.4", "0", "0.001", "rs8", "1000"],            # no Z: the standard error is 0
        ["2", "800", "C", "T", "0.5", "0.25", "0.001", "rs8", "1000"],         # rs8 again, and sound: this row is kept
        ["X", "100", "A", "G", "0.1", "0.2", "0.6", "rs9", "1000"],
        ["1", "200", "C", "T", "0.1", "0.2", "0.6", "rs2", "NA"]])             # no N; and rs2 is written already
    with (tmp_path / "study.tsv").open("a") as file:
        file.write("1\t300\tA\tC\n")  # a row that ends early: read with its last columns empty, so without an N
    table, done = formatted(tmp_path / "study.tsv", *hapmap3)
    assert table.to_csv(sep="\t", index=False, lineterminator="\n") == ("SNP\tA1\tA2\tN\tZ\nrs1\tG\tA\t1000\t2\nrs2\tG\tA\t900.5\t-2.14285714286\nrs4\tT\tG\t1000\t1\n"
                                                                        "rs5\tA\tG\t1000\t0.5\nrs7\tC\tA\t1000\t4\nrs8\tC\tT\t1000\t2\n")
    assert {key: done[key] for key in ("rows", "not_hapmap3", "other_alleles", "without_n", "without_z", "listed_again", "written")} == {
        "rows": 14, "not_hapmap3": 2, "other_alleles": 2, "without_n": 2, "without_z": 1, "listed_again": 1, "written": 6}
    assert (done["n_from"], done["z_from"], done["columns"]["snp"]) == ("n", "beta / standard_error", "rs_id")


@pytest.mark.parametrize("header, row, n, roles, line", [
    (["SNP", "A1", "A2", "Z", "BETA", "SE"], ["rs1", "A", "G", "3.5", "9", "1"], 5000, None, "rs1\tA\tG\t5000\t3.5"),  # a Z stands before an effect; the study's N for a file without one
    (["variant_id", "effect_allele", "other_allele", "odds_ratio", "standard_error"], ["1:100:A:G", "A", "G", "2", "0.5"], 123456.7, None, "rs1\tA\tG\t123456.7\t1.38629436112"),
    (["MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "TotalSampleSize"], ["chr1:100", "a", "g", "0.3", "0.1", "777"], None, {"a1": "Allele1", "a2": "Allele2", "n": "TotalSampleSize"},
     "rs1\tA\tG\t777\t3"),  # columns named for their roles
    (["rsid", "effect_allele", "other_allele", "z_score", "N", "Neff"], ["rs1", "A", "G", "-1e1", "100", "50"], 9, {"n": "Neff"}, "rs1\tA\tG\t50\t-10"),
])
def test_columns_play_their_roles_by_name_or_as_told(tmp_path, hapmap3, header, row, n, roles, line):
    write(tmp_path / "study.tsv", header, [row])
    table, _ = formatted(tmp_path / "study.tsv", *hapmap3, n=n, roles=roles)
    assert table.to_csv(sep="\t", index=False, header=False, lineterminator="\n") == line + "\n"


@pytest.mark.parametrize("text", ["SNP,A1,A2,Z,N\nrs1,A,G,1.5,10\n", "SNP A1  A2 Z N\nrs1 A G   1.5 10\n"])
def test_a_file_may_be_delimited_by_commas_or_blanks_and_zipped(tmp_path, hapmap3, text):
    with gzip.open(tmp_path / "study.out", "wt") as file:  # zipped, whatever its name
        file.write(text)
    assert formatted(tmp_path / "study.out", *hapmap3)[0].values.tolist() == [["rs1", "A", "G", "10", "1.5"]]


@pytest.mark.parametrize("header, n, roles, problem", [
    (["SNP", "Z", "N"], None, None, "no alleles"), (["SNP", "A1", "A2", "BETA", "N"], None, None, "no Z, and no effect"),
    (["SNP", "A1", "A2", "Z"], None, None, "no N among"), (["SNP", "A1", "A2", "Z"], 10, {"n": "Weight"}, "n=Weight: no such role, or no such column"),
    (["SNP", "A1", "A2", "Z"], 10, {"effect": "Z"}, "effect=Z: no such role")])
def test_a_file_that_cannot_be_formatted_says_why(tmp_path, hapmap3, header, n, roles, problem):
    write(tmp_path / "study.tsv", header, [["rs1"] * len(header)])
    with pytest.raises(ValueError, match=problem):
        formatted(tmp_path / "study.tsv", *hapmap3, n=n, roles=roles)


def test_a_role_goes_to_the_first_name_the_release_prefers():
    assert named(["ID", "rsid", "CHR", "BP", "A1", "A2", "OR", "SE", "N", "effective_sample_size"]) == {"snp": "rsid", "chr": "CHR", "pos": "BP", "a1": "A1", "a2": "A2", "or": "OR", "se": "SE", "n": "N"}


def test_a_study_s_folder_in_the_catalog_and_the_file_in_it():
    assert fetch.folder("GCST006414") == fetch.CATALOG + "GCST006001-GCST007000/GCST006414/"
    assert fetch.folder("GCST90292609") == fetch.CATALOG + "GCST90292001-GCST90293000/GCST90292609/" and fetch.folder("GCST90271000").endswith("GCST90270001-GCST90271000/GCST90271000/")
    names = ["GCST1_buildGRCh37.tsv.gz", "GCST1_buildGRCh37.tsv.gz-meta.yaml", "README.txt", "md5sum.txt"]
    assert fetch.chosen(names) == "GCST1_buildGRCh37.tsv.gz" and fetch.chosen(names + ["other.tbl.gz"], "other.tbl.gz") == "other.tbl.gz"
    for listing, file, problem in ((names + ["other.tbl.gz"], "", "2 files of statistics"), (names, "absent.tsv", "no absent.tsv in the folder"), (["README.txt"], "", "none")):
        with pytest.raises(ValueError, match=problem):
            fetch.chosen(listing, file)
    assert fetch.located("https://example.org/study.txt.gz") == ("https://example.org/study.txt.gz", "")  # a URL is taken as it is


def test_a_study_s_file_is_found_in_its_folder_with_the_checksum_stated(monkeypatch):
    pages = {fetch.folder("GCST1001"): '<a href="../">Parent</a><a href="?C=N;O=D">Name</a><a href="harmonised/">harmonised/</a><a href="GCST1001.tsv.gz">x</a><a href="md5sum.txt">m</a>',
             fetch.folder("GCST1001") + "md5sum.txt": "0123abcd  GCST1001.tsv.gz\nffff  GCST1001.tsv.gz-meta.yaml\n",
             fetch.folder("GCST1001") + "harmonised/": '<a href="a.f.tsv.gz">a</a><a href="a.h.tsv.gz">b</a><a href="a.h.tsv.gz-meta.yaml">c</a><a href="readme.txt">d</a>'}
    monkeypatch.setattr(fetch, "read", lambda url: pages[url].encode())
    assert fetch.located("GCST1001") == (fetch.folder("GCST1001") + "GCST1001.tsv.gz", "0123abcd")
    assert fetch.located("GCST1001", "harmonised/a.h.tsv.gz") == (fetch.folder("GCST1001") + "harmonised/a.h.tsv.gz", "")  # no checksum stated there
    with pytest.raises(ValueError, match="2 files of statistics"):
        fetch.located("GCST1001", "harmonised/")


def test_a_download_is_kept_and_checked(tmp_path):
    (tmp_path / "source.txt").write_text("SNP A1 A2\n")
    got = fetch.download((tmp_path / "source.txt").as_uri(), tmp_path / "kept" / "source.txt")
    assert (tmp_path / "kept" / "source.txt").read_text() == "SNP A1 A2\n" and got["bytes"] == 10 and len(got["sha256"]) == 64 and len(got["md5"]) == 32
    assert fetch.download("https://example.invalid/never-asked", tmp_path / "kept" / "source.txt", got["md5"])["sha256"] == got["sha256"]  # there already: not fetched again
    with pytest.raises(ValueError, match="where its folder states"):
        fetch.download("https://example.invalid/never-asked", tmp_path / "kept" / "source.txt", "0" * 32)
    assert not (tmp_path / "kept" / "source.txt").exists()  # a file that is not the one stated is not kept
