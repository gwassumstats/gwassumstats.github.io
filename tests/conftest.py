"""What the tests share: the real release and passphrase, from the environment and never from a file, and a made-up
miniature release, with what a repository would add to it, that tests the build anywhere."""
import json
import os
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
SITE = ROOT / "site"
DATA = SITE / "data"

INVENTORY = ["Trait Name", "Trait Identifier", "Origin", "Reference", "N", "M", "Population", "Signed SumStat", "Allele Info", "Cross-Trait LDSC"]
H2 = ["trait_identifier", "h2", "h2_se", "h2_z", "intercept"]
RG = ["trait1", "trait2", "genetic_correlation", "genetic_correlation_se", "genetic_correlation_z", "genetic_correlation_p", "genetic_covariance", "genetic_covariance_se"]
NAN = ["nan"] * 6


@pytest.fixture
def release():
    """The release folder the site is built from."""
    folder = os.environ.get("SUMSTATS_RELEASE")
    if not folder or not Path(folder).is_dir():
        pytest.skip("set SUMSTATS_RELEASE to the release folder")
    return Path(folder)


@pytest.fixture
def passphrase():
    return os.environ.get("SUMSTATS_PASSPHRASE") or pytest.skip("set SUMSTATS_PASSPHRASE")


def write(path: Path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join("\t".join(map(str, row)) for row in [header, *rows]) + "\n", encoding="utf-8")


PAPERS = {  # as PubMed would record the papers of the made-up release
    "1": {"authors": ["Smith J", "Doe A Jr", "Roe B", "Poe C"], "title": "Height of 1,000 people", "journal": "Nature", "year": "2010", "volume": "1", "issue": "2", "pages": "10-9", "doi": "10.1000/one"},
    "2": {"authors": ["{The Example Consortium}"], "title": "Loci & genes for Alzheimer's disease", "journal": "Nat Genet", "year": "2021", "volume": "53", "issue": "", "pages": "", "doi": "10.1000/two"},
    "3": {"authors": ["Kim A", "Park B"], "title": "Trait D in East Asians", "journal": "Nat Genet", "year": "2018", "volume": "50", "issue": "1", "pages": "100-110", "doi": "10.1000/three"},
    "4": {"authors": ["Kim C"], "title": "A biobank", "journal": "Nature", "year": "2018", "volume": "9", "issue": "", "pages": "e12", "doi": ""},
    "111": {"authors": ["Ortín Vela S", "Li X"], "title": "A new trait", "journal": "Nat Genet", "year": "2026", "volume": "58", "issue": "1", "pages": "1-9", "doi": "10.1000/new"},
    "222": {"authors": ["Ito K"], "title": "Traits E and F", "journal": "Nat Commun", "year": "2024", "volume": "15", "issue": "", "pages": "1234", "doi": "10.1000/ito"},
}
NAMES = [["A_unsigned", "Smith_2010_Height", "drafted"], ["B", "TEC_2021_Alzheimers", "drafted"], ["C", "Lee_2022_TraitC", "check: no paper found"],
         ["GCST0001", "OrtinVela_2025_TraitNew", "drafted"], ["D_EAS", "Kim_2018_TraitDEAS", "drafted"], ["GCST0002", "Ito_2024_TraitEEAS", "by hand"],
         ["GCST0003", "Ito_2024_TraitFEAS", "drafted"]]


def data_of(release: Path) -> Path:
    """The folder that plays data/ for a made-up release."""
    return release.parent / "data"


@pytest.fixture
def mini(tmp_path):
    """A release folder of seven studies, laid out as a real one: four EUR (one unsigned, one added in 2026) and
    three EAS (two added). The pair tables name the pairs in no particular order. Beside it, data_of(release): the
    papers of the studies and their short names, as data/ holds them for the real release."""
    release = tmp_path / "release_2026_20260907"
    data = data_of(release)
    write(data / "references.tsv", ["reference", "study", "pmid", "doi", "how", "paper"], [
        ["Smith et al., 2010 Nature", "", "1", "10.1000/one", "author, year and journal: the only paper", "Smith J 2010 Nature: Height of 1,000 people"],
        ["NG_Lee_2021_x", "", "2", "10.1000/two", "by hand: a consortium, found by its title", "{The Example Consortium} 2021 Nat Genet: Loci & genes"],
        ["Lee et al. 2022 Cell", "", "9", "", "author, year and journal: of 3 papers, check: none stands out by its title", "Lee Q 2022 Cell: Something else"],
        ["Kim et al. 2018 Nat Genet", "", "3;4", "10.1000/three;", "by hand: the paper, and the resource it names", "Kim A 2018 Nat Genet: Trait D in East Asians"]])
    (data / "papers.json").write_text(json.dumps(PAPERS, ensure_ascii=False), encoding="utf-8")
    write(data / "consortia.tsv", ["group", "short"], [["# a comment line", ""], ["The Example Consortium", "TEC"]])
    write(data / "names.tsv", ["id", "name", "how"], NAMES)
    inventory = [
        ["Height (old)", "A_unsigned", "Old 2015", "Smith et al., 2010 Nature ", 1000, 900, "EUR", "NO", "YES", "NO"],
        ["Alzheimer’s disease", "B", "New 2024", "NG_Lee_2021_x", 2000, 910, "EUR", "YES", "YES", "YES"],
        ["Trait C", "C", "New 2024", "Lee et al. 2022 Cell", 3000, 920, "EUR", "YES", "YES", "YES"],
        ["Trait new", "GCST0001", "New 2026", "Ortín et al. 2025 Nat Genet", 4000, 930, "EUR", "YES", "YES", "YES"],
        ["Trait D", "D_EAS", "New 2021", "Kim et al. 2018 Nat Genet", 5000, 940, "EAS", "YES", "YES", "YES"],
        ["Trait E", "GCST0002", "New 2026", "Ito et al. 2024 Nat Commun", 6000, 950, "EAS", "YES", "YES", "YES"],
        ["Trait F", "GCST0003", "New 2026", "Ito et al. 2024 Nat Commun", 7000, 960, "EAS", "YES", "YES", "YES"],
    ]
    write(release / "summary_statistics_2026.txt", INVENTORY, inventory)
    write(release / "signed_summary_statistics_2026.txt", INVENTORY, [row for row in inventory if row[7] == "YES"])
    write(release / "sldsc_heritabilities_EUR_2026.txt", H2, [["A_unsigned", "0.2", "0.02", "10.0", "1.01"], ["B", "-0.0058", "0.01", "-0.58", "0.99"],
                                                             ["C", "0.393", "0.0674", "5.83086", "1.0297"], ["GCST0001", "0.1461", "0.009", "16.23333", "0.7683"]])
    write(release / "sldsc_heritabilities_EAS_2026.txt", H2, [["D_EAS", "0.1", "0.01", "10.0", "1.0"], ["GCST0002", "0.3", "0.03", "10.0", "1.1"], ["GCST0003", "0.05", "0.02", "2.5", "1.2"]])
    write(release / "non_redundent_traits_EUR_2026.txt", H2, [["GCST0001", "0.1461", "0.009", "16.23333", "0.7683"], ["C", "0.393", "0.0674", "5.83086", "1.0297"]])
    write(release / "non_redundent_traits_EAS_2026.txt", H2, [["GCST0002", "0.3", "0.03", "10.0", "1.1"]])
    write(release / "ldsc_genetic_correlations_trait_pairs_EUR_2026.txt", RG, [
        ["GCST0001", "B", "-0.048", "0.0826", "-0.5809", "0.5613", "-0.0037", "0.0065"],  # the later study first
        ["C", "GCST0001", "nan", "nan", "nan (nan) (rg out of bounds)", "nan (nan) (rg out of bounds)", "-0.0023", "0.0007"],
        ["B", "C", "7.6455e-06", "0.096", "7.9634e-05", "0.", "1.3498e-06", "0.0169"]])
    write(release / "ldsc_genetic_correlations_trait_pairs_EAS_2026.txt", RG, [
        ["D_EAS", "GCST0002", "0.5", "0.05", "10.", "5.e-06", "0.01", "0.001"], ["GCST0003", "GCST0002", *NAN], ["D_EAS", "GCST0003", "1.1027", "0.0168", "65.6", "0.", "0.02", "0.002"]])
    write(release / "_audit" / "sumstats_file_index.tsv", ["Trait Identifier", "Release", "Population", "H2 Sumstats", "RG Sumstats", "M in summary", "M in RG input"], [
        ["A_unsigned", "2024", "EUR", "/lab/old/sumstats/A_unsigned.sumstats", "", 900, ""],
        ["B", "2024", "EUR", "/lab/old/sumstats/B.sumstats", "/lab/old/sumstats/B.sumstats", 910, 910],
        ["C", "2024", "EUR", "/lab/old/sumstats/C.sumstats", "/lab/old/sumstats/C.sumstats", 920, 920],
        ["GCST0001", "2026 addition", "EUR", "/lab/new/sumstats/EUR/GCST0001.sumstats", "/lab/new/sumstats_rg/EUR/GCST0001.sumstats", 930, 926],
        ["D_EAS", "2024", "EAS", "/lab/old/sumstats/D_EAS.sumstats", "/lab/old/sumstats/D_EAS.sumstats", 940, 940],
        ["GCST0002", "2026 addition", "EAS", "/lab/new/sumstats/EAS/GCST0002.sumstats", "/lab/new/sumstats_rg/EAS/GCST0002.sumstats", 950, 950],
        ["GCST0003", "2026 addition", "EAS", "/lab/new/sumstats/EAS/GCST0003.sumstats", "/lab/new/sumstats_rg/EAS/renamed.txt", 960, 958]])
    write(release / "_audit" / "new_study_provenance.tsv", ["Trait Identifier", "Target Paper Identifier", "Reference", "PMID", "Paper Link", "Catalog Study Link", "Origin in source master"], [
        ["GCST0001", "NG_Ortin_2025_x", "Ortín et al. 2025 Nat Genet", "111", "https://pubmed.ncbi.nlm.nih.gov/111/", "https://www.ebi.ac.uk/gwas/studies/GCST0001", "a list"],
        ["GCST0002", "NC_Ito_2024_y", "Ito et al. 2024 Nat Commun", "222", "https://pubmed.ncbi.nlm.nih.gov/222/", "https://www.ebi.ac.uk/gwas/studies/GCST0002", "a list"],
        ["GCST0003", "NC_Ito_2024_y", "Ito et al. 2024 Nat Commun", "222", "https://pubmed.ncbi.nlm.nih.gov/222/", "https://www.ebi.ac.uk/gwas/studies/GCST0003", "a list"]])
    write(release / "_audit" / "inputs" / "non_redundent_traits_EUR_2024.txt", H2, [["B", "0.3", "0.01", "30.0", "1.0"], ["C", "0.393", "0.0674", "5.83086", "1.0297"]])
    (release / "README.md").write_text("A made-up release.\n")
    return release
