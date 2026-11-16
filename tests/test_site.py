"""The committed site: its data against the release it is built from, and nothing private in the repository.

The release folder and the passphrase come from SUMSTATS_RELEASE and SUMSTATS_PASSPHRASE; a test that needs one that
is not set is skipped, and says so. Run with both before publishing.
"""
import csv
import functools
import getpass
import json
import math
import os
import pwd
import re
import subprocess
import tomllib
from pathlib import Path

import pytest
from cryptography.exceptions import InvalidTag

import numpy as np

from build import crypt, merge
from build.__main__ import BLOB, build
from build.release import tree
from conftest import DATA, ROOT, SITE
from papers import cite
from qc import labels

SITE_DATA = json.loads((DATA / "studies.json").read_text())
STUDIES, RELEASE = SITE_DATA["studies"], SITE_DATA["release"]
BY_ID = {d["id"]: d for d in STUDIES}
PAPERS = json.loads((DATA / "papers.json").read_text())
COLUMNS = {"rg": "genetic_correlation", "se": "genetic_correlation_se", "z": "genetic_correlation_z", "p": "genetic_correlation_p",
           "gcov": "genetic_covariance", "gcov_se": "genetic_covariance_se"}  # site name -> the release's
# What is never published: a path on O2 or on someone's computer, and a user name there (two to four letters, then
# two to four digits). This file names none either: it is committed too.
PATH = re.compile(r"/n/[a-z0-9_]+|/(?:home|Users)/\w")
NAME = r"[a-z]{2,4}[0-9]{2,4}"
USER = re.compile(rf"(?<![A-Za-z0-9_#.-]){NAME}(?![A-Za-z0-9_-])")
NOT_USERS = {"base64", "sha256", "win32"}
DOI = re.compile(r"\b10\.\d{4,9}/[^\s\"<>]+")  # a paper's DOI can end in the shape of a user name: it is not one


def correlations(anc):
    return json.loads((DATA / f"rg_{anc}.json").read_text())


def by_pair(t):
    """{the two IDs: the pair's values}: the lists hold the pairs i < j of `ids`, row by row."""
    order = [frozenset((a, b)) for i, a in enumerate(t["ids"]) for b in t["ids"][i + 1:]]
    return dict(zip(order, zip(*(t[c] for c in COLUMNS), strict=True), strict=True))


def rows(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def number(text):
    """A cell of a results table as the site holds it: a number, or None where the release has no estimate."""
    try:
        return None if math.isnan(float(text)) else float(text)
    except ValueError:
        return None


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def committed():
    """Every file that is committed or would be: {path: text}."""
    return {name: (ROOT / name).read_text(encoding="utf-8", errors="replace")
            for name in git("ls-files", "--cached", "--others", "--exclude-standard").splitlines() if (ROOT / name).is_file()}


def people():
    """O2 user names known here: whoever runs the tests and, where the release is at hand, the owners of its files."""
    folder, names = Path(os.environ.get("SUMSTATS_RELEASE") or ROOT), {getpass.getuser()}
    for uid in {path.stat().st_uid for path in [folder, *folder.parents, *(folder.rglob("*") if folder != ROOT else [])]}:
        try:
            names.add(pwd.getpwuid(uid).pw_name)
        except KeyError:  # an account that is gone
            pass
    return {name for name in names if re.fullmatch(NAME, name)}  # not root, nor the account of a build machine


# ---------------------------------------------------------------- the data


def test_counts_are_the_release_s():
    count = lambda **want: sum(all(d.get(k) == v for k, v in want.items()) for d in STUDIES)
    assert len(STUDIES) == len(BY_ID) == 602
    assert (count(anc="EUR"), count(anc="EAS")) == (546, 56)
    assert [sum("panel" in d for d in STUDIES if d["anc"] == anc) for anc in ("EUR", "EAS")] == [85, 17]
    assert count(signed=True) == count(cross=True) == 576 and count(release=2026) == 205 and count(release=2024) == 397
    for anc, signed, failed in (("EUR", 520, 821), ("EAS", 56, 55)):
        t = correlations(anc)
        assert t["ids"] == [d["id"] for d in STUDIES if d["anc"] == anc and d["cross"]] and len(t["ids"]) == signed
        assert {len(t[c]) for c in COLUMNS} == {signed * (signed - 1) // 2} and t["rg"].count(None) == failed
        assert RELEASE[anc] | {"worst": None} == RELEASE[anc] | {"worst": None, "studies": count(anc=anc), "signed": signed, "pairs": len(t["rg"]), "failed": failed}


def test_a_panel_is_what_about_says():
    for anc in RELEASE["ancestries"]:
        panel = sorted((d for d in STUDIES if d["anc"] == anc and "panel" in d), key=lambda d: d["panel"])
        estimates = by_pair(correlations(anc))
        within = [estimates[frozenset((a["id"], b["id"]))][0] for i, a in enumerate(panel) for b in panel[i + 1:]]
        assert [d["panel"] for d in panel] == list(range(1, len(panel) + 1))
        assert all(d["h2_z"] > 6 and d["cross"] for d in panel)
        assert all(rg is not None and rg ** 2 <= 0.1 for rg in within)


def test_a_heatmap_is_in_the_order_of_its_clustering():
    for anc in RELEASE["ancestries"]:
        t = correlations(anc)
        n, estimates = len(t["ids"]), by_pair(t)
        rg = np.array([[0 if a == b else estimates[frozenset((a, b))][0] or 0 for b in t["ids"]] for a in t["ids"]])  # no estimate: 0
        distance = 1 - np.minimum(np.abs(rg), 1)
        panel = [i for i, name in enumerate(t["ids"]) if "panel" in BY_ID[name]]
        assert t["heat"] == {"all": tree(distance), "panel": tree(distance[np.ix_(panel, panel)], panel)}
        for members, heat in ((list(range(n)), t["heat"]["all"]), (panel, t["heat"]["panel"])):
            heights = [height for *_, height in heat["merges"]]
            assert sorted(heat["order"]) == members and len(heights) == len(members) - 1 and heights == sorted(heights)


def test_every_study_has_a_short_name():
    names = [d["name"] for d in STUDIES]
    assert all(merge.NAME.fullmatch(name) for name in names) and len({name.lower() for name in names}) == len(names) == 602
    kept = cite.table(merge.DATA / "names.tsv")
    assert {row["id"]: row["name"] for row in kept} == {d["id"]: d["name"] for d in STUDIES} and all(row["how"] for row in kept)
    assert BY_ID["PASS_BMI1"]["name"] == "Speliotes_2010_BMI" and BY_ID["GCST90651074"]["name"] == "Kany_2025_AorticStenosis"  # the owner's examples
    assert sum(d["name"].endswith("EAS") for d in STUDIES) == sum(d["anc"] == "EAS" for d in STUDIES) == 56  # an East Asian study says so, and no other
    # a study of one sex says so, where the GWAS Catalog knows and the trait is one of both sexes; and where the release's own trait name or ID does
    assert BY_ID["GCST90566420"]["name"] == "Gu_2025_FastingGlucoseFemaleEAS" and BY_ID["GCST90503313"]["name"] == "Pagadala_2025_TotalTestosteroneMale"
    assert BY_ID["GCST90566419"]["name"] == "Gu_2025_GDMEAS" and BY_ID["GCST90483500"]["name"] == "Venkatesh_2025_PCOS" and BY_ID["GCST90301704"]["name"] == "Yang_2024_Preeclampsia"
    assert not [d["name"] for d in STUDIES if re.search(r"\(female\)|_Female$", d["trait"] + d["id"]) and "Female" not in d["name"]]
    # two names say what the file is, not what the release calls it
    assert BY_ID["PASS_PCT_Vuckovic2020"]["name"] == "Vuckovic_2020_Plateletcrit" and BY_ID["UKB_460K.repro_AgeFirstBirth_Female"]["name"] == "Loh_2018_AgeFirstBirthFemale"


def test_every_study_cites_its_paper():
    """A paper for each of the 602: the release's own for a study it gives a PMID, the one of data/references.tsv for the others."""
    assert all(d.get("cite") for d in STUDIES) and all(d["cite"] == [d["pmid"]] for d in STUDIES if "pmid" in d)
    cited = {pmid for d in STUDIES for pmid in d["cite"]}
    assert cited == set(PAPERS) == set(SITE_DATA["papers"]) and SITE_DATA["papers"] == {pmid: paper["label"] for pmid, paper in PAPERS.items()}
    keys = [paper["key"] for paper in PAPERS.values()]
    assert len(set(keys)) == len(keys) and all(re.fullmatch(r"[A-Za-z]+\d{4}[a-z]?", key) for key in keys)
    assert all(paper["bib"].startswith(f"@article{{{paper['key']},\n  author = {{") and paper["bib"].endswith(f"  pmid = {{{pmid}}}\n}}\n") for pmid, paper in PAPERS.items())
    assert all(paper["title"] and paper["journal"] and paper["authors"] and re.fullmatch(r"(19|20)\d\d", paper["year"]) for paper in PAPERS.values())
    # a study's short name begins with whom and when its paper is cited by, which is also the paper's key
    assert not [d["id"] for d in STUDIES if not PAPERS[d["cite"][0]]["key"].startswith(cite.plain("".join(d["name"].split("_")[:2])))]
    # UK Biobank's studies cite the paper of their statistics, then the resource; the Pan-UK Biobank's, that project's paper
    assert all(d["cite"] == ["29892013", "30305743"] for d in STUDIES if d["id"].startswith("UKB_460K")) and all(d["cite"] == ["40968291"] for d in STUDIES if d["id"].startswith("PANUKBB_"))
    references = cite.table(merge.DATA / "references.tsv")
    assert {row["reference"] for row in references} == {d["ref"] for d in STUDIES if "pmid" not in d} and not [row for row in references if "check" in row["how"] or not row["pmid"]]


def test_the_checks_catch_the_known_mislabels():
    """Four studies of the Pan-UK Biobank whose label and genetic correlations disagree: the label check finds each, by
    itself, from the trait names and the genetic correlations alone."""
    tables = {anc: correlations(anc) for anc in RELEASE["ancestries"]}
    signed = {d["id"]: {"trait": d["trait"], "h2_z": d["h2_z"]} for d in STUDIES if d["cross"]}
    found = labels.check(signed, tables, tomllib.loads((merge.DATA / "checks.toml").read_text())["label"])
    ldl = "PANUKBB_LDLDirectAdjustedByMedication"
    assert any(level == "fail" and "FEV1" in text for level, text in found[ldl])  # its best match is FEV1/FVC
    assert any(level == "fail" for level, _ in found["PANUKBB_IndirectBilirubin"])
    assert found["PANUKBB_OperationTonsillectomyAdenoids"] and found["PANUKBB_AtrialFibrillationAndFlutter"]
    # what the site says of them, after a person looked: three hold another trait; the fourth is what it says, and the
    # study it disagreed with is the one at fault, its Z not going with its alleles
    assert [BY_ID[study]["qc"] for study in (ldl, "PANUKBB_IndirectBilirubin", "PANUKBB_OperationTonsillectomyAdenoids")] == ["fail"] * 3
    assert BY_ID["PANUKBB_AtrialFibrillationAndFlutter"]["qc"] == "pass" and BY_ID["PASS_AtrialFibrillation_Nielsen2018"]["qc"] == "fail"


def test_what_the_checks_say_of_the_release():
    counts = {level: sum(d["qc"] == level for d in STUDIES) for level in ("pass", "warn", "fail")}
    assert RELEASE["qc"] == counts and sum(counts.values()) == 602 and counts["pass"] > 500  # the thresholds are set so that a study with nothing wrong passes
    for d in STUDIES:
        levels = [level for level, _ in d.get("notes", [])]
        assert set(levels) <= {"fail", "warn", "note"} and levels == sorted(levels, key=["fail", "warn", "note"].index)  # the worst first
        assert d["qc"] == ("fail" if "fail" in levels else "warn" if "warn" in levels else "pass")
        assert all(study in BY_ID and study != d["id"] for _, text in d.get("notes", []) for study in re.findall(r"\{([^}]+)\}", text))  # a finding names studies as {ID}
        assert 0 < d["hm3"] <= 1 and d["chi2"] > 0 and d["lambda"] > 0
    # the largest GWAS of height passes; so does a classic whose statistics were corrected for inflation, with a note
    assert BY_ID["PASS_Height_Yengo2022"]["qc"] == "pass" and BY_ID["PASS_BMI1"]["qc"] == "pass" and BY_ID["PASS_BMI1"]["notes"][0][0] == "note"
    # the two labels that a person found wrong warn, and say what the file holds
    for study, holds in (("PASS_PCT_Vuckovic2020", "Plateletcrit"), ("UKB_460K.repro_AgeFirstBirth_Female", "Age first birth")):
        assert BY_ID[study]["qc"] == "warn" and any(text.startswith("the release calls it") and holds in text for _, text in BY_ID[study]["notes"])
    metrics = cite.table(merge.DATA / "metrics.tsv")
    assert {row["id"] for row in metrics if row["file"] == "h2"} == set(BY_ID) and {row["id"] for row in metrics if row["file"] == "rg"} < set(BY_ID)
    assert not [d["id"] for d in STUDIES for _, text in d.get("notes", []) if "was not measured" in text or "SNPs in the file; the inventory says" in text]


def test_a_problem_is_reported_through_a_form_that_warns_it_is_public():
    form, script = (ROOT / ".github" / "ISSUE_TEMPLATE" / "study.yml").read_text(), (SITE / "app.js").read_text()
    assert "id: study" in form and "public" in form and 'template: "study.yml"' in script and "study: `${d.name} (${d.id})`" in script  # the page fills in the study
    nominate = (ROOT / ".github" / "ISSUE_TEMPLATE" / "candidate.yml").read_text()
    assert "public" in nominate and "issues/new?template=candidate.yml" in (SITE / "index.html").read_text()  # and a study is nominated through another


def test_the_data_is_small():
    sizes = {path.name: path.stat().st_size for path in DATA.iterdir()}
    assert set(sizes) == {"studies.json", "papers.json", "candidates.json", "rg_EUR.json", "rg_EAS.json", BLOB}
    assert sum(sizes.values()) < 20e6  # the plan: well under 100 MB


def test_the_site_data_is_a_build_of_the_release(release, tmp_path):
    build(release, tmp_path, "not the passphrase", 1000)
    differ = [name for name in ("studies.json", "papers.json", "candidates.json", "rg_EUR.json", "rg_EAS.json") if (tmp_path / name).read_bytes() != (DATA / name).read_bytes()]
    assert not differ, "rebuild the site: python -m build"


def test_every_value_is_the_release_s(release):
    """Read again from the release's own tables, without the build's code."""
    one = lambda pattern: functools.reduce(lambda _, path: path, sorted(release.glob(pattern)))
    yes = lambda cell: {"YES": True, "NO": False}[cell]
    inventory = rows(one("summary_statistics_*.txt"))
    assert [(r["Trait Identifier"], r["Trait Name"].strip(), r["Population"], int(r["N"]), int(r["M"]), r["Origin"].strip(), r["Reference"].strip(),
             yes(r["Signed SumStat"]), yes(r["Allele Info"]), yes(r["Cross-Trait LDSC"])) for r in inventory] \
        == [(d["id"], d["trait"], d["anc"], d["n"], d["m"], d["origin"], d["ref"], d["signed"], d["alleles"], d["cross"]) for d in STUDIES]
    for anc in ("EUR", "EAS"):
        mine = [d for d in STUDIES if d["anc"] == anc]
        h2 = rows(one(f"sldsc_heritabilities_{anc}_*.txt"))
        assert {r["trait_identifier"]: tuple(float(r[c]) for c in ("h2", "h2_se", "h2_z", "intercept")) for r in h2} \
            == {d["id"]: (d["h2"], d["h2_se"], d["h2_z"], d["intercept"]) for d in mine}
        assert [r["trait_identifier"] for r in rows(one(f"non_redundent_traits_{anc}_*.txt"))] \
            == [d["id"] for d in sorted((d for d in mine if "panel" in d), key=lambda d: d["panel"])]
        table, t = rows(one(f"ldsc_genetic_correlations_trait_pairs_{anc}_*.txt")), correlations(anc)
        theirs = {frozenset((r["trait1"], r["trait2"])): tuple(number(r[c]) for c in COLUMNS.values()) for r in table}
        ours = by_pair(t)
        assert len(table) == len(theirs) == len(ours)
        assert not [sorted(pair) for pair in theirs if theirs[pair] != ours.get(pair)][:5]
        reasons = {frozenset((r["trait1"], r["trait2"])): " ".join(r["genetic_correlation_z"].split("(")[-1].rstrip(")").split())
                   for r in table if "out of bounds" in r["genetic_correlation_z"]}
        assert reasons == {pair: t["notes"][str(k)] for k, pair in enumerate(ours) if str(k) in t["notes"]} and len(reasons) == len(t["notes"])
    index = rows(release / "_audit" / "sumstats_file_index.tsv")
    assert [(r["Trait Identifier"], int(r["Release"][:4]), int(r["M in RG input"]) if r["M in RG input"] else None) for r in index] \
        == [(d["id"], d["release"], d.get("m_rg")) for d in STUDIES]
    papers = rows(release / "_audit" / "new_study_provenance.tsv")
    assert {r["Trait Identifier"]: (r["PMID"], r["Paper Link"], r["Catalog Study Link"], r["Target Paper Identifier"]) for r in papers} \
        == {d["id"]: (d["pmid"], d["paper"], d["catalog"], d["paper_id"]) for d in STUDIES if "pmid" in d}
    assert len(papers) == RELEASE["all"]["added"] and len({r["Target Paper Identifier"] for r in papers}) == RELEASE["papers"]


def test_the_texts_quote_the_release(release):
    notes = re.sub(r"<[^>]+>|\s+", " ", (SITE / "index.html").read_text())
    kept = {r["Check"]: int(r["Details"]) for r in rows(release / "_audit" / "validation_results.tsv") if r["Check"].endswith("complete subgraph")}
    eur, eas = kept["EUR complete subgraph"], kept["EAS complete subgraph"]
    assert re.sub(r"\s+", " ", f"{RELEASE['EUR']['signed'] - eur} EUR studies and {RELEASE['EAS']['signed'] - eas} EAS study go, {eur} and {eas} stay.") in re.sub(r"\s+", " ", notes)


def test_the_numbers_in_the_texts_are_the_data_s():
    """The page fills them in from the data when it loads; what is written in must agree, for a reader of the file."""
    html = (SITE / "index.html").read_text()
    written = re.findall(r'data-v="([\w.]+)">([^<]*)<', html)
    value = lambda path: functools.reduce(lambda o, k: o[k], path.split("."), RELEASE)
    shown = lambda v: f"{v:,}" if isinstance(v, int) else BY_ID[v]["name"] if v in BY_ID else v.replace("-", "/")  # a study by its short name
    assert len(written) > 40
    assert not [(path, text) for path, text in written if text != shown(value(path))]
    named = re.findall(r'data-study="([^"]+)">([^<]*)<', html)
    assert len(named) >= 2 and not [(study, text) for study, text in named if text != BY_ID[study]["name"]]
    corrected = BY_ID["GCST90503316"]
    assert re.search(rf'<span data-study="GCST90503316">{corrected["name"]}</span> uses a corrected\s+estimate \(h² {corrected["h2"]}, SE {corrected["h2_se"]}\)', html)
    assert f"{RELEASE['year']} release" in html and "Gaspard Kerner" in html and "Ben Strober" in html


# ---------------------------------------------------------------- what is private


def test_the_blob_opens_with_the_passphrase_only(passphrase):
    blob = json.loads((DATA / BLOB).read_text())
    assert (blob["kdf"], blob["cipher"]) == ("PBKDF2-HMAC-SHA256", "AES-256-GCM") and blob["iterations"] >= 600_000
    assert [len(crypt.base64.b64decode(blob[k])) for k in ("salt", "iv")] == [16, 12]
    private = json.loads(crypt.decrypt(blob, passphrase))
    assert set(private["files"]) == set(BY_ID) and all(private["folders"][folder] and name for folder, name in private["files"].values())
    for wrong in (passphrase + " ", passphrase[:-1], passphrase.upper(), ""):
        with pytest.raises(InvalidTag):
            crypt.decrypt(blob, wrong)


def test_the_blob_holds_the_file_index_and_every_file_is_there(release, passphrase, tmp_path):
    private = json.loads(crypt.decrypt(json.loads((DATA / BLOB).read_text()), passphrase))
    build(release, tmp_path, "not the passphrase", 1000)
    assert private == json.loads(crypt.decrypt(json.loads((tmp_path / BLOB).read_text()), "not the passphrase")), "rebuild the site: python -m build"
    where = lambda file: f"{private['folders'][file[0]]}/{file[1]}"
    index = rows(release / "_audit" / "sumstats_file_index.tsv")
    assert [(r["Trait Identifier"], r["H2 Sumstats"]) for r in index] == [(d["id"], where(private["files"][d["id"]])) for d in STUDIES]
    copies = {r["RG Sumstats"] for r in index} - {r["H2 Sumstats"] for r in index} - {""}  # what rg was estimated on, where it is another file
    assert copies and not {path.rpartition("/")[0] for path in copies} & set(private["folders"])  # one path per study: these are not listed
    assert private["release"] == str(release) and all((release / name).is_file() for name in private["tables"])
    head = lambda path: Path(path).open().readline().rstrip("\n").split("\t")  # About: "What is in a file?"
    wrong = [r["Trait Identifier"] for r in index for path in {r["H2 Sumstats"], r["RG Sumstats"]} - {""}
             if head(path) != (["SNP", "A1", "A2", "N", "Z"] if BY_ID[r["Trait Identifier"]]["alleles"] else ["SNP", "N", "Z"])]
    assert not wrong


def test_no_path_and_no_user_name_is_committed():
    known, files = people(), committed()
    assert len(files) > 10
    del files[f"site/data/{BLOB}"]  # ciphertext: a few of its letters, digits and slashes can look like a path or a name
    paths = {name: PATH.findall(text)[:3] for name, text in files.items() if PATH.search(text)}
    shaped = {name: sorted(set(USER.findall(DOI.sub("", text))) - NOT_USERS)[:5] for name, text in files.items()}
    named = [name for name, text in files.items() if any(person in text for person in known)]
    assert not paths and not {name: found for name, found in shaped.items() if found} and not named


def test_the_passphrase_and_the_paths_are_nowhere_but_in_the_blob(passphrase):
    files = committed()
    private = json.loads(crypt.decrypt(json.loads(files[f"site/data/{BLOB}"]), passphrase))
    release = Path(private["release"])
    secret = [str(release), release.name, str(Path(*release.parts[:3])), *private["folders"]]  # the folder, its name, the root it is under
    leaks = [name for name, text in files.items() if passphrase.lower() in text.lower() or any(s in text for s in secret)]
    assert not leaks  # the names of the files only: a failure must not print the passphrase


def test_the_history_is_clean():
    """Nothing private in an earlier commit or in a commit message either."""
    if not (ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    # every line a commit added or removed, and every message; --full-index: a short blob hash can take the shape of a user name
    everything = git("log", "--all", "-p", "--full-index", "--format=%B")
    readable = git("log", "--all", "-p", "--full-index", "--format=%B", "--", ".", f":(exclude)site/data/{BLOB}")  # without the ciphertext
    assert not PATH.findall(readable)[:3]
    assert not sorted(set(USER.findall(DOI.sub("", readable))) - NOT_USERS)[:5]
    assert not [k for k, person in enumerate(sorted(people())) if person in readable]
    passphrase = os.environ.get("SUMSTATS_PASSPHRASE")
    leaked = bool(passphrase) and passphrase.lower() in everything.lower()
    assert not leaked  # said without printing it


def test_a_row_is_a_study():
    """Never a "dataset": in the page's text, and in the script but for the old address ?dataset=ID and the DOM's own word."""
    html, script = ((SITE / name).read_text().lower() for name in ("index.html", "app.js"))
    assert "dataset" not in html and "dataset" not in re.sub(r'\.dataset\b|"dataset"', "", script)
    assert 'q.set("study", q.get("dataset"))' in script  # the old address still opens the study


def test_the_page_loads_nothing_from_elsewhere():
    """The passphrase is typed into this page: no script, style or request of another site."""
    html, script, style = ((SITE / name).read_text() for name in ("index.html", "app.js", "style.css"))
    assert all(rule in html for rule in ("default-src 'none'", "script-src 'self'", "style-src 'self'", "connect-src 'self'", "form-action 'none'"))
    loaded = re.findall(r'<(?:script|link|img|iframe|source)\b[^>]*?\b(?:src|href)="([^"]+)"', html)
    assert len(loaded) == 4 and not [url for url in loaded if re.match(r"(?:[a-z]+:)?//", url)]  # the icon (data:), the style, the inventory, the script
    assert re.findall(r"https?://[\w./-]+", script) == ["https://pubmed.ncbi.nlm.nih.gov/", "https://doi.org/", "https://github.com/gwassumstats/gwassumstats.github.io/issues/new", "https://www.ebi.ac.uk/gwas/studies/",
            "http://www.w3.org/2000/svg"]  # four links and a name: no request
    assert "url(" not in style and "@import" not in style and " style=" not in html
