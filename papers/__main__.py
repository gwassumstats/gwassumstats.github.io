"""Finds the paper of every study of a release in PubMed, and keeps what it found where it can be read and corrected.

  uv run python -m papers [RELEASE]

RELEASE defaults to $SUMSTATS_RELEASE. A study that the release gives a PMID keeps it. For the others, each reference
of the inventory ("Speliotes et al., 2010 Nat Genet", "NG_Kanai_2018_BBJ") is looked up once: first author, year and
journal; of several papers, the one whose title names the study's trait. What is found goes to

  data/references.tsv   reference (and study, where one reference stands for several papers) -> PMID, DOI, and how it
                        was matched. A row that says "check" needs a look. Correct a row by hand and write "by hand"
                        in `how`: a later run leaves the rows that are there alone, and looks up only new references.
  data/papers.json      the PubMed record of every paper: authors, title, journal, year, volume, pages, DOI

The build merges both into the site's data. Rows that say "check" are left out until someone has looked.
"""
import argparse
import csv
import json
import os
import re
import urllib.error
from pathlib import Path

from build import release as rel
from papers import cite, pubmed

ROOT = Path(__file__).parents[1]
REFERENCES, PAPERS = ROOT / "data" / "references.tsv", ROOT / "data" / "papers.json"
FIELDS = ["reference", "study", "pmid", "doi", "how", "paper"]
JOURNALS = {  # as the inventory writes a journal -> as PubMed abbreviates it
    "nat genet": "Nat Genet", "nat gent": "Nat Genet", "ng": "Nat Genet", "nat commun": "Nat Commun", "nat comm": "Nat Commun", "nc": "Nat Commun",
    "nature": "Nature", "science": "Science", "cell": "Cell", "lancet": "Lancet", "cell genomics": "Cell Genom", "cell genom": "Cell Genom",
    "ajhg": "Am J Hum Genet", "mol psych": "Mol Psychiatry", "nat neuro": "Nat Neurosci", "plos genet": "PLoS Genet", "pg": "PLoS Genet",
    "nat aging": "Nat Aging", "nat hum behav": "Nat Hum Behav", "ejhg": "Eur J Hum Genet", "the lancet oncology": "Lancet Oncol",
    "the lancet psychiatry": "Lancet Psychiatry", "lancet neurolo": "Lancet Neurol", "translational psychiatry": "Transl Psychiatry",
    "am j psychiatry": "Am J Psychiatry", "american journal of psychiatry": "Am J Psychiatry", "annals of the rheumatic diseases": "Ann Rheum Dis",
    "human molecular genetics": "Hum Mol Genet", "hummolgenet": "Hum Mol Genet", "diabetes": "Diabetes", "eneuro": "eNeuro", "elife": "Elife",
    "carcinogenesis": "Carcinogenesis"}
GWAS = ("genome-wide", "genome wide", "gwas", "loci", "locus", "association", "genetic", "variant", "heritab", "polygenic", "genomic")
DULL = {"count", "level", "levels", "volume", "disease", "diseases", "disorder", "average", "total", "percentage", "with", "from", "type", "rate", "score", "study"}


def parse(reference: str) -> tuple[str, str, str]:
    """First author, year and PubMed's journal of a reference of the inventory, each "" if it names none."""
    coded = re.fullmatch(r"([A-Za-z]+)_([A-Za-z'-]+)_(\d{4})(?:_.*)?", reference)  # NG_Kanai_2018_BBJ
    if coded:
        return coded[2], coded[3], JOURNALS.get(coded[1].lower(), "")
    year = re.search(r"\b(?:19|20)\d{2}\b", reference)
    named = re.match(r"(.+?)\s+et\.? al\.?,?", reference)  # Speliotes et al., 2010 Nat Genet
    author = named[1] if named else ""
    rest = re.sub(r"\b(?:19|20)\d{2}\b", "", reference[len(named[0]):] if named else reference.partition(",")[2])
    return author, year[0] if year else "", JOURNALS.get(re.sub(r"[^a-z ]", "", rest.lower()).strip(), "")


def stems(text: str) -> set[str]:
    """The words of a trait or a title that say what it is about, cut to five letters: "Osteoarthritis in Hip" -> osteo."""
    return {word[:5] for word in re.findall(r"[a-z]{4,}", text.lower()) if word not in DULL}


def choose(traits: list[str], found: dict[str, dict]) -> tuple[str, str]:
    """Of the papers found, the one whose title shares the most words with the traits, and how sure that is."""
    about = set().union(*map(stems, traits))
    shared = {pmid: len(about & stems(paper["title"])) for pmid, paper in found.items()}
    best = max(shared, key=shared.get)
    looks = any(word in found[best]["title"].lower() for word in GWAS)
    if len(found) == 1:
        return best, "the only paper" + (": its title names the trait" if shared[best] else ": a genetic study" if looks else "; check: its title does not name the trait")
    alone = sorted(shared.values())[-2] < shared[best]
    return best, f"of {len(found)} papers, " + ("the one whose title names the trait" if alone and shared[best] else "check: none stands out by its title")


def resolve(reference: str, studies: list[dict]) -> list[dict]:
    """The rows of data/references.tsv for a reference: one for all its studies or, if they are of several papers, one a study."""
    author, year, journal = parse(reference)
    if not author:
        return [{"reference": reference, "study": "", "pmid": "", "doi": "", "how": "check: no first author to look up", "paper": ""}]
    term = f"{author}[1au]" + (f" AND {year}[dp]" if year else "") + (f' AND "{journal}"[ta]' if journal else "")
    by = "author, year and journal" if year and journal else "author and year" if year else "author and journal"
    papers = lambda term: {pmid: paper for pmid, paper in pubmed.records(pubmed.search(term)).items()
                           if not re.match(r"(author |publisher )?(correction|erratum|corrigendum)", paper["title"].lower())}
    found = papers(term)
    if not found and journal and year:  # the journal may be another one: Nalls 2019
        found, by = papers(f"{author}[1au] AND {year}[dp]"), "check: author and year, another journal"
    if not found:
        return [{"reference": reference, "study": "", "pmid": "", "doi": "", "how": f"check: nothing for {term}", "paper": ""}]
    picks = {s["id"]: choose([s["trait"]], found) for s in studies}
    whole = choose([s["trait"] for s in studies], found)
    row = lambda study, pmid, how: {"reference": reference, "study": study, "pmid": pmid, "doi": found[pmid]["doi"], "how": f"{by}: {how}" if "check" not in by else f"{by}; {how}",
                                    "paper": label(found[pmid])}
    if len({pmid for pmid, _ in picks.values()}) == 1 or any("check" in how for _, how in picks.values()):
        return [row("", *whole)]
    return [row(study, *pick) for study, pick in picks.items()]  # Jansen et al. 2019 Nat Genet: Alzheimer's disease, and insomnia


def label(paper: dict) -> str:
    return f"{paper['authors'][0] if paper['authors'] else '?'} {paper['year']} {paper['journal']}: {paper['title']}"


def main():
    parser = argparse.ArgumentParser(prog="python -m papers", description="Finds the paper of every study of a release in PubMed.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    d = rel.studies(Path(args.release))
    rows = cite.table(REFERENCES) if REFERENCES.exists() else []
    known = {row["reference"] for row in rows}
    open_references = {}
    for study, trait, ref in d.loc[d["pmid"].isna(), ["trait", "ref"]].itertuples(name=None):
        if ref not in known:
            open_references.setdefault(ref, []).append({"id": study, "trait": trait})
    for k, (reference, group) in enumerate(sorted(open_references.items())):
        try:
            new = resolve(reference, group)
        except (urllib.error.URLError, TimeoutError) as error:
            new = [{"reference": reference, "study": "", "pmid": "", "doi": "", "how": f"check: PubMed did not answer ({error})", "paper": ""}]
        rows += new
        print(f"[{k + 1}/{len(open_references)}] {reference}: {'; '.join(row['how'] for row in new)}", flush=True)
    rows.sort(key=lambda row: (row["reference"].lower(), row["study"]))
    with REFERENCES.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, FIELDS, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None)
        writer.writeheader()
        writer.writerows(rows)
    pmids = sorted({pmid for row in rows for pmid in row["pmid"].split(";") if pmid} | set(d["pmid"].dropna()), key=int)
    PAPERS.write_text(json.dumps(pubmed.records(pmids), ensure_ascii=False, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    look = [row for row in rows if "check" in row["how"]]
    print(f"{len(rows)} rows for {len({row['reference'] for row in rows})} references: {len(rows) - len(look)} matched, {len(look)} to check; {len(pmids)} papers")


if __name__ == "__main__":
    main()
