"""How a paper is cited: who and which year for its label, its BibTeX entry and the entry's key."""
import csv
import re
import unicodedata
from collections import Counter
from pathlib import Path

SUFFIXES = {"Jr", "Sr", "2nd", "3rd", "4th", "II", "III", "IV"}
NAMED = re.compile(r"([A-Z][\w' -]*?)\s+et\.? al")           # Speliotes et al., 2010 Nat Genet
CODED = re.compile(r"[A-Za-z]+_([A-Za-z][A-Za-z'-]*)_\d{4}(?:_.*)?")  # NG_Kanai_2018_BBJ
YEAR = re.compile(r"\b(?:19|20)\d{2}\b")
PREPRINT = re.compile(r"[a-z]{2,3}rxiv", re.I)


def table(path: Path) -> list[dict]:
    """The rows of a tab-delimited file of data/, without its comment lines."""
    with path.open(encoding="utf-8", newline="") as file:
        return [row for row in csv.DictReader(file, delimiter="\t", quoting=csv.QUOTE_NONE) if not next(iter(row.values())).startswith("#")]


def person(name: str) -> tuple[str, str, str]:
    """A PubMed author, "Taylor HA Jr", as surname, initials and suffix; a name of another shape whole, as the surname."""
    *front, last = name.split()
    suffix = ""
    if last in SUFFIXES and len(front) > 1:
        suffix, last = last, front.pop()
    return (" ".join(front), last, suffix) if front and last.isalpha() and last.isupper() and len(last) <= 5 else (name, "", "")


def who(paper: dict, references: list[str], groups: dict[str, str]) -> tuple[str, bool]:
    """Who a paper is cited by, and whether "et al." follows: its first author's surname. A paper that PubMed files
    under a group (in braces) is cited by the group's short name (data/consortia.tsv), or else by the person that the
    inventory's reference names, or else by the group's name."""
    first = paper["authors"][0] if paper["authors"] else ""
    if not first.startswith("{"):
        return person(first)[0] or "Anonymous", len(paper["authors"]) > 1
    if first.strip("{}") in groups:
        return groups[first.strip("{}")], False
    named = next((found[1] for found in (NAMED.match(r) or CODED.fullmatch(r) for r in references) if found), None)
    return (named, True) if named else (first.strip("{}"), False)


def year(paper: dict, references: list[str]) -> str:
    """The year a paper is cited by: the one the inventory gives it, so that a paper keeps the year the lab knows it by
    (online in 2025, in print in 2026: 2025). PubMed's year if the inventory names none, names a preprint, or names
    several and PubMed's among them; of several others, the most frequent, then the earliest."""
    named = Counter(found[0] for found in (YEAR.search(r) for r in references if not PREPRINT.search(r)) if found)
    if not named or paper["year"] in named:
        return paper["year"]
    return min(named, key=lambda y: (-named[y], y))


def label(paper: dict, references: list[str], groups: dict[str, str]) -> str:
    """Speliotes et al. 2010 Nat Genet; PGC 2014 Nature."""
    name, others = who(paper, references, groups)
    return f"{name}{' et al.' if others else ''} {year(paper, references)} {paper['journal']}"


def plain(text: str) -> str:
    """Letters and digits only, without accents: Ortín Vela -> OrtinVela."""
    return re.sub(r"[^A-Za-z0-9]", "", unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode())


def keys(cited: dict[str, tuple[str, str]]) -> dict[str, str]:
    """A BibTeX key for each paper, {PMID: key}, from who and which year it is cited by ({PMID: (who, year)}):
    Speliotes2010; where several papers share one, a letter in the order of their PMIDs: Jansen2019a, Jansen2019b."""
    base = {pmid: plain(name) + when for pmid, (name, when) in cited.items()}
    shared = Counter(base.values())
    order = {pmid: k for key in shared for k, pmid in enumerate(sorted((p for p in base if base[p] == key), key=int))}
    return {pmid: key + ("abcdefghijklmnopqrstuvwxyz"[order[pmid]] if shared[key] > 1 else "") for pmid, key in base.items()}


def pages(text: str) -> str:
    """PubMed's pages as BibTeX writes a range: 937-48 -> 937--948; 1214-1231.e11 -> 1214--1231.e11."""
    short = re.fullmatch(r"(\d+)-(\d+)", text)
    if short and len(short[2]) < len(short[1]):
        return f"{short[1]}--{short[1][:-len(short[2])]}{short[2]}"
    return text.replace("-", "--", 1)


def bibtex(key: str, pmid: str, paper: dict) -> str:
    """The paper as a BibTeX entry. An author is "Surname, Suffix, I. N."; a group stays whole, in braces."""
    def author(name: str) -> str:
        if name.startswith("{"):
            return name
        surname, initials, suffix = person(name)
        return ", ".join(part for part in (surname, suffix, " ".join(f"{letter}." for letter in initials)) if part)

    text = lambda value: re.sub(r"([&%$#_])", r"\\\1", value.replace("{", "").replace("}", ""))
    fields = {"author": " and ".join(map(author, paper["authors"])), "title": "{" + text(paper["title"]) + "}", "journal": text(paper["journal"]),
              "year": paper["year"], "volume": paper["volume"], "number": paper["issue"], "pages": pages(paper["pages"]), "doi": paper["doi"], "pmid": pmid}
    return f"@article{{{key},\n" + ",\n".join(f"  {name} = {{{value}}}" for name, value in fields.items() if value) + "\n}\n"
