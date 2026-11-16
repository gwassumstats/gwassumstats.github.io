"""What the repository adds to a release, from data/: the paper of each study, and its short name."""
import json
import re
from pathlib import Path

import pandas as pd

from papers import cite

DATA = Path(__file__).parents[1] / "data"
NAMED = 3  # the authors of a paper that the site names before "et al."
NAME = re.compile(r"[A-Za-z]+_\d{4}_[A-Za-z0-9]+")  # a short name: FirstAuthor_Year_Trait


def sources(d: pd.DataFrame, data: Path) -> tuple[list[list[str]], dict, dict, dict]:
    """The papers each study cites, as PMIDs, the paper of its statistics first; PubMed's record of each paper; for
    each paper, how the inventory writes it for the studies whose statistics are its own; and the groups' short names.

    A study that the release gives a PMID cites that paper. Another cites what data/references.tsv says of its
    reference, for this study or for all. A row there that still says "check" is left out: its studies cite nothing yet."""
    rows = {(row["reference"], row["study"]): [pmid for pmid in row["pmid"].split(";") if pmid]
            for row in cite.table(data / "references.tsv") if "check" not in row["how"]}
    cited = [[pmid] if pd.notna(pmid) else rows.get((ref, study)) or rows.get((ref, "")) or [] for study, pmid, ref in zip(d.index, d["pmid"], d["ref"])]
    records = json.loads((data / "papers.json").read_text(encoding="utf-8"))
    missing = sorted({pmid for pmids in cited for pmid in pmids} - set(records))
    if missing:
        raise ValueError(f"data/papers.json lacks the paper {missing[0]}: run python -m papers")
    references = {pmid: [] for pmids in cited for pmid in pmids}
    for pmids, ref in zip(cited, d["ref"]):
        if pmids:
            references[pmids[0]].append(ref)
    return cited, records, references, {row["group"]: row["short"] for row in cite.table(data / "consortia.tsv")}


def citations(d: pd.DataFrame, data: Path) -> tuple[list[list[str]], dict[str, dict]]:
    """The papers each study cites (see sources), and what the site shows of every paper cited, {PMID: ...}, by PMID:
    its label, the first authors, the bibliographic record, and its BibTeX entry under a key of its own."""
    cited, records, references, groups = sources(d, data)
    keys = cite.keys({pmid: (cite.who(records[pmid], refs, groups)[0], cite.year(records[pmid], refs)) for pmid, refs in references.items()})
    papers = {}
    for pmid in sorted(references, key=int):
        paper = records[pmid]
        papers[pmid] = {"key": keys[pmid], "label": cite.label(paper, references[pmid], groups), "authors": paper["authors"][:NAMED], "more": max(0, len(paper["authors"]) - NAMED),
                        **{field: paper[field] for field in ("title", "journal", "year", "volume", "issue", "pages", "doi")}, "bib": cite.bibtex(keys[pmid], pmid, paper)}
    return cited, papers


def names(d: pd.DataFrame, data: Path) -> pd.Series:
    """The short name of each study, from data/names.tsv. Every study has one, of the form FirstAuthor_Year_Trait,
    and no two studies the same, whatever the case of its letters."""
    rows = cite.table(data / "names.tsv")
    given = pd.Series({row["id"]: row["name"] for row in rows}, dtype=object)
    if len(given) != len(rows):
        raise ValueError("data/names.tsv names a study twice")
    lacking = d.index.difference(given.index)
    if len(lacking):
        raise ValueError(f"data/names.tsv lacks the short name of {lacking[0]}: run python -m names")
    named = given.reindex(d.index)
    odd, twice = named[~named.str.fullmatch(NAME.pattern)], named[named.str.lower().duplicated()]
    if len(odd):
        raise ValueError(f"data/names.tsv: {odd.iloc[0]} is not FirstAuthor_Year_Trait")
    if len(twice):
        raise ValueError(f"data/names.tsv: {twice.iloc[0]} names two studies")
    return named
