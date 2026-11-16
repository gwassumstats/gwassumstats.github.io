"""Drafts a short name for every study that has none: FirstAuthor_Year_Trait, as Speliotes_2010_BMI.

  uv run python -m names [RELEASE]

RELEASE defaults to $SUMSTATS_RELEASE. The first author and the year are those the study's paper is cited by
(papers/cite.py; run python -m papers first). The trait is the release's trait name, shortened with
data/abbreviations.tsv. A study of women or of men only says so, Female or Male, unless its trait already does: the
GWAS Catalog's record of the study tells (data/catalog.tsv, if python -m candidates wrote it). An East Asian study
ends in EAS.

The names are kept in data/names.tsv: study ID, short name, and how the name was made. The build merges them into the
site's data after checking that every study has one, of that form, and no two the same. A row that says "check" was a
guess. Correct any name by hand and write "by hand" in `how`: a later run keeps the rows that are there, and drafts
names for new studies only.

A short name is what the site calls a study. It names no file: the study's ID stays its key.
"""
import argparse
import csv
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

from build import merge, release as rel
from papers import cite

DATA = merge.DATA
SMALL = {"of", "the", "and", "in", "at", "for", "to", "by", "with", "a", "an", "or", "on", "level", "levels"}  # words a trait's short form leaves out


def trait(name: str, phrases: dict[str, str]) -> str:
    """A trait name as one word. A name that is a phrase of the list is what stands for it. A name that ends in its
    own code is that code: Hemoglobin (HGB) -> HGB. In any other, the phrases of the list are replaced, the longest
    first, and the other words kept, each with a capital: Fasting Glucose -> FastingGlucose."""
    text = name.replace("’", "'")
    low = text.lower()
    if low in phrases:
        return cite.plain(phrases[low])
    code = re.fullmatch(r".*\((\w*[A-Z0-9]\w*)\)", text)
    if code and f"({code[1].lower()})" not in phrases:
        return cite.plain(code[1])
    spans = []  # start, end and short form of each phrase found, none inside another
    for phrase in sorted(phrases, key=len, reverse=True):
        for found in re.finditer(rf"(?<![a-z0-9]){re.escape(phrase)}(?![a-z0-9])", low):
            if all(found.end() <= start or found.start() >= end for start, end, _ in spans):
                spans.append((found.start(), found.end(), phrases[phrase]))
    words = lambda piece: [word if re.search(r"[A-Z0-9]", word[1:]) else word.capitalize() for word in re.findall(r"[A-Za-z0-9]+(?:'[A-Za-z]+)*", piece) if word.lower() not in SMALL]
    out, at = [], 0
    for start, end, short in sorted(spans):
        out += [*words(text[at:start]), short]
        at = end
    return cite.plain("".join([*out, *words(text[at:])]))


def one_sex(record: dict) -> str:
    """Whose statistics a study holds, by the GWAS Catalog's record of it: "Female", "Male", or "" for both or unknown.
    The submitters' word goes first; without it, a sample of female or of male cases and controls tells."""
    told = {"F": "Female", "M": "Male", "combined": ""}.get(record.get("sex", ""))
    if told is not None:
        return told
    sample = record.get("sample", "").lower()
    return "Female" if re.search(r"\b(female|women)", sample) else "Male" if re.search(r"\b(male|men)\b|\bmales\b", sample) else ""


def draft(studies: list[dict], cited: dict[str, tuple[str, str] | None], phrases: dict[str, str], catalog: dict[str, dict]) -> dict[str, tuple[str, str]]:
    """{study ID: (short name, how it was made)}. `studies`: id, trait, anc, n and ref of each; `cited`: who and which
    year each study's paper is cited by, or None for a study without a paper; `catalog`: the GWAS Catalog's `sex` and
    `sample` of a study, where it knows them."""
    by_id, parts, how = {s["id"]: s for s in studies}, {}, {}
    for s in studies:
        paper = cited.get(s["id"])
        named = cite.NAMED.match(s["ref"]) or cite.CODED.fullmatch(s["ref"])
        who, year = paper or (named[1] if named else "Unknown", (cite.YEAR.search(s["ref"]) or ["0000"])[0])
        short, sex = trait(s["trait"], phrases), one_sex(catalog.get(s["id"], {}))
        says = bool(sex) and not re.search("Female|Male", short)  # a name never hides that a study is of one sex
        parts[s["id"]] = [cite.plain(who), year, short + (sex if says else "") + ("EAS" if s["anc"] == "EAS" else "")]
        how[s["id"]] = ("drafted" + (f": of {'women' if sex == 'Female' else 'men'} only, by the GWAS Catalog" if says else "")) if paper else "check: no paper found"
    same = defaultdict(list)
    for study, name in parts.items():
        same["_".join(name).lower()].append(study)
    for group in (g for g in same.values() if len(g) > 1):  # studies of one paper and one trait: nothing tells them apart but their size
        for study in sorted(group, key=lambda study: -by_id[study]["n"])[1:]:
            parts[study][2] += f"N{round(by_id[study]['n'] / 1000)}k"
            how[study] = "check: told apart by its sample size only"
    return {study: ("_".join(name), how[study]) for study, name in parts.items()}


def main():
    parser = argparse.ArgumentParser(prog="python -m names", description="Drafts a short name for every study of a release that has none.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    d = rel.studies(Path(args.release))
    cited, records, references, groups = merge.sources(d, DATA)
    by = {pmid: (cite.who(records[pmid], refs, groups)[0], cite.year(records[pmid], refs)) for pmid, refs in references.items()}
    studies = [{"id": study, **row} for study, row in d[["trait", "anc", "n", "ref"]].to_dict("index").items()]
    phrases = {row["phrase"]: row["short"] or "" for row in cite.table(DATA / "abbreviations.tsv")}
    catalog = {row["id"]: row for row in cite.table(DATA / "catalog.tsv")} if (DATA / "catalog.tsv").exists() else {}
    new = draft(studies, {study: by[pmids[0]] if pmids else None for study, pmids in zip(d.index, cited)}, phrases, catalog)  # all drafted together, so that they differ

    path = DATA / "names.tsv"
    rows = cite.table(path) if path.exists() else []
    kept = {row["id"] for row in rows}
    rows += [{"id": study, "name": name, "how": how} for study, (name, how) in new.items() if study not in kept]
    order = {study: k for k, study in enumerate(d.index)}
    rows.sort(key=lambda row: order.get(row["id"], len(order)))
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, ["id", "name", "how"], delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None)
        writer.writeheader()
        writer.writerows(rows)
    twice = [name for name, n in Counter(row["name"].lower() for row in rows).items() if n > 1]
    odd = [row["name"] for row in rows if not merge.NAME.fullmatch(row["name"])]
    print(f"{len(rows)} names, {len(rows) - len(kept)} drafted now; {sum('check' in row['how'] for row in rows)} to check; "
          f"{len(twice)} twice: {twice[:8]}; {len(odd)} of another form: {odd[:8]}")


if __name__ == "__main__":
    main()
