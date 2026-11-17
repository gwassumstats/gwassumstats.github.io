"""Drafts the complex trait of every study that has none: a name, an EFO term where one fits, and a broad category.

  uv run python -m traits [RELEASE]

RELEASE defaults to $SUMSTATS_RELEASE. Studies of one trait share it whatever their paper, ancestry, sex or
adjustment: the trait part of a study's short name, without those qualifiers (EAS, Female, AdjBMI, MTAG...), is the
key that groups them. A trait is named by hand (data/trait_aliases.tsv, which also joins keys that are one trait), or
else after the most frequent spelling of it among its studies' names in the release. Its EFO term is the term of
that name in the GWAS Catalog's table of trait mappings (downloaded to .cache/), or else the term the Catalog maps
its studies to, or else the one term the Catalog maps a trait so worded to. Its category is the first rule of data/categories.tsv that its name, or the Catalog's parent term of
its EFO term, matches.

The result is data/traits.tsv: study ID, trait, EFO term, category, and how they were found. A row that says "check"
was a guess. Correct a row by hand and write "by hand" in `how`: a later run keeps the rows that are there.
"""
import argparse
import csv
import os
import re
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from build import merge, release as rel
from papers import cite

DATA = merge.DATA
MAPPINGS = "https://www.ebi.ac.uk/gwas/api/search/downloads/trait_mappings"
QUALIFIER = re.compile(r"(EAS|Female|Male|MTAG|AdjBMI|UnadjBMI|AdjMed|Combined|NoUKB|UKB|N\d+k)$")  # what a short name adds to its trait


def key(name: str) -> str:
    """The trait of a short name, without what tells its study from others of the trait: Kanai_2018_SBPEAS -> SBP."""
    trait = name.split("_", 2)[2]
    while QUALIFIER.search(trait) and QUALIFIER.sub("", trait):
        trait = QUALIFIER.sub("", trait)
    return trait


TYPOS = {"cigarattes": "cigarettes", "duation": "duration", "epillepsy": "epilepsy", "hyperthension": "hypertension", "gastroesophagal": "gastroesophageal", "brainste": "brainstem",
         "even": "ever", "fibrilation": "fibrillation", "apen": "open", "alzheimers": "Alzheimer's", "alzheimer's": "Alzheimer's", "crohn's": "Crohn's", "raynaud's": "Raynaud's"}  # as the release spells them
DROPPED = re.compile(r"\s*\([^)]*\)|\b(levels?|direct|goldmann correlated|(combined )?automated (manual )?reading|adjusted by medication|of the human)\b", re.I)


def spelled(names: list[str]) -> str:
    """One name for a trait, from the names its studies have in the release: the most frequent spelling, then the
    longest, without what is in brackets and without the words of a measurement's variant, as a sentence."""
    def clean(name: str) -> str:
        words = re.sub(r"(?<=[a-z])(?=[A-Z][a-z])", " ", DROPPED.sub(" ", name.replace("’", "'"))).split()
        words = [TYPOS.get(word.lower(), word if re.search(r"[A-Z0-9]", word[1:]) else word.lower()) for word in words]  # an abbreviation keeps its capitals
        return (" ".join(words)[:1].upper() + " ".join(words)[1:]).strip()

    counted = Counter(clean(name) for name in names)
    return max(counted, key=lambda name: (counted[name], len(name), name))


def category(parent: str, name: str, rules: list[tuple[str, str]]) -> str:
    """A trait's broad category: the first rule whose pattern its name or its EFO parent term contains."""
    return next((broad for pattern, broad in rules if re.search(pattern, f"{name} | {parent}", re.I)), "other")


def draft(studies: pd.DataFrame, catalog: pd.DataFrame, mappings: pd.DataFrame, rules: list[tuple[str, str]], aliases: dict[str, tuple[str, str]]) -> pd.DataFrame:
    """One row per study: trait, efo, category, how. `studies`: name and trait (its name in the release) by ID;
    `catalog`: the GWAS Catalog's `efo` of a study it knows; `mappings`: the Catalog's table of trait mappings;
    `aliases`: {key: (trait, EFO label)} for the traits named by hand."""
    keys = studies["name"].map(key)
    named = {k: aliases[k][0] if k in aliases else spelled(list(group)) for k, group in studies["trait"].groupby(keys)}
    terms = mappings.assign(efo=mappings["EFO URI"].str.rsplit("/", n=1).str[1]).drop_duplicates("efo").set_index("efo")  # EFO ID -> its label and parent
    by_label = {label.lower(): efo for efo, label in terms["EFO term"].items()}
    one = mappings.assign(text=mappings["Disease trait"].str.lower(), efo=mappings["EFO URI"].str.rsplit("/", n=1).str[1]).groupby("text")["efo"].agg(lambda efo: efo.iloc[0] if efo.nunique() == 1 else "")
    by_trait = one[one != ""].to_dict()  # a trait as a paper words it -> its EFO term, where the Catalog maps those words to one term only
    trait = keys.map(named)
    rows = {}
    for name, group in studies.groupby(trait, sort=False):
        given = Counter(efo for efo in catalog["efo"].reindex(group.index).dropna() if efo and "," not in efo)
        labels = [label.lower() for k in set(keys[group.index]) if (label := aliases.get(k, ("", ""))[1])] + [name.lower(), f"{name.lower()} measurement"]
        looked = next((by_label[label] for label in labels if label in by_label), "")
        worded = next((by_trait[text] for text in [name.lower(), *group["trait"].str.lower()] if text in by_trait), "")
        efo, how = ((looked, "the EFO term of this name in the GWAS Catalog") if looked else (given.most_common(1)[0][0], "the GWAS Catalog's term for its studies") if given
                    else (worded, "the EFO term that the GWAS Catalog maps a trait of this name to") if worded else ("", "check: no EFO term found"))
        broad = category(terms.at[efo, "Parent term"] if efo in terms.index else "", name, rules)
        for study in group.index:
            rows[study] = {"trait": name, "efo": efo, "category": broad, "how": how + ("; check: no category" if broad == "other" else "")}
    return pd.DataFrame.from_dict(rows, orient="index").reindex(studies.index)


def main():
    parser = argparse.ArgumentParser(prog="python -m traits", description="Drafts the complex trait of every study of a release that has none.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--cache", type=Path, default=Path(__file__).parents[1] / ".cache" / "gwas-catalog", help="where the Catalog's tables are kept")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    d = rel.studies(Path(args.release))
    studies = d[["trait"]].assign(name=merge.names(d, DATA))
    reviewed = {row["id"]: row["label"] for row in cite.table(DATA / "labels.tsv")}  # a label of the release that is wrong: the trait is the file's
    studies["trait"] = [reviewed.get(study, trait) for study, trait in studies["trait"].items()]
    kept = args.cache / "trait_mappings.tsv"
    if not kept.exists():
        kept.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(urllib.request.Request(MAPPINGS, headers={"User-Agent": "gwassumstats"}), timeout=600) as response:
            kept.write_bytes(response.read())
    mappings = pd.read_csv(kept, sep="\t", dtype=str, keep_default_na=False)
    catalog = pd.read_csv(DATA / "catalog.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("id")
    rules = [(row["pattern"], row["category"]) for row in cite.table(DATA / "categories.tsv")]
    aliases = {row["key"]: (row["trait"], row["efo_label"] or "") for row in cite.table(DATA / "trait_aliases.tsv")}
    new = draft(studies, catalog, mappings, rules, aliases)

    path = DATA / "traits.tsv"
    rows = cite.table(path) if path.exists() else []
    known = {row["id"] for row in rows}
    rows += [{"id": study, **row} for study, row in new.to_dict("index").items() if study not in known]
    order = {study: k for k, study in enumerate(d.index)}
    rows.sort(key=lambda row: order.get(row["id"], len(order)))
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, ["id", "trait", "efo", "category", "how"], delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, quotechar=None)
        writer.writeheader()
        writer.writerows(rows)
    traits = defaultdict(set)
    for row in rows:
        traits[row["trait"]].add((row["efo"], row["category"]))
    torn = sorted(trait for trait, said in traits.items() if len(said) > 1)
    print(f"{len(rows)} studies, {len(rows) - len(known)} drafted now; {len(traits)} traits; {sum('check' in row['how'] for row in rows)} rows to check; "
          f"{len(torn)} traits with two EFO terms or categories: {torn[:6]}")


if __name__ == "__main__":
    main()
