"""The GWAS Catalog's tables of studies and ancestries, and the rule that picks candidates from them."""
import datetime
import re
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

DOWNLOADS = {"studies.tsv": "https://www.ebi.ac.uk/gwas/api/search/downloads/studies/v1.0.3.1",
             "ancestry.tsv": "https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry"}
STEPS = ["with summary statistics", "in a paper of at most {max_studies} such studies", "of European or East Asian ancestry only",
         "of at least {min_n:,} individuals", "in a journal that StatGenPapers reads"]  # the rule, one condition at a time


def download(cache: Path, refresh: bool = False) -> str:
    """Fetches the Catalog's tables into `cache`, unless they are there; returns the day they were fetched."""
    cache.mkdir(parents=True, exist_ok=True)
    for name, url in DOWNLOADS.items():
        if refresh or not (cache / name).exists():
            request = urllib.request.Request(url, headers={"User-Agent": "gwassumstats"})
            with urllib.request.urlopen(request, timeout=600) as response, (cache / name).open("wb") as file:
                while chunk := response.read(1 << 20):
                    file.write(chunk)
    return f"{datetime.date.fromtimestamp(min((cache / name).stat().st_mtime for name in DOWNLOADS)):%Y-%m-%d}"


def load(cache: Path) -> pd.DataFrame:
    """One row per study of the Catalog, by accession: its paper, its trait, whether its summary statistics are
    public, the trait and EFO term it is mapped to, and from the first stage of its ancestry table the number of
    individuals (`n`), their ancestries (`ancestry`), and the number of each ("European 53293; East Asian 1937":
    `ancestries`)."""
    read = lambda name: pd.read_csv(cache / name, sep="\t", dtype=str, keep_default_na=False, quoting=3)
    studies, stage = read("studies.tsv"), read("ancestry.tsv")[lambda t: t["STAGE"] == "initial"]
    groups = stage.assign(group=stage["BROAD ANCESTRAL CATEGORY"].str.split(r",\s*")).explode("group")[lambda t: t["group"] != "NR"]  # NR: not reported
    table = pd.DataFrame({
        "id": studies["STUDY ACCESSION"], "pmid": studies["PUBMED ID"], "author": studies["FIRST AUTHOR"], "published": studies["DATE"],
        "added": studies["DATE ADDED TO CATALOG"], "journal": studies["JOURNAL"], "title": studies["STUDY"], "trait": studies["DISEASE/TRAIT"],
        "sample": studies["INITIAL SAMPLE SIZE"], "cohort": studies["COHORT"].str.replace(r"\s*\|\s*", ", ", regex=True),
        "sumstats": studies["FULL SUMMARY STATISTICS"] == "yes", "location": studies["SUMMARY STATS LOCATION"],
        "mapped": studies["MAPPED_TRAIT"], "efo": studies["MAPPED_TRAIT_URI"].str.replace(r"[^,\s]*/", "", regex=True)}).set_index("id")
    counted = stage.assign(n=pd.to_numeric(stage["NUMBER OF INDIVDUALS"], errors="coerce").fillna(0).astype(int)).groupby(["STUDY ACCESSION", "BROAD ANCESTRAL CATEGORY"], sort=False)["n"].sum().reset_index()
    table["n"] = counted.groupby("STUDY ACCESSION")["n"].sum().reindex(table.index).fillna(0).astype(int)
    table["ancestry"] = groups.groupby("STUDY ACCESSION")["group"].agg(lambda g: ", ".join(sorted(set(g)))).reindex(table.index).fillna("")
    table["ancestries"] = (counted["BROAD ANCESTRAL CATEGORY"] + " " + counted["n"].astype(str)).groupby(counted["STUDY ACCESSION"]).agg("; ".join).reindex(table.index).fillna("")
    return table


def fetch(url: str, tries: int = 5) -> str:
    """The text at `url`; a server that refuses is asked again, a little later each time."""
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "gwassumstats"}), timeout=120) as response:
                return response.read().decode("utf-8", "replace")
        except (urllib.error.URLError, TimeoutError):
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


def sexes(studies: pd.DataFrame, cache: Path) -> pd.Series:
    """Whose statistics each study holds, as its submitters told the Catalog: "F", "M", "combined", or "" if they did
    not say. It is in the metadata file beside the study's summary statistics, fetched once into `cache`."""
    def one(study: str, location: str) -> str:
        kept = cache / f"{study}.yaml"
        if not kept.exists():
            folder = location.replace("http://", "https://").rstrip("/") + "/"
            names = re.findall(r'href="([^"/]+-meta\.yaml)"', fetch(folder))
            kept.write_text(fetch(folder + names[0]) if names else "", encoding="utf-8")
        found = re.search(r"^sex:\s*(\w+)", kept.read_text(encoding="utf-8"), re.M)
        return found[1] if found and found[1] != "NR" else ""

    cache.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(2) as pool:  # the Catalog's server refuses many requests at once
        return pd.Series(list(pool.map(one, studies.index, studies["location"])), index=studies.index)


def conditions(studies: pd.DataFrame, rule: dict) -> pd.DataFrame:
    """Whether each study meets each condition of the rule: a column per step of STEPS, in their order."""
    per_paper = studies.groupby("pmid")["sumstats"].transform("sum")
    return pd.DataFrame({0: studies["sumstats"], 1: per_paper <= rule["max_studies"], 2: studies["ancestry"].isin(rule["ancestries"]),
                         3: studies["n"] >= rule["min_n"], 4: studies["journal"].isin(rule["journals"])})


def papers(studies: pd.DataFrame) -> list[dict]:
    """Studies grouped by paper, the latest paper first: who, where, when, and the studies with trait, size and ancestry."""
    out = []
    for pmid, rows in studies.groupby("pmid", sort=False):
        first = rows.iloc[0]
        out.append({"pmid": pmid, "author": first["author"], "journal": first["journal"], "published": first["published"], "added": rows["added"].max(),
                    "title": first["title"], "studies": [{"id": i, "trait": r["trait"], "n": int(r["n"]), "ancestry": r["ancestry"]} for i, r in rows.iterrows()]})
    return sorted(out, key=lambda p: (p["published"], p["pmid"]), reverse=True)
