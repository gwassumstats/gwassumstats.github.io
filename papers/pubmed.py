"""PubMed, through NCBI's E-utilities: a search, and the bibliographic record of a paper."""
import json
import re
import time
import urllib.parse
import urllib.request

API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PAUSE = 0.4  # seconds between requests: NCBI allows three a second without a key
MONTHS = {name: number for number, name in enumerate("Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split(), 1)}


def date(text: str) -> str:
    """A date as PubMed writes it, "2010 Nov 14", "2010 Nov", "2019 Jan-Feb" or "2010", as 2010-11-14, 2010-11, 2019-01
    or 2010: as precise as PubMed gives it, and no more."""
    found = re.match(r"(\d{4})(?: ([A-Z][a-z]{2})[a-z]*(?:-[A-Za-z]+)?(?: (\d{1,2}))?)?", text)
    if not found:
        return ""
    year, month, day = found.groups()
    if month not in MONTHS:  # a year alone, or a season
        return year
    return f"{year}-{MONTHS[month]:02d}" + (f"-{int(day):02d}" if day else "")


def earliest(*dates: str) -> str:
    """The earliest of dates of unequal precision; of two that may be the same day, the more precise: 2019-01 and
    2019-01-14 give 2019-01-14, 2019-01 and 2019-02-02 give 2019-01."""
    return min((d for d in dates if d), key=lambda d: d + "-99-99"[len(d) - 4:], default="")


def ask(service: str, **query) -> dict:
    time.sleep(PAUSE)
    url = f"{API}/{service}.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "tool": "gwassumstats", **query})
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.load(response)


def search(term: str, most: int = 20) -> list[str]:
    """The PMIDs of the papers that match a PubMed query, the best match first."""
    return ask("esearch", term=term, retmax=most, sort="relevance")["esearchresult"].get("idlist", [])


def records(pmids: list[str]) -> dict[str, dict]:
    """{PMID: the paper}: authors (as PubMed spells them: "Speliotes EK"; a consortium in braces, without the address
    PubMed may add to it), title, journal, year (of the issue), volume, issue, pages, doi, and `published`: the day
    the paper first appeared, online or in print, as precisely as PubMed dates it."""
    out = {}
    for start in range(0, len(pmids), 100):
        result = ask("esummary", id=",".join(pmids[start:start + 100]))["result"]
        for pmid in result["uids"]:
            r = result[pmid]
            doi = next((i["value"] for i in r.get("articleids", []) if i["idtype"] == "doi"), "")
            name = lambda a: "{" + a["name"].partition(". Electronic address")[0] + "}" if a.get("authtype") == "CollectiveName" else a["name"]
            out[pmid] = {"authors": [name(a) for a in r.get("authors", [])], "title": r.get("title", "").rstrip("."), "journal": r.get("source", ""),
                         "year": r.get("pubdate", "")[:4], "volume": r.get("volume", ""), "issue": r.get("issue", ""), "pages": r.get("pages", ""), "doi": doi,
                         "published": earliest(date(r.get("pubdate", "")), date(r.get("epubdate", "")))}
    return out
