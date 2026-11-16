"""PubMed, through NCBI's E-utilities: a search, and the bibliographic record of a paper."""
import json
import time
import urllib.parse
import urllib.request

API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PAUSE = 0.4  # seconds between requests: NCBI allows three a second without a key


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
    PubMed may add to it), title, journal, year, volume, issue, pages, doi."""
    out = {}
    for start in range(0, len(pmids), 100):
        result = ask("esummary", id=",".join(pmids[start:start + 100]))["result"]
        for pmid in result["uids"]:
            r = result[pmid]
            doi = next((i["value"] for i in r.get("articleids", []) if i["idtype"] == "doi"), "")
            name = lambda a: "{" + a["name"].partition(". Electronic address")[0] + "}" if a.get("authtype") == "CollectiveName" else a["name"]
            out[pmid] = {"authors": [name(a) for a in r.get("authors", [])], "title": r.get("title", "").rstrip("."), "journal": r.get("source", ""),
                         "year": r.get("pubdate", "")[:4], "volume": r.get("volume", ""), "issue": r.get("issue", ""), "pages": r.get("pages", ""), "doi": doi}
    return out
