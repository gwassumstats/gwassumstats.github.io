"""Where the summary statistics of a study are, and their download: a study's folder in the GWAS Catalog, or a URL."""
import hashlib
import http.client
import re
import time
import urllib.error
import urllib.request
from pathlib import Path

CATALOG = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
AGENT = {"User-Agent": "gwassumstats"}
STATISTICS = (".tsv", ".txt", ".csv", ".tbl", ".out", ".gz", ".bgz", ".zip")  # how a file of statistics may end; the rest of a folder is about it


def folder(accession: str) -> str:
    """The folder of a study in the GWAS Catalog, which keeps its studies by the thousand."""
    number, width = int(accession[4:]), len(accession) - 4
    first = (number - 1) // 1000 * 1000 + 1
    return f"{CATALOG}GCST{first:0{width}d}-GCST{first + 999:0{width}d}/{accession}/"


def read(url: str, tries: int = 5) -> bytes:
    """What is at `url`; a server that does not answer is asked again, a little later each time."""
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=AGENT), timeout=120) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError):
            if attempt == tries - 1:
                raise
            time.sleep(5 * (attempt + 1))


def chosen(names: list[str], file: str = "") -> str:
    """The file of statistics among the names of a folder: the one named, or the only one there is."""
    files = [name for name in names if name.lower().endswith(STATISTICS) and not name.lower().endswith(("md5sum.txt", "readme.txt"))]
    if file:
        if file not in names:
            raise ValueError(f"no {file} in the folder, which has {', '.join(names)}")
        return file
    if len(files) != 1:
        raise ValueError(f"the folder has {len(files)} files of statistics ({', '.join(files) or 'none'}): name one")
    return files[0]


def located(source: str, file: str = "") -> tuple[str, str]:
    """The URL of a source's file, and the MD5 that its folder states for it ("" if none). `source`: an accession of
    the GWAS Catalog, or a URL, taken as it is. `file`: for an accession, the file's name in the study's folder, with
    the folder inside it if it is in one ("harmonised/x.tsv.gz"; "harmonised/" for that folder's only file); nothing
    for the study's only file."""
    if not re.fullmatch(r"GCST\d+", source):
        return source, ""
    inside, _, named = file.rpartition("/")
    base = folder(source) + (f"{inside}/" if inside else "")
    names = re.findall(r'href="([^"?/][^"/]*)"', read(base).decode())
    name = chosen(names, named)
    stated = dict(reversed(line.split()) for line in read(base + "md5sum.txt").decode().splitlines() if len(line.split()) == 2) if "md5sum.txt" in names else {}
    return base + name, stated.get(name, "")


def download(url: str, to: Path, md5: str = "") -> dict:
    """Fetches `url` into the file `to`, unless it is there already, and says what it is: its URL, its size and its
    checksums. A file whose MD5 is not the one stated is an error, and is not kept. A response cut short by the
    server shows in that checksum, where the folder states one."""
    for attempt in range(0 if to.exists() else 3):  # a transfer that breaks off is begun again, twice
        to.parent.mkdir(parents=True, exist_ok=True)
        part = to.with_name(to.name + ".part")
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=AGENT), timeout=600) as response, part.open("wb") as file:
                while block := response.read(1 << 22):
                    file.write(block)
        except (urllib.error.URLError, TimeoutError, ConnectionError, http.client.HTTPException):
            if attempt == 2:
                raise
            time.sleep(30)
        else:
            part.replace(to)
            break
    sums = hashlib.sha256(), hashlib.md5()
    with to.open("rb") as file:
        while block := file.read(1 << 22):
            for digest in sums:
                digest.update(block)
    if md5 and sums[1].hexdigest() != md5:
        to.unlink()
        raise ValueError(f"{url}: MD5 {sums[1].hexdigest()}, where its folder states {md5}")
    return {"url": url, "bytes": to.stat().st_size, "sha256": sums[0].hexdigest(), "md5": sums[1].hexdigest(), "md5_stated": md5}
