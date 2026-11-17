"""A source's file, downloaded and written in the release's format."""
from pathlib import Path

from sources import fetch
from sources.format import formatted, reference


def made(source: str, to: Path, snplist: str, bim_prefix: str, n: float | None = None, file: str = "", roles: dict[str, str] | None = None, own_n: bool = True) -> dict:
    """Downloads a source's file beside `to`'s folder (sources/<source>/), writes it in the release's format to `to`,
    and returns its record: the source (URL, size, checksums) and what the formatting did. `source`, `file`: see
    fetch.located; `snplist`, `bim_prefix`, `n`, `roles`, `own_n`: see format.reference and format.formatted."""
    url, md5 = fetch.located(source, file)
    kept = to.parents[1] / "sources" / (source if url != source else "by-url") / url.rsplit("/", 1)[1]
    got = fetch.download(url, kept, md5)
    table, done = formatted(kept, *reference(snplist, bim_prefix), n=n, roles=roles, own_n=own_n)
    to.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(to, sep="\t", index=False, lineterminator="\n")
    return {"source": got, "formatting": done}
