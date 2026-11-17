"""Builds the site's data from a release folder.

  SUMSTATS_PASSPHRASE=... uv run python -m build RELEASE

RELEASE defaults to $SUMSTATS_RELEASE; without $SUMSTATS_PASSPHRASE the build asks for the passphrase. Neither is
ever written to the repository. The build writes, in site/data/:

  studies.json         the inventory with h², each study's short name, the paper it cites and
                       what the quality checks say of it, and the release's counts       public
  papers.json          every paper cited: authors, title, journal, DOI, BibTeX entry     public
  candidates.json      the candidates for the next release, and the rule's backtest      public
  rg_{ancestry}.json   every genetic correlation of an ancestry                          public
  paths.enc.json       where each study's file is, encrypted (crypt.py)                  readable after login

The papers are those of data/references.tsv and data/papers.json (python -m papers); the short names, those of
data/names.tsv (python -m names). A study without a short name stops the build. The quality checks (qc/checks.py)
judge each study from its file as measured (data/metrics.tsv: python -m qc.scan), its results and its label. The
candidates are those of data/candidates.json (python -m candidates), as they are.

A release that has not changed rewrites the same files, and keeps the encrypted one as it is, so readers stay logged in.
"""
import argparse
import getpass
import gzip
import json
import os
from pathlib import Path

from cryptography.exceptions import InvalidTag

from build import crypt, merge, release as rel

SITE_DATA = Path(__file__).parents[1] / "site" / "data"
BLOB = "paths.enc.json"


def dumps(value) -> str:
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"))


def holds(path: Path, plain: bytes, passphrase: str, iterations: int) -> bool:
    """Whether the blob at `path` already encrypts `plain` under this passphrase."""
    try:
        blob = json.loads(path.read_text())
        return blob["iterations"] == iterations and crypt.decrypt(blob, passphrase) == plain
    except (FileNotFoundError, InvalidTag, KeyError, ValueError):
        return False


def build(release: Path, out: Path, passphrase: str, iterations: int = crypt.ITERATIONS, data: Path = merge.DATA) -> list[Path]:
    """Writes the site's data for the release at `release` into `out`, with what `data` adds to it; returns the files."""
    d = rel.studies(release)
    rg = {anc: rel.correlations(release, d, anc) for anc in d["anc"].unique()}
    cited, papers = merge.citations(d, data)
    judged = merge.quality(d, rg, data)
    public = d.assign(name=merge.names(d, data)).join(judged[["qc", "hm3", "chi2", "lambda"]])[["name", *rel.PUBLIC, "qc", "hm3", "chi2", "lambda"]].reset_index()
    records = public.astype(object).where(public.notna(), None).to_dict("records")
    records = [{k: v for k, v in r.items() if v is not None} | ({"cite": pmids} if pmids else {}) | ({"notes": notes} if notes else {})
               for r, pmids, notes in zip(records, cited, judged["notes"])]
    summary = rel.summary(release, d, rg) | {"qc": {level: int((judged["qc"] == level).sum()) for level in ("pass", "warn", "fail")}}
    lines = lambda items: ",\n".join(items)  # a study, a paper a line, for readable diffs
    cited_as = {pmid: {"label": p["label"], "published": p["published"]} for pmid, p in papers.items()}  # what the table of studies needs of a paper
    texts = {"studies.json": f'{{\n"release":{dumps(summary)},\n"papers":{dumps(cited_as)},\n"studies":[\n'
                             + lines(map(dumps, records)) + "\n]}\n",
             "papers.json": "{\n" + lines(f"{dumps(pmid)}:{dumps(paper)}" for pmid, paper in papers.items()) + "\n}\n",
             "candidates.json": dumps(json.loads((data / "candidates.json").read_text(encoding="utf-8"))) + "\n"}
    for anc, t in rg.items():
        texts[f"rg_{anc}.json"] = "{\n" + ",\n".join(f"{dumps(k)}:{dumps(v)}" for k, v in t.items()) + "\n}\n"
    out.mkdir(parents=True, exist_ok=True)
    for name, text in texts.items():
        (out / name).write_text(text, encoding="utf-8")
    plain = dumps(rel.private(release, d)).encode()
    if not holds(out / BLOB, plain, passphrase, iterations):
        (out / BLOB).write_text(json.dumps(crypt.encrypt(plain, passphrase, iterations), indent=1) + "\n")
    return [out / name for name in [*texts, BLOB]]


def main():
    parser = argparse.ArgumentParser(prog="python -m build", description="Builds the site's data from a release folder.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--out", type=Path, default=SITE_DATA, help="where to write (default: site/data)")
    parser.add_argument("--data", type=Path, default=merge.DATA, help="what the repository adds to the release (default: data)")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    passphrase = os.environ.get("SUMSTATS_PASSPHRASE") or getpass.getpass("Passphrase: ")
    if not passphrase:
        parser.error("the passphrase is empty")
    for path in build(Path(os.path.abspath(args.release)), args.out, passphrase, data=args.data):
        size = path.stat().st_size
        print(f"{path.name:<16} {size / 1e3:>8,.0f} kB  ({len(gzip.compress(path.read_bytes(), 6)) / 1e3:,.0f} kB gzipped)")


if __name__ == "__main__":
    main()
