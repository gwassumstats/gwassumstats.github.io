"""The next version of the collection, in next/ of the folder of the corrected copies: every study's file as it now
stands, the corrected copy where there is one and else the release's, with an inventory.

  SUMSTATS_CORRECTIONS=<the folder of the corrected copies> uv run python -m corrections.next [RELEASE]

RELEASE defaults to $SUMSTATS_RELEASE and is only read. Run it after `python -m build`: the inventory is the site's
(site/data/studies.json), with the study's file. It writes, under <folder>/next/ (the folder may be purged: a run
writes it all again):

  sumstats/<ID>.sumstats   a link to the study's file: its corrected copy, or the release's file; no file is copied
  inventory.tsv            a line a study: short name, ID, trait, ancestry, N, effective N, h² and on the liability
                           scale, QC, whether it is recommended or superseded, whether it includes UK Biobank, what was
                           corrected, and the paths of its link and of the file it leads to
  README.md                what the folder is
"""
import argparse
import json
import os
from pathlib import Path

import pandas as pd

from build import overlay, release as rel

ROOT = Path(__file__).parents[1]
COLUMNS = ["short_name", "id", "trait", "release_name", "ancestry", "n", "cases", "controls", "effective_n", "h2", "h2_se", "h2_z", "h2_liability", "h2_liability_se",
           "qc", "recommended", "superseded_by", "uk_biobank", "corrected", "reference", "pmid", "path", "file"]
README = """# The next version of the collection

Written by `python -m corrections.next` of the GWAS Sumstats repository, from the release and the corrected copies
beside this folder. Nothing here is copied: `sumstats/` holds a link per study to its file, the corrected copy where
there is one, else the release's file. `inventory.tsv` lists every study with what the site shows of it. The folder
may be purged: the command writes it again.
"""


def inventory(studies: list[dict], files: dict[str, str], out: Path) -> pd.DataFrame:
    """A row a study, from the site's records of the studies and the file each now stands in (`files`)."""
    names = {d["id"]: d["name"] for d in studies}
    row = lambda d: {"short_name": d["name"], "id": d["id"], "trait": d["of"], "release_name": d["trait"], "ancestry": d["anc"], "n": d["n"],
                     "cases": d.get("cases"), "controls": d.get("controls"), "effective_n": d.get("neff") if d.get("kind") == "binary" else d["n"],
                     "h2": d["h2"], "h2_se": d["h2_se"], "h2_z": d["h2_z"], "h2_liability": d.get("h2l"), "h2_liability_se": d.get("h2l_se"), "qc": d["qc"],
                     "recommended": "yes" if d.get("rec") else "", "superseded_by": names.get(d.get("sup"), ""),
                     "uk_biobank": {True: "yes", False: "no"}.get(d.get("ukb"), ""),
                     "corrected": "; ".join(kind for kind, _, _ in d.get("fixed", [])), "reference": d.get("label", d.get("ref", "")), "pmid": ";".join(d.get("cite", [])),
                     "path": str(out / "sumstats" / f"{d['id']}.sumstats"), "file": files[d["id"]]}
    table = pd.DataFrame([row(d) for d in studies], columns=COLUMNS)
    return table.astype({column: "Int64" for column in ("n", "cases", "controls", "effective_n")})


def linked(files: dict[str, str], folder: Path):
    """A link `folder`/<ID>.sumstats to each study's file; a link of a study that is gone is removed."""
    folder.mkdir(parents=True, exist_ok=True)
    for old in folder.iterdir():
        if old.is_symlink() and old.stem not in files:
            old.unlink()
    for study, file in files.items():
        link = folder / f"{study}.sumstats"
        if link.is_symlink() and os.readlink(link) == file:
            continue
        if link.is_symlink():
            link.unlink()
        link.symlink_to(file)


def write(release: Path, copies: Path, site: Path = ROOT / "site" / "data" / "studies.json", data: Path = ROOT / "data") -> Path:
    """Writes next/ under `copies`; returns its folder."""
    built = json.loads(site.read_text(encoding="utf-8"))
    theirs = rel.studies(release)
    d, *_ = overlay.laid(theirs, {anc: rel.pairs_of(release, anc) for anc in theirs["anc"].unique()}, data, copies)
    studies = built["studies"]
    if missing := [s["id"] for s in studies if s["id"] not in d.index] or [study for study in d.index if study not in {s["id"] for s in studies}]:
        raise ValueError(f"the site's studies are not the release's ({missing[0]}): run python -m build first")
    for s in studies:  # the label the inventory cites, as the site shows it
        s["label"] = built["papers"][s["cite"][0]]["label"] if s.get("cite") else s["ref"]
    files = d["h2_file"].to_dict()
    if absent := [study for study, file in files.items() if not Path(file).is_file()]:
        raise ValueError(f"the file of {absent[0]} is not there: {files[absent[0]]}")
    out = copies / "next"
    linked(files, out / "sumstats")
    table = inventory(studies, files, out)
    table.to_csv(out / "inventory.tsv", sep="\t", index=False)
    (out / "README.md").write_text(README, encoding="utf-8")
    return out


def main():
    parser = argparse.ArgumentParser(prog="python -m corrections.next", description="Writes the next version of the collection: a link to each study's file as it stands, and an inventory.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder, read only (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--out", default=os.environ.get("SUMSTATS_CORRECTIONS"), help="the folder of the corrected copies (default: $SUMSTATS_CORRECTIONS)")
    args = parser.parse_args()
    if not args.release or not args.out:
        parser.error("name the release folder and the folder of the corrected copies, or set SUMSTATS_RELEASE and SUMSTATS_CORRECTIONS")
    out = write(Path(args.release).resolve(), Path(args.out).resolve())
    table = pd.read_csv(out / "inventory.tsv", sep="\t", dtype=str, keep_default_na=False)
    print(f"next/: {len(table)} studies, {table['file'].str.startswith(str(out.parent / 'sumstats') + '/').sum()} by their corrected copy; "
          f"{(table['recommended'] == 'yes').sum()} recommended; {(table['uk_biobank'] == 'yes').sum()} with UK Biobank")


if __name__ == "__main__":
    main()
