"""Writes a corrected copy of every study that data/corrections.tsv says to take from its source again, to flip, to
clean or to give another N.

  SUMSTATS_CORRECTIONS=<a folder outside the repository> uv run --group qc python -m corrections [RELEASE] \\
      [--snplist FILE --bim-prefix PREFIX] [--account NAME]

RELEASE defaults to $SUMSTATS_RELEASE and is only read: its file index says where each study's file is. A study with
a `source` is downloaded and formatted as the release formats (sources/), which needs FILE, the HapMap3 SNPs with
their alleles, and PREFIX, the reference panel's .bim files: a Slurm task each. The other studies share one task.
Under the folder named by $SUMSTATS_CORRECTIONS the command writes:

  sumstats/<ID>.sumstats   the corrected copy, in the release's format: nothing but what a fix changes is rewritten
  provenance/<ID>.json     where it came from, the checksums of source and copy, the rows of each, and what changed
  sources/<SOURCE>/        the files as downloaded, kept for a later run
  corrections.tsv          a line a copy: study, fixes, rows, checksum

The folder may be purged: a run writes everything again, from the release's files and from the sources. Nothing here
enters the repository, and nothing is written to the release or to any other folder.
"""
import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

from corrections import fix
from papers import cite
from sources.make import made
from sources.slurm import submitted

ROOT = Path(__file__).parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        while chunk := file.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def planned(table: list[dict]) -> dict[str, list[dict]]:
    """{study: its fixes that change the file, in the table's order}; a relabel changes no file. A study's source,
    if it has one, is where its other fixes start from: it comes first."""
    plan = {}
    for row in table:
        if row["fix"] in (*fix.FIXES, "source"):
            plan.setdefault(row["id"], []).append(row)
        elif row["fix"] != "relabel":
            raise ValueError(f"data/corrections.tsv: {row['fix']} is no fix; expected one of {', '.join(['source', *fix.FIXES, 'relabel'])}")
    if late := [study for study, fixes in plan.items() if any(row["fix"] == "source" for row in fixes[1:])]:
        raise ValueError(f"data/corrections.tsv: the source of {late[0]} comes after another fix of it, or twice")
    return plan


def correct(plan: dict[str, list[dict]], sources: dict[str, str], out: Path, code: str, snplist: str = "", bim_prefix: str = "", sizes: dict[str, float] | None = None) -> list[dict]:
    """One Slurm task: writes each study's corrected copy and its record under `out`; returns a line a study.
    `sources`: the release's file of each study; `snplist`, `bim_prefix`: the reference that a study taken from its
    source is formatted with; `sizes`: {accession: N}, the GWAS Catalog's, for a source's file that has no N."""
    lines = []
    for study, fixes in plan.items():
        source, copy, record = Path(sources[study]), out / "sumstats" / f"{study}.sumstats", out / "provenance" / f"{study}.json"
        rows = before = fix.read(source)
        changed, taken = {}, {}
        copy.parent.mkdir(parents=True, exist_ok=True)
        record.parent.mkdir(parents=True, exist_ok=True)
        if fixes[0]["fix"] == "source":  # the study's statistics from their source, in place of the release's file
            said = fix.sourced(fixes[0]["value"])
            taken = made(said["source"], copy, snplist, bim_prefix, said["n"] if said["n"] is not None else (sizes or {}).get(said["source"]), said["file"], said["roles"], own_n=said["n"] is None)
            rows = fix.read(copy)
            changed = {"rows_of_the_source": taken["formatting"]["rows"], "snps_also_in_the_release_file": int(rows["SNP"].isin(before["SNP"]).sum())}
        for row in fixes[bool(taken):]:
            rows, what = fix.FIXES[row["fix"]](rows, row["value"])
            changed |= what
        rows.to_csv(copy, sep="\t", index=False, lineterminator="\n")
        if list(fix.read(copy).columns) != list(before.columns) or not fix.finite(rows["Z"]).all() or rows["SNP"].duplicated().any() and any(r["fix"] == "clean" for r in fixes):
            raise ValueError(f"{study}: the corrected copy is not what was meant")
        line = {"id": study, "fixes": "; ".join(row["fix"] + (f" {row['value']}" if row["value"] else "") for row in fixes), "rows_before": len(before), "rows": len(rows), "sha256": sha256(copy)}
        theirs = {"path": str(source), "sha256": sha256(source), "rows": len(before)}  # the release's file: the source of a copy, unless the copy has a source of its own
        record.write_text(json.dumps({"study": study, "fixes": [{"fix": row["fix"], "value": row["value"], "why": row["why"]} for row in fixes], "changed": changed,
                                      **({"source": taken["source"], "formatting": taken["formatting"], "in_place_of": theirs} if taken else {"source": theirs}),
                                      "copy": {"path": str(copy), "sha256": line["sha256"], "rows": len(rows)}, "made": f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}", "code": code}, indent=1) + "\n")
        lines.append(line | {key: value for key, value in changed.items()})
        print(line, changed, flush=True)
    return lines


def main():
    parser = argparse.ArgumentParser(prog="python -m corrections", description="Writes corrected copies of the studies of data/corrections.tsv, outside the repository.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder, read only (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--out", default=os.environ.get("SUMSTATS_CORRECTIONS"), help="where to write (default: $SUMSTATS_CORRECTIONS)")
    parser.add_argument("--snplist", default="", help="the HapMap3 SNPs and their alleles (SNP, A1, A2): for the studies taken from their source")
    parser.add_argument("--bim-prefix", default="", help="the reference panel's .bim files, up to the chromosome's number: for the same")
    parser.add_argument("--account", help="the Slurm account to charge")
    parser.add_argument("--name", default=f"{datetime.date.today():%Y%m%d}_corrections", help="the name of the Slurm tasks and of their logs")
    args = parser.parse_args()
    if not args.release or not args.out:
        parser.error("name the release folder and the folder to write to, or set SUMSTATS_RELEASE and SUMSTATS_CORRECTIONS")
    release, out = Path(args.release).resolve(), Path(args.out).resolve()
    if out.is_relative_to(ROOT) or out.is_relative_to(release) or release.is_relative_to(out):
        parser.error("the corrected copies do not go inside the repository, nor inside the release")
    index = pd.read_csv(release / "_audit" / "sumstats_file_index.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("Trait Identifier")["H2 Sumstats"]
    plan = planned(cite.table(ROOT / "data" / "corrections.tsv"))
    if unknown := sorted(set(plan) - set(index.index)):
        parser.error(f"data/corrections.tsv names a study that the release lacks: {unknown[0]}")
    from_source = [study for study, fixes in plan.items() if fixes[0]["fix"] == "source"]
    if from_source and not (args.snplist and args.bim_prefix):
        parser.error(f"{len(from_source)} studies are taken from their source: name the reference they are formatted with, --snplist and --bim-prefix")
    sizes = {row["accession"]: float(row["n"]) for row in cite.table(ROOT / "data" / "samples.tsv")}
    code = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    parts = [{study: plan[study] for study in plan if study not in from_source}, *({study: plan[study]} for study in from_source)]  # a task for each study from its source, one for the rest
    parts = [part for part in parts if part]
    k = len(parts)
    lines = [line for task in submitted(args.name, out, args.account, correct, parts, [{study: index[study] for study in part} for part in parts], [out] * k, [code] * k,
                                        [args.snplist] * k, [args.bim_prefix] * k, [sizes] * k) for line in task]
    pd.DataFrame(lines).to_csv(out / "corrections.tsv", sep="\t", index=False)
    missing = [study for study in plan if not (out / "sumstats" / f"{study}.sumstats").exists() or not (out / "provenance" / f"{study}.json").exists()]
    print(f"{len(lines)} copies written, {len(missing)} missing")
    if missing:
        sys.exit(f"missing: {', '.join(missing)}")


if __name__ == "__main__":
    main()
