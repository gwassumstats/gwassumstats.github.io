"""Writes a corrected copy of every study that data/corrections.tsv says to flip, to clean or to give another N.

  SUMSTATS_CORRECTIONS=<a folder outside the repository> uv run --group qc python -m corrections [RELEASE] [--account NAME]

RELEASE defaults to $SUMSTATS_RELEASE and is only read: its file index says where each study's file is. Under the
folder named by $SUMSTATS_CORRECTIONS the command writes, in one Slurm task:

  sumstats/<ID>.sumstats   the corrected copy, in the release's format: nothing but what the fix changes is rewritten
  provenance/<ID>.json     where it came from, the checksums of source and copy, the rows of each, and what changed
  corrections.tsv          a line a copy: study, fixes, rows, checksum

The folder may be purged: a run writes everything again from the release's files. Nothing here enters the repository,
and nothing is written to the release or to any other folder.
"""
import argparse
import datetime
import getpass
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

from corrections import fix
from papers import cite

ROOT = Path(__file__).parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as file:
        while chunk := file.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def planned(table: list[dict]) -> dict[str, list[dict]]:
    """{study: its fixes that change the file, in the table's order}; a relabel changes no file."""
    plan = {}
    for row in table:
        if row["fix"] in fix.FIXES:
            plan.setdefault(row["id"], []).append(row)
        elif row["fix"] != "relabel":
            raise ValueError(f"data/corrections.tsv: {row['fix']} is no fix; expected one of {', '.join([*fix.FIXES, 'relabel'])}")
    return plan


def correct(plan: dict[str, list[dict]], sources: dict[str, str], out: Path, code: str) -> list[dict]:
    """One Slurm task: writes each study's corrected copy and its record under `out`; returns a line a study."""
    lines = []
    for study, fixes in plan.items():
        source, copy, record = Path(sources[study]), out / "sumstats" / f"{study}.sumstats", out / "provenance" / f"{study}.json"
        rows = before = fix.read(source)
        changed = {}
        for row in fixes:
            rows, what = fix.FIXES[row["fix"]](rows, row["value"])
            changed |= what
        copy.parent.mkdir(parents=True, exist_ok=True)
        record.parent.mkdir(parents=True, exist_ok=True)
        rows.to_csv(copy, sep="\t", index=False, lineterminator="\n")
        if list(fix.read(copy).columns) != list(before.columns) or not fix.finite(rows["Z"]).all() or rows["SNP"].duplicated().any() and any(r["fix"] == "clean" for r in fixes):
            raise ValueError(f"{study}: the corrected copy is not what was meant")
        line = {"id": study, "fixes": "; ".join(row["fix"] + (f" {row['value']}" if row["value"] else "") for row in fixes), "rows_before": len(before), "rows": len(rows), "sha256": sha256(copy)}
        record.write_text(json.dumps({"study": study, "fixes": [{"fix": row["fix"], "value": row["value"], "why": row["why"]} for row in fixes], "changed": changed,
                                      "source": {"path": str(source), "sha256": sha256(source), "rows": len(before)}, "copy": {"path": str(copy), "sha256": line["sha256"], "rows": len(rows)},
                                      "made": f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}", "code": code}, indent=1) + "\n")
        lines.append(line | {key: value for key, value in changed.items()})
        print(line, changed, flush=True)
    return lines


def main():
    import submitit  # only this command needs it: the group "qc"

    parser = argparse.ArgumentParser(prog="python -m corrections", description="Writes corrected copies of the studies of data/corrections.tsv, outside the repository.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder, read only (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--out", default=os.environ.get("SUMSTATS_CORRECTIONS"), help="where to write (default: $SUMSTATS_CORRECTIONS)")
    parser.add_argument("--account", help="the Slurm account to charge")
    parser.add_argument("--name", default=f"{datetime.date.today():%Y%m%d}_corrections", help="the name of the Slurm task and of its logs")
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
    logs = out / "submitit-logs"
    folder, record = logs / args.name, logs / f"{args.name}.log"
    queued = subprocess.run(["squeue", "-h", "-u", getpass.getuser(), "-n", args.name], capture_output=True, text=True).stdout.strip()
    if folder.exists() or record.exists() or queued:
        parser.error(f"{args.name} has run already, or is running: choose another --name")
    logs.mkdir(parents=True, exist_ok=True)
    code = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    executor = submitit.AutoExecutor(folder=folder)
    executor.update_parameters(name=args.name, slurm_partition="short", cpus_per_task=1, mem_gb=8, timeout_min=60, slurm_srun_args=["--export=ALL"],
                               slurm_additional_parameters={"account": args.account} if args.account else {})
    jobs = executor.map_array(correct, [plan], [{study: index[study] for study in plan}], [out], [code])

    def note(text):
        with record.open("a") as log:
            log.write(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {text}\n")
        print(text, flush=True)

    note(f"submitted {jobs[0].job_id}: {len(plan)} studies; short, 1 CPU, 8 GB, 60 min; expected: sumstats/ and provenance/ with {len(plan)} files each, corrections.tsv")
    submitit.helpers.monitor_jobs(jobs)
    lines = [line for job in jobs for line in job.result()]  # raises what the task raised
    pd.DataFrame(lines).to_csv(out / "corrections.tsv", sep="\t", index=False)
    missing = [study for study in plan if not (out / "sumstats" / f"{study}.sumstats").exists() or not (out / "provenance" / f"{study}.json").exists()]
    note(f"finished: {len(lines)} copies written, {len(missing)} missing")
    if missing:
        sys.exit(f"missing: {', '.join(missing)}")


if __name__ == "__main__":
    main()
