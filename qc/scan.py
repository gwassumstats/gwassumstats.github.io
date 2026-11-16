"""Measures every sumstats file of a release, on Slurm: the per-SNP numbers behind the site's quality checks.

  uv run --group qc python -m qc.scan RELEASE --hm3 EUR='PATTERN' --hm3 EAS='PATTERN' --logs FOLDER [--account NAME]

RELEASE is the release folder (default: $SUMSTATS_RELEASE); its file index says where each study's files are. A
PATTERN names the LD score weight files of an ancestry, one per chromosome: their SNPs are the HapMap3 SNPs that LDSC
regresses on. FOLDER takes the Slurm logs, which hold paths: keep it outside the repository. The result,
data/metrics.tsv, holds none: a row per file, named by its study's ID.
"""
import argparse
import datetime
import getpass
import glob
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = Path(__file__).parents[1] / "data" / "metrics.tsv"
MEDIAN_CHI2 = 0.4549364  # of a chi-square with one degree of freedom: lambda GC is the median of Z² over it
TASKS = 12  # a file takes a second or two: many files a task


def measure(path: str, hm3: pd.Index) -> dict:
    """The numbers of one file: its columns and rows, what is wrong with its SNPs and alleles, its N and its Z."""
    table = pd.read_csv(path, sep="\t", dtype={"SNP": str, "A1": str, "A2": str})
    n, z = (pd.to_numeric(table[c], errors="coerce") if c in table else pd.Series(dtype=float) for c in ("N", "Z"))
    z_missing = int((~np.isfinite(z)).sum())
    z = z[np.isfinite(z)]  # an infinite Z is counted, and left out of the rest
    chi2 = z ** 2
    out = {"columns": " ".join(table.columns), "rows": len(table), "duplicates": int(table["SNP"].duplicated().sum()),
           "not_rs": int((~table["SNP"].str.fullmatch(r"rs\d+", na=False)).sum()), "in_hm3": int(hm3.isin(table["SNP"]).sum()), "hm3": len(hm3)}  # each SNP once
    if {"A1", "A2"} <= set(table.columns):
        a1, a2 = table["A1"].str.upper(), table["A2"].str.upper()
        out["bad_alleles"] = int((~a1.isin(list("ACGT")) | ~a2.isin(list("ACGT")) | (a1 == a2)).sum())
        out["ambiguous"] = int((a1 + a2).isin(["AT", "TA", "CG", "GC"]).sum())  # strand-ambiguous
    return out | {"n_missing": int((~(n > 0)).sum()), "n_min": n.min(), "n_median": n.median(), "n_max": n.max(),
                  "z_missing": z_missing, "z_mean": z.mean(), "z_sd": z.std(), "z_max": z.abs().max(),
                  "z_positive": (z > 0).mean(), "chi2_mean": chi2.mean(), "lambda_gc": chi2.median() / MEDIAN_CHI2}


def scan(files: list[tuple[str, str, str, str]], patterns: dict[str, str]) -> list[dict]:
    """One Slurm task: measures (study ID, "h2" or "rg", ancestry, path) one by one, against its ancestry's HapMap3 SNPs."""
    hm3 = {anc: pd.Index(pd.concat(pd.read_csv(path, sep="\t", usecols=["SNP"]) for path in sorted(glob.glob(pattern)))["SNP"])
           for anc, pattern in patterns.items() if anc in {f[2] for f in files}}
    rows = []
    for k, (study, kind, anc, path) in enumerate(files):
        start = time.perf_counter()
        rows.append({"id": study, "file": kind} | measure(path, hm3[anc]))
        print(f"[{k + 1}/{len(files)}] {study} {kind}: {time.perf_counter() - start:.1f} s", flush=True)
    return rows


def files_of(release: Path) -> list[tuple[str, str, str, str]]:
    """Every file of the release, once: a study's h² file and, if it is another file, its rg file."""
    index = pd.read_csv(release / "_audit" / "sumstats_file_index.tsv", sep="\t", dtype=str, keep_default_na=False)
    h2 = index.assign(kind="h2", path=index["H2 Sumstats"])
    rg = index.assign(kind="rg", path=index["RG Sumstats"])[lambda t: (t["path"] != "") & (t["path"] != t["H2 Sumstats"])]
    return list(pd.concat([h2, rg])[["Trait Identifier", "kind", "Population", "path"]].itertuples(index=False, name=None))


def main():
    import submitit  # only the scan needs it: the group "qc"

    parser = argparse.ArgumentParser(prog="python -m qc.scan", description="Measures every sumstats file of a release, on Slurm.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--hm3", action="append", required=True, metavar="ANCESTRY=PATTERN", help="the LD score weight files of an ancestry")
    parser.add_argument("--logs", type=Path, required=True, help="a folder outside the repository, for the Slurm logs")
    parser.add_argument("--account", help="the Slurm account to charge")
    parser.add_argument("--name", default=f"{datetime.date.today():%Y%m%d}_qc_scan", help="the name of the array and of its logs")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    patterns = dict(pattern.split("=", 1) for pattern in args.hm3)
    files = files_of(Path(args.release))
    missing = sorted({anc for _, _, anc, _ in files} - {anc for anc, pattern in patterns.items() if glob.glob(pattern)})
    if missing:
        parser.error(f"no --hm3 files for {', '.join(missing)}")
    folder, record = args.logs / args.name, args.logs / f"{args.name}.log"
    queued = subprocess.run(["squeue", "-h", "-u", getpass.getuser(), "-n", args.name], capture_output=True, text=True).stdout.strip()
    if folder.exists() or record.exists() or queued:
        parser.error(f"{args.name} has run already, or is running: choose another --name")
    args.logs.mkdir(parents=True, exist_ok=True)

    executor = submitit.AutoExecutor(folder=folder)
    executor.update_parameters(name=args.name, slurm_partition="short", cpus_per_task=1, mem_gb=8, timeout_min=60, slurm_srun_args=["--export=ALL"],
                               slurm_additional_parameters={"account": args.account} if args.account else {})
    chunks = [files[k::TASKS] for k in range(TASKS)]
    jobs = executor.map_array(scan, chunks, [patterns] * TASKS)

    def note(text):
        with record.open("a") as log:
            log.write(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {text}\n")
        print(text, flush=True)

    note(f"submitted {len(jobs)} tasks for {len(files)} files: jobs {', '.join(job.job_id for job in jobs)}; short, 1 CPU, 8 GB, 60 min; result: {METRICS}")
    submitit.helpers.monitor_jobs(jobs)
    rows, failures = [], []
    for chunk, job in zip(chunks, jobs, strict=True):
        try:
            rows += job.result()
        except Exception as error:  # every task is asked, so that every failure is named
            failures.append(f"{job.job_id} ({len(chunk)} files from {chunk[0][0]}): {error}")
    note(f"finished: {len(jobs) - len(failures)} tasks completed, {len(failures)} failed, {len(rows)} files measured")
    if failures:
        raise SystemExit("\n".join(failures))
    metrics = pd.DataFrame(rows).set_index(["id", "file"]).loc[[(study, kind) for study, kind, _, _ in files]]  # in the release's order; a missing file is an error
    metrics.round(4).to_csv(METRICS, sep="\t")
    note(f"wrote {len(metrics)} rows to {METRICS}")


if __name__ == "__main__":
    main()
