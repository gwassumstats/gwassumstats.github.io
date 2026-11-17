"""Checks that sources/ formats a file as the release does, by making studies of the release again from their source.

  SUMSTATS_CORRECTIONS=<a folder outside the repository> uv run --group qc python -m sources ACCESSION ... \\
      --snplist FILE --bim-prefix PREFIX [--account NAME]

ACCESSION: a study that the release formatted from its file in the GWAS Catalog (the 2026 additions; RELEASE, by
default $SUMSTATS_RELEASE, is only read). FILE lists the HapMap3 SNPs with their alleles (SNP, A1, A2) and PREFIX
names the reference panel's .bim files, one per chromosome: those the release formatted with. Each study gets a Slurm
task, which downloads its file again, formats it and compares the result with the release's file. Under the folder:

  sources/<ACCESSION>/        the file as downloaded, kept for a later run
  check/<ACCESSION>.sumstats  the file made again; check/<ACCESSION>.json: its source, what the formatting did, and
                              how it compares with the release's
  check.tsv                   a line a study: identical or not, and if not, where the two differ
"""
import argparse
import datetime
import json
import os
from pathlib import Path

import pandas as pd

from build import release as rel
from sources.make import made
from sources.slurm import submitted

ROOT = Path(__file__).parents[1]


def again(accession: str, theirs: str, n: float, out: Path, snplist: str, bim_prefix: str) -> dict:
    """One Slurm task: makes a study of the release again and compares the two files, as text."""
    ours = out / "check" / f"{accession}.sumstats"
    record = made(accession, ours, snplist, bim_prefix, n)
    a, b = (pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False) for path in (ours, theirs))
    both = a.merge(b, on="SNP", suffixes=("", "_release"))
    line = {"study": accession, "identical": ours.read_bytes() == Path(theirs).read_bytes(), "rows": len(a), "rows_release": len(b), "in_both": len(both),
            "same_alleles": int(((both["A1"] == both["A1_release"]) & (both["A2"] == both["A2_release"])).sum()), "same_n": int((both["N"] == both["N_release"]).sum()),
            "same_z": int((both["Z"] == both["Z_release"]).sum()), "largest_z_difference": float((both["Z"].astype(float) - both["Z_release"].astype(float)).abs().max()),
            "n_from": record["formatting"]["n_from"], "z_from": record["formatting"]["z_from"]}
    (out / "check" / f"{accession}.json").write_text(json.dumps(record | {"against_the_release": line}, indent=1) + "\n")
    print(line, flush=True)
    return line


def main():
    parser = argparse.ArgumentParser(prog="python -m sources", description="Makes studies of the release again from the GWAS Catalog, and compares.")
    parser.add_argument("accessions", nargs="+", metavar="ACCESSION", help="a study that the release formatted from the GWAS Catalog's file")
    parser.add_argument("--release", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder, read only (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--out", default=os.environ.get("SUMSTATS_CORRECTIONS"), help="where to write (default: $SUMSTATS_CORRECTIONS)")
    parser.add_argument("--snplist", required=True, help="the HapMap3 SNPs and their alleles: SNP, A1, A2")
    parser.add_argument("--bim-prefix", required=True, help="the reference panel's .bim files, up to the chromosome's number")
    parser.add_argument("--account", help="the Slurm account to charge")
    parser.add_argument("--name", default=f"{datetime.date.today():%Y%m%d}_sources_check", help="the name of the Slurm tasks and of their logs")
    args = parser.parse_args()
    if not args.release or not args.out:
        parser.error("name the release folder and the folder to write to, or set SUMSTATS_RELEASE and SUMSTATS_CORRECTIONS")
    release, out = Path(args.release).resolve(), Path(args.out).resolve()
    if out.is_relative_to(ROOT) or out.is_relative_to(release) or release.is_relative_to(out):
        parser.error("nothing is written inside the repository, nor inside the release")
    files = pd.read_csv(release / "_audit" / "sumstats_file_index.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("Trait Identifier")["H2 Sumstats"]
    sizes = rel.studies(release)["n"]
    if unknown := [accession for accession in args.accessions if accession not in files.index]:
        parser.error(f"the release has no study {unknown[0]}")
    k = len(args.accessions)
    lines = submitted(args.name, out, args.account, again, args.accessions, [files[a] for a in args.accessions], [float(sizes[a]) for a in args.accessions], [out] * k, [args.snplist] * k, [args.bim_prefix] * k)
    pd.DataFrame(lines).to_csv(out / "check.tsv", sep="\t", index=False)
    print(f"{sum(line['identical'] for line in lines)} of {len(lines)} files are the release's, byte for byte")


if __name__ == "__main__":
    main()
