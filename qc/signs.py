"""The sign test: for every pair of signed studies, at the SNPs that are strong in both, how often their Z have the
same sign once the alleles are aligned. On Slurm.

  uv run --group qc python -m qc.signs RELEASE --alleles FILE --logs FOLDER [--account NAME]

Two studies of one trait agree at nearly every such SNP; two of opposite traits, at nearly none. A study whose Z does
not go with its alleles agrees with every other at about half: its sign carries nothing. RELEASE is the release folder
(default: $SUMSTATS_RELEASE); FILE lists the HapMap3 SNPs with their alleles (SNP, A1, A2): each file's Z is turned to
the first allele of the list, so that every pair is compared at once. FOLDER takes the Slurm logs: keep it outside
the repository. The result, data/signs.tsv, holds no path: a line a pair, with the SNPs strong in both and those of
the same sign; of each study's pairs, the clearest, the largest, and those with a study of its trait.
"""
import argparse
import datetime
import os
from pathlib import Path

import numpy as np
import pandas as pd

from papers import cite
from sources.slurm import submitted

ROOT = Path(__file__).parents[1]
SIGNS = ROOT / "data" / "signs.tsv"
STRONG, LEAST, KEPT = 4, 50, 3  # |Z| above which a SNP is strong; the SNPs a pair needs to count; the pairs kept of a study, of each kind
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def oriented(path: str, alleles: pd.DataFrame) -> pd.Series:
    """Of a file's SNPs with |Z| above STRONG, by SNP: +1 if the first allele of the list raises the trait there, -1 if
    it lowers it. A SNP whose alleles are not the list's two, on either strand, or are so either way round, is left out."""
    table = pd.read_csv(path, sep="\t", usecols=["SNP", "A1", "A2", "Z"], dtype={"SNP": str, "A1": str, "A2": str})
    z = pd.to_numeric(table["Z"], errors="coerce")
    strong = table.assign(Z=z)[np.isfinite(z) & (z.abs() > STRONG) & table["SNP"].isin(alleles.index)].drop_duplicates("SNP")
    a1, a2, listed = strong["A1"].str.upper(), strong["A2"].str.upper(), alleles.loc[strong["SNP"]]
    r1, r2, c1, c2 = listed["A1"].to_numpy(), listed["A2"].to_numpy(), a1.str.translate(COMPLEMENT), a2.str.translate(COMPLEMENT)
    as_listed, turned = (a1 == r1) & (a2 == r2) | (c1 == r1) & (c2 == r2), (a1 == r2) & (a2 == r1) | (c1 == r2) & (c2 == r1)
    sign = np.sign(strong["Z"]) * np.where(as_listed, 1, -1)
    return pd.Series(sign.to_numpy(), index=strong["SNP"].to_numpy(), dtype=int)[(as_listed != turned).to_numpy()]


def compared(files: dict[str, str], snplist: str) -> pd.DataFrame:
    """One Slurm task: every pair of files with at least LEAST SNPs strong in both, as a, b, the number of those SNPs
    (n) and the number at which the two Z have the same sign (same)."""
    alleles = pd.read_csv(snplist, sep="\t", dtype=str).drop_duplicates("SNP").set_index("SNP")
    place = pd.Series(np.arange(len(alleles)), index=alleles.index)
    signs = np.zeros((len(files), len(alleles)), dtype=np.int8)
    for k, (study, path) in enumerate(files.items()):
        found = oriented(path, alleles)
        signs[k, place[found.index].to_numpy()] = found.to_numpy()
        print(f"[{k + 1}/{len(files)}] {study}: {len(found):,} strong SNPs", flush=True)
    both, product = np.zeros((len(files), len(files))), np.zeros((len(files), len(files)))
    for start in range(0, signs.shape[1], 100_000):  # a block of SNPs at a time: the whole table as numbers would not fit
        block = signs[:, start:start + 100_000].astype(np.float32)
        both += np.abs(block) @ np.abs(block).T
        product += block @ block.T
    a, b = np.triu_indices(len(files), 1)
    pairs = pd.DataFrame({"a": np.array(list(files))[a], "b": np.array(list(files))[b], "n": both[a, b].round().astype(int), "same": ((both + product) / 2)[a, b].round().astype(int)})
    return pairs[pairs["n"] >= LEAST].reset_index(drop=True)


def chosen(pairs: pd.DataFrame, traits: dict[str, str]) -> pd.DataFrame:
    """Of every study's pairs, those worth keeping: the KEPT clearest (the share of the same sign farthest from a
    half; of two alike, the pair of more SNPs), the KEPT of most SNPs, and every pair with a study of its trait."""
    both = pd.concat([pairs, pairs.rename(columns={"a": "b", "b": "a"})], ignore_index=True)  # each pair from either study
    both = both.assign(clear=(both["same"] / both["n"] - 0.5).abs(), kin=both["a"].map(traits) == both["b"].map(traits))
    keep = pd.concat([both.sort_values(["clear", "n"], ascending=False).groupby("a").head(KEPT), both.sort_values("n", ascending=False).groupby("a").head(KEPT), both[both["kin"]]])
    keys = {frozenset(pair) for pair in zip(keep["a"], keep["b"])}
    return pairs[[frozenset(pair) in keys for pair in zip(pairs["a"], pairs["b"])]].reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(prog="python -m qc.signs", description="Compares the signs of every pair of signed studies at the SNPs strong in both, on Slurm.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--alleles", required=True, metavar="FILE", help="the HapMap3 SNPs with their alleles: SNP, A1, A2")
    parser.add_argument("--logs", type=Path, required=True, help="a folder outside the repository, for the Slurm logs")
    parser.add_argument("--account", help="the Slurm account to charge")
    parser.add_argument("--name", default=f"{datetime.date.today():%Y%m%d}_qc_signs", help="the name of the Slurm task and of its logs")
    args = parser.parse_args()
    if not args.release:
        parser.error("name the release folder, or set SUMSTATS_RELEASE")
    release = Path(args.release)
    index = pd.read_csv(release / "_audit" / "sumstats_file_index.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("Trait Identifier")["H2 Sumstats"]
    inventory = pd.read_csv(next(release.glob("signed_summary_statistics_*.txt")), sep="\t", dtype=str, keep_default_na=False)
    signed = [study.strip() for study, alleles in zip(inventory["Trait Identifier"], inventory["Allele Info"]) if alleles.strip() == "YES"]
    pairs = submitted(args.name, args.logs.resolve(), args.account, compared, [{study: index[study] for study in signed}], [args.alleles], memory=16, minutes=120)[0]
    kept = chosen(pairs, {row["id"]: row["trait"] for row in cite.table(ROOT / "data" / "traits.tsv")})
    kept.to_csv(SIGNS, sep="\t", index=False)
    print(f"{len(signed)} signed studies: {len(pairs):,} pairs with at least {LEAST} SNPs strong in both; {len(kept):,} kept in {SIGNS.name}")


if __name__ == "__main__":
    main()
