"""Formats a file of summary statistics as the release does: SNP, A1, A2, N and Z, on the HapMap3 SNPs.

The rules are those of the release's own scripts for its 2026 additions, so that a file made here can stand beside
the release's:

  1. A row's SNP is the rsID found in its SNP column, if that is a HapMap3 SNP; otherwise the HapMap3 SNP at its
     chromosome and position (GRCh37), which may also be written into the SNP column ("1:752566:G:A").
  2. Its alleles must be the two of that SNP in the reference, in either order, on either strand.
  3. A1 is the effect allele and A2 the other, as the file has them. Z is the file's Z; without one, the effect over
     its standard error; without an effect, the logarithm of the odds ratio over the standard error.
  4. N is the file's, SNP by SNP; without one, the N given for the study. It is above 0.
  5. Of several rows of a SNP, the first that passes all of this is kept.

A file may be delimited by tabs or commas, as in the release, or by blanks. A column plays a role by its name, in the
release's order of preference (ROLES), unless it is named for the role. A row with more fields than the header is
left out, as in the release; one that ends early is read with its last columns empty, where the release leaves it out.
"""
import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

ROLES = {"snp": ["snp", "rs_id", "rsid", "variant_id", "variant", "markername", "marker", "id"],
         "chr": ["chromosome", "chr", "chrom", "#chrom"],
         "pos": ["base_pair_location", "bp", "pos", "position", "base_pair", "basepair"],
         "a1": ["effect_allele", "a1", "tested_allele", "alt", "alt_allele"],
         "a2": ["other_allele", "a2", "reference_allele", "ref_allele", "ref"],
         "n": ["n", "nmiss", "sample_size", "n_total", "total_sample_size", "effective_sample_size"],
         "z": ["z", "z_score", "zscore"],
         "beta": ["beta", "effect", "effect_size", "estimate"],
         "se": ["standard_error", "se", "stderr", "beta_se"],
         "or": ["odds_ratio", "or"]}
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CHUNK = 500_000  # rows read at a time


def reference(snplist: str, bim_prefix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The HapMap3 SNPs that the release formats to: their two alleles, by rsID (`snplist`: SNP, A1, A2); and where
    each is (chr, pos, SNP), from the reference panel's .bim files, one per chromosome (`bim_prefix`)."""
    alleles = pd.read_csv(snplist, sep="\t", dtype=str)
    alleles = alleles.assign(SNP=alleles["SNP"].str.lower(), A1=alleles["A1"].str.upper(), A2=alleles["A2"].str.upper()).set_index("SNP")
    bim = pd.concat(pd.read_csv(f"{bim_prefix}{chromosome}.bim", sep=r"\s+", header=None, usecols=[0, 1, 3], names=["chr", "SNP", "pos"], dtype=str) for chromosome in range(1, 23))
    bim = bim.assign(SNP=bim["SNP"].str.lower(), chr=chromosomes(bim["chr"]), pos=positions(bim["pos"]))
    return alleles, bim[bim["SNP"].isin(alleles.index)].reset_index(drop=True)


def number(text: pd.Series) -> pd.Series:
    """Each value as a number; what is not a finite number is missing."""
    return pd.to_numeric(text, errors="coerce").astype(float).replace([np.inf, -np.inf], np.nan)


def chromosomes(text: pd.Series) -> pd.Series:
    """Chromosomes as 1 to 22 ("CHR01" is 1); any other is ""."""
    named = text.str.strip().str.upper().str.removeprefix("CHR").str.lstrip("0")
    return named.where(named.isin([str(k) for k in range(1, 23)]), "")


def positions(text: pd.Series) -> pd.Series:
    """Positions as whole numbers ("752566.0" is 752566); what is not a number is ""."""
    return np.trunc(number(text)).astype("Int64").astype("string").fillna("").astype(object)


def named(header: list[str], roles: dict[str, str] | None = None) -> dict[str, str]:
    """The column that plays each role: the one named for it in `roles`, or the first of ROLES that the file has."""
    lower = {name.strip().lower(): name for name in reversed(header)}  # of two columns of one name, the first
    wrong = [f"{role}={column}" for role, column in (roles or {}).items() if role not in ROLES or column.strip().lower() not in lower]
    if wrong:
        raise ValueError(f"{wrong[0]}: no such role, or no such column among {', '.join(header)}")
    found = {role: lower[candidate] for role, candidates in ROLES.items() if (candidate := next((c for c in candidates if c in lower), None))}
    return found | {role: lower[column.strip().lower()] for role, column in (roles or {}).items()}


def matching(a1: pd.Series, a2: pd.Series, r1: pd.Series, r2: pd.Series) -> pd.Series:
    """Whether a row's alleles are the two of the reference, in either order, on either strand."""
    letters = a1.isin(list("ACGT")) & a2.isin(list("ACGT")) & r1.isin(list("ACGT")) & r2.isin(list("ACGT"))
    as_is = (a1 == r1) & (a2 == r2) | (a1 == r2) & (a2 == r1)
    c1, c2 = a1.str.translate(COMPLEMENT), a2.str.translate(COMPLEMENT)
    return letters & (as_is | (c1 == r1) & (c2 == r2) | (c1 == r2) & (c2 == r1))


def rows(chunk: pd.DataFrame, columns: dict[str, str], alleles: pd.DataFrame, where: pd.DataFrame, n: float | None) -> tuple[pd.DataFrame, dict]:
    """The rows of a part of a file that pass the rules, in the file's order, as SNP, A1, A2, N and Z (N and Z as
    numbers); and how many rows fell at each rule."""
    text = lambda role, at=slice(None): chunk.loc[at, columns[role]].str.strip() if role in columns else pd.Series("", index=chunk.loc[at].index)
    raw = text("snp")
    rs = raw.str.extract(r"(?i)(rs[0-9]+)", expand=False).str.lower()
    known = rs.isin(alleles.index)
    chromosome, position = chromosomes(text("chr")), positions(text("pos"))
    written = raw.str.extract(r"(?i)(?:chr)?([0-9]{1,2})[:_\-]([0-9]+)").fillna("")  # a position written into the SNP column, read if a row has no other
    lacking = (chromosome == "") | (position == "")
    chromosome, position = chromosome.where(~lacking, chromosomes(written[0])), position.where(~lacking, positions(written[1]))
    elsewhere = ~known & (chromosome != "") & (position != "")
    placed = pd.DataFrame({"row": chunk.index[elsewhere], "chr": chromosome[elsewhere], "pos": position[elsewhere]}).merge(where, on=["chr", "pos"])
    found = pd.concat([pd.DataFrame({"row": chunk.index[known], "SNP": rs[known]}), placed[["row", "SNP"]]]).sort_values("row", kind="stable", ignore_index=True)
    of = lambda role: text(role, found["row"]).reset_index(drop=True)  # a column, for the rows found
    found = found.assign(A1=of("a1").str.upper(), A2=of("a2").str.upper())
    fits = matching(found["A1"], found["A2"], alleles["A1"].reindex(found["SNP"]).reset_index(drop=True), alleles["A2"].reindex(found["SNP"]).reset_index(drop=True))
    fits &= ~found["row"].where(fits).duplicated(keep=False)  # of several SNPs at a row's position, the only one whose alleles fit
    counts = {"rows": len(chunk), "not_hapmap3": len(chunk) - found["row"].nunique(), "other_alleles": found["row"].nunique() - int(fits.sum())}
    found = found[fits].reset_index(drop=True)
    z = number(of("z")) if "z" in columns else pd.Series(np.nan, index=found.index)
    if "se" in columns:
        error = number(of("se")).where(lambda e: e > 0)
        if "beta" in columns:
            z = z.fillna(number(of("beta")) / error)
        if "or" in columns:
            z = z.fillna(np.log(number(of("or")).where(lambda ratio: ratio > 0)) / error)
    found = found.assign(N=number(of("n")) if "n" in columns else n, Z=z)
    sized, signed = found["N"] > 0, np.isfinite(found["Z"])
    return found.loc[sized & signed, ["SNP", "A1", "A2", "N", "Z"]], counts | {"without_n": int((~sized).sum()), "without_z": int((sized & ~signed).sum())}


def formatted(source: Path, alleles: pd.DataFrame, where: pd.DataFrame, n: float | None = None, roles: dict[str, str] | None = None, own_n: bool = True) -> tuple[pd.DataFrame, dict]:
    """A file of summary statistics in the release's format, as text ready to write: SNP, A1, A2, N and Z. And what was
    done: the columns used, the rows read, the rows that fell at each rule, and the rows written. `n`: the N of the
    study, for a file without one, or, if not `own_n`, in place of the file's; `roles`: {role: column}, for a file
    whose columns are not named as ROLES expects."""
    with open(source, "rb") as file:
        zipped = file.read(2) == b"\x1f\x8b"
    with (gzip.open if zipped else open)(source, "rt", encoding="utf-8", errors="replace", newline="") as file:
        first = file.readline().rstrip("\r\n")
    delimiter = "\t" if "\t" in first else "," if "," in first else r"\s+"
    header = first.split() if delimiter == r"\s+" else next(csv.reader([first], delimiter=delimiter))
    columns = {role: column for role, column in named(header, roles).items() if own_n or role != "n"}
    if not {"a1", "a2"} <= set(columns):
        raise ValueError(f"no alleles among the columns {', '.join(header)}")
    if "z" not in columns and not ({"beta", "se"} <= set(columns) or {"or", "se"} <= set(columns)):
        raise ValueError(f"no Z, and no effect or odds ratio with a standard error, among the columns {', '.join(header)}")
    if "n" not in columns and n is None:
        raise ValueError(f"no N among the columns {', '.join(header)}, and none given for the study")
    kept, counts = [], {}
    for chunk in pd.read_csv(source, sep=delimiter, dtype=str, keep_default_na=False, usecols=list(set(columns.values())), chunksize=CHUNK, on_bad_lines="skip",
                             encoding_errors="replace", compression="gzip" if zipped else None):
        part, counted = rows(chunk, columns, alleles, where, n)
        kept.append(part)
        counts = {key: counts.get(key, 0) + value for key, value in counted.items()}
    table = pd.concat(kept, ignore_index=True)
    once = table.drop_duplicates("SNP")
    done = {"columns": columns, "n_from": columns.get("n", "the N given for the study"),
            "z_from": columns["z"] if "z" in columns else " / ".join([columns["beta"] if "beta" in columns else f"log {columns['or']}", columns["se"]]),
            **counts, "listed_again": len(table) - len(once), "written": len(once)}
    return once.assign(N=once["N"].map("{:.10g}".format), Z=once["Z"].map("{:.12g}".format)), done
