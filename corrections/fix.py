"""The fixes a file can get, on its rows as text: nothing but what a fix changes is rewritten."""
import math

import pandas as pd


def read(path) -> pd.DataFrame:
    """A sumstats file as text, cell by cell: its numbers keep the digits they were written with."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def finite(z: pd.Series) -> pd.Series:
    """Whether each Z is a number: not empty, not NA, not infinite."""
    return pd.to_numeric(z, errors="coerce").map(math.isfinite)


def flip(rows: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Z with its sign turned, digit for digit: for a file whose Z was signed for A2."""
    turned = rows["Z"].map(lambda z: z[1:] if z.startswith("-") else z if float(z) == 0 else f"-{z}")
    return rows.assign(Z=turned), {"z_turned": int((turned != rows["Z"]).sum())}


def clean(rows: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Without the rows whose Z is not a number, then without the later rows of a SNP listed twice: the first stays,
    as LDSC's munge keeps it."""
    numbers = rows[finite(rows["Z"])]
    once = numbers.drop_duplicates("SNP")
    return once, {"z_not_a_number": len(rows) - len(numbers), "listed_twice": len(numbers) - len(once)}


def rescale(rows: pd.DataFrame, n: float) -> tuple[pd.DataFrame, dict]:
    """Each SNP's N scaled so that the largest N becomes `n`: a SNP that every cohort has gets the paper's count of
    cases and controls, and a SNP with fewer cohorts keeps its share. For a file whose N is an effective N."""
    largest = pd.to_numeric(rows["N"]).max()
    return rows.assign(N=rows["N"].map(lambda value: f"{float(value) * n / largest:.3f}")), {"n_scaled_by": round(n / largest, 6), "n_largest_before": float(largest)}


def sourced(value: str) -> dict:
    """What a `source` fix says, as in "GCST006906; a1=Allele1; a2=Allele2; n=446696": where the study's statistics
    are (`source`: an accession of the GWAS Catalog, or a URL); the file's name in the study's folder, if the folder
    has several (`file`); the N of every SNP, in place of the file's (`n`, a number); and the columns that play the
    roles of sources/format.py, where their names do not say so (`roles`; n=Column for the column that holds N)."""
    source, *settings = [part.strip() for part in value.split(";") if part.strip()]
    if lone := [setting for setting in settings if "=" not in setting]:
        raise ValueError(f"data/corrections.tsv: {lone[0]} of the source {source} is not key=value")
    given = {key.strip(): column.strip() for key, column in (setting.split("=", 1) for setting in settings)}
    try:
        n = float(given.get("n", ""))
        del given["n"]
    except ValueError:  # no N given, or the name of the column that holds it
        n = None
    return {"source": source, "file": given.pop("file", ""), "n": n, "roles": given}


FIXES = {"flip": lambda rows, value: flip(rows), "clean": lambda rows, value: clean(rows), "n": lambda rows, value: rescale(rows, float(value))}
