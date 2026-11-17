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
    """Each SNP's N scaled so that the most frequent N becomes `n`: a SNP that every cohort has gets the paper's N,
    and a SNP with fewer cohorts keeps its share."""
    most = float(rows["N"].mode()[0])
    return rows.assign(N=rows["N"].map(lambda value: f"{float(value) * n / most:.3f}")), {"n_scaled_by": round(n / most, 6), "n_most_frequent_before": most}


FIXES = {"flip": lambda rows, value: flip(rows), "clean": lambda rows, value: clean(rows), "n": lambda rows, value: rescale(rows, float(value))}
