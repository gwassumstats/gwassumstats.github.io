"""The effective N of a study, and the study to use of each trait and ancestry.

A study of a binary trait has the power of its cases and controls, not of their sum: its effective N is
4 / (1 / cases + 1 / controls), which is 4 v (1 - v) N for a share v of cases among N individuals. A quantitative
study's is its N. Of the studies of one trait and one ancestry that pass the quality checks, the one of the largest
effective N is recommended, or of the largest h² z where the effective N of one of them is not known; an older study
of that trait and ancestry, and of the same sex, is superseded by it.
"""
import math

import pandas as pd

SEXES = ("Female", "Male")  # what a short name says of a study of one sex, where its trait is of both


def effective(kind: str, n: float, cases: float | None = None, controls: float | None = None, largest: float | None = None) -> float | None:
    """The effective N of a study of `n` individuals: `n` if its trait is quantitative; 4 v (1 - v) times its
    individuals if it is binary, v the share of cases among `cases` and `controls`; None if its cases are not known,
    or its kind. The cases and controls recorded may be those of a larger sample than the file holds, of which the
    file is then taken to have the same share of cases: the individuals counted are theirs, but no more than the
    study's `n`, or than the `largest` N of its file's SNPs if that is above it."""
    if kind == "quantitative":
        return float(n)
    if kind != "binary" or not cases or not controls or math.isnan(cases) or math.isnan(controls):
        return None
    total, share = cases + controls, cases / (cases + controls)
    return 4 * share * (1 - share) * min(total, max(n, largest or 0))


def sex(name: str) -> str:
    """The sex that a study's short name says it is of: Female, Male, or "" for both."""
    return next((word for word in SEXES if word in name.split("_", 2)[-1]), "")


def older(a: str, b: str) -> bool:
    """Whether the day `a` is before the day `b`, as far as both are known: 2018-07 is not before 2018-07-15."""
    k = min(len(a), len(b))
    return bool(k) and a[:k] < b[:k]


def chosen(d: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Which study is recommended of each trait and ancestry, and which study supersedes each. `d`, by study: `of`
    (its trait), `anc`, `qc`, `n`, `neff` (its effective N, missing if not known), `h2_z`, `name` and `published` (the
    day its paper was published, or its month; "" if none is known).

    The first: by study, "neff" or "z" for a recommended study, after what it was chosen by, and None for the others.
    Of a trait and an ancestry, the studies that pass are ranked by effective N, the largest first, then by h² z. If
    the effective N of one of them is not known, as where a binary study's cases are not public, they are all ranked
    by h² z instead, then by N: an effective N is no measure against a study that has none, and N alone overstates
    the power of a study of few cases among many controls; h² z is the power that LDSC measured. The first is
    recommended.

    The second: by study, the ID of the study that supersedes it, or None. Every other study of the trait and
    ancestry whose paper is older than the recommended study's, and which is of the same sex, is superseded by it: a
    study of one sex, or a newer one, is another analysis, not an older version."""
    recommended, superseded = pd.Series(None, index=d.index, dtype=object), pd.Series(None, index=d.index, dtype=object)
    for _, group in d.groupby(["of", "anc"], sort=False):
        passing = group[group["qc"] == "pass"]
        if passing.empty:
            continue
        known = passing["neff"].notna().all()
        best = passing.sort_values(["neff", "h2_z"] if known else ["h2_z", "n"], ascending=False, kind="stable").index[0]
        recommended[best] = "neff" if known else "z"
        when, whose = d.at[best, "published"], sex(d.at[best, "name"])
        for study, row in group.iterrows():
            if study != best and older(row["published"], when) and sex(row["name"]) == whose:
                superseded[study] = best
    return recommended, superseded
