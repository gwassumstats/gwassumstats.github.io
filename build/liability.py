"""h² on the liability scale: of a binary study, with the population prevalence of its trait that a paper of the
collection states (data/prevalence.tsv).

LDSC estimates h² on the observed scale, of a sample with a share v of cases. On the liability scale (Lee et al. 2011,
PMID 21376301) it is h² K² (1 - K)² / (φ(t)² v (1 - v)), where K is the population prevalence, t the threshold that
a share K of a standard normal lies above, and φ its density. LDSC divides by the N of the file: an N above the
individuals of its cases and controls, as where it counts some twice, lowers h² as many times, which is undone; an N
below them is taken to be of a part of the sample, with its share of cases, as for the effective N
(build/recommend.py). A study some of whose cases are proxy cases, relatives of a case, has none: its cases are not
those of the trait.
"""
import math
from pathlib import Path
from statistics import NormalDist

from papers import cite


def factor(prevalence: float, cases: float, controls: float, n: float | None = None) -> float:
    """What h² on the observed scale is multiplied by on the liability scale, for a study of `n` individuals (the N
    LDSC used) whose sample has these cases and controls."""
    v, normal = cases / (cases + controls), NormalDist()
    density = normal.pdf(normal.inv_cdf(1 - prevalence))
    return prevalence ** 2 * (1 - prevalence) ** 2 / (density ** 2 * v * (1 - v)) * max(1.0, (n or 0) / (cases + controls))


def prevalences(data: Path) -> dict[tuple[str, str], dict]:
    """{(trait, ancestry): {k, source, quote}} of data/prevalence.tsv: the population prevalence that a paper states,
    in its words. One row a trait and ancestry."""
    table = {}
    for row in cite.table(data / "prevalence.tsv"):
        key, k = (row["trait"], row["anc"]), float(row["prevalence"])
        if key in table:
            raise ValueError(f"data/prevalence.tsv: {key[0]} ({key[1]}) has two rows")
        if not 0 < k < 1:
            raise ValueError(f"data/prevalence.tsv: the prevalence of {key[0]} ({key[1]}) is not between 0 and 1")
        table[key] = {"k": k, "source": row["source"], "quote": row["quote"]}
    return table


def scaled(kind: str, proxy: bool, cases: float, controls: float, n: float, h2: float, se: float, prevalence: float | None) -> tuple[float, float] | tuple[None, None]:
    """A study's h² and its standard error on the liability scale, or None where it has none: its trait is not
    binary, some of its cases are proxy cases, its cases are not known, or no prevalence is."""
    if kind != "binary" or proxy or prevalence is None or not cases or not controls or math.isnan(cases) or math.isnan(controls):
        return None, None
    f = factor(prevalence, cases, controls, n)
    return h2 * f, se * f
