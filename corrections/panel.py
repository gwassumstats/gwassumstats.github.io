"""The release's rule for its panels of non-redundant traits, to apply it again once studies are corrected.

From the release's own script: of the studies that have genetic correlations, those in pairs that LDSC could not
estimate are left out one at a time, the one in most such pairs first, until every pair that remains was estimated.
Of the rest, the studies with an h² z above 6 are candidates. The candidate with the highest z is taken, and it and
every candidate whose rg with it has a square above 0.1 leave; and so on until none is left.
"""
MIN_Z, MAX_RG2 = 6, 0.1


def selected(z: dict[str, float], rg: dict[frozenset, float | None]) -> list[str]:
    """The panel, in the order its studies are taken. `z`: the h² z of each study; `rg`: the genetic correlation of
    each pair of them, None for a pair that could not be estimated."""
    active = {study for pair in rg for study in pair}
    failed = [pair for pair, value in rg.items() if value is None]
    while bad := [pair for pair in failed if pair <= active]:
        degree = {}
        for pair in bad:
            for study in pair:
                degree[study] = degree.get(study, 0) + 1
        active.remove(max(degree, key=lambda study: (degree[study], study)))  # in most pairs that failed; of several, the last by name
    remain = {study: z[study] for study in sorted(active) if study in z and z[study] > MIN_Z}
    panel = []
    while remain:
        best = max(remain, key=remain.get)
        panel.append(best)
        remain = {study: value for study, value in remain.items() if study != best and rg[frozenset((study, best))] ** 2 <= MAX_RG2}
    return panel
