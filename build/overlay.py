"""The corrections laid over a release: what the site shows of a study that has a corrected copy.

The release is never changed. A study that data/corrections.tsv corrects has a copy outside the repository, and LDSC
was run on the copies (corrections/): data/corrected_h2.tsv, data/corrected_rg.tsv. Here the corrected numbers take
the place of the release's, which are kept beside them; a study whose copy has genetic correlations where the
release had none becomes a signed study; and the panels of non-redundant traits are chosen again, by the release's
rule (corrections/panel.py).
"""
import math
from pathlib import Path

import pandas as pd

from build.release import RG_COLUMNS
from corrections import fix
from corrections.__main__ import planned
from corrections.panel import selected
from papers import cite

KEPT = ["n", "m", "m_rg", "signed", "cross", "h2", "h2_se", "h2_z", "intercept", "panel"]  # what a correction may change of a study, kept as the release has it


def finite(text) -> bool:
    try:
        return math.isfinite(float(text))
    except (TypeError, ValueError):
        return False


def laid(d: pd.DataFrame, pairs: dict[str, pd.DataFrame], data: Path, copies: Path | None = None) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], dict[str, dict], dict[str, dict]]:
    """The studies and the tables of pairs with the corrections laid over them; `fixed`, {study: its corrections as
    [fix, what it gives, why]}; and `was`, {study: the release's values of what changed}, with, under each ancestry,
    {pair as "a b": the release's [rg, se]} for the pairs estimated again. `copies`: the folder of the corrected
    copies, whose files then take the place of the release's in the studies' paths. Without data/corrected_h2.tsv
    nothing is laid: the corrections are listed only."""
    rows = [row for row in cite.table(data / "corrections.tsv") if row["id"] in d.index] if (data / "corrections.tsv").exists() else []
    fixed = {}
    for row in rows:
        fixed.setdefault(row["id"], []).append([row["fix"], row["value"], row["why"]])
    if not (data / "corrected_h2.tsv").exists():
        return d, pairs, fixed, {}
    plan = planned(rows)
    h2 = pd.read_csv(data / "corrected_h2.tsv", sep="\t", dtype=str).set_index("trait_identifier")
    if lacking := [study for study in plan if study not in h2.index]:
        raise ValueError(f"data/corrected_h2.tsv lacks {lacking[0]}: run python -m corrections.ldsc")
    rg = pd.read_csv(data / "corrected_rg.tsv", sep="\t", dtype=str, keep_default_na=False)
    sizes = {row["accession"]: float(row["n"]) for row in cite.table(data / "samples.tsv")} if (data / "samples.tsv").exists() else {}
    new = d.copy()
    for study, fixes in plan.items():
        for column in ("h2", "h2_se", "h2_z", "intercept"):
            new.at[study, column] = float(h2.at[study, column])
        new.at[study, "m"], new.at[study, "m_rg"] = int(h2.at[study, "m"]), pd.NA  # the SNPs of the copy; the release's count for rg is not the copy's
        for row in fixes:  # the N that a correction gives: its own, or that of its source in the GWAS Catalog
            said = fix.sourced(row["value"]) if row["fix"] == "source" else {}
            n = float(row["value"]) if row["fix"] == "n" else said.get("n") or sizes.get(said.get("source"))
            if n:
                new.at[study, "n"] = round(n)
    correlated = set(rg["trait1"]) | set(rg["trait2"])
    new["cross"] = new["cross"] | new.index.isin([study for study in plan if study in correlated])
    new["signed"] = new["signed"] | new["cross"]
    tables, was = {}, {}
    for anc, table in pairs.items():
        mine = rg[rg["trait1"].map(new["anc"]) == anc]
        again = {frozenset(pair) for pair in zip(mine["trait1"], mine["trait2"])}
        replaced = table[[frozenset(pair) in again for pair in zip(table["trait1"], table["trait2"])]]
        tables[anc] = pd.concat([table.drop(replaced.index), mine.set_axis(table.columns, axis=1)], ignore_index=True)
        was[anc] = {f"{a} {b}": [float(r) if finite(r) else None, float(se) if finite(se) else None]
                    for a, b, r, se in zip(replaced["trait1"], replaced["trait2"], replaced[RG_COLUMNS["rg"]], replaced[RG_COLUMNS["se"]])}
        if not any(new.at[study, "anc"] == anc for study in plan):  # nothing of this ancestry is corrected: its panel stands
            continue
        estimated = tables[anc].iloc[:, 2:].map(finite).all(axis=1)
        chosen = selected(new["h2_z"][new["h2_z"].map(finite)].to_dict(),
                          {frozenset(pair): float(r) if good else None for pair, r, good in zip(zip(tables[anc]["trait1"], tables[anc]["trait2"]), tables[anc][RG_COLUMNS["rg"]], estimated)})
        new.loc[new["anc"] == anc, "panel"] = pd.Series(range(1, len(chosen) + 1), index=chosen, dtype="Int64").reindex(new.index[new["anc"] == anc])
    if copies is not None:
        for study in plan:
            new.at[study, "h2_file"] = str(copies / "sumstats" / f"{study}.sumstats")
            new.at[study, "rg_file"] = new.at[study, "h2_file"] if new.at[study, "cross"] else None
    alike = lambda a, b: pd.isna(a) and pd.isna(b) or not pd.isna(a) and not pd.isna(b) and a == b
    for study in new.index:  # of a corrected study, what changed; of another, its place in the panel if it entered or left
        before = {column: d.at[study, column] for column in (KEPT if study in fixed else ["panel"]) if not alike(d.at[study, column], new.at[study, column])}
        if study not in fixed and pd.isna(d.at[study, "panel"]) == pd.isna(new.at[study, "panel"]):
            continue
        if before:
            was[study] = {column: None if pd.isna(value) else value.item() if hasattr(value, "item") else value for column, value in before.items()}
    return new, tables, fixed, was
