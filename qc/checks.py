"""Judges every study: pass, warn or fail, and why.

What is judged: the study's files, measured SNP by SNP (qc/scan.py -> data/metrics.tsv); its LDSC results, from the
release; its label, against its genetic correlations (qc/labels.py); and its row of the inventory against the GWAS
Catalog's record (data/catalog.tsv). The thresholds are in data/checks.toml. A finding is a fail (do not use the
study as it is), a warn (use it with care) or a note (worth knowing; the study passes). Its text names another study
as {ID}: the site shows that study's short name there.
"""
from itertools import combinations

import pandas as pd

from papers import cite
from qc import labels

LEVELS = ["fail", "warn", "note"]  # the worst first; a study with notes only passes


def files(studies: pd.DataFrame, metrics: pd.DataFrame, rules: dict):
    """Findings on the files of each study: (study, level, text). Of a study's copy for rg, only what its file lacks."""
    for (study, kind), m in metrics.iterrows():
        s, of = studies.loc[study], "its copy for rg: " if kind == "rg" else ""
        columns = "SNP A1 A2 N Z" if s["alleles"] else "SNP N Z"  # the inventory says whether a study has alleles
        if m["columns"] != columns:
            yield study, "fail", f"{of}its columns are {m['columns']}, not {columns}"
        for count, what in ((m["duplicates"], "listed twice"), (m["z_missing"], "with a Z that is not a number")):
            if count:
                yield study, "fail" if count > rules["max_bad"] * m["rows"] else "warn", f"{of}{count:,.0f} SNPs {what}"
        for count, level, what in ((m["not_rs"], "warn", "without an rsID"), (m["bad_alleles"], "fail", "with alleles that are not two of A, C, G, T"),
                                   (m["n_missing"], "warn", "without an N")):
            if count > 0:
                yield study, level, f"{of}{count:,.0f} SNPs {what}"
        if kind == "rg":
            if m["ambiguous"] > 0:
                yield study, "warn", f"{of}{m['ambiguous']:,.0f} strand-ambiguous SNPs"
            continue  # what follows is said once a study, of its file
        if m["rows"] != s["m"]:
            yield study, "warn", f"{m['rows']:,} SNPs in the file; the inventory says {s['m']:,}"
        if abs(m["n_median"] / s["n"] - 1) > rules["n_within"]:
            yield study, "note", f"the median N of its SNPs is {m['n_median']:,.0f}; the inventory says {s['n']:,}"
        if m["in_hm3"] < rules["min_snps"]:
            yield study, "fail", f"only {m['in_hm3']:,} of the {m['hm3']:,} HapMap3 SNPs that LDSC regresses on: too few for LDSC"
        elif m["in_hm3"] < rules["min_coverage"] * m["hm3"]:
            yield study, "warn", f"{m['in_hm3'] / m['hm3']:.0%} of the HapMap3 SNPs that LDSC regresses on"


def copies(metrics: pd.DataFrame):
    """A warn for each study whose file is, by every number measured of it, the file of another study."""
    files = metrics.xs("h2", level="file")
    same = files.groupby([column for column in files.columns if column != "columns"], dropna=False).groups
    for group in (list(group) for group in same.values() if len(group) > 1):
        for study in group:
            yield study, "warn", "by every number measured, its file is that of " + " and ".join(f"{{{other}}}" for other in group if other != study) + ": one file under two IDs"


def statistics(studies: pd.DataFrame, metrics: pd.DataFrame, rules: dict):
    """Findings on the Z scores of each study and on its LDSC results: (study, level, text)."""
    low, high = rules["signed_positive"]
    for study, s in studies.join(metrics.xs("h2", level="file"), how="inner").iterrows():
        holds_abs = s["z_positive"] > rules["unsigned_positive"]
        if s["signed"] and holds_abs:
            yield study, "fail", f"signed in the inventory, but {s['z_positive']:.0%} of its Z are positive: it holds |Z|"
        elif not s["signed"] and not holds_abs:
            yield study, "warn", f"unsigned in the inventory, but {s['z_positive']:.0%} of its Z are positive: they look signed"
        elif s["signed"] and not low <= s["z_positive"] <= high:
            yield study, "warn", f"{s['z_positive']:.0%} of its Z are positive"
        if s["chi2_mean"] < rules["min_mean_chi2"]:
            yield study, "warn", f"mean χ² {s['chi2_mean']:.2f}: no signal above chance" + ("; the statistics are deflated" if s["chi2_mean"] < rules["deflated_chi2"] else "")
        elif s["intercept"] >= rules["high_intercept"] and (s["intercept"] - 1) / (s["chi2_mean"] - 1) > rules["max_ratio"]:
            yield study, "warn", f"LDSC intercept {s['intercept']:.2f} with a mean χ² of {s['chi2_mean']:.2f}: most of the inflation is not polygenic"
        if s["intercept"] < rules["min_intercept"]:
            yield study, "note", f"LDSC intercept {s['intercept']:.2f}: the statistics look corrected for inflation (genomic control), which lowers h²"
        if s["h2"] <= 0:
            yield study, "fail", f"h² {s['h2']:.4f}: not above 0, so nothing for LDSC to work with"
        elif s["h2_z"] < rules["min_h2_z"]:
            yield study, "warn", f"h² z {s['h2_z']:.2f}: h² is within {rules['min_h2_z']} standard errors of 0"
        if s["h2"] > 1:
            yield study, "warn", f"h² {s['h2']:.2f}: above 1"


def twins(studies: pd.DataFrame, tables: dict, rules: dict):
    """A note for each study that is nearly identical to others: how many, and the largest of them."""
    near = {}
    for t in tables.values():
        for (a, b), rg, se in zip(combinations(t["ids"], 2), t["rg"], t["se"], strict=True):
            if rg is not None and rg >= rules["twin_rg"] and se <= rules["twin_se"]:
                near.setdefault(a, []).append(b)
                near.setdefault(b, []).append(a)
    for study, others in near.items():
        largest = max(others, key=lambda other: studies.at[other, "n"])
        larger = studies.at[largest, "n"] > studies.at[study, "n"]
        yield study, "note", (f"nearly identical to {len(others)} other {'study' if len(others) == 1 else 'studies'} (rg of at least {rules['twin_rg']})"
                              + (f": {{{largest}}}" if len(others) == 1 else f"; the largest is {{{largest}}}" if larger else "; none is larger")
                              + (f", N {studies.at[largest, 'n']:,}" if larger else ""))


def metadata(studies: pd.DataFrame, catalog: pd.DataFrame):
    """Findings on the inventory against the GWAS Catalog, for the studies that are in it: (study, level, text)."""
    names = {"EUR": "European", "EAS": "East Asian"}
    letters = lambda text: cite.plain(text).lower()
    for study, c in catalog.iterrows():
        s = studies.loc[study]
        if int(c["n"]) != s["n"]:
            yield study, "warn", f"N {s['n']:,}; the GWAS Catalog says {int(c['n']):,}"
        if c["ancestry"] != names.get(s["anc"], s["anc"]):
            yield study, "fail", f"ancestry {s['anc']}; the GWAS Catalog says {c['ancestry'] or 'nothing'}"
        if letters(c["author"].rsplit(" ", 1)[0]) not in letters(s["ref"]) or c["published"][:4] not in s["ref"]:
            yield study, "warn", f"reference {s['ref']}; the GWAS Catalog says {c['author']} {c['published'][:4]}"


def judge(studies: pd.DataFrame, tables: dict, metrics: pd.DataFrame, catalog: pd.DataFrame, reviewed: dict[str, tuple[str, str]],
          by_hand: list[tuple[str, str, str]], rules: dict) -> pd.DataFrame:
    """One row per study, in their order: `qc` (pass, warn or fail) and `notes`, its findings as [level, text], the
    worst first. `metrics`: by study and file ("h2", "rg"); `catalog`: by study; `reviewed`: {study: (the trait its
    file holds, how that is known)} for the studies whose label in the release a person found wrong; `by_hand`: what a
    person found of a study, as (study, level, text), counted like a finding of the checks. A study that was not
    measured warns."""
    measured = metrics.index.get_level_values("id")
    found = [*files(studies, metrics[measured.isin(studies.index)], rules["file"]), *copies(metrics[measured.isin(studies.index)]), *statistics(studies, metrics, rules["statistics"]),
             *metadata(studies, catalog[catalog.index.isin(studies.index)]), *twins(studies, tables, rules["label"])]
    found += [(study, "warn", "its file was not measured") for study in studies.index.difference(measured)]
    found += [(study, "warn", f"the release calls it {studies.at[study, 'trait']}; its file holds {label}: {why}") for study, (label, why) in reviewed.items() if study in studies.index]
    if wrong := [level for _, level, _ in by_hand if level not in LEVELS]:
        raise ValueError(f"a finding made by hand is a {wrong[0]}: expected one of {', '.join(LEVELS)}")
    found += [finding for finding in by_hand if finding[0] in studies.index]
    signed = studies[studies["cross"]].assign(trait=lambda t: [reviewed[study][0] if study in reviewed else trait for study, trait in t["trait"].items()])
    found += [(study, level, text) for study, findings in labels.check(signed[["trait", "h2_z"]].to_dict("index"), tables, rules["label"]).items() for level, text in findings]
    notes = {study: [] for study in studies.index}
    for study, level, text in sorted(found, key=lambda finding: LEVELS.index(finding[1])):  # a stable sort: within a level, in the order found
        notes[study].append([level, text])
    worst = lambda findings: next((level for level, _ in findings if level != "note"), "pass")
    return pd.DataFrame({"qc": [worst(notes[study]) for study in studies.index], "notes": [notes[study] for study in studies.index]}, index=studies.index)
