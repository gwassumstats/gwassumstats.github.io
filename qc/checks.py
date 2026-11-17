"""Judges every study: pass, warn or fail, and why.

What is judged: the study's files, measured SNP by SNP (qc/scan.py -> data/metrics.tsv); its LDSC results, from the
release; its label, against its genetic correlations (qc/labels.py); its row of the inventory against the GWAS
Catalog's record (data/catalog.tsv); and the N of its file against the samples the Catalog records (data/samples.tsv). The thresholds are in data/checks.toml. A finding is a fail (do not use the
study as it is), a warn (use it with care) or a note (worth knowing; the study passes). Its text names another study
as {ID}: the site shows that study's short name there.
"""
import re
from itertools import combinations

import pandas as pd

from papers import cite
from qc import labels

LEVELS = ["fail", "warn", "note"]  # the worst first; a study with notes only passes
ANCESTRIES = {"EUR": "European", "EAS": "East Asian"}  # as the GWAS Catalog names them


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
        if s["signed"] and m.get("list_reversed", 0) > rules["list_reversed"] * m["rows"]:
            yield study, "warn", (f"its alleles are those of the HapMap3 list, the other way round, at {m['list_reversed'] / m['rows']:.0%} of its SNPs: the mark of a formatting that set the alleles "
                                  "without regard to the effect allele. Where such a file was compared with its source, its Z went with its alleles at half of the SNPs only")
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


def recorded(row: pd.Series, ancestry: str) -> dict:
    """A sample that the GWAS Catalog records, as a study of `ancestry` sees it: its individuals (`n`); those who may
    be of that ancestry (`own`: every group that the Catalog gives it, alone or with others, or gives none); and its
    cases and controls of that ancestry, counted from the Catalog's words for the sample."""
    groups = [group.rsplit(" ", 1) for group in row["ancestries"].split("; ") if group]
    own = sum(int(count) for named, count in groups if {ancestry, "NR"} & set(named.split(", ")))  # NR: not reported
    counted = {"cases": 0, "controls": 0}
    for count, words, kind in re.findall(r"(\d[\d,]*) ([^,]*?)\b(cases|controls)\b", row["sample"]):
        if own == int(row["n"]) or ancestry in words:
            counted[kind] += int(count.replace(",", ""))
    return {"accession": row["accession"], "n": int(row["n"]), "own": own, **counted}


def sizes(studies: pd.DataFrame, metrics: pd.DataFrame, samples: pd.DataFrame, cited: dict[str, str], rules: dict):
    """Findings on the N of each study's file against the samples that the GWAS Catalog records: (study, level, text).

    The release counts individuals: cases and controls together. A file is compared by the largest N of its SNPs, that
    of a SNP every cohort has, and by their median. `samples`: the Catalog's records (data/samples.tsv). A study that
    is itself a record of the Catalog is compared with that record; another with every record of its paper (`cited`:
    {study: PMID}), for the checks cannot tell which of them it is."""
    accessions = set(samples["accession"])
    plain = lambda number: f"{number:,.2f}".removesuffix(".00")
    low, high = rules["n_of_paper"]
    for study, m in metrics.xs("h2", level="file").iterrows():
        if study not in studies.index:
            continue
        ancestry = ANCESTRIES.get(studies.at[study, "anc"], studies.at[study, "anc"])
        rows, whose = (samples[samples["accession"] == study], "it") if study in accessions else (samples[samples["pmid"] == cited.get(study, "")], "its paper")
        records = [record for record in (recorded(row, ancestry) for _, row in rows.iterrows()) if record["n"]]
        largest, median = m["n_max"], m["n_median"]
        totals = sorted({size for record in records for size in (record["n"], record["own"]) if size})
        same = lambda size: abs(largest / size - 1) <= rules["same_n"]
        derived = not all(float(m[column]).is_integer() for column in ("n_min", "n_median", "n_max")) and not any(map(same, totals))  # not whole numbers, nor the count of a sample
        mixed = [record for record in records if record["own"] < record["n"] and same(record["n"])]
        effective = [(4 / (1 / record["cases"] + 1 / record["controls"]) / halved, halved, record) for record in records if record["cases"] and record["controls"] for halved in (1, 2)]
        fits = [(largest / size, halved, record) for size, halved, record in effective if derived and rules["effective"][0] <= largest / size <= rules["effective"][1]]
        if mixed and largest > (1 + rules["other_ancestry"]) * max(record["own"] for record in records):
            yield study, "warn", (f"its N, {largest:,.0f} at most, is that of a sample of several ancestries: {mixed[0]['accession']} has {mixed[0]['n']:,} individuals, at most "
                                  f"{mixed[0]['own']:,} of them {ancestry}, and no sample recorded for {whose} has more of these. The file looks like the analysis of all ancestries together")
        elif fits:
            ratio, halved, record = min(fits, key=lambda fit: (round(abs(fit[0] - 1), 6), fit[2]["own"] < fit[2]["n"]))  # the nearest; of two alike, the sample of one ancestry
            total = record["cases"] + record["controls"]
            yield study, "warn" if total / largest - 1 > rules["n_within"] else "note", (
                f"its N looks like {'half an effective N, as PGC files carry (NEFFDIV2)' if halved == 2 else 'an effective N'}: {plain(largest)} at most, where the {record['cases']:,} cases and "
                f"{record['controls']:,} controls of {record['accession']} give {largest / ratio:,.0f}, and {total:,} in all. h² is {total / largest:.2f} times what their count gives")
        elif totals and not any(low <= median / size <= high for size in totals):
            yield study, "note", (f"the median N of its SNPs, {median:,.0f}, is {'above' if median > totals[-1] else 'far from'} every sample the GWAS Catalog records for {whose} ({totals[0]:,}"
                                  + (f" to {totals[-1]:,}" if len(totals) > 1 else "") + "): "
                                  + ("it may count samples that the statistics do not hold, such as a replication stage" if median > totals[-1] else "a part of the study, or an effective N"))
        elif derived:
            value, where = next((m[column], where) for column, where in (("n_max", "at most"), ("n_median", "for most SNPs"), ("n_min", "at least")) if not float(m[column]).is_integer())
            yield study, "note", f"its N is not whole numbers ({plain(value)} {where}): the weights of a meta-analysis, or an effective N"


def signs(studies: pd.DataFrame, pairs: pd.DataFrame, traits: dict[str, str], names: dict[str, str], faulty: set[str], rules: dict):
    """Findings from the sign test (qc/signs.py -> data/signs.tsv): (study, level, text). `pairs`: a, b, the SNPs
    with |Z| > 4 in both (n) and those at which the two Z have the same sign (same). Two studies of one trait
    (`traits`) and one ancestry agree at nearly all such SNPs: a pair that agrees at some only, or at none, warns
    both, unless one is `faulty`, found at fault by a person, which explains it, or the two are not of the same sex
    by their short `names`, where a trait may differ. Across ancestries two sound studies agree less, and are not
    judged. A study with many such SNPs in a pair, and no pair that agrees or disagrees clearly, warns too: its sign
    carries nothing."""
    pairs = pairs[pairs["a"].isin(studies.index) & pairs["b"].isin(studies.index)]
    pairs = pairs.assign(share=pairs["same"] / pairs["n"])
    sex = lambda study: next((word for word in ("Female", "Male") if word in names.get(study, "")), "")
    for a, b, n, share in pairs[["a", "b", "n", "share"]].itertuples(index=False):
        if traits.get(a) != traits.get(b) or studies.at[a, "anc"] != studies.at[b, "anc"] or {a, b} & faulty or sex(a) != sex(b) or share >= rules["agree"]:
            continue
        for study, other in ((a, b), (b, a)):
            if share <= 1 - rules["agree"]:
                yield study, "warn", (f"at the {n:,} SNPs with |Z| > 4 in both, its Z has the sign opposite to that of {{{other}}}, a study of its trait, for {1 - share:.0%}: "
                                      "one of the two is signed for the other allele, or measures the trait the other way round")
            else:
                yield study, "warn", (f"at the {n:,} SNPs with |Z| > 4 in both, its Z has the sign of that of {{{other}}}, a study of its trait, for {share:.0%} only: "
                                      "two studies of one trait agree at nearly all, unless the Z of one does not go with its alleles")
    both = pd.concat([pairs, pairs.rename(columns={"a": "b", "b": "a"})])
    for study, own in both[~both["b"].isin(faulty)].groupby("a"):  # a pair with a study at fault says nothing of the other
        if study not in faulty and own["n"].max() >= rules["clear_snps"] and ((own["share"] - 0.5).abs() < rules["unclear"]).all():
            yield study, "warn", (f"its sign is clear against no other study: of its {len(own)} pairs listed, one of {own['n'].max():,} SNPs with |Z| > 4 in both, the same sign at {own['share'].min():.0%} to "
                                  f"{own['share'].max():.0%}. A study whose Z does not go with its alleles looks so")


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
    letters = lambda text: cite.plain(text).lower()
    for study, c in catalog.iterrows():
        s = studies.loc[study]
        if int(c["n"]) != s["n"]:
            yield study, "warn", f"N {s['n']:,}; the GWAS Catalog says {int(c['n']):,}"
        if c["ancestry"] != ANCESTRIES.get(s["anc"], s["anc"]):
            yield study, "fail", f"ancestry {s['anc']}; the GWAS Catalog says {c['ancestry'] or 'nothing'}"
        if letters(c["author"].rsplit(" ", 1)[0]) not in letters(s["ref"]) or c["published"][:4] not in s["ref"]:
            yield study, "warn", f"reference {s['ref']}; the GWAS Catalog says {c['author']} {c['published'][:4]}"


def judge(studies: pd.DataFrame, tables: dict, metrics: pd.DataFrame, catalog: pd.DataFrame, reviewed: dict[str, tuple[str, str]],
          by_hand: list[tuple[str, str, str]], rules: dict, samples: pd.DataFrame | None = None, cited: dict[str, str] | None = None,
          compared: pd.DataFrame | None = None, traits: dict[str, str] | None = None, names: dict[str, str] | None = None) -> pd.DataFrame:
    """One row per study, in their order: `qc` (pass, warn or fail) and `notes`, its findings as [level, text], the
    worst first. `metrics`: by study and file ("h2", "rg"); `catalog`: by study; `reviewed`: {study: (the trait its
    file holds, how that is known)} for the studies whose label in the release a person found wrong; `by_hand`: what a
    person found of a study, as (study, level, text), counted like a finding of the checks. A study that was not
    measured warns."""
    measured = metrics.index.get_level_values("id")
    found = [*files(studies, metrics[measured.isin(studies.index)], rules["file"]), *copies(metrics[measured.isin(studies.index)]), *statistics(studies, metrics, rules["statistics"]),
             *metadata(studies, catalog[catalog.index.isin(studies.index)]), *twins(studies, tables, rules["label"])]
    found += [(study, "warn", "its file was not measured") for study in studies.index.difference(measured)]
    if samples is not None:
        found += sizes(studies, metrics, samples, cited or {}, rules["file"])
    found += [(study, "warn", f"the release calls it {studies.at[study, 'trait']}; its file holds {label}: {why}") for study, (label, why) in reviewed.items() if study in studies.index]
    if wrong := [level for _, level, _ in by_hand if level not in LEVELS]:
        raise ValueError(f"a finding made by hand is a {wrong[0]}: expected one of {', '.join(LEVELS)}")
    found += [finding for finding in by_hand if finding[0] in studies.index]
    signed = studies[studies["cross"]].assign(trait=lambda t: [reviewed[study][0] if study in reviewed else trait for study, trait in t["trait"].items()])
    faulty = {study for study, level, _ in by_hand if level == "fail"}  # a study that a person found at fault: its partners are not blamed for disagreeing with it
    if compared is not None:
        found += signs(studies, compared, traits or {}, names or {}, faulty, rules["signs"])
    found += [(study, level, text) for study, findings in labels.check(signed[["trait", "h2_z"]].to_dict("index"), tables, rules["label"], faulty).items() for level, text in findings]
    notes = {study: [] for study in studies.index}
    for study, level, text in sorted(found, key=lambda finding: LEVELS.index(finding[1])):  # a stable sort: within a level, in the order found
        notes[study].append([level, text])
    worst = lambda findings: next((level for level, _ in findings if level != "note"), "pass")
    return pd.DataFrame({"qc": [worst(notes[study]) for study in studies.index], "notes": [notes[study] for study in studies.index]}, index=studies.index)
