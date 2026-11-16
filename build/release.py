"""Reads a release folder of the collection into what the site shows.

A release folder holds the inventory (one row per study), the S-LDSC heritabilities and the nonredundant panel of
each ancestry, the cross-trait LDSC results of every pair of signed studies of an ancestry, and under _audit/ the
papers of the studies it adds and the index of the per-SNP files. Its file names end in the release year.
"""
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd

# site name -> column of the release's pair tables
RG_COLUMNS = {"rg": "genetic_correlation", "se": "genetic_correlation_se", "z": "genetic_correlation_z",
              "p": "genetic_correlation_p", "gcov": "genetic_covariance", "gcov_se": "genetic_covariance_se"}
# a study's public fields; the path of its file (h2_file) and of the copy rg was estimated on (rg_file) are not among them
PUBLIC = ["trait", "anc", "n", "m", "origin", "ref", "signed", "alleles", "cross", "release", "h2", "h2_se", "h2_z", "intercept",
          "panel", "m_rg", "pmid", "paper", "catalog", "paper_id"]


def one(folder: Path, pattern: str) -> Path:
    """The only file in `folder` that matches `pattern`."""
    found = sorted(folder.glob(pattern))
    if len(found) != 1:
        raise FileNotFoundError(f"expected one file {pattern} in the release folder, found {len(found)}")
    return found[0]


def table(path: Path) -> pd.DataFrame:
    """A tab-delimited file with a header row: every cell as text, without spaces around it."""
    cells = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, quoting=csv.QUOTE_NONE, encoding="utf-8")
    return cells.apply(lambda column: column.str.strip())


def per_ancestry(release: Path, pattern: str, ancestries) -> pd.DataFrame:
    """The release's files of one kind, one per ancestry, stacked; `anc` is the ancestry of a row's file."""
    return pd.concat([table(one(release, pattern.format(anc))).assign(anc=anc) for anc in ancestries], ignore_index=True)


def by(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    """`frame` with the rows named by a column, which must name each once."""
    if frame[column].duplicated().any():
        raise ValueError(f"{column}: {frame.loc[frame[column].duplicated(), column].iloc[0]} occurs twice")
    return frame.set_index(column)


def flag(column: pd.Series) -> pd.Series:
    if not column.isin(["YES", "NO"]).all():
        raise ValueError(f"{column.name}: expected YES or NO")
    return column == "YES"


def studies(release: Path) -> pd.DataFrame:
    """One row per study, by ID, in the inventory's order: the inventory, the study's h², its rank in the
    nonredundant panel of its ancestry, the release it came with, its paper if the release names one, and the paths
    of its file and of the copy that rg was estimated on (`h2_file`, `rg_file`: never public)."""
    inventory = table(one(release, "summary_statistics_*.txt"))
    d = pd.DataFrame({
        "id": inventory["Trait Identifier"], "trait": inventory["Trait Name"], "anc": inventory["Population"],
        "n": inventory["N"].astype(int), "m": inventory["M"].astype(int), "origin": inventory["Origin"], "ref": inventory["Reference"],
        "signed": flag(inventory["Signed SumStat"]), "alleles": flag(inventory["Allele Info"]), "cross": flag(inventory["Cross-Trait LDSC"]),
    }).pipe(by, "id")
    if list(table(one(release, "signed_summary_statistics_*.txt"))["Trait Identifier"]) != list(d.index[d["signed"]]):
        raise ValueError("the signed inventory is not the inventory's signed studies")
    ancestries = d["anc"].unique()

    h2 = per_ancestry(release, "sldsc_heritabilities_{}_*.txt", ancestries).pipe(by, "trait_identifier")
    panel = per_ancestry(release, "non_redundent_traits_{}_*.txt", ancestries).pipe(by, "trait_identifier")
    if set(h2.index) != set(d.index) or not (h2["anc"] == d["anc"].reindex(h2.index)).all():
        raise ValueError("the heritability files do not cover the inventory, ancestry by ancestry")
    if not (panel["anc"] == d["anc"].reindex(panel.index)).all():
        raise ValueError("a nonredundant panel lists a study of another ancestry, or an unknown one")
    d = d.join(h2[["h2", "h2_se", "h2_z", "intercept"]].astype(float))
    d["panel"] = (panel.groupby("anc").cumcount() + 1).astype("Int64")  # the panel files are in selection order

    files = table(release / "_audit" / "sumstats_file_index.tsv").pipe(by, "Trait Identifier")
    if list(files.index) != list(d.index):
        raise ValueError("the file index does not list the inventory's studies")
    d["release"] = files["Release"].str.extract(r"(\d{4})", expand=False).astype(int)  # "2024", "2026 addition"
    d["m_rg"] = pd.to_numeric(files["M in RG input"].where(files["M in RG input"] != "")).astype("Int64")
    d["h2_file"], d["rg_file"] = files["H2 Sumstats"], files["RG Sumstats"].where(files["RG Sumstats"] != "")
    if d["h2_file"].eq("").any() or (d["rg_file"].notna() != d["cross"]).any():
        raise ValueError("the file index lacks an h² file, or its rg files are not those of the cross-trait LDSC studies")

    papers = table(release / "_audit" / "new_study_provenance.tsv").pipe(by, "Trait Identifier")
    if not papers.index.isin(d.index).all():
        raise ValueError("the provenance table names a study that the inventory lacks")
    d["pmid"], d["paper"], d["catalog"], d["paper_id"] = (papers[c] for c in ("PMID", "Paper Link", "Catalog Study Link", "Target Paper Identifier"))
    if not pd.concat([d["paper"], d["catalog"]]).dropna().str.fullmatch(r"https://\S+").all():
        raise ValueError("a paper or catalog link is not an https address")
    return d


def correlations(release: Path, d: pd.DataFrame, anc: str) -> dict:
    """The cross-trait LDSC results of an ancestry: `ids`, its studies in the inventory's order, and one list per
    column over the pairs i < j of `ids`, row by row: (0,1), (0,2), ..., (1,2), ... A pair without an estimate
    holds null; `notes` gives LDSC's reason where the release keeps one, by the pair's position. `heat` orders the
    studies, by their place in `ids`, for a heatmap: hierarchical clustering on the distance 1 - |rg| (see tree)."""
    ids = d.index[(d["anc"] == anc) & d["cross"]]
    pairs = table(one(release, f"ldsc_genetic_correlations_trait_pairs_{anc}_*.txt"))
    n = len(ids)
    i, j = np.sort([ids.get_indexer_for(pairs["trait1"]), ids.get_indexer_for(pairs["trait2"])], axis=0)
    position = i * (2 * n - i - 1) // 2 + j - i - 1
    if i.min() < 0 or (i == j).any() or len(position) != n * (n - 1) // 2 or len(np.unique(position)) != len(position):
        raise ValueError(f"the {anc} pair table is not every pair of the {n} {anc} cross-trait LDSC studies, once")
    out = {"ids": list(ids)}
    for name, column in RG_COLUMNS.items():
        values = np.full(len(position), np.nan)
        values[position] = pd.to_numeric(pairs[column], errors="coerce").to_numpy(dtype=float, na_value=np.nan)  # "nan" and the like
        out[name] = np.where(np.isnan(values), None, values).tolist()
    reason = pairs["genetic_correlation_z"].str.extract(r"\((\w+)\s+out of bounds\)", expand=False)  # "nan (nan) (rg out of bounds)"
    out["notes"] = {str(k): f"{r} out of bounds" for k, r in sorted(zip(position[reason.notna().to_numpy()].tolist(), reason.dropna()))}
    # The order of the heatmaps: of every study of the table, and of the nonredundant panel, each clustered on its own
    rg = np.zeros((n, n))
    rg[np.triu_indices(n, 1)] = np.nan_to_num(np.array(out["rg"], dtype=float))  # a pair without an estimate counts as rg = 0
    distance = 1 - np.minimum(np.abs(rg + rg.T), 1)  # an estimate beyond ±1 counts as ±1
    panel = np.flatnonzero(d.loc[ids, "panel"].notna().to_numpy())
    out["heat"] = {"all": tree(distance), "panel": tree(distance[np.ix_(panel, panel)], panel.tolist())}
    return out


def tree(distance: np.ndarray, names=None) -> dict:
    """Average-linkage clustering of a symmetric matrix of distances. `order`: the leaves (as `names`, by default their
    row numbers) from one end of the tree to the other. `merges`: first to last, [node, node, height], where a node
    is a leaf's place in `order` or, from the number of leaves on, an earlier merge."""
    n = len(distance)
    d = np.array(distance, dtype=float)
    np.fill_diagonal(d, np.inf)
    size, node, leaves, merges = np.ones(n), list(range(n)), [[i] for i in range(n)], []
    for step in range(n - 1):
        a, b = divmod(int(d.argmin()), n)  # the two closest clusters; of equally close pairs, the first
        merges.append((node[a], node[b], float(d[a, b])))
        d[a, :] = d[:, a] = (size[a] * d[a] + size[b] * d[b]) / (size[a] + size[b])  # the mean distance between members
        d[b, :] = d[:, b] = d[a, a] = np.inf
        size[a], node[a], leaves[a] = size[a] + size[b], n + step, leaves[a] + leaves[b]
    order = leaves[a] if n > 1 else list(range(n))
    place = {leaf: k for k, leaf in enumerate(order)}
    return {"order": [order[k] if names is None else names[order[k]] for k in range(n)],
            "merges": [[place.get(i, i), place.get(j, j), round(height, 4)] for i, j, height in merges]}


def summary(release: Path, d: pd.DataFrame, rg: dict) -> dict:
    """The counts the release notes quote: of the whole release (`all`) and of each ancestry, now and before."""
    year = int(re.search(r"_(\d{4})\.txt$", one(release, "summary_statistics_*.txt").name)[1])
    packaged = re.search(r"_(\d{4})(\d{2})(\d{2})$", release.name)
    new = d["release"] == year

    def count(rows):
        return {"studies": len(rows), "before": int((~new[rows.index]).sum()), "added": int(new[rows.index].sum()),
                "signed": int(rows["cross"].sum()), "signed_before": int((rows["cross"] & ~new[rows.index]).sum())}

    out = {"year": year, "packaged": "-".join(packaged.groups()) if packaged else None,
           "papers": int(d.loc[new, "paper_id"].nunique()), "ancestries": list(rg),
           "all": {**count(d), "unsigned": int((~d["cross"]).sum()), "unsigned_added": int((~d["cross"] & new).sum()),
                   "no_alleles": int((~d["alleles"]).sum())}}
    for anc, t in rg.items():
        rows, n = d[d["anc"] == anc], len(t["ids"])
        failed = np.isnan(np.array(t["rg"], dtype=float))
        i, j = np.triu_indices(n, 1)  # the pairs, in the order of the lists
        per_study = np.bincount(np.concatenate([i[failed], j[failed]]), minlength=n)
        panel = rows.index[rows["panel"].notna()]
        out[anc] = count(rows) | {"pairs": len(failed), "failed": int(failed.sum()),
                                  "worst": {"id": t["ids"][per_study.argmax()], "failed": int(per_study.max()), "pairs": n - 1},
                                  "panel": len(panel), "panel_added": int(new[panel].sum())}
        out[anc]["pairs_before"] = out[anc]["signed_before"] * (out[anc]["signed_before"] - 1) // 2
        previous = sorted((release / "_audit" / "inputs").glob(f"non_redundent_traits_{anc}_*.txt"))
        if len(previous) == 1:  # the panel of the release before, kept with the release's inputs
            before = table(previous[0])["trait_identifier"]
            out[anc] |= {"panel_before": len(before), "panel_kept": int(panel.isin(before).sum())}
    return out


def private(release: Path, d: pd.DataFrame) -> dict:
    """What only a logged-in reader sees: the release folder with its tables, and where each study's file is, as
    [folder number, file name]. It is the file h² was estimated on; the copy without strand-ambiguous SNPs that rg
    was estimated on, where the release keeps one, is not listed."""
    folder, _, name = (d["h2_file"].str.rpartition("/")[part] for part in range(3))
    folders = pd.unique(folder).tolist()
    return {"release": str(release), "tables": sorted(p.name for p in release.iterdir() if p.is_file()), "folders": folders,
            "files": dict(zip(d.index, zip(folder.map(folders.index).tolist(), name)))}
