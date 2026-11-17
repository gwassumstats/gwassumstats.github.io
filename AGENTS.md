# GWAS Sumstats

A catalog of the Price lab's GWAS summary statistics, at <https://gwassumstats.github.io> (GitHub org `gwassumstats`).
Read `README.md` before changing anything. Sibling of StatGenPapers (`../statgenpapers.github.io`): follow its
conventions unless this file says otherwise.

Status (2026-11-16): of `notebook/plans/2026-11-16-portal.md`, phase 1 is built (the catalog: Sumstats,
Heritability, Genetic correlation, Release notes, About, a page per study, the login that shows the O2 paths), and of
phases 2 and 3 the papers, the short names, the cart (ticked studies export as TSV and BibTeX), the quality checks
with the label check and the link to report a problem, and the candidates for the next release with the backtest of
their rule. Open: fixes to the studies that the checks flag wait for the owner (the review is in the notebook).

## Hard constraints

1. Never publish, in a committed file, a commit message or a page:
   1. per-SNP data or raw files;
   2. the passphrase;
   3. an O2 path, except encrypted in `site/data/paths.enc.json`;
   4. an O2 user name or home directory;
   5. content of the lab's mailing list: emails, or who flagged what.
2. The passphrase and the release folder's path come from the environment (`SUMSTATS_PASSPHRASE`, `SUMSTATS_RELEASE`)
   or a prompt. Never write them into a file of the repository, a test or a page.
3. `site/data/` is written by `python -m build` only. The files stay on O2: the repository holds metadata, h² and rg.
   `data/` is drafted by scripts and corrected by hand; a script never overwrites a row that is there.
4. The site loads nothing from other sites (no analytics, fonts or scripts): the passphrase is typed into it.
5. Before a push, run `uv run pytest` with both variables set, and `dev/site_check.py`. A skipped test did not run.
   What a Slurm job writes besides `data/metrics.tsv` (its logs hold paths) stays outside the repository, under
   `notebook/results/`. Cancel a Slurm job by its explicit ID only.
6. Push workflow files over SSH (the `gh` token lacks the `workflow` scope).

## Owner's preferences

1. Desktop only: no phone layout, no phone checks.
2. Short, plain text everywhere: About is an FAQ; docs stay brief.
3. Small UI details (size, spacing, motion): pick the best-practice option yourself; show a screenshot of visible changes.
4. h² and rg of every study are public, restricted sources (23andMe) included. Only paths are encrypted.
5. Committed site data stays well under 100 MB: compact files, loaded when needed, rather than dropped content.
6. An estimate shows as "estimate ± SE", both to two significant digits, in one column that sorts by the estimate; z to
   two significant digits too. No bars in tables. rg is red above 0, blue below.
7. A study is called by its short name everywhere: tables, pickers, heatmaps, tooltips, titles, notes. Its ID shows
   only on its page and in exports; its trait name in the release only on its page, as "Release name". The search
   finds all three.

## Terminology

- *Sumstats:* GWAS summary statistics files.
- *Study:* one GWAS summary-statistics file of the collection, keyed by its Trait Identifier (its ID). Never "dataset".
- *Short name:* what the site calls a study: `FirstAuthor_Year_Trait` (`Speliotes_2010_BMI`), in `data/names.tsv`.
  A label only: this repository never renames, moves or writes a file on O2, nor changes an ID. A future release may
  adopt short names as file names for new studies; that is the release owner's step.
- *Paper, reference:* the paper a study cites (`data/references.tsv`, PubMed's record in `data/papers.json`); its
  *label* is who and which year it is cited by, and its journal: `Speliotes et al. 2010 Nat Genet`. A paper's day
  of publication is the earliest of PubMed's online and print dates, as precise as PubMed gives it (`2010-10-10`,
  or `2018-07`); the Reference column and the BibTeX export sort by it.
- *Cart:* the studies a reader ticked, kept in the browser; it exports their names, IDs and paths (TSV) and their
  papers (BibTeX).
- *QC, checks:* what `qc/checks.py` says of a study: *pass*, *warn* or *fail*, with its *findings*, each a fail, a
  warn or a *note* (worth knowing; the study passes). Thresholds: `data/checks.toml`.
- *Label check:* the part of QC that compares a study's label (its trait name in the release) with its genetic
  correlations: a *twin* is a study with rg of 0.95 or more; a label is *confirmed* when a study of the same label
  agrees. It flags, with the best-matching studies as evidence; it never relabels.
- *Trait:* the complex trait a study is of (`data/traits.tsv`): one per study, shared by studies of the trait
  whatever their paper, ancestry, sex or adjustment; with the EFO term that fits, if one does, and a broad *category*.
- *Candidate:* a GWAS of the GWAS Catalog with public sumstats, added after the last release's *window* (the days
  of publication it took: 2024/03/01 to 2026/08/26 for 2026), that passes the *rule* of `data/candidates.toml`.
- *Backtest:* the rule applied to the window of the last release, against the studies that release added: what it
  finds (recall), misses, and finds beyond.
- *Reviewed label:* a row of `data/labels.tsv`: what a study's file holds where the release's trait name is wrong,
  found by a person. Label fixes beyond that wait for the owner's approval.
- *Release:* the collection as published on O2 in one *release folder*; the site shows one (2026: 602 studies).
- *Origin:* the update that added a study (`Finucane2015`, ..., `New 2026`); "Release" on a study's page.
- *Signed:* a study whose Z has a direction. Only signed studies have rg.
- *File, path:* a study's per-SNP file, which h² was estimated on, and where it is on O2: one path per study. For a
  study added in 2026, rg was estimated on a copy without strand-ambiguous SNPs; the site does not list the copies.
- *Pair:* two signed studies of one ancestry. Its rg can be *not estimated* (LDSC gave none).
- *Panel:* the nonredundant panel of an ancestry (85 EUR, 17 EAS): h² z > 6, and rg² of at most 0.1 between any two.
- *Strongest:* a study's genetic correlations with |z| ≥ 3, by |rg|; a study's page lists ten.
- *Heatmap order:* hierarchical clustering, average linkage on 1 − |rg|, computed by the build.
- *Passphrase:* the lab's shared login. *Blob:* `site/data/paths.enc.json`, what it decrypts: the paths, nothing else.
