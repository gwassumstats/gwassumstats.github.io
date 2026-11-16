# GWAS Sumstats

A catalog of the Price lab's GWAS summary statistics, at <https://gwassumstats.github.io> (GitHub org `gwassumstats`).
Read `README.md` before changing anything. Sibling of StatGenPapers (`../statgenpapers.github.io`): follow its
conventions unless this file says otherwise.

Status (2026-11-16): phase 1 of `notebook/plans/2026-11-16-portal.md` is built, the catalog: Sumstats, Heritability,
Genetic correlation, Release notes, About, a page per dataset, and the login that shows the O2 paths. Next: phase 2,
tools for members (cart, citations, label check); phase 3, the next release (candidates, backtest, QC).

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
4. The site loads nothing from other sites (no analytics, fonts or scripts): the passphrase is typed into it.
5. Before a push, run `uv run pytest` with both variables set, and `dev/site_check.py`. A skipped test did not run.
6. Push workflow files over SSH (the `gh` token lacks the `workflow` scope).

## Owner's preferences

1. Desktop only: no phone layout, no phone checks.
2. Short, plain text everywhere: About is an FAQ; docs stay brief.
3. Small UI details (size, spacing, motion): pick the best-practice option yourself; show a screenshot of visible changes.
4. h² and rg of every dataset are public, restricted sources (23andMe) included. Only paths are encrypted.
5. Committed site data stays well under 100 MB: compact files, loaded when needed, rather than dropped content.

## Terminology

- *Sumstats:* GWAS summary statistics files.
- *Dataset:* one GWAS of the collection: an ID, its row of the inventory, an h² file and, if signed, an rg file.
- *Release:* the collection as published on O2 in one *release folder*; the site shows one (2026: 602 datasets).
- *Origin:* the update that added a dataset (`Finucane2015`, ..., `New 2026`); the site's column "Release".
- *Signed:* a dataset whose Z has a direction. Only signed datasets have rg.
- *h² file, rg file:* a dataset's per-SNP file for S-LDSC and for cross-trait LDSC; one file up to 2024, two from 2026.
- *Pair:* two signed datasets of one ancestry. Its rg can be *not estimated* (LDSC gave none).
- *Panel:* the nonredundant panel of an ancestry (85 EUR, 17 EAS): h² z > 6, and rg² of at most 0.1 between any two.
- *Strongest:* a dataset's genetic correlations with |z| ≥ 3, by |rg|; a dataset's page lists ten.
- *Heatmap order:* hierarchical clustering, average linkage on 1 − |rg|, computed by the build.
- *Passphrase:* the lab's shared login. *Blob:* `site/data/paths.enc.json`, what it decrypts: the paths, nothing else.
