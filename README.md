# GWAS Sumstats

<https://gwassumstats.github.io> is a catalog of the Price lab's collection of GWAS summary statistics: 602 studies
in LDSC format, with their heritability and genetic correlations. The files stay on the O2 cluster; lab members log in
with the lab's passphrase to see where they are. It costs $0.

## How it works

1. **Release.** The collection is published on O2 as a release folder: the inventory, the S-LDSC heritabilities, the
   cross-trait LDSC genetic correlations, the nonredundant panels and the index of the per-SNP files.
2. **What the repository adds** (`data/`): the paper of each study and its short name, both drafted by a script and
   kept as small tables that are read and corrected by hand; and the numbers of each file, measured SNP by SNP on
   Slurm, from which the quality checks judge the study: pass, warn or fail.
3. **Build** (`build/`, run on O2): `python -m build` reads the release folder and `data/`, and writes `site/data/`,
   which is committed. Nothing is built on GitHub, and the site is never edited by hand to change a number.
4. **Site.** Plain HTML, CSS and JavaScript in `site/`, for desktop screens. It loads nothing from other sites. A push
   that touches `site/` deploys it (`.github/workflows/deploy.yml`).

| File in `site/data/` | Content | Read by |
|---|---|---|
| `studies.json` | The inventory with h², each study's short name, paper, quality checks, cases and controls, and whether it is the study to use of its trait; the release's counts | Anyone |
| `papers.json` | Every paper cited: authors, title, journal, DOI, BibTeX entry | Anyone, when first needed |
| `candidates.json` | The candidates for the next release, and the backtest of the rule that finds them | Anyone, when first needed |
| `rg_EUR.json`, `rg_EAS.json` | Every genetic correlation of an ancestry; the order of its heatmaps | Anyone, when first needed |
| `paths.enc.json` | The release folder and each study's file on O2, encrypted | Whoever has the passphrase |

| File in `data/` | Content | Written by |
|---|---|---|
| `references.tsv` | A reference of the inventory → PMID, DOI, the day the paper was published, and how it was matched | `python -m papers`, then by hand |
| `papers.json` | PubMed's record of every paper cited, with the day it first appeared, online or in print | `python -m papers` |
| `consortia.tsv` | A group that PubMed names as author → the short name it is cited by | By hand |
| `names.tsv` | A study's ID → its short name, and how the name was made | `python -m names`, then by hand |
| `abbreviations.tsv` | A phrase of a trait name → what stands for it in a short name | By hand |
| `traits.tsv` | A study's ID → its complex trait, the EFO term that fits, a broad category, and how they were found | `python -m traits`, then by hand |
| `trait_aliases.tsv`, `categories.tsv` | Traits named by hand, with the label of their EFO term; the rule that gives a trait its category | By hand |
| `metrics.tsv` | What was measured of each file: rows, SNPs listed twice, HapMap3 SNPs, N, Z, mean χ² | `python -m qc.scan` |
| `signs.tsv` | The sign test: pairs of signed studies, the SNPs strong in both, and those of the same sign | `python -m qc.signs` |
| `metrics_corrected.tsv`, `signs_corrected.tsv` | The same of the corrected copies: each copy as measured, and the sign test with the copies in place of the release's files | The same, with `--copies` |
| `checks.toml` | The thresholds of the quality checks | By hand |
| `labels.tsv` | A study whose trait name in the release a person found wrong → what its file holds, and how that is known | By hand |
| `findings.tsv` | What else a person found of a study, shown and counted like a finding of the checks; of the release's file, or of the study whatever its file | By hand |
| `corrections.tsv` | What is corrected of a study, and why: a source to take it from again, a sign to turn, rows to drop, an N, a label | By hand |
| `corrected_h2.tsv`, `corrected_rg.tsv` | What LDSC finds on the corrected copies: h², and rg of each with every other study of its ancestry | `python -m corrections.ldsc` |
| `candidates.toml` | The rule that finds candidates in the GWAS Catalog, and the window of the last release | By hand |
| `candidates.json` | The candidates, the papers set aside, and the rule's backtest | `python -m candidates` |
| `catalog.tsv` | The GWAS Catalog's record of each study it knows: paper, N, ancestry, sex, cohorts | `python -m candidates` |
| `samples.tsv` | Every sample the GWAS Catalog records for a paper cited, and every record that `cases.tsv` names: its individuals, by ancestry, and its cases and controls | `python -m candidates` |
| `cases.tsv` | A study's ID → whether its trait is binary or quantitative, its cases and controls, the record they are from, and how they were found | `python -m cases`, then by hand |

**Quality checks** (`qc/checks.py`). Every study goes through the same checks: its file as measured, its LDSC
results, its label against its genetic correlations (`qc/labels.py`), its row against the GWAS Catalog's record, and
the N of its file against the samples the Catalog records for its paper. The release counts individuals, cases and
controls together: an N that is an effective N, or half of one, or that of a sample of several ancestries, is named.
The sign test (`qc/signs.py`) compares every pair of signed studies at the SNPs with |Z| > 4 in both: two studies of one
trait and ancestry have the same sign at nearly all of them (99% or more in 302 of 309 pairs of the 2026 release), and
a study whose Z does not go with its alleles agrees with every other at about half. The scan also counts the SNPs
whose alleles are those of the HapMap3 list in its order, or the other way round: a file that has them the other way
round throughout had its alleles set from the list, without regard to the effect allele.
A fail means do not use the study as it is; a warn, use it with care; a note is worth knowing, and the study passes.
The site shows each finding on the study's page. A corrected study is judged twice: as released, and as corrected,
which is the verdict the site's tables show: its copy as measured, with the copy's h² and genetic correlations and the
sign test of the copies; a study with a corrected label, under the trait its file holds, where it warns for the name
it keeps. What a person found of the release's file is no finding of the corrected study, unless
`data/findings.tsv` says it holds of the study whatever its file. The label check flags a study whose label its genetic correlations
contradict; it does not relabel it. A person who has looked adds the study to `data/labels.tsv`.

**Candidates** (`candidates/`). GWAS with public summary statistics that the GWAS Catalog added after the last
release's window, and that pass a rule read off that release: one ancestry, European or East Asian; at least 50,000
individuals; a paper of at most 100 studies, in a journal that StatGenPapers reads. The backtest applies the rule to
the window of the last release and counts what it finds of the studies that the release added: 202 of 205 for 2026.
The site lists the candidates, the papers set aside and the backtest; a form on GitHub takes nominations.

**Corrections** (`corrections/`). The release is never changed. A study that the review found at fault gets a
corrected copy outside this repository, in the folder named by `SUMSTATS_CORRECTIONS`: `python -m corrections` reads
the release's file and `data/corrections.tsv`, and writes the copy with a record of its source, the checksums of both
and what changed. That folder may be purged: a run writes everything again.

**The overlay** (`build/overlay.py`). The site shows a corrected study by its copy: the build lays the corrected
h², N, SNPs and genetic correlations over the release's and keeps the release's beside them (`was`), names what was
corrected and why (`fixed`), gives genetic correlations to a study that its source signs, and chooses the panel of
an ancestry again, by the release's rule, where one of its studies was corrected. The checks and the Release notes
stay those of the release as it is.

**Sources** (`sources/`). A study whose file cannot be mended is taken from its source again: its folder in the GWAS
Catalog, or a URL. The file is downloaded, checked against the checksum its folder states, and formatted by the rules
of the release's own scripts (`sources/format.py`): SNP, A1, A2, N and Z on the HapMap3 SNPs, A1 the effect allele.
`python -m sources ACCESSION ...` checks the rules: it makes studies of the release again from the Catalog and
compares them with the release's files, byte for byte. Nine of the 2026 additions, one of each layout of columns the
release met, were made again on 2026-11-16: every one is the release's file.

**LDSC on the corrected copies** (`corrections/ldsc.py`, `corrections/panel.py`). `python -m corrections.ldsc` runs
S-LDSC and LDSC on each corrected copy with the release's options and reference files, against the file of every
other study of its ancestry, and writes what the logs say to `data/corrected_h2.tsv` and `data/corrected_rg.tsv`.
`corrections/panel.py` is the release's rule for its panels of non-redundant traits, read off its script: from the
release's own tables it gives the release's two panels, study for study.

**Cases, controls and effective N** (`cases/`, `build/recommend.py`). A study of a binary trait has the power of
its cases and controls, not of their sum: its effective N is 4 / (1/cases + 1/controls). `python -m cases` drafts
`data/cases.tsv` from the records of each study's sample: the GWAS Catalog's record of the study, or of the source
its corrected copy was taken from; else the record of its paper whose sample is the study's by its N; else, for a
study of the Pan-UK Biobank, the phenotype of its N in that project's manifest. The rest is by hand, each row with
the record or the paper its numbers are from. A binary study whose cases no public record gives has no effective N.
Where the record is of a larger sample than the file holds, the file is taken to have its share of cases.

**The study to use of a trait** (`build/recommend.py`). Of the studies of one trait and one ancestry that pass the
checks, as the site shows them, the one of the largest effective N is recommended; h² z settles a tie. If the
effective N of one of them is not known, effective N cannot rank it, and its N would overstate few cases among many
controls: h² z, the power that LDSC measured, decides between the study of the largest effective N and each study
without one, then N. An older study of the same trait, ancestry and sex is superseded by the recommended one; a study
of one sex, or a newer one, is another analysis and is left unmarked. A trait of which nothing passes has no
recommended study.

**Short names.** A study is called `FirstAuthor_Year_Trait` (`Speliotes_2010_BMI`) everywhere on the site: who and
which year its paper is cited by, and its trait, abbreviated. The name is a label only. A study's ID in the release
stays its key, names its file, and shows on its page and in the exports.

**Login.** `paths.enc.json` is encrypted with AES-256-GCM under a key derived from the passphrase (PBKDF2-HMAC-SHA256,
1,000,000 rounds, a random salt). The browser derives the key, decrypts the paths, and keeps the key in `localStorage`
until logout. The passphrase is short: this keeps the paths from the public, not from someone set on guessing it. The
files themselves are protected by O2.

## Common tasks

1. **Rebuild**, on O2, after the release changed or for a new one:

   ```bash
   export SUMSTATS_RELEASE=<release folder>
   uv run python -m papers    # the papers of new references, from PubMed: data/references.tsv, data/papers.json
   uv run python -m names     # short names for new studies: data/names.tsv
   uv run python -m traits    # the trait, EFO term and category of new studies: data/traits.tsv
   uv run python -m cases     # binary or quantitative, cases and controls of new studies: data/cases.tsv
   SUMSTATS_PASSPHRASE=... SUMSTATS_CORRECTIONS=<the folder of the corrected copies> uv run python -m build
   ```

   Without `SUMSTATS_PASSPHRASE` the build asks for it. `SUMSTATS_CORRECTIONS` gives a corrected study the path of its
   copy, behind the login; without it such a study keeps the path of the release's file. The folder's path shows on the site's Release notes after
   login. Then run the checks and commit `data/` and `site/data/`. An unchanged release rewrites nothing, and readers
   stay logged in.
2. **Review papers and names.** A row of `data/references.tsv`, `data/names.tsv`, `data/traits.tsv` or
   `data/cases.tsv` whose `how` says "check" was a guess, or could not be drafted. Correct the row, write "by hand" and why in `how`, and rebuild: the scripts keep the rows that are there.
   A study whose reference still says "check" shows no paper; a study without a short name stops the build.
3. **Measure the files**, on O2, when files were added or changed. The scan runs on Slurm (12 tasks, a few minutes)
   and writes `data/metrics.tsv`; the sign test (one task, a quarter of an hour) writes `data/signs.tsv`; then rebuild.

   ```bash
   uv run --group qc python -m qc.scan --hm3 EUR='<weights of EUR, one file per chromosome>' \
       --hm3 EAS='<weights of EAS>' --alleles <the HapMap3 SNPs with their alleles> --logs <a folder outside the repository>
   uv run --group qc python -m qc.signs --alleles <the same list> --logs <the same folder>
   ```

   The weights are those LDSC regresses with (`weights.hm3_noMHC.*.l2.ldscore.gz`). The logs hold paths: never commit
   them. Cancel a Slurm job by its ID only.
4. **Review a study the checks flag.** Write what you find, with its evidence, in `data/findings.tsv`: it shows on
   the study's page with what the checks found. If its label in the release is wrong and the right one is known, add
   it to `data/labels.tsv`: the label check then reads the corrected label, and the study warns with the evidence
   until the release is corrected. Then rebuild.
5. **Write the corrected copies**, on O2 (a Slurm task for each study taken from its source, one for the rest; minutes each):

   ```bash
   SUMSTATS_CORRECTIONS=<a folder outside the repository> uv run --group qc python -m corrections \
       --snplist <the HapMap3 SNPs with their alleles> --bim-prefix <the reference panel's .bim files>
   ```

   The two files are those the release formatted with; they are needed for the studies taken from their source.
   Then LDSC on the copies, as the release ran it (a Slurm task for each; the settings file, kept outside the
   repository, names LDSC and its reference files: see `corrections/ldsc.py`):

   ```bash
   SUMSTATS_CORRECTIONS=<the same folder> uv run --group qc python -m corrections.ldsc --settings <file>
   ```

   Then measure the copies and run the sign test with them, so that each corrected study is judged as corrected:
   the commands of step 3 with `--copies <the same folder>` write `data/metrics_corrected.tsv` and
   `data/signs_corrected.tsv`. Then rebuild.

6. **Look for candidates**, from any machine: `uv run python -m candidates --refresh` downloads the GWAS Catalog's
   tables (to `.cache/`, not committed) and writes `data/candidates.json`, `data/catalog.tsv` and
   `data/samples.tsv`; then rebuild.
7. **A new release:** set its window in `data/candidates.toml`, run the steps above, then update the Release notes
   in `site/index.html`: the tests say which written-in numbers no longer match.
8. **Change the passphrase:** rebuild with the new one. Every reader logs in again.

## Checks

```bash
export SUMSTATS_RELEASE=<release folder> SUMSTATS_PASSPHRASE=... SUMSTATS_CORRECTIONS=<folder of the copies>   # on O2; never in a file of the repository
uv run pytest                                                      # the build; the data against the release; nothing private committed
uv run --with playwright python dev/site_check.py                  # the site in a headless browser, with login
python3 -m http.server -d site 8000                                # preview
```

Without the variables, the tests that need them are skipped. Run with all three before a push.

## Setup

1. **GitHub:** Settings → Pages → Source: GitHub Actions. The workflow asks for its own write permission.
2. **Browser check:** once, `uv run --with playwright python -m playwright install chromium`.

## Layout

| Path | What |
|---|---|
| `build/` | The build: `release.py` reads a release folder, `merge.py` adds what `data/` holds, `crypt.py` encrypts the paths |
| `papers/` | Finds each study's paper in PubMed (`python -m papers`); `cite.py`: how a paper is cited, its BibTeX entry |
| `names/` | Drafts the short names (`python -m names`) |
| `traits/` | Drafts each study's complex trait, EFO term and category (`python -m traits`) |
| `qc/` | Quality checks: `scan.py` measures the files on Slurm, `checks.py` judges each study, `labels.py` checks its label |
| `corrections/` | Writes corrected copies of the studies of `data/corrections.tsv`, outside the repository, with provenance |
| `candidates/` | Finds candidates for the next release in the GWAS Catalog, and backtests the rule (`python -m candidates`) |
| `data/` | What the repository adds to a release: papers and short names, reviewable |
| `site/` | The website; `site/data/` is written by the build |
| `dev/` | Local only: `site_check.py` |
| `tests/` | The build on a made-up release; the committed data against the real one; what must stay private |
