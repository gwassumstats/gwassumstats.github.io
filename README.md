# GWAS Sumstats

<https://gwassumstats.github.io> is a catalog of the Price lab's collection of GWAS summary statistics: 602 studies
in LDSC format, with their heritability and genetic correlations. The files stay on the O2 cluster; lab members log in
with the lab's passphrase to see where they are. It costs $0.

## How it works

1. **Release.** The collection is published on O2 as a release folder: the inventory, the S-LDSC heritabilities, the
   cross-trait LDSC genetic correlations, the nonredundant panels and the index of the per-SNP files.
2. **Build** (`build/`, run on O2): `python -m build` reads the release folder and writes `site/data/`, which is
   committed. Nothing is built on GitHub, and the site is never edited by hand to change a number.
3. **Site.** Plain HTML, CSS and JavaScript in `site/`, for desktop screens. It loads nothing from other sites. A push
   that touches `site/` deploys it (`.github/workflows/deploy.yml`).

| File in `site/data/` | Content | Read by |
|---|---|---|
| `studies.json` | The inventory with h², and the release's counts | Anyone |
| `rg_EUR.json`, `rg_EAS.json` | Every genetic correlation of an ancestry; the order of its heatmaps | Anyone, when first needed |
| `paths.enc.json` | The release folder and each study's file on O2, encrypted | Whoever has the passphrase |

**Login.** `paths.enc.json` is encrypted with AES-256-GCM under a key derived from the passphrase (PBKDF2-HMAC-SHA256,
1,000,000 rounds, a random salt). The browser derives the key, decrypts the paths, and keeps the key in `localStorage`
until logout. The passphrase is short: this keeps the paths from the public, not from someone set on guessing it. The
files themselves are protected by O2.

## Common tasks

1. **Rebuild**, on O2, after the release changed or for a new one:

   ```bash
   SUMSTATS_PASSPHRASE=... uv run python -m build <release folder>
   ```

   Without `SUMSTATS_PASSPHRASE` the build asks for it; `SUMSTATS_RELEASE` can name the folder instead. The folder's
   path shows on the site's Release notes after login. Then run the checks and commit `site/data/`. An unchanged
   release rewrites nothing, and readers stay logged in.
2. **A new release:** rebuild, then update the Release notes in `site/index.html`: the tests say which written-in
   numbers no longer match.
3. **Change the passphrase:** rebuild with the new one. Every reader logs in again.

## Checks

```bash
export SUMSTATS_RELEASE=<release folder> SUMSTATS_PASSPHRASE=...   # on O2; never in a file of the repository
uv run pytest                                                      # the build; the data against the release; nothing private committed
uv run --with playwright python dev/site_check.py                  # the site in a headless browser, with login
python3 -m http.server -d site 8000                                # preview
```

Without the two variables, the tests that need them are skipped. Run with both before a push.

## Setup

1. **GitHub:** Settings → Pages → Source: GitHub Actions. The workflow asks for its own write permission.
2. **Browser check:** once, `uv run --with playwright python -m playwright install chromium`.

## Layout

| Path | What |
|---|---|
| `build/` | The build: `release.py` reads a release folder, `crypt.py` encrypts the paths |
| `site/` | The website; `site/data/` is written by the build |
| `dev/` | Local only: `site_check.py` |
| `tests/` | The build on a made-up release; the committed data against the real one; what must stay private |
