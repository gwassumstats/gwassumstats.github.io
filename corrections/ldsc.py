"""LDSC for the corrected copies, run as the release ran it: the heritability of each (S-LDSC, baselineLD v2.2) and its
genetic correlation with every other study of its ancestry.

  SUMSTATS_CORRECTIONS=<the folder of the corrected copies> uv run --group qc python -m corrections.ldsc [RELEASE] \\
      --settings FILE [--account NAME]

RELEASE defaults to $SUMSTATS_RELEASE and is only read. FILE, kept outside the repository, says where LDSC and its
reference files are, as the release's own scripts name them:

  python = "<the Python 2.7 that runs LDSC>"
  ldsc = "<ldsc.py>"
  [EUR]                    # and [EAS]
  baseline = "<prefix of the baselineLD v2.2 LD scores>"      # --ref-ld-chr, for h²
  frequencies = "<prefix of the allele frequencies>"          # --frqfile-chr, for h²
  ldscores = "<prefix of the LD scores>"                      # --ref-ld-chr, for rg
  weights = "<prefix of the regression weights>"              # --w-ld-chr, for both

Each corrected study gets a Slurm task: one S-LDSC run, and one LDSC run of its copy against the file of every other
study of its ancestry that has genetic correlations; another corrected study is taken by its corrected copy. A study
has genetic correlations if the release has them for it, or if its copy was taken from its source, where the effect
allele is known: a file that the release held unsigned can be signed in its source. The logs stay under the folder
(ldsc/h2/, ldsc/rg/). A log stands as long as the files it is of have not changed: what it holds is not run again.
So the command can be given again, under another --name, after a task failed or after more studies were corrected:
only what is missing is run, and nothing at all if every log stands. What the logs say goes to the repository, without
a path, once no study lacks LDSC:

  data/corrected_h2.tsv   a line a corrected study: h², its standard error, z and the intercept, as the release's table
  data/corrected_rg.tsv   a line a pair with a corrected study, as the release's table; a pair that LDSC could not
                          estimate has nan throughout
"""
import argparse
import datetime
import math
import os
import re
import subprocess
import tomllib
from pathlib import Path

import pandas as pd

from build import release as rel
from corrections.__main__ import planned
from papers import cite
from sources.slurm import submitted

ROOT = Path(__file__).parents[1]
H2 = ["trait_identifier", "h2", "h2_se", "h2_z", "intercept"]
RG = ["trait1", "trait2", "genetic_correlation", "genetic_correlation_se", "genetic_correlation_z", "genetic_correlation_p", "genetic_covariance", "genetic_covariance_se"]
PAIR = {"rg": r"Genetic Correlation:\s*([^\s]+)\s*\(([^)]+)\)", "z": r"Z-score:\s*([^\s]+)", "p": r"P:\s*([^\s]+)", "gcov": r"Total Observed scale gencov:\s*([^\s]+)\s*\(([^)]+)\)"}


def finite(text: str) -> bool:
    """Whether a value of a log is a finite number."""
    try:
        return math.isfinite(float(text))
    except ValueError:
        return False


def heritability(log: str) -> dict:
    """What an S-LDSC log says of a study, as the release's table has it: h², its standard error, their ratio to five
    decimals, and the intercept."""
    if not (found := re.search(r"Total Observed scale h2: ([^ ]+) \(([^)]+)\)", log)):
        raise ValueError("S-LDSC did not finish: its log has no h²")
    h2, se = found.groups()
    return {"h2": h2, "h2_se": se, "h2_z": str(round(float(h2) / float(se), 5)), "intercept": re.search(r"Intercept: ([^ ]+) \(([^)]+)\)", log)[1]}


def listed(log: str) -> list[str]:
    """The files of an LDSC log of genetic correlations, as its call names them: the study's, then the others'."""
    found = re.search(r"--rg (\S+)", log)
    return found[1].split(",") if found else []


def correlations(log: str) -> dict[str, list[str]]:
    """What an LDSC log of one study against others says of each pair, by the other's file: rg, its standard error, z,
    p, the genetic covariance and its standard error, as LDSC wrote them. A pair that LDSC could not estimate, or
    whose numbers are not all finite, has nan throughout, as in the release. A log that LDSC did not finish says
    nothing."""
    others, parts = listed(log)[1:], re.split(r"Computing rg for phenotype \d+/\d+", log)[1:]
    if "Summary of Genetic Correlation Results" not in log or len(parts) != len(others):
        return {}
    pairs = {}
    for other, part in zip(others, parts):
        found = [re.search(pattern, part) for pattern in PAIR.values()]
        values = [value for match in found if match for value in match.groups()]
        pairs[other] = values if all(found) and all(map(finite, values)) else ["nan"] * 6
    return pairs


def kept(prefix: Path, file: str) -> str:
    """The text of the LDSC log at `prefix`, if it still stands: it is there, and later than the file it is of."""
    log = Path(f"{prefix}.log")
    return log.read_text(errors="replace") if log.exists() and log.stat().st_mtime > Path(file).stat().st_mtime else ""


def standing(prefix: Path, file: str) -> dict[str, list[str]]:
    """The pairs of a study's log of genetic correlations that still stand, by the other study's file: the log is
    whole, of the study's file as it is now, and the other's file has not changed since."""
    text = kept(prefix, file)
    if listed(text)[:1] != [file]:
        return {}
    done = Path(f"{prefix}.log").stat().st_mtime
    return {other: values for other, values in correlations(text).items() if Path(other).exists() and Path(other).stat().st_mtime < done}


def lacking(study: str, file: str, partners: dict[str, str], out: Path, covered: frozenset = frozenset()) -> bool:
    """Whether LDSC has yet to run for a study: its S-LDSC log does not stand, or its log of genetic correlations
    lacks the pair of a partner as the files are now. Partners that are `covered`, corrected studies, do not count:
    their own logs hold the pair."""
    pairs = standing(out / "ldsc" / "rg" / study, file)
    return "Total Observed scale h2" not in kept(out / "ldsc" / "h2" / study, file) or any(path not in pairs and other not in covered for other, path in partners.items())


def found(study: str, file: str, partners: dict[str, str], out: Path) -> dict:
    """What a study's logs say as they stand: its h², and each partner's pair, or None for a pair that they do not
    hold as the files are now."""
    pairs = standing(out / "ldsc" / "rg" / study, file)
    return {"study": study, "h2": heritability(kept(out / "ldsc" / "h2" / study, file)), "rg": {other: pairs.get(path) for other, path in partners.items()}}


def rerun(study: str, ancestry: str, file: str, partners: dict[str, str], out: Path, settings: dict, covered: frozenset = frozenset()) -> dict:
    """One Slurm task: S-LDSC on a corrected copy, then LDSC of it against each partner's file, with the release's
    options; returns what the logs say. What a log that still stands holds is not run again."""
    run, of = [settings["python"], settings["ldsc"]], settings[ancestry]
    h2, rg = out / "ldsc" / "h2" / study, out / "ldsc" / "rg" / study
    h2.parent.mkdir(parents=True, exist_ok=True)
    rg.parent.mkdir(parents=True, exist_ok=True)
    if "Total Observed scale h2" not in kept(h2, file):
        subprocess.run([*run, "--h2", file, "--ref-ld-chr", of["baseline"], "--w-ld-chr", of["weights"], "--overlap-annot", "--print-coefficients", "--frqfile-chr", of["frequencies"], "--out", str(h2)])
    if lacking(study, file, partners, out, covered):
        subprocess.run([*run, "--rg", ",".join([file, *partners.values()]), "--ref-ld-chr", of["ldscores"], "--w-ld-chr", of["weights"], "--out", str(rg)])
    if lacking(study, file, partners, out, covered):
        raise ValueError(f"LDSC did not finish for {study}: see its logs")
    return found(study, file, partners, out)


def main():
    parser = argparse.ArgumentParser(prog="python -m corrections.ldsc", description="Runs LDSC on the corrected copies as the release ran it, and keeps what it found.")
    parser.add_argument("release", nargs="?", default=os.environ.get("SUMSTATS_RELEASE"), help="the release folder, read only (default: $SUMSTATS_RELEASE)")
    parser.add_argument("--out", default=os.environ.get("SUMSTATS_CORRECTIONS"), help="the folder of the corrected copies (default: $SUMSTATS_CORRECTIONS)")
    parser.add_argument("--settings", type=Path, required=True, help="a file outside the repository: where LDSC and its reference files are")
    parser.add_argument("--only", nargs="+", metavar="ID", help="run LDSC for these corrected studies only, where others are being run elsewhere")
    parser.add_argument("--account", help="the Slurm account to charge")
    parser.add_argument("--name", default=f"{datetime.date.today():%Y%m%d}_corrected_ldsc", help="the name of the Slurm tasks and of their logs")
    args = parser.parse_args()
    if not args.release or not args.out:
        parser.error("name the release folder and the folder of the corrected copies, or set SUMSTATS_RELEASE and SUMSTATS_CORRECTIONS")
    release, out = Path(args.release).resolve(), Path(args.out).resolve()
    settings = tomllib.loads(args.settings.read_text())
    d = rel.studies(release)
    index = pd.read_csv(release / "_audit" / "sumstats_file_index.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("Trait Identifier")
    theirs = index["RG Sumstats"].where(index["RG Sumstats"] != "", index["H2 Sumstats"])  # what the release estimated rg on
    plan = planned(cite.table(ROOT / "data" / "corrections.tsv"))
    corrected = {study: str(out / "sumstats" / f"{study}.sumstats") for study in plan}
    if missing := [study for study, copy in corrected.items() if not Path(copy).is_file()]:
        parser.error(f"no corrected copy of {missing[0]} in {out.name}/sumstats: run python -m corrections first")
    files = theirs.to_dict() | corrected
    correlated = d["cross"] | d.index.isin([study for study, fixes in plan.items() if fixes[0]["fix"] == "source"])  # the release has rg for it, or its copy is from its source
    partners = [{other: files[other] for other in d.index[(d["anc"] == d.at[study, "anc"]) & correlated] if other != study} if correlated[study] else {} for study in corrected]
    covered = {study: frozenset(corrected) - {study} for study in corrected}
    behind = lambda: [k for k, study in enumerate(corrected) if lacking(study, corrected[study], partners[k], out, covered[study])]
    if todo := [k for k in behind() if not args.only or list(corrected)[k] in args.only]:
        studies = [list(corrected)[k] for k in todo]
        submitted(args.name, out / "submitit-logs", args.account, rerun, studies, [d.at[study, "anc"] for study in studies], [corrected[study] for study in studies], [partners[k] for k in todo], [out] * len(todo),
                  [settings] * len(todo), [covered[study] for study in studies], memory=16, minutes=120, strict=False)
    if left := [list(corrected)[k] for k in behind()]:
        parser.exit(1, f"LDSC has yet to finish for {len(left)} studies ({', '.join(left[:3])}{', ...' if len(left) > 3 else ''}): see their logs, or give the command again; nothing was written\n")
    tasks = [found(study, corrected[study], partners[k], out) for k, study in enumerate(corrected)]  # from the logs as they stand, whoever ran them
    pd.DataFrame([[task["study"], *task["h2"].values()] for task in tasks], columns=H2).to_csv(ROOT / "data" / "corrected_h2.tsv", sep="\t", index=False)
    pairs = {}  # a pair of two corrected studies may be estimated from each: the first stands
    for task in tasks:
        for other, values in task["rg"].items():
            if values:
                pairs.setdefault(frozenset((task["study"], other)), [task["study"], other, *values])
    if lost := [f"{study} and {other}" for study, others in zip(corrected, partners) for other in others if frozenset((study, other)) not in pairs]:
        parser.exit(1, f"no log holds the pair of {lost[0]} as their files are now ({len(lost)} such pairs): nothing was written\n")
    pd.DataFrame(pairs.values(), columns=RG).to_csv(ROOT / "data" / "corrected_rg.tsv", sep="\t", index=False)
    print(f"{len(corrected)} corrected studies: {len(pairs)} pairs, {sum(row[2] == 'nan' for row in pairs.values())} of them not estimated")


if __name__ == "__main__":
    main()
