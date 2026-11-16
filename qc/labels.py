"""The label check: does a study's label agree with what it is genetically correlated with?

Two studies carry the same label if their trait names, without the words that do not make a label, are the same
letters, or one begins with the other's and that one has at least two words, or they differ by a slip of the pen.
A study is confirmed if a study with the same label agrees with it. Three findings:

  fail  It is nearly identical to a study with an unrelated label, which studies of that label confirm, and
        nothing confirms its own label: it most likely carries that other trait.
  warn  It is nearly identical to a study with an unrelated label, and nothing confirms its own.
  warn  A study with the same label disagrees with it, and nothing confirms it: one of the two is not what it says.

A finding names the other study as {ID}, and gives its label, the trait name of the release.
"""
import difflib
import re
from itertools import combinations


def words(trait: str, ignore) -> list[str]:
    """The words of a trait name that make its label, in the singular: "LDL Direct Adjusted By Medication" -> ldl, direct."""
    spaced = re.sub(r"(?<=[a-z])(?=[A-Z][a-z])", " ", trait).lower().replace("haem", "hem")  # DirectBilirubin; haemoglobin
    kept = [word for word in re.findall(r"[a-z0-9]+", spaced) if word not in ignore]
    return [word[:-1] if len(word) > 4 and word.endswith("s") and not word.endswith("ss") else word for word in kept]


def same(a: list[str], b: list[str]) -> bool:
    """Whether two labels, as words, are one label."""
    (few, x), (_, y) = sorted(((len(w), "".join(w)) for w in (a, b)), key=lambda label: len(label[1]))
    if x == y:
        return bool(x)
    if few < 2:  # one word is a label of its own: Thyroid is not Thyroid stimulating hormone
        return False
    if y.startswith(x):
        return not y[len(x)].isdigit()  # Hair Balding Pattern is not Hair Balding Pattern 2: a number makes another label
    if len(x) < 15 or re.sub(r"\D", "", x) != re.sub(r"\D", "", y[:len(x) + 1]):  # a slip of the pen, never of a number
        return False
    return max(difflib.SequenceMatcher(None, x, y[:k]).ratio() for k in (len(x), len(x) + 1)) >= 0.95


def related(a: list[str], b: list[str], pairs) -> bool:
    """Whether two labels share a word, or are a pair that may be nearly identical."""
    if {word[:6] for word in a if len(word) >= 4} & {word[:6] for word in b if len(word) >= 4}:
        return True
    x, y = " ".join(a), " ".join(b)
    return any((p in x and q in y) or (q in x and p in y) for p, q in pairs)


def check(studies: dict[str, dict], tables: dict[str, dict], rules: dict) -> dict[str, list[tuple[str, str]]]:
    """{study ID: its findings, each (fail or warn, what was found)}. `studies`: {ID: trait, h2_z}; `tables`: per
    ancestry, `ids` and the lists `rg` and `se` over the pairs of ids, row by row; `rules`: the [label] thresholds."""
    ignore = set(rules["ignore"])
    label = {i: words(d["trait"], ignore) for i, d in studies.items()}
    pairs = [tuple(" ".join(words(part, ignore)) for part in pair) for pair in rules["related"]]  # in the labels' own spelling
    apart = {frozenset(pair) for pair in rules.get("distinct", [])}  # two studies that share a label and are known to measure different things
    name = lambda i: f"{{{i}}} ({studies[i]['trait']})"
    confirmed, conflicts, twins = set(), [], []
    for t in tables.values():
        for (a, b), rg, se in zip(combinations(t["ids"], 2), t["rg"], t["se"], strict=True):
            if rg is None:
                continue
            if same(label[a], label[b]):
                if abs(rg) >= rules["agree_rg"]:
                    confirmed |= {a, b}
                elif abs(rg) + 2 * se < rules["conflict_rg"] and min(studies[a]["h2_z"], studies[b]["h2_z"]) > rules["conflict_h2_z"] and frozenset((a, b)) not in apart:
                    conflicts.append((a, b, rg))
            elif rg >= rules["twin_rg"] and se <= rules["twin_se"] and not related(label[a], label[b], pairs):
                twins.append((a, b, rg))
    found = {}
    for a, b, rg in twins:  # the one whose label nothing confirms carries, it seems, the other's trait
        for mine, other in ((a, b), (b, a)):
            if mine not in confirmed:
                level = "fail" if other in confirmed else "warn"
                found.setdefault(mine, []).append((level, f"label in doubt: rg {rg:.2f} with {name(other)}, a study of another trait by its label; "
                                                   + ("other studies confirm that label, none this study's" if level == "fail" else "one of the two labels is wrong")))
    failed = {i for i, findings in found.items() if any(level == "fail" for level, _ in findings)}
    for a, b, rg in conflicts:  # a failed study explains the disagreement: the other is not blamed for it
        for mine, other in ((a, b), (b, a)):
            if mine not in confirmed and (mine in failed or other not in failed):
                found.setdefault(mine, []).append(("warn", f"label in doubt: rg {rg:.2f} only with {name(other)}, a study of the same trait by its label"))
    return found
