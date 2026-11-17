"""Browser check of the site, served from site/ on this machine.

A headless browser opens every tab: the tables, the filters and the sorting, the ticks and the files they export, the
charts, a study's page with its paper and what the checks found, the candidates, and the login, which must show the
paths only for the right passphrase, keep them over a reload and drop them at logout. What the page shows is compared with site/data, read here independently.

  uv run --with playwright python -m playwright install chromium   # once
  SUMSTATS_PASSPHRASE=... uv run --with playwright python dev/site_check.py

With SITE_CHECK_SHOTS=folder, it saves a picture of each view there, for review.
"""
import functools
import http.server
import json
import os
import re
import sys
import threading
import time
import urllib.parse
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).parents[1]
SITE = ROOT / "site"
sys.path.insert(0, str(ROOT))
from build import crypt  # noqa: E402

PASSPHRASE = os.environ.get("SUMSTATS_PASSPHRASE") or sys.exit("set SUMSTATS_PASSPHRASE: the check logs in")
SHOTS = os.environ.get("SITE_CHECK_SHOTS")
DATA = json.loads((SITE / "data" / "studies.json").read_text())
STUDIES, RELEASE, LABELS, TRAITS = DATA["studies"], DATA["release"], DATA["papers"], DATA["traits"]
BY_ID = {d["id"]: d for d in STUDIES}
PAPERS = json.loads((SITE / "data" / "papers.json").read_text())
CANDIDATES = json.loads((SITE / "data" / "candidates.json").read_text())
RG = {anc: json.loads((SITE / "data" / f"rg_{anc}.json").read_text()) for anc in RELEASE["ancestries"]}
PATHS = json.loads(crypt.decrypt(json.loads((SITE / "data" / "paths.enc.json").read_text()), PASSPHRASE))
PUBMED = "https://pubmed.ncbi.nlm.nih.gov/"
failed = False


class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def end_headers(self):  # never an older copy of the site
        self.send_header("Cache-Control", "no-store")
        super().end_headers()


def check(label, ok, why=None):
    """PASS or FAIL; a failure also prints `why()`, the state it found."""
    global failed
    failed |= not ok
    print("PASS" if ok else "FAIL", label, flush=True)
    if not ok and why:
        try:
            print("    ", why(), flush=True)
        except Exception as e:  # the check goes on
            print("    ", repr(e), flush=True)


def estimate(a, b):
    """The pair's values in site/data: the lists hold the pairs i < j of `ids`, row by row."""
    t = RG[BY_ID[a]["anc"]]
    i, j = sorted((t["ids"].index(a), t["ids"].index(b)))
    k = i * (2 * len(t["ids"]) - i - 1) // 2 + j - i - 1
    return {c: t[c][k] for c in ("rg", "se", "z", "p", "gcov", "gcov_se")} | {"note": t["notes"].get(str(k))}


def neighbours(a):
    """Every other study of a's ancestry with its estimate against a, the strongest first."""
    pairs = [{"id": b, **estimate(a, b)} for b in RG[BY_ID[a]["anc"]]["ids"] if b != a]
    return sorted(pairs, key=lambda r: -1 if r["rg"] is None else abs(r["rg"]), reverse=True)


def file_path(study):
    folder, name = PATHS["files"][study]
    return f"{PATHS['folders'][folder]}/{name}"


def sig(x):
    """A number as the site writes it: two significant digits, in full, with a true minus sign. From the number's exact
    value, a half rounded away from 0, as the browser does."""
    if x is None:
        return "–"
    if x == 0:
        return "0"
    exact = Decimal(x)
    rounded = exact.quantize(Decimal(1).scaleb(exact.adjusted() - 1), rounding=ROUND_HALF_UP)
    return f"{rounded:.{max(0, 1 - rounded.adjusted())}f}".replace("-", "−")


def plus_minus(value, se):
    return "–" if value is None else f"{sig(value)} ± {sig(se)}"


def label(d):
    """A study's paper as the table names it."""
    return LABELS[d["cite"][0]]["label"] if "cite" in d else d["ref"]


def published(d):
    """The day a study's paper was published: 2010-10-10, or 2018-07 where PubMed gives no day."""
    return LABELS[d["cite"][0]]["published"] if "cite" in d else ""


def said(finding):
    """A finding of the checks as the site words it: another study, {ID} in the data, by its short name."""
    return re.sub(r"\{([^}]+)\}", lambda found: BY_ID[found[1]]["name"], finding)


def text(d):
    """What the search reads of a study: its short name, its ID and name in the release, its trait with its term and
    category, its release, its paper."""
    return "\n".join([d["name"], d["id"], d["trait"], d["of"], *TRAITS[d["of"]], d["origin"], d["ref"], *(LABELS[pmid]["label"] for pmid in d.get("cite", []))]).lower()


def table_of(studies, paths):
    """The file of the studies ticked: a line each, with its trait, its cases, controls and effective N, what the checks say, whether it is corrected (by a copy, or its label
    alone), whether it is the study to use of its trait or which supersedes it, and the path after login."""
    head = ["short_name", "id", "release_name", "trait", "ancestry", "n", "cases", "controls", "effective_n", "h2", "h2_se", "h2_liability", "h2_liability_se", "reference", "published", "pmid", "qc",
            "corrected", "recommended", "superseded_by", "uk_biobank", *(["path"] if paths else [])]
    corrected = lambda d: "copy" if any(kind != "relabel" for kind, _, _ in d.get("fixed", [])) else "label" if "fixed" in d else ""
    lines = [[d["name"], d["id"], d["trait"], d["of"], d["anc"], str(d["n"]), str(d.get("cases", "")), str(d.get("controls", "")), str(d.get("neff", "") if d.get("kind") == "binary" else d["n"]),
              *(str(d.get(key, "")) for key in ("h2", "h2_se", "h2l", "h2l_se")), label(d), published(d), ";".join(d.get("cite", [])), d["qc"], corrected(d), "yes" if "rec" in d else "",
              BY_ID[d["sup"]]["name"] if "sup" in d else "", {True: "yes", False: "no"}.get(d.get("ukb"), ""), *([file_path(d["id"])] if paths else [])]
             for d in studies]
    return "".join("\t".join(line) + "\n" for line in [head, *lines])


def references_of(studies):
    """The BibTeX file of the studies ticked: an entry a paper, the earliest published first, under a comment that names its studies."""
    cited = {}
    for d in studies:
        for pmid in d.get("cite", []):
            cited.setdefault(pmid, []).append(f"{d['name']} ({d['id']})")
    return "\n".join(f"% {', '.join(cited[pmid])}\n{PAPERS[pmid]['bib']}" for pmid in sorted(cited, key=lambda pmid: (PAPERS[pmid]["published"], PAPERS[pmid]["key"])))


def run(pw, url, name, **context):
    """Every check, in one browser window."""
    browser = pw.chromium.launch()
    ctx = browser.new_context(**context)
    ctx.grant_permissions(["clipboard-read", "clipboard-write"])
    page = ctx.new_page()
    errors, elsewhere = [], []
    page.on("console", lambda m: m.type in ("error", "warning") and errors.append(m.text))
    page.on("pageerror", lambda e: errors.append(str(e)))
    page.on("requestfailed", lambda r: errors.append(f"request failed: {r.url}"))
    page.on("request", lambda r: r.url.startswith(url) or r.url.startswith(("data:", "blob:")) or elsewhere.append(r.url))
    count = lambda: page.inner_text("#count")
    bodies = lambda table: page.locator(f"#{table} tbody")
    first = lambda table: bodies(table).first.locator("td.study").inner_text()
    fits = lambda: page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    private = lambda: [f for f in [PATHS["release"], *PATHS["folders"]] if f in page.content()]
    heads = lambda table: [" ".join(text.split()) for text in page.locator(f"#{table} th").all_inner_texts()]
    cells = lambda table: [" ".join(text.split()) for text in bodies(table).first.locator("td").all_inner_texts()]  # of its first row
    kept = lambda: json.loads(page.evaluate("localStorage.getItem('cart')") or "[]")

    def shot(view, full=True):
        if SHOTS:
            page.screenshot(path=f"{SHOTS}/{name}-{view}.png", full_page=full)

    def show(query, ready):
        page.goto(url + query)
        expect(page.locator(ready).first).to_be_visible(timeout=20000)

    def downloaded(link):
        """The name and the text of the file that a click on the cart's link hands over."""
        with page.expect_download() as download:
            page.locator("#cart a", has_text=link).click()
        return download.value.suggested_filename, Path(download.value.path()).read_text(encoding="utf-8")

    # the studies
    show("", "#sumstats-table tbody")
    check(f"[{name}] the table lists every study", bodies("sumstats-table").count() == len(STUDIES) == 602 and count() == "602 studies", count)
    d = STUDIES[0]
    check(f"[{name}] its columns: a study by its short name; h² ± SE and h² z to two significant digits; what the checks say; no ID, no trait, no Signed, no Release",
          heads("sumstats-table") == ["", "Study", "Ancestry", "N", "M", "h² ± SE", "h² z", "QC", "UKB", "Reference"]
          and cells("sumstats-table") == ["", d["name"], d["anc"], f"{d['n']:,}", f"{d['m']:,}", plus_minus(d["h2"], d["h2_se"]), sig(d["h2_z"]), d["qc"], {True: "yes", False: "no"}.get(d.get("ukb"), "–"), label(d)],
          lambda: (heads("sumstats-table"), cells("sumstats-table")))
    links = page.locator("#sumstats-table td.ref a").evaluate_all("list => list.map((a) => [a.textContent, a.href, a.target, a.rel])")
    check(f"[{name}] every study's reference names its paper and leads to PubMed, in a new tab",
          links == [[label(d), f"{PUBMED}{d['cite'][0]}/", "_blank", "noopener noreferrer"] for d in STUDIES], lambda: links[:3])
    shot("sumstats", full=False)
    ok = True
    for sought in ("height uk", d["id"], d["name"].lower(), "body mass index", "bycroft"):  # words; an ID; a short name; a name in the release; a paper
        page.fill("#q", sought)
        found = [s for s in STUDIES if all(w in text(s) for w in sought.lower().split())]
        ok &= 0 < len(found) < 602 and bodies("sumstats-table").count() == len(found) and count() == f"{len(found)} of 602 studies"
    check(f"[{name}] the search finds every word, in short name, ID, name in the release, trait, release or reference", ok, count)
    page.fill("#q", "")
    page.select_option("#f-anc", "EAS")
    eas = [d for d in STUDIES if d["anc"] == "EAS"]
    ok = count() == "56 of 602 studies"
    for box, keep in (("z6", lambda d: d["h2_z"] > 6), ("signed", lambda d: d["signed"]), ("panel", lambda d: "panel" in d)):
        page.check(f"#f-{box}")
        ok &= count() == f"{sum(map(keep, eas))} of 602 studies"
        page.uncheck(f"#f-{box}")
    page.check("#f-panel")
    ok &= count() == "17 of 602 studies"
    page.select_option("#f-anc", "EUR")
    ok &= count() == f"{sum('panel' in d for d in STUDIES if d['anc'] == 'EUR')} of 602 studies"
    page.uncheck("#f-panel")
    page.select_option("#f-anc", "")
    page.check("#f-rec")
    ok &= count() == f"{RELEASE['recommended']['studies']} of 602 studies" and set(page.locator("#sumstats-table td .qc").all_inner_texts()) == {"pass"}
    page.uncheck("#f-rec")
    for value, keep in (("yes", lambda d: d.get("ukb") is True), ("no", lambda d: d.get("ukb") is False), ("unknown", lambda d: "ukb" not in d)):
        page.select_option("#f-ukb", value)
        ok &= count() == f"{sum(map(keep, STUDIES))} of 602 studies" and RELEASE["ukb"][value] == sum(map(keep, STUDIES))
    page.select_option("#f-ukb", "")
    check(f"[{name}] the filters count as the data does: ancestry, h² z > 6, signed, nonredundant, recommended, UK Biobank", ok and count() == "602 studies", count)
    ok = True
    for level in ("pass", "warn", "fail"):
        page.select_option("#f-qc", level)
        ok &= count() == f"{RELEASE['qc'][level]} of 602 studies" and set(page.locator("#sumstats-table td .qc").all_inner_texts()) == {level}
    page.select_option("#f-qc", "")
    page.locator("#sumstats-table th", has_text="QC").locator("button").click()  # the worst first
    failing = [d for d in STUDIES if d["qc"] == "fail"]
    ok &= page.locator("#sumstats-table td.study").all_inner_texts()[:len(failing) + 1] == [d["name"] for d in failing] + [next(d["name"] for d in STUDIES if d["qc"] == "warn")]
    why = bodies("sumstats-table").first.locator(".qc").get_attribute("title")
    ok &= why == "\n".join(f"{level}: {said(text)}" for level, text in failing[0]["notes"] if level != "note")
    for _ in range(2):
        page.locator("#sumstats-table th", has_text="QC").locator("button").click()
    check(f"[{name}] the checks' verdict filters and sorts the studies, the worst first, and says why under the pointer", ok and first("sumstats-table") == STUDIES[0]["name"],
          lambda: (count(), why))
    by_n = sorted(STUDIES, key=lambda d: d["n"])
    sort_n = page.locator("#sumstats-table th:nth-child(4) button")  # N
    order = []
    for _ in range(3):
        sort_n.click()
        order.append(first("sumstats-table"))
    check(f"[{name}] a column sorts largest first, then smallest first, then as released",
          order == [by_n[-1]["name"], by_n[0]["name"], STUDIES[0]["name"]], lambda: order)
    by_paper = [d["name"] for d in sorted(STUDIES, key=lambda d: (published(d), label(d), d["name"]))]  # the earliest paper first; then first author; then short name
    sort_reference, order = page.locator("#sumstats-table th", has_text="Reference").locator("button"), []
    for _ in range(3):
        sort_reference.click()
        order.append(page.locator("#sumstats-table td.study").all_inner_texts())
    days = [published(BY_ID[next(d["id"] for d in STUDIES if d["name"] == study)]) for study in order[0]]
    check(f"[{name}] the Reference column sorts by the day the paper was published, the earliest first, then the latest first, then as released",
          order == [by_paper, by_paper[::-1], [d["name"] for d in STUDIES]] and days == sorted(days) and days[0] < "2011" and days[-1] > "2026", lambda: (order[0][:4], days[:4], days[-3:]))
    check(f"[{name}] logged out, no path is in the page", not private() and page.locator(".paths").count() == 0, private)

    # the cart: ticks, and the files they export
    boxes = page.locator("#sumstats-table td.tick input")
    boxes.nth(0).check()
    boxes.nth(2).check()
    mine = [STUDIES[0], STUDIES[2]]
    all_box = page.locator("#sumstats-table th.tick input")
    check(f"[{name}] ticking studies counts them and offers the exports",
          page.text_content("#cart") == " · 2 ticked: download IDs (TSV) or references (BibTeX) · untick all" and kept() == [d["id"] for d in mine]
          and all_box.evaluate("box => box.indeterminate && !box.checked"), lambda: page.text_content("#cart"))
    shot("cart", full=False)
    got = downloaded("IDs (TSV)")
    check(f"[{name}] the studies ticked download as a table: short name, ID, name in the release, ancestry, N, reference, PMID; no path",
          got == ("gwas-sumstats-studies.tsv", table_of(mine, paths=False)), lambda: got)
    got = downloaded("references (BibTeX)")
    check(f"[{name}] their papers download as BibTeX, an entry a paper, under the studies that cite it",
          got == ("gwas-sumstats-references.bib", references_of(mine)) and got[1].count("@article{") == len({d["cite"][0] for d in mine}), lambda: got[1][:300])
    page.fill("#q", "teslovich")
    all_box.check()
    mine = [d for d in STUDIES if d in mine or "teslovich" in text(d)]
    ok = len(mine) == 5 and bodies("sumstats-table").count() == 3 and all_box.is_checked() and kept() == [d["id"] for d in mine]  # in the order ticked: the two, then the three
    page.fill("#q", "")
    page.check("#f-ticked")
    ok &= count() == "5 of 602 studies" and page.locator("#sumstats-table td.study").all_inner_texts() == [d["name"] for d in mine]
    page.reload()
    expect(page.locator("#cart")).to_contain_text("5 ticked", timeout=20000)
    ok &= boxes.evaluate_all("list => list.filter((box) => box.checked).length") == 5
    page.click("nav a[data-tab=h2]")
    expect(page.locator("#h2-table tbody").first).to_be_visible()
    ok &= page.locator("#h2-table td.tick input").evaluate_all("list => list.filter((box) => box.checked).length") == 5
    page.click("nav a[data-tab=sumstats]")
    check(f"[{name}] the header's box ticks the studies listed; the ticks stay over a reload and in the other table, and filter the studies", ok,
          lambda: (kept(), count()))

    # login
    page.click("#login-open")
    page.fill("#login-pass", PASSPHRASE[::-1] + "?")
    page.click("#login button[type=submit]")
    expect(page.locator("#login-note")).to_have_text("Wrong passphrase.", timeout=20000)
    shot("login-wrong", full=False)
    check(f"[{name}] a wrong passphrase is refused and nothing is kept",
          not private() and page.evaluate("localStorage.getItem('key')") is None and page.inner_text("#login-open") == "login")
    page.fill("#login-pass", PASSPHRASE)
    t0 = time.time()
    page.click("#login button[type=submit]")
    expect(page.locator("#sumstats-table .paths").first).to_be_visible(timeout=30000)
    print(f"NOTE [{name}] login took {time.time() - t0:.1f} s", flush=True)
    shown = page.locator("#sumstats-table tbody").first.locator(".path").all_inner_texts()
    check(f"[{name}] the passphrase shows each study's path",
          page.locator("#sumstats-table .paths").count() == 602 and shown == [file_path(STUDIES[0]["id"])] and page.inner_text("#login-open") == "logout"
          and page.locator("#login").is_hidden() and page.input_value("#login-pass") == "", lambda: shown)
    shot("sumstats-logged-in", full=False)
    page.locator("#sumstats-table tbody").first.locator(".paths a").click()
    check(f"[{name}] copy puts the path on the clipboard", page.evaluate("navigator.clipboard.readText()") == file_path(STUDIES[0]["id"]))
    got = downloaded("IDs and paths (TSV)")
    check(f"[{name}] after login, the table of the studies ticked has their paths", got == ("gwas-sumstats-studies.tsv", table_of(mine, paths=True)), lambda: got[1][:300])
    page.reload()
    expect(page.locator("#sumstats-table .paths").first).to_be_visible(timeout=30000)
    check(f"[{name}] a reload stays logged in", page.inner_text("#login-open") == "logout")
    added = next(d["id"] for d in STUDIES if d["release"] == RELEASE["year"] and d["cross"])
    show(f"?study={added}", "#view-study .paths")
    shown = page.locator("#view-study .path").all_inner_texts()
    line = page.locator("#view-study .row", has=page.locator(".paths"))
    check(f"[{name}] an added study's page has one line, Path, with its file and copy",
          shown == [file_path(added)] and line.locator(".label").inner_text() == "Path" and line.locator("a").all_inner_texts() == ["copy"]
          and "h² file" not in page.inner_text("#view-study") and "rg file" not in page.inner_text("#view-study"), lambda: shown)
    shot("study-logged-in")
    show("?tab=h2", "#h2-table .paths")
    check(f"[{name}] the heritability table shows the paths too", page.locator("#h2-table .paths").count() == 602 and fits())
    show("?tab=notes", "#o2-release")
    check(f"[{name}] the release notes say where the release is", PATHS["release"] in page.inner_text("#o2-release") and fits())
    shot("notes-logged-in")
    show("?tab=about", "#o2-folders")
    folders = page.locator("#o2-folders li").all_inner_texts()
    check(f"[{name}] About lists the folders of the files, each with its number of files",
          folders == [f"{folder}: {count:,} files" for folder, count in ((f, sum(file[0] == k for file in PATHS["files"].values())) for k, f in enumerate(PATHS["folders"]))]
          and fits(), lambda: folders)
    shot("about-logged-in")
    page.click("#login-open")
    expect(page.locator("#o2-folders")).to_be_hidden()
    show("", "#sumstats-table tbody")
    check(f"[{name}] logout hides the paths and forgets the key; the export has no path again",
          not private() and page.locator(".paths").count() == 0 and page.evaluate("localStorage.getItem('key')") is None and page.inner_text("#login-open") == "login"
          and downloaded("IDs (TSV)")[1] == table_of(mine, paths=False), private)
    page.locator("#cart a", has_text="untick all").click()
    check(f"[{name}] untick all empties the cart", page.locator("#cart").is_hidden() and kept() == [] and boxes.evaluate_all("list => list.every((box) => !box.checked)"))

    # heritability
    show("?tab=h2", "#h2-plot svg")
    top = max(STUDIES, key=lambda d: d["h2"])
    bins = [round(d["h2"] * 1e4) // 200 for d in STUDIES]  # of 0.02
    mode = max(sorted(set(bins)), key=bins.count)
    check(f"[{name}] heritability: a column per bin of h², and the studies by h²",
          page.inner_text("#h2-caption") == f"h² of the 602 studies shown, in bins of 0.02: most often {mode * 0.02:.2f} to {(mode + 1) * 0.02:.2f}."
          and first("h2-table") == top["name"] and page.locator("#h2-plot .bar").count() == len(set(bins)) and fits(),
          lambda: page.inner_text("#h2-caption"))
    check(f"[{name}] its table: h² ± SE in one column, sorted by h², without bars",
          heads("h2-table") == ["", "Study", "Ancestry", "h² ± SE ▼", "Liability h² ± SE", "h² z", "Intercept", "QC"] and page.locator(".meter").count() == 0
          and cells("h2-table") == ["", top["name"], top["anc"], plus_minus(top["h2"], top["h2_se"]), plus_minus(top.get("h2l"), top.get("h2l_se")), sig(top["h2_z"]), f"{top['intercept']:.4f}", top["qc"]],
          lambda: (heads("h2-table"), cells("h2-table")))
    shot("h2", full=False)
    box = page.locator("#h2-plot svg").bounding_box()
    page.mouse.move(box["x"] + box["width"] * 0.12, box["y"] + box["height"] * 0.5)
    check(f"[{name}] a column says how many studies it stands for", bool(re.fullmatch(r"\d+ (study|studies) with h² from [\d.]+ to [\d.]+", page.inner_text("#tip"))), lambda: page.inner_text("#tip"))
    shot("h2-hover", full=False)
    page.mouse.move(0, 0)

    # genetic correlation
    eur, eas_ids = RG["EUR"]["ids"], RG["EAS"]["ids"]
    noted = int(next(iter(RG["EUR"]["notes"])))  # a pair whose missing estimate has a reason
    rows = [i for i in range(len(eur)) if i * (2 * len(eur) - i - 1) // 2 <= noted]
    i = rows[-1]
    j = noted - i * (2 * len(eur) - i - 1) // 2 + i + 1
    pairs = [(eur[0], eur[-1]), (eur[-1], eur[1]), (eas_ids[3], eas_ids[-2]), (eur[i], eur[j])]
    ok, seen = True, []
    for a, b in pairs:
        e = estimate(a, b)
        show(f"?tab=rg&a={a}&b={b}", "#pair-result .fields")
        seen.append(page.locator("#pair-result .row").first.locator(".value").inner_text())
        ok &= seen[-1] == (plus_minus(e["rg"], e["se"]) if e["rg"] is not None else f"not estimated: LDSC gave none for this pair ({e['note']})")
        ok &= page.input_value("#rg-a") == BY_ID[a]["name"] and page.input_value("#rg-b") == BY_ID[b]["name"]
        ok &= page.inner_text("#pair-result h3") == f"{BY_ID[a]['name']} and {BY_ID[b]['name']}"
    check(f"[{name}] a pair's genetic correlation is the one in the data, in either order, in EUR and EAS, or says why there is none", ok, lambda: seen)
    a = "UKB_460K.body_HEIGHTz"
    near = neighbours(a)
    strong = [r for r in near if r["rg"] is not None and abs(r["z"]) >= 3]
    show(f"?tab=rg&a={a}&b={strong[1]['id']}", "#neighbours-table tbody")
    listed = lambda: page.locator("#neighbours-table td.study").all_inner_texts()
    ok = listed() == [BY_ID[r["id"]]["name"] for r in strong[:10]] and fits()
    ok &= heads("neighbours-table") == ["Study", "rg ± SE", "z", "p"] and page.inner_text("#neighbours-head").startswith(f"{BY_ID[a]['name']}: {len(strong)} of ")
    ok &= cells("neighbours-table")[1:3] == [plus_minus(strong[0]["rg"], strong[0]["se"]), sig(strong[0]["z"])]
    shot("rg")
    page.click("#neighbours-more")
    ok &= bodies("neighbours-table").count() == len(strong)
    page.uncheck("#rg-strict")
    page.click("#neighbours-more")
    ok &= bodies("neighbours-table").count() == len(near) == len(eur) - 1
    page.locator("#neighbours-table th:nth-child(3) button").click()  # z
    by_z = max((r for r in near if r["z"] is not None), key=lambda r: r["z"])
    ok &= listed()[0] == BY_ID[by_z["id"]]["name"]
    check(f"[{name}] a study's genetic correlations: the strongest ten with |z| ≥ 3, then all, then every pair, sortable", ok, lambda: listed()[:10])
    show("?tab=rg", "#heat svg")
    ok = True
    for sought in ("schizo", "UKB_460K.body_HEIGHTz", "body mass index"):  # a short name; an ID; a name in the release
        page.fill("#rg-a", sought)
        options = page.locator("#rg-a-list a")
        expect(options.first).to_be_visible()
        matches = [d for d in STUDIES if d["cross"] and all(w in text(d) for w in sought.lower().split())]
        ok &= options.all_inner_texts() == [f"{d['name']} {d['anc']}" for d in matches[:60]] and len(matches) > 0
    page.fill("#rg-a", "schizo")
    matches = [d for d in STUDIES if d["cross"] and "schizo" in text(d)]
    shot("rg-picker", full=False)
    page.keyboard.press("ArrowDown")
    page.keyboard.press("ArrowDown")
    page.keyboard.press("Enter")
    expect(page.locator("#neighbours-table tbody").first).to_be_visible(timeout=20000)
    check(f"[{name}] typing a short name, an ID or a name in the release lists the signed studies that match, by short name; arrows and Enter pick one",
          ok and page.evaluate("new URLSearchParams(location.search).get('a')") == matches[1]["id"] and page.input_value("#rg-a") == matches[1]["name"],
          lambda: (page.locator("#rg-a-list a").all_inner_texts()[:5], page.url))
    order = lambda anc, kind: [BY_ID[RG[anc]["ids"][i]] for i in RG[anc]["heat"][kind]["order"]]  # the build's clustering

    def named(rows):
        """The heatmap names these studies beside its rows, in this order (a long name cut short)."""
        shown = page.locator("#heat text").all_text_contents()
        return len(shown) == len(rows) and all(d["name"].startswith(s.rstrip("…")) for s, d in zip(shown, rows))

    def tip_for(a, b):
        e = estimate(a["id"], b["id"])
        return ("not estimated" if e["rg"] is None else plus_minus(e["rg"], e["se"])) + f"\n{a['name']} × {b['name']}"

    def cell(rows, i, j):
        """The middle of the heatmap's cell in row i, column j, scrolled to the middle of the window."""
        def middle():
            box = page.locator("#heat image").bounding_box()
            return box["x"] + (j + 0.5) * box["width"] / len(rows), box["y"] + (i + 0.5) * box["height"] / len(rows)
        page.evaluate("y => scrollBy(0, y - innerHeight / 2)", middle()[1])
        return middle()

    show("?tab=rg", "#heat image")
    ok = named(order("EUR", "panel")) and page.locator("#heat .tree").count() == 1 and fits()
    ends = page.locator("#heat-caption .steps i").evaluate_all("list => [list[0], list.at(-1)].map((i) => getComputedStyle(i).backgroundColor.match(/\\d+/g).map(Number))")
    ok &= ends[0][2] > ends[0][0] and ends[1][0] > ends[1][2] and "Blue below 0, red above" in page.inner_text("#heat-caption")  # blue for -, red for +
    shot("rg-heat")
    page.locator("#heat-tabs a", has_text="EAS nonredundant").click()
    expect(page.locator("#heat text")).to_have_count(17, timeout=20000)
    rows = order("EAS", "panel")
    ok &= named(rows)
    page.mouse.move(*cell(rows, 4, 14))
    tip = page.inner_text("#tip")
    ok &= tip == tip_for(rows[4], rows[14])
    shot("rg-heat-eas", full=False)
    page.mouse.click(*cell(rows, 4, 14))
    ok &= page.evaluate("location.search") == f"?tab=rg&a={rows[4]['id']}&b={rows[14]['id']}" and named(rows)
    check(f"[{name}] a panel's heatmap is in the order of the build's clustering; a cell names its pair, by short names, and looks it up", ok, lambda: (tip, page.url))
    show("?tab=rg", "#heat image")
    ok = True
    for anc in RG:
        page.locator("#heat-tabs a", has_text=f"{anc} all").click()
        rows = order(anc, "all")
        expect(page.locator("#heat svg")).to_have_attribute("aria-label", re.compile(f"the {len(rows)} signed {anc} studies"), timeout=20000)
        ok &= len(rows) == len(RG[anc]["ids"]) and fits() and (named(rows) if len(rows) <= 100 else page.locator("#heat text").count() == 0)
        for i, j in ((0, len(rows) - 1), (len(rows) - 2, 1), (len(rows) // 3, len(rows) // 2)):
            page.mouse.move(*cell(rows, i, j))
            tip = page.inner_text("#tip")
            ok &= tip == tip_for(rows[i], rows[j])
        shot(f"rg-heat-all-{anc}")
        page.mouse.click(*cell(rows, i, j))  # the pair is looked up, and the heatmap stays
        expect(page.locator("#pair-result .fields")).to_be_visible(timeout=20000)
        ok &= page.evaluate("location.search") == f"?tab=rg&a={rows[i]['id']}&b={rows[j]['id']}" and page.inner_text("#heat-tabs .current").startswith(f"{anc} all")
    page.mouse.move(0, 0)
    check(f"[{name}] the heatmap of all signed studies of an ancestry, clustered too: a cell names its pair", ok, lambda: tip)

    # a study's page, the texts, the address
    show(f"?study={added}", "#view-study .pairs tbody")
    d, strong = BY_ID[added], [r for r in neighbours(added) if r["rg"] is not None and abs(r["z"]) >= 3]
    fields = dict(zip(page.locator("#view-study .label").all_inner_texts(), page.locator("#view-study .value").all_inner_texts()))
    links = page.locator("#view-study .fields a").evaluate_all("list => list.map((a) => a.href)")
    check(f"[{name}] a study's page: its short name as title; its ID and name in the release; its fields; its strongest genetic correlations",
          page.inner_text("#view-study h2") == d["name"] and fields["ID"] == added and fields["Release name"] == d["trait"] and fields["N"] == f"{d['n']:,}"
          and fields["Release"] == d["origin"] and fields["Signed"] == "yes" and fields["h² ± SE"] == plus_minus(d["h2"], d["h2_se"]) and fields["h² z"] == sig(d["h2_z"])
          and fields["Path"] == "login to see it" and d["catalog"] in links
          and page.locator("#view-study .pairs td.study").all_inner_texts() == [BY_ID[r["id"]]["name"] for r in strong[:10]]
          and page.title() == f"{d['name']} · GWAS Sumstats" and fits(), lambda: fields)
    shot("study")
    both = next(d for d in STUDIES if len(d.get("cite", [])) == 2)  # the paper of its statistics, and the resource it names
    show(f"?study={both['id']}", "#view-study .citation a")
    expect(page.locator("#view-study .citation", has_text=PAPERS[both["cite"][0]]["title"])).to_be_visible(timeout=20000)
    cited = page.locator("#view-study .citation").all_inner_texts()
    links = page.locator("#view-study .citation a").evaluate_all("list => list.map((a) => [a.textContent, a.href])")

    def in_full(pmid):
        p = PAPERS[pmid]
        authors = ", ".join(author.strip("{}") for author in p["authors"]) + (", et al" if p["more"] else "")
        issue = p["volume"] + (f"({p['issue']})" if p["issue"] else "")
        return (f"{authors}. {p['title']}. {p['journal']} {p['year']}" + (f";{issue}" if issue else "") + (f":{p['pages']}" if p["pages"] else "")
                + f". Published {p['published'].replace('-', '/')}. PubMed" + (" · DOI" if p["doi"] else ""))

    check(f"[{name}] a study's page cites its papers in full, with their PubMed and DOI links",
          cited == [in_full(pmid) for pmid in both["cite"]]
          and links == [pair for pmid in both["cite"] for pair in (["PubMed", f"{PUBMED}{pmid}/"], ["DOI", f"https://doi.org/{PAPERS[pmid]['doi']}"])], lambda: (cited, links))
    shot("study-papers", full=False)
    page.locator("#view-study > .check input").check()
    ok = kept() == [both["id"]]
    show("", "#sumstats-table tbody")
    ok &= page.text_content("#cart").startswith(" · 1 ticked") and page.locator("#sumstats-table td.tick input:checked").count() == 1
    page.locator("#cart a", has_text="untick all").click()
    check(f"[{name}] a study is ticked on its page too", ok and kept() == [], kept)
    doubted = BY_ID["PANUKBB_LDLDirectAdjustedByMedication"]  # a label that its genetic correlations contradict; corrected since, so that its page has both verdicts
    show(f"?study={doubted['id']}", "#view-study .finding")

    def told(label, notes):
        """The block of a verdict on the page, and its findings if they are these, each naming other studies by their short names, linked."""
        block = page.locator("#view-study .row", has=page.locator(".label", has_text=re.compile(f"^{label}$"))).locator(".value")
        findings = [" ".join(line.split()) for line in block.locator(".finding").all_inner_texts()]
        named_there = block.locator(".finding a").evaluate_all("list => list.map((a) => [a.textContent, a.getAttribute('href')])")
        right = findings == [f"{level} {said(text)}" for level, text in notes] and named_there == [[BY_ID[study]["name"], f"?study={study}"] for _, text in notes for study in re.findall(r"\{([^}]+)\}", text)]
        return block, findings if right else []

    block, findings = told("QC as corrected", doubted["notes"])
    before, released = told("QC as released", doubted["was"]["notes"])
    report = block.locator("a", has_text="report a problem").get_attribute("href")
    check(f"[{name}] a study's page says what the checks found, the worst first, naming other studies by their short names, with what was measured of its file",
          block.locator("div div").first.inner_text() == "warn · report a problem" and before.locator("div div").first.inner_text() == "fail" and before.locator(".hint").count() == 0
          and block.locator(".hint").inner_text() == f"Its file: {int(doubted['hm3'] * 100)}% of the HapMap3 SNPs that LDSC regresses on; mean χ² {doubted['chi2']:.2f}; λGC {doubted['lambda']:.2f}."
          and findings and findings[0].startswith("warn the release calls it ") and released and released[0].startswith("fail ") and any("FEV1FVC" in finding for finding in released), lambda: (findings, released))
    check(f"[{name}] its link to report a problem opens the repository's form with the study filled in, and no path",
          report.startswith("https://github.com/gwassumstats/gwassumstats.github.io/issues/new?template=study.yml&")
          and dict(urllib.parse.parse_qsl(report.split("?", 1)[1])) == {"template": "study.yml", "title": f"{doubted['name']}: ", "study": f"{doubted['name']} ({doubted['id']})"},
          lambda: report)
    shot("study-qc", full=False)
    # the corrections: studies shown, and checked, as corrected, with what the release has beside
    corrected = [d for d in STUDIES if "fixed" in d]
    show("", "#sumstats-table tbody")
    page.select_option("#f-qc", "corrected")
    tags = page.locator("#sumstats-table tbody tr:first-child td:nth-child(8)").all_inner_texts()
    check(f"[{name}] the corrected studies are listed on demand, each with what the checks say of it as corrected",
          tags == [f"{d['qc']} corrected" for d in corrected] and len(corrected) == RELEASE["corrected"]["studies"] == 37 and page.inner_text("#count") == "37 of 602 studies"
          and len({d["qc"] for d in corrected}) == 3, lambda: tags[:3])
    shot("corrected", full=False)
    page.select_option("#f-qc", "")
    row = lambda label: " ".join(page.locator("#view-study .row", has=page.locator(".label", has_text=re.compile(f"^{re.escape(label)}$"))).locator(".value").inner_text().split())
    d = BY_ID["PASS_AtrialFibrillation_Nielsen2018"]
    show(f"?study={d['id']}", "#view-study .fix")
    told, w = row("Corrected"), d["was"]
    check(f"[{name}] a corrected study's page says what was done and why, and what the release has; the checks of the copy and of the release's file are told apart",
          told.startswith("taken again from its source, GCST006414, and formatted as the release formats: the alleles of the release's file do not go with its Z")
          and told.endswith(f"The numbers on this page are of the corrected copy. In the release: h² {plus_minus(w['h2'], w['h2_se'])}, z {sig(w['h2_z'])}; {w['m']:,} SNPs.")
          and page.locator("#view-study .fix a", has_text="GCST006414").get_attribute("href") == "https://www.ebi.ac.uk/gwas/studies/GCST006414"
          and row("Nonredundant panel") == f"yes: {d['panel']} of {sum('panel' in s for s in STUDIES if s['anc'] == 'EUR')} in the EUR panel, by h² z; in the release's panel: yes"
          and row("h² ± SE") == plus_minus(d["h2"], d["h2_se"]) and row("QC as corrected").startswith("pass · report a problem note nearly identical to 1 other study")
          and "Its corrected copy: " in row("QC as corrected") and row("QC as released").startswith("fail fail its Z does not go with its alleles")
          and page.locator("#view-study .label", has_text=re.compile("^QC$")).count() == 0, lambda: (told, row("QC as corrected"), row("QC as released")))
    shot("study-corrected", full=False)
    other, t = "PANUKBB_AtrialFibrillationAndFlutter", RG["EUR"]
    i, j = sorted((t["ids"].index(d["id"]), t["ids"].index(other)))
    k = i * (2 * len(t["ids"]) - i - 1) // 2 + j - i - 1
    show(f"?tab=rg&a={d['id']}&b={other}", "#pair-result .fields")
    before = " ".join(page.locator("#pair-result .row", has_text="In the release").inner_text().split())
    check(f"[{name}] a pair estimated again on a corrected copy shows the release's estimate beside the new one",
          before == f"In the release rg {plus_minus(*t['was'][str(k)])}: estimated again on a corrected copy" and t["rg"][k] > 0.9 > 0 > t["was"][str(k)][0], lambda: before)
    relabelled = BY_ID["PANUKBB_IndirectBilirubin"]
    show(f"?study={relabelled['id']}", "#view-study .fix")
    check(f"[{name}] a study whose label alone is corrected says what its file holds: it warns under that trait, where it failed under the release's",
          row("Corrected").startswith("its file holds another trait: Estimated glomerular filtration rate (creatinine, cystatin C): its N, 401,570,") and "h2" not in relabelled["was"]
          and "The numbers on this page" not in row("Corrected") and row("QC as corrected").startswith("warn · report a problem warn the release calls it Indirect Bilirubin; its file holds Estimated glomerular filtration rate")
          and "Its file: " in row("QC as corrected") and row("QC as released").startswith("fail fail its N, 401,570, is that of the Pan-UK Biobank's eGFR"), lambda: (row("QC as corrected"), row("QC as released")))
    untouched = next(d for d in STUDIES if "fixed" not in d and d["qc"] == "warn")
    show(f"?study={untouched['id']}", "#view-study .fields")
    check(f"[{name}] a study without a correction has one verdict", row("QC").startswith("warn") and page.locator("#view-study .label", has_text=re.compile("^QC as ")).count() == 0, lambda: row("QC"))
    # cases, controls and effective N; the study to use of a trait
    def settled(study):
        """A study's page, once its genetic correlations have loaded: leaving it earlier would break off their request."""
        show(f"?study={study}", "#view-study .fields")
        expect(page.locator("#view-study .section-head", has_text="Loading the genetic correlations")).to_have_count(0, timeout=20000)

    veins = BY_ID["PASS_VaricoseVeins_Ahmed2022"]
    settled(veins["id"])
    counted = row("Cases and controls")
    check(f"[{name}] a binary study's page gives its cases and controls, its effective N and where they are from",
          counted == f"{veins['cases']:,} and {veins['controls']:,}: effective N {veins['neff']:,} · the paper The paper's 22,473 cases and 379,183 controls of UK Biobank; they sum to its N."
          and page.locator("#view-study .row", has_text="Cases and controls").locator("a", has_text="the paper").get_attribute("href") == "https://pubmed.ncbi.nlm.nih.gov/35654884/"
          and veins["neff"] == round(4 / (1 / veins["cases"] + 1 / veins["controls"])) and row("Recommended") == "yes: the one study of this trait and ancestry that passes the checks", lambda: (counted, row("Recommended")))
    shot("study-cases", full=False)
    stenosis = BY_ID["GCST90651074"]
    settled(stenosis["id"])
    check(f"[{name}] cases that are those of a record of the GWAS Catalog link to it",
          row("Cases and controls") == f"{stenosis['cases']:,} and {stenosis['controls']:,}: effective N {stenosis['neff']:,} · GCST90651074"
          and page.locator("#view-study .row", has_text="Cases and controls").locator("a").get_attribute("href") == "https://www.ebi.ac.uk/gwas/studies/GCST90651074", lambda: row("Cases and controls"))
    # h² on the liability scale, with the prevalence a paper states; or why there is none. UK Biobank, and how it is known
    adhd, proxy = BY_ID["PASS_ADHD_Demontis2018"], BY_ID["PASS_Alzheimers_Jansen2019"]
    settled(adhd["id"])
    liability, biobank = row("Liability-scale h²"), row("UK Biobank")
    settled(proxy["id"])
    check(f"[{name}] a binary study's page gives h² on the liability scale with the prevalence, quoted from its paper; a study of proxy cases has none; whether a study includes UK Biobank, and how it is known",
          liability.startswith(f"{plus_minus(adhd['h2l'], adhd['h2l_se'])} At a population prevalence of 5%: “") and "Demontis et al. 2023" in liability
          and row("Liability-scale h²") == "none: some of its cases are proxy cases, relatives of a case, and not cases of the trait"
          and biobank == "no · the paper Its paper's full text does not name UK Biobank." and row("UK Biobank").startswith("yes"), lambda: (liability, biobank, row("Liability-scale h²"), row("UK Biobank")))
    told = {}
    for study in ("UKB_460K.cancer_PROSTATE", "PASS_ProstateCancer", "PASS_BMI1", "PASS_Celiac", "UKB_460K.body_HEIGHTz"):
        settled(study)
        told[study] = (row("Recommended"), row("Cases and controls") if page.locator("#view-study .label", has_text="Cases and controls").count() else None)
        if study == "UKB_460K.cancer_PROSTATE":
            shot("study-cases-not-public", full=False)
    check(f"[{name}] a study whose cases are not public says so, and h² z sets it against the largest of the others; an older study is superseded; a trait of which nothing passes has no study to use",
          told["UKB_460K.cancer_PROSTATE"] == ("no: the study to use is Schumacher_2018_ProstateCancer", "not known: no effective N A yes-or-no trait by its name. Its cases are not public: the GWAS Catalog's records of the "
                                               "paper count individuals, and no other public record of them was found.")
          and told["PASS_ProstateCancer"][0] == "yes: of the 3 studies of this trait and ancestry that pass the checks, the largest h² z among the one of the largest effective N and the 1 whose cases are not known, which effective N cannot rank"
          and told["PASS_BMI1"] == ("no: superseded by Yengo_2018_BMI", None) and told["PASS_Celiac"][0] == "no: no study of this trait and ancestry passes the checks"
          and told["UKB_460K.body_HEIGHTz"] == ("no: superseded by Yengo_2022_Height", None), lambda: told)
    # the candidates for the next release, and the rule's backtest
    show("?tab=candidates", "#cand-steps tbody tr")
    table = lambda name: [[" ".join(cell.split()) for cell in row] for row in page.locator(f"#{name} tbody tr").evaluate_all("rows => rows.map((row) => [...row.cells].map((cell) => cell.innerText))")]
    slash, test = lambda date: date.replace("-", "/"), CANDIDATES["backtest"]
    paper = lambda p: f"{p['author']}, {p['journal']}, {slash(p['published'])}"
    listed = [[paper(p) if k == 0 else "", s["id"], s["trait"], f"{s['n']:,}", s["ancestry"]] for p in CANDIDATES["candidates"] for k, s in enumerate(p["studies"])]
    lead, aside = page.inner_text("#cand-lead"), page.inner_text("#cand-aside")
    links = page.locator("#cand-list a").evaluate_all("list => list.map((a) => a.href)")
    check(f"[{name}] candidates: the studies that the GWAS Catalog added since the release and that pass the rule, by paper, with their links; the papers set aside",
          table("cand-list") == listed and lead.startswith(f"{len(listed)} studies in {len(CANDIDATES['candidates'])} paper") and slash(CANDIDATES["until"]) in lead and slash(CANDIDATES["catalog"]["fetched"]) in lead
          and all(f"https://www.ebi.ac.uk/gwas/studies/{row[1]}" in links for row in listed) and f"{PUBMED}{CANDIDATES['candidates'][0]['pmid']}/" in links
          and all(paper(p) in aside and group["why"] in aside for group in CANDIDATES["aside"] for p in group["papers"])
          and page.locator("#view-candidates a", has_text="Nominate a study").get_attribute("href") == "https://github.com/gwassumstats/gwassumstats.github.io/issues/new?template=candidate.yml"
          and page.title() == "Candidates · GWAS Sumstats" and fits(), lambda: (table("cand-list"), lead, aside))
    last, said_test, missed = test["steps"][-1], page.inner_text("#cand-test"), page.locator("#cand-missed a").evaluate_all("list => list.map((a) => [a.textContent, a.getAttribute('href')])")
    check(f"[{name}] the backtest: what each condition of the rule leaves, what it finds of the studies that the release added, and the ones it misses by their short names",
          table("cand-steps") == [[("… and " if k else "") + step["rule"], f"{step['studies']:,}", f"{step['papers']:,}", f"{step['found']:,}", f"{step['found_papers']:,}"] for k, step in enumerate(test["steps"])]
          and f"It finds {last['found']} of the {test['taken']} studies ({100 * last['found'] / test['taken']:.1f}%) and {last['found_papers']} of the {test['taken_papers']} papers." in said_test
          and missed == [[BY_ID[m["id"]]["name"], f"?study={m['id']}"] for m in test["missed"]] and all(m["why"] in page.inner_text("#cand-missed") for m in test["missed"])
          and len(page.locator("#cand-rule li").all_inner_texts()) == 5, lambda: (table("cand-steps"), said_test, missed))
    shot("candidates")
    page.locator("#cand-beyond a").click()
    check(f"[{name}] what the rule finds beyond the release: counted, and the papers listed on demand",
          f"{test['beyond']['studies']} studies that the release did not add" in page.inner_text("#cand-beyond") and page.locator("#cand-more li").count() == len(test["beyond"]["new_papers"])
          and page.locator("#cand-more").is_visible() and page.inner_text("#cand-beyond a") == "hide the papers", lambda: page.inner_text("#cand-beyond"))
    # the studies by trait
    show("?tab=traits", "#traits-list tbody tr")
    of = {}
    for study in STUDIES:
        of.setdefault(study["of"], []).append(study)
    listed = sorted(of, key=lambda trait: (TRAITS[trait][1], trait.lower()))
    row = lambda trait: [trait, TRAITS[trait][1], TRAITS[trait][0] or "–", *(str(sum(study["anc"] == anc for study in of[trait])) for anc in RELEASE["ancestries"]),
                         ", ".join(study["name"] for anc in RELEASE["ancestries"] for study in of[trait] if "rec" in study and study["anc"] == anc) or "–", ", ".join(study["name"] for study in of[trait])]
    rows, lead = table("traits-list"), page.inner_text("#traits-lead")
    old = page.locator("#traits-list a.superseded", has_text="Speliotes_2010_BMI")
    check(f"[{name}] traits: every trait with its category, its term, the study to use of each ancestry and its studies by their short names, in the order of the categories",
          [r[:2] for r in rows] == [[trait, TRAITS[trait][1]] for trait in listed] and rows[0] == row(listed[0]) and row("Body mass index") in rows and len(rows) == RELEASE["traits"]["all"] == 335
          and row("Body mass index")[5] == "Yengo_2018_BMI, Akiyama_2017_BMIEAS" and row("Celiac disease")[5] == "–" and old.get_attribute("title") == "superseded by Yengo_2018_BMI"
          and page.locator("#traits-list a.superseded").count() == RELEASE["recommended"]["superseded"]
          and lead.startswith("The 602 studies are of 335 traits: 327 have a study of EUR ancestry, 49 have a study of EAS ancestry. 316 have a term") and page.title() == "Traits · GWAS Sumstats"
          and page.locator("#traits-list a", has_text="EFO_0004340").get_attribute("href") == "https://www.ebi.ac.uk/ols4/ontologies/efo/classes?short_form=EFO_0004340", lambda: (rows[:2], lead))
    page.fill("#traits-q", "sleep duration")
    check(f"[{name}] the search narrows the traits by their names, terms, categories and studies",
          [r[0] for r in table("traits-list")] == sorted((trait for trait in of if all(w in "\n".join([trait, *TRAITS[trait], *(study["name"] for study in of[trait])]).lower() for w in ("sleep", "duration"))),
                                                         key=lambda trait: (TRAITS[trait][1], trait.lower())) and 0 < page.locator("#traits-list tbody tr").count() < 10, lambda: table("traits-list"))
    shot("traits")
    page.fill("#traits-q", "no such trait")
    check(f"[{name}] a search that finds nothing says so", page.locator("#traits-none").is_visible() and page.locator("#traits-list").is_hidden())
    page.fill("#traits-q", "")
    page.locator("#traits-list a", has_text=of["Body mass index"][0]["name"]).first.click()
    expect(page.locator("#view-study .fields")).to_be_visible()
    naming = " ".join(page.locator("#view-study .row", has_text="Trait").first.inner_text().split())
    check(f"[{name}] a study's page names its trait, with its category and its term, and leads back to the trait's row",
          naming == "Trait Body mass index (anthropometric; EFO_0004340)" and page.url.endswith(f"?study={of['Body mass index'][0]['id']}"), lambda: naming)
    page.locator("#view-study .row", has_text="Trait").first.locator("a", has_text="Body mass index").click()
    expect(page.locator("#traits-list tbody tr").first).to_be_visible()
    check(f"[{name}] and from there the traits are searched for it", page.input_value("#traits-q") == "Body mass index" and "Body mass index" in [r[0] for r in table("traits-list")] and page.url.endswith("?tab=traits"),
          lambda: page.input_value("#traits-q"))
    unsigned = next(d for d in STUDIES if not d["cross"])
    show(f"?study={unsigned['id']}", "#view-study .section-head")
    check(f"[{name}] an unsigned study's page says it has no genetic correlations", "unsigned" in page.inner_text("#view-study .section-head"))
    for tab in ("notes", "about"):
        show(f"?tab={tab}", f"#view-{tab} h2")
        values = page.locator(f"#view-{tab} [data-v]").all_inner_texts()
        named_here = page.locator(f"#view-{tab} [data-study] a, #view-{tab} [data-v] a").evaluate_all("list => list.map((a) => [a.textContent, a.getAttribute('href')])")
        check(f"[{name}] {tab}: its numbers are filled in, a study it names is a link by its short name, and nothing private shows",
              all(values) and fits() and not private() and page.locator("#o2-release").is_hidden() and page.locator("#o2-folders").is_hidden()
              and len(named_here) > 0 and all(BY_ID[href.removeprefix("?study=")]["name"] == text for text, href in named_here), lambda: (values, named_here))
        shot(tab)
    show("", "#sumstats-table tbody")
    page.locator("#sumstats-table td.study a").first.click()
    expect(page.locator("#view-study h2")).to_have_text(STUDIES[0]["name"])
    page.click("nav a[data-tab=h2]")
    expect(page.locator("#h2-plot svg")).to_be_visible()
    page.go_back()
    expect(page.locator("#view-study h2")).to_have_text(STUDIES[0]["name"])
    page.go_back()
    expect(page.locator("#sumstats-table tbody").first).to_be_visible()
    check(f"[{name}] links change the address, and Back returns", page.evaluate("location.search") == "")
    show(f"?dataset={added}", "#view-study h2")
    check(f"[{name}] an old link, ?dataset=ID, opens the study's page as ?study=ID",
          page.evaluate("location.search") == f"?study={added}" and page.inner_text("#view-study h2") == BY_ID[added]["name"], lambda: page.url)
    check(f"[{name}] no error, and nothing is loaded from elsewhere", not errors and not elsewhere, lambda: (errors, elsewhere))
    browser.close()


with sync_playwright() as pw:
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), functools.partial(Quiet, directory=SITE))
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{server.server_port}/"
    if SHOTS:
        Path(SHOTS).mkdir(parents=True, exist_ok=True)
    run(pw, url, "light", viewport={"width": 1280, "height": 900})
    if SHOTS:  # the dark colours, for a look
        run(pw, url, "dark", viewport={"width": 1280, "height": 900}, color_scheme="dark")
    server.shutdown()
print("site check:", "FAILED" if failed else "passed")
sys.exit(1 if failed else 0)
