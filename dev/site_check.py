"""Browser check of the site, served from site/ on this machine.

A headless browser opens every tab: the tables, the filters and the sorting, the charts, a dataset's page, and the
login, which must show the paths only for the right passphrase, keep them over a reload and drop them at logout. What
the page shows is compared with site/data, read here independently.

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
from pathlib import Path

from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).parents[1]
SITE = ROOT / "site"
sys.path.insert(0, str(ROOT))
from build import crypt  # noqa: E402

PASSPHRASE = os.environ.get("SUMSTATS_PASSPHRASE") or sys.exit("set SUMSTATS_PASSPHRASE: the check logs in")
SHOTS = os.environ.get("SITE_CHECK_SHOTS")
DATA = json.loads((SITE / "data" / "datasets.json").read_text())
DATASETS, RELEASE = DATA["datasets"], DATA["release"]
BY_ID = {d["id"]: d for d in DATASETS}
RG = {anc: json.loads((SITE / "data" / f"rg_{anc}.json").read_text()) for anc in RELEASE["ancestries"]}
PATHS = json.loads(crypt.decrypt(json.loads((SITE / "data" / "paths.enc.json").read_text()), PASSPHRASE))
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
    """Every other dataset of a's ancestry with its estimate against a, the strongest first."""
    pairs = [{"id": b, **estimate(a, b)} for b in RG[BY_ID[a]["anc"]]["ids"] if b != a]
    return sorted(pairs, key=lambda r: -1 if r["rg"] is None else abs(r["rg"]), reverse=True)


def file_path(dataset, kind):
    file = PATHS["files"][dataset][kind]
    return file and f"{PATHS['folders'][file[0]]}/{file[1]}"


def text(d):
    return "\n".join([d["trait"], d["id"], d["origin"], d["ref"]]).lower()


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
    page.on("request", lambda r: r.url.startswith(url) or r.url.startswith("data:") or elsewhere.append(r.url))
    count = lambda: page.inner_text("#count")
    bodies = lambda table: page.locator(f"#{table} tbody")
    first = lambda table: bodies(table).first.locator("td.trait").inner_text()
    fits = lambda: page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    private = lambda: [f for f in [PATHS["release"], *PATHS["folders"]] if f in page.content()]

    def shot(view, full=True):
        if SHOTS:
            page.screenshot(path=f"{SHOTS}/{name}-{view}.png", full_page=full)

    def show(query, ready):
        page.goto(url + query)
        expect(page.locator(ready).first).to_be_visible(timeout=20000)

    # the datasets
    show("", "#sumstats-table tbody")
    check(f"[{name}] the table lists every dataset", bodies("sumstats-table").count() == len(DATASETS) == 602 and count() == "602 datasets", count)
    shot("sumstats", full=False)
    page.fill("#q", "height uk")
    found = [d for d in DATASETS if all(w in text(d) for w in ("height", "uk"))]
    check(f"[{name}] the search finds every word, in trait, ID, release or reference",
          0 < len(found) < 602 and bodies("sumstats-table").count() == len(found) and count() == f"{len(found)} of 602 datasets", count)
    page.fill("#q", "")
    page.select_option("#f-anc", "EAS")
    eas = [d for d in DATASETS if d["anc"] == "EAS"]
    ok = count() == "56 of 602 datasets"
    for box, keep in (("z6", lambda d: d["h2_z"] > 6), ("signed", lambda d: d["signed"]), ("panel", lambda d: "panel" in d)):
        page.check(f"#f-{box}")
        ok &= count() == f"{sum(map(keep, eas))} of 602 datasets"
        page.uncheck(f"#f-{box}")
    page.check("#f-panel")
    ok &= count() == "17 of 602 datasets"
    page.select_option("#f-anc", "EUR")
    ok &= count() == "85 of 602 datasets"
    page.uncheck("#f-panel")
    page.select_option("#f-anc", "")
    check(f"[{name}] the filters count as the data does: ancestry, h² z > 6, signed, nonredundant", ok and count() == "602 datasets", count)
    by_n = sorted(DATASETS, key=lambda d: d["n"])
    sort_n = page.locator("#sumstats-table th:nth-child(4) button")  # N
    order = []
    for _ in range(3):
        sort_n.click()
        order.append(first("sumstats-table"))
    check(f"[{name}] a column sorts largest first, then smallest first, then as released",
          order == [by_n[-1]["trait"], by_n[0]["trait"], DATASETS[0]["trait"]], lambda: order)
    check(f"[{name}] logged out, no path is in the page", not private() and page.locator(".paths").count() == 0, private)

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
    check(f"[{name}] the passphrase shows each dataset's files",
          page.locator("#sumstats-table .paths").count() == 602 and shown == [file_path(DATASETS[0]["id"], 0)] and page.inner_text("#login-open") == "logout"
          and page.locator("#login").is_hidden() and page.input_value("#login-pass") == "", lambda: shown)
    shot("sumstats-logged-in", full=False)
    page.locator("#sumstats-table tbody").first.locator(".paths a").click()
    check(f"[{name}] copy puts the path on the clipboard", page.evaluate("navigator.clipboard.readText()") == file_path(DATASETS[0]["id"], 0))
    page.reload()
    expect(page.locator("#sumstats-table .paths").first).to_be_visible(timeout=30000)
    check(f"[{name}] a reload stays logged in", page.inner_text("#login-open") == "logout")
    added = next(d["id"] for d in DATASETS if d["release"] == RELEASE["year"] and d["cross"])
    show(f"?dataset={added}", "#view-dataset .paths")
    shown = page.locator("#view-dataset .path").all_inner_texts()
    check(f"[{name}] an added dataset's page lists its h² file and its rg file",
          shown == [file_path(added, 0), file_path(added, 1)] and shown[0] != shown[1], lambda: shown)
    shot("dataset-logged-in")
    show("?tab=h2", "#h2-table .paths")
    check(f"[{name}] the heritability table lists the files too", page.locator("#h2-table .paths").count() == 602 and fits())
    show("?tab=notes", "#o2-release")
    check(f"[{name}] the release notes say where the release is", PATHS["release"] in page.inner_text("#o2-release") and fits())
    shot("notes-logged-in")
    show("?tab=about", "#o2-folders")
    check(f"[{name}] About lists the folders of the files", page.locator("#o2-folders li").count() == len(PATHS["folders"]) and fits(),
          lambda: page.inner_text("#o2-folders"))
    shot("about-logged-in")
    page.click("#login-open")
    expect(page.locator("#o2-folders")).to_be_hidden()
    show("", "#sumstats-table tbody")
    check(f"[{name}] logout hides the paths and forgets the key",
          not private() and page.locator(".paths").count() == 0 and page.evaluate("localStorage.getItem('key')") is None and page.inner_text("#login-open") == "login", private)

    # heritability
    show("?tab=h2", "#h2-plot svg")
    top = max(DATASETS, key=lambda d: d["h2"])
    bins = [round(d["h2"] * 1e4) // 200 for d in DATASETS]  # of 0.02
    mode = max(sorted(set(bins)), key=bins.count)
    check(f"[{name}] heritability: a column per bin of h², and the datasets by h²",
          page.inner_text("#h2-caption") == f"h² of the 602 datasets shown, in bins of 0.02: most often {mode * 0.02:.2f} to {(mode + 1) * 0.02:.2f}."
          and first("h2-table") == top["trait"] and page.locator("#h2-plot .bar").count() == len(set(bins)) and fits(),
          lambda: page.inner_text("#h2-caption"))
    shot("h2", full=False)
    box = page.locator("#h2-plot svg").bounding_box()
    page.mouse.move(box["x"] + box["width"] * 0.12, box["y"] + box["height"] * 0.5)
    check(f"[{name}] a column says how many datasets it stands for", "dataset" in page.inner_text("#tip") and "h² from" in page.inner_text("#tip"), lambda: page.inner_text("#tip"))
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
        ok &= seen[-1] == (f"{e['rg']:.4f} (SE {e['se']:.4f})" if e["rg"] is not None else f"not estimated: LDSC gave none for this pair ({e['note']})")
        ok &= page.input_value("#rg-a") == BY_ID[a]["trait"] and page.input_value("#rg-b") == BY_ID[b]["trait"]
    check(f"[{name}] a pair's genetic correlation is the one in the data, in either order, in EUR and EAS, or says why there is none", ok, lambda: seen)
    a = "UKB_460K.body_HEIGHTz"
    near = neighbours(a)
    strong = [r for r in near if r["rg"] is not None and abs(r["z"]) >= 3]
    show(f"?tab=rg&a={a}&b={strong[1]['id']}", "#neighbours-table tbody")
    listed = lambda: page.locator("#neighbours-table td.trait").all_inner_texts()
    ok = listed() == [BY_ID[r["id"]]["trait"] for r in strong[:10]] and fits()
    shot("rg")
    page.click("#neighbours-more")
    ok &= bodies("neighbours-table").count() == len(strong)
    page.uncheck("#rg-strict")
    page.click("#neighbours-more")
    ok &= bodies("neighbours-table").count() == len(near) == len(eur) - 1
    page.locator("#neighbours-table th:nth-child(5) button").click()  # z
    by_z = max((r for r in near if r["z"] is not None), key=lambda r: r["z"])
    ok &= listed()[0] == BY_ID[by_z["id"]]["trait"]
    check(f"[{name}] a dataset's genetic correlations: the strongest ten with |z| ≥ 3, then all, then every pair, sortable", ok, lambda: listed()[:10])
    show("?tab=rg", "#heat svg")
    page.fill("#rg-a", "schizo")
    options = page.locator("#rg-a-list a")
    expect(options.first).to_be_visible()
    matches = [d for d in DATASETS if d["cross"] and "schizo" in text(d)]
    ok = options.count() == len(matches)
    shot("rg-picker", full=False)
    page.keyboard.press("ArrowDown")
    page.keyboard.press("ArrowDown")
    page.keyboard.press("Enter")
    expect(page.locator("#neighbours-table tbody").first).to_be_visible(timeout=20000)
    check(f"[{name}] typing in a field lists the signed datasets that match; arrows and Enter pick one",
          ok and page.evaluate("new URLSearchParams(location.search).get('a')") == matches[1]["id"] and page.input_value("#rg-a") == matches[1]["trait"],
          lambda: (options.count(), len(matches), page.url))
    order = lambda anc, kind: [BY_ID[RG[anc]["ids"][i]] for i in RG[anc]["heat"][kind]["order"]]  # the build's clustering

    def named(rows):
        """The heatmap names these datasets beside its rows, in this order (a long name cut short)."""
        shown = page.locator("#heat text").all_text_contents()
        return len(shown) == len(rows) and all(d["trait"].startswith(s.rstrip("…")) for s, d in zip(shown, rows))

    def tip_for(a, b):
        e = estimate(a["id"], b["id"])
        return ("not estimated" if e["rg"] is None else f"{e['rg']:.4f} (SE {e['se']:.4f})") + f"\n{a['trait']} × {b['trait']}"

    def cell(rows, i, j):
        """The middle of the heatmap's cell in row i, column j, scrolled to the middle of the window."""
        def middle():
            box = page.locator("#heat image").bounding_box()
            return box["x"] + (j + 0.5) * box["width"] / len(rows), box["y"] + (i + 0.5) * box["height"] / len(rows)
        page.evaluate("y => scrollBy(0, y - innerHeight / 2)", middle()[1])
        return middle()

    show("?tab=rg", "#heat image")
    ok = named(order("EUR", "panel")) and page.locator("#heat .tree").count() == 1 and fits()
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
    check(f"[{name}] a panel's heatmap is in the order of the build's clustering; a cell names its pair and looks it up", ok, lambda: (tip, page.url))
    show("?tab=rg", "#heat image")
    ok = True
    for anc in RG:
        page.locator("#heat-tabs a", has_text=f"{anc} all").click()
        rows = order(anc, "all")
        expect(page.locator("#heat svg")).to_have_attribute("aria-label", re.compile(f"the {len(rows)} signed {anc} datasets"), timeout=20000)
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
    check(f"[{name}] the heatmap of all signed datasets of an ancestry, clustered too: a cell names its pair", ok, lambda: tip)

    # a dataset's page, the texts, the address
    show(f"?dataset={added}", "#view-dataset .pairs tbody")
    d, strong = BY_ID[added], [r for r in neighbours(added) if r["rg"] is not None and abs(r["z"]) >= 3]
    fields = dict(zip(page.locator("#view-dataset .label").all_inner_texts(), page.locator("#view-dataset .value").all_inner_texts()))
    links = page.locator("#view-dataset .fields a").evaluate_all("list => list.map((a) => a.href)")
    check(f"[{name}] a dataset's page: its fields, its paper and GWAS Catalog links, its strongest genetic correlations",
          page.inner_text("#view-dataset h2") == d["trait"] and fields["ID"] == added and fields["N"] == f"{d['n']:,}" and fields["Release"] == d["origin"]
          and fields["h²"] == f"{d['h2']:.4f} (SE {d['h2_se']:.4f}), z {d['h2_z']:.2f}" and d["paper"] in links and d["catalog"] in links
          and page.locator("#view-dataset .pairs td.trait").all_inner_texts() == [BY_ID[r["id"]]["trait"] for r in strong[:10]]
          and page.title() == f"{d['trait']} · GWAS Sumstats" and fits(), lambda: fields)
    shot("dataset")
    unsigned = next(d for d in DATASETS if not d["cross"])
    show(f"?dataset={unsigned['id']}", "#view-dataset .section-head")
    check(f"[{name}] an unsigned dataset's page says it has no genetic correlations", "unsigned" in page.inner_text("#view-dataset .section-head"))
    for tab in ("notes", "about"):
        show(f"?tab={tab}", f"#view-{tab} h2")
        values = page.locator(f"#view-{tab} [data-v]").all_inner_texts()
        check(f"[{name}] {tab}: its numbers are filled in, and nothing private shows", all(values) and fits() and not private()
              and page.locator("#o2-release").is_hidden() and page.locator("#o2-folders").is_hidden(), lambda: values)
        shot(tab)
    show("", "#sumstats-table tbody")
    page.locator("#sumstats-table td.trait a").first.click()
    expect(page.locator("#view-dataset h2")).to_have_text(DATASETS[0]["trait"])
    page.click("nav a[data-tab=h2]")
    expect(page.locator("#h2-plot svg")).to_be_visible()
    page.go_back()
    expect(page.locator("#view-dataset h2")).to_have_text(DATASETS[0]["trait"])
    page.go_back()
    expect(page.locator("#sumstats-table tbody").first).to_be_visible()
    check(f"[{name}] links change the address, and Back returns", page.evaluate("location.search") == "")
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
