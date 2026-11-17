// GWAS Sumstats: a catalog read from data/*.json (written by build/). Where the files are on O2 is in
// data/paths.enc.json, encrypted under the lab's passphrase; it is decrypted here, and nothing leaves this browser.
// A study is called by its short name (d.name) everywhere; its ID in the release shows on its page and in exports.

const TITLE = "GWAS Sumstats";
const TABS = { sumstats: "", h2: "Heritability", rg: "Genetic correlation", candidates: "Candidates", notes: "Release notes", about: "About" };  // tab → its page title
const STRONG_Z = 3;  // a genetic correlation is listed among a study's strongest if its |z| is at least this
const TOP = 10;      // the strongest shown before "show all"
const PANEL_MAX = Math.sqrt(0.1);  // a panel's heatmap: its colours end at the panels' bound, rg² = 0.1
const PUBMED = "https://pubmed.ncbi.nlm.nih.gov/", DOI = "https://doi.org/";  // where a paper's links lead
const ISSUES = "https://github.com/gwassumstats/gwassumstats.github.io/issues/new";  // where a problem is reported
const CATALOG = "https://www.ebi.ac.uk/gwas/studies/";  // where the GWAS Catalog shows a study
const QC = ["pass", "warn", "fail"];  // what the quality checks say of a study, the best first
const $ = (id) => document.getElementById(id);
const el = (tag, cls, text) => Object.assign(document.createElement(tag), { className: cls || "", textContent: text ?? "" });
const svg = (tag, attrs = {}) => {
  const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [name, value] of Object.entries(attrs)) node.setAttribute(name, value);
  return node;
};

const state = {
  loaded: false, studies: [], byId: new Map(), release: null,
  papers: {},         // PMID → a paper's label, "Speliotes et al. 2010 Nat Genet", and the day it was published
  view: "sumstats",   // a tab, or "study": the page of state.study
  study: null,
  filter: { q: "", anc: "", qc: "", z6: false, signed: false, panel: false, ticked: false },  // of the studies under "sumstats" and "heritability"
  cart: new Set(),    // the IDs of the studies ticked, for the exports; kept in this browser
  sort: { sumstats: null, h2: { key: "h2", dir: -1 }, neighbours: null },  // null: the release's order; for neighbours, the strongest first
  rg: { a: null, b: null, strict: true, all: false, heat: null },  // the pair looked up, the neighbours listed, the heatmap drawn ({ anc, kind })
  paths: null,        // where the files are, once logged in
};
const tables = new Map();  // ancestry → promise of its genetic correlations
const rows = { sumstats: new Map(), h2: new Map() };  // study id → its <tbody>, per table
let blob;  // promise of the encrypted paths
let records;  // promise of the papers: authors, title, journal, BibTeX entry
let found;  // promise of the candidates for the next release, and of the rule's backtest

// ---------------------------------------------------------------- helpers

// A short message at the bottom, gone after 5 s.
function toast(message) {
  $("toast").replaceChildren(message);
  clearTimeout(toast.timer);
  toast.timer = setTimeout(() => $("toast").replaceChildren(), 5000);  // empty, not hidden: screen readers announce each new message
}

// A link that acts in the page.
const link = (handler) => (event) => {
  event.preventDefault();
  handler(event);
};
const action = (text, handler, cls = "dashed") => {
  const a = Object.assign(el("a", cls, text), { href: "#" });
  a.addEventListener("click", link(handler));
  return a;
};

// localStorage, which may be off: read with one argument, write with two (null removes)
function store(key, value) {
  try {
    if (value === undefined) return localStorage.getItem(key);
    value === null ? localStorage.removeItem(key) : localStorage.setItem(key, value);
  } catch {}
  return null;
}

const int = (n) => n.toLocaleString("en-US");
const fixed = (digits) => (x) => (x == null ? "–" : x.toFixed(digits));
// A number to two significant digits, written out in full: 0.146 → 0.15, 16.2 → 16, −0.00824 → −0.0082
function sig(x) {
  if (x == null) return "–";
  if (x === 0) return "0";
  const [mantissa, exponent] = x.toExponential(1).split("e");  // rounded here: 0.0996 → 1.0e-1
  return Number(`${mantissa}e${exponent}`).toFixed(Math.max(0, 1 - Number(exponent))).replace("-", "−");
}
// An estimate and its standard error, each to two significant digits: as text, and as a table cell whose ± signs
// stand one under the other
const plusMinus = (value, se) => (value == null ? "–" : `${sig(value)} ± ${sig(se)}`);
const plusMinusCell = (value, se) => (value == null ? ["–"] : [el("span", "est", sig(value)), " ± ", el("span", "se", sig(se))]);
const plusMinusLabel = (name) => () => [el("span", "est", name), " ± ", el("span", "se", "SE")];
const pValue = (p) => (p == null ? "–" : p === 0 ? "<1e-300" : p < 1e-4 ? p.toExponential(2) : p.toFixed(4));  // LDSC prints 0 below about 1e-308
const yesNo = (flag) => (flag ? "yes" : "no");
const howMany = (n, one, many = `${one}s`) => `${int(n)} ${n === 1 ? one : many}`;
const words = (text) => text.toLowerCase().split(/\s+/).filter(Boolean);
const httpUrl = (url) => (/^https:\/\/\S+$/.test(url ?? "") ? url : null);
const rgLabel = () => {
  const label = el("span", "", "r");
  label.append(el("sub", "", "g"));
  return label;
};
const rgPlusMinusLabel = () => {
  const [name, ...rest] = plusMinusLabel("")();
  name.append(rgLabel());
  return [name, ...rest];
};
// A long name in parts (PASS_Height_Yengo2022, NG_Kanai_2018_BBJ): it may wrap after an underscore
const parts = (text) => text.split("_").flatMap((part, k, all) => (k < all.length - 1 ? [`${part}_`, el("wbr")] : [part]));
const outside = (url, ...content) => {
  const a = Object.assign(el("a"), { href: url, target: "_blank", rel: "noopener noreferrer" });
  a.append(...content);
  return a;
};
const studyLink = (d) => {
  const a = Object.assign(el("a"), { href: `?study=${encodeURIComponent(d.id)}` });
  a.append(...parts(d.name));
  return a;
};
// A finding of the quality checks names another study as {ID}: its parts, that study as its short name, linked or not
const finding = (text, linked) => text.split(/\{([^}]+)\}/).map((part, k) => (k % 2 && state.byId.has(part) ? (linked ? studyLink(state.byId.get(part)) : state.byId.get(part).name) : part));
// What the checks say of a study, in a word; the fails and warns show when the pointer rests on it
const verdict = (d) => {
  const word = el("span", `qc ${d.qc}`, d.qc), reasons = (d.notes ?? []).filter(([level]) => level !== "note");
  if (reasons.length) word.title = reasons.map(([level, text]) => `${level}: ${finding(text).join("")}`).join("\n");
  return word;
};
// A study's paper, by its label, linked to PubMed; the inventory's own words for a study whose paper is not known
const reference = (d) => (d.cite ? outside(`${PUBMED}${d.cite[0]}/`, d.label) : parts(d.ref));
// The order of two studies by their papers: the earlier published first; then by label, which begins with the first
// author; then by short name. Plain comparisons: a day is written 2010-10-10, or 2018-07 where PubMed gives no day.
const plainly = (x, y) => (x < y ? -1 : x > y ? 1 : 0);
const byPublication = (a, b) => plainly(a.published, b.published) || plainly(a.label, b.label) || plainly(a.name, b.name);

const load = (name) => fetch(`data/${name}`, { cache: "no-cache" }).then((r) => (r.ok ? r.json() : Promise.reject(new Error(`${name}: ${r.status}`))));
const papers = () => (records ??= load("papers.json").catch((e) => { records = undefined; throw e; }));  // loaded when first needed

// Hands the reader a file made here
function download(name, text, type) {
  const a = Object.assign(el("a"), { href: URL.createObjectURL(new Blob([text], { type })), download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 10000);
}

// ---------------------------------------------------------------- genetic correlations

// An ancestry's table, loaded when first needed: `ids`, and one list per column over the pairs i < j of `ids`, row by
// row: (0,1), (0,2), …, (1,2), … (build/release.py)
function correlations(anc) {
  if (!tables.has(anc)) {
    const table = load(`rg_${anc}.json`).then((t) => Object.assign(t, { n: t.ids.length, at: new Map(t.ids.map((id, i) => [id, i])) }));
    tables.set(anc, table);
    table.catch(() => tables.delete(anc));  // a failed load is tried again at the next need
  }
  return tables.get(anc);
}

const position = (t, i, j) => (i > j ? position(t, j, i) : i * (2 * t.n - i - 1) / 2 + j - i - 1);  // of the pair of studies i and j in the lists

function estimate(t, i, j) {
  const k = position(t, i, j);
  return { rg: t.rg[k], se: t.se[k], z: t.z[k], p: t.p[k], gcov: t.gcov[k], gcov_se: t.gcov_se[k], note: t.notes[k] };
}

// Every other study of d's ancestry with its estimate against d
function neighbours(t, d) {
  const i = t.at.get(d.id);
  return t.ids.flatMap((id, j) => (j === i ? [] : [{ d: state.byId.get(id), ...estimate(t, i, j) }]));
}
const strong = (pairs) => pairs.filter((r) => r.rg != null && Math.abs(r.z) >= STRONG_Z);
const strongest = (x, y) => (y.rg == null ? -1 : Math.abs(y.rg)) - (x.rg == null ? -1 : Math.abs(x.rg));  // by |rg|; pairs without an estimate last

// ---------------------------------------------------------------- login: the paths

const bytes = (base64) => Uint8Array.from(atob(base64), (c) => c.charCodeAt(0));
const encrypted = () => (blob ??= load("paths.enc.json").catch((e) => { blob = undefined; throw e; }));

// PBKDF2 turns the passphrase into the AES key; its many rounds are what make guessing slow (build/crypt.py)
async function deriveKey(passphrase, { salt, iterations }) {
  const material = await crypto.subtle.importKey("raw", new TextEncoder().encode(passphrase), "PBKDF2", false, ["deriveBits"]);
  return new Uint8Array(await crypto.subtle.deriveBits({ name: "PBKDF2", hash: "SHA-256", salt: bytes(salt), iterations }, material, 256));
}

// Throws an OperationError for a wrong key
async function decrypt({ iv, data }, key) {
  const aes = await crypto.subtle.importKey("raw", key, "AES-GCM", false, ["decrypt"]);
  return JSON.parse(new TextDecoder().decode(await crypto.subtle.decrypt({ name: "AES-GCM", iv: bytes(iv) }, aes, bytes(data))));
}

function setPaths(paths) {
  state.paths = paths;
  for (const cache of Object.values(rows)) cache.clear();  // the rows gain or lose their paths
  renderCart();  // and the export its paths
  render();
}

async function login(event) {
  event.preventDefault();
  const passphrase = $("login-pass").value.trim(), note = $("login-note"), submit = $("login").querySelector("[type=submit]");
  if (!passphrase) return $("login-pass").focus();
  note.classList.remove("error");
  note.textContent = "Checking…";
  submit.disabled = true;
  try {
    const sealed = await encrypted(), key = await deriveKey(passphrase, sealed);
    const paths = await decrypt(sealed, key);
    store("key", JSON.stringify({ salt: sealed.salt, key: btoa(String.fromCharCode(...key)) }));  // this browser stays logged in
    closeLogin();
    setPaths(paths);
  } catch (e) {
    note.textContent = e.name === "OperationError" ? "Wrong passphrase." : "Could not check the passphrase; reload and try again.";
    note.classList.add("error");
    $("login-pass").select();
  } finally {
    submit.disabled = false;
  }
}

// The key saved at login opens the paths again; a key from before the paths were re-encrypted does not, and is dropped.
async function restore() {
  let saved = null;
  try { saved = JSON.parse(store("key")); } catch {}
  if (!saved?.key) return state.paths && setPaths(null);
  const sealed = await encrypted().catch(() => null);
  if (!sealed) return;  // no connection: still logged in, the paths show at the next load
  try {
    if (saved.salt !== sealed.salt) throw new Error("the paths were encrypted again");
    setPaths(await decrypt(sealed, bytes(saved.key)));
  } catch {
    store("key", null);
    if (state.paths) setPaths(null);
  }
}

function logout() {
  store("key", null);
  $("login-pass").value = "";
  setPaths(null);
}

function openLogin(open = $("login").hidden) {
  if (!open) return closeLogin();
  $("login").hidden = false;
  $("login-pass").focus();
  renderBar();
}

function closeLogin() {
  $("login").hidden = true;
  $("login-pass").value = "";
  $("login-note").textContent = "";
  $("login-note").classList.remove("error");
  renderBar();
}

// Where a study's file is on O2; null when logged out
function filePath(d) {
  const file = state.paths?.files[d.id];
  return file ? `${state.paths.folders[file[0]]}/${file[1]}` : null;
}

function copy(text) {
  navigator.clipboard.writeText(text).then(() => toast("Copied."), () => toast("Could not copy; select the path instead."));
}

// A study's path with "copy", after `label`; null when logged out
function pathLine(d, label = "") {
  const path = filePath(d);
  if (!path) return null;
  const line = el("div", "paths");
  line.append(label, el("span", "path", path), action("copy", () => copy(path)));
  return line;
}

// ---------------------------------------------------------------- the studies: filters, tables

function shown() {
  const f = state.filter, sought = words(f.q);
  return state.studies.filter((d) => (!f.anc || d.anc === f.anc) && (!f.qc || d.qc === f.qc) && (!f.z6 || d.h2_z > 6) && (!f.signed || d.signed) && (!f.panel || d.panel)
    && (!f.ticked || state.cart.has(d.id)) && sought.every((w) => d.text.includes(w)));
}

const compare = (x, y) => (typeof x === "string" ? x.localeCompare(y, "en", { numeric: true, sensitivity: "base" }) : (x ?? -Infinity) - (y ?? -Infinity));
// `list` by a column, in the direction asked: by `order`, the column's own, or else by its values; ties, and no sort
// at all, keep the order of `list`
const sorted = (list, sort, value = (item) => item[sort.key], order = (a, b) => compare(value(a), value(b))) => (sort ? [...list].sort((a, b) => sort.dir * order(a, b)) : list);

// ---------------------------------------------------------------- the cart: studies ticked, and what is exported of them

const ticked = () => state.studies.filter((d) => state.cart.has(d.id));
// The ticks kept in this browser, without those of studies that are gone
function readCart() {
  let ids = [];
  try { ids = JSON.parse(store("cart")) ?? []; } catch {}
  state.cart = new Set(Array.isArray(ids) ? ids.filter((id) => state.byId.has(id)) : []);
}
let tickAll;  // the box in the header of the table shown: ticks every study listed

function tick(studies, on) {
  for (const d of studies) on ? state.cart.add(d.id) : state.cart.delete(d.id);
  store("cart", JSON.stringify([...state.cart]));
  for (const cache of Object.values(rows)) {  // the same study's box in the other table
    for (const d of studies) {
      const box = cache.get(d.id)?.querySelector(".tick input");
      if (box) box.checked = on;
    }
  }
  if (state.filter.ticked && (state.view === "sumstats" || state.view === "h2")) renderList();  // the rows listed change
  else renderCart();
}

// A study's box. Ticking leaves the table as it is, so the keyboard stays on the row.
function tickBox(d) {
  const box = Object.assign(el("input"), { type: "checkbox", checked: state.cart.has(d.id) });
  box.setAttribute("aria-label", `Tick ${d.name}`);
  box.addEventListener("change", () => tick([d], box.checked));
  return box;
}

function tickAllBox() {
  tickAll = Object.assign(el("input"), { type: "checkbox", title: "Tick every study listed" });
  tickAll.setAttribute("aria-label", "Tick every study listed");
  tickAll.addEventListener("change", () => tick(shown(), tickAll.checked));
  return tickAll;
}

// How many are ticked, and the exports; the header's box says whether the studies listed are
function renderCart() {
  const n = state.cart.size, cart = $("cart");
  cart.hidden = !n;
  if (n) {
    cart.replaceChildren(` · ${int(n)} ticked: download `, action(state.paths ? "IDs and paths (TSV)" : "IDs (TSV)", exportStudies), " or ",
      action("references (BibTeX)", exportReferences), " · ", action("untick all", () => tick(ticked(), false)));
  }
  if (!tickAll?.isConnected) return;
  const list = shown(), on = list.filter((d) => state.cart.has(d.id)).length;
  tickAll.checked = on > 0 && on === list.length;
  tickAll.indeterminate = on > 0 && on < list.length;
}

// The studies ticked as a table: short name, ID, and what a script needs to find and cite each; the path after login
function exportStudies() {
  const columns = [["short_name", (d) => d.name], ["id", (d) => d.id], ["release_name", (d) => d.trait], ["ancestry", (d) => d.anc], ["n", (d) => d.n],
    ["reference", (d) => d.label], ["published", (d) => d.published], ["pmid", (d) => (d.cite ?? []).join(";")], ...(state.paths ? [["path", filePath]] : [])];
  const cell = (value) => String(value ?? "").replace(/[\t\r\n]+/g, " ");
  const lines = [columns.map(([name]) => name), ...ticked().map((d) => columns.map(([, value]) => cell(value(d))))];
  download("gwas-sumstats-studies.tsv", `${lines.map((line) => line.join("\t")).join("\n")}\n`, "text/tab-separated-values");
}

// The papers of the studies ticked, as BibTeX: an entry a paper, under a comment that names its studies
async function exportReferences() {
  let all;
  try {
    all = await papers();
  } catch {
    return toast("The references did not load; check the connection and try again.");
  }
  const cited = new Map(), named = (d) => `${d.name} (${d.id})`;
  for (const d of ticked()) for (const pmid of d.cite ?? []) cited.set(pmid, [...(cited.get(pmid) ?? []), named(d)]);
  const entries = [...cited].sort(([a], [b]) => plainly(all[a].published, all[b].published) || plainly(all[a].key, all[b].key))  // the earliest paper first
    .map(([pmid, studies]) => `% ${studies.join(", ")}\n${all[pmid].bib}`);
  const without = ticked().filter((d) => !d.cite).map(named);
  download("gwas-sumstats-references.bib", entries.join("\n") + (without.length ? `\n% no paper known: ${without.join(", ")}\n` : ""), "application/x-bibtex");
}

// A column: key (the field it shows and sorts by), label, title (what the label means), num (a format: right-aligned,
// largest first), show (the cell's content, if not the field itself), cls, first (-1: a click sorts the largest first),
// fixed (it does not sort).
const TICK = { key: "ticked", label: tickAllBox, cls: "tick", fixed: true, show: tickBox };
const STUDY = { key: "name", label: "Study", title: "Short name: first author, year and trait", cls: "study", show: studyLink };
const ANCESTRY = { key: "anc", label: "Ancestry", cls: "short" };
const H2 = { key: "h2", label: plusMinusLabel("h²"), title: "Heritability and its standard error: S-LDSC, observed scale", cls: "pm", first: -1, show: (d) => plusMinusCell(d.h2, d.h2_se) };
const H2_Z = { key: "h2_z", label: "h² z", title: "h² over its standard error", num: sig };
const CHECKS = { key: "rank", label: "QC", title: "Quality checks: pass, warn or fail; a study's page says why", cls: "short", first: -1, show: verdict };
const COLUMNS = {
  sumstats: [TICK, STUDY, ANCESTRY,
    { key: "n", label: "N", title: "Sample size", num: int },
    { key: "m", label: "M", title: "SNPs in the file", num: int },
    H2, H2_Z, CHECKS,
    { key: "label", label: "Reference", title: "The study's paper; it sorts by the day the paper was published", cls: "ref", show: reference, order: byPublication }],
  h2: [TICK, STUDY, ANCESTRY, H2, H2_Z,
    { key: "intercept", label: "Intercept", title: "Intercept of the LD score regression", num: fixed(4) }, CHECKS],
  neighbours: [{ ...STUDY, show: (r) => studyLink(r.d) },
    { key: "rg", label: rgPlusMinusLabel, title: "Genetic correlation and its standard error", cls: "pm", first: -1, show: (r) => plusMinusCell(r.rg, r.se) },
    { key: "z", label: "z", num: sig },
    { key: "p", label: "p", num: pValue }],
};

const labelOf = (c) => [typeof c.label === "function" ? c.label() : c.label].flat();

// The header row. With `onSort`, a click on a column sorts by it, a second turns the order round, a third returns to
// the first order.
function header(columns, sort, onSort) {
  const tr = el("tr");
  for (const c of columns) {
    const th = el("th", c.cls || (c.num ? "num" : "")), first = c.first ?? (c.num ? -1 : 1);
    th.scope = "col";
    if (onSort && !c.fixed) {
      const button = el("button");
      button.type = "button";
      button.append(...labelOf(c));
      button.title = `${c.title ? `${c.title}. ` : ""}Sort by this column.`;
      if (sort?.key === c.key) {
        th.setAttribute("aria-sort", sort.dir > 0 ? "ascending" : "descending");
        button.append(sort.dir > 0 ? " ▲" : " ▼");
      }
      button.addEventListener("click", () => onSort(sort?.key !== c.key ? { key: c.key, dir: first } : sort.dir === first ? { key: c.key, dir: -first } : null));
      th.append(button);
    } else {
      th.append(...labelOf(c));
      if (c.title) th.title = c.title;
    }
    tr.append(th);
  }
  const thead = el("thead");
  thead.append(tr);
  return thead;
}

// One <tbody> per row: its cells and, under them, `below` (a study's path)
function body(columns, item, below) {
  const tbody = el("tbody"), tr = el("tr");
  for (const c of columns) {
    const td = el("td", c.cls || (c.num ? "num" : ""));
    td.append(...[c.show ? c.show(item) : c.num ? c.num(item[c.key]) : item[c.key]].flat());
    tr.append(td);
  }
  tbody.append(tr);
  if (below) {
    const sub = el("tr", "sub"), td = el("td");
    td.colSpan = columns.length;
    td.append(below);
    sub.append(td);
    tbody.append(sub);
  }
  return tbody;
}

// A table's only row when it has nothing to list
function nothing(columns, text) {
  const tbody = el("tbody"), tr = el("tr"), td = el("td", "empty", text);
  td.colSpan = columns.length;
  tr.append(td);
  tbody.append(tr);
  return tbody;
}

function renderList() {
  const view = state.view, columns = COLUMNS[view], list = shown(), all = state.studies.length;
  $("count").textContent = list.length === all ? howMany(all, "study", "studies") : `${int(list.length)} of ${howMany(all, "study", "studies")}`;
  const cache = rows[view];
  const by = columns.find((c) => c.key === state.sort[view]?.key)?.order;
  const bodies = sorted(list, state.sort[view], undefined, by).map((d) => cache.get(d.id) ?? cache.set(d.id, body(columns, d, pathLine(d))).get(d.id));
  const head = header(columns, state.sort[view], (sort) => { state.sort[view] = sort; renderList(); });
  $(`${view}-table`).replaceChildren(head, ...(bodies.length ? bodies : [nothing(columns, "No study matches.")]));
  renderCart();
  if (view === "h2") drawHistogram(list);
}

// ---------------------------------------------------------------- charts

function showTip(event, ...content) {
  const tip = $("tip");
  tip.replaceChildren(...content);
  tip.hidden = false;
  const { width, height } = tip.getBoundingClientRect();
  tip.style.left = `${Math.max(4, Math.min(event.clientX + 14, innerWidth - width - 4))}px`;
  tip.style.top = `${event.clientY + 18 + height > innerHeight ? event.clientY - height - 10 : event.clientY + 18}px`;
}
const hideTip = () => { $("tip").hidden = true; };

// How many studies have each h²: a column per bin of 0.02, in a frame that stays as the filters change
function drawHistogram(list) {
  const figure = $("h2-plot"), width = figure.clientWidth, height = 190, left = 30, right = 8, top = 8, bottom = 22;
  if (!width) return;
  const bin = (d) => Math.floor(Math.round(d.h2 * 1e4) / 200);  // h² has four decimals: in whole numbers it bins exactly
  const all = state.studies.map(bin), first = Math.min(...all), bins = Math.max(...all) - first + 1, w = (width - left - right) / bins;
  const counts = new Array(bins).fill(0);
  for (const d of list) counts[bin(d) - first]++;
  const peak = Math.max(1, ...counts), tick = [1, 2, 5, 10, 20, 50, 100, 200].find((t) => peak / t <= 5) ?? 500;
  const y = (count) => top + (height - top - bottom) * (1 - count / (Math.ceil(peak / tick) * tick));
  const range = (k) => `${((first + k) * 0.02).toFixed(2)} to ${((first + k + 1) * 0.02).toFixed(2)}`;
  const chart = svg("svg", { width, height, viewBox: `0 0 ${width} ${height}`, role: "img" });
  for (let count = 0; count <= peak; count += tick) {  // a hairline and its count, every `tick` studies
    chart.append(svg("line", { class: "grid", x1: left, x2: width - right, y1: y(count) + 0.5, y2: y(count) + 0.5 }));
    chart.append(Object.assign(svg("text", { x: left - 6, y: y(count), "text-anchor": "end", "dominant-baseline": "central" }), { textContent: count }));
  }
  for (let k = -first; k < bins; k += 10) {  // h² under the columns, every 0.2
    chart.append(Object.assign(svg("text", { x: left + k * w, y: height - 6, "text-anchor": "middle" }), { textContent: ((first + k) * 0.02).toFixed(1) }));
  }
  const bars = counts.map((count, k) => {
    const thick = Math.min(24, w - 2), x0 = left + k * w + (w - thick) / 2, r = Math.min(4, thick / 2, y(0) - y(count));  // a 2px gap; the top corners rounded
    const bar = svg("path", { class: "bar", d: `M${x0},${y(0)}V${y(count) + r}q0,${-r} ${r},${-r}h${thick - 2 * r}q${r},0 ${r},${r}V${y(0)}z` });
    if (count) chart.append(bar);
    return bar;
  });
  chart.addEventListener("pointermove", (event) => {
    const k = Math.floor((event.clientX - chart.getBoundingClientRect().left - left) / w);
    bars.forEach((bar, b) => bar.classList.toggle("on", b === k));
    if (k < 0 || k >= bins) return hideTip();
    showTip(event, el("b", "", howMany(counts[k], "study", "studies")), ` with h² from ${range(k)}`);
  });
  chart.addEventListener("pointerleave", () => { hideTip(); bars.forEach((bar) => bar.classList.remove("on")); });
  const label = list.length ? `h² of the ${howMany(list.length, "study", "studies")} shown, in bins of 0.02: most often ${range(counts.indexOf(peak))}.` : "";
  chart.setAttribute("aria-label", label || "No study matches.");
  figure.replaceChildren(chart);
  $("h2-caption").textContent = label;
}

// The colours of a genetic correlation, as [red, green, blue]: `of(v)` for v from -1 (--neg, blue) through 0 (--subtle)
// to 1 (--pos, red), and `self` for a study against itself. Read from the style sheet, so they follow the light or
// dark mode.
function shades() {
  const css = getComputedStyle(document.documentElement);
  const [neg, zero, pos, self] = ["--neg", "--subtle", "--pos", "--border"].map((name) => {
    const hex = css.getPropertyValue(name).trim();
    return [1, 3, 5].map((k) => parseInt(hex.slice(k, k + 2), 16));
  });
  return { self, of: (v) => zero.map((c, k) => Math.round(c + ((v < 0 ? neg : pos)[k] - c) * Math.min(1, Math.abs(v)))) };
}

const heatShown = () => state.rg.heat ?? { anc: state.rg.a?.anc ?? state.release.ancestries[0], kind: "panel" };

// A heatmap of the genetic correlations among the studies of a nonredundant panel ("panel") or among all signed
// studies of an ancestry ("all"): red above 0, blue below. The build ordered the studies by hierarchical clustering
// (build/release.py); the tree of that clustering is on the left, the leaves beside their rows.
async function drawHeat() {
  const figure = $("heat"), shown = heatShown(), { anc, kind } = shown;
  $("heat-tabs").replaceChildren(...state.release.ancestries.flatMap((name) => [["panel", `${name} nonredundant (${state.release[name].panel})`], ["all", `${name} all (${state.release[name].signed})`]]
    .map(([k, label]) => {
      const current = name === anc && k === kind, tab = action(label, () => { state.rg.heat = { anc: name, kind: k }; drawHeat(); }, current ? "current" : "");
      if (current) tab.setAttribute("aria-current", "true");
      return tab;
    })).flatMap((tab, k) => (k ? [" | ", tab] : [tab])));
  const drawn = `${anc} ${kind} ${figure.clientWidth}`;
  if (figure.dataset.drawn === drawn) return;  // as it stands
  let t;
  try {
    t = await correlations(anc);
  } catch {
    return figure.replaceChildren(el("p", "empty", "The genetic correlations did not load; check the connection and reload."));
  }
  const width = figure.clientWidth, now = heatShown();
  if (state.view !== "rg" || now.anc !== anc || now.kind !== kind || !width) return;  // another heatmap was asked for meanwhile
  figure.dataset.drawn = drawn;
  const { order, merges } = t.heat[kind], n = order.length, studies = order.map((i) => state.byId.get(t.ids[i]));
  const limit = kind === "panel" ? PANEL_MAX : 1, colours = shades();
  const TREE = 90, names = n <= 100 ? 270 : 0;  // the tree on the left; the studies' names on the right, where a row is tall enough to read one
  const cell = Math.max(2, Math.min(28, Math.floor((width - TREE - names) / n))), size = cell * n;

  const canvas = Object.assign(document.createElement("canvas"), { width: n, height: n }), pixels = canvas.getContext("2d").createImageData(n, n);
  for (let i = 0; i < n; i++) {  // a pixel per pair; the chart stretches the picture
    for (let j = 0; j < n; j++) {
      const [r, g, b] = i === j ? colours.self : colours.of((t.rg[position(t, order[i], order[j])] ?? 0) / limit);
      pixels.data.set([r, g, b, 255], 4 * (i * n + j));
    }
  }
  canvas.getContext("2d").putImageData(pixels, 0, 0);
  const chart = svg("svg", { class: "heat", width: TREE + size + names, height: size, viewBox: `0 0 ${TREE + size + names} ${size}`, role: "img",
    "aria-label": `Genetic correlations among ${kind === "panel" ? `the ${n} studies of the ${anc} nonredundant panel` : `the ${n} signed ${anc} studies`}, ordered by hierarchical clustering` });
  chart.append(svg("image", { x: TREE, y: 0, width: size, height: size, preserveAspectRatio: "none", href: canvas.toDataURL() }));

  const heights = merges.map((m) => m[2]), low = Math.min(...heights), high = Math.max(...heights);
  const reach = (height) => TREE - 4 - (TREE - 8) * (high > low ? (height - low) / (high - low) : 0);  // the root at the left edge
  const nodes = order.map((_, k) => [TREE, (k + 0.5) * cell]);  // where each leaf, then each merge, stands
  const branches = merges.map(([i, j, height]) => {
    const x = reach(height);
    nodes.push([x, (nodes[i][1] + nodes[j][1]) / 2]);
    return `M${nodes[i][0].toFixed(1)},${nodes[i][1].toFixed(1)}H${x.toFixed(1)}V${nodes[j][1].toFixed(1)}H${nodes[j][0].toFixed(1)}`;
  });
  chart.append(svg("path", { class: "tree", d: branches.join("") }));
  if (names) {
    const font = Math.min(12, cell), fits = Math.floor((names - 8) / (font * 0.6));
    for (const [k, d] of studies.entries()) {
      const name = svg("text", { x: TREE + size + 8, y: (k + 0.5) * cell, "dominant-baseline": "central" });
      name.textContent = d.name.length > fits ? `${d.name.slice(0, fits - 1)}…` : d.name;
      name.style.fontSize = `${font}px`;
      chart.append(name);
    }
  }
  const frame = chart.appendChild(svg("rect", { class: "on", width: cell, height: cell, visibility: "hidden" }));  // around the cell under the pointer

  const under = (event) => {
    const box = chart.getBoundingClientRect(), scale = box.width / (TREE + size + names);
    const j = Math.floor(((event.clientX - box.left) / scale - TREE) / cell), i = Math.floor((event.clientY - box.top) / scale / cell);
    return i >= 0 && j >= 0 && i < n && j < n ? [i, j] : null;
  };
  const leave = () => { hideTip(); frame.setAttribute("visibility", "hidden"); };
  chart.addEventListener("pointermove", (event) => {
    const here = under(event);
    if (!here) return leave();
    const [i, j] = here, e = i === j ? null : estimate(t, order[i], order[j]);
    for (const [name, value] of [["x", TREE + j * cell], ["y", i * cell], ["visibility", "visible"]]) frame.setAttribute(name, value);
    showTip(event, ...(e ? [el("b", "", e.rg == null ? "not estimated" : plusMinus(e.rg, e.se)), el("br"), `${studies[i].name} × ${studies[j].name}`] : [studies[i].name]));
  });
  chart.addEventListener("pointerleave", leave);
  chart.addEventListener("click", (event) => {  // a cell's pair, looked up above
    const here = under(event);
    if (!here || here[0] === here[1]) return;
    leave();
    state.rg.heat = shown;  // the heatmap stays, whatever the pair
    go({ tab: "rg", a: studies[here[0]].id, b: studies[here[1]].id });
    $("pair").scrollIntoView({ block: "start" });
  });
  figure.replaceChildren(chart);

  const scale = el("span", "scale"), steps = el("span", "steps");
  steps.append(...[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1].map((v) => { const swatch = el("i"); swatch.style.background = `rgb(${colours.of(v).join(" ")})`; return swatch; }));
  scale.append(rgLabel(), ` −${limit.toFixed(2)}`, steps, `+${limit.toFixed(2)}`);
  $("heat-caption").replaceChildren(scale, ` Blue below 0, red above${kind === "panel" ? "; no pair of a panel is beyond these" : ""}. A pair without an estimate shows as 0.`
    + " Rows and columns are in the order of the tree: hierarchical clustering on 1 − |", rgLabel(), "|, average linkage. A click on a cell looks its pair up.");
}

// ---------------------------------------------------------------- genetic correlation: a pair, a study's neighbours

// A text field that picks a study: typing lists the matches of `among()` under it; a click, or the arrows and
// Enter, picks one; an emptied field picks none.
function picker(input, list, among, picked, onPick) {
  let found = [], active = -1;
  const close = () => {
    list.hidden = true;
    input.setAttribute("aria-expanded", "false");
    input.removeAttribute("aria-activedescendant");
    active = -1;
  };
  const pick = (d) => {
    close();
    if (d) input.value = d.name;
    onPick(d);
  };
  const mark = (k) => {
    active = k;
    [...list.children].forEach((option, i) => option.classList.toggle("active", i === k));
    k < 0 ? input.removeAttribute("aria-activedescendant") : input.setAttribute("aria-activedescendant", list.children[k].id);
    list.children[k]?.scrollIntoView({ block: "nearest" });
  };
  const open = () => {
    const sought = words(input.value);
    found = among().filter((d) => sought.every((w) => d.text.includes(w))).slice(0, 60);
    list.replaceChildren(...(found.length ? found.map((d, k) => {
      const option = Object.assign(el("a", "", d.name), { href: `?study=${encodeURIComponent(d.id)}`, id: `${list.id}-${k}` });
      option.setAttribute("role", "option");
      option.append(" ", el("span", "hint", d.anc));
      option.addEventListener("mousedown", (event) => event.preventDefault());  // the field keeps the focus, so the click lands before the list closes
      option.addEventListener("click", link(() => pick(d)));
      return option;
    }) : [el("span", "hint", "No signed study matches.")]));
    list.hidden = false;
    input.setAttribute("aria-expanded", "true");
    active = -1;
  };
  input.addEventListener("input", () => (input.value.trim() ? open() : pick(null)));
  input.addEventListener("focus", () => { input.select(); open(); });
  input.addEventListener("blur", () => { close(); input.value = picked()?.name ?? ""; });  // what was typed and not picked goes
  input.addEventListener("keydown", (event) => {
    if (event.key === "Escape") return close();
    if (event.key === "Enter") {
      event.preventDefault();
      if (!list.hidden && found.length) pick(found[Math.max(0, active)]);
    } else if ((event.key === "ArrowDown" || event.key === "ArrowUp") && found.length) {
      event.preventDefault();
      if (list.hidden) open();
      mark((active + (event.key === "ArrowDown" ? 1 : found.length - 1) + (active < 0 && event.key === "ArrowUp" ? 1 : 0)) % found.length);
    }
  });
}

function field(label, ...value) {
  const row = el("div", "row"), name = el("span", "label"), content = el("div", "value");
  name.append(...[label].flat());
  content.append(...value);
  row.append(name, content);
  return row;
}

function pairResult(a, b, e) {
  const head = el("h3"), fields = el("div", "fields");
  head.append(studyLink(a), " and ", studyLink(b));
  if (e.rg == null) fields.append(field(rgLabel(), `not estimated: LDSC gave none for this pair${e.note ? ` (${e.note})` : ""}`));
  else fields.append(field(rgPlusMinusLabel(), plusMinus(e.rg, e.se)), field("z", sig(e.z)), field("p", pValue(e.p)));
  if (e.gcov != null) fields.append(field("Genetic covariance ± SE", plusMinus(e.gcov, e.gcov_se)));
  $("pair-result").replaceChildren(head, fields);
}

async function renderRg() {
  const { a, b } = state.rg;
  for (const [input, d] of [[$("rg-a"), a], [$("rg-b"), b]]) if (document.activeElement !== input) input.value = d?.name ?? "";
  $("rg-b").disabled = !a;
  $("pair-result").hidden = !(a && b);
  $("neighbours").hidden = !a;
  drawHeat();
  if (!a) return;
  const pending = setTimeout(() => { $("neighbours-head").textContent = "Loading the genetic correlations…"; }, 150);  // said only if it takes a moment
  let t;
  try {
    t = await correlations(a.anc);
  } catch {
    $("neighbours-head").textContent = "The genetic correlations did not load; check the connection and reload.";
    return;
  } finally {
    clearTimeout(pending);
  }
  if (state.view !== "rg" || state.rg.a !== a || state.rg.b !== b) return;  // another pair was asked for meanwhile
  if (b) pairResult(a, b, estimate(t, t.at.get(a.id), t.at.get(b.id)));
  const pairs = neighbours(t, a), missing = pairs.filter((r) => r.rg == null).length, listed = state.rg.strict ? strong(pairs) : pairs;
  const sort = state.sort.neighbours, columns = COLUMNS.neighbours;
  const ordered = sort ? sorted(listed, sort, (r) => (sort.key === "name" ? r.d.name : r[sort.key])) : [...listed].sort(strongest);
  const head = el("span");
  head.append(studyLink(a), `: ${state.rg.strict ? `${int(listed.length)} of ` : ""}${howMany(pairs.length - missing, "genetic correlation")}`
    + `${missing ? `, ${int(missing)} not estimated` : ""}${sort ? "" : ", the strongest first"}`);
  $("neighbours-head").replaceChildren(head);
  $("rg-strict").checked = state.rg.strict;
  const bodies = (state.rg.all ? ordered : ordered.slice(0, TOP)).map((r) => body(columns, r));
  $("neighbours-table").replaceChildren(header(columns, sort, (next) => { state.sort.neighbours = next; renderRg(); }),
    ...(bodies.length ? bodies : [nothing(columns, `No genetic correlation with |z| ≥ ${STRONG_Z}.`)]));
  const more = $("neighbours-more");
  more.hidden = ordered.length <= TOP;
  more.textContent = state.rg.all ? `show the first ${TOP}` : `show all ${int(ordered.length)}`;
}

// ---------------------------------------------------------------- a study's page

// A paper in full: authors, title, journal, and its links
function citation(pmid, paper) {
  const line = el("div", "citation"), issue = `${paper.volume}${paper.issue ? `(${paper.issue})` : ""}`;
  line.append(`${paper.authors.map((name) => name.replace(/^\{|\}$/g, "")).join(", ")}${paper.more ? ", et al" : ""}. ${paper.title}. `
    + `${paper.journal} ${paper.year}${issue ? `;${issue}` : ""}${paper.pages ? `:${paper.pages}` : ""}. Published ${day(paper.published)}. `,
    outside(`${PUBMED}${pmid}/`, "PubMed"), ...(paper.doi ? [" · ", outside(`${DOI}${paper.doi}`, "DOI")] : []));
  return line;
}

// What the quality checks found of a study, the worst first; what was measured of its file; where to report a problem
function checked(d) {
  const block = el("div"), head = el("div"), measured = el("div", "hint");
  const report = outside(`${ISSUES}?${new URLSearchParams({ template: "study.yml", title: `${d.name}: `, study: `${d.name} (${d.id})` })}`, "report a problem");
  head.append(verdict(d), " · ", report);
  block.append(head, ...(d.notes ?? []).map(([level, text]) => {
    const line = el("div", "finding");
    line.append(el("span", `qc ${level}`, level), " ", ...finding(text, true));
    return line;
  }));
  if (d.hm3 != null) {
    measured.textContent = `Its file: ${Math.floor(d.hm3 * 100)}% of the HapMap3 SNPs that LDSC regresses on; mean χ² ${d.chi2.toFixed(2)}; λGC ${d.lambda.toFixed(2)}.`;
    block.append(measured);
  }
  return block;
}

async function renderStudy() {
  const d = state.study, page = $("view-study"), fields = el("div", "fields"), cited = el("div"), ticking = el("label", "check");
  cited.append(...(d.cite ? d.cite.map((pmid) => { const line = el("div", "citation"); line.append(outside(`${PUBMED}${pmid}/`, state.papers[pmid].label)); return line; }) : [d.ref]));
  if (d.cite) papers().then((all) => { if (state.study === d) cited.replaceChildren(...d.cite.map((pmid) => citation(pmid, all[pmid]))); }, () => {});  // the labels stay if it fails
  ticking.append(tickBox(d), " ticked for export");
  fields.append(
    field("ID", ...parts(d.id)), field("Release name", d.trait), field("Ancestry", d.anc), field("N", int(d.n)),
    field("M", d.m_rg != null && d.m_rg !== d.m ? `${int(d.m)} SNPs; ${int(d.m_rg)} without the strand-ambiguous ones, for rg` : `${int(d.m)} SNPs`),
    field("Release", d.origin), field("Reference", cited),
    ...(httpUrl(d.catalog) ? [field("GWAS Catalog", outside(d.catalog, d.id))] : []),
    field("Signed", yesNo(d.signed)), field("Allele info", yesNo(d.alleles)), field("Cross-trait LDSC", yesNo(d.cross)),
    field("Nonredundant panel", d.panel ? `yes: ${int(d.panel)} of ${int(state.release[d.anc].panel)} in the ${d.anc} panel, by h² z` : "no"),
    field("h² ± SE", plusMinus(d.h2, d.h2_se)), field("h² z", sig(d.h2_z)), field("Intercept", fixed(4)(d.intercept)),
    field("QC", checked(d)),
    field("Path", pathLine(d) ?? action("login to see it", () => openLogin(true))));
  const head = el("h2"), section = el("div");
  head.append(...parts(d.name));
  page.replaceChildren(head, ticking, fields, section);
  if (!d.cross) return section.append(el("p", "section-head", "No genetic correlations: the sumstats are unsigned."));
  section.append(el("p", "section-head", "Loading the genetic correlations…"));
  let t;
  try {
    t = await correlations(d.anc);
  } catch {
    return section.replaceChildren(el("p", "section-head", "The genetic correlations did not load; check the connection and reload."));
  }
  if (state.study !== d || !section.isConnected) return;  // another page was opened meanwhile
  const pairs = neighbours(t, d), top = strong(pairs).sort(strongest).slice(0, TOP);
  const title = el("p", "section-head"), all = Object.assign(el("a", "", `all ${int(pairs.length)} pairs`), { href: `?tab=rg&a=${encodeURIComponent(d.id)}` });
  title.append(top.length ? `Strongest genetic correlations with |z| ≥ ${STRONG_Z}` : `No genetic correlation with |z| ≥ ${STRONG_Z}`, " · ", all);
  const table = el("table", "data pairs"), scroll = el("div", "scroll");
  table.append(header(COLUMNS.neighbours), ...top.map((r) => body(COLUMNS.neighbours, r)));
  scroll.append(table);
  section.replaceChildren(title, ...(top.length ? [scroll] : []));
}

// ---------------------------------------------------------------- candidates for the next release

const day = (date) => date.replaceAll("-", "/");
const paperOf = (p) => outside(`${PUBMED}${p.pmid}/`, `${p.author}, ${p.journal}, ${day(p.published)}`);
const cells = (tag, row) => {  // a table row: text, or what to put in the cell; a number to the right
  const tr = el("tr");
  tr.append(...row.map((content) => {
    const cell = el(tag, typeof content === "number" ? "num" : "");
    cell.append(...[typeof content === "number" ? int(content) : content].flat());
    return cell;
  }));
  return tr;
};
const grid = (table, head, rows) => {
  const thead = el("thead"), tbody = el("tbody");
  thead.append(cells("th", head));
  tbody.append(...rows.map((row) => cells("td", row)));
  table.replaceChildren(thead, tbody);
};

// The candidates the GWAS Catalog holds for the next release, the rule that finds them, and what the rule finds of
// the studies that the last release added (candidates/, data/candidates.json)
async function renderCandidates() {
  let c;
  try {
    c = await (found ??= load("candidates.json").catch((e) => { found = undefined; throw e; }));
  } catch {
    $("cand-lead").textContent = "The candidates did not load; check the connection and reload.";
    return;
  }
  if (state.view !== "candidates") return;
  const year = state.release.year, test = c.backtest, last = test.steps.at(-1), beyond = test.beyond;
  const studies = c.candidates.flatMap((p) => p.studies.map((s, k) => [k ? "" : paperOf(p), outside(`${CATALOG}${s.id}`, s.id), s.trait, s.n, s.ancestry]));
  $("cand-lead").textContent = `${howMany(studies.length, "study", "studies")} in ${howMany(c.candidates.length, "paper")}: GWAS with public summary statistics that the GWAS Catalog added after `
    + `${day(c.until)}, where the ${year} release stopped, and that pass the rule below. The Catalog was read on ${day(c.catalog.fetched)}; its latest paper was then of ${day(c.catalog.published_to)}.`;
  $("cand-list").hidden = !studies.length;
  grid($("cand-list"), ["Paper", "Study", "Trait", "N", "Ancestry"], studies);
  const aside = c.aside.flatMap((group) => group.papers.map((p) => [paperOf(p), ` (${howMany(p.studies, "study", "studies")}: ${group.why})`]));
  $("cand-aside").replaceChildren(...(aside.length ? ["Set aside: ", ...aside.flatMap((item, k) => [k ? "; " : "", ...item]), "."] : []),
    c.other ? ` ${howMany(c.other, "more study", "more studies")} of the Catalog ${c.other === 1 ? "is" : "are"} of another ancestry or too small.` : "");
  $("cand-rule").replaceChildren(...["Its summary statistics are public in the Catalog.", `Its paper has at most ${int(c.rule.max_studies)} such studies.`,
    `Its individuals are of one ancestry: ${c.rule.ancestries.join(" or ")}.`, `It has at least ${int(c.rule.min_n)} individuals.`,
    `Its paper is in a journal that StatGenPapers reads: ${c.rule.journals.join(", ")}.`].map((text) => el("li", "", text)));
  $("cand-test").textContent = `The backtest: the same rule, on the GWAS published from ${day(test.since)} to ${day(test.until)}, against the ${int(test.taken)} studies of `
    + `${int(test.taken_papers)} papers that the ${year} release added. It finds ${int(last.found)} of the ${int(test.taken)} studies (${(100 * last.found / test.taken).toFixed(1)}%) and `
    + `${int(last.found_papers)} of the ${int(test.taken_papers)} papers.`;
  grid($("cand-steps"), ["GWAS published in those days", "Studies", "Papers", `Of the ${int(test.taken)} added`, `Of their ${int(test.taken_papers)} papers`],
    test.steps.map((step, k) => [`${k ? "… and " : ""}${step.rule}`, step.studies, step.papers, step.found, step.found_papers]));
  const why = new Map();  // reason → the studies missed for it
  for (const m of test.missed) why.set(m.why, [...(why.get(m.why) ?? []), state.byId.has(m.id) ? studyLink(state.byId.get(m.id)) : m.id]);
  $("cand-missed").replaceChildren(...(why.size ? ["Missed: ", ...[...why].flatMap(([reason, missed], k) => [k ? "; " : "", ...missed.flatMap((study, j) => [j ? ", " : "", study]), ` (${reason})`]), "."]
    : ["It misses none."]));
  const more = $("cand-more"), toggle = action("", () => { more.hidden = !more.hidden; label(); }), label = () => { toggle.textContent = more.hidden ? "show the papers" : "hide the papers"; };
  label();
  $("cand-beyond").replaceChildren(`Beyond those, the rule finds ${howMany(beyond.studies, "study", "studies")} that the release did not add: ${int(beyond.in_taken_papers)} in papers it took others from, `
    + `${int(beyond.studies - beyond.in_taken_papers)} in ${howMany(beyond.new_papers.length, "other paper")}. `, toggle);
  more.replaceChildren(...beyond.new_papers.map((p) => {
    const item = el("li");
    item.append(paperOf(p), `: ${p.title} (${howMany(p.studies.length, "study", "studies")})`);
    return item;
  }));
}

// ---------------------------------------------------------------- the bar, the tabs, the address

function renderBar() {
  $("login-open").textContent = state.paths ? "logout" : "login";
  if (state.paths) $("login-open").removeAttribute("aria-expanded");  // "logout" opens nothing
  else $("login-open").setAttribute("aria-expanded", !$("login").hidden);
}

// What only a logged-in reader sees in the texts: the release folder, and the folders of the files
function renderPrivate() {
  const p = state.paths, release = $("o2-release"), folders = $("o2-folders");
  release.hidden = folders.hidden = !p;
  if (!p) return release.replaceChildren(), folders.replaceChildren();
  const code = (text) => el("code", "", text), files = Object.values(p.files);
  const note = el("p", "o2");
  note.append("On O2, the release is in ", code(p.release), ", with ", ...p.tables.flatMap((name, k) => [k ? ", " : "", code(name)]), ".");
  release.replaceChildren(el("h2", "", "Where is it on O2?"), note);
  const list = el("ol", "o2");
  list.append(...p.folders.map((folder, k) => {
    const item = el("li");
    item.append(code(folder), `: ${howMany(files.filter((file) => file[0] === k).length, "file")}`);
    return item;
  }));
  folders.replaceChildren(el("p", "", "On O2, a file is named after its study's ID, in one of these folders:"), list);
}

function render() {
  renderBar();
  if (!state.loaded) return;
  if (state.view === "sumstats" || state.view === "h2") renderList();
  else if (state.view === "rg") renderRg();
  else if (state.view === "candidates") renderCandidates();
  else if (state.view === "study") renderStudy();
  renderPrivate();
}

// The address says what to show: ?tab=h2 | rg | notes | about, ?study=ID, or nothing for the first tab; under
// ?tab=rg, a and b name the pair looked up.
function route() {
  const q = new URLSearchParams(location.search), signed = (id) => (state.byId.get(id)?.cross ? state.byId.get(id) : null);
  if (q.has("dataset")) {  // a link from before a row was called a study
    q.set("study", q.get("dataset"));
    q.delete("dataset");
    history.replaceState(null, "", `?${q}`);
  }
  state.study = state.byId.get(q.get("study")) ?? null;
  state.view = state.study || (q.has("study") && !state.loaded) ? "study" : Object.hasOwn(TABS, q.get("tab")) ? q.get("tab") : "sumstats";
  const a = signed(q.get("a")), b = signed(q.get("b"));
  if (a !== state.rg.a) state.rg.all = false;
  if (a && a.anc !== state.rg.heat?.anc) state.rg.heat = null;  // the heatmap follows the first study's ancestry
  Object.assign(state.rg, { a, b: a && b && b !== a && b.anc === a.anc ? b : null });
  for (const name of [...Object.keys(TABS), "study"]) $(`view-${name}`).hidden = name !== state.view;
  $("listview").hidden = state.view !== "sumstats" && state.view !== "h2";
  for (const tab of document.querySelectorAll("[data-tab]")) {
    const current = tab.dataset.tab === state.view;
    tab.classList.toggle("current", current);
    current ? tab.setAttribute("aria-current", "page") : tab.removeAttribute("aria-current");
  }
  document.title = [state.study?.name ?? TABS[state.view], TITLE].filter(Boolean).join(" · ");
  hideTip();
  render();
}

// Every navigation adds a history entry, so Back returns to the previous tab or study.
function go(params) {
  const search = String(new URLSearchParams(params));
  if (search !== location.search.slice(1)) history.pushState(null, "", search ? `?${search}` : location.pathname);
  route();
}

// ---------------------------------------------------------------- start

function start(data) {
  state.release = data.release;
  state.studies = data.studies;
  state.papers = data.papers;
  for (const d of data.studies) {
    d.label = d.cite ? state.papers[d.cite[0]].label : d.ref;  // its paper, as the table lists it
    d.published = d.cite ? state.papers[d.cite[0]].published : "";  // and the day that paper was published, by which the table sorts it
    d.rank = QC.indexOf(d.qc);  // the worst sorts first
    d.text = [d.name, d.id, d.trait, d.origin, d.ref, ...(d.cite ?? []).map((pmid) => state.papers[pmid].label)].join("\n").toLowerCase();  // what the search reads
    state.byId.set(d.id, d);
  }
  state.loaded = true;
  readCart();
  $("f-anc").append(...data.release.ancestries.map((anc) => el("option", "", anc)));
  for (const node of document.querySelectorAll("[data-v]")) {  // the numbers in the texts, from the release itself; a study by its name
    const value = node.dataset.v.split(".").reduce((o, k) => o?.[k], data.release);
    if (state.byId.has(value)) node.replaceChildren(studyLink(state.byId.get(value)));
    else if (value != null) node.textContent = typeof value === "number" ? int(value) : String(value).replace(/^(\d{4})-(\d\d)-(\d\d)$/, "$1/$2/$3");
  }
  for (const node of document.querySelectorAll("[data-study]")) {  // a study the texts name
    if (state.byId.has(node.dataset.study)) node.replaceChildren(studyLink(state.byId.get(node.dataset.study)));
  }
  route();
  restore();
}

function listen() {
  const f = state.filter;
  $("filters").addEventListener("submit", (event) => event.preventDefault());
  $("q").addEventListener("input", () => { f.q = $("q").value; renderList(); });
  for (const name of ["anc", "qc"]) $(`f-${name}`).addEventListener("change", () => { f[name] = $(`f-${name}`).value; renderList(); });
  for (const name of ["z6", "signed", "panel", "ticked"]) $(`f-${name}`).addEventListener("change", () => { f[name] = $(`f-${name}`).checked; renderList(); });

  $("pair").addEventListener("submit", (event) => event.preventDefault());
  const signed = () => state.studies.filter((d) => d.cross);
  picker($("rg-a"), $("rg-a-list"), signed, () => state.rg.a, (d) => go(d ? { tab: "rg", a: d.id, ...(state.rg.b && state.rg.b.anc === d.anc && state.rg.b !== d ? { b: state.rg.b.id } : {}) } : { tab: "rg" }));
  picker($("rg-b"), $("rg-b-list"), () => signed().filter((d) => d.anc === state.rg.a?.anc && d !== state.rg.a), () => state.rg.b, (d) => go({ tab: "rg", a: state.rg.a.id, ...(d ? { b: d.id } : {}) }));
  $("rg-strict").addEventListener("change", () => { Object.assign(state.rg, { strict: $("rg-strict").checked, all: false }); renderRg(); });
  $("neighbours-more").addEventListener("click", link(() => { state.rg.all = !state.rg.all; renderRg(); }));

  $("login-open").addEventListener("click", link(() => (state.paths ? logout() : openLogin())));
  $("login-cancel").addEventListener("click", () => { closeLogin(); $("login-open").focus({ preventScroll: true }); });
  $("login").addEventListener("submit", login);
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && !$("login").hidden) { closeLogin(); $("login-open").focus({ preventScroll: true }); }
  });
  addEventListener("storage", (event) => {  // a login, a logout or a tick in another tab of this browser applies here too
    if (event.key === "key") restore();
    if (event.key === "cart" && state.loaded) { readCart(); for (const cache of Object.values(rows)) cache.clear(); render(); }
  });

  // A link to another tab or study of this page changes the page in place; Ctrl, Cmd or Shift opens it as the browser would.
  document.addEventListener("click", (event) => {
    const a = event.target.closest?.("a[href]");
    if (!a || event.defaultPrevented || event.button || event.ctrlKey || event.metaKey || event.shiftKey || a.target || a.download || a.getAttribute("href").startsWith("#")) return;
    const url = new URL(a.href);
    if (url.origin !== location.origin || url.pathname !== location.pathname) return;
    event.preventDefault();
    go(Object.fromEntries(url.searchParams));
    scrollTo(0, 0);
  });
  addEventListener("popstate", route);
  let width = innerWidth;
  addEventListener("resize", () => {  // the charts are drawn to the column's width
    if (innerWidth === width || !state.loaded) return;
    width = innerWidth;
    if (state.view === "h2") drawHistogram(shown());
    if (state.view === "rg") drawHeat();
  });
  matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {  // the heatmap's picture has its colours painted in
    delete $("heat").dataset.drawn;
    if (state.loaded && state.view === "rg") drawHeat();
  });
}

listen();
route();  // the tab the address names shows at once; its data follows
load("studies.json").then(start, () => { $("failed").hidden = false; });
