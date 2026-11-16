// GWAS Sumstats: a catalog read from data/*.json (written by build/). Where the files are on O2 is in
// data/paths.enc.json, encrypted under the lab's passphrase; it is decrypted here, and nothing leaves this browser.

const TITLE = "GWAS Sumstats";
const TABS = { sumstats: "", h2: "Heritability", rg: "Genetic correlation", notes: "Release notes", about: "About" };  // tab → its page title
const STRONG_Z = 3;  // a genetic correlation is listed among a dataset's strongest if its |z| is at least this
const TOP = 10;      // the strongest shown before "show all"
const PANEL_MAX = Math.sqrt(0.1);  // a panel's heatmap: its colours end at the panels' bound, rg² = 0.1
const $ = (id) => document.getElementById(id);
const el = (tag, cls, text) => Object.assign(document.createElement(tag), { className: cls || "", textContent: text ?? "" });
const svg = (tag, attrs = {}) => {
  const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [name, value] of Object.entries(attrs)) node.setAttribute(name, value);
  return node;
};

const state = {
  loaded: false, datasets: [], byId: new Map(), release: null,
  view: "sumstats",   // a tab, or "dataset": the page of state.dataset
  dataset: null,
  filter: { q: "", anc: "", z6: false, signed: false, panel: false },  // of the datasets under "sumstats" and "heritability"
  sort: { sumstats: null, h2: { key: "h2", dir: -1 }, neighbours: null },  // null: the release's order; for neighbours, the strongest first
  rg: { a: null, b: null, strict: true, all: false, heat: null },  // the pair looked up, the neighbours listed, the heatmap drawn ({ anc, kind })
  paths: null,        // where the files are, once logged in
};
const tables = new Map();  // ancestry → promise of its genetic correlations
const rows = { sumstats: new Map(), h2: new Map() };  // dataset id → its <tbody>, per table
let blob;  // promise of the encrypted paths

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
const pValue = (p) => (p == null ? "–" : p === 0 ? "<1e-300" : p < 1e-4 ? p.toExponential(2) : p.toFixed(4));  // LDSC prints 0 below about 1e-308
const yesNo = (flag) => (flag ? "yes" : "no");
const plural = (n, noun) => `${int(n)} ${noun}${n === 1 ? "" : "s"}`;
const words = (text) => text.toLowerCase().split(/\s+/).filter(Boolean);
const httpUrl = (url) => (/^https:\/\/\S+$/.test(url ?? "") ? url : null);
const rgLabel = () => {
  const label = el("span", "", "r");
  label.append(el("sub", "", "g"));
  return label;
};
// A long name in parts (PASS_Height_Yengo2022, NG_Kanai_2018_BBJ): it may wrap after an underscore
const parts = (text) => text.split("_").flatMap((part, k, all) => (k < all.length - 1 ? [`${part}_`, el("wbr")] : [part]));
const outside = (url, ...content) => {
  const a = Object.assign(el("a"), { href: url, target: "_blank", rel: "noopener noreferrer" });
  a.append(...content);
  return a;
};
const datasetLink = (d) => Object.assign(el("a", "", d.trait), { href: `?dataset=${encodeURIComponent(d.id)}` });
const reference = (d) => (httpUrl(d.paper) ? outside(d.paper, ...parts(d.ref)) : parts(d.ref));

const load = (name) => fetch(`data/${name}`, { cache: "no-cache" }).then((r) => (r.ok ? r.json() : Promise.reject(new Error(`${name}: ${r.status}`))));

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

const position = (t, i, j) => (i > j ? position(t, j, i) : i * (2 * t.n - i - 1) / 2 + j - i - 1);  // of the pair of datasets i and j in the lists

function estimate(t, i, j) {
  const k = position(t, i, j);
  return { rg: t.rg[k], se: t.se[k], z: t.z[k], p: t.p[k], gcov: t.gcov[k], gcov_se: t.gcov_se[k], note: t.notes[k] };
}

// Every other dataset of d's ancestry with its estimate against d
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
  for (const cache of Object.values(rows)) cache.clear();  // the rows gain or lose their files
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

// A dataset's h² file (kind 0) or rg file (kind 1) on O2; null when logged out or when it has none
function filePath(d, kind) {
  const file = state.paths?.files[d.id]?.[kind];
  return file ? `${state.paths.folders[file[0]]}/${file[1]}` : null;
}

function copy(text) {
  navigator.clipboard.writeText(text).then(() => toast("Copied."), () => toast("Could not copy; select the path instead."));
}

// A dataset's files, each with "copy"; null when logged out
function fileList(d) {
  const h2 = filePath(d, 0), rg = filePath(d, 1);
  if (!h2) return null;
  const list = el("div", "paths");
  for (const [label, path] of rg === h2 ? [["h² and rg file", h2]] : [["h² file", h2], ["rg file", rg]]) {
    if (!path) continue;
    const row = el("div");
    row.append(`${label}: `, el("span", "path", path), action("copy", () => copy(path)));
    list.append(row);
  }
  return list;
}

// ---------------------------------------------------------------- the datasets: filters, tables

function shown() {
  const f = state.filter, sought = words(f.q);
  return state.datasets.filter((d) => (!f.anc || d.anc === f.anc) && (!f.z6 || d.h2_z > 6) && (!f.signed || d.signed) && (!f.panel || d.panel)
    && sought.every((w) => d.text.includes(w)));
}

const compare = (x, y) => (typeof x === "string" ? x.localeCompare(y, "en", { numeric: true, sensitivity: "base" }) : (x ?? -Infinity) - (y ?? -Infinity));
// `list` by a column, in the direction asked; ties, and no sort at all, keep the order of `list`
const sorted = (list, sort, value = (item) => item[sort.key]) => (sort ? [...list].sort((a, b) => sort.dir * compare(value(a), value(b))) : list);

function meter(h2) {  // h² as a share of 1; an estimate outside 0 to 1 shows as the nearer end
  const track = el("span", "meter"), fill = el("i");
  fill.style.width = `${Math.max(0, Math.min(1, h2)) * 100}%`;
  track.append(fill);
  return track;
}

// A column: key (the field it shows and sorts by), label, title (what the label means), num (a format: right-aligned,
// largest first), show (the cell's content, if not the field itself), cls.
const TRAIT = { key: "trait", label: "Trait", cls: "trait", show: datasetLink };
const ID = { key: "id", label: "ID", cls: "id", show: (d) => parts(d.id) };
const ANCESTRY = { key: "anc", label: "Ancestry", cls: "short" };
const COLUMNS = {
  sumstats: [TRAIT, ID, ANCESTRY,
    { key: "n", label: "N", title: "Sample size", num: int },
    { key: "m", label: "M", title: "SNPs in the h² file", num: int },
    { key: "h2", label: "h²", title: "Heritability: S-LDSC, observed scale", num: fixed(4) },
    { key: "h2_z", label: "h² z", title: "h² over its standard error", num: fixed(2) },
    { key: "signed", label: "Signed", title: "Signed sumstats: used for genetic correlation", cls: "short", show: (d) => yesNo(d.signed) },
    { key: "origin", label: "Release", title: "The update that added the dataset", cls: "short" },
    { key: "ref", label: "Reference", cls: "ref", show: reference }],
  h2: [TRAIT, ID, ANCESTRY,
    { key: "h2", label: "h²", title: "Heritability: S-LDSC, observed scale", num: fixed(4), show: (d) => [meter(d.h2), fixed(4)(d.h2)] },
    { key: "h2_se", label: "SE", title: "Standard error of h²", num: fixed(4) },
    { key: "h2_z", label: "z", title: "h² over its standard error", num: fixed(2) },
    { key: "intercept", label: "Intercept", title: "Intercept of the LD score regression", num: fixed(4) }],
  neighbours: [{ ...TRAIT, show: (r) => datasetLink(r.d) }, { ...ID, show: (r) => parts(r.d.id) },
    { key: "rg", label: rgLabel, title: "Genetic correlation", num: fixed(4) },
    { key: "se", label: "SE", title: "Standard error of the genetic correlation", num: fixed(4) },
    { key: "z", label: "z", num: fixed(2) },
    { key: "p", label: "p", num: pValue }],
};

const labelOf = (c) => (typeof c.label === "function" ? c.label() : c.label);

// The header row. With `onSort`, a click on a column sorts by it, a second turns the order round, a third returns to
// the first order.
function header(columns, sort, onSort) {
  const tr = el("tr");
  for (const c of columns) {
    const th = el("th", c.cls || (c.num ? "num" : "")), first = c.num ? -1 : 1;
    th.scope = "col";
    if (onSort) {
      const button = el("button");
      button.type = "button";
      button.append(labelOf(c));
      button.title = `${c.title ? `${c.title}. ` : ""}Sort by this column.`;
      if (sort?.key === c.key) {
        th.setAttribute("aria-sort", sort.dir > 0 ? "ascending" : "descending");
        button.append(sort.dir > 0 ? " ▲" : " ▼");
      }
      button.addEventListener("click", () => onSort(sort?.key !== c.key ? { key: c.key, dir: first } : sort.dir === first ? { key: c.key, dir: -first } : null));
      th.append(button);
    } else {
      th.append(labelOf(c));
      th.title = c.title ?? "";
    }
    tr.append(th);
  }
  const thead = el("thead");
  thead.append(tr);
  return thead;
}

// One <tbody> per row: its cells and, under them, `below` (a dataset's files)
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
  const view = state.view, columns = COLUMNS[view], list = shown(), all = state.datasets.length;
  $("count").textContent = list.length === all ? plural(all, "dataset") : `${int(list.length)} of ${plural(all, "dataset")}`;
  const cache = rows[view];
  const bodies = sorted(list, state.sort[view]).map((d) => cache.get(d.id) ?? cache.set(d.id, body(columns, d, fileList(d))).get(d.id));
  const head = header(columns, state.sort[view], (sort) => { state.sort[view] = sort; renderList(); });
  $(`${view}-table`).replaceChildren(head, ...(bodies.length ? bodies : [nothing(columns, "No dataset matches.")]));
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

// How many datasets have each h²: a column per bin of 0.02, in a frame that stays as the filters change
function drawHistogram(list) {
  const figure = $("h2-plot"), width = figure.clientWidth, height = 190, left = 30, right = 8, top = 8, bottom = 22;
  if (!width) return;
  const bin = (d) => Math.floor(Math.round(d.h2 * 1e4) / 200);  // h² has four decimals: in whole numbers it bins exactly
  const all = state.datasets.map(bin), first = Math.min(...all), bins = Math.max(...all) - first + 1, w = (width - left - right) / bins;
  const counts = new Array(bins).fill(0);
  for (const d of list) counts[bin(d) - first]++;
  const peak = Math.max(1, ...counts), tick = [1, 2, 5, 10, 20, 50, 100, 200].find((t) => peak / t <= 5) ?? 500;
  const y = (count) => top + (height - top - bottom) * (1 - count / (Math.ceil(peak / tick) * tick));
  const range = (k) => `${((first + k) * 0.02).toFixed(2)} to ${((first + k + 1) * 0.02).toFixed(2)}`;
  const chart = svg("svg", { width, height, viewBox: `0 0 ${width} ${height}`, role: "img" });
  for (let count = 0; count <= peak; count += tick) {  // a hairline and its count, every `tick` datasets
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
    showTip(event, el("b", "", plural(counts[k], "dataset")), ` with h² from ${range(k)}`);
  });
  chart.addEventListener("pointerleave", () => { hideTip(); bars.forEach((bar) => bar.classList.remove("on")); });
  const label = list.length ? `h² of the ${plural(list.length, "dataset")} shown, in bins of 0.02: most often ${range(counts.indexOf(peak))}.` : "";
  chart.setAttribute("aria-label", label || "No dataset matches.");
  figure.replaceChildren(chart);
  $("h2-caption").textContent = label;
}

// The colours of a genetic correlation, as [red, green, blue]: `of(v)` for v from -1 (--neg) through 0 (--subtle) to 1
// (--mark), and `self` for a dataset against itself. Read from the style sheet, so they follow the light or dark mode.
function shades() {
  const css = getComputedStyle(document.documentElement);
  const [neg, zero, pos, self] = ["--neg", "--subtle", "--mark", "--border"].map((name) => {
    const hex = css.getPropertyValue(name).trim();
    return [1, 3, 5].map((k) => parseInt(hex.slice(k, k + 2), 16));
  });
  return { self, of: (v) => zero.map((c, k) => Math.round(c + ((v < 0 ? neg : pos)[k] - c) * Math.min(1, Math.abs(v)))) };
}

const heatShown = () => state.rg.heat ?? { anc: state.rg.a?.anc ?? state.release.ancestries[0], kind: "panel" };

// A heatmap of the genetic correlations among the datasets of a nonredundant panel ("panel") or among all signed
// datasets of an ancestry ("all"): blue above 0, red below. The build ordered the datasets by hierarchical clustering
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
  const { order, merges } = t.heat[kind], n = order.length, datasets = order.map((i) => state.byId.get(t.ids[i]));
  const limit = kind === "panel" ? PANEL_MAX : 1, colours = shades();
  const TREE = 90, names = n <= 100 ? 230 : 0;  // the tree on the left; the traits' names on the right, where a row is tall enough to read one
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
    "aria-label": `Genetic correlations among ${kind === "panel" ? `the ${n} datasets of the ${anc} nonredundant panel` : `the ${n} signed ${anc} datasets`}, ordered by hierarchical clustering` });
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
    for (const [k, d] of datasets.entries()) {
      const name = svg("text", { x: TREE + size + 8, y: (k + 0.5) * cell, "dominant-baseline": "central" });
      name.textContent = d.trait.length > fits ? `${d.trait.slice(0, fits - 1).trimEnd()}…` : d.trait;
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
    showTip(event, ...(e ? [el("b", "", e.rg == null ? "not estimated" : fixed(4)(e.rg)), e.rg == null ? "" : ` (SE ${fixed(4)(e.se)})`, el("br"), `${datasets[i].trait} × ${datasets[j].trait}`] : [datasets[i].trait]));
  });
  chart.addEventListener("pointerleave", leave);
  chart.addEventListener("click", (event) => {  // a cell's pair, looked up above
    const here = under(event);
    if (!here || here[0] === here[1]) return;
    leave();
    state.rg.heat = shown;  // the heatmap stays, whatever the pair
    go({ tab: "rg", a: datasets[here[0]].id, b: datasets[here[1]].id });
    $("pair").scrollIntoView({ block: "start" });
  });
  figure.replaceChildren(chart);

  const scale = el("span", "scale"), steps = el("span", "steps");
  steps.append(...[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1].map((v) => { const swatch = el("i"); swatch.style.background = `rgb(${colours.of(v).join(" ")})`; return swatch; }));
  scale.append(rgLabel(), ` −${limit.toFixed(2)}`, steps, `+${limit.toFixed(2)}`);
  $("heat-caption").replaceChildren(scale, ` Red below 0, blue above${kind === "panel" ? "; no pair of a panel is beyond these" : ""}. A pair without an estimate shows as 0.`
    + " Rows and columns are in the order of the tree: hierarchical clustering on 1 − |", rgLabel(), "|, average linkage. A click on a cell looks its pair up.");
}

// ---------------------------------------------------------------- genetic correlation: a pair, a dataset's neighbours

// A text field that picks a dataset: typing lists the matches of `among()` under it; a click, or the arrows and
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
    if (d) input.value = d.trait;
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
      const option = Object.assign(el("a", "", d.trait), { href: `?dataset=${encodeURIComponent(d.id)}`, id: `${list.id}-${k}` });
      option.setAttribute("role", "option");
      option.append(" ", el("span", "hint", `${d.id} · ${d.anc}`));
      option.addEventListener("mousedown", (event) => event.preventDefault());  // the field keeps the focus, so the click lands before the list closes
      option.addEventListener("click", link(() => pick(d)));
      return option;
    }) : [el("span", "hint", "No signed dataset matches.")]));
    list.hidden = false;
    input.setAttribute("aria-expanded", "true");
    active = -1;
  };
  input.addEventListener("input", () => (input.value.trim() ? open() : pick(null)));
  input.addEventListener("focus", () => { input.select(); open(); });
  input.addEventListener("blur", () => { close(); input.value = picked()?.trait ?? ""; });  // what was typed and not picked goes
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
  name.append(label);
  content.append(...value);
  row.append(name, content);
  return row;
}
const withSe = (value, se) => (value == null ? "–" : `${fixed(4)(value)} (SE ${fixed(4)(se)})`);

function pairResult(a, b, e) {
  const head = el("h3"), fields = el("div", "fields");
  head.append(datasetLink(a), " and ", datasetLink(b));
  if (e.rg == null) fields.append(field(rgLabel(), `not estimated: LDSC gave none for this pair${e.note ? ` (${e.note})` : ""}`));
  else fields.append(field(rgLabel(), withSe(e.rg, e.se)), field("z", fixed(2)(e.z)), field("p", pValue(e.p)));
  if (e.gcov != null) fields.append(field("Genetic covariance", withSe(e.gcov, e.gcov_se)));
  $("pair-result").replaceChildren(head, fields);
}

async function renderRg() {
  const { a, b } = state.rg;
  for (const [input, d] of [[$("rg-a"), a], [$("rg-b"), b]]) if (document.activeElement !== input) input.value = d?.trait ?? "";
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
  const ordered = sort ? sorted(listed, sort, (r) => (sort.key === "trait" ? r.d.trait : sort.key === "id" ? r.d.id : r[sort.key])) : [...listed].sort(strongest);
  const head = el("span");
  head.append(datasetLink(a), ` (${a.id}): ${state.rg.strict ? `${int(listed.length)} of ` : ""}${plural(pairs.length - missing, "genetic correlation")}`
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

// ---------------------------------------------------------------- a dataset's page

async function renderDataset() {
  const d = state.dataset, page = $("view-dataset"), fields = el("div", "fields"), files = fileList(d);
  fields.append(
    field("ID", ...parts(d.id)), field("Ancestry", d.anc), field("N", int(d.n)),
    field("M", d.m_rg != null && d.m_rg !== d.m ? `${int(d.m)} SNPs in the h² file, ${int(d.m_rg)} in the rg file` : `${int(d.m)} SNPs`),
    field("Release", d.origin), field("Reference", ...[reference(d)].flat(), d.pmid ? ` (PMID ${d.pmid})` : ""),
    ...(d.paper_id ? [field("Short name", ...parts(d.paper_id))] : []),
    ...(httpUrl(d.catalog) ? [field("GWAS Catalog", outside(d.catalog, d.id))] : []),
    field("Signed", yesNo(d.signed)), field("Allele info", yesNo(d.alleles)), field("Cross-trait LDSC", yesNo(d.cross)),
    field("Nonredundant panel", d.panel ? `yes: ${int(d.panel)} of ${int(state.release[d.anc].panel)} in the ${d.anc} panel, by h² z` : "no"),
    field("h²", `${withSe(d.h2, d.h2_se)}, z ${fixed(2)(d.h2_z)}`), field("Intercept", fixed(4)(d.intercept)),
    field("Files on O2", files ?? action("login to see them", () => openLogin(true))));
  const head = el("h2", "", d.trait), section = el("div");
  page.replaceChildren(head, fields, section);
  if (!d.cross) return section.append(el("p", "section-head", "No genetic correlations: the sumstats are unsigned."));
  section.append(el("p", "section-head", "Loading the genetic correlations…"));
  let t;
  try {
    t = await correlations(d.anc);
  } catch {
    return section.replaceChildren(el("p", "section-head", "The genetic correlations did not load; check the connection and reload."));
  }
  if (state.dataset !== d || !section.isConnected) return;  // another page was opened meanwhile
  const pairs = neighbours(t, d), top = strong(pairs).sort(strongest).slice(0, TOP);
  const title = el("p", "section-head"), all = Object.assign(el("a", "", `all ${int(pairs.length)} pairs`), { href: `?tab=rg&a=${encodeURIComponent(d.id)}` });
  title.append(top.length ? `Strongest genetic correlations with |z| ≥ ${STRONG_Z}` : `No genetic correlation with |z| ≥ ${STRONG_Z}`, " · ", all);
  const table = el("table", "data pairs"), scroll = el("div", "scroll");
  table.append(header(COLUMNS.neighbours), ...top.map((r) => body(COLUMNS.neighbours, r)));
  scroll.append(table);
  section.replaceChildren(title, ...(top.length ? [scroll] : []));
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
    const item = el("li"), count = (kind) => files.filter((file) => file[kind]?.[0] === k).length;
    const kinds = [[count(0), "h²"], [count(1), "rg"]].filter(([n]) => n);
    item.append(code(folder), `: ${kinds.map(([n, kind]) => `${int(n)} ${kind} files`).join(", ")}`);
    return item;
  }));
  folders.replaceChildren(el("p", "", "On O2, a file is named after its dataset's ID, in one of these folders:"), list);
}

function render() {
  renderBar();
  if (!state.loaded) return;
  if (state.view === "sumstats" || state.view === "h2") renderList();
  else if (state.view === "rg") renderRg();
  else if (state.view === "dataset") renderDataset();
  renderPrivate();
}

// The address says what to show: ?tab=h2 | rg | notes | about, ?dataset=ID, or nothing for the first tab; under
// ?tab=rg, a and b name the pair looked up.
function route() {
  const q = new URLSearchParams(location.search), signed = (id) => (state.byId.get(id)?.cross ? state.byId.get(id) : null);
  state.dataset = state.byId.get(q.get("dataset")) ?? null;
  state.view = state.dataset || (q.has("dataset") && !state.loaded) ? "dataset" : Object.hasOwn(TABS, q.get("tab")) ? q.get("tab") : "sumstats";
  const a = signed(q.get("a")), b = signed(q.get("b"));
  if (a !== state.rg.a) state.rg.all = false;
  if (a && a.anc !== state.rg.heat?.anc) state.rg.heat = null;  // the heatmap follows the first dataset's ancestry
  Object.assign(state.rg, { a, b: a && b && b !== a && b.anc === a.anc ? b : null });
  for (const name of [...Object.keys(TABS), "dataset"]) $(`view-${name}`).hidden = name !== state.view;
  $("listview").hidden = state.view !== "sumstats" && state.view !== "h2";
  for (const tab of document.querySelectorAll("[data-tab]")) {
    const current = tab.dataset.tab === state.view;
    tab.classList.toggle("current", current);
    current ? tab.setAttribute("aria-current", "page") : tab.removeAttribute("aria-current");
  }
  document.title = [state.dataset?.trait ?? TABS[state.view], TITLE].filter(Boolean).join(" · ");
  hideTip();
  render();
}

// Every navigation adds a history entry, so Back returns to the previous tab or dataset.
function go(params) {
  const search = String(new URLSearchParams(params));
  if (search !== location.search.slice(1)) history.pushState(null, "", search ? `?${search}` : location.pathname);
  route();
}

// ---------------------------------------------------------------- start

function start(data) {
  state.release = data.release;
  state.datasets = data.datasets;
  for (const d of data.datasets) {
    d.text = [d.trait, d.id, d.origin, d.ref].join("\n").toLowerCase();  // what the search reads
    state.byId.set(d.id, d);
  }
  state.loaded = true;
  $("f-anc").append(...data.release.ancestries.map((anc) => el("option", "", anc)));
  for (const node of document.querySelectorAll("[data-v]")) {  // the numbers in the texts, from the release itself
    const value = node.dataset.v.split(".").reduce((o, k) => o?.[k], data.release);
    if (value != null) node.textContent = typeof value === "number" ? int(value) : String(value).replace(/^(\d{4})-(\d\d)-(\d\d)$/, "$1/$2/$3");
  }
  route();
  restore();
}

function listen() {
  const f = state.filter;
  $("filters").addEventListener("submit", (event) => event.preventDefault());
  $("q").addEventListener("input", () => { f.q = $("q").value; renderList(); });
  $("f-anc").addEventListener("change", () => { f.anc = $("f-anc").value; renderList(); });
  for (const name of ["z6", "signed", "panel"]) $(`f-${name}`).addEventListener("change", () => { f[name] = $(`f-${name}`).checked; renderList(); });

  $("pair").addEventListener("submit", (event) => event.preventDefault());
  const signed = () => state.datasets.filter((d) => d.cross);
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
  addEventListener("storage", (event) => { if (event.key === "key") restore(); });  // a login or logout in another tab of this browser applies here too

  // A link to another tab or dataset of this page changes the page in place; Ctrl, Cmd or Shift opens it as the browser would.
  document.addEventListener("click", (event) => {
    const a = event.target.closest?.("a[href]");
    if (!a || event.defaultPrevented || event.button || event.ctrlKey || event.metaKey || event.shiftKey || a.target || a.getAttribute("href").startsWith("#")) return;
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
load("datasets.json").then(start, () => { $("failed").hidden = false; });
