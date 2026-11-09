// make_screenshots.mjs - reproducible screenshots of the BaySimChip Shiny app (e.g. for Fig. 2)
//
// Drives a headless Chrome through the DevTools protocol (Node >= 22: built-in
// fetch and WebSocket, no npm packages needed). The app must be running, e.g.
//   Rscript -e "shiny::runApp('BaySimChipR', port = 8765)"
//
// Usage: node tools/make_screenshots.mjs <output-folder>
//   env: BSC_URL (default http://localhost:8765), CHROME (path to chrome.exe / msedge.exe)

import { spawn } from 'node:child_process';
import { mkdirSync, writeFileSync } from 'node:fs';
import path from 'node:path';

const APP = process.env.BSC_URL || 'http://localhost:8765';
const OUT = path.resolve(process.argv[2] || 'screenshots');
const CHROME = process.env.CHROME || 'C:/Program Files/Google/Chrome/Application/chrome.exe';
const PORT = 9333;
const W = 1600, H = 1000, SCALE = 2;   // CSS pixels and device pixel ratio -> 3200 x 2000 PNG

mkdirSync(OUT, { recursive: true });
const chrome = spawn(CHROME, [
  '--headless=new', '--disable-gpu', '--hide-scrollbars', '--no-first-run',
  `--remote-debugging-port=${PORT}`, `--user-data-dir=${path.join(OUT, '.chrome-profile')}`,
  `--window-size=${W},${H}`, 'about:blank'], { stdio: 'ignore' });
process.on('exit', () => chrome.kill());   // also when a step fails

const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
async function devtools(url, method = 'GET') {
  for (let i = 0; i < 100; i++) {
    try { const r = await fetch(url, { method }); return await r.json(); } catch (e) { await sleep(200); }
  }
  throw new Error('Chrome DevTools endpoint not reachable');
}

const target = await devtools(`http://127.0.0.1:${PORT}/json/new?about:blank`, 'PUT');
const ws = new WebSocket(target.webSocketDebuggerUrl);
await new Promise((r) => ws.addEventListener('open', r));
let msgId = 0;
const pending = new Map();
ws.addEventListener('message', (ev) => {
  const m = JSON.parse(ev.data);
  if (m.method === 'Runtime.exceptionThrown') {
    const d = m.params.exceptionDetails;
    console.log('  [page exception] ' + ((d.exception && d.exception.description) || d.text).split('\n')[0]);
  }
  if (m.method === 'Runtime.consoleAPICalled' && ['error', 'warning'].includes(m.params.type)) {
    console.log('  [page ' + m.params.type + '] ' + m.params.args.map(x => x.value ?? x.description).join(' ').slice(0, 200));
  }
  if (m.id && pending.has(m.id)) {
    const { resolve, reject } = pending.get(m.id);
    pending.delete(m.id);
    if (m.error) reject(new Error(JSON.stringify(m.error))); else resolve(m.result);
  }
});
const send = (method, params = {}) => new Promise((resolve, reject) => {
  const i = ++msgId;
  pending.set(i, { resolve, reject });
  ws.send(JSON.stringify({ id: i, method, params }));
});
async function evaluate(expr) {
  const r = await send('Runtime.evaluate', { expression: expr, awaitPromise: true, returnByValue: true });
  if (r.exceptionDetails) {
    throw new Error((r.exceptionDetails.exception && r.exceptionDetails.exception.description) || r.exceptionDetails.text);
  }
  return r.result.value;
}
async function waitFor(expr, timeout = 30000) {
  const t0 = Date.now();
  while (Date.now() - t0 < timeout) {
    try { if (await evaluate(expr)) return; } catch (e) { /* page not ready yet */ }
    await sleep(250);
  }
  throw new Error('timeout waiting for: ' + expr);
}
async function idle(extra = 700) {
  await sleep(300);
  await waitFor("!document.documentElement.classList.contains('shiny-busy')", 60000);
  await sleep(extra);
}
async function viewport(height) {
  await send('Emulation.setDeviceMetricsOverride', { width: W, height, deviceScaleFactor: SCALE, mobile: false });
  await sleep(800);
}
// height (CSS px) that shows the visible tab down to its last card
async function pageHeight(min = 900) {
  const h = await evaluate(`Math.ceil(Math.max(...Array.from(document.querySelectorAll(
    '.tab-pane.active .card, .tab-pane.active .accordion')).filter(c => c.offsetParent !== null)
    .map(c => c.getBoundingClientRect().bottom + window.scrollY)) + 16)`);
  return Math.max(min, h);
}
async function shot(name) {
  await evaluate("document.querySelectorAll('.shiny-notification').forEach(n => n.remove()); true");
  // redraw the model diagrams of the visible tab (labels at their final positions)
  await evaluate(`document.querySelectorAll('.tab-pane.active [id^="graph"]').forEach(g => {
    const net = g.chart; if (!net) return;
    const pos = net.getPositions();
    for (const id in pos) net.moveNode(id, pos[id].x, pos[id].y);   // forces edge geometry update
    net.redraw(); net.fit({animation: false}); net.redraw(); }); true`);
  await sleep(800);
  const labels = await evaluate(`Array.from(document.querySelectorAll('.tab-pane.active [id^="graph"]'))
    .filter(g => g.chart).map(g => Object.entries(g.chart.body.edges)
      .map(([id, e]) => id + '@' + Math.round(e.labelModule.size.left) + ',' + Math.round(e.labelModule.size.top))
      .join(' ')).join(' | ')`);
  console.log('  edge labels: ' + labels);
  // card positions (CSS px) for cropping the screenshot later
  const cards = await evaluate(`Array.from(document.querySelectorAll('.navbar, .tab-pane.active .card'))
    .filter(c => c.classList.contains('navbar') || c.offsetParent !== null).map(c => {
      const r = c.getBoundingClientRect(), h = c.querySelector('.card-header');
      return { title: c.classList.contains('navbar') ? 'navbar' : (h ? h.innerText.trim().split('\\n')[0] : ''),
               nested: !!c.parentElement.closest('.card'),
               x: Math.round(r.left), y: Math.round(r.top), w: Math.round(r.width), h: Math.round(r.height) };
    })`);
  writeFileSync(path.join(OUT, name.replace(/\.png$/, '.json')), JSON.stringify({ scale: SCALE, cards }, null, 1));
  const r = await send('Page.captureScreenshot', { format: 'png', fromSurface: true });
  writeFileSync(path.join(OUT, name), Buffer.from(r.data, 'base64'));
  console.log('saved', path.join(OUT, name));
}

await send('Page.enable');
await send('Runtime.enable');
await send('Emulation.setDeviceMetricsOverride', { width: W, height: H, deviceScaleFactor: SCALE, mobile: false });
await send('Page.navigate', { url: APP });
await waitFor('window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected()', 60000);
await idle(1500);
await evaluate(`window.__b = {
  click(txt) {
    const b = Array.from(document.querySelectorAll('button, a')).filter(x => x.offsetParent !== null)
      .find(x => x.textContent.trim() === txt);
    if (!b) throw new Error('no visible button ' + txt); b.click(); return true; },
  sel(id, v) {
    const s = document.getElementById(id); if (!s) throw new Error('no select ' + id);
    if (s.selectize) s.selectize.setValue(v); else { s.value = v; $(s).trigger('change'); } return true; },
  val(id, v) {
    const e = document.getElementById(id); if (!e) throw new Error('no input ' + id);
    $(e).val(v).trigger('change'); return true; },
  tab(txt) {
    const a = Array.from(document.querySelectorAll('a.nav-link')).find(x => x.textContent.trim() === txt);
    if (!a) throw new Error('no tab ' + txt); a.click(); return true; }
}; true`);

// part of the project name shown in the header once an example is loaded
const exampleTitle = {
  CNBio_NSB: 'non-specific binding', CNBio_evaporation: 'with evaporation',
  DCF_gut_liver_DD_100uM: 'daily dose 100', DCF_gut_liver_SD_1000uM: 'single dose 1000'
};
async function loadExample(name) {
  for (let attempt = 1; attempt <= 4; attempt++) {
    await evaluate(`__b.tab('Model structure')`); await idle();
    await evaluate(`__b.sel('structure-example', '${name}')`); await sleep(600);
    await evaluate(`__b.click('Load example')`); await idle(2000);
    const nav = (await evaluate(`document.querySelector('.navbar').innerText`)).replace(/\s+/g, ' ');
    console.log(`load ${name} (attempt ${attempt}): ${nav.slice(-70)}`);
    if (nav.toLowerCase().includes(exampleTitle[name])) return;
    const modal = await evaluate(`(document.querySelector('.modal.show') || {}).innerText || ''`);
    if (modal) console.log('  modal open: ' + modal.replace(/\s+/g, ' ').slice(0, 150));
    const notes = await evaluate(`Array.from(document.querySelectorAll('.shiny-notification'))
      .map(n => n.innerText.replace(/\\s+/g, ' ')).join(' | ')`);
    if (notes) console.log('  notifications: ' + notes.slice(0, 200));
    const state = await evaluate(`JSON.stringify({ connected: Shiny.shinyapp.isConnected(),
      example: Shiny.shinyapp.$inputValues['structure-example'] })`);
    console.log('  state: ' + state);
  }
  throw new Error('could not load example ' + name);
}

// (A) Model structure: NSB model (fluid + solid compartments), solid compartment selected
await viewport(1000);
await loadExample('CNBio_NSB');
await evaluate(`__b.sel('structure-editor-element', 'compartment|NSB_liver')`); await idle(1500);
await viewport(await pageHeight());
await shot('Fig2A_model_structure.png');

// (B) Experimental conditions: diclofenac gut-liver, daily dosing with media exchange
await loadExample('DCF_gut_liver_DD_100uM');
await evaluate(`__b.tab('Experimental conditions')`); await idle(1200);
await evaluate(`__b.sel('conditions-editor-element', 'compartment|Liver')`); await idle(1500);
await viewport(await pageHeight());
await shot('Fig2B_conditions_liver.png');
await evaluate(`__b.sel('conditions-editor-element', 'compartment|Apical')`); await idle(1500);
await viewport(await pageHeight());
await shot('Fig2B_conditions_apical.png');
await viewport(1000);

// (C) Simulation: same virtual experiment with and without evaporation (overlay)
await loadExample('CNBio_evaporation');
await evaluate(`__b.tab('Simulation')`); await idle(1200);
await evaluate(`__b.sel('simulation-comps', ['Liver'])`); await idle(600);
await evaluate(`__b.val('simulation-run_label', 'with evaporation')`); await sleep(400);
await evaluate(`__b.click('Simulate')`); await idle(2500);
await evaluate(`__b.tab('Experimental conditions')`); await idle(1000);
for (const comp of ['Apical', 'Liver']) {
  await evaluate(`__b.sel('conditions-editor-element', 'compartment|${comp}')`); await idle(1000);
  await evaluate(`__b.val('conditions-kev', 0)`); await idle(1000);
}
await evaluate(`__b.tab('Simulation')`); await idle(1000);
await evaluate(`__b.val('simulation-run_label', 'without evaporation')`); await sleep(400);
await evaluate(`__b.click('Simulate')`); await idle(3000);
await viewport(await pageHeight());
await shot('Fig2C_simulation.png');

// (D) Fitting and analysis: diclofenac gut-liver single dose, CL and PS estimated
await loadExample('DCF_gut_liver_SD_1000uM');
await evaluate(`__b.tab('Fitting and analysis')`); await idle(1200);
await evaluate(`__b.click('Fit model')`); await idle(5000);
await viewport(await pageHeight());
await shot('Fig2D_fitting.png');

ws.close();
chrome.kill();
process.exit(0);
