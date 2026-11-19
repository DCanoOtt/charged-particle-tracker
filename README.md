# Charged Particle Tracker

Single-file web app (SRIM-like): multilayer targets, charged-particle transport by condensed-history Monte Carlo
(energy loss from NIST/ICRU tables, straggling, δ-rays, scattering, simple nuclear-reaction model). Units: metric; energies in MeV, lengths in cm internally.

## Layout
| Path | Content |
|---|---|
| `data/` | Generated reference data: `xray.js` (atomic X-ray data, `tools/gen_xray.py`, needs `pip install xraylib`), `elements.js`, `densities.js`, `stopping-tables.js` (PSTAR, ASTAR, ICRU 73/90 from libdedx). Regenerate with `tools/parse_libdedx.py`, then apply `tools/repair_stab.py` (27 ICRU 73 heavy-ion rows of the libdedx data are shifted, mislabelled or duplicated; `tests/test_stab.js` checks the result). |
| `src/engine/` | Physics and Monte Carlo as ES modules (explicit `import`/`export`; `build.py` joins them in name order into the Web Worker source, after checking that every top-level name is unique; `ESBUILD=… node tests/test_engine.js --esm` bundles them with esbuild instead and checks that the results are the same): `10-constants` · `11-materials` · `12-species` · `13-stopping-model` (Bethe/LS, ZBL) · `14-stopping-tables` (table interpolation, data-set choice) · `15-electrons` · `16-nuclear-reactions` · `20-tables` (per-material tables, scattering, hard events) · `21-transport` (event loop, tallies) · `30-exports`. |
| `src/ui/` | Interface, concatenated in name order: `00-bootstrap` (worker) · `01-presets` · `02-formula` · `03-state` · `10-chart` (SVG renderer) · `20-helpers` · `21-panels` · `22-layer-stack` · `23-material-editor` · `24-run` · `30-summary` · `31-plots` (all plot definitions) · `40-zoom` · `41-export` · `43-endf` · `44-alpha-library` · `45-plot-endf` · `50-help` · `90-init`; plus `style.css`, `body.html`. |
| `tests/test_app.py` | Browser regression tests (density, molar mass, transport, flux smoothness, results reset, setup round trip, decay/activation/γ/measurement windows). |
| `tests/test_stab.js`, `test_decay.js`, `test_endf.js`, `test_data.js`, `test_engine.js` | Node tests: audit of the stopping-power tables, decay-data parsers and chains (ENDF-6, ENSDF, X-rays), ENDF cross sections, audit of elements / densities / half-lives / ICRP table, transport engine (seeded runs, reference stopping powers and ranges, energy balance, golden output `tests/golden/engine.json`; `node tests/test_engine.js --update` rewrites it after an intended physics change). |
| `src/help/` | The text of the instructions (HTML), one file per chapter / section of chapter 11; `build.py` joins them into the constant `HELP_HTML`. |
| `src/types.d.ts`, `jsconfig.json` | Types of the central data structures (decay entry, X-ray data, chart options…). Opt-in checking: a module with `// @ts-check` on its first line is checked by `npx tsc -p jsconfig.json` (now `48-decay-xray.js`, `48-decay-ensdf.js`). |
| `build.py` | Assembles `dist/index.html` (complete standalone document), `dist/artifact.html` (bare fragment for the Claude Artifact tool, which adds its own document skeleton) and `dist/engine.js`, `dist/ui.js`; runs `node --check`. |

## Workflow
```
python3 build.py --review   # also dist/review.txt: all the code in one text file (generated tables replaced by placeholders) for a reviewer / AI assistant
python3 build.py            # build + syntax check (--debug also writes dist/index.debug.html with an inline source map: the browser's debugger shows the original files and line numbers)
npx eslint dist/*.js        # optional static check (eslint.config.js)
npx prettier --write "src/**/*.js"   # formatting (.prettierrc.json); the source is kept in this format
python3 tests/test_app.py   # needs Playwright + Chromium
for t in stab decay endf data engine; do node tests/test_$t.js; done
```
All files share one global scope per bundle (engine, UI); the number prefix is the load order. Every module starts with a banner line `/* ===== src/ui/49-gamma.js ===== */` in `dist/*.js` and in the page. Big topics are split over several files (`48-activation*.js`, `49-gamma*.js`): the parts are cut at function boundaries and keep the load order. Add new code to the module that owns the topic, not to the end of a file.

## Model notes
* Electronic stopping: table data set selectable (auto / NIST+ICRU / semi-empirical); log-log Hermite interpolation, √E below the table, Bethe above.
* Nuclear stopping: ZBL; hard recoils are explicit events, the soft part is continuous with variance.
* Step length ≤ 10 % of the range, with a random first step in each layer so that the primaries' step grids do not coincide.
* Ion flux: track-length weighted, each step spread uniformly over its energy interval (exact bin overlap).
* Gas density: ideal gas P·M/(R·T) unless the user disables it.

## Supported cross-section libraries
* **Charged-particle induced reactions (p, d, t, ³He, α):** the program reads TENDL files and also **JENDL-5** and **ENDF/B-VIII.0** files in ENDF-6 format. Both evaluations contain charged-particle sublibraries for only a few isotopes (the light ones). JENDL-5 can be downloaded from <https://wwwndc.jaea.go.jp/jendl/j5/j5.html> and ENDF/B-VIII.0 from <https://www.nndc.bnl.gov/endf-releases/?version=B-VIII.0>. Select the folder with the files in the *Reactions* tab; the file name identifies the projectile and the isotope (for example `d-003_Li_006.endf` for ENDF/B-VIII.0 or `d_003-Li-006.dat` for JENDL-5).
* **Neutron detector model (Configure n-detector & DAQ):** the n+¹H and n+¹²C data are built into the program (`src/ui/51-nd-data.js`, ENDF/B-VIII.0 ACE files 1001.800nc and 6012.800nc at 293.6 K, lib80x); `tools/make_nd_data.js` regenerates them. Above 20 MeV the data are approximate (see the in-program help).
* **Pile-up (Configure n-detector & DAQ → DAQ):** pulse width τ and a paralysable / non-paralysable model give the pile-up level and the count rate after pile-up in the counting-statistics table.
* **PPAC detector with a U-238 target (Configure n-detector & DAQ → Detector type):** the U-238 neutron-induced fission cross section is built into the program (`src/ui/51-nd-data-u.js`, JENDL-5 `U-238-n.jendl5`, 300 K, up to 195 MeV, constant above); regenerate it with `python3 tools/make_nd_u238.py U-238-n.jendl5`. JENDL-5 can be downloaded from <https://wwwndc.jaea.go.jp/jendl/j5/j5.html>.

## Alpha reaction data (TENDL, ENDF-6)
* `src/ui/43-endf.js` reads MF 3 (cross sections) and the product yields of MF 6; distributions are skipped by record structure. `44-alpha-library.js` handles the folder (File System Access API, directory input as fallback), loads the isotopes of the target and folds the flux; `45-plot-endf.js` is the tab.
* The engine tallies the path length of the primary ions per layer and energy bin (`result.Lh`, 1000 bins up to the maximum beam energy); folding with σ(E) is done in the UI, so any channel can be inspected without a new run.
* `tests/make_endf.py` writes synthetic ENDF-6 files for the tests (not evaluated data); `tests/test_endf.js` tests the parser.
