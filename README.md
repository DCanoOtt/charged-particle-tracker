# Charged Particle Tracker

Single-file web app (SRIM-like): multilayer targets, charged-particle transport by condensed-history Monte Carlo
(energy loss from NIST/ICRU tables, straggling, δ-rays, scattering, simple nuclear-reaction model). Units: metric; energies in MeV, lengths in cm internally.

## Quick start
* **Web / offline:** open `dist/index.html` (built with `python3 build.py`) in any current Chromium-based browser or Firefox. It is one self-contained file; nothing is uploaded anywhere. Settings are kept in the browser (localStorage); use *Save setup* for a portable file.
* **GitHub Pages:** the workflow `.github/workflows/pages.yml` publishes the build on every push to `main` (enable *Settings → Pages → Source: GitHub Actions*).
* **Desktop app (Windows, macOS, Linux):** Electron wrapper in `desktop/`.
  ```
  python3 build.py
  cd desktop && npm install && npm start      # run
  npm run dist                                # installers in desktop/release/
  ```
  Pushing a tag (`git tag v0.1.0 && git push --tags`) runs `.github/workflows/release.yml`, which builds the installers on the three systems and attaches them to a GitHub release. The installers are not code-signed (Windows SmartScreen / macOS Gatekeeper will warn).

## Optional data you supply
* TENDL / ENDF-6 charged-particle reaction libraries (`TENDL-a`, `-p`, `-d`, `-t`, `-h`) and ENDF-6 decay data (MF 8 / MT 457) are chosen from your disk in the app and read locally; they are not part of this repository.

## Data and references bundled in the code
Stopping-power tables (NIST PSTAR/ASTAR, ICRU 73/90) via libdedx; decay data from ICRP-107 with NUBASE2020 half-lives (radioactivedecay 0.6.1) and mendeleev 1.3.0. Check the licences of these sources before choosing the licence of the repository (see `tools/` and `tools_gen_*.py` for how the tables were generated). The *Help* tab of the app lists the references.

## Status
Validated with synthetic data and regression tests (`tests/`); not yet compared systematically against SRIM/PSTAR/ASTAR or real TENDL/ENDF decay libraries.

## Layout
| Path | Content |
|---|---|
| `data/` | Generated reference data: `elements.js`, `densities.js`, `stopping-tables.js` (PSTAR, ASTAR, ICRU 73/90 from libdedx). Regenerate with `tools/parse_libdedx.py`. |
| `src/engine/` | Physics and Monte Carlo, concatenated in name order into the Web Worker source: `10-constants` · `11-materials` · `12-species` · `13-stopping-model` (Bethe/LS, ZBL) · `14-stopping-tables` (table interpolation, data-set choice) · `15-electrons` · `16-nuclear-reactions` · `20-tables` (per-material tables, scattering, hard events) · `21-transport` (event loop, tallies) · `30-exports`. |
| `src/ui/` | Interface, concatenated in name order: `00-bootstrap` (worker) · `01-presets` · `02-formula` · `03-state` · `10-chart` (SVG renderer) · `20-helpers` · `21-panels` · `22-layer-stack` · `23-material-editor` · `24-run` · `30-summary` · `31-plots` (all plot definitions) · `40-zoom` · `41-export` · `43-endf` · `44-alpha-library` · `45-plot-endf` · `50-help` · `90-init`; plus `style.css`, `body.html`. |
| `tests/test_app.py` | Regression tests (density, molar mass, transport, flux smoothness, results reset, setup round trip). |
| `build.py` | Assembles `dist/index.html` (and `dist/engine.js`, `dist/ui.js`), runs `node --check`. |

## Workflow
```
python3 build.py            # build + syntax check
npx eslint dist/*.js        # optional static check (eslint.config.js)
python3 tests/test_app.py   # needs Playwright + Chromium
```
All files share one global scope per bundle (engine, UI); the number prefix is the load order. Add new code to the module that owns the topic, not to the end of a file.

## Model notes
* Electronic stopping: table data set selectable (auto / NIST+ICRU / semi-empirical); log-log Hermite interpolation, √E below the table, Bethe above.
* Nuclear stopping: ZBL; hard recoils are explicit events, the soft part is continuous with variance.
* Step length ≤ 10 % of the range, with a random first step in each layer so that the primaries' step grids do not coincide.
* Ion flux: track-length weighted, each step spread uniformly over its energy interval (exact bin overlap).
* Gas density: ideal gas P·M/(R·T) unless the user disables it.

## Alpha reaction data (TENDL, ENDF-6)
* `src/ui/43-endf.js` reads MF 3 (cross sections) and the product yields of MF 6; distributions are skipped by record structure. `44-alpha-library.js` handles the folder (File System Access API, directory input as fallback), loads the isotopes of the target and folds the flux; `45-plot-endf.js` is the tab.
* The engine tallies the path length of the primary ions per layer and energy bin (`result.Lh`, 1000 bins up to the maximum beam energy); folding with σ(E) is done in the UI, so any channel can be inspected without a new run.
* `tests/make_endf.py` writes synthetic ENDF-6 files for the tests (not evaluated data); `tests/test_endf.js` tests the parser.
