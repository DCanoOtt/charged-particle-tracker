/* ====== material presets ====== */
const el = (sym, amt) => ({ kind: 'el', sym, amt });
const PRESETS = [
  // solids
  { g: 'Solids', name: 'Silicon', phase: 'solid', mode: 'atoms', rows: [el('Si', 1)] },
  { g: 'Solids', name: 'Aluminium', phase: 'solid', mode: 'atoms', rows: [el('Al', 1)] },
  { g: 'Solids', name: 'Titanium', phase: 'solid', mode: 'atoms', rows: [el('Ti', 1)] },
  { g: 'Solids', name: 'Iron', phase: 'solid', mode: 'atoms', rows: [el('Fe', 1)] },
  { g: 'Solids', name: 'Nickel', phase: 'solid', mode: 'atoms', rows: [el('Ni', 1)] },
  { g: 'Solids', name: 'Copper', phase: 'solid', mode: 'atoms', rows: [el('Cu', 1)] },
  { g: 'Solids', name: 'Tungsten', phase: 'solid', mode: 'atoms', rows: [el('W', 1)] },
  { g: 'Solids', name: 'Gold', phase: 'solid', mode: 'atoms', rows: [el('Au', 1)] },
  { g: 'Solids', name: 'Lead', phase: 'solid', mode: 'atoms', rows: [el('Pb', 1)] },
  { g: 'Solids', name: 'Carbon (graphite)', phase: 'solid', mode: 'atoms', rows: [el('C', 1)], rho: 2.2, I: 78 },
  { g: 'Solids', name: 'Beryllium', phase: 'solid', mode: 'atoms', rows: [el('Be', 1)] },
  { g: 'Solids', name: 'Lithium fluoride (LiF)', phase: 'solid', mode: 'atoms', rows: [el('Li', 1), el('F', 1)], rho: 2.635, I: 94 },
  { g: 'Solids', name: 'Silicon dioxide (SiO2)', phase: 'solid', mode: 'atoms', rows: [el('Si', 1), el('O', 2)], rho: 2.32, I: 139.2 },
  { g: 'Solids', name: 'Aluminium oxide (Al2O3)', phase: 'solid', mode: 'atoms', rows: [el('Al', 2), el('O', 3)], rho: 3.97, I: 145.2 },
  { g: 'Solids', name: 'Silicon nitride (Si3N4)', phase: 'solid', mode: 'atoms', rows: [el('Si', 3), el('N', 4)], rho: 3.17 },
  { g: 'Solids', name: 'Sodium iodide (NaI)', phase: 'solid', mode: 'atoms', rows: [el('Na', 1), el('I', 1)], rho: 3.667, I: 452 },
  { g: 'Solids', name: 'Caesium iodide (CsI)', phase: 'solid', mode: 'atoms', rows: [el('Cs', 1), el('I', 1)], rho: 4.51, I: 553.1 },
  { g: 'Solids', name: 'BGO (Bi4Ge3O12)', phase: 'solid', mode: 'atoms', rows: [el('Bi', 4), el('Ge', 3), el('O', 12)], rho: 7.13, I: 534.1 },
  { g: 'Solids', name: 'Polyethylene (CH2)', phase: 'solid', mode: 'atoms', rows: [el('C', 1), el('H', 2)], rho: 0.94, I: 57.4 },
  { g: 'Solids', name: 'PMMA (C5H8O2)', phase: 'solid', mode: 'atoms', rows: [el('C', 5), el('H', 8), el('O', 2)], rho: 1.19, I: 74 },
  { g: 'Solids', name: 'Polystyrene (C8H8)', phase: 'solid', mode: 'atoms', rows: [el('C', 8), el('H', 8)], rho: 1.06, I: 68.7 },
  { g: 'Solids', name: 'Kapton (C22H10N2O5)', phase: 'solid', mode: 'atoms', rows: [el('H', 10), el('C', 22), el('N', 2), el('O', 5)], rho: 1.42, I: 79.6 },
  { g: 'Solids', name: 'Mylar (C10H8O4)', phase: 'solid', mode: 'atoms', rows: [el('C', 10), el('H', 8), el('O', 4)], rho: 1.38, I: 78.7 },
  { g: 'Solids', name: 'Teflon (C2F4)', phase: 'solid', mode: 'atoms', rows: [el('C', 2), el('F', 4)], rho: 2.2, I: 99.1 },
  { g: 'Solids', name: 'Stainless steel 304', phase: 'solid', mode: 'mass', rows: [el('Fe', 69.5), el('Cr', 19), el('Ni', 9.5), el('Mn', 2)], rho: 8.0 },
  // liquids
  { g: 'Liquids', name: 'Water', phase: 'liquid', mode: 'atoms', rows: [el('H', 2), el('O', 1)], rho: 1.0, I: 78 },
  { g: 'Liquids', name: 'Heavy water (D2O)', phase: 'liquid', mode: 'atoms', rows: [el('H', 2), el('O', 1)], rho: 1.107, I: 78, iso: { H: [[1, 0], [2, 100]] } },
  { g: 'Liquids', name: 'Liquid argon', phase: 'liquid', mode: 'atoms', rows: [el('Ar', 1)], rho: 1.396 },
  { g: 'Liquids', name: 'Liquid xenon', phase: 'liquid', mode: 'atoms', rows: [el('Xe', 1)], rho: 2.94 },
  { g: 'Liquids', name: 'Liquid hydrogen', phase: 'liquid', mode: 'atoms', rows: [el('H', 2)], rho: 0.0708, I: 21.8 },
  { g: 'Liquids', name: 'Ethanol (C2H5OH)', phase: 'liquid', mode: 'atoms', rows: [el('C', 2), el('H', 6), el('O', 1)], rho: 0.789, I: 62.9 },
  // gases
  { g: 'Gases', name: 'Air (dry)', phase: 'gas', mode: 'mass', rows: [el('N', 75.5268), el('O', 23.1781), el('Ar', 1.2827), el('C', 0.0124)], molMass: 28.96, I: 85.7, gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Hydrogen gas (H2)', phase: 'gas', mode: 'atoms', rows: [el('H', 2)], I: 19.2, gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Helium gas', phase: 'gas', mode: 'atoms', rows: [el('He', 1)], gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Nitrogen gas (N2)', phase: 'gas', mode: 'atoms', rows: [el('N', 2)], I: 82, gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Argon gas', phase: 'gas', mode: 'atoms', rows: [el('Ar', 1)], gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Neon gas', phase: 'gas', mode: 'atoms', rows: [el('Ne', 1)], gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Xenon gas', phase: 'gas', mode: 'atoms', rows: [el('Xe', 1)], gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Methane (CH4)', phase: 'gas', mode: 'atoms', rows: [el('C', 1), el('H', 4)], I: 41.7, gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Isobutane (C4H10)', phase: 'gas', mode: 'atoms', rows: [el('C', 4), el('H', 10)], I: 48.3, gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Carbon dioxide (CO2)', phase: 'gas', mode: 'atoms', rows: [el('C', 1), el('O', 2)], I: 85, gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'P10 (Ar 90% / CH4 10%)', phase: 'gas', mode: 'mass', rows: [el('Ar', 95.73), el('C', 3.198), el('H', 1.074)], molMass: 37.557, gas: { P: 101.325, T: 293.15 } },
  { g: 'Gases', name: 'Deuterium gas (D2)', phase: 'gas', mode: 'atoms', rows: [el('H', 2)], I: 19.2, iso: { H: [[1, 0], [2, 100]] }, gas: { P: 101.325, T: 293.15 } },
];
PRESETS.forEach(p => { if (p.rho == null && p.phase !== 'gas') delete p.rho; });
let USERLIB = [];       // user-saved materials (defs)
function libMap() { const m = {}; PRESETS.forEach(p => m[p.name] = p); USERLIB.forEach(p => m[p.name] = p); return m; }
function presetDef(name) { const p = libMap()[name]; const d = clone(p); delete d.g; return d; }
