/* Types of the central data structures (documentation + editor/tsc support; not part of the build).
   Opt in per file with a first-line `// @ts-check` and JSDoc (@param / @returns), then check with:  npx tsc -p jsconfig.json  */

/** Nuclide id: 10·(1000·Z + A) + isomer index (0 = ground state). */
type NuclideId = number;

/** One decay mode of a nuclide. */
interface DecayBranch {
  /** fraction of the decays that go this way (the fractions of an entry add up to 1) */
  f: number;
  /** daughter id (see NuclideId); 0 = no daughter (spontaneous fission) */
  d: NuclideId;
  /** label, e.g. 'β−', 'EC/β+', 'IT', 'α' */
  lab: string;
  /** β+ share of an EC/β+ branch, where the data give it (ENSDF E records) */
  pp?: number;
}

/** Photon line: [energy (keV), photons per decay, optional flag (1 = annihilation)]. */
type PhotonLine = [number, number] | [number, number, number];

/** Decay data of one nuclide, as the readers (internal table, ENDF-6, ENSDF) deliver them. */
interface DecayEntry {
  st: 'stable' | 'rad';
  /** half-life (s) */
  T: number | null;
  br: DecayBranch[];
  /** discrete γ lines */
  gam?: PhotonLine[];
  /** total γ intensity per decay */
  gtot?: number;
  /** X-ray lines [energy keV, photons per decay] and their sum */
  xr?: [number, number][];
  xtot?: number;
  /** annihilation photons per decay */
  ann?: number;
  /** continuous γ spectrum: energies (keV) and photons per decay and keV */
  gc?: { E: number[]; y: number[] };
}

/** Atomic data of one element (data/xray.js). */
interface XrayElement {
  /** K, L1, L2, L3 edge energies (keV) */
  e: number[];
  /** fluorescence yields ωK, ωL1, ωL2, ωL3 */
  w: number[];
  /** Coster–Kronig probabilities f12, f13, f23 */
  f: number[];
  /** radiative yields of K→L2 and K→L3 per K vacancy */
  kl: number[];
  /** lines of each shell: [energy (keV), photons per vacancy] */
  K: [number, number][]; L1: [number, number][]; L2: [number, number][]; L3: [number, number][];
}
declare const XRAY: Record<number, XrayElement>;

/** Vacancies per decay in the K and L shells. */
interface Vacancies { K: number; L1: number; L2: number; L3: number }

/** A curve of renderChart. */
interface ChartSeries {
  name: string;
  type?: 'line' | 'hist' | 'pts' | 'bar';
  x?: ArrayLike<number>; y: ArrayLike<number>; edges?: ArrayLike<number>;
  color: string; w?: number; dash?: string; hidden?: boolean;
}
/** Options of renderChart (all numbers as strings or numbers; '' = automatic). */
interface ChartOptions {
  series: ChartSeries[];
  xlabel: string; ylabel: string; title?: string; stats?: string[];
  xlog?: boolean; ylog?: boolean;
  xmin?: string | number; xmax?: string | number; ymin?: string | number; ymax?: string | number;
  vlines?: { x: number; label?: string; color?: string; dash?: string; w?: number; op?: number; rot?: boolean }[];
  bands?: { x0: number; x1: number; label?: string }[];
}
/** Geometry returned by renderChart (pixel margins and data ranges), used by zoom, tooltips and exports. */
interface ChartGeom { W: number; H: number; m: { l: number; r: number; t: number; b: number }; iw: number; ih: number; x0: number; x1: number; y0: number; y1: number; xlog: boolean; ylog: boolean }
