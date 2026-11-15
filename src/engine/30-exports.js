/* global ELTAB, ELDENS */
// Entry of the engine: build.py joins the modules (imports and exports removed; every top-level name is unique) and the page uses the object PT.
import { IVAL, SYMS, SYMZ, isoMass } from './10-constants.js';
import { buildMaterial } from './11-materials.js';
import { makeSpecies } from './12-species.js';
import { tmaxHeavy } from './13-stopping-model.js';
import { sigmaReaction } from './16-nuclear-reactions.js';
import { simulate, stopTable } from './21-transport.js';
export const PT = { buildMaterial, makeSpecies, simulate, stopTable, ELTAB, ELDENS, IVAL, SYMS, SYMZ, isoMass, sigmaReaction, tmaxHeavy };
