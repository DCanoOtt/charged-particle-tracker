/* ====== Neutron kinematics, time-of-flight folding and the spectrum reconstructed from the time of flight (pure functions, no DOM) ====== */
const NEU_MC2 = 939.56542,
  NEU_C = 0.299792458; // neutron rest energy (MeV), speed of light (m/ns)
// time of flight (ns) and |dE/dt| (MeV/ns) of a neutron of kinetic energy E (MeV) over L metres; sigma_E of a time resolution st (ns, 1 sigma)
function neuKin(E, L) {
  const g = 1 + E / NEU_MC2,
    b = Math.sqrt(1 - 1 / (g * g)),
    t = L / (b * NEU_C);
  return { t, dEdt: (NEU_MC2 * g * g * g * b * b) / t };
}
const neuSigE = (E, L, st) => neuKin(E, L).dEdt * st;
// Gaussian smearing of a spectrum y (per MeV, nbo bins up to Eo) with the energy-dependent width of the time resolution; the integral is kept
function neuFold(y, nbo, Eo, L, st) {
  const dE = Eo / nbo,
    out = new Float64Array(nbo);
  for (let i = 0; i < nbo; i++) {
    if (!(y[i] > 0)) continue;
    const sg = neuSigE((i + 0.5) * dE, L, st);
    if (!(sg > dE / 2)) {
      out[i] += y[i];
      continue;
    }
    const j0 = Math.max(0, Math.floor(i - (5 * sg) / dE)),
      j1 = Math.min(nbo - 1, Math.ceil(i + (5 * sg) / dE)),
      w = [];
    let sum = 0;
    for (let j = j0; j <= j1; j++) {
      const v = Math.exp((-0.5 * ((j - i) * dE) ** 2) / (sg * sg));
      w.push(v);
      sum += v;
    }
    for (let j = j0; j <= j1; j++) out[j] += (y[i] * w[j - j0]) / sum;
  }
  return out;
}
// normal distribution function Φ(x) (erfc by the Chebyshev fit of Numerical Recipes, relative error < 1.2e-7)
function neuPhi(x) {
  const z = Math.abs(x) / Math.SQRT2,
    t = 1 / (1 + 0.5 * z),
    ec =
      t *
      Math.exp(
        -z * z -
          1.26551223 +
          t *
            (1.00002368 +
              t *
                (0.37409196 +
                  t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277))))))))
      );
  return x >= 0 ? 1 - 0.5 * ec : 0.5 * ec;
}
// kinetic energy (MeV) of a neutron that needs the time t (ns) for L metres; Infinity when t <= L/c
function neuEfromT(t, L) {
  const b = L / (NEU_C * t);
  return b >= 1 ? Infinity : b <= 0 ? 0 : NEU_MC2 * (1 / Math.sqrt(1 - b * b) - 1);
}
// Expected spectrum reconstructed from the time of flight. cnt = expected detected neutrons in the nbo bins (0 - Eo) of the true energy;
// every neutron gets the time t(E) over L metres + a Gaussian time error (st = 1 sigma, ns) and is put in the bin of E' = E(t') of nb bins (0 - Eo).
// Per source bin, 8 energies are taken and the probability of every reconstructed bin is the integral of the Gaussian between the times of its edges.
// Returns { mu: expected counts per reconstructed bin, lost: counts reconstructed above Eo }
function neuTofRecon(cnt, Eo, L, st, nb) {
  const nbo = cnt.length,
    dE = Eo / nbo,
    dO = Eo / nb,
    mu = new Float64Array(nb),
    te = new Float64Array(nb + 1); // times of the edges of the reconstructed bins (ns), decreasing with the bin index
  let tot = 0,
    got = 0;
  te[0] = Infinity;
  for (let i = 1; i <= nb; i++) te[i] = neuKin(i * dO, L).t;
  const NS = 8;
  for (let j = 0; j < nbo; j++) {
    if (!(cnt[j] > 0)) continue;
    tot += cnt[j];
    for (let k = 0; k < NS; k++) {
      const w = cnt[j] / NS,
        t0 = neuKin(Math.max((j + (k + 0.5) / NS) * dE, 1e-9), L).t;
      if (!(st > 0)) {
        const i = Math.min(nb, Math.floor(((j + (k + 0.5) / NS) * dE) / dO));
        if (i < nb) {
          mu[i] += w;
          got += w;
        }
        continue;
      }
      const Elo = neuEfromT(t0 + 6 * st, L),
        Ehi = neuEfromT(t0 - 6 * st, L),
        i0 = Math.max(0, Math.floor(Elo / dO)),
        i1 = Math.min(nb - 1, Number.isFinite(Ehi) ? Math.ceil(Ehi / dO) : nb - 1);
      for (let i = i0; i <= i1; i++) {
        const a = Number.isFinite(te[i]) ? neuPhi((te[i] - t0) / st) : 1,
          b = neuPhi((te[i + 1] - t0) / st),
          p = a - b;
        if (p > 0) {
          mu[i] += w * p;
          got += w * p;
        }
      }
    }
  }
  return { mu, lost: tot - got };
}
// Poisson random number with mean m (Knuth for small m, normal approximation above 40)
function neuPoisson(m, rng) {
  if (!(m > 0)) return 0;
  if (m < 40) {
    const L = Math.exp(-m);
    let k = 0,
      p = 1;
    do {
      k++;
      p *= rng();
    } while (p > L);
    return k - 1;
  }
  const u = Math.max(rng(), 1e-300),
    z = Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng());
  return Math.max(0, Math.round(m + Math.sqrt(m) * z));
}
