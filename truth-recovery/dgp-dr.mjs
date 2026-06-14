// dgp-dr.mjs — STANDALONE seeded known-truth dose-response DGP.
//
// Models the TWO-STAGE pooling path of metasprint-dose-response:
//   Each study i has a true per-unit-dose log-RR slope  slope_i ~ N(muSlope, tau2).
//   The study runs a binomial dose-response trial: a reference arm at dose 0 plus
//   several dose arms, events a ~ Bin(n, p0 * exp(slope_i * dose)).
//   From those counts we fit a within-study weighted log-linear slope + its SE,
//   and emit one row in the app's two-stage format: { slope, slopeSE } (effect scale).
//
// The estimand is the TRUE POOLED MEAN SLOPE  muSlope  (and the between-study tau2).
// computeMetaAnalysis is then handed these per-study slopes as effectType 'MD'
// (additive scale — slopes are already on the log-RR-per-unit-dose scale).

// ---- deterministic RNG (mulberry32) ----
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function randn(rng) {
  // Box-Muller
  let u = 0, v = 0;
  while (u === 0) u = rng();
  while (v === 0) v = rng();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
function binom(rng, n, p) {
  if (p <= 0) return 0;
  if (p >= 1) return n;
  let k = 0;
  for (let i = 0; i < n; i++) if (rng() < p) k++;
  return k;
}

const Z975 = 1.959963984540054; // qnorm(0.975)

/**
 * Generate one synthetic study's two-stage row by fitting a within-study
 * weighted log-linear dose-response slope from binomial counts.
 *
 * @returns {{slope, slopeSE, lowerCI, upperCI, effectEstimate, effectType, trueSlope}|null}
 */
function makeStudy(rng, { muSlope, tau2, doses, p0, n }) {
  const slope_i = muSlope + (tau2 > 0 ? Math.sqrt(tau2) * randn(rng) : 0);
  // reference at dose 0 plus the supplied dose arms
  const allDoses = [0, ...doses];
  const pts = [];
  for (const d of allDoses) {
    const p = Math.min(0.999, Math.max(0.001, p0 * Math.exp(slope_i * d)));
    let a = binom(rng, n, p);
    a = Math.min(n - 1, Math.max(1, a)); // avoid 0/n log-RR blow-ups (Haldane-ish)
    const phat = a / n;
    // log-RR vs reference handled by WLS on logit-free log-risk: use log(phat)
    const logp = Math.log(phat);
    // variance of log(phat) ~ (1-phat)/(n*phat)
    const v = (1 - phat) / (n * phat);
    pts.push({ dose: d, y: logp, v });
  }
  // within-study WLS slope of log-risk on dose  ->  per-unit-dose log-RR slope
  let sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;
  for (const pt of pts) {
    const w = 1 / pt.v;
    sw += w; swx += w * pt.dose; swy += w * pt.y;
    swxx += w * pt.dose * pt.dose; swxy += w * pt.dose * pt.y;
  }
  const det = sw * swxx - swx * swx;
  if (Math.abs(det) < 1e-15) return null;
  const slope = (sw * swxy - swx * swy) / det;
  const seSlope = Math.sqrt(sw / det);
  if (!isFinite(slope) || !isFinite(seSlope) || seSlope <= 0) return null;

  // Emit in the app's effect-with-CI format (additive 'MD' scale: the slope itself)
  return {
    trueSlope: slope_i,
    slope, slopeSE: seSlope,
    effectEstimate: slope,
    lowerCI: slope - Z975 * seSlope,
    upperCI: slope + Z975 * seSlope,
    effectType: 'MD',
  };
}

/**
 * One simulated meta-analysis: k studies, true pooled mean slope muSlope,
 * between-study heterogeneity tau2.
 */
// NOTE on dose scale: the count model is p = p0*exp(slope*dose). To stay on the
// non-saturated part of the curve (where a log-linear within-study slope is a clean
// estimand) we keep slope*maxDose well below ~0.7. With muSlope=0.10 the default
// doses [1,2,4] give p0*exp(0.4)=0.30 at the top dose — no clamping, clean recovery.
function makeMetaAnalysis(seed, { k = 12, muSlope = 0.10, tau2 = 0.0,
  doses = [1, 2, 4], p0 = 0.20, n = 400 } = {}) {
  const rng = mulberry32(seed);
  const studies = [];
  let tries = 0;
  while (studies.length < k && tries < k * 10) {
    tries++;
    const s = makeStudy(rng, { muSlope, tau2, doses, p0, n });
    if (s) studies.push(s);
  }
  return { studies, truth: { muSlope, tau2, k } };
}

export { mulberry32, randn, makeStudy, makeMetaAnalysis };
