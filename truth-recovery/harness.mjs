// harness.mjs — wires the repo's OWN computeMetaAnalysis (two-stage slope pooling)
// to the known-truth DGP and measures:
//   (1) coverage of the true pooled mean slope by the pooled CI,
//   (2) bias of the pooled slope point estimate,
//   (3) whether the DL tau2 estimate collapses to 0 under known heterogeneity,
//   (4) whether REML recovers tau2.
//
// Run:  node truth-recovery/harness.mjs

import { computeMetaAnalysis, estimateREML, normalQuantile } from './engine.mjs';
import { makeMetaAnalysis } from './dgp-dr.mjs';

function runScenario({ label, tau2, k = 12, reps = 600, useHKSJ = false, muSlope = 0.10 }) {
  let covered = 0, biasSum = 0, tau2Sum = 0, tau2RemlSum = 0, tau2Zero = 0, valid = 0;
  let widthSum = 0;
  for (let r = 0; r < reps; r++) {
    const { studies, truth } = makeMetaAnalysis(1000 + r * 7919, { k, muSlope, tau2 });
    const res = computeMetaAnalysis(studies, 0.95, { hksj: useHKSJ });
    if (!res) continue;
    valid++;
    // pooled is on the additive ('MD') slope scale (isRatio=false)
    const lo = res.pooledLo, hi = res.pooledHi, est = res.pooled;
    if (truth.muSlope >= lo && truth.muSlope <= hi) covered++;
    biasSum += est - truth.muSlope;
    widthSum += hi - lo;
    tau2Sum += res.tau2;
    tau2RemlSum += res.tau2REML;
    if (res.tau2 < 1e-9) tau2Zero++;
  }
  return {
    label, tau2True: tau2, k, reps: valid,
    coverage: covered / valid,
    bias: biasSum / valid,
    meanWidth: widthSum / valid,
    tau2_DL_mean: tau2Sum / valid,
    tau2_REML_mean: tau2RemlSum / valid,
    fracTau2Zero: tau2Zero / valid,
    method: useHKSJ ? 'DL+HKSJ' : 'DL',
  };
}

function fmt(x, d = 4) { return Number(x).toFixed(d); }

const scenarios = [
  { label: 'Homogeneous (tau2=0)',        tau2: 0.0,    muSlope: 0.10 },
  { label: 'Mild het (tau2=0.0005)',      tau2: 0.0005, muSlope: 0.10 },
  { label: 'Moderate het (tau2=0.002)',   tau2: 0.002,  muSlope: 0.10 },
  { label: 'Strong het (tau2=0.005)',     tau2: 0.005,  muSlope: 0.10 },
];

console.log('=== metasprint-dose-response: two-stage slope pooling truth-recovery ===');
console.log('Estimand: true pooled mean per-unit-dose log-RR slope (muSlope=0.10), k=12, reps=600');
console.log('');
const results = [];
for (const sc of scenarios) {
  const dl = runScenario({ ...sc, useHKSJ: false });
  results.push(dl);
  console.log(`[${dl.label}]  method=${dl.method}`);
  console.log(`   coverage(95% target) = ${fmt(dl.coverage, 3)}   bias = ${fmt(dl.bias)}   meanCIwidth = ${fmt(dl.meanWidth)}`);
  console.log(`   tau2_true = ${fmt(dl.tau2True, 5)}   tau2_DL_mean = ${fmt(dl.tau2_DL_mean, 5)}   tau2_REML_mean = ${fmt(dl.tau2_REML_mean, 5)}   frac(tau2_DL==0) = ${fmt(dl.fracTau2Zero, 3)}`);
}

// HKSJ on strongest heterogeneity, to compare coverage repair
console.log('');
const hk = runScenario({ label: 'Strong het + HKSJ', tau2: 0.005, muSlope: 0.10, useHKSJ: true });
console.log(`[${hk.label}]  method=${hk.method}`);
console.log(`   coverage(95% target) = ${fmt(hk.coverage, 3)}   bias = ${fmt(hk.bias)}   meanCIwidth = ${fmt(hk.meanWidth)}`);

export { runScenario };
