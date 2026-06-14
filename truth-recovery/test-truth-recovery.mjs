// test-truth-recovery.mjs — node --test
// Validates the metasprint-dose-response TWO-STAGE slope-pooling engine against
// a known-truth binomial dose-response DGP.
//
// Run:  node --test truth-recovery/test-truth-recovery.mjs

import { test } from 'node:test';
import assert from 'node:assert/strict';
import { runScenario } from './harness.mjs';
import { computeMetaAnalysis, estimateREML } from './engine.mjs';
import { makeMetaAnalysis } from './dgp-dr.mjs';

// 1. Point estimate of the pooled mean slope is unbiased (homogeneous case).
test('pooled slope is unbiased under homogeneity', () => {
  const r = runScenario({ label: 'homog', tau2: 0.0, reps: 400, muSlope: 0.10 });
  assert.ok(Math.abs(r.bias) < 0.01, `bias ${r.bias} should be ~0`);
});

// 2. CRITICAL (dose-response-pro bug check): the DL tau2 estimator must NOT
//    collapse to 0 under genuine between-study heterogeneity. It must track truth.
test('DL tau2 recovers true heterogeneity (no silent collapse to 0)', () => {
  const r = runScenario({ label: 'het', tau2: 0.005, reps: 400, muSlope: 0.10 });
  // mean DL tau2 must be within 30% of the true value, NOT ~0
  assert.ok(r.tau2_DL_mean > 0.0035 && r.tau2_DL_mean < 0.0065,
    `DL tau2 mean ${r.tau2_DL_mean} should be near true 0.005, not collapsed`);
  // and it must not collapse to exactly 0 in most reps
  assert.ok(r.fracTau2Zero < 0.05,
    `frac(tau2==0) ${r.fracTau2Zero} too high — would indicate silent collapse`);
});

// 3. Coverage of the pooled CI is near-nominal under homogeneity (DGP/engine sane).
test('coverage near nominal 0.95 under homogeneity', () => {
  const r = runScenario({ label: 'homog', tau2: 0.0, reps: 500, muSlope: 0.10 });
  assert.ok(r.coverage > 0.92 && r.coverage <= 0.99,
    `coverage ${r.coverage} should be ~0.95`);
});

// 4. Under heterogeneity, DL coverage degrades (known DL under-coverage) but is not
//    catastrophic, and HKSJ moves coverage UP toward nominal (repair direction).
test('HKSJ improves coverage over plain DL under heterogeneity', () => {
  const dl = runScenario({ label: 'het-DL', tau2: 0.005, reps: 500, muSlope: 0.10, useHKSJ: false });
  const hk = runScenario({ label: 'het-HKSJ', tau2: 0.005, reps: 500, muSlope: 0.10, useHKSJ: true });
  assert.ok(dl.coverage > 0.80, `DL coverage ${dl.coverage} not catastrophic`);
  assert.ok(hk.coverage >= dl.coverage,
    `HKSJ coverage ${hk.coverage} should be >= DL ${dl.coverage}`);
});

// 5. REML tau2 also recovers truth (alternative estimator sanity).
test('REML tau2 recovers true heterogeneity', () => {
  const { studies } = makeMetaAnalysis(7, { k: 20, muSlope: 0.10, tau2: 0.004 });
  const data = studies.map(s => ({ yi: s.slope, vi: s.slopeSE * s.slopeSE }));
  const t2 = estimateREML(data);
  assert.ok(t2 > 0.0015 && t2 < 0.0090, `REML tau2 ${t2} should be near 0.004`);
});
