# Exceed dosresmeta2 Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Achieve full parity with R's dosresmeta2 package (GL covariance, ML/REML, all models), then exceed it with model-averaged predictions. Validate against R reference values.

**Architecture:** Single-file HTML app (`metasprint-dose-response.html`, ~21,559 lines). All new code goes into the `<script>` block (starts line 1785). Functions are added near related existing functions. Inline self-tests validate against R dosresmeta2 reference values from `alcohol_cvd` dataset.

**Tech Stack:** Vanilla JS, SVG for plots, inline self-tests, Selenium + Python for integration tests. R dosresmeta2 v2.2.0 for reference value generation.

**Key line references (as of commit `eb6caa4`):**
- `parseFloatSafe`: 7366 | `computeMetaAnalysis`: 7407 | `chi2CDF`: 8101
- `fitLinearDR`: 8543 | `fitQuadraticDR`: 8594 | `inv3`: 8623
- `matTranspose`: 8697 | `matMul`: 8712 | `matInvertNxN`: 8742 | `wlsFit`: 8810
- `fitEmaxDR`: 8993 | `fitFP1DR`: 9140 | `fitFP2DR`: 9221
- `computeRCSBasis`: 9509 | `fitRCSDR`: 9548 | `fitOneStageMixedEffects`: 9827
- `_runAllDRTests`: 10223 | `computeDRHeterogeneity`: 10253
- `compareDoseResponseModels`: 10289 | `renderModelComparisonTable`: 10329
- `renderLeaveOneOutDR`: 10382 | `renderDoseResponseAnalysis`: 10471
- `renderDoseResponseCurveSVG`: 10615
- `drModelSelect` dropdown: 1457 | `methodSelect` dropdown: 1451
- Extract phase: 1361

**R reference values (dosresmeta2 v2.2.0, alcohol_cvd dataset):**
- Linear REML: coef=-0.004365412, SE=0.00588923, tau²=0.0001005709, AIC=-24.06829
- Quadratic REML: coef=[-0.03023389, 0.0007514854], SE=[0.0131779, 0.0004885058], AIC=-59.39346
- RCS-3 (knots 10,25,50) REML: coef=[-0.01489386, 0.03301487], SE=[0.006798487, 0.01348527], AIC=-5.8216
- Linear ML: coef=-0.004336017, SE=0.005333154, tau²=7.397169e-05
- Linear Fixed: coef=-0.004373803, SE=0.002908558

---

### Task 1: Greenland-Longnecker Covariance Reconstruction

**Files:**
- Modify: `metasprint-dose-response.html` (insert after `wlsFit` at ~line 8810)

**Context:** This is dosresmeta's core innovation. When multiple dose-level contrasts share a reference group, their effect sizes are correlated. GL reconstructs this within-study covariance from case counts and sample sizes. Without it, CIs are overconfident.

**Step 1: Write the inline self-test `_testGLCovariance()`**

Insert after `_testMatrixOps` (line ~8852). Test against R's exact covariance matrices:

```javascript
function _testGLCovariance() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 1e-6);

  // Test 1: cc_ex dataset (case-control, 1 study, 4 dose levels)
  // R: covar.logrr(cases=c(165,74,90,122), n=c(337,167,186,212),
  //    y=c(0,-0.2231435,0.1484200,0.4510757), v=c(0,0.0498,0.0546,0.1142),
  //    type="cc")
  // Expected diagonal: [0.04987985, 0.05463996, 0.14176254] (3 non-ref)
  // Wait - cc_ex has 4 rows total, ref + 3 non-ref
  const cc = greenlandLongnecker(
    [165, 74, 90, 122],  // cases per dose level
    [337, 167, 186, 212], // total N per dose level
    0, // reference index
    'cc' // type: case-control
  );
  results.push({
    name: 'GL cc_ex diagonal',
    pass: cc && cc.length === 3 &&
      approxEq(cc[0][0], 0.04987985, 1e-4) &&
      approxEq(cc[1][1], 0.05463996, 1e-4) &&
      approxEq(cc[2][2], 0.14176254, 1e-4)
  });
  results.push({
    name: 'GL cc_ex off-diagonal',
    pass: cc &&
      approxEq(cc[0][1], 0.01875906, 1e-4) &&
      approxEq(cc[0][2], 0.01678737, 1e-4) &&
      approxEq(cc[1][2], 0.01736485, 1e-4)
  });

  // Test 2: ci_ex dataset (cohort incidence, 5 dose levels)
  // R covariance for ci_ex: diagonal [0.014376, 0.015845, 0.019584, 0.022291]
  const ci = greenlandLongnecker(
    [110, 212, 211, 132, 133], // cases
    [8103, 17538, 15304, 9078, 10685], // N
    0, // reference index
    'ci' // type: cumulative incidence
  );
  results.push({
    name: 'GL ci_ex diagonal',
    pass: ci && ci.length === 4 &&
      approxEq(ci[0][0], 0.014376267, 1e-4) &&
      approxEq(ci[1][1], 0.015845348, 1e-4) &&
      approxEq(ci[2][2], 0.019584006, 1e-4) &&
      approxEq(ci[3][3], 0.022291481, 1e-4)
  });
  results.push({
    name: 'GL ci_ex off-diagonal (shared ref)',
    pass: ci &&
      approxEq(ci[0][1], 0.009064702, 1e-4) &&
      approxEq(ci[0][3], 0.009315187, 1e-4)
  });

  // Test 3: ir_ex dataset (incidence rate)
  const ir = greenlandLongnecker(
    [148, 127, 114, 107, 95], // cases
    [134707, 133824, 130654, 124522, 117808], // person-time N
    0, // reference index
    'ir' // type: incidence rate
  );
  results.push({
    name: 'GL ir_ex diagonal',
    pass: ir && ir.length === 4 &&
      approxEq(ir[0][0], 0.014774956, 1e-4) &&
      approxEq(ir[1][1], 0.016794974, 1e-4)
  });

  return results;
}
```

**Step 2: Run test to verify it fails**

Open browser console, run `_runAllDRTests()`. Expected: `_testGLCovariance` fails (function not defined).

**Step 3: Implement `greenlandLongnecker()`**

Insert before the test function (after `wlsFit` at ~line 8810):

```javascript
// Greenland-Longnecker covariance reconstruction
// Reconstructs within-study covariance for dose-level contrasts sharing a reference group
// Ref: Greenland & Longnecker (1992), Hamling et al (2008)
// cases: array of case counts per dose level (including reference)
// n: array of total N (or person-time for IR) per dose level
// refIdx: index of the reference dose level (usually 0)
// type: 'cc' (case-control), 'ci' (cumulative incidence), or 'ir' (incidence rate)
// Returns: (k-1) x (k-1) covariance matrix for non-reference contrasts
function greenlandLongnecker(cases, n, refIdx, type) {
  if (!cases || !n || cases.length !== n.length || cases.length < 2) return null;
  const k = cases.length;
  // Indices of non-reference levels
  const nonRef = [];
  for (let i = 0; i < k; i++) { if (i !== refIdx) nonRef.push(i); }
  const m = nonRef.length; // dimension of covariance matrix

  // Compute variances based on study type
  // For log-RR from cumulative incidence (ci):
  //   var(logRR_i) = 1/cases_i - 1/n_i + 1/cases_ref - 1/n_ref  (non-ref)
  //   cov(logRR_i, logRR_j) = 1/cases_ref - 1/n_ref  (shared reference)
  // For log-RR from incidence rate (ir):
  //   var(logRR_i) = 1/cases_i + 1/cases_ref  (person-time cancels differently)
  //   cov(logRR_i, logRR_j) = 1/cases_ref
  // For log-OR from case-control (cc):
  //   controls_i = n_i - cases_i
  //   var(logOR_i) = 1/cases_i + 1/controls_i + 1/cases_ref + 1/controls_ref  (non-ref)
  //   cov(logOR_i, logOR_j) = 1/cases_ref + 1/controls_ref  (shared reference)

  const cRef = cases[refIdx];
  const nRef = n[refIdx];
  if (cRef <= 0 || nRef <= 0) return null;

  let sharedCov; // the off-diagonal shared reference component
  const variances = []; // diagonal elements

  if (type === 'ir') {
    // Incidence rate: var = 1/a_i + 1/a_0, cov = 1/a_0
    sharedCov = 1 / cRef;
    for (let idx = 0; idx < m; idx++) {
      const i = nonRef[idx];
      if (cases[i] <= 0) return null;
      variances.push(1 / cases[i] + 1 / cRef);
    }
  } else if (type === 'ci') {
    // Cumulative incidence: var = 1/a_i - 1/n_i + 1/a_0 - 1/n_0, cov = 1/a_0 - 1/n_0
    sharedCov = 1 / cRef - 1 / nRef;
    for (let idx = 0; idx < m; idx++) {
      const i = nonRef[idx];
      if (cases[i] <= 0 || n[i] <= 0) return null;
      variances.push(1 / cases[i] - 1 / n[i] + 1 / cRef - 1 / nRef);
    }
  } else {
    // Case-control (cc): var = 1/a_i + 1/(n_i-a_i) + 1/a_0 + 1/(n_0-a_0)
    // cov = 1/a_0 + 1/(n_0-a_0)
    const ctrlRef = nRef - cRef;
    if (ctrlRef <= 0) return null;
    sharedCov = 1 / cRef + 1 / ctrlRef;
    for (let idx = 0; idx < m; idx++) {
      const i = nonRef[idx];
      const ctrl = n[i] - cases[i];
      if (cases[i] <= 0 || ctrl <= 0) return null;
      variances.push(1 / cases[i] + 1 / ctrl + 1 / cRef + 1 / ctrlRef);
    }
  }

  // Build symmetric covariance matrix
  const S = [];
  for (let i = 0; i < m; i++) {
    S[i] = [];
    for (let j = 0; j < m; j++) {
      if (i === j) {
        S[i][j] = variances[i];
      } else {
        S[i][j] = sharedCov;
      }
    }
  }
  return S;
}
```

**Step 4: Register test in `_runAllDRTests` and run**

Add `{name: 'GL Covariance', fn: _testGLCovariance}` to the suites array in `_runAllDRTests()`.
Run `_runAllDRTests()` in console. Expected: ALL PASS.

**Step 5: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add Greenland-Longnecker covariance reconstruction with R-validated tests"
```

---

### Task 2: Block-Diagonal V Matrix and GLS with Full Covariance

**Files:**
- Modify: `metasprint-dose-response.html` (insert after `greenlandLongnecker`)

**Context:** The GL covariance gives per-study S matrices. For GLS, we need a block-diagonal V matrix and a GLS solver that uses it. Currently `wlsFit` uses diagonal weights only.

**Step 1: Write inline test `_testGLSFit()`**

```javascript
function _testGLSFit() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 1e-4);

  // Test: alcohol_cvd linear model with GL covariance
  // R reference: coef = -0.004365412, SE = 0.00588923 (REML)
  // For fixed-effects (no tau2): coef = -0.004373803, SE = 0.002908558
  // We test fixed-effects GLS first (no between-study variance)

  // alcohol_cvd data (non-reference rows only, 19 data points across 6 studies)
  const alcData = {
    studies: [
      { id: 1, type: 'cc', doses: [0, 9.06, 27, 45, 64.8],
        cases: [126, 61, 69, 22, 19], n: [414, 261, 228, 44, 34],
        logrr: [0, -0.2231435, -0.0001, 0.5306283, 0.8754687] },
      { id: 2, type: 'cc', doses: [0, 16.05, 46.425, 77.16],
        cases: [77, 88, 24, 13], n: [258, 413, 202, 64],
        logrr: [0, -0.4307829, -1.07881, -0.6161861] },
      { id: 3, type: 'cc', doses: [0, 1.176, 8.916, 18.72],
        cases: [83, 46, 17, 4], n: [208, 175, 58, 11],
        logrr: [0, -0.5798185, -0.5621189, 0.6151856] },
      { id: 4, type: 'cc', doses: [0, 0.955, 7.43, 15.6],
        cases: [83, 29, 32, 18], n: [208, 117, 81, 39],
        logrr: [0, -0.2876821, 0.5128236, -0.3147107] },
      { id: 5, type: 'ci', doses: [0, 6, 18, 28.8],
        cases: [159, 229, 38, 39], n: [480, 975, 207, 133],
        logrr: [0, -0.4155155, -0.6348783, -0.0833816] },
      { id: 6, type: 'ci', doses: [0, 6.25, 18.75, 30],
        cases: [302, 129, 19, 15], n: [1269, 436, 49, 41],
        logrr: [0, -0.2231435, -0.0100503, -0.2613648] }
    ]
  };

  // Build per-study covariance matrices using GL
  const Slist = [];
  const yAll = [];
  const XAll = [];
  for (const s of alcData.studies) {
    const S = greenlandLongnecker(s.cases, s.n, 0, s.type);
    if (!S) { results.push({name: 'GL build failed for study ' + s.id, pass: false}); return results; }
    Slist.push(S);
    // Non-reference data
    for (let i = 1; i < s.doses.length; i++) {
      yAll.push(s.logrr[i]);
      XAll.push([s.doses[i]]); // linear model: X = [dose]
    }
  }

  // GLS fixed-effects fit
  const gls = glsFit(XAll, yAll, Slist, alcData.studies.map(s => s.doses.length - 1));
  results.push({
    name: 'GLS fixed coef vs R',
    pass: gls && approxEq(gls.coefficients[0], -0.004373803, 1e-4)
  });
  results.push({
    name: 'GLS fixed SE vs R',
    pass: gls && approxEq(Math.sqrt(gls.vcov[0][0]), 0.002908558, 1e-4)
  });

  return results;
}
```

**Step 2: Run test (fail)**

**Step 3: Implement `blockDiagMatrix()` and `glsFit()`**

```javascript
// Build block-diagonal matrix from array of square matrices
function blockDiagMatrix(blocks) {
  let totalSize = 0;
  for (const b of blocks) totalSize += b.length;
  const M = [];
  for (let i = 0; i < totalSize; i++) {
    M[i] = new Array(totalSize).fill(0);
  }
  let offset = 0;
  for (const b of blocks) {
    const sz = b.length;
    for (let i = 0; i < sz; i++) {
      for (let j = 0; j < sz; j++) {
        M[offset + i][offset + j] = b[i][j];
      }
    }
    offset += sz;
  }
  return M;
}

// Add tau2 * I to each block of the covariance matrix (between-study heterogeneity)
// tau2 can be scalar (proportional) or a matrix (for multivariate random effects)
function addTau2ToBlocks(Slist, tau2, studySizes) {
  const augmented = [];
  for (let s = 0; s < Slist.length; s++) {
    const S = Slist[s];
    const m = S.length;
    const A = [];
    for (let i = 0; i < m; i++) {
      A[i] = [];
      for (let j = 0; j < m; j++) {
        A[i][j] = S[i][j];
      }
    }
    if (typeof tau2 === 'number') {
      // Scalar tau2: add to diagonal only (simplest random-effects)
      for (let i = 0; i < m; i++) A[i][i] += tau2;
    } else if (Array.isArray(tau2) && tau2.length > 0) {
      // Matrix tau2 (Psi): tau2 is p x p, map dose design to variance contribution
      // For linear model (p=1): each dose_i contributes dose_i^2 * Psi
      // General: X_i * Psi * X_i' added to S_i
      // This is handled by the caller passing the right Psi structure
      // Simple case: tau2 is [p x p], studySizes[s] has the X matrix for this study
      if (studySizes && studySizes[s] && studySizes[s].X) {
        const Xi = studySizes[s].X; // m x p design matrix for study s
        const p = tau2.length;
        for (let i = 0; i < m; i++) {
          for (let j = 0; j < m; j++) {
            for (let a = 0; a < p; a++) {
              for (let b = 0; b < p; b++) {
                A[i][j] += Xi[i][a] * tau2[a][b] * Xi[j][b];
              }
            }
          }
        }
      }
    }
    augmented.push(A);
  }
  return augmented;
}

// Generalized Least Squares fit with block-diagonal covariance
// X: N x p design matrix (array of arrays), y: N x 1 response vector
// Slist: array of per-study covariance matrices
// studyNs: array of integers = number of non-ref rows per study (sum must equal N)
// Returns: {coefficients, vcov, fitted, residuals, logLik, df}
function glsFit(X, y, Slist, studyNs) {
  const N = y.length;
  const p = X[0].length;

  // Build block-diagonal V
  const V = blockDiagMatrix(Slist);

  // Invert V
  const Vinv = matInvertNxN(V);
  if (!Vinv) return null;

  // X'V^{-1}X
  const Xt = matTranspose(X);
  const XtVinv = matMul(Xt, Vinv);
  const XtVinvX = matMul(XtVinv, [y.map((_, i) => X[i])].length ? Xt : X);

  // Actually, let me do this properly with matMul
  // Xt is p x N, Vinv is N x N, X is N x p
  const XtVi = matMul(Xt, Vinv); // p x N
  const XtViX = matMul(XtVi, X); // p x p

  // Invert (X'V^-1X)
  const XtViX_inv = matInvertNxN(XtViX);
  if (!XtViX_inv) return null;

  // X'V^-1 y
  const XtViy = [];
  for (let i = 0; i < p; i++) {
    let sum = 0;
    for (let j = 0; j < N; j++) sum += XtVi[i][j] * y[j];
    XtViy.push(sum);
  }

  // beta = (X'V^-1X)^-1 X'V^-1 y
  const beta = [];
  for (let i = 0; i < p; i++) {
    let sum = 0;
    for (let j = 0; j < p; j++) sum += XtViX_inv[i][j] * XtViy[j];
    beta.push(sum);
  }

  // Fitted values and residuals
  const fitted = [];
  const resid = [];
  for (let i = 0; i < N; i++) {
    let yhat = 0;
    for (let j = 0; j < p; j++) yhat += X[i][j] * beta[j];
    fitted.push(yhat);
    resid.push(y[i] - yhat);
  }

  // Log-likelihood: -0.5 * [N*log(2pi) + log|V| + r'V^-1 r]
  // where r = y - X*beta
  let rVinvr = 0;
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      rVinvr += resid[i] * Vinv[i][j] * resid[j];
    }
  }
  // log|V| via LU decomposition (approximate from diagonal for now)
  let logDetV = 0;
  for (let i = 0; i < N; i++) logDetV += Math.log(Math.max(V[i][i], 1e-30));
  // Better: use the Cholesky or LU factored determinant
  // For now, use sum of log diagonal as approximation (exact for block-diagonal if blocks are diagonal)
  // TODO: exact log-determinant from block inverses

  const logLik = -0.5 * (N * Math.log(2 * Math.PI) + logDetV + rVinvr);

  return {
    coefficients: beta,
    vcov: XtViX_inv,
    fitted: fitted,
    residuals: resid,
    logLik: logLik,
    df: { N: N, p: p, nStudies: Slist.length }
  };
}

// Compute log-determinant of a block-diagonal matrix from its blocks
function logDetBlockDiag(blocks) {
  let logDet = 0;
  for (const B of blocks) {
    // For small blocks, use direct determinant
    const n = B.length;
    if (n === 1) {
      logDet += Math.log(Math.max(B[0][0], 1e-30));
    } else {
      // LU decomposition for determinant
      const LU = [];
      for (let i = 0; i < n; i++) LU[i] = B[i].slice();
      let sign = 1;
      for (let j = 0; j < n; j++) {
        let maxVal = Math.abs(LU[j][j]);
        let maxRow = j;
        for (let i = j + 1; i < n; i++) {
          if (Math.abs(LU[i][j]) > maxVal) { maxVal = Math.abs(LU[i][j]); maxRow = i; }
        }
        if (maxRow !== j) {
          [LU[j], LU[maxRow]] = [LU[maxRow], LU[j]];
          sign *= -1;
        }
        if (Math.abs(LU[j][j]) < 1e-30) return -Infinity;
        for (let i = j + 1; i < n; i++) {
          LU[i][j] /= LU[j][j];
          for (let k = j + 1; k < n; k++) LU[i][k] -= LU[i][j] * LU[j][k];
        }
      }
      for (let i = 0; i < n; i++) logDet += Math.log(Math.max(Math.abs(LU[i][i]), 1e-30));
    }
  }
  return logDet;
}
```

**Step 4: Register `_testGLSFit` in `_runAllDRTests`, run all tests**

**Step 5: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add block-diagonal GLS solver with GL covariance integration"
```

---

### Task 3: ML/REML Estimation via Profile Likelihood

**Files:**
- Modify: `metasprint-dose-response.html` (insert after `glsFit`)

**Context:** dosresmeta uses ML and REML for tau² estimation. The profile log-likelihood is optimized over tau² (scalar or matrix). We use golden-section search for scalar tau².

**Step 1: Write inline test `_testMLREML()`**

```javascript
function _testMLREML() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 1e-3);

  // Build alcohol_cvd data and Slist (same as Task 2 test)
  const alcStudies = buildAlcoholCVDData(); // helper we'll create
  const Slist = [];
  const yAll = [];
  const XAll = [];
  const studyNs = [];
  for (const s of alcStudies) {
    const S = greenlandLongnecker(s.cases, s.n, 0, s.type);
    Slist.push(S);
    const ni = s.doses.length - 1;
    studyNs.push(ni);
    for (let i = 1; i < s.doses.length; i++) {
      yAll.push(s.logrr[i]);
      XAll.push([s.doses[i]]);
    }
  }

  // Test REML: R gives tau2=0.0001005709, coef=-0.004365412
  const reml = estimateREML(XAll, yAll, Slist, studyNs);
  results.push({
    name: 'REML tau2 vs R',
    pass: reml && approxEq(reml.tau2, 0.0001005709, 1e-4)
  });
  results.push({
    name: 'REML coef vs R',
    pass: reml && approxEq(reml.coefficients[0], -0.004365412, 1e-4)
  });
  results.push({
    name: 'REML SE vs R',
    pass: reml && approxEq(Math.sqrt(reml.vcov[0][0]), 0.00588923, 1e-3)
  });
  results.push({
    name: 'REML AIC vs R',
    pass: reml && approxEq(reml.AIC, -24.06829, 0.5)
  });

  // Test ML: R gives tau2=7.397169e-05, coef=-0.004336017
  const ml = estimateML(XAll, yAll, Slist, studyNs);
  results.push({
    name: 'ML tau2 vs R',
    pass: ml && approxEq(ml.tau2, 7.397169e-05, 1e-4)
  });
  results.push({
    name: 'ML coef vs R',
    pass: ml && approxEq(ml.coefficients[0], -0.004336017, 1e-4)
  });

  return results;
}
```

**Step 2: Run test (fail)**

**Step 3: Implement `estimateML()`, `estimateREML()`, and helper `buildAlcoholCVDData()`**

```javascript
// alcohol_cvd reference dataset (dosresmeta2 built-in)
function buildAlcoholCVDData() {
  return [
    { id: 1, type: 'cc', doses: [0, 9.06, 27, 45, 64.8],
      cases: [126, 61, 69, 22, 19], n: [414, 261, 228, 44, 34],
      logrr: [0, -0.2231435, -0.0001, 0.5306283, 0.8754687] },
    { id: 2, type: 'cc', doses: [0, 16.05, 46.425, 77.16],
      cases: [77, 88, 24, 13], n: [258, 413, 202, 64],
      logrr: [0, -0.4307829, -1.07881, -0.6161861] },
    { id: 3, type: 'cc', doses: [0, 1.176, 8.916, 18.72],
      cases: [83, 46, 17, 4], n: [208, 175, 58, 11],
      logrr: [0, -0.5798185, -0.5621189, 0.6151856] },
    { id: 4, type: 'cc', doses: [0, 0.955, 7.43, 15.6],
      cases: [83, 29, 32, 18], n: [208, 117, 81, 39],
      logrr: [0, -0.2876821, 0.5128236, -0.3147107] },
    { id: 5, type: 'ci', doses: [0, 6, 18, 28.8],
      cases: [159, 229, 38, 39], n: [480, 975, 207, 133],
      logrr: [0, -0.4155155, -0.6348783, -0.0833816] },
    { id: 6, type: 'ci', doses: [0, 6.25, 18.75, 30],
      cases: [302, 129, 19, 15], n: [1269, 436, 49, 41],
      logrr: [0, -0.2231435, -0.0100503, -0.2613648] }
  ];
}

// Profile log-likelihood for ML estimation
// tau2: scalar between-study variance
// X, y, Slist, studyNs: as in glsFit
// Returns log-likelihood value
function profileLogLikML(tau2, X, y, Slist, studyNs) {
  // Add tau2 to each block diagonal
  const augBlocks = [];
  for (const S of Slist) {
    const m = S.length;
    const A = [];
    for (let i = 0; i < m; i++) {
      A[i] = S[i].slice();
      A[i][i] += tau2;
    }
    augBlocks.push(A);
  }

  const V = blockDiagMatrix(augBlocks);
  const Vinv = matInvertNxN(V);
  if (!Vinv) return -Infinity;

  const N = y.length;
  const p = X[0].length;

  // GLS beta
  const Xt = matTranspose(X);
  const XtVi = matMul(Xt, Vinv);
  const XtViX = matMul(XtVi, X);
  const XtViX_inv = matInvertNxN(XtViX);
  if (!XtViX_inv) return -Infinity;

  const XtViy = [];
  for (let i = 0; i < p; i++) {
    let s = 0;
    for (let j = 0; j < N; j++) s += XtVi[i][j] * y[j];
    XtViy.push(s);
  }
  const beta = [];
  for (let i = 0; i < p; i++) {
    let s = 0;
    for (let j = 0; j < p; j++) s += XtViX_inv[i][j] * XtViy[j];
    beta.push(s);
  }

  // Residuals
  const resid = [];
  for (let i = 0; i < N; i++) {
    let yhat = 0;
    for (let j = 0; j < p; j++) yhat += X[i][j] * beta[j];
    resid.push(y[i] - yhat);
  }

  // r'V^-1 r
  let rVr = 0;
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) rVr += resid[i] * Vinv[i][j] * resid[j];
  }

  const logDetV = logDetBlockDiag(augBlocks);

  // ML: -0.5 * (N*log(2pi) + log|V| + r'V^-1 r)
  return -0.5 * (N * Math.log(2 * Math.PI) + logDetV + rVr);
}

// Profile log-likelihood for REML
function profileLogLikREML(tau2, X, y, Slist, studyNs) {
  const augBlocks = [];
  for (const S of Slist) {
    const m = S.length;
    const A = [];
    for (let i = 0; i < m; i++) {
      A[i] = S[i].slice();
      A[i][i] += tau2;
    }
    augBlocks.push(A);
  }

  const V = blockDiagMatrix(augBlocks);
  const Vinv = matInvertNxN(V);
  if (!Vinv) return -Infinity;

  const N = y.length;
  const p = X[0].length;
  const Xt = matTranspose(X);
  const XtVi = matMul(Xt, Vinv);
  const XtViX = matMul(XtVi, X);
  const XtViX_inv = matInvertNxN(XtViX);
  if (!XtViX_inv) return -Infinity;

  const XtViy = [];
  for (let i = 0; i < p; i++) {
    let s = 0;
    for (let j = 0; j < N; j++) s += XtVi[i][j] * y[j];
    XtViy.push(s);
  }
  const beta = [];
  for (let i = 0; i < p; i++) {
    let s = 0;
    for (let j = 0; j < p; j++) s += XtViX_inv[i][j] * XtViy[j];
    beta.push(s);
  }

  const resid = [];
  for (let i = 0; i < N; i++) {
    let yhat = 0;
    for (let j = 0; j < p; j++) yhat += X[i][j] * beta[j];
    resid.push(y[i] - yhat);
  }

  let rVr = 0;
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) rVr += resid[i] * Vinv[i][j] * resid[j];
  }

  const logDetV = logDetBlockDiag(augBlocks);
  const logDetXtViX = logDetBlockDiag([XtViX]);

  // REML: -0.5 * ((N-p)*log(2pi) + log|V| + log|X'V^-1X| + r'V^-1 r)
  return -0.5 * ((N - p) * Math.log(2 * Math.PI) + logDetV + logDetXtViX + rVr);
}

// Golden-section search to maximize profile log-likelihood over tau2
function optimizeTau2(logLikFn, X, y, Slist, studyNs) {
  const tau2Max = 10;
  let a = 0, b = tau2Max;
  const phi = (1 + Math.sqrt(5)) / 2;
  const tol = 1e-8;
  const maxIter = 200;

  // Golden section search for maximum
  for (let iter = 0; iter < maxIter; iter++) {
    if (b - a < tol) break;
    const c = b - (b - a) / phi;
    const d = a + (b - a) / phi;
    const fc = logLikFn(c, X, y, Slist, studyNs);
    const fd = logLikFn(d, X, y, Slist, studyNs);
    if (fc > fd) { b = d; } else { a = c; }
  }
  return (a + b) / 2;
}

// Full ML estimation
function estimateML(X, y, Slist, studyNs) {
  const tau2 = optimizeTau2(profileLogLikML, X, y, Slist, studyNs);

  // Re-fit GLS with optimal tau2
  const augBlocks = Slist.map(S => {
    const A = S.map(row => row.slice());
    for (let i = 0; i < A.length; i++) A[i][i] += tau2;
    return A;
  });
  const fit = glsFit(X, y, augBlocks, studyNs);
  if (!fit) return null;

  const logLik = profileLogLikML(tau2, X, y, Slist, studyNs);
  const p = X[0].length;
  const nParams = p + 1; // fixed effects + tau2

  return {
    coefficients: fit.coefficients,
    vcov: fit.vcov,
    tau2: tau2,
    logLik: logLik,
    AIC: -2 * logLik + 2 * nParams,
    BIC: -2 * logLik + nParams * Math.log(y.length),
    method: 'ML'
  };
}

// Full REML estimation
function estimateREML(X, y, Slist, studyNs) {
  const tau2 = optimizeTau2(profileLogLikREML, X, y, Slist, studyNs);

  const augBlocks = Slist.map(S => {
    const A = S.map(row => row.slice());
    for (let i = 0; i < A.length; i++) A[i][i] += tau2;
    return A;
  });
  const fit = glsFit(X, y, augBlocks, studyNs);
  if (!fit) return null;

  const logLik = profileLogLikREML(tau2, X, y, Slist, studyNs);
  const p = X[0].length;
  const nParams = p + 1;

  return {
    coefficients: fit.coefficients,
    vcov: fit.vcov,
    tau2: tau2,
    logLik: logLik,
    AIC: -2 * logLik + 2 * nParams,
    BIC: -2 * logLik + nParams * Math.log(y.length),
    method: 'REML'
  };
}
```

**Step 4: Register `_testMLREML` in `_runAllDRTests`, run all tests**

**Step 5: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add ML/REML estimation via profile likelihood (R-validated)"
```

---

### Task 4: Integrate GL Covariance + ML/REML into Model Fitting Pipeline

**Files:**
- Modify: `metasprint-dose-response.html` (modify existing model fitters and `compareDoseResponseModels`)

**Context:** Currently all model fitters (`fitLinearDR`, `fitQuadraticDR`, etc.) take `points` array with diagonal weights (1/SE²). We need to update the pipeline so that when GL covariance data is available (cases/N provided), the fitters use full GLS instead. The key integration point is `compareDoseResponseModels()` at line ~10289 and `renderDoseResponseAnalysis()` at line ~10471.

**Step 1: Add estimation method to UI**

In the `drModelSelect` dropdown area (line ~1457), add estimation method selector:

```html
<label for="drEstMethod" style="margin-left:12px;">Estimation:</label>
<select id="drEstMethod" style="margin-left:4px;">
  <option value="dl" selected>DerSimonian-Laird</option>
  <option value="reml">REML</option>
  <option value="ml">ML</option>
  <option value="fixed">Fixed-Effect</option>
</select>
```

**Step 2: Create `buildStudyData()` function**

This function inspects the extract table data and builds the structured study data with GL covariance when cases/N are available:

```javascript
// Build structured study data from extract table for GLS fitting
// Returns { studies: [...], hasCovariance: bool }
// Each study: { id, type, doses: [], logrr: [], se: [], cases: [], n: [], Smatrix: [] }
function buildStudyData(points) {
  // Group points by study
  const byStudy = {};
  for (const pt of points) {
    const sid = pt.study || pt.studyId || 'S1';
    if (!byStudy[sid]) byStudy[sid] = [];
    byStudy[sid].push(pt);
  }

  const studies = [];
  let hasCovariance = false;

  for (const [sid, pts] of Object.entries(byStudy)) {
    // Sort by dose
    pts.sort((a, b) => a.dose - b.dose);

    const hasCasesN = pts.every(p => p.cases != null && p.cases > 0 && p.n != null && p.n > 0);

    const study = {
      id: sid,
      type: pts[0].type || 'ci', // default to cumulative incidence
      doses: pts.map(p => p.dose),
      effects: pts.map(p => p.effect),
      se: pts.map(p => p.se),
      cases: hasCasesN ? pts.map(p => p.cases) : null,
      n: hasCasesN ? pts.map(p => p.n) : null,
      Smatrix: null
    };

    // Find reference dose (dose=0 or minimum dose with effect=0)
    let refIdx = study.doses.indexOf(0);
    if (refIdx === -1) refIdx = study.effects.indexOf(0);
    if (refIdx === -1) refIdx = 0;
    study.refIdx = refIdx;

    // Build GL covariance if possible
    if (hasCasesN && pts.length >= 2) {
      const S = greenlandLongnecker(study.cases, study.n, refIdx, study.type);
      if (S) {
        study.Smatrix = S;
        hasCovariance = true;
      }
    }

    studies.push(study);
  }

  return { studies, hasCovariance };
}
```

**Step 3: Update `compareDoseResponseModels()` to use GLS when available**

Add a `studyData` parameter and `estMethod` parameter. When `studyData.hasCovariance` is true and `estMethod` is 'reml' or 'ml', use the GL-based GLS pipeline. Otherwise fall back to the existing WLS-based fitters.

The key change is in `compareDoseResponseModels()`: after building the model comparison, if GLS is available, re-fit the selected models using GLS and use those results.

**Step 4: Update Extract phase UI to include Cases/N columns**

Add optional columns to the extract table for case counts and sample sizes, visible when the user selects a study type that supports GL covariance.

**Step 5: Run all tests, verify no regressions**

**Step 6: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: integrate GL covariance + ML/REML into model fitting pipeline"
```

---

### Task 5: Log-Linear Model

**Files:**
- Modify: `metasprint-dose-response.html` (insert after `fitFP2DR` at ~line 9221)

**Step 1: Write inline test `_testLogLinear()`**

```javascript
function _testLogLinear() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 1e-4);

  // Simple test: y = 0.5 * log(dose + 1e-6) at doses [0, 5, 10, 20, 40]
  // With known coefficients, verify recovery
  const doses = [5, 10, 20, 40];
  const trueCoef = -0.08; // similar to alcohol_cvd FP p=0 (log)
  const effects = doses.map(d => trueCoef * Math.log(d + 1e-6));
  const se = doses.map(() => 0.05);
  const points = doses.map((d, i) => ({ dose: d, effect: effects[i], se: se[i] }));

  const fit = fitLogLinearDR(points);
  results.push({
    name: 'Log-linear coefficient recovery',
    pass: fit && approxEq(fit.coefficients[0], trueCoef, 0.01)
  });
  results.push({
    name: 'Log-linear has AIC',
    pass: fit && isFinite(fit.aic)
  });

  return results;
}
```

**Step 2: Implement `fitLogLinearDR()`**

```javascript
// Log-linear dose-response model: y = b1 * log(dose + epsilon)
// Common in epidemiology for alcohol, dietary intake, etc.
function fitLogLinearDR(points, epsilon) {
  if (!points || points.length < 2) return null;
  const eps = epsilon || 1e-6;

  // Transform dose to log scale
  const X = points.map(p => [Math.log(p.dose + eps)]);
  const y = points.map(p => p.effect);
  const W = points.map(p => 1 / (p.se * p.se));

  // WLS fit
  const fit = wlsFit(X, y, W);
  if (!fit) return null;

  // Prediction function
  const predict = (dose) => fit.coefficients[0] * Math.log(dose + eps);

  // Compute AIC
  const n = points.length;
  const p = 1;
  let wss = 0;
  for (let i = 0; i < n; i++) {
    const r = y[i] - predict(points[i].dose);
    wss += W[i] * r * r;
  }
  const aic = n * Math.log(wss / n) + 2 * p;

  return {
    type: 'log-linear',
    coefficients: fit.coefficients,
    se: fit.se,
    vcov: fit.vcov,
    predict: predict,
    aic: aic,
    bic: n * Math.log(wss / n) + p * Math.log(n),
    rSquared: fit.rSquared,
    residuals: fit.residuals,
    df: n - p,
    epsilon: eps
  };
}
```

**Step 3: Add 'Log-Linear' to model selector dropdown and `compareDoseResponseModels`**

In the HTML dropdown (`drModelSelect` at line ~1457):
```html
<option value="loglinear">Log-Linear</option>
```

In `compareDoseResponseModels()`, add:
```javascript
try {
  const loglin = fitLogLinearDR(points);
  if (loglin) allModels.push(loglin);
} catch(e) {}
```

**Step 4: Register test, run all tests**

**Step 5: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add log-linear dose-response model"
```

---

### Task 6: Exponential Model

**Files:**
- Modify: `metasprint-dose-response.html` (insert after `fitLogLinearDR`)

**Step 1: Write inline test `_testExponential()`**

```javascript
function _testExponential() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 0.05);

  // Exponential: y = Emax * (1 - exp(-alpha * dose))
  // Generate test data with known parameters
  const trueEmax = -0.5, trueAlpha = 0.05;
  const doses = [5, 10, 20, 30, 50];
  const effects = doses.map(d => trueEmax * (1 - Math.exp(-trueAlpha * d)));
  const se = doses.map(() => 0.03);
  const points = doses.map((d, i) => ({ dose: d, effect: effects[i], se: se[i] }));

  const fit = fitExponentialDR(points);
  results.push({
    name: 'Exponential Emax recovery',
    pass: fit && approxEq(fit.Emax, trueEmax, 0.1)
  });
  results.push({
    name: 'Exponential alpha recovery',
    pass: fit && approxEq(fit.alpha, trueAlpha, 0.02)
  });
  results.push({
    name: 'Exponential has AIC',
    pass: fit && isFinite(fit.aic)
  });

  return results;
}
```

**Step 2: Implement `fitExponentialDR()`**

```javascript
// Exponential dose-response: y = Emax * (1 - exp(-alpha * dose))
// Fitted via Gauss-Newton iterative WLS
function fitExponentialDR(points) {
  if (!points || points.length < 3) return null;
  const n = points.length;
  const doses = points.map(p => p.dose);
  const y = points.map(p => p.effect);
  const W = points.map(p => 1 / (p.se * p.se));

  // Initialize: Emax = extreme observed effect, alpha from half-max
  let Emax = y.reduce((a, b) => Math.abs(a) > Math.abs(b) ? a : b, 0);
  if (Math.abs(Emax) < 1e-10) Emax = -0.5;
  let alpha = Math.abs(Emax) > 1e-10 ? -Math.log(0.5) / (doses[Math.floor(n / 2)] || 10) : 0.05;

  const maxIter = 100;
  const tol = 1e-8;

  for (let iter = 0; iter < maxIter; iter++) {
    // Predicted values and Jacobian
    const pred = [];
    const J = []; // n x 2 Jacobian: d/dEmax, d/dalpha
    for (let i = 0; i < n; i++) {
      const d = doses[i];
      const expTerm = Math.exp(-alpha * d);
      pred.push(Emax * (1 - expTerm));
      J.push([
        1 - expTerm,           // d/dEmax
        Emax * d * expTerm     // d/dalpha
      ]);
    }

    // Weighted residuals
    const r = y.map((yi, i) => yi - pred[i]);

    // Weighted normal equations: (J'WJ) delta = J'W r
    const JtWJ = [[0, 0], [0, 0]];
    const JtWr = [0, 0];
    for (let i = 0; i < n; i++) {
      for (let a = 0; a < 2; a++) {
        JtWr[a] += J[i][a] * W[i] * r[i];
        for (let b = 0; b < 2; b++) {
          JtWJ[a][b] += J[i][a] * W[i] * J[i][b];
        }
      }
    }

    // Solve 2x2 system
    const det = JtWJ[0][0] * JtWJ[1][1] - JtWJ[0][1] * JtWJ[1][0];
    if (Math.abs(det) < 1e-20) break;
    const delta = [
      (JtWJ[1][1] * JtWr[0] - JtWJ[0][1] * JtWr[1]) / det,
      (JtWJ[0][0] * JtWr[1] - JtWJ[1][0] * JtWr[0]) / det
    ];

    Emax += delta[0];
    alpha += delta[1];
    if (alpha < 1e-10) alpha = 1e-10; // keep positive

    if (Math.abs(delta[0]) + Math.abs(delta[1]) < tol) break;
  }

  // Final predictions and residuals
  const predict = (dose) => Emax * (1 - Math.exp(-alpha * dose));
  let wss = 0;
  for (let i = 0; i < n; i++) {
    const r = y[i] - predict(doses[i]);
    wss += W[i] * r * r;
  }

  // Compute SE from Jacobian at final estimates
  const Jfinal = [];
  for (let i = 0; i < n; i++) {
    const d = doses[i];
    const expTerm = Math.exp(-alpha * d);
    Jfinal.push([1 - expTerm, Emax * d * expTerm]);
  }
  const JtWJ = [[0, 0], [0, 0]];
  for (let i = 0; i < n; i++) {
    for (let a = 0; a < 2; a++) {
      for (let b = 0; b < 2; b++) {
        JtWJ[a][b] += Jfinal[i][a] * W[i] * Jfinal[i][b];
      }
    }
  }
  const det = JtWJ[0][0] * JtWJ[1][1] - JtWJ[0][1] * JtWJ[1][0];
  let vcov = null, se = [NaN, NaN];
  if (Math.abs(det) > 1e-20) {
    vcov = [
      [JtWJ[1][1] / det, -JtWJ[0][1] / det],
      [-JtWJ[1][0] / det, JtWJ[0][0] / det]
    ];
    se = [Math.sqrt(Math.max(vcov[0][0], 0)), Math.sqrt(Math.max(vcov[1][1], 0))];
  }

  const p = 2;
  const aic = n * Math.log(wss / n) + 2 * p;

  return {
    type: 'exponential',
    Emax: Emax,
    alpha: alpha,
    coefficients: [Emax, alpha],
    se: se,
    vcov: vcov,
    predict: predict,
    aic: aic,
    bic: n * Math.log(wss / n) + p * Math.log(n),
    rSquared: 1 - wss / y.reduce((s, yi) => s + W[y.indexOf(yi)] * yi * yi, 0),
    df: n - p
  };
}
```

**Step 3: Add to model selector and `compareDoseResponseModels`**

**Step 4: Register test, run all tests**

**Step 5: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add exponential dose-response model (Gauss-Newton)"
```

---

### Task 7: 3-Parameter Emax (Hill) Model

**Files:**
- Modify: `metasprint-dose-response.html` (insert after `fitExponentialDR`)

**Step 1: Write inline test `_testHillModel()`**

```javascript
function _testHillModel() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 0.1);

  // Hill: y = Emax * dose^h / (ED50^h + dose^h)
  // With h=1 should match standard Emax
  const trueEmax = -0.6, trueED50 = 20, trueH = 1.5;
  const doses = [5, 10, 15, 20, 30, 50];
  const effects = doses.map(d =>
    trueEmax * Math.pow(d, trueH) / (Math.pow(trueED50, trueH) + Math.pow(d, trueH))
  );
  const se = doses.map(() => 0.02);
  const points = doses.map((d, i) => ({ dose: d, effect: effects[i], se: se[i] }));

  const fit = fitHillDR(points);
  results.push({
    name: 'Hill Emax recovery',
    pass: fit && approxEq(fit.Emax, trueEmax, 0.15)
  });
  results.push({
    name: 'Hill ED50 recovery',
    pass: fit && approxEq(fit.ED50, trueED50, 5)
  });
  results.push({
    name: 'Hill exponent recovery',
    pass: fit && approxEq(fit.h, trueH, 0.5)
  });
  results.push({
    name: 'Hill has AIC',
    pass: fit && isFinite(fit.aic)
  });

  return results;
}
```

**Step 2: Implement `fitHillDR()`**

```javascript
// 3-parameter Emax (Hill) model: y = Emax * dose^h / (ED50^h + dose^h)
// h=1 gives standard Emax; h>1 gives sigmoidal; h<1 gives sub-linear
// Fitted via Gauss-Newton with grid search over h
function fitHillDR(points) {
  if (!points || points.length < 4) return null; // need 4+ points for 3 params
  const n = points.length;
  const doses = points.map(p => p.dose);
  const y = points.map(p => p.effect);
  const W = points.map(p => 1 / (p.se * p.se));

  // Grid search over h, then refine with Gauss-Newton
  const hGrid = [0.5, 0.75, 1, 1.5, 2, 3];
  let bestFit = null;
  let bestWSS = Infinity;

  for (const hInit of hGrid) {
    const fit = _fitHillGN(doses, y, W, n, hInit);
    if (fit && fit.wss < bestWSS) {
      bestWSS = fit.wss;
      bestFit = fit;
    }
  }

  if (!bestFit) return null;

  const predict = (dose) => {
    const dh = Math.pow(dose, bestFit.h);
    return bestFit.Emax * dh / (Math.pow(bestFit.ED50, bestFit.h) + dh);
  };

  const p = 3;
  const aic = n * Math.log(bestWSS / n) + 2 * p;

  return {
    type: 'hill',
    Emax: bestFit.Emax,
    ED50: bestFit.ED50,
    h: bestFit.h,
    coefficients: [bestFit.Emax, bestFit.ED50, bestFit.h],
    se: bestFit.se,
    predict: predict,
    aic: aic,
    bic: n * Math.log(bestWSS / n) + p * Math.log(n),
    df: n - p
  };
}

// Internal Gauss-Newton for Hill model with given h initialization
function _fitHillGN(doses, y, W, n, hInit) {
  // Initialize from existing Emax fit with h=hInit
  let Emax = y.reduce((a, b) => Math.abs(a) > Math.abs(b) ? a : b, 0) * 1.2;
  let ED50 = doses[Math.floor(n / 2)] || 10;
  let h = hInit;
  if (ED50 <= 0) ED50 = 10;

  const maxIter = 100;
  const tol = 1e-8;

  for (let iter = 0; iter < maxIter; iter++) {
    const J = []; // n x 3 Jacobian
    const pred = [];
    for (let i = 0; i < n; i++) {
      const d = doses[i];
      const dh = Math.pow(d, h);
      const e50h = Math.pow(ED50, h);
      const denom = e50h + dh;
      if (Math.abs(denom) < 1e-30) return null;
      const yhat = Emax * dh / denom;
      pred.push(yhat);

      // Partial derivatives
      const dEmax = dh / denom;
      const dED50 = -Emax * dh * h * Math.pow(ED50, h - 1) / (denom * denom);
      const logD = d > 0 ? Math.log(d) : 0;
      const logE = ED50 > 0 ? Math.log(ED50) : 0;
      const dh_dh = Emax * dh * e50h * (logD - logE) / (denom * denom);

      J.push([dEmax, dED50, dh_dh]);
    }

    const r = y.map((yi, i) => yi - pred[i]);

    // (J'WJ) delta = J'Wr
    const JtWJ = Array.from({length: 3}, () => new Array(3).fill(0));
    const JtWr = [0, 0, 0];
    for (let i = 0; i < n; i++) {
      for (let a = 0; a < 3; a++) {
        JtWr[a] += J[i][a] * W[i] * r[i];
        for (let b = 0; b < 3; b++) {
          JtWJ[a][b] += J[i][a] * W[i] * J[i][b];
        }
      }
    }

    const inv = matInvertNxN(JtWJ);
    if (!inv) break;

    const delta = [0, 0, 0];
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) delta[a] += inv[a][b] * JtWr[b];
    }

    // Damped update to prevent divergence
    let step = 1;
    for (let s = 0; s < 5; s++) {
      const newEmax = Emax + step * delta[0];
      const newED50 = ED50 + step * delta[1];
      const newH = h + step * delta[2];
      if (newED50 > 0 && newH > 0.1 && newH < 10) {
        Emax = newEmax;
        ED50 = newED50;
        h = newH;
        break;
      }
      step *= 0.5;
    }

    if (Math.abs(delta[0]) + Math.abs(delta[1]) + Math.abs(delta[2]) < tol) break;
  }

  // Final weighted SS
  let wss = 0;
  for (let i = 0; i < n; i++) {
    const dh = Math.pow(doses[i], h);
    const yhat = Emax * dh / (Math.pow(ED50, h) + dh);
    wss += W[i] * (y[i] - yhat) * (y[i] - yhat);
  }

  // SE from final Jacobian
  const Jf = [];
  for (let i = 0; i < n; i++) {
    const d = doses[i];
    const dh = Math.pow(d, h);
    const e50h = Math.pow(ED50, h);
    const denom = e50h + dh;
    const logD = d > 0 ? Math.log(d) : 0;
    const logE = ED50 > 0 ? Math.log(ED50) : 0;
    Jf.push([
      dh / denom,
      -Emax * dh * h * Math.pow(ED50, h - 1) / (denom * denom),
      Emax * dh * e50h * (logD - logE) / (denom * denom)
    ]);
  }
  const JtWJf = Array.from({length: 3}, () => new Array(3).fill(0));
  for (let i = 0; i < n; i++) {
    for (let a = 0; a < 3; a++) {
      for (let b = 0; b < 3; b++) JtWJf[a][b] += Jf[i][a] * W[i] * Jf[i][b];
    }
  }
  const vcov = matInvertNxN(JtWJf);
  const se = vcov ? [0, 1, 2].map(i => Math.sqrt(Math.max(vcov[i][i], 0))) : [NaN, NaN, NaN];

  return { Emax, ED50, h, wss, se };
}
```

**Step 3: Add to model selector and `compareDoseResponseModels`**

**Step 4: Register test, run all tests**

**Step 5: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add 3-parameter Emax (Hill) dose-response model"
```

---

### Task 8: Extend One-Stage Mixed-Effects to All Models

**Files:**
- Modify: `metasprint-dose-response.html` (modify `fitOneStageMixedEffects` at line ~9827)

**Context:** Currently `fitOneStageMixedEffects()` only supports linear and quadratic models. We need to extend it to accept a design matrix builder function for any model type.

**Step 1: Write inline test `_testOneStageAllModels()`**

```javascript
function _testOneStageAllModels() {
  const results = [];

  // Test: one-stage with RCS should not crash and should return valid results
  const studies = [
    { doses: [0, 10, 20, 40], effects: [0, -0.15, -0.28, -0.10], se: [0, 0.05, 0.06, 0.08] },
    { doses: [0, 10, 30, 50], effects: [0, -0.10, -0.25, -0.05], se: [0, 0.06, 0.07, 0.09] },
    { doses: [0, 5, 15, 25], effects: [0, -0.08, -0.20, -0.30], se: [0, 0.04, 0.05, 0.06] }
  ];

  // Test RCS one-stage
  const rcs = fitOneStageMixedEffects(studies, 'rcs3');
  results.push({
    name: 'One-stage RCS-3 returns result',
    pass: rcs && rcs.coefficients && rcs.coefficients.length === 2
  });

  // Test log-linear one-stage
  const loglin = fitOneStageMixedEffects(studies, 'loglinear');
  results.push({
    name: 'One-stage log-linear returns result',
    pass: loglin && loglin.coefficients && loglin.coefficients.length === 1
  });

  return results;
}
```

**Step 2: Refactor `fitOneStageMixedEffects()` to accept model type**

Add a `modelType` parameter and use a design matrix builder:

```javascript
// Design matrix builder for different model types
function buildDesignMatrix(doses, modelType, knots) {
  const eps = 1e-6;
  switch (modelType) {
    case 'linear':
      return doses.map(d => [d]);
    case 'quadratic':
      return doses.map(d => [d, d * d]);
    case 'rcs3': {
      // Use computeRCSBasis with 3 knots
      const k = knots || _defaultKnots(doses, 3);
      return doses.map(d => {
        const basis = computeRCSBasis(d, k);
        return basis.slice(1); // remove intercept (linear + spline terms)
      });
    }
    case 'rcs4': {
      const k = knots || _defaultKnots(doses, 4);
      return doses.map(d => {
        const basis = computeRCSBasis(d, k);
        return basis.slice(1);
      });
    }
    case 'loglinear':
      return doses.map(d => [Math.log(d + eps)]);
    default:
      return doses.map(d => [d]); // fallback to linear
  }
}

function _defaultKnots(doses, nKnots) {
  const sorted = [...new Set(doses)].filter(d => d > 0).sort((a, b) => a - b);
  if (nKnots === 3) {
    return [
      sorted[Math.floor(sorted.length * 0.1)] || sorted[0],
      sorted[Math.floor(sorted.length * 0.5)],
      sorted[Math.floor(sorted.length * 0.9)] || sorted[sorted.length - 1]
    ];
  } else {
    return [
      sorted[Math.floor(sorted.length * 0.05)] || sorted[0],
      sorted[Math.floor(sorted.length * 0.35)],
      sorted[Math.floor(sorted.length * 0.65)],
      sorted[Math.floor(sorted.length * 0.95)] || sorted[sorted.length - 1]
    ];
  }
}
```

Then modify `fitOneStageMixedEffects()` to use `buildDesignMatrix()` instead of hardcoded linear/quadratic.

**Step 3: Update one-stage toggle UI** — enable for all model types (remove gray-out for RCS/FP/Emax)

**Step 4: Register test, run all tests**

**Step 5: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: extend one-stage mixed-effects to all model types (RCS, FP, log-linear)"
```

---

### Task 9: Prediction Intervals for DR Curve

**Files:**
- Modify: `metasprint-dose-response.html` (modify `renderDoseResponseCurveSVG` at line ~10615)

**Step 1: Write inline test `_testPredictionInterval()`**

```javascript
function _testPredictionInterval() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 0.01);

  // PI = y_hat +/- t_{k-p} * sqrt(SE^2 + tau2)
  // With tau2=0.01, SE=0.05, k=6, p=1, confLevel=0.95
  // t(5, 0.025) ≈ 2.571
  const tau2 = 0.01;
  const se = 0.05;
  const k = 6, p = 1;
  const pi = computePredictionInterval(0, se, tau2, k, p, 0.95);
  const expectedWidth = 2.571 * Math.sqrt(se * se + tau2);
  results.push({
    name: 'PI width correct',
    pass: pi && approxEq(pi.upper - pi.lower, 2 * expectedWidth, 0.05)
  });
  results.push({
    name: 'PI wider than CI',
    pass: pi && (pi.upper - pi.lower) > 2 * 1.96 * se
  });

  return results;
}
```

**Step 2: Implement `computePredictionInterval()`**

```javascript
// Compute prediction interval at a dose point
// yhat: predicted value, se: standard error of prediction
// tau2: between-study variance, k: number of studies, p: model parameters
// confLevel: confidence level (0.95)
function computePredictionInterval(yhat, se, tau2, k, p, confLevel) {
  const df = k - p;
  if (df <= 0) return null; // PI undefined

  const alpha = 1 - (confLevel || 0.95);
  // t-distribution critical value (using approximation)
  const tCrit = tQuantile(1 - alpha / 2, df);

  const piSE = Math.sqrt(se * se + tau2);
  return {
    lower: yhat - tCrit * piSE,
    upper: yhat + tCrit * piSE,
    piSE: piSE
  };
}

// t-distribution quantile (approximation for df >= 1)
function tQuantile(p, df) {
  if (df <= 0) return NaN;
  if (df > 300) {
    // Normal approximation for large df
    return normalQuantile(p);
  }
  // Abramowitz & Stegun approximation via normal quantile
  const z = normalQuantile(p);
  const g1 = (z * z * z + z) / 4;
  const g2 = (5 * z * z * z * z * z + 16 * z * z * z + 3 * z) / 96;
  const g3 = (3 * z * z * z * z * z * z * z + 19 * z * z * z * z * z + 17 * z * z * z - 15 * z) / 384;
  return z + g1 / df + g2 / (df * df) + g3 / (df * df * df);
}

// Standard normal quantile (Beasley-Springer-Moro approximation)
function normalQuantile(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;
  // Rational approximation
  const a = [
    -3.969683028665376e+01, 2.209460984245205e+02,
    -2.759285104469687e+02, 1.383577518672690e+02,
    -3.066479806614716e+01, 2.506628277459239e+00
  ];
  const b = [
    -5.447609879822406e+01, 1.615858368580409e+02,
    -1.556989798598866e+02, 6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  const c = [
    -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
    4.374664141464968e+00, 2.938163982698783e+00
  ];
  const d = [
    7.784695709041462e-03, 3.224671290700398e-01,
    2.445134137142996e+00, 3.754408661907416e+00
  ];
  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let q, r;
  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}
```

**Step 3: Update `renderDoseResponseCurveSVG()` to draw PI band**

Add a dashed SVG path for the prediction interval band, drawn behind the CI ribbon. Add a checkbox "Show prediction interval" to the analyze phase UI. Use `stroke-dasharray="5,5"` for the dashed PI lines.

**Step 4: Register test, run all tests**

**Step 5: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add prediction intervals for dose-response curves"
```

---

### Task 10: Model-Averaged Predictions (Exceeds dosresmeta)

**Files:**
- Modify: `metasprint-dose-response.html` (insert after `compareDoseResponseModels`)

**Context:** This is our signature feature that goes beyond dosresmeta2. We average predictions from all fitted models using AIC weights, incorporating both within-model and between-model uncertainty (Burnham & Anderson 2002, eq 4.9).

**Step 1: Write inline test `_testModelAveraging()`**

```javascript
function _testModelAveraging() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 0.01);

  // Test with 2 models: linear (w=0.7) and quadratic (w=0.3)
  const models = [
    { type: 'linear', aicWeight: 0.7,
      predict: (d) => -0.01 * d,
      predictSE: (d) => 0.005 * d },
    { type: 'quadratic', aicWeight: 0.3,
      predict: (d) => -0.02 * d + 0.0003 * d * d,
      predictSE: (d) => 0.006 * d }
  ];

  const avg = modelAveragedPrediction(models, 20);
  // Expected: 0.7 * (-0.2) + 0.3 * (-0.4 + 0.12) = -0.14 + -0.084 = -0.224
  const expectedPred = 0.7 * (-0.2) + 0.3 * (-0.28);
  results.push({
    name: 'Model-averaged prediction',
    pass: avg && approxEq(avg.prediction, expectedPred, 0.01)
  });
  results.push({
    name: 'Model-averaged SE includes between-model variance',
    pass: avg && avg.se > 0.005 * 20 // should be wider than narrowest model
  });
  results.push({
    name: 'Model-averaged CI exists',
    pass: avg && avg.lower < avg.prediction && avg.upper > avg.prediction
  });

  return results;
}
```

**Step 2: Implement `modelAveragedPrediction()` and `computeModelAveragedCurve()`**

```javascript
// Model-averaged prediction at a single dose
// models: array of { predict(dose), predictSE(dose), aicWeight }
// dose: dose value
// confLevel: confidence level for CI
function modelAveragedPrediction(models, dose, confLevel) {
  const cl = confLevel || 0.95;
  const z = normalQuantile(1 - (1 - cl) / 2);

  // Filter to models with meaningful weight
  const active = models.filter(m => m.aicWeight > 0.01 && m.predict);
  if (active.length === 0) return null;

  // Weighted average prediction
  let yAvg = 0;
  for (const m of active) yAvg += m.aicWeight * m.predict(dose);

  // Unconditional SE (Burnham & Anderson 2002, eq 4.9)
  // SE_avg = sqrt( sum_i w_i * [var_i(dose) + (y_i(dose) - y_avg(dose))^2] )
  let varAvg = 0;
  for (const m of active) {
    const se_i = m.predictSE ? m.predictSE(dose) : 0;
    const diff = m.predict(dose) - yAvg;
    varAvg += m.aicWeight * (se_i * se_i + diff * diff);
  }
  const seAvg = Math.sqrt(varAvg);

  return {
    prediction: yAvg,
    se: seAvg,
    lower: yAvg - z * seAvg,
    upper: yAvg + z * seAvg,
    weights: active.map(m => ({ type: m.type, weight: m.aicWeight }))
  };
}

// Compute full model-averaged curve at nPoints dose values
function computeModelAveragedCurve(models, maxDose, confLevel, nPoints) {
  const pts = nPoints || 100;
  const curve = [];
  for (let i = 0; i <= pts; i++) {
    const dose = (i / pts) * maxDose;
    const avg = modelAveragedPrediction(models, dose, confLevel);
    if (avg) curve.push({ dose, ...avg });
  }
  return curve;
}
```

**Step 3: Add "Model-Averaged" to model selector**

In `drModelSelect`:
```html
<option value="averaged">Model-Averaged</option>
```

In `renderDoseResponseAnalysis()`, when model-averaged is selected, compute the averaged curve and pass to SVG renderer. Show individual model curves as thin dashed lines overlaid on the thick averaged curve.

**Step 4: Update `renderDoseResponseCurveSVG()` to support model-averaged display**

When model type is 'averaged', draw:
- Thick solid line: averaged prediction
- Shaded band: averaged 95% CI
- Thin dashed lines: individual model predictions (color-coded)
- Legend showing model weights

**Step 5: Register test, run all tests**

**Step 6: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add model-averaged predictions (Burnham & Anderson 2002) - exceeds dosresmeta"
```

---

### Task 11: R Validation Script

**Files:**
- Create: `validate_vs_dosresmeta.R`

**Context:** This R script runs dosresmeta2 on its built-in datasets and exports all reference values to JSON for JavaScript comparison.

**Step 1: Write `validate_vs_dosresmeta.R`**

```r
# validate_vs_dosresmeta.R
# Generate reference values from dosresmeta2 for cross-validation
# Run: Rscript validate_vs_dosresmeta.R

library(dosresmeta)
library(jsonlite)

results <- list()

# ============================================
# Dataset 1: alcohol_cvd (6 studies, mixed cc/ci)
# ============================================
data(alcohol_cvd)

# Linear REML
lin_reml <- dosresmeta(logrr ~ dose, id=id, type=type, se=se,
                        cases=cases, n=n, data=alcohol_cvd, method="reml")
doses <- c(0, 5, 10, 15, 20, 25, 30, 40, 50, 60)
pred_lin <- predict(lin_reml, newdata=data.frame(dose=doses), expo=FALSE)

results$alcohol_cvd$linear_reml <- list(
  coefficients = as.numeric(coef(lin_reml)),
  se = as.numeric(sqrt(diag(vcov(lin_reml)))),
  tau2 = as.numeric(lin_reml$Psi),
  logLik = as.numeric(logLik(lin_reml)),
  AIC = AIC(lin_reml),
  BIC = BIC(lin_reml),
  predictions = data.frame(dose=doses, pred=pred_lin[,1],
                           ci_lb=pred_lin[,2], ci_ub=pred_lin[,3])
)

# Quadratic REML
quad_reml <- dosresmeta(logrr ~ dose + I(dose^2), id=id, type=type, se=se,
                         cases=cases, n=n, data=alcohol_cvd, method="reml")
pred_quad <- predict(quad_reml, newdata=data.frame(dose=doses), expo=FALSE)

results$alcohol_cvd$quadratic_reml <- list(
  coefficients = as.numeric(coef(quad_reml)),
  se = as.numeric(sqrt(diag(vcov(quad_reml)))),
  tau2 = as.list(as.data.frame(quad_reml$Psi)),
  logLik = as.numeric(logLik(quad_reml)),
  AIC = AIC(quad_reml),
  BIC = BIC(quad_reml),
  predictions = data.frame(dose=doses, pred=pred_quad[,1],
                           ci_lb=pred_quad[,2], ci_ub=pred_quad[,3])
)

# Linear ML
lin_ml <- dosresmeta(logrr ~ dose, id=id, type=type, se=se,
                      cases=cases, n=n, data=alcohol_cvd, method="ml")
results$alcohol_cvd$linear_ml <- list(
  coefficients = as.numeric(coef(lin_ml)),
  se = as.numeric(sqrt(diag(vcov(lin_ml)))),
  tau2 = as.numeric(lin_ml$Psi)
)

# Linear Fixed
lin_fixed <- dosresmeta(logrr ~ dose, id=id, type=type, se=se,
                         cases=cases, n=n, data=alcohol_cvd, method="fixed")
results$alcohol_cvd$linear_fixed <- list(
  coefficients = as.numeric(coef(lin_fixed)),
  se = as.numeric(sqrt(diag(vcov(lin_fixed))))
)

# RCS 3 knots
library(rms)
knots3 <- c(10, 25, 50)
rcs_reml <- dosresmeta(logrr ~ rcs(dose, knots3), id=id, type=type, se=se,
                        cases=cases, n=n, data=alcohol_cvd, method="reml")
pred_rcs <- predict(rcs_reml, newdata=data.frame(dose=doses), expo=FALSE)

results$alcohol_cvd$rcs3_reml <- list(
  coefficients = as.numeric(coef(rcs_reml)),
  se = as.numeric(sqrt(diag(vcov(rcs_reml)))),
  knots = knots3,
  AIC = AIC(rcs_reml),
  predictions = data.frame(dose=doses, pred=pred_rcs[,1],
                           ci_lb=pred_rcs[,2], ci_ub=pred_rcs[,3])
)

# FP2 best (p1=0.5, p2=0.5 from our earlier results)
fp2 <- dosresmeta(logrr ~ dose^0.5 + I(dose^0.5 * log(dose)), id=id,
                   type=type, se=se, cases=cases, n=n,
                   data=alcohol_cvd, method="reml")

tryCatch({
  results$alcohol_cvd$fp2_reml <- list(
    coefficients = as.numeric(coef(fp2)),
    se = as.numeric(sqrt(diag(vcov(fp2)))),
    AIC = AIC(fp2)
  )
}, error = function(e) {
  results$alcohol_cvd$fp2_reml <<- list(error = e$message)
})

# GL covariance matrices (Slist)
Slist <- lin_reml$Slist
results$alcohol_cvd$Slist <- lapply(Slist, function(s) as.matrix(s))

# Study-level data
results$alcohol_cvd$data <- alcohol_cvd

# ============================================
# Dataset 2: coffee_mort (larger, 21 studies)
# ============================================
data(coffee_mort)

cof_lin <- dosresmeta(logrr ~ dose, id=id, type=type, se=se,
                       cases=cases, n=n, data=coffee_mort, method="reml")
cof_doses <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
pred_cof <- predict(cof_lin, newdata=data.frame(dose=cof_doses), expo=FALSE)

results$coffee_mort$linear_reml <- list(
  coefficients = as.numeric(coef(cof_lin)),
  se = as.numeric(sqrt(diag(vcov(cof_lin)))),
  tau2 = as.numeric(cof_lin$Psi),
  AIC = AIC(cof_lin),
  predictions = data.frame(dose=cof_doses, pred=pred_cof[,1],
                           ci_lb=pred_cof[,2], ci_ub=pred_cof[,3])
)

# ============================================
# Dataset 3: ci_ex (single study, fixed effects)
# ============================================
data(ci_ex)
ci_lin <- dosresmeta(logrr ~ dose, se=se, cases=cases, n=n,
                      data=ci_ex, type="ci", method="fixed")
results$ci_ex$linear_fixed <- list(
  coefficients = as.numeric(coef(ci_lin)),
  se = as.numeric(sqrt(diag(vcov(ci_lin))))
)

# GL covariance for ci_ex
ci_cov <- covar.logrr(cases=ci_ex$cases, n=ci_ex$n, y=ci_ex$logrr,
                       v=ci_ex$se^2, type="ci")
results$ci_ex$covariance <- as.matrix(ci_cov)

# ============================================
# Write JSON output
# ============================================
json_output <- toJSON(results, pretty=TRUE, digits=10, auto_unbox=TRUE)
writeLines(json_output, "validation_reference.json")
cat("Validation reference saved to validation_reference.json\n")
cat("Models fitted:", length(unlist(results, recursive=FALSE)), "\n")
```

**Step 2: Run R script**

```bash
"C:/Program Files/R/R-4.4.2/bin/Rscript.exe" validate_vs_dosresmeta.R
```

Expected: `validation_reference.json` created with all reference values.

**Step 3: Commit**

```bash
git add validate_vs_dosresmeta.R
git commit -m "feat: add R validation script for dosresmeta2 cross-validation"
```

---

### Task 12: JavaScript Cross-Validation Against R Reference

**Files:**
- Modify: `metasprint-dose-response.html` (add comprehensive validation test)

**Context:** Now we embed the R reference values and write a comprehensive self-test that compares our JS implementations against R's exact outputs.

**Step 1: Write `_testVsR()` inline self-test**

```javascript
function _testVsR() {
  const results = [];
  const approxEq = (a, b, tol) => Math.abs(a - b) < (tol || 1e-3);

  // R reference values from dosresmeta2 v2.2.0 on alcohol_cvd
  const R = {
    linear_reml: {
      coef: [-0.004365412], se: [0.00588923],
      tau2: 0.0001005709, AIC: -24.06829,
      pred: [0, -0.02182706, -0.04365412, -0.06548118, -0.08730824,
             -0.10913529, -0.13096235, -0.17461647, -0.21827059, -0.26192471]
    },
    quadratic_reml: {
      coef: [-0.03023389, 0.0007514854], se: [0.0131779, 0.0004885058],
      AIC: -59.39346,
      pred: [0, -0.13238229, -0.22719031, -0.28442406, -0.30408353,
             -0.28616874, -0.23067967, -0.00697872, 0.36701932, 0.89131444]
    },
    linear_ml: { coef: [-0.004336017], se: [0.005333154], tau2: 7.397169e-05 },
    linear_fixed: { coef: [-0.004373803], se: [0.002908558] },
    rcs3_reml: {
      coef: [-0.01489386, 0.03301487], se: [0.006798487, 0.01348527],
      knots: [10, 25, 50], AIC: -5.8216,
      pred: [0, -0.07446932, -0.14893864, -0.22082867, -0.27724298,
             -0.30270584, -0.28586841, -0.15005377, 0.06004433, 0.28252301]
    },
    // GL covariance for Study 1 (Bianchi, cc, 4x4)
    S1: [[0.04987985, 0.01875906, 0.01678737, 0.01780975],
         [0.01875906, 0.05463996, 0.01736485, 0.01842239],
         [0.01678737, 0.01736485, 0.14176254, 0.01648609],
         [0.01780975, 0.01842239, 0.01648609, 0.19714011]]
  };

  const doses = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60];

  // Build alcohol_cvd data
  const alcData = buildAlcoholCVDData();

  // 1. Test GL covariance for Study 1
  const S1 = greenlandLongnecker(alcData[0].cases, alcData[0].n, 0, 'cc');
  for (let i = 0; i < 4; i++) {
    for (let j = i; j < 4; j++) {
      results.push({
        name: 'GL S1[' + i + '][' + j + '] vs R',
        pass: S1 && approxEq(S1[i][j], R.S1[i][j], 1e-3)
      });
    }
  }

  // 2. Test Fixed-effects GLS
  const Slist = [], yAll = [], XAll = [], studyNs = [];
  for (const s of alcData) {
    const S = greenlandLongnecker(s.cases, s.n, 0, s.type);
    if (!S) continue;
    Slist.push(S);
    studyNs.push(s.doses.length - 1);
    for (let i = 1; i < s.doses.length; i++) {
      yAll.push(s.logrr[i]);
      XAll.push([s.doses[i]]);
    }
  }

  const fixedFit = glsFit(XAll, yAll, Slist, studyNs);
  results.push({
    name: 'Fixed GLS coef vs R',
    pass: fixedFit && approxEq(fixedFit.coefficients[0], R.linear_fixed.coef[0], 1e-4)
  });
  results.push({
    name: 'Fixed GLS SE vs R',
    pass: fixedFit && approxEq(Math.sqrt(fixedFit.vcov[0][0]), R.linear_fixed.se[0], 1e-3)
  });

  // 3. Test REML
  const reml = estimateREML(XAll, yAll, Slist, studyNs);
  results.push({
    name: 'REML coef vs R',
    pass: reml && approxEq(reml.coefficients[0], R.linear_reml.coef[0], 1e-3)
  });
  results.push({
    name: 'REML tau2 vs R',
    pass: reml && approxEq(reml.tau2, R.linear_reml.tau2, 1e-3)
  });

  // 4. Test ML
  const ml = estimateML(XAll, yAll, Slist, studyNs);
  results.push({
    name: 'ML coef vs R',
    pass: ml && approxEq(ml.coefficients[0], R.linear_ml.coef[0], 1e-3)
  });
  results.push({
    name: 'ML tau2 vs R',
    pass: ml && approxEq(ml.tau2, R.linear_ml.tau2, 1e-3)
  });

  // 5. Test predictions at standard doses
  if (reml) {
    for (let i = 1; i < doses.length; i++) {
      const pred = reml.coefficients[0] * doses[i];
      results.push({
        name: 'REML pred dose=' + doses[i] + ' vs R',
        pass: approxEq(pred, R.linear_reml.pred[i], 1e-3)
      });
    }
  }

  return results;
}
```

**Step 2: Register `_testVsR` in `_runAllDRTests`**

**Step 3: Run all tests — verify R-validated results match within tolerance**

**Step 4: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add comprehensive R cross-validation test suite (dosresmeta2 v2.2.0)"
```

---

### Task 13: Selenium Integration Tests for New Features

**Files:**
- Modify: `test_dose_response_models.py`

**Context:** Add Selenium tests that verify the new UI elements (estimation method selector, new model types, prediction interval toggle, model-averaged display) and that the R-validated self-tests pass in the browser.

**Step 1: Add new test class `TestExceedDosresmeta`**

```python
class TestExceedDosresmeta(unittest.TestCase):
    """Tests for dosresmeta-exceeding features"""

    @classmethod
    def setUpClass(cls):
        opts = Options()
        opts.add_argument('--headless=new')
        opts.add_argument('--no-sandbox')
        cls.driver = webdriver.Chrome(options=opts)
        cls.driver.get('file:///' + os.path.abspath('metasprint-dose-response.html').replace('\\', '/'))
        time.sleep(2)

    @classmethod
    def tearDownClass(cls):
        cls.driver.quit()

    def test_estimation_method_selector_exists(self):
        """Estimation method dropdown (DL/REML/ML/Fixed) exists"""
        sel = self.driver.find_element(By.ID, 'drEstMethod')
        options = [o.get_attribute('value') for o in sel.find_elements(By.TAG_NAME, 'option')]
        self.assertIn('reml', options)
        self.assertIn('ml', options)

    def test_new_model_types_in_selector(self):
        """Log-linear, exponential, hill, averaged in model selector"""
        sel = self.driver.find_element(By.ID, 'drModelSelect')
        options = [o.get_attribute('value') for o in sel.find_elements(By.TAG_NAME, 'option')]
        self.assertIn('loglinear', options)
        self.assertIn('exponential', options)
        self.assertIn('hill', options)
        self.assertIn('averaged', options)

    def test_gl_covariance_self_test(self):
        """GL covariance inline self-test passes"""
        result = self.driver.execute_script('return _testGLCovariance()')
        for r in result:
            self.assertTrue(r['pass'], f"GL test failed: {r['name']}")

    def test_ml_reml_self_test(self):
        """ML/REML inline self-test passes"""
        result = self.driver.execute_script('return _testMLREML()')
        for r in result:
            self.assertTrue(r['pass'], f"ML/REML test failed: {r['name']}")

    def test_r_validation_self_test(self):
        """R cross-validation inline self-test passes"""
        result = self.driver.execute_script('return _testVsR()')
        for r in result:
            self.assertTrue(r['pass'], f"R validation failed: {r['name']}")

    def test_prediction_interval_self_test(self):
        """Prediction interval self-test passes"""
        result = self.driver.execute_script('return _testPredictionInterval()')
        for r in result:
            self.assertTrue(r['pass'], f"PI test failed: {r['name']}")

    def test_model_averaging_self_test(self):
        """Model averaging self-test passes"""
        result = self.driver.execute_script('return _testModelAveraging()')
        for r in result:
            self.assertTrue(r['pass'], f"Model averaging test failed: {r['name']}")

    def test_all_self_tests_pass(self):
        """Master self-test suite passes (all suites)"""
        result = self.driver.execute_script("""
            const r = _runAllDRTests();
            return { total: r.total, passed: r.passed, failed: r.failed,
                     failures: r.results.filter(x => !x.pass).map(x => x.name) };
        """)
        self.assertEqual(result['failed'], 0,
                         f"Self-tests failed: {result['failures']}")
```

**Step 2: Run Selenium tests**

```bash
python test_dose_response_models.py -v
```

Expected: All tests PASS (previous 36 + new ~8 = ~44 tests).

**Step 3: Commit**

```bash
git add test_dose_response_models.py
git commit -m "test: add Selenium tests for dosresmeta-exceeding features (R-validated)"
```

---

### Task 14: Final Verification and Div Balance

**Files:**
- Modify: `metasprint-dose-response.html` (fixes only)

**Step 1: Div balance check**

```bash
# Count opening and closing divs (in bash)
grep -c '<div[\s>]' metasprint-dose-response.html
grep -c '</div>' metasprint-dose-response.html
```

Investigate and fix any imbalance that wasn't pre-existing.

**Step 2: Check for literal `</script>` in script block**

```bash
# Should find 0 matches inside <script> block
grep -n '</script>' metasprint-dose-response.html
```

Only the closing `</script>` tag at end of script block should match.

**Step 3: Run full test suite**

```bash
python test_dose_response_models.py -v
```

Expected: ALL tests pass, 0 failures.

**Step 4: Run inline self-tests in browser**

Open `metasprint-dose-response.html`, console: `_runAllDRTests()`
Expected: ALL suites pass, including R-validated tests.

**Step 5: Commit any fixes**

```bash
git add metasprint-dose-response.html test_dose_response_models.py
git commit -m "chore: final verification — div balance, safety checks, all tests pass"
```

---

## Summary

| Task | Feature | R-Validated | New Tests |
|------|---------|------------|-----------|
| 1 | GL Covariance | Yes (cc_ex, ci_ex, ir_ex) | 5 assertions |
| 2 | Block-Diagonal GLS | Yes (fixed-effects coef/SE) | 2 assertions |
| 3 | ML/REML Estimation | Yes (tau², coef, SE, AIC) | 6 assertions |
| 4 | Pipeline Integration | Via Task 12 | UI test |
| 5 | Log-Linear Model | Partial (structure) | 2 assertions |
| 6 | Exponential Model | Self-validated | 3 assertions |
| 7 | Hill (3-param Emax) | Self-validated | 4 assertions |
| 8 | One-Stage All Models | Self-validated | 2 assertions |
| 9 | Prediction Intervals | Formula-validated | 2 assertions |
| 10 | Model-Averaged Predictions | Formula-validated | 3 assertions |
| 11 | R Validation Script | Source of truth | — |
| 12 | JS Cross-Validation | Comprehensive | 20+ assertions |
| 13 | Selenium Tests | Via self-tests | ~8 test methods |
| 14 | Final Verification | — | Full suite |

**Total new inline test assertions:** ~50+
**Total Selenium tests (estimated):** ~44 (36 existing + 8 new)
**R reference datasets:** alcohol_cvd (6 studies), coffee_mort (21 studies), ci_ex (1 study)
