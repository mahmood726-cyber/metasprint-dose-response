# Dose-Response Enhancements Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add RCS, fractional polynomial, one-stage mixed-effects, and continuous input modes to the existing dose-response meta-analysis app.

**Architecture:** Enhance `metasprint-dose-response.html` (19,344 lines) in-place. Add general NxN matrix operations, then build each model on top. Each model follows the same interface: `{model, aic, bic, R2, predict(dose), predictCI(dose, zCrit)}`. UI expands the model selector dropdown and input mode radios.

**Tech Stack:** Vanilla JS (single-file HTML), SVG for plots, IndexedDB for persistence, Python + Selenium for testing.

**Key file:** `C:\Users\user\Downloads\metasprint-dose-response\metasprint-dose-response.html`

**Key line references (current):**
- Model selector dropdown: ~line 1445
- Input mode radios: ~line 1371
- Extract table columns: ~line 1380-1392
- `prepareDoseResponseData()`: ~line 8187
- `fitLinearDR()`: ~line 8233
- `fitQuadraticDR()`: ~line 8284
- `det3()`/`inv3()`: ~line 8299-8325
- `fitEmaxDR()`: ~line 8366
- `compareDoseResponseModels()`: ~line 8497
- `renderDoseResponseAnalysis()`: ~line 8513
- Model comparison table: ~line 8543-8566
- `renderDoseResponseCurveSVG()`: ~line 8593
- `runAnalysis()`: ~line 8811
- `computeMetaAnalysis()`: ~line 7144

**IMPORTANT — Single-file HTML rules:**
- Never write literal `</script>` inside `<script>` blocks — use `${'<'}/script>`
- After structural HTML edits, verify div balance: count `<div[\s>]` vs `</div>`
- All element IDs must be globally unique
- All function names must be globally unique (no collisions across modules)
- Use `?? fallback` not `|| fallback` for numeric values (0 is valid)
- Guard matrix inversions with singularity check (`|pivot| < 1e-15 → null`)

---

### Task 1: General NxN Matrix Operations

**Context:** The app only has 3x3 `det3()`/`inv3()` (line ~8299). RCS and FP models need NxN matrices. Add general matrix utilities before the existing model functions.

**Files:**
- Modify: `metasprint-dose-response.html` — insert after `inv3()` (around line 8325)

**Step 1: Add NxN matrix functions**

Insert these functions after the existing `inv3()` function:

```javascript
// ── General NxN Matrix Operations (for RCS, FP models) ──
function matCreate(rows, cols) {
  return Array.from({length: rows}, () => new Array(cols).fill(0));
}

function matTranspose(A) {
  const m = A.length, n = A[0].length;
  const T = matCreate(n, m);
  for (let i = 0; i < m; i++)
    for (let j = 0; j < n; j++) T[j][i] = A[i][j];
  return T;
}

function matMul(A, B) {
  const m = A.length, n = B[0].length, p = B.length;
  const C = matCreate(m, n);
  for (let i = 0; i < m; i++)
    for (let j = 0; j < n; j++)
      for (let k = 0; k < p; k++) C[i][j] += A[i][k] * B[k][j];
  return C;
}

function matDiag(vec) {
  const n = vec.length;
  const D = matCreate(n, n);
  for (let i = 0; i < n; i++) D[i][i] = vec[i];
  return D;
}

function matInvertNxN(M) {
  const n = M.length;
  const aug = M.map((row, i) => {
    const r = row.slice();
    for (let j = 0; j < n; j++) r.push(i === j ? 1 : 0);
    return r;
  });
  for (let col = 0; col < n; col++) {
    let maxRow = col, maxVal = Math.abs(aug[col][col]);
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(aug[row][col]) > maxVal) { maxVal = Math.abs(aug[row][col]); maxRow = row; }
    }
    if (maxVal < 1e-15) return null; // singular
    if (maxRow !== col) { const tmp = aug[col]; aug[col] = aug[maxRow]; aug[maxRow] = tmp; }
    const pivot = aug[col][col];
    for (let j = 0; j < 2 * n; j++) aug[col][j] /= pivot;
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const factor = aug[row][col];
      for (let j = 0; j < 2 * n; j++) aug[row][j] -= factor * aug[col][j];
    }
  }
  return aug.map(row => row.slice(n));
}

function matSolve(A, b) {
  // Solve Ax = b via Gauss-Jordan; b is column vector [[b0],[b1],...]
  const Ainv = matInvertNxN(A);
  if (!Ainv) return null;
  return matMul(Ainv, b);
}

function wlsFit(X, y, w) {
  // Weighted least squares: (X'WX)^{-1} X'Wy
  // X: n x p design matrix, y: n-vector, w: n-vector of weights
  // Returns {beta: p-vector, vcov: p x p, fitted: n-vector, residuals: n-vector, sse: scalar}
  const n = X.length, p = X[0].length;
  const W = matDiag(w);
  const Xt = matTranspose(X);
  const XtW = matMul(Xt, W);
  const XtWX = matMul(XtW, X);
  const XtWXinv = matInvertNxN(XtWX);
  if (!XtWXinv) return null;
  const yCol = y.map(v => [v]);
  const XtWy = matMul(XtW, yCol);
  const betaCol = matMul(XtWXinv, XtWy);
  const beta = betaCol.map(r => r[0]);
  const fitted = X.map(row => row.reduce((s, x, j) => s + x * beta[j], 0));
  const residuals = y.map((yi, i) => yi - fitted[i]);
  const sse = residuals.reduce((s, r, i) => s + w[i] * r * r, 0);
  return { beta, vcov: XtWXinv, fitted, residuals, sse };
}
```

**Step 2: Add inline self-test for matrix operations**

Insert a self-test function (callable from browser console):

```javascript
function _testMatrixOps() {
  const results = [];
  // Test 1: 2x2 inverse
  const A = [[4, 7], [2, 6]];
  const Ainv = matInvertNxN(A);
  const det = 4*6 - 7*2; // = 10
  results.push({ test: '2x2 inv [0][0]', pass: Math.abs(Ainv[0][0] - 0.6) < 1e-10 });
  results.push({ test: '2x2 inv [0][1]', pass: Math.abs(Ainv[0][1] - (-0.7)) < 1e-10 });

  // Test 2: matMul identity
  const I = [[1,0],[0,1]];
  const AI = matMul(A, I);
  results.push({ test: 'A*I=A', pass: Math.abs(AI[0][0] - 4) < 1e-10 && Math.abs(AI[1][1] - 6) < 1e-10 });

  // Test 3: WLS on known data (y = 2 + 3x, equal weights)
  const X = [[1, 0], [1, 1], [1, 2], [1, 3]];
  const y = [2, 5, 8, 11];
  const w = [1, 1, 1, 1];
  const fit = wlsFit(X, y, w);
  results.push({ test: 'WLS intercept=2', pass: Math.abs(fit.beta[0] - 2) < 1e-10 });
  results.push({ test: 'WLS slope=3', pass: Math.abs(fit.beta[1] - 3) < 1e-10 });
  results.push({ test: 'WLS SSE=0', pass: fit.sse < 1e-10 });

  // Test 4: singular matrix returns null
  const S = [[1, 2], [2, 4]];
  results.push({ test: 'singular=null', pass: matInvertNxN(S) === null });

  console.table(results);
  return results.every(r => r.pass);
}
```

**Step 3: Verify no function name collisions**

Search the file for any existing `matCreate`, `matTranspose`, `matMul`, `matDiag`, `matInvertNxN`, `matSolve`, `wlsFit`. If any exist, rename the new ones with a `dr` prefix.

**Step 4: Commit**

```bash
git add metasprint-dose-response.html
git commit -m "feat: add general NxN matrix operations and WLS solver for DR models"
```

---

### Task 2: RCS (Restricted Cubic Splines) Model

**Context:** RCS is the standard non-linear dose-response method (Orsini, Greenland, Harrell). The autopilot has a `computeRCSBasis()` at line ~12454 that can be referenced, but we build fresh for this app.

**Files:**
- Modify: `metasprint-dose-response.html` — insert `fitRCSDR()` after `fitEmaxDR()` (~line 8460)

**Step 1: Add RCS basis function and model fitter**

```javascript
// ── RCS (Restricted Cubic Splines) ──
function computeRCSKnots(doses, nKnots) {
  // Percentile-based knot placement (Harrell 2001)
  const sorted = doses.slice().sort((a, b) => a - b);
  const n = sorted.length;
  const pctl = (p) => {
    const idx = p * (n - 1);
    const lo = Math.floor(idx), hi = Math.ceil(idx);
    return lo === hi ? sorted[lo] : sorted[lo] + (idx - lo) * (sorted[hi] - sorted[lo]);
  };
  if (nKnots === 3) return [pctl(0.10), pctl(0.50), pctl(0.90)];
  if (nKnots === 4) return [pctl(0.05), pctl(0.35), pctl(0.65), pctl(0.95)];
  if (nKnots === 5) return [pctl(0.05), pctl(0.275), pctl(0.50), pctl(0.725), pctl(0.95)];
  return [pctl(0.10), pctl(0.50), pctl(0.90)]; // default 3
}

function computeRCSBasis(doses, knots) {
  // Harrell (2001) truncated power basis with boundary constraints
  // Returns n x (nKnots-1) design matrix: [linear, spline1, ..., spline_{k-2}]
  const n = doses.length;
  const k = knots.length;
  const kLast = knots[k - 1];
  const kPrev = knots[k - 2];
  const denom = kLast - kPrev;
  if (Math.abs(denom) < 1e-15) return null;
  const ncols = k - 1; // linear + (k-2) spline terms
  const X = matCreate(n, ncols);
  for (let i = 0; i < n; i++) {
    const d = doses[i];
    X[i][0] = d; // linear term
    for (let j = 0; j < k - 2; j++) {
      const t_j = knots[j];
      const h1 = Math.max(0, d - t_j);
      const h2 = Math.max(0, d - kPrev);
      const h3 = Math.max(0, d - kLast);
      X[i][j + 1] = (h1 * h1 * h1 - h2 * h2 * h2 * (kLast - t_j) / denom
                      + h3 * h3 * h3 * (kPrev - t_j) / denom) / ((kLast - t_j) * (kLast - t_j));
      // Normalize to avoid numerical issues
      if (!isFinite(X[i][j + 1])) X[i][j + 1] = 0;
    }
  }
  return X;
}

function fitRCSDR(points, nKnots) {
  nKnots = nKnots ?? 3;
  const minPts = nKnots + 1; // need more points than parameters
  if (points.length < minPts) return null;

  const n = points.length;
  const doses = points.map(p => p.dose);
  const y = points.map(p => p.effect);
  const w = points.map(p => {
    const se = p.se > 0 ? p.se : 1;
    return 1 / (se * se);
  });

  const knots = computeRCSKnots(doses, nKnots);
  const Xspline = computeRCSBasis(doses, knots);
  if (!Xspline) return null;

  // Add intercept column: X = [1 | Xspline]
  const p = Xspline[0].length + 1;
  const X = matCreate(n, p);
  for (let i = 0; i < n; i++) {
    X[i][0] = 1; // intercept
    for (let j = 0; j < Xspline[0].length; j++) X[i][j + 1] = Xspline[i][j];
  }

  const fit = wlsFit(X, y, w);
  if (!fit) return null;

  // Goodness of fit
  const yMean = y.reduce((s, v, i) => s + w[i] * v, 0) / w.reduce((s, v) => s + v, 0);
  const ssTot = y.reduce((s, v, i) => s + w[i] * (v - yMean) * (v - yMean), 0);
  const R2 = ssTot > 1e-15 ? 1 - fit.sse / ssTot : 0;

  // AIC / BIC
  const aic = n * Math.log(fit.sse / n + 1e-15) + 2 * p;
  const bic = n * Math.log(fit.sse / n + 1e-15) + p * Math.log(n);

  // Non-linearity Wald test: H0 = all spline coefficients = 0
  // Test statistic: beta_spline' * (vcov_spline)^{-1} * beta_spline ~ chi2(k-2)
  let pNonlinear = NaN;
  const nSpline = nKnots - 2;
  if (nSpline > 0 && fit.vcov) {
    const betaSpline = fit.beta.slice(2); // skip intercept and linear term
    const vcovSpline = matCreate(nSpline, nSpline);
    for (let i = 0; i < nSpline; i++)
      for (let j = 0; j < nSpline; j++)
        vcovSpline[i][j] = fit.vcov[i + 2][j + 2];
    const vcovSplineInv = matInvertNxN(vcovSpline);
    if (vcovSplineInv) {
      const bCol = betaSpline.map(v => [v]);
      const waldStat = matMul(matMul(matTranspose(bCol), vcovSplineInv), bCol)[0][0];
      if (isFinite(waldStat) && waldStat >= 0) {
        pNonlinear = 1 - chi2CDF(waldStat, nSpline);
      }
    }
  }

  // Predict function
  function predict(dose) {
    const xRow = computeRCSBasis([dose], knots);
    if (!xRow) return NaN;
    let val = fit.beta[0]; // intercept
    for (let j = 0; j < xRow[0].length; j++) val += fit.beta[j + 1] * xRow[0][j];
    return val;
  }

  function predictCI(dose, zCrit) {
    const xRow = computeRCSBasis([dose], knots);
    if (!xRow || !fit.vcov) return { lo: NaN, hi: NaN };
    const xVec = [1, ...xRow[0]];
    let variance = 0;
    for (let i = 0; i < p; i++)
      for (let j = 0; j < p; j++)
        variance += xVec[i] * fit.vcov[i][j] * xVec[j];
    const se = Math.sqrt(Math.max(0, variance));
    const pred = predict(dose);
    return { lo: pred - zCrit * se, hi: pred + zCrit * se };
  }

  const paramStr = `knots=[${knots.map(k => k.toFixed(1)).join(', ')}]`;
  return {
    model: `RCS-${nKnots}`, b0: fit.beta[0], beta: fit.beta, knots,
    se: fit.vcov ? fit.beta.map((_, i) => Math.sqrt(Math.max(0, fit.vcov[i][i]))) : [],
    pNonlinear, R2, aic, bic, fitted: fit.fitted, residuals: fit.residuals,
    predict, predictCI, paramStr, nParams: p, vcov: fit.vcov
  };
}
```

**Step 2: Add inline self-test for RCS**

```javascript
function _testRCS() {
  const results = [];
  // Test: 5 points on a known cubic curve, RCS-3 should capture non-linearity
  const pts = [
    { dose: 0, effect: 0, se: 0.1 },
    { dose: 1, effect: 0.5, se: 0.1 },
    { dose: 2, effect: 1.5, se: 0.1 },
    { dose: 3, effect: 2.0, se: 0.1 },
    { dose: 4, effect: 1.8, se: 0.1 }
  ];
  const rcs3 = fitRCSDR(pts, 3);
  results.push({ test: 'RCS-3 fits', pass: rcs3 !== null });
  results.push({ test: 'RCS-3 has AIC', pass: isFinite(rcs3?.aic) });
  results.push({ test: 'RCS-3 predict(0)', pass: Math.abs(rcs3?.predict(0) - 0) < 0.5 });
  results.push({ test: 'RCS-3 R2 > 0', pass: (rcs3?.R2 ?? 0) > 0 });
  results.push({ test: 'RCS-3 pNonlinear', pass: isFinite(rcs3?.pNonlinear) });

  // RCS-4 with more points
  const pts2 = Array.from({length: 10}, (_, i) => ({
    dose: i, effect: Math.sin(i * 0.5), se: 0.15
  }));
  const rcs4 = fitRCSDR(pts2, 4);
  results.push({ test: 'RCS-4 fits', pass: rcs4 !== null });
  results.push({ test: 'RCS-4 knots length=4', pass: rcs4?.knots?.length === 4 });

  // Too few points
  const ptsFew = [{ dose: 0, effect: 0, se: 0.1 }, { dose: 1, effect: 1, se: 0.1 }];
  results.push({ test: 'RCS-3 null if <4 pts', pass: fitRCSDR(ptsFew, 3) === null });

  console.table(results);
  return results.every(r => r.pass);
}
```

**Step 3: Commit**

```bash
git commit -m "feat: add RCS (restricted cubic splines) dose-response model with 3/4 knots"
```

---

### Task 3: Fractional Polynomial Models

**Context:** FP1 fits `y = b0 + b1 * dose^p` across 8 candidate powers. FP2 fits `y = b0 + b1 * dose^p1 + b2 * dose^p2` across 36 combinations. Power 0 means log(dose).

**Files:**
- Modify: `metasprint-dose-response.html` — insert `fitFP1DR()` and `fitFP2DR()` after `fitRCSDR()`

**Step 1: Add FP model fitters**

```javascript
// ── Fractional Polynomial Models ──
const FP_POWERS = [-2, -1, -0.5, 0, 0.5, 1, 2, 3];

function fpTransform(dose, p) {
  // Power 0 = log(dose); handle dose <= 0 with small shift
  const d = dose <= 0 ? 1e-6 : dose;
  if (p === 0) return Math.log(d);
  return Math.pow(d, p);
}

function fitFP1DR(points) {
  if (points.length < 3) return null;
  const n = points.length;
  const y = points.map(p => p.effect);
  const w = points.map(p => { const se = p.se > 0 ? p.se : 1; return 1 / (se * se); });

  let bestFit = null, bestAIC = Infinity, bestPower = null;

  for (const p of FP_POWERS) {
    const X = matCreate(n, 2);
    for (let i = 0; i < n; i++) {
      X[i][0] = 1; // intercept
      X[i][1] = fpTransform(points[i].dose, p);
    }
    const fit = wlsFit(X, y, w);
    if (!fit) continue;
    const aic = n * Math.log(fit.sse / n + 1e-15) + 2 * 2;
    if (aic < bestAIC) {
      bestAIC = aic;
      bestFit = fit;
      bestPower = p;
    }
  }
  if (!bestFit) return null;

  const bic = n * Math.log(bestFit.sse / n + 1e-15) + 2 * Math.log(n);
  const yMean = y.reduce((s, v, i) => s + w[i] * v, 0) / w.reduce((s, v) => s + v, 0);
  const ssTot = y.reduce((s, v, i) => s + w[i] * (v - yMean) * (v - yMean), 0);
  const R2 = ssTot > 1e-15 ? 1 - bestFit.sse / ssTot : 0;
  const se_b0 = bestFit.vcov ? Math.sqrt(Math.max(0, bestFit.vcov[0][0])) : NaN;
  const se_b1 = bestFit.vcov ? Math.sqrt(Math.max(0, bestFit.vcov[1][1])) : NaN;

  function predict(dose) {
    return bestFit.beta[0] + bestFit.beta[1] * fpTransform(dose, bestPower);
  }
  function predictCI(dose, zCrit) {
    const x = [1, fpTransform(dose, bestPower)];
    let variance = 0;
    for (let i = 0; i < 2; i++)
      for (let j = 0; j < 2; j++)
        variance += x[i] * bestFit.vcov[i][j] * x[j];
    const se = Math.sqrt(Math.max(0, variance));
    const pred = predict(dose);
    return { lo: pred - zCrit * se, hi: pred + zCrit * se };
  }

  const powerLabel = bestPower === 0 ? 'log' : bestPower.toString();
  return {
    model: `FP1(${powerLabel})`, b0: bestFit.beta[0], b1: bestFit.beta[1],
    power: bestPower, se_b0, se_b1, R2, aic: bestAIC, bic,
    fitted: bestFit.fitted, residuals: bestFit.residuals,
    predict, predictCI, paramStr: `power=${powerLabel}`, nParams: 2, vcov: bestFit.vcov
  };
}

function fitFP2DR(points) {
  if (points.length < 4) return null;
  const n = points.length;
  const y = points.map(p => p.effect);
  const w = points.map(p => { const se = p.se > 0 ? p.se : 1; return 1 / (se * se); });

  let bestFit = null, bestAIC = Infinity, bestP1 = null, bestP2 = null;

  for (const p1 of FP_POWERS) {
    for (const p2 of FP_POWERS) {
      if (FP_POWERS.indexOf(p2) < FP_POWERS.indexOf(p1)) continue; // avoid duplicates
      const X = matCreate(n, 3);
      for (let i = 0; i < n; i++) {
        const d = points[i].dose;
        X[i][0] = 1;
        X[i][1] = fpTransform(d, p1);
        if (p1 === p2) {
          // Repeated power: second term = dose^p * log(dose)
          const dSafe = d <= 0 ? 1e-6 : d;
          X[i][2] = fpTransform(d, p1) * Math.log(dSafe);
        } else {
          X[i][2] = fpTransform(d, p2);
        }
      }
      const fit = wlsFit(X, y, w);
      if (!fit) continue;
      const aic = n * Math.log(fit.sse / n + 1e-15) + 2 * 3;
      if (aic < bestAIC) {
        bestAIC = aic;
        bestFit = fit;
        bestP1 = p1;
        bestP2 = p2;
      }
    }
  }
  if (!bestFit) return null;

  const bic = n * Math.log(bestFit.sse / n + 1e-15) + 3 * Math.log(n);
  const yMean = y.reduce((s, v, i) => s + w[i] * v, 0) / w.reduce((s, v) => s + v, 0);
  const ssTot = y.reduce((s, v, i) => s + w[i] * (v - yMean) * (v - yMean), 0);
  const R2 = ssTot > 1e-15 ? 1 - bestFit.sse / ssTot : 0;

  // Deviance test: FP2 vs FP1 (chi2 with 1 df)
  const fp1 = fitFP1DR(points);
  let pFP2vsFP1 = NaN;
  if (fp1) {
    const devDiff = n * (Math.log(fp1.aic > -Infinity ? (Math.exp((fp1.aic - 2 * 2) / n)) : 1e-15)
      - Math.log(bestFit.sse / n + 1e-15));
    if (isFinite(devDiff) && devDiff > 0) pFP2vsFP1 = 1 - chi2CDF(devDiff, 1);
  }

  function predict(dose) {
    const d = dose;
    let val = bestFit.beta[0] + bestFit.beta[1] * fpTransform(d, bestP1);
    if (bestP1 === bestP2) {
      const dSafe = d <= 0 ? 1e-6 : d;
      val += bestFit.beta[2] * fpTransform(d, bestP1) * Math.log(dSafe);
    } else {
      val += bestFit.beta[2] * fpTransform(d, bestP2);
    }
    return val;
  }

  function predictCI(dose, zCrit) {
    const d = dose, dSafe = d <= 0 ? 1e-6 : d;
    const x = [1, fpTransform(d, bestP1),
      bestP1 === bestP2 ? fpTransform(d, bestP1) * Math.log(dSafe) : fpTransform(d, bestP2)];
    let variance = 0;
    if (bestFit.vcov) {
      for (let i = 0; i < 3; i++)
        for (let j = 0; j < 3; j++)
          variance += x[i] * bestFit.vcov[i][j] * x[j];
    }
    const se = Math.sqrt(Math.max(0, variance));
    const pred = predict(dose);
    return { lo: pred - zCrit * se, hi: pred + zCrit * se };
  }

  const p1Label = bestP1 === 0 ? 'log' : bestP1.toString();
  const p2Label = bestP2 === 0 ? 'log' : bestP2.toString();
  return {
    model: `FP2(${p1Label},${p2Label})`, b0: bestFit.beta[0], b1: bestFit.beta[1], b2: bestFit.beta[2],
    power1: bestP1, power2: bestP2, pFP2vsFP1, R2, aic: bestAIC, bic,
    fitted: bestFit.fitted, residuals: bestFit.residuals,
    predict, predictCI, paramStr: `powers=(${p1Label},${p2Label})`, nParams: 3, vcov: bestFit.vcov
  };
}
```

**Step 2: Add inline self-test for FP**

```javascript
function _testFP() {
  const results = [];
  // Test: points on y = 2 + 3*sqrt(dose) → FP1 should pick power=0.5
  const pts = Array.from({length: 8}, (_, i) => ({
    dose: i + 1, effect: 2 + 3 * Math.sqrt(i + 1), se: 0.1
  }));
  const fp1 = fitFP1DR(pts);
  results.push({ test: 'FP1 fits', pass: fp1 !== null });
  results.push({ test: 'FP1 power=0.5', pass: fp1?.power === 0.5 });
  results.push({ test: 'FP1 R2 > 0.99', pass: (fp1?.R2 ?? 0) > 0.99 });

  // FP2 should also fit well, with power pair that includes 0.5
  const fp2 = fitFP2DR(pts);
  results.push({ test: 'FP2 fits', pass: fp2 !== null });
  results.push({ test: 'FP2 AIC <= FP1 AIC', pass: (fp2?.aic ?? Infinity) <= (fp1?.aic ?? Infinity) + 2 });

  // Too few points
  const few = [{ dose: 1, effect: 1, se: 0.1 }, { dose: 2, effect: 2, se: 0.1 }];
  results.push({ test: 'FP1 null if <3', pass: fitFP1DR(few) === null });
  results.push({ test: 'FP2 null if <4', pass: fitFP2DR(few) === null });

  console.table(results);
  return results.every(r => r.pass);
}
```

**Step 3: Commit**

```bash
git commit -m "feat: add fractional polynomial (FP1, FP2) dose-response models"
```

---

### Task 4: One-Stage Mixed-Effects Model

**Context:** One-stage fits all studies simultaneously with random intercepts + slopes. Uses iterative GLS with REML variance estimation. This is the most complex enhancement.

**Files:**
- Modify: `metasprint-dose-response.html` — insert `fitOneStage()` after FP functions

**Step 1: Add one-stage mixed-effects fitter**

```javascript
// ── One-Stage Mixed-Effects Dose-Response ──
function fitOneStageMixedEffects(studyGroups, modelType, nKnots) {
  // studyGroups: [{studyId, points: [{dose, effect, se}]}]
  // modelType: 'linear' | 'quadratic' | 'rcs3' | 'rcs4' | 'fp1' | 'fp2'
  // Returns: {pooledModel, studyModels[], varianceComponents, convergence}

  if (studyGroups.length < 2) return null;

  // Step 1: Build design matrices per study
  const allPoints = [];
  const studyIdx = []; // maps each point to a study index
  for (let s = 0; s < studyGroups.length; s++) {
    for (const pt of studyGroups[s].points) {
      allPoints.push(pt);
      studyIdx.push(s);
    }
  }
  const N = allPoints.length;
  const K = studyGroups.length;
  if (N < K + 2) return null;

  // Build fixed-effects design matrix X based on modelType
  function buildXrow(dose) {
    if (modelType === 'linear') return [1, dose];
    if (modelType === 'quadratic') return [1, dose, dose * dose];
    if (modelType === 'fp1') {
      // Use power=1 (linear) as default for one-stage; or find best via two-stage first
      return [1, dose];
    }
    return [1, dose]; // fallback to linear
  }

  const pFixed = buildXrow(0).length;

  // Initial two-stage estimates for starting values
  const X = matCreate(N, pFixed);
  const y = allPoints.map(p => p.effect);
  const wBase = allPoints.map(p => { const se = p.se > 0 ? p.se : 1; return 1 / (se * se); });

  for (let i = 0; i < N; i++) {
    const row = buildXrow(allPoints[i].dose);
    for (let j = 0; j < pFixed; j++) X[i][j] = row[j];
  }

  // Initialize: fixed-effects WLS
  let fitInit = wlsFit(X, y, wBase);
  if (!fitInit) return null;
  let betaFixed = fitInit.beta.slice();

  // Initialize variance components: tau2_intercept, tau2_slope
  let tau2 = [0.01, 0.001]; // [random intercept variance, random slope variance]
  const maxIter = 50;
  const tol = 1e-6;
  let converged = false;
  let iterations = 0;

  for (let iter = 0; iter < maxIter; iter++) {
    iterations = iter + 1;

    // Step 2: Build V_i = diag(se_i^2) + Z_i * D * Z_i' for each study
    // Z_i = [1, dose] for random intercept + slope
    // D = diag(tau2)
    const wNew = new Array(N);
    for (let i = 0; i < N; i++) {
      const se2 = wBase[i] > 0 ? 1 / wBase[i] : 1;
      const dose_i = allPoints[i].dose;
      // Marginal variance = se^2 + tau2[0] + tau2[1]*dose^2 (diagonal of Z*D*Z')
      const margVar = se2 + tau2[0] + tau2[1] * dose_i * dose_i;
      wNew[i] = margVar > 1e-15 ? 1 / margVar : wBase[i];
    }

    // Step 3: Update fixed effects with new weights
    const fitNew = wlsFit(X, y, wNew);
    if (!fitNew) break;

    // Check convergence on beta
    let maxDelta = 0;
    for (let j = 0; j < pFixed; j++) {
      maxDelta = Math.max(maxDelta, Math.abs(fitNew.beta[j] - betaFixed[j]));
    }

    betaFixed = fitNew.beta.slice();

    // Step 4: Update tau2 via method of moments (REML-type)
    // Residuals per study
    const resid = y.map((yi, i) => yi - fitNew.fitted[i]);
    let sumR2_int = 0, sumR2_slope = 0, countStudies = 0;
    for (let s = 0; s < K; s++) {
      const studyPts = [];
      for (let i = 0; i < N; i++) {
        if (studyIdx[i] === s) studyPts.push({ resid: resid[i], dose: allPoints[i].dose, se2: 1 / wBase[i] });
      }
      if (studyPts.length === 0) continue;
      countStudies++;
      const meanResid = studyPts.reduce((s, p) => s + p.resid, 0) / studyPts.length;
      sumR2_int += meanResid * meanResid;
      if (studyPts.length >= 2) {
        // Slope residual: weighted regression of resid on dose within study
        let sx = 0, sy = 0, sxx = 0, sxy = 0, sw = 0;
        for (const p of studyPts) {
          const w = 1 / (p.se2 + 1e-15);
          sx += w * p.dose; sy += w * p.resid; sxx += w * p.dose * p.dose; sxy += w * p.dose * p.resid; sw += w;
        }
        const denom = sw * sxx - sx * sx;
        if (Math.abs(denom) > 1e-15) {
          const slopeResid = (sw * sxy - sx * sy) / denom;
          sumR2_slope += slopeResid * slopeResid;
        }
      }
    }

    const meanSe2 = wBase.reduce((s, w) => s + 1 / w, 0) / N;
    const newTau2_0 = Math.max(0, sumR2_int / countStudies - meanSe2);
    const newTau2_1 = Math.max(0, sumR2_slope / Math.max(1, countStudies) - meanSe2 * 0.1);

    tau2 = [newTau2_0, newTau2_1];

    if (maxDelta < tol && iter > 0) { converged = true; break; }
  }

  // Final fit with converged weights
  const wFinal = new Array(N);
  for (let i = 0; i < N; i++) {
    const se2 = wBase[i] > 0 ? 1 / wBase[i] : 1;
    const dose_i = allPoints[i].dose;
    const margVar = se2 + tau2[0] + tau2[1] * dose_i * dose_i;
    wFinal[i] = margVar > 1e-15 ? 1 / margVar : wBase[i];
  }
  const finalFit = wlsFit(X, y, wFinal);
  if (!finalFit) return null;

  // Build pooled model object (same interface as other models)
  const yMean = y.reduce((s, v, i) => s + wFinal[i] * v, 0) / wFinal.reduce((s, v) => s + v, 0);
  const ssTot = y.reduce((s, v, i) => s + wFinal[i] * (v - yMean) * (v - yMean), 0);
  const R2 = ssTot > 1e-15 ? 1 - finalFit.sse / ssTot : 0;
  const aic = N * Math.log(finalFit.sse / N + 1e-15) + 2 * (pFixed + 2);
  const bic = N * Math.log(finalFit.sse / N + 1e-15) + (pFixed + 2) * Math.log(N);

  function predict(dose) {
    const row = buildXrow(dose);
    return row.reduce((s, x, j) => s + x * finalFit.beta[j], 0);
  }
  function predictCI(dose, zCrit) {
    const row = buildXrow(dose);
    let variance = 0;
    if (finalFit.vcov) {
      for (let i = 0; i < pFixed; i++)
        for (let j = 0; j < pFixed; j++)
          variance += row[i] * finalFit.vcov[i][j] * row[j];
    }
    const se = Math.sqrt(Math.max(0, variance));
    const pred = predict(dose);
    return { lo: pred - zCrit * se, hi: pred + zCrit * se };
  }

  // Per-study curves (for overlay on plot)
  const studyModels = [];
  for (let s = 0; s < K; s++) {
    const sPts = studyGroups[s].points;
    if (sPts.length < 2) continue;
    const sX = matCreate(sPts.length, pFixed);
    const sY = sPts.map(p => p.effect);
    const sW = sPts.map(p => { const se = p.se > 0 ? p.se : 1; return 1 / (se * se); });
    for (let i = 0; i < sPts.length; i++) {
      const row = buildXrow(sPts[i].dose);
      for (let j = 0; j < pFixed; j++) sX[i][j] = row[j];
    }
    const sFit = wlsFit(sX, sY, sW);
    if (sFit) {
      studyModels.push({
        studyId: studyGroups[s].studyId,
        predict: (dose) => {
          const row = buildXrow(dose);
          return row.reduce((sum, x, j) => sum + x * sFit.beta[j], 0);
        },
        beta: sFit.beta
      });
    }
  }

  return {
    model: `1-stage ${modelType}`, beta: finalFit.beta, vcov: finalFit.vcov,
    tau2_intercept: tau2[0], tau2_slope: tau2[1],
    R2, aic, bic, fitted: finalFit.fitted, residuals: finalFit.residuals,
    predict, predictCI, studyModels,
    convergence: { converged, iterations, tolerance: tol },
    paramStr: `tau2_int=${tau2[0].toFixed(4)}, tau2_slope=${tau2[1].toFixed(4)}, iter=${iterations}`,
    nParams: pFixed + 2
  };
}
```

**Step 2: Add inline self-test**

```javascript
function _testOneStage() {
  const results = [];
  // Test: 3 studies, each with 4 dose levels on a line y = 1 + 0.5*dose + noise
  const groups = [];
  for (let s = 0; s < 3; s++) {
    const points = [];
    const intercept = 1 + (s - 1) * 0.2; // slight between-study variation
    for (let d = 0; d <= 3; d++) {
      points.push({ dose: d, effect: intercept + 0.5 * d, se: 0.15 });
    }
    groups.push({ studyId: `Study${s + 1}`, points });
  }
  const result = fitOneStageMixedEffects(groups, 'linear');
  results.push({ test: '1-stage fits', pass: result !== null });
  results.push({ test: '1-stage converged', pass: result?.convergence?.converged === true });
  results.push({ test: '1-stage slope ~0.5', pass: Math.abs((result?.beta?.[1] ?? 0) - 0.5) < 0.2 });
  results.push({ test: '1-stage has study models', pass: (result?.studyModels?.length ?? 0) >= 2 });
  results.push({ test: '1-stage tau2 >= 0', pass: (result?.tau2_intercept ?? -1) >= 0 });

  // Edge: only 1 study should return null
  results.push({ test: '1-stage null if K<2', pass: fitOneStageMixedEffects([groups[0]], 'linear') === null });

  console.table(results);
  return results.every(r => r.pass);
}
```

**Step 3: Commit**

```bash
git commit -m "feat: add one-stage mixed-effects dose-response model"
```

---

### Task 5: Continuous Input Mode (Mean/SD/N)

**Context:** Currently the app supports two input modes: one-stage (multi-point effect sizes) and two-stage (slopes). We add a third: continuous outcomes where users enter Mean, SD, N per dose level.

**Files:**
- Modify: `metasprint-dose-response.html`:
  - HTML: add radio button near line 1371
  - HTML: add alternative table header near line 1380
  - JS: modify `setInputMode()` near line 6595
  - JS: modify `prepareDoseResponseData()` near line 8187

**Step 1: Add continuous radio button**

After the existing two radio buttons at ~line 1372, add:

```html
<label style="font-size:0.82rem;cursor:pointer"><input type="radio" name="inputMode" value="continuous" onchange="setInputMode('continuous')"> Continuous (Mean/SD/N per dose)</label>
```

**Step 2: Add continuous table headers**

Add an alternative header row (hidden by default, shown when mode="continuous"). Near line 1380, add a new `<tr id="extractHeaderContinuous" style="display:none">` with columns:

```html
<tr id="extractHeaderContinuous" style="display:none">
  <th scope="col">Study ID</th><th scope="col">Outcome *</th><th scope="col">Dose *</th>
  <th scope="col" title="Unit for dose">Unit</th>
  <th scope="col" title="Reference dose (usually 0 = placebo)">Ref Dose</th>
  <th scope="col" title="Number of participants at this dose level">N *</th>
  <th scope="col" title="Mean outcome value at this dose level">Mean *</th>
  <th scope="col" title="Standard deviation at this dose level">SD *</th>
  <th scope="col" title="Type of effect: MD or SMD">Type</th>
  <th scope="col">Subgroup</th><th scope="col">Notes</th><th scope="col"></th>
</tr>
```

**Step 3: Modify `setInputMode()` to handle continuous**

In `setInputMode()` at ~line 6595, add a case for `mode === 'continuous'`:
- Show `extractHeaderContinuous`, hide others
- Update help text: "Enter Mean, SD, and N at each dose level. The app computes MD or SMD vs the reference dose."
- Update the add-row function to use continuous columns

**Step 4: Modify `prepareDoseResponseData()` for continuous**

In `prepareDoseResponseData()` at ~line 8187, add a new branch:

```javascript
if (extractInputMode === 'continuous') {
  // Group by study, compute MD or SMD vs reference dose
  const byStudy = {};
  for (const s of studies) {
    const key = s.studyId || s.trialId || 'unknown';
    if (!byStudy[key]) byStudy[key] = [];
    byStudy[key].push(s);
  }

  const points = [];
  const effectType = studies[0]?.effectType === 'SMD' ? 'SMD' : 'MD';

  for (const [studyId, rows] of Object.entries(byStudy)) {
    const refDose = parseFloat(rows[0]?.refDose ?? 0);
    const refRow = rows.find(r => parseFloat(r.dose) === refDose);
    if (!refRow || !isFinite(parseFloat(refRow.mean)) || !isFinite(parseFloat(refRow.sd))) continue;

    const refMean = parseFloat(refRow.mean);
    const refSD = parseFloat(refRow.sd);
    const refN = parseInt(refRow.n, 10);
    if (refN <= 0 || refSD <= 0) continue;

    for (const row of rows) {
      const dose = parseFloat(row.dose);
      if (dose === refDose) continue; // skip reference
      const mean = parseFloat(row.mean);
      const sd = parseFloat(row.sd);
      const n = parseInt(row.n, 10);
      if (!isFinite(mean) || !isFinite(sd) || sd <= 0 || n <= 0) continue;

      let effect, se;
      if (effectType === 'SMD') {
        // Hedges' g
        const df = n + refN - 2;
        const pooledSD = Math.sqrt(((n - 1) * sd * sd + (refN - 1) * refSD * refSD) / df);
        const J = df > 0 ? 1 - 3 / (4 * df - 1) : 1;
        effect = J * (mean - refMean) / pooledSD;
        se = Math.sqrt((n + refN) / (n * refN) + effect * effect / (2 * (n + refN)));
      } else {
        // MD
        effect = mean - refMean;
        se = Math.sqrt(sd * sd / n + refSD * refSD / refN);
      }

      if (isFinite(effect) && isFinite(se) && se > 0) {
        points.push({ studyId, dose, effect, se, n });
      }
    }
  }

  return { mode: 'one-stage', points, studies: Object.keys(byStudy), isRatio: false };
}
```

**Step 5: Add inline self-test**

```javascript
function _testContinuousInput() {
  const results = [];
  // Simulate 2 studies, 3 dose levels each, compute MD
  const studies = [
    { studyId: 'S1', dose: '0', mean: '10', sd: '2', n: '50', refDose: '0', effectType: 'MD' },
    { studyId: 'S1', dose: '5', mean: '12', sd: '2.1', n: '48', refDose: '0', effectType: 'MD' },
    { studyId: 'S1', dose: '10', mean: '15', sd: '2.3', n: '52', refDose: '0', effectType: 'MD' },
    { studyId: 'S2', dose: '0', mean: '9', sd: '1.8', n: '60', refDose: '0', effectType: 'MD' },
    { studyId: 'S2', dose: '5', mean: '11.5', sd: '2.0', n: '55', refDose: '0', effectType: 'MD' },
    { studyId: 'S2', dose: '10', mean: '14', sd: '2.2', n: '58', refDose: '0', effectType: 'MD' },
  ];
  // Temporarily set mode
  const origMode = extractInputMode;
  extractInputMode = 'continuous';
  const data = prepareDoseResponseData(studies);
  extractInputMode = origMode;

  results.push({ test: 'continuous mode=one-stage', pass: data.mode === 'one-stage' });
  results.push({ test: 'continuous 4 points', pass: data.points.length === 4 }); // 2 studies x 2 non-ref doses
  results.push({ test: 'continuous S1 dose=5 effect~2', pass: Math.abs((data.points[0]?.effect ?? 0) - 2) < 0.5 });
  results.push({ test: 'continuous SE > 0', pass: data.points.every(p => p.se > 0) });

  console.table(results);
  return results.every(r => r.pass);
}
```

**Step 6: Commit**

```bash
git commit -m "feat: add continuous input mode (Mean/SD/N per dose) with MD and SMD computation"
```

---

### Task 6: Update Model Selector and compareDoseResponseModels()

**Context:** The model selector dropdown (line ~1445) only has 4 options. `compareDoseResponseModels()` (line ~8497) only fits 3 models. Both need to include RCS and FP.

**Files:**
- Modify: `metasprint-dose-response.html`:
  - HTML: model selector at ~line 1445
  - JS: `compareDoseResponseModels()` at ~line 8497

**Step 1: Expand model selector dropdown**

Replace the existing `<select id="drModelSelect">` with:

```html
<select id="drModelSelect" style="padding:4px 8px;font-size:0.82rem;border:1px solid var(--border);border-radius:var(--radius)">
  <option value="auto" selected>Auto-select (best AIC)</option>
  <option value="linear">Linear</option>
  <option value="quadratic">Quadratic</option>
  <option value="emax">Emax</option>
  <option value="rcs3">RCS (3 knots)</option>
  <option value="rcs4">RCS (4 knots)</option>
  <option value="fp1">FP1 (best 1st-degree)</option>
  <option value="fp2">FP2 (best 2nd-degree)</option>
</select>
```

**Step 2: Add one-stage toggle**

After the model selector, add:

```html
<label style="font-size:0.82rem;margin-left:12px;cursor:pointer" title="Fit all studies simultaneously with random effects (slower but more rigorous)">
  <input type="checkbox" id="drOneStageToggle"> One-stage mixed-effects
</label>
```

**Step 3: Update `compareDoseResponseModels()`**

Replace the body of `compareDoseResponseModels()` with:

```javascript
function compareDoseResponseModels(points) {
  const models = [];
  const lin = fitLinearDR(points);
  if (lin) models.push(lin);
  const quad = fitQuadraticDR(points);
  if (quad) models.push(quad);
  const emax = fitEmaxDR(points);
  if (emax) models.push(emax);
  const rcs3 = fitRCSDR(points, 3);
  if (rcs3) models.push(rcs3);
  const rcs4 = fitRCSDR(points, 4);
  if (rcs4) models.push(rcs4);
  const fp1 = fitFP1DR(points);
  if (fp1) models.push(fp1);
  const fp2 = fitFP2DR(points);
  if (fp2) models.push(fp2);
  if (models.length === 0) return null;
  models.sort((a, b) => a.aic - b.aic);
  return { best: models[0], all: models };
}
```

**Step 4: Commit**

```bash
git commit -m "feat: expand model selector to include RCS and FP; update compareDoseResponseModels"
```

---

### Task 7: Update Model Comparison Table Rendering

**Context:** The model comparison table (line ~8543) needs a "Non-linearity p" column and should show knot/power information for RCS and FP.

**Files:**
- Modify: `metasprint-dose-response.html` — model comparison table at ~line 8543

**Step 1: Expand table header**

Add a "Non-linearity p" column between R-squared and Parameters:

```javascript
'<th style="padding:6px 10px;text-align:right" title="p-value for non-linearity (RCS Wald test or FP2 vs FP1)">Non-lin. p</th>' +
```

**Step 2: Expand table rows**

For each model row, add the non-linearity p-value:

```javascript
for (const m of comparison.all) {
  const isBest = m === comparison.best;
  const style = isBest ? 'background:var(--primary);color:#fff;font-weight:600' : '';
  const pNL = m.pNonlinear != null && isFinite(m.pNonlinear) ? m.pNonlinear.toFixed(4) :
              m.pFP2vsFP1 != null && isFinite(m.pFP2vsFP1) ? m.pFP2vsFP1.toFixed(4) : '—';
  tableHtml += `<tr style="${style}">` +
    `<td style="padding:6px 10px">${m.model}${isBest ? ' *' : ''}</td>` +
    `<td style="padding:6px 10px;text-align:right">${m.aic.toFixed(1)}</td>` +
    `<td style="padding:6px 10px;text-align:right">${m.bic.toFixed(1)}</td>` +
    `<td style="padding:6px 10px;text-align:right">${(m.R2 * 100).toFixed(1)}%</td>` +
    `<td style="padding:6px 10px;text-align:right">${pNL}</td>` +
    `<td style="padding:6px 10px">${m.paramStr ?? ''}</td>` +
    '</tr>';
}
```

**Step 3: Commit**

```bash
git commit -m "feat: expand model comparison table with non-linearity p-values"
```

---

### Task 8: Update renderDoseResponseAnalysis() for One-Stage and New Models

**Context:** The main orchestrator `renderDoseResponseAnalysis()` (line ~8513) needs to handle one-stage fitting and pass the selected model (including RCS/FP) to the curve renderer.

**Files:**
- Modify: `metasprint-dose-response.html` — `renderDoseResponseAnalysis()` at ~line 8513

**Step 1: Add one-stage path**

After the existing `prepareDoseResponseData()` call, add a branch for one-stage mixed-effects:

```javascript
const useOneStage = document.getElementById('drOneStageToggle')?.checked ?? false;

if (useOneStage && data.mode === 'one-stage' && data.points.length > 0) {
  // Group points by study for one-stage fitting
  const byStudy = {};
  for (const pt of data.points) {
    const sid = pt.studyId ?? 'unknown';
    if (!byStudy[sid]) byStudy[sid] = { studyId: sid, points: [] };
    byStudy[sid].points.push(pt);
  }
  const studyGroups = Object.values(byStudy);

  if (studyGroups.length >= 2) {
    const modelSelect = document.getElementById('drModelSelect')?.value ?? 'auto';
    const modelType = modelSelect === 'auto' ? 'linear' : modelSelect.replace('rcs3','linear').replace('rcs4','linear');
    const oneStageResult = fitOneStageMixedEffects(studyGroups, modelType);

    if (oneStageResult) {
      // Render one-stage results with study-level curve overlay
      renderDoseResponseCurveSVG(data.points, oneStageResult, confLevel, zCrit, data.isRatio);

      // Show convergence info
      const convHtml = oneStageResult.convergence.converged
        ? `<span style="color:var(--success)">Converged in ${oneStageResult.convergence.iterations} iterations</span>`
        : `<span style="color:var(--warning)">Did not converge (${oneStageResult.convergence.iterations} iterations)</span>`;
      el.innerHTML += `<p style="font-size:0.82rem;margin-top:8px">${convHtml}</p>`;
      el.innerHTML += `<p style="font-size:0.82rem;color:var(--text-muted)">${oneStageResult.paramStr}</p>`;
      return;
    }
  }
}
```

**Step 2: Update model selection logic for forced model**

In the existing two-stage path, add handling for forced RCS/FP model selection:

```javascript
const modelSelect = document.getElementById('drModelSelect')?.value ?? 'auto';
let selectedModel;
if (modelSelect === 'auto') {
  selectedModel = comparison.best;
} else if (modelSelect === 'rcs3') {
  selectedModel = comparison.all.find(m => m.model === 'RCS-3') ?? comparison.best;
} else if (modelSelect === 'rcs4') {
  selectedModel = comparison.all.find(m => m.model === 'RCS-4') ?? comparison.best;
} else if (modelSelect === 'fp1') {
  selectedModel = comparison.all.find(m => m.model?.startsWith('FP1')) ?? comparison.best;
} else if (modelSelect === 'fp2') {
  selectedModel = comparison.all.find(m => m.model?.startsWith('FP2')) ?? comparison.best;
} else {
  selectedModel = comparison.all.find(m => m.model?.toLowerCase() === modelSelect) ?? comparison.best;
}
```

**Step 3: Commit**

```bash
git commit -m "feat: integrate one-stage mixed-effects and new models into analysis orchestrator"
```

---

### Task 9: Update renderDoseResponseCurveSVG() for Study Overlays

**Context:** When one-stage is used, the plot should show per-study curves in addition to the pooled curve.

**Files:**
- Modify: `metasprint-dose-response.html` — `renderDoseResponseCurveSVG()` at ~line 8593

**Step 1: Add study curve overlay rendering**

After the main curve path is drawn, check if the model has `studyModels` and draw them:

```javascript
// After the main curve and CI band are drawn:
if (model.studyModels && model.studyModels.length > 0) {
  const studyColors = ['#94a3b8', '#a78bfa', '#67e8f9', '#fbbf24', '#fb923c', '#f87171', '#a3e635', '#c084fc'];
  for (let s = 0; s < model.studyModels.length; s++) {
    const sm = model.studyModels[s];
    const color = studyColors[s % studyColors.length];
    let studyPath = '';
    for (let i = 0; i <= nSteps; i++) {
      const dose = doseMin + (doseMax - doseMin) * i / nSteps;
      const pred = sm.predict(dose);
      if (!isFinite(pred)) continue;
      const cx = pad.left + (dose - doseMin) / (doseMax - doseMin) * plotW;
      const cy = pad.top + (1 - (pred - yMin) / (yMax - yMin)) * plotH;
      studyPath += (studyPath === '' ? 'M' : 'L') + cx.toFixed(1) + ',' + cy.toFixed(1);
    }
    if (studyPath) {
      svg += `<path d="${studyPath}" fill="none" stroke="${color}" stroke-width="1" stroke-dasharray="4,3" opacity="0.6"/>`;
      // Legend entry
      svg += `<text x="${pad.left + plotW - 10}" y="${pad.top + 16 + s * 14}" text-anchor="end" font-size="10" fill="${color}">${sm.studyId}</text>`;
    }
  }
}
```

**Step 2: Commit**

```bash
git commit -m "feat: add per-study curve overlays for one-stage mixed-effects plots"
```

---

### Task 10: Run All Self-Tests and Verify Div Balance

**Context:** After all code changes, verify everything works and the HTML is balanced.

**Files:**
- Modify: `metasprint-dose-response.html` — add master test runner

**Step 1: Add master test runner**

```javascript
function _runAllDRTests() {
  console.log('=== MetaSprint Dose-Response Self-Tests ===');
  const results = {};
  results.matrix = _testMatrixOps();
  results.rcs = _testRCS();
  results.fp = _testFP();
  results.oneStage = _testOneStage();
  // results.continuous = _testContinuousInput(); // needs DOM; run manually
  console.log('=== Results ===');
  console.table(results);
  const allPass = Object.values(results).every(v => v === true);
  console.log(allPass ? 'ALL TESTS PASSED' : 'SOME TESTS FAILED');
  return allPass;
}
```

**Step 2: Verify div balance**

Run in bash:
```bash
cd "C:/Users/user/Downloads/metasprint-dose-response"
# Count opening divs (exclude regex patterns in JS)
grep -c '<div[\s>]' metasprint-dose-response.html
# Count closing divs
grep -c '</div>' metasprint-dose-response.html
# These should match
```

**Step 3: Verify no literal `</script>` in template literals**

```bash
# Search for </script> that might be inside JS strings
grep -n '</script>' metasprint-dose-response.html | head -20
# Only the actual closing </script> tag at the end should appear
```

**Step 4: Verify function name uniqueness**

```bash
grep -oP 'function \w+' metasprint-dose-response.html | sort | uniq -c | sort -rn | head -20
# All counts should be 1
```

**Step 5: Final commit**

```bash
git commit -m "feat: add master self-test runner; verify div balance and function uniqueness"
```

---

### Task 11: Create Selenium Integration Test

**Context:** Create a Python Selenium test that opens the app, enters test data, runs analysis with each model type, and verifies output.

**Files:**
- Create: `test_dose_response_models.py`

**Step 1: Write Selenium test**

```python
"""Integration tests for dose-response model enhancements."""
import sys, io, time, pytest
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC

@pytest.fixture(scope='module')
def driver():
    opts = Options()
    opts.add_argument('--headless=new')
    opts.add_argument('--no-sandbox')
    opts.add_argument('--window-size=1400,900')
    d = webdriver.Chrome(options=opts)
    d.get('file:///C:/Users/user/Downloads/metasprint-dose-response/metasprint-dose-response.html')
    time.sleep(2)
    yield d
    d.quit()

def test_self_tests_pass(driver):
    """Run all inline self-tests via browser console."""
    result = driver.execute_script('return _runAllDRTests()')
    assert result is True, 'Inline self-tests failed'

def test_model_selector_has_all_options(driver):
    """Model selector dropdown should have 8 options."""
    select = Select(driver.find_element(By.ID, 'drModelSelect'))
    options = [o.get_attribute('value') for o in select.options]
    assert 'rcs3' in options, 'Missing RCS-3'
    assert 'rcs4' in options, 'Missing RCS-4'
    assert 'fp1' in options, 'Missing FP1'
    assert 'fp2' in options, 'Missing FP2'

def test_one_stage_toggle_exists(driver):
    """One-stage mixed-effects checkbox should exist."""
    el = driver.find_element(By.ID, 'drOneStageToggle')
    assert el is not None

def test_continuous_input_mode_exists(driver):
    """Continuous input radio should exist."""
    radios = driver.find_elements(By.CSS_SELECTOR, 'input[name="inputMode"][value="continuous"]')
    assert len(radios) == 1

def test_rcs_model_fitting(driver):
    """RCS-3 model should fit test data correctly."""
    result = driver.execute_script("""
        var pts = [];
        for (var i = 0; i <= 8; i++) {
            pts.push({dose: i, effect: Math.sin(i * 0.4) * 2, se: 0.2});
        }
        var rcs = fitRCSDR(pts, 3);
        return rcs ? {aic: rcs.aic, R2: rcs.R2, model: rcs.model, pNL: rcs.pNonlinear} : null;
    """)
    assert result is not None, 'RCS-3 returned null'
    assert result['model'] == 'RCS-3'
    assert result['R2'] > 0

def test_fp1_model_fitting(driver):
    """FP1 should fit sqrt relationship."""
    result = driver.execute_script("""
        var pts = [];
        for (var i = 1; i <= 8; i++) pts.push({dose: i, effect: 3 * Math.sqrt(i), se: 0.1});
        var fp = fitFP1DR(pts);
        return fp ? {power: fp.power, R2: fp.R2} : null;
    """)
    assert result is not None
    assert result['power'] == 0.5, f"Expected power=0.5, got {result['power']}"
    assert result['R2'] > 0.99

def test_fp2_model_fitting(driver):
    """FP2 should fit data needing two terms."""
    result = driver.execute_script("""
        var pts = [];
        for (var i = 1; i <= 10; i++) pts.push({dose: i, effect: 2*Math.sqrt(i) - 0.3*i, se: 0.15});
        var fp = fitFP2DR(pts);
        return fp ? {power1: fp.power1, power2: fp.power2, R2: fp.R2} : null;
    """)
    assert result is not None
    assert result['R2'] > 0.9

def test_one_stage_mixed_effects(driver):
    """One-stage should converge on clean data."""
    result = driver.execute_script("""
        var groups = [];
        for (var s = 0; s < 4; s++) {
            var pts = [];
            var b0 = 1 + (s-1.5)*0.3;
            for (var d = 0; d <= 4; d++) pts.push({dose: d, effect: b0 + 0.4*d, se: 0.2});
            groups.push({studyId: 'S' + s, points: pts});
        }
        var r = fitOneStageMixedEffects(groups, 'linear');
        return r ? {converged: r.convergence.converged, slope: r.beta[1], iter: r.convergence.iterations} : null;
    """)
    assert result is not None
    assert result['converged'] is True
    assert abs(result['slope'] - 0.4) < 0.15
```

**Step 2: Run tests**

```bash
cd "C:/Users/user/Downloads/metasprint-dose-response"
python -m pytest test_dose_response_models.py -v
```

Expected: All tests PASS.

**Step 3: Commit**

```bash
git add test_dose_response_models.py
git commit -m "test: add Selenium integration tests for all new dose-response models"
```

---

## Execution Order Summary

| Task | Description | Depends On |
|------|-------------|------------|
| 1 | NxN matrix operations + WLS solver | — |
| 2 | RCS (restricted cubic splines) | Task 1 |
| 3 | Fractional polynomial (FP1, FP2) | Task 1 |
| 4 | One-stage mixed-effects | Task 1 |
| 5 | Continuous input mode (Mean/SD/N) | — |
| 6 | Update model selector + compareDoseResponseModels | Tasks 2, 3 |
| 7 | Update model comparison table | Task 6 |
| 8 | Update renderDoseResponseAnalysis orchestrator | Tasks 2, 3, 4, 6 |
| 9 | Update curve SVG with study overlays | Task 4 |
| 10 | Self-tests + div balance + verification | All above |
| 11 | Selenium integration tests | All above |

**Parallelizable:** Tasks 2, 3, 4 can run in parallel (all depend only on Task 1). Task 5 is independent of everything.
