// engine.mjs — VERBATIM-extracted pure functions from metasprint-dose-response.html
// ADDITIVE: only module-scope const LOG_FLOOR + export {} appended. Bodies unmodified.
const LOG_FLOOR = 1e-10;

  function safeLog(v) { return Math.log(Math.max(LOG_FLOOR, v)); }

  function normalCDF(x) {

    const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;

    const p = 0.3275911;

    const sign = x < 0 ? -1 : 1;

    x = Math.abs(x) / Math.sqrt(2);

    const t = 1 / (1 + p * x);

    const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

    return 0.5 * (1 + sign * y);

  }

  function normalQuantile(p) {

    // Rational approximation (Abramowitz & Stegun 26.2.23)

    if (p <= 0) return -Infinity;

    if (p >= 1) return Infinity;

    if (p === 0.5) return 0;

    const a = p < 0.5 ? p : 1 - p;

    const t = Math.sqrt(-2 * Math.log(a));

    const c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;

    const d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;

    let z = t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t);

    return p < 0.5 ? -z : z;

  }

  function computeMetaAnalysis(studies, confLevel, opts) {

    confLevel = confLevel ?? 0.95;

    const alpha = 1 - confLevel;

    const zCrit = normalQuantile(1 - alpha / 2);

    const useHKSJ = opts?.hksj ?? false; // Knapp-Hartung/Sidik-Jonkman modification



    const valid = studies.filter(s =>

      s.effectEstimate !== null && s.lowerCI !== null && s.upperCI !== null

    );

    if (valid.length === 0) return null;



    const isRatio = ['OR', 'RR', 'HR'].includes(valid[0].effectType);



    // P1-2: Studies report 95% CIs regardless of analysis confLevel

    const studyCiZ = normalQuantile(0.975);



    const data = valid.map(s => {

      // P0-1: Guard against log(0) for ratio measures

      if (isRatio && (s.effectEstimate <= 0 || s.lowerCI <= 0 || s.upperCI <= 0)) return null;

      const yi = isRatio ? Math.log(s.effectEstimate) : s.effectEstimate;

      const lo = isRatio ? Math.log(s.lowerCI) : s.lowerCI;

      const hi = isRatio ? Math.log(s.upperCI) : s.upperCI;

      const sei = (hi - lo) / (2 * studyCiZ);

      // P0-7: Guard against zero CI width (sei=0 causes division by zero)

      if (sei <= 0 || !isFinite(sei)) return null;

      return { ...s, yi, sei, vi: sei * sei, wi: 1 / (sei * sei) };

    }).filter(d => d !== null);



    if (data.length === 0) return null;



    // Fixed-effect estimate

    const sumW = data.reduce((a, d) => a + d.wi, 0);

    const muFE = data.reduce((a, d) => a + d.wi * d.yi, 0) / sumW;



    // Q statistic

    const Q = data.reduce((a, d) => a + d.wi * (d.yi - muFE) ** 2, 0);

    const df = data.length - 1;



    // DerSimonian-Laird tau-squared

    // Guard: C=0 when all weights are equal (degenerate); fall back to tau2=0

    const C = sumW - data.reduce((a, d) => a + d.wi * d.wi, 0) / sumW;

    const tau2 = df > 0 && C > 1e-15 ? Math.max(0, (Q - df) / C) : 0;



    // Random-effects weights

    const reData = data.map(d => {

      const wi_re = 1 / (d.vi + tau2);

      return { ...d, wi_re };

    });

    const sumW_re = reData.reduce((a, d) => a + d.wi_re, 0);

    const muRE = reData.reduce((a, d) => a + d.wi_re * d.yi, 0) / sumW_re;

    const seRE = Math.sqrt(1 / sumW_re);



    // I-squared (guard: Q=0 when all effects identical → 0/0; I2 is 0% by definition)

    const I2 = df > 0 && Q > 0 ? Math.max(0, (Q - df) / Q * 100) : (df > 0 ? 0 : null);



    // HKSJ/Knapp-Hartung modification (Hartung & Knapp 2001, Sidik & Jonkman 2002)

    // Replaces z-based CI with t-based CI using adjusted variance

    let seCI = seRE;  // SE used for confidence intervals

    let critVal = zCrit;  // Critical value for CIs

    if (useHKSJ && df > 0) {

      // q* = (1/(k-1)) * sum(w*_i * (y_i - mu_RE)^2)

      const qStar = reData.reduce((a, d) => a + d.wi_re * (d.yi - muRE) ** 2, 0) / df;

      // Apply max(1, q*) to prevent HKSJ from narrowing CIs vs DL

      const qAdj = Math.max(1, qStar);

      seCI = Math.sqrt(qAdj / sumW_re);

      critVal = tQuantile(1 - alpha / 2, df);

    }



    // z-test (or t-test under HKSJ)

    const z = muRE / seCI;

    const pValue = useHKSJ && df > 0

      ? 2 * (1 - tCDFfn(Math.abs(z), df))

      : 2 * (1 - normalCDF(Math.abs(z)));



    // Back-transform

    const pooled = isRatio ? Math.exp(muRE) : muRE;

    const pooledLo = isRatio ? Math.exp(muRE - critVal * seCI) : muRE - critVal * seCI;

    const pooledHi = isRatio ? Math.exp(muRE + critVal * seCI) : muRE + critVal * seCI;



    // Per-study weights (%)

    const totalW = reData.reduce((a, d) => a + d.wi_re, 0);

    const studyResults = reData.map(d => ({

      ...d,

      weightPct: (d.wi_re / totalW * 100).toFixed(1),

      display: isRatio ? Math.exp(d.yi) : d.yi,

      // Study-level CIs: always 95% (studyCiZ), independent of analysis confLevel

      displayLo: isRatio ? Math.exp(d.yi - studyCiZ * d.sei) : d.yi - studyCiZ * d.sei,

      displayHi: isRatio ? Math.exp(d.yi + studyCiZ * d.sei) : d.yi + studyCiZ * d.sei

    }));



    // Q-test p-value (chi-squared with df degrees of freedom)

    const QpValue = df > 0 ? 1 - chi2CDF(Q, df) : 1;



    // Prediction interval (t-distribution based)

    // Cochrane Handbook v6.5 (Jan 2025): uses k-1 df (updated from k-2)

    // Ref: Higgins et al. 2009, Riley et al. 2011

    let piLo = null, piHi = null;

    if (data.length >= 3) {

      const piDf = Math.max(1, data.length - 1);  // k-1 per Cochrane Handbook v6.5 (Jan 2025)

      const tCrit = tQuantile(1 - alpha / 2, piDf);

      const piSE = Math.sqrt(tau2 + seRE * seRE);

      piLo = isRatio ? Math.exp(muRE - tCrit * piSE) : muRE - tCrit * piSE;

      piHi = isRatio ? Math.exp(muRE + tCrit * piSE) : muRE + tCrit * piSE;

    }



    // REML-based I² (Cochrane Handbook v6.5, Jan 2025): tau²_REML/(tau²_REML + v_typical)

    // v_typical = (k-1)*sumW / (sumW² - sum(wi²)) — typical within-study variance

    // P0-2 fix: compute REML tau2 here (not externally) so I2_REML uses the correct value

    const sumW2 = data.reduce((s, d) => s + d.wi * d.wi, 0);

    const vTypical = df > 0 && (sumW * sumW - sumW2) > 0

      ? df * sumW / (sumW * sumW - sumW2) : null;

    const tau2REML = data.length >= 3 ? estimateREML(data) : tau2;

    const I2_REML = vTypical !== null ? Math.max(0, tau2REML / (tau2REML + vTypical) * 100) : I2;



    return {

      pooled, pooledLo, pooledHi, tau2, tau2REML, I2, I2_REML, Q, QpValue, df, pValue,

      k: data.length, isRatio, studyResults,

      muRE, seRE, seCI, muFE, confLevel, zCrit, piLo, piHi,

      method: useHKSJ ? 'DL+HKSJ' : 'DL'

    };

  }

  function tQuantile(p, df) {

    if (df <= 0) return normalQuantile(p);

    if (df === 1) return Math.tan(Math.PI * (p - 0.5)); // Cauchy exact

    if (df === 2) {  // Exact formula for df=2

      const a = 2 * p - 1;

      // Guard: when p≈0 or p≈1, a*a rounds to 1 → division by zero

      const denom = 1 - a * a;

      if (denom < 1e-15) return a > 0 ? 1e15 : -1e15;

      return a * Math.sqrt(2 / denom);

    }

    if (df >= 200) return normalQuantile(p); // normal approximation

    // Hybrid: Newton-Raphson with bisection fallback for robustness

    const sign = p >= 0.5 ? 1 : -1;

    const pp = p >= 0.5 ? p : 1 - p;  // work with upper tail

    // Initial guess from normal quantile, corrected for heavy tails

    let x = normalQuantile(pp);

    // Cornish-Fisher correction for small df

    if (df < 30) {

      const g1 = 1 / (4 * df);

      x = x + (x * x * x + x) * g1;

    }

    // Newton-Raphson with clamped steps

    let converged = false;

    for (let i = 0; i < 30; i++) {

      const cdf = tCDFfn(x, df);

      const pdf = Math.pow(1 + x * x / df, -(df + 1) / 2) / (Math.sqrt(df) * betaFn(0.5, df / 2));

      if (pdf < 1e-15) break;

      const step = (cdf - pp) / pdf;

      const clampedStep = Math.abs(step) > Math.abs(x) * 0.5 + 1

        ? Math.sign(step) * (Math.abs(x) * 0.5 + 1) : step;

      x -= clampedStep;

      if (Math.abs(step) < 1e-10) { converged = true; break; }

    }

    // Bisection fallback if Newton didn't converge

    if (!converged) {

      let lo = normalQuantile(pp);

      let hi = Math.max(lo * 3, 50);  // generous upper bound

      // Ensure bracket: tCDF(hi) > pp

      while (tCDFfn(hi, df) < pp && hi < 1e6) hi *= 2;

      for (let i = 0; i < 80; i++) {

        const mid = (lo + hi) / 2;

        if (tCDFfn(mid, df) < pp) lo = mid; else hi = mid;

        if (hi - lo < 1e-10) break;

      }

      x = (lo + hi) / 2;

    }

    return sign * x;

  }

  function tCDFfn(t, df) {

    const x = df / (df + t * t);

    const p = 0.5 * regIncBeta(df / 2, 0.5, x);

    return t >= 0 ? 1 - p : p;

  }

  function regIncBeta(a, b, x) {

    if (x <= 0) return 0;

    if (x >= 1) return 1;

    const lnBeta = lnGamma(a) + lnGamma(b) - lnGamma(a + b);

    const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - lnBeta);

    // Lentz's continued fraction

    let f = 1e-30, c = 1e-30, d = 0;

    for (let m = 0; m <= 200; m++) {

      let num;

      if (m === 0) num = 1;

      else if (m % 2 === 0) {

        const k = m / 2;

        num = k * (b - k) * x / ((a + 2 * k - 1) * (a + 2 * k));

      } else {

        const k = (m - 1) / 2;

        num = -((a + k) * (a + b + k) * x) / ((a + 2 * k) * (a + 2 * k + 1));

      }

      d = 1 + num * d; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1 / d;

      c = 1 + num / c; if (Math.abs(c) < 1e-30) c = 1e-30;

      f *= c * d;

      if (Math.abs(c * d - 1) < 1e-10) break;

    }

    return front * f / a;

  }

  function lnGamma(z) {

    const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,

              -1.231739572450155, 0.001208650973866179, -0.000005395239384953];

    let x = z, y = z, tmp = x + 5.5;

    tmp -= (x + 0.5) * Math.log(tmp);

    let ser = 1.000000000190015;

    for (let j = 0; j < 6; j++) ser += c[j] / ++y;

    return -tmp + Math.log(2.5066282746310005 * ser / x);

  }

  function betaFn(a, b) {

    return Math.exp(lnGamma(a) + lnGamma(b) - lnGamma(a + b));

  }

  function chi2CDF(x, df) {

    if (x <= 0 || df <= 0) return 0;

    return regIncGamma(df / 2, x / 2);

  }

  function regIncGamma(a, x) {

    if (x < 0 || a <= 0) return 0;

    if (x === 0) return 0;

    if (x < a + 1) {

      // Series expansion

      let sum = 1 / a, term = 1 / a;

      for (let n = 1; n < 200; n++) {

        term *= x / (a + n);

        sum += term;

        if (Math.abs(term) < 1e-10 * Math.abs(sum)) break;

      }

      return sum * Math.exp(-x + a * Math.log(x) - lnGamma(a));

    } else {

      // Continued fraction (Numerical Recipes, Press et al.)

      // Computes Q(a,x) = 1 - P(a,x) via Lentz's modified method

      let b = x + 1 - a;

      let c = 1e30;

      let d = 1 / b;

      let h = d;

      for (let i = 1; i <= 200; i++) {

        const an = -i * (i - a);

        b += 2;

        d = an * d + b; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1 / d;

        c = b + an / c; if (Math.abs(c) < 1e-30) c = 1e-30;

        const del = c * d;

        h *= del;

        if (Math.abs(del - 1) < 1e-10) break;

      }

      return 1 - Math.exp(-x + a * Math.log(x) - lnGamma(a)) * h;

    }

  }

  function estimateREML(studyData, maxIter, tol) {

    maxIter = maxIter ?? 50;

    tol = tol ?? 1e-5;

    const k = studyData.length;

    if (k < 2) return 0;



    // Start from DL estimate

    const ws = studyData.map(d => 1 / d.vi);

    const sumW = ws.reduce((a, w) => a + w, 0);

    const muFE = ws.reduce((a, w, i) => a + w * studyData[i].yi, 0) / sumW;

    const Q = ws.reduce((a, w, i) => a + w * (studyData[i].yi - muFE) ** 2, 0);

    const C = sumW - ws.reduce((a, w) => a + w * w, 0) / sumW;

    // Guard: C=0 when all weights are equal (degenerate case); fall back to tau2=0

    let tau2 = C > 1e-15 ? Math.max(0, (Q - (k - 1)) / C) : 0;



    for (let iter = 0; iter < maxIter; iter++) {

      const w = studyData.map(d => 1 / (d.vi + tau2));

      const sW = w.reduce((a, b) => a + b, 0);

      const mu = w.reduce((s, wi, i) => s + wi * studyData[i].yi, 0) / sW;



      // REML score (Viechtbauer 2005, eq. 12; +1/sW term is the REML bias correction

      // that distinguishes REML from ML — accounts for uncertainty in estimating mu)

      const num = w.reduce((s, wi, i) =>

        s + wi * wi * ((studyData[i].yi - mu) ** 2 - studyData[i].vi), 0);

      const sW2 = w.reduce((s, wi) => s + wi * wi, 0);

      // Guard: sW2 or sW underflow to 0 when tau2 is extremely large (all weights ≈ 0)

      if (sW2 < 1e-30 || sW < 1e-30) break;



      const tau2New = Math.max(0, num / sW2 + 1 / sW);

      if (Math.abs(tau2New - tau2) < tol) { tau2 = tau2New; break; }

      tau2 = tau2New;

    }

    return tau2;

  }

  function fitLinearDR(points) {

    if (points.length < 2) return null;

    const n = points.length;

    let sw = 0, swx = 0, swy = 0, swxx = 0, swxy = 0;

    for (const p of points) {

      const w = 1 / (p.se * p.se);

      sw += w; swx += w * p.dose; swy += w * p.effect;

      swxx += w * p.dose * p.dose; swxy += w * p.dose * p.effect;

    }

    const det = sw * swxx - swx * swx;

    if (Math.abs(det) < 1e-15) return null;

    const b0 = (swxx * swy - swx * swxy) / det;

    const b1 = (sw * swxy - swx * swy) / det;

    const se_b0 = Math.sqrt(swxx / det);

    const se_b1 = Math.sqrt(sw / det);

    const cov_01 = -swx / det; // off-diagonal of (X'WX)^{-1}

    // Goodness of fit

    let ssRes = 0, ssTot = 0;

    const yBar = swy / sw;

    const fitted = [];

    const residuals = [];

    for (const p of points) {

      const w = 1 / (p.se * p.se);

      const yHat = b0 + b1 * p.dose;

      fitted.push(yHat);

      residuals.push(p.effect - yHat);

      ssRes += w * (p.effect - yHat) * (p.effect - yHat);

      ssTot += w * (p.effect - yBar) * (p.effect - yBar);

    }

    const R2 = ssTot > 0 ? 1 - ssRes / ssTot : 0;

    const df = n - 2;

    const pValue = df > 0

      ? 2 * (1 - normalCDF(Math.abs(b1 / se_b1)))

      : NaN; // DR-10: guard zero df

    const k = 2; // parameters

    const aic = n * Math.log(ssRes / n + 1e-15) + 2 * k;

    const bic = n * Math.log(ssRes / n + 1e-15) + k * Math.log(n);

    return { model: 'Linear', b0, b1, se_b0, se_b1, pValue, R2, aic, bic, fitted, residuals, nPoints: n,

      predict: (d) => b0 + b1 * d,

      predictCI: (d, z) => {

        const yh = b0 + b1 * d;

        // DR-1 fix: include cov(b0,b1) cross-term

        const seh = Math.sqrt(se_b0 * se_b0 + 2 * d * cov_01 + se_b1 * se_b1 * d * d);

        return { y: yh, lo: yh - z * seh, hi: yh + z * seh };

      }

    };

  }

export { computeMetaAnalysis, fitLinearDR, estimateREML, normalQuantile, normalCDF, tQuantile, chi2CDF, safeLog };
