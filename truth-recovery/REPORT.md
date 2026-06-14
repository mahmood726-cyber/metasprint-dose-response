# Truth-Recovery Validation — metasprint-dose-response

## Verdict: GENUINE ENGINE — VALIDATED. No `dose-response-pro` tau2 bug. (Honest negative: known DL under-coverage under heterogeneity.)

## What was tested
The engine is a single 2.3 MB self-contained HTML app (`metasprint-dose-response.html`, ~25k lines inline JS).
**No WebAssembly** — all dose-response and pooling math is pure JS, fully testable headless.

We extracted, **verbatim**, the pure statistical core into `engine.mjs`:
- `computeMetaAnalysis` — random-effects pooling (DerSimonian-Laird, optional HKSJ), used by the app's
  **two-stage** mode to pool per-study dose-response slopes.
- `estimateREML` — REML tau2 estimator (EM iteration).
- `fitLinearDR` — within-study weighted-least-squares linear dose-response fit.
- Numeric helpers (`tQuantile`, `chi2CDF`, `normalCDF/Quantile`, incomplete beta/gamma) — all verbatim.

`dgp-dr.mjs` is a standalone seeded **binomial** dose-response DGP: each study draws a true per-unit-dose
log-RR slope `slope_i ~ N(muSlope, tau2)`, runs a reference arm + dose arms with counts
`a ~ Bin(n, p0*exp(slope_i*dose))`, and a within-study weighted log-linear fit emits the app's two-stage
`{slope, slopeSE}` row. **Estimand = true pooled mean slope `muSlope` and true between-study `tau2`.**

## The specific bug we checked for (sibling `dose-response-pro`)
That sibling's DL tau2 denominator was dimensionally wrong (`sumTrV - df`: sum of small point-variances
minus an integer -> always negative -> DL returns 0 -> model silently runs fixed-effect -> slope-CI
under-coverage collapse).

**This engine does NOT have that bug.** Line ~15128 uses the correct DL form:
`C = sumW - sum(wi^2)/sumW;  tau2 = max(0, (Q-df)/C)`. REML (line ~16558) uses the same correct `C`.
The harness confirms recovery empirically.

## Results (muSlope=0.10, k=12, 600 reps/scenario)

| Scenario | true tau2 | DL tau2_hat (mean) | REML tau2_hat | bias | 95% CI coverage | frac(tau2_hat=0) |
|---|---|---|---|---|---|---|
| Homogeneous   | 0.00000 | 0.00012 | 0.00012 | -0.0006 | 0.968 | 0.597 |
| Mild het      | 0.00050 | 0.00052 | 0.00052 | +0.0005 | 0.913 | 0.188 |
| Moderate het  | 0.00200 | 0.00200 | 0.00201 | +0.0016 | 0.903 | 0.008 |
| Strong het    | 0.00500 | 0.00504 | 0.00505 | +0.0022 | 0.900 | 0.000 |
| Strong + HKSJ | 0.00500 | -       | -       | +0.0022 | 0.927 | -    |

## Findings
1. **tau2 recovery is essentially exact** — DL and REML means track true tau2 across all levels
   (0.0005->0.005). Bug-signature (collapse to ~0 under real heterogeneity) is ABSENT; `frac(tau2_hat=0)`
   falls toward 0 as true heterogeneity rises, exactly as it should.
2. **Pooled slope is unbiased** (|bias| <= 0.0022 ~ 2% of the true slope).
3. **Coverage near-nominal under homogeneity (0.968)**; expected well-documented DL under-coverage under
   heterogeneity (~0.90 at k=12). Property of DL+z, not a defect. The engine already ships HKSJ, which
   moves coverage back up (0.900->0.927). Honest negative, as-measured.

## Recommendation
**Ship validation; no code fix required.** The dose-response pooling core is correct and the DL tau2
estimator is dimensionally sound (unlike `dose-response-pro`). For few-study / heterogeneous analyses,
prefer the already-implemented HKSJ option (and/or REML tau2); consider defaulting HKSJ for the two-stage
pool when k is small, since plain DL+z under-covers at ~k=12.

## Reproduce
```
node truth-recovery/harness.mjs
node --test truth-recovery/test-truth-recovery.mjs   # 5/5 pass
```
