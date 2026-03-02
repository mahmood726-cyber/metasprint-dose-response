# Exceed dosresmeta2 — Design Document

**Date**: 2026-03-02
**App**: `metasprint-dose-response.html` (~21,559 lines)
**Goal**: Achieve full parity with R's dosresmeta2 package, then exceed it with model-averaged predictions. Validate all results against R reference values.

## Context

The app currently has Linear, Quadratic, Emax, RCS (3/4 knot), FP1, FP2 models with AIC weights, LOO influence, Q_DR heterogeneity, and CSV import. However, it lacks the **Greenland-Longnecker covariance reconstruction** (dosresmeta's core innovation), proper **ML/REML estimation**, one-stage support beyond linear/quadratic, and prediction intervals.

### What We Already Exceed dosresmeta On
- AIC model averaging weights (Burnham & Anderson 2002)
- Automated model comparison table (8 models side-by-side)
- Leave-one-out DR influence analysis
- Auto-select best model by AIC
- DR Universe (landscape visualization, dose-ranging trial detection)

### Critical Gaps vs dosresmeta2
1. Greenland-Longnecker covariance reconstruction (THE core gap)
2. ML/REML estimation (we use method-of-moments proxy)
3. One-stage mixed-effects limited to linear/quadratic
4. Missing models: log-linear, exponential, 3-param Emax (Hill)
5. No prediction intervals on DR curve
6. No model-averaged predictions

## Section 1: Greenland-Longnecker Covariance Reconstruction

### What
Reconstruct within-study covariance between dose-level contrasts that share a reference group. Without this, CIs are overconfident because correlated contrasts are treated as independent.

### Implementation
- **Function**: `greenlandLongnecker(cases, n, dose, referenceIndex)`
  - Input: per-dose-level case counts, sample sizes, doses, and which level is reference
  - Output: covariance matrix S for the (k-1) non-reference contrasts
- **Formula** (for log-RR):
  - `var(logRR_i) = 1/cases_i - 1/n_i + 1/cases_ref - 1/n_ref`
  - `cov(logRR_i, logRR_j) = 1/cases_ref - 1/n_ref` (shared reference contribution)
- **For continuous outcomes (MD/SMD)**:
  - `cov(MD_i, MD_j) = var_ref / n_ref` (shared reference group variance)
- **Hamling alternative**: For case-control data where marginal totals are known
  - Reconstruct pseudo-cell counts from adjusted ORs + marginal totals
  - `hamlingCovariance(or, ci_lower, ci_upper, cases_total, controls_total)`
- **Integration**: Block-diagonal V matrix in GLS: `V = blockDiag(S_1, S_2, ..., S_k)`
  - Each study contributes a block; between-study blocks are zero
  - GLS: `β = (X'V⁻¹X)⁻¹ X'V⁻¹y` using full V instead of diagonal weights

### Data Flow
1. User enters case counts + sample sizes per dose level (or the app infers from effect/SE)
2. `greenlandLongnecker()` builds per-study covariance matrix
3. `blockDiagV()` assembles full V matrix
4. All model fitting functions updated to use GLS with full V

### UI
- New optional columns in Extract table: "Cases" and "N" per dose level
- Auto-detect: if cases/N provided, use GL covariance; otherwise fall back to diagonal (current behavior)
- Info tooltip explaining GL covariance and why it matters

## Section 2: ML/REML Estimation

### What
Replace method-of-moments tau² with proper ML/REML for valid log-likelihood, AIC/BIC, and likelihood ratio tests.

### Implementation
- **Profile log-likelihood** (for tau² scalar, assuming proportional heterogeneity):
  - `logL(τ²) = -0.5 * [log|V + τ²I| + (y - Xβ)'(V + τ²I)⁻¹(y - Xβ)]`
  - REML adds: `+ 0.5 * log|X'(V + τ²I)⁻¹X|`
- **Optimization**: Golden-section search on τ² ∈ [0, τ²_max]
  - τ²_max = 10 × max(diagonal of V)
  - Tolerance: 1e-8, max 100 iterations
  - Falls back to DL if optimization fails
- **Output**: τ²_ML or τ²_REML, log-likelihood, AIC = -2logL + 2p, BIC = -2logL + p·log(n)
- **LRT**: For nested models (e.g., linear vs quadratic): χ² = -2(logL_reduced - logL_full)

### UI
- Estimation method selector: "DL (default)" / "ML" / "REML"
- When ML/REML selected, AIC/BIC in model comparison table uses proper likelihood values
- LRT column added to model comparison when ML/REML active

## Section 3: One-Stage Mixed-Effects for All Models

### What
Extend `fitOneStageMixedEffects()` from linear/quadratic only to support RCS, FP1, FP2, Emax, and the new models.

### Implementation
- **Current state**: `fitOneStageMixedEffects()` builds design matrix X for linear (1 col) or quadratic (2 cols)
- **Extension**: Accept a `designMatrixFn(doses)` callback that returns X columns for any model type
  - Linear: `[dose]`
  - Quadratic: `[dose, dose²]`
  - RCS-3: `[dose, spline1(dose)]` (2 basis functions after boundary constraints)
  - RCS-4: `[dose, spline1(dose), spline2(dose)]`
  - FP1: `[dose^p]` (best power from FP1 grid search)
  - FP2: `[dose^p1, dose^p2]` or `[dose^p, dose^p·log(dose)]` when p1=p2
  - Emax: linearized via iterative GLS (Gauss-Newton on Emax, ED50)
  - Log-linear: `[log(dose + ε)]`
  - Exponential: linearized via iterative GLS
  - Hill: linearized via iterative GLS (Emax, ED50, h)
- **Random effects**: Per-study random intercept + random slope on first dose term
- **Convergence**: Same IGLS with REML, max 50 iterations, tol 1e-6

### UI
- One-stage toggle works for all model types (currently grayed out for RCS/FP/Emax)
- Per-study curves overlaid on pooled curve for any model

## Section 4: Additional Models

### What
Three new dose-response model types to match dosresmeta2's flexibility.

### Log-Linear
- **Formula**: `y = β₁ · log(dose + c)` where c = small constant (1e-6 or user-set)
- **Fitting**: WLS on transformed dose axis
- **Use case**: Common in epidemiology (alcohol-CVD, radiation exposure)

### Exponential
- **Formula**: `y = Emax · (1 - exp(-α · dose))`
- **Fitting**: Gauss-Newton iterative WLS (2 parameters: Emax, α)
- **Initialization**: Emax₀ = max observed effect, α₀ = -log(0.5) / median(dose)
- **Use case**: Saturation curves where effect plateaus

### 3-Parameter Emax (Hill)
- **Formula**: `y = Emax · dose^h / (ED50^h + dose^h)`
- **Fitting**: Gauss-Newton iterative WLS (3 parameters: Emax, ED50, h)
- **Initialization**: h₀ = 1 (reduces to standard Emax), then grid search h ∈ {0.5, 1, 1.5, 2, 3}
- **Use case**: Sigmoidal dose-response with variable steepness

### UI
- Three new entries in model selector dropdown
- Model comparison table shows all 11 models when data supports them
- Hill exponent h displayed in results summary

## Section 5: Prediction Intervals

### What
Display prediction intervals on the DR curve showing the range of plausible effects in a new study.

### Implementation
- **Formula**: `PI(d) = ŷ(d) ± t_{k-p, α/2} × √(SE²(d) + τ²)`
  - SE²(d) = variance of the pooled estimate at dose d
  - τ² = between-study variance (from DL, ML, or REML)
  - t critical value with k-p degrees of freedom (p = number of model parameters)
- **For non-linear models**: SE(d) varies along the curve, computed from `var(ŷ(d)) = x(d)' · Cov(β) · x(d)`
- **Edge case**: When k ≤ p+1, PI is undefined (show warning)

### UI
- Dashed band on DR curve SVG (distinct from solid 95% CI ribbon)
- Toggle: "Show prediction interval" checkbox
- Legend entry: "95% PI (new study range)"

## Section 6: Model-Averaged Predictions (Exceeds dosresmeta)

### What
Weighted average of DR curves from all fitted models, using AIC weights. This accounts for model selection uncertainty — dosresmeta2 does NOT have this.

### Implementation
- **Formula**: `ŷ_avg(d) = Σ w_i · ŷ_i(d)` where w_i = AIC weight of model i
- **SE**: `SE_avg(d) = √(Σ w_i · [SE²_i(d) + (ŷ_i(d) - ŷ_avg(d))²])` (Burnham & Anderson 2002, eq 4.9)
  - Includes both within-model uncertainty and between-model uncertainty
- **CI**: `ŷ_avg(d) ± z_{α/2} · SE_avg(d)`
- **Prediction**: Only average models with w_i > 0.01 (skip negligible contributors)

### UI
- New option in model selector: "Model-Averaged"
- DR curve shows thick line = averaged, thin dashed lines = individual model curves
- Tooltip on curve shows contributing models and their weights
- Model comparison table highlights contributing models

## Section 7: R Validation Suite

### What
Automated cross-validation against R dosresmeta2 using built-in datasets.

### Implementation

#### R Script: `validate_vs_dosresmeta.R`
- Uses dosresmeta2 built-in datasets: `alcohol_cvd`, `cc_ex`, `ci_ex`
- For each dataset, fits: linear, quadratic, RCS-3, RCS-4, FP1, FP2
- Extracts: coefficients, SEs, AIC, BIC, tau², log-likelihood, predicted values at 10 dose points
- Outputs: `validation_reference.json` with all reference values

#### JavaScript Validation: `_testVsR()` inline self-test
- Loads `validation_reference.json` (embedded or fetched)
- Re-fits same models on same data
- Compares: coefficients (tol 1e-4), SEs (tol 1e-3), AIC (tol 0.1), predictions (tol 1e-3)
- Reports: PASS/FAIL per model per dataset

#### Selenium Tests: `TestRValidation` class
- Runs `Rscript validate_vs_dosresmeta.R` to generate reference JSON
- Loads app, injects reference data, triggers comparison
- Asserts all models within tolerance

### Datasets
1. **alcohol_cvd**: Continuous dose (g/day), log-RR outcome, 8 studies, GL covariance
2. **cc_ex**: Case-control, OR outcome, Hamling covariance
3. **ci_ex**: Cohort incidence, RR outcome, GL covariance

## Safety / Constraints

- All matrix inversions check `|det| < 1e-15` → return null
- NaN/Inf checks on all computed values
- Zero-dose handling: `dose + 1e-6` for log/negative powers
- Convergence failure → graceful fallback to two-stage DL with warning
- SE > 0 required for all study-level estimates
- GL covariance: cases must be > 0 and ≤ N for all dose levels
- Hamling: OR must be > 0, CI must bracket OR
- All changes are additive — existing functionality must not break
- Div balance verification after all HTML edits
- No literal `</script>` inside script blocks
