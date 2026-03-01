# Dose-Response Enhancements Design

**Date**: 2026-03-01
**App**: `metasprint-dose-response.html` (19,344 lines)
**Approach**: Enhance existing app in-place with 4 major features

## Context

The existing app has Linear, Quadratic, and Emax dose-response models with two input modes (one-stage multi-point, two-stage slope) and DL/HKSJ/FE pooling. Four enhancements are needed to reach dosresmeta2-level capability.

## Enhancement 1: RCS (Restricted Cubic Splines)

### What
Add restricted cubic spline models (3-knot and 4-knot) to the model selector.

### Implementation
- **Knot placement**: Percentile-based defaults (3 knots: 10th/50th/90th; 4 knots: 5th/35th/65th/95th)
- **Basis function**: Harrell (2001) truncated power basis with boundary constraints
- **Fitting**: Weighted least squares on the RCS design matrix
- **Non-linearity test**: Wald test H0 = all spline coefficients = 0
- **Output**: Coefficients, SEs, p-values, AIC/BIC, R-squared, non-linearity p-value
- **Reference**: The autopilot has `computeRCSBasis()` at line ~12454 that can be ported

### UI
- Two new entries in model selector dropdown: "RCS (3 knots)" and "RCS (4 knots)"
- Collapsible knot configuration panel (auto vs manual knot positions)
- Dose-response curve rendered with smooth spline + 95% CI ribbon

## Enhancement 2: Fractional Polynomial Models

### What
First-degree (FP1) and second-degree (FP2) fractional polynomial models.

### Implementation
- **Power set**: {-2, -1, -0.5, 0, 0.5, 1, 2, 3} where 0 = log(dose)
- **FP1**: y = b0 + b1 * dose^p (8 candidate models)
- **FP2**: y = b0 + b1 * dose^p1 + b2 * dose^p2 (36 candidate combinations)
  - When p1 = p2: y = b0 + b1 * dose^p + b2 * dose^p * log(dose)
- **Selection**: Best deviance/AIC among all candidates
- **Fitting**: Closed-form WLS for each power combination
- **Comparison**: Deviance tests (FP2 vs FP1 vs Linear)

### UI
- Two new entries in model selector: "FP1 (best)" and "FP2 (best)"
- Model comparison table shows best FP1 and FP2 powers alongside other models
- Display selected powers in results summary

## Enhancement 3: One-Stage Mixed-Effects

### What
Random-effects dose-response model fitting all studies simultaneously.

### Implementation
- **Fixed effects**: Dose terms from selected model (linear/RCS/polynomial/FP)
- **Random effects**: Per-study random intercepts + random slopes (unstructured covariance)
- **Estimation**: Iterative GLS with REML variance component estimation
  - Initialize with two-stage estimates
  - Iterate: update fixed effects (GLS), update variance components (REML), check convergence
  - Max 50 iterations, tolerance 1e-6
- **Output**: Pooled curve + per-study curves, variance components, likelihood-based CIs
- **Fallback**: If convergence fails (singular Hessian, >50 iterations), fall back to two-stage with warning

### UI
- Toggle in Analyze phase: "Two-stage (default)" vs "One-stage mixed-effects"
- When one-stage selected, show per-study curves overlaid on pooled curve
- Convergence indicator (iterations, final tolerance)

## Enhancement 4: Continuous Input Mode

### What
Mean/SD/N per dose level as a third input mode for continuous outcomes.

### Implementation
- **Input fields**: Dose, N, Mean, SD per dose level per study
- **Effect computation**:
  - MD: mean_i - mean_ref, SE = sqrt(SD_i^2/N_i + SD_ref^2/N_ref)
  - SMD (Hedges' g): with pooled SD and small-sample correction J
- **Covariance**: Within-study contrasts share the reference group, so covariance = var(ref)
- **Integration**: Computed effects feed into all existing models

### UI
- New radio button in Extract phase: "Continuous (Mean/SD/N)"
- When selected, table columns change to: Study, Dose, N, Mean, SD, Unit
- Auto-computes effect sizes when user enters data (live preview)

## Model Selector (Final)

Dropdown options after all enhancements:
1. Auto-select (best AIC)
2. Linear
3. Quadratic
4. Emax
5. RCS (3 knots)
6. RCS (4 knots)
7. FP1 (best first-degree)
8. FP2 (best second-degree)

## Model Comparison Table (Final)

| Model | AIC | BIC | R-squared | Non-linearity p | Parameters |
|-------|-----|-----|-----------|-----------------|------------|
| Linear | ... | ... | ... | — | b1 |
| Quadratic | ... | ... | ... | p(b2) | b1, b2 |
| Emax | ... | ... | ... | — | Emax, ED50 |
| RCS-3 | ... | ... | ... | p(spline) | k1, k2 |
| RCS-4 | ... | ... | ... | p(spline) | k1, k2, k3 |
| FP1 (p=X) | ... | ... | ... | — | b1, power |
| FP2 (p1=X, p2=Y) | ... | ... | ... | p(FP2 vs FP1) | b1, b2, powers |

Best model highlighted. Auto-select uses AIC with parsimony preference (linear wins ties).

## Safety Guards

- All matrix inversions check for singularity (|pivot| < 1e-15 -> null)
- NaN/Inf checks on all computed values
- Zero-dose handling: add small epsilon (1e-6) for log/negative powers
- Convergence failure -> graceful fallback with user warning
- SE > 0 required for all study-level estimates

## Testing Strategy

- Unit tests for each model against R dosresmeta2 reference values
- Edge cases: k=1, k=2, all-same-dose, zero-dose-only, extreme heterogeneity
- Div balance verification after all HTML edits
- localStorage key uniqueness check
