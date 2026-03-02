# Top 5 Advanced Features Design

**Goal:** Add dose-finding metrics, Bayesian DR, DR publication bias, dose-specific forest plot, and breakpoint bootstrap CI.

## Feature 1: Dose-Finding & Safe-Dose Estimation
- computeDoseMetrics(model, maxDose, clinicalThreshold): bisection root-finding
- MED, ED50, ED90, therapeutic window, NOAEL
- Annotated markers on SVG + summary card

## Feature 2: Bayesian DR (Laplace Approximation)
- fitBayesianDR(points, modelType, priorSD): ML mode + inverse Hessian
- Credible interval bands on SVG (darker than CI)
- Prior: N(0, 10^2) on coefficients; posterior summary table

## Feature 3: DR-Specific Publication Bias
- renderDRFunnelPlot(points, model): dose on X, 1/SE on Y, residual coloring
- testDoseReportingBias(points): chi-squared on dose distribution uniformity
- Dose-arm Egger-like regression on residuals

## Feature 4: Dose-Specific Forest Plot
- renderDoseForestPlot(points, confLevel): rows = study x dose arm
- Grouped by dose bin (5-8 groups), pooled diamond per bin
- Container: #doseForestContainer

## Feature 5: Breakpoint Bootstrap CI
- bootstrapBreakpoints(points, modelType, nBoot, confLevel): 500 resamples
- 95% CI for plateau and inflection doses
- Shaded bands on SVG breakpoint markers
