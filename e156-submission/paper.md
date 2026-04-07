Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

MetaSprint Dose-Response: Browser-Based Nonlinear Dose-Response Meta-Analysis with R Cross-Validation

Can a browser-based application implement dose-response meta-analysis with nonlinear models validated against established R packages without requiring installation? MetaSprint Dose-Response is a single HTML file of 25,800 lines implementing eight dose-response models including linear, quadratic, Emax, fractional polynomial, and restricted cubic splines with Greenland-Longnecker covariance reconstruction for correlated dose contrasts. The engine provides ML and REML estimation via profile likelihood with golden-section optimization, one-stage mixed-effects modeling, AIC-weighted model averaging across fitted curves, and Bayesian Laplace approximation with credible interval bands. Cross-validation against R dosresmeta2 v2.2.0 confirmed agreement within tolerance of 1e-4 for coefficients, standard errors, tau-squared, AIC, and predictions across linear, quadratic, and spline models. Leave-one-out sensitivity analysis across six studies and bootstrap breakpoint confidence intervals with 500 cluster resamples provide robustness assessment for dose-finding decisions. The platform delivers publication-quality dose-response curves directly from a structured 40-day sprint workflow. However, a limitation is that the current implementation supports only two-level study-dose clustering without three-level hierarchical extensions.

Outside Notes

Type: methods
Primary estimand: Dose-response coefficients
App: MetaSprint Dose-Response v1.0
Data: R dosresmeta2 v2.2.0 reference datasets (alcohol_cvd, coffee_mort)
Code: https://github.com/mahmood726-cyber/metasprint-dose-response
Version: 1.0
Validation: DRAFT

References

1. Crippa A, Orsini N. Dose-response meta-analysis of differences in means. BMC Med Res Methodol. 2016;16:91.
2. Greenland S, Longnecker MP. Methods for trend estimation from summarized dose-response data, with applications to meta-analysis. Am J Epidemiol. 1992;135(11):1301-1309.
3. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
