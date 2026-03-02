# validate_vs_dosresmeta.R
# Generate reference values from dosresmeta2 v2.2.0 for cross-validation
# Run: Rscript validate_vs_dosresmeta.R

library(dosresmeta)
library(jsonlite)

results <- list()

# ============================================
# Dataset 1: alcohol_cvd (6 studies, mixed cc/ci)
# ============================================
data(alcohol_cvd)

# --- Linear REML ---
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
  predictions = data.frame(dose=doses, pred=pred_lin$pred,
                           ci_lb=pred_lin$ci.lb, ci_ub=pred_lin$ci.ub)
)

# --- Quadratic REML ---
quad_reml <- dosresmeta(logrr ~ dose + I(dose^2), id=id, type=type, se=se,
                         cases=cases, n=n, data=alcohol_cvd, method="reml")
pred_quad <- predict(quad_reml, newdata=data.frame(dose=doses), expo=FALSE)

results$alcohol_cvd$quadratic_reml <- list(
  coefficients = as.numeric(coef(quad_reml)),
  se = as.numeric(sqrt(diag(vcov(quad_reml)))),
  tau2_matrix = as.matrix(quad_reml$Psi),
  logLik = as.numeric(logLik(quad_reml)),
  AIC = AIC(quad_reml),
  BIC = BIC(quad_reml),
  predictions = data.frame(dose=doses, pred=pred_quad$pred,
                           ci_lb=pred_quad$ci.lb, ci_ub=pred_quad$ci.ub)
)

# --- Linear ML ---
lin_ml <- dosresmeta(logrr ~ dose, id=id, type=type, se=se,
                      cases=cases, n=n, data=alcohol_cvd, method="ml")
results$alcohol_cvd$linear_ml <- list(
  coefficients = as.numeric(coef(lin_ml)),
  se = as.numeric(sqrt(diag(vcov(lin_ml)))),
  tau2 = as.numeric(lin_ml$Psi),
  logLik = as.numeric(logLik(lin_ml)),
  AIC = AIC(lin_ml)
)

# --- Linear Fixed ---
lin_fixed <- dosresmeta(logrr ~ dose, id=id, type=type, se=se,
                         cases=cases, n=n, data=alcohol_cvd, method="fixed")
results$alcohol_cvd$linear_fixed <- list(
  coefficients = as.numeric(coef(lin_fixed)),
  se = as.numeric(sqrt(diag(vcov(lin_fixed))))
)

# --- GL Covariance matrices (Slist) ---
Slist <- lin_reml$Slist
results$alcohol_cvd$Slist <- lapply(Slist, function(s) as.matrix(s))

# --- Study-level data ---
results$alcohol_cvd$data <- alcohol_cvd

# --- RCS-3 REML ---
library(splines)
rcs3_reml <- dosresmeta(logrr ~ rcs(dose, knots=c(10, 25, 50)), id=id, type=type, se=se,
                         cases=cases, n=n, data=alcohol_cvd, method="reml")
pred_rcs3 <- predict(rcs3_reml, newdata=data.frame(dose=doses), expo=FALSE)

results$alcohol_cvd$rcs3_reml <- list(
  coefficients = as.numeric(coef(rcs3_reml)),
  se = as.numeric(sqrt(diag(vcov(rcs3_reml)))),
  tau2_matrix = as.matrix(rcs3_reml$Psi),
  logLik = as.numeric(logLik(rcs3_reml)),
  AIC = AIC(rcs3_reml),
  BIC = BIC(rcs3_reml),
  knots = c(10, 25, 50),
  predictions = data.frame(dose=doses, pred=pred_rcs3$pred,
                           ci_lb=pred_rcs3$ci.lb, ci_ub=pred_rcs3$ci.ub)
)

# --- Quadratic ML ---
quad_ml <- dosresmeta(logrr ~ dose + I(dose^2), id=id, type=type, se=se,
                       cases=cases, n=n, data=alcohol_cvd, method="ml")
results$alcohol_cvd$quadratic_ml <- list(
  coefficients = as.numeric(coef(quad_ml)),
  se = as.numeric(sqrt(diag(vcov(quad_ml)))),
  tau2_matrix = as.matrix(quad_ml$Psi),
  AIC = AIC(quad_ml)
)

# --- Quadratic Fixed ---
quad_fixed <- dosresmeta(logrr ~ dose + I(dose^2), id=id, type=type, se=se,
                          cases=cases, n=n, data=alcohol_cvd, method="fixed")
results$alcohol_cvd$quadratic_fixed <- list(
  coefficients = as.numeric(coef(quad_fixed)),
  se = as.numeric(sqrt(diag(vcov(quad_fixed))))
)

# --- Leave-one-out (Linear REML): drop each study, refit ---
loo_results <- list()
study_ids <- unique(alcohol_cvd$id)
for (sid in study_ids) {
  loo_data <- alcohol_cvd[alcohol_cvd$id != sid, ]
  tryCatch({
    loo_fit <- dosresmeta(logrr ~ dose, id=id, type=type, se=se,
                           cases=cases, n=n, data=loo_data, method="reml")
    loo_results[[as.character(sid)]] <- list(
      coefficients = as.numeric(coef(loo_fit)),
      se = as.numeric(sqrt(diag(vcov(loo_fit)))),
      AIC = AIC(loo_fit)
    )
  }, error = function(e) {
    loo_results[[as.character(sid)]] <<- list(error = e$message)
  })
}
results$alcohol_cvd$leave_one_out <- loo_results

# ============================================
# Dataset 2: coffee_mort (21+ studies)
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
  predictions = data.frame(dose=cof_doses, pred=pred_cof$pred,
                           ci_lb=pred_cof$ci.lb, ci_ub=pred_cof$ci.ub)
)

# --- Coffee quadratic ---
cof_quad <- dosresmeta(logrr ~ dose + I(dose^2), id=id, type=type, se=se,
                        cases=cases, n=n, data=coffee_mort, method="reml")
cof_pred_q <- predict(cof_quad, newdata=data.frame(dose=cof_doses), expo=FALSE)

results$coffee_mort$quadratic_reml <- list(
  coefficients = as.numeric(coef(cof_quad)),
  se = as.numeric(sqrt(diag(vcov(cof_quad)))),
  AIC = AIC(cof_quad),
  predictions = data.frame(dose=cof_doses, pred=cof_pred_q$pred,
                           ci_lb=cof_pred_q$ci.lb, ci_ub=cof_pred_q$ci.ub)
)

# ============================================
# Dataset 3: ci_ex (single study, fixed effects)
# ============================================
data(ci_ex)

ci_cov <- covar.logrr(cases=ci_ex$cases, n=ci_ex$n, y=ci_ex$logrr,
                       v=ci_ex$se^2, type="ci")
results$ci_ex$covariance <- as.matrix(ci_cov)
results$ci_ex$data <- ci_ex

# ============================================
# Dataset 4: cc_ex (single case-control study)
# ============================================
data(cc_ex)

cc_cov <- covar.logrr(cases=cc_ex$case, n=cc_ex$n, y=cc_ex$logrr,
                       v=cc_ex$logrr, type="cc",
                       data=cc_ex)
results$cc_ex$data <- cc_ex

tryCatch({
  results$cc_ex$covariance <- as.matrix(cc_cov)
}, error = function(e) {
  results$cc_ex$covariance_error <<- e$message
})

# ============================================
# Write JSON output
# ============================================
json_output <- toJSON(results, pretty=TRUE, digits=10, auto_unbox=TRUE)
writeLines(json_output, "validation_reference.json")
cat("Validation reference saved to validation_reference.json\n")
cat("Datasets: alcohol_cvd, coffee_mort, ci_ex, cc_ex\n")
cat("Models: linear (REML/ML/fixed), quadratic (REML/ML/fixed), RCS-3 (REML)\n")
cat("Also: leave-one-out, GL covariance matrices\n")
