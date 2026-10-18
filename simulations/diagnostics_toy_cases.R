##############################################################################
# diagnostics_toy_cases.R
#
# Helpers for the diagnostics write-up (diagnostics_toy_cases.qmd):
#   toy_panel()            the factor-model generator of augmented_comparison_simulation.R,
#                          also returning the TRUE factor loadings
#   run_toy_case()         simulate a panel and fit dsc_placebo() (RMSPE scale, diagnostics;
#                          augment = "ridge" for the ridge-augmented estimator)
#   plot_true_loadings()   donors' true loadings, their convex hull, and the treated unit
#   diagnostics_table()    the four diagnostics and the assessment, as a one-row tibble
# Requires synthetic_control_functions_claude.R to be sourced first.
##############################################################################
library(tidyverse)

toy_panel <- function(treated, factors, n_donors = 40, T0 = 40, T1 = 12, rep_seed) {
  set.seed(rep_seed)
  TT <- T0 + T1; N <- n_donors + 1; K <- 2
  f <- sapply(seq_len(K), function(k) {
    u <- if (factors == "stationary") rnorm(TT, 0, 0.05) else rnorm(TT, 0.005, 0.03)
    if (factors == "stationary") as.numeric(stats::filter(u, 0.8, method = "recursive")) else cumsum(u)
  })
  lam <- matrix(rnorm(N * K), N, K)
  if (treated == "edge")    lam[1, ] <- apply(lam[-1, , drop = FALSE], 2, max)
  if (treated == "outside") lam[1, ] <- apply(lam[-1, , drop = FALSE], 2, max) + 1
  e <- sapply(seq_len(N), function(i) as.numeric(stats::filter(rnorm(TT, 0, 0.02 * sqrt(1 - 0.5^2)), 0.5, method = "recursive")))
  alpha <- rnorm(N, 0, 0.3)
  logy <- sapply(seq_len(N), function(i) 5 + alpha[i] + 0.01 * seq_len(TT) + f %*% lam[i, ] + e[, i])
  units <- sprintf("u%02d", seq_len(N))
  list(data = tibble(unit = rep(units, each = TT), period = rep(seq_len(TT), N), y = exp(as.vector(logy))),
       loadings = tibble(unit = units, lambda1 = lam[, 1], lambda2 = lam[, 2], treated = units == "u01"),
       T0 = T0)
}

run_toy_case <- function(treated, factors, rep_seed, T0 = 40, augment = "none") {
  p <- toy_panel(treated, factors, T0 = T0, rep_seed = rep_seed)
  fit <- dsc_placebo(as.data.frame(p$data), "unit", "period", "y", "u01", treatment_period = c(period = T0 + 1),
                     scale = "rmspe", warn = FALSE, augment = augment)
  list(treated = treated, factors = factors, seed = rep_seed, panel = p, fit = fit)
}

plot_true_loadings <- function(case) {
  donors <- case$panel$loadings |> filter(!treated)
  hull <- donors |> slice(grDevices::chull(lambda1, lambda2))
  ggplot() +
    geom_polygon(data = hull, aes(lambda1, lambda2), fill = "grey90", colour = "grey60") +
    geom_point(data = donors, aes(lambda1, lambda2), colour = "grey35", size = 1.8) +
    geom_point(data = filter(case$panel$loadings, treated), aes(lambda1, lambda2), shape = 17, size = 4, colour = "firebrick") +
    coord_equal() + theme_minimal(base_size = 11) +
    labs(title = "True factor loadings", subtitle = "Grey: donors and their convex hull; red triangle: treated unit",
         x = "Loading on factor 1", y = "Loading on factor 2")
}

diagnostics_table <- function(case) {
  dg <- case$fit$diagnostics
  e <- case$fit$effects |> filter(post)
  tibble(case = paste(case$factors, case$treated),
         `fit: placebos worse` = scales::percent(dg$fit_rank$share_placebos_worse, 1),
         `fit flag` = dg$fit_rank$flag,
         factors_est = dg$hull$n_factors,
         `hull: donors further out` = if (is.null(dg$hull$share_donors_beyond)) NA_character_ else scales::percent(dg$hull$share_donors_beyond, 1),
         `hull flag` = dg$hull$flag,
         `implied drift (%)` = round(dg$extrapolation$max_post_drift_pct %||% NA_real_, 1),
         `factor AR(1)` = if (is.null(dg$persistence$factor_ar1)) NA_real_ else round(max(dg$persistence$factor_ar1), 2),
         trending = dg$persistence$trending,
         `cumulative error (%)` = round(100 * last(e$cum_effect) / sum(e$observed), 1),
         `interval covers 0` = last(e$cum_eff_low_pw) <= 0 & last(e$cum_eff_high_pw) >= 0,
         assessment = dg$assessment)
}

