# firpo_possebom_comparison.R
#
# How do dsc_placebo()'s intervals relate to the confidence sets of Firpo and
# Possebom (2018, Journal of Causal Inference)? Applied to West Germany
# (scpi::scpi_germany, 16 donors, treatment from 1991) with the same demeaned
# synthetic control estimator throughout, so only the inference differs.
#
# Firpo-Possebom (FP) constant-effect confidence set: test the sharp null
# "the effect is c in every post-treatment period" by subtracting c from the
# treated unit's post-treatment outcomes, recomputing every unit's post/pre
# RMSPE ratio, and ranking the treated unit; the set is {c : p(c) > gamma}.
# FP keep the treated unit in the placebos' donor pools (following Imbens and
# Rubin), so the placebo statistics change with c.
#
#   A  FP: RMSPE ratio, treated unit kept in the placebo donor pools
#   B  RMSPE ratio, treated unit removed from the placebo donor pools
#      (as in dsc_placebo()): placebo statistics do not depend on c
#   C  dsc_placebo() simultaneous band (90%): the constant effects inside the
#      band in every post-treatment period
#
# Run from the project root: source("simulations/firpo_possebom_comparison.R")
# Uses future + furrr::future_map (10 workers).

suppressPackageStartupMessages({ library(tidyverse); library(future); library(furrr) })
source("synthetic_control_functions/synthetic_control_functions_claude.R")
out_dir <- file.path("simulations", "results_firpo_possebom"); dir.create(out_dir, showWarnings = FALSE)

data("scpi_germany", package = "scpi")
ge <- as.data.frame(scpi_germany)[, c("country", "year", "gdp")]
Y <- ge |> pivot_wider(names_from = country, values_from = gdp) |> arrange(year)
years <- Y$year; Y <- as.matrix(Y[, -1]); post <- years >= 1991
treated <- "West Germany"; units <- colnames(Y); gamma <- 0.10

fit_dsc <- function(X, yv) {                        # as inside dsc_placebo()
  sc <- mean(abs(yv[!post]))
  Xp <- X[!post, , drop = FALSE] / sc; yp <- yv[!post] / sc
  w <- CVXR::Variable(ncol(Xp)); b0 <- CVXR::Variable(1)
  prob <- CVXR::Problem(CVXR::Minimize(CVXR::sum_squares(yp - b0 - Xp %*% w)), list(w >= 0, sum(w) == 1))
  CVXR::psolve(prob)
  if (!CVXR::status(prob) %in% c("optimal", "optimal_inaccurate")) try(CVXR::psolve(prob, solver = "CLARABEL", max_iter = 10000), silent = TRUE)
  wv <- as.numeric(CVXR::value(w)); wv[abs(wv) < 1e-8] <- 0
  as.numeric(CVXR::value(b0)) * sc + as.numeric(X %*% wv)
}
ratio <- function(gap) sqrt(mean(gap[post]^2)) / sqrt(mean(gap[!post]^2))

# treated unit's own fit uses pre-treatment data only, so it does not depend on c
gap_T <- Y[, treated] - fit_dsc(Y[, setdiff(units, treated), drop = FALSE], Y[, treated])

# B: placebos without the treated unit (independent of c)
ratio_P_excl <- sapply(setdiff(units, treated), function(j)
  ratio(Y[, j] - fit_dsc(Y[, setdiff(units, c(j, treated)), drop = FALSE], Y[, j])))

# A: placebos with the treated unit in the pool, its post outcomes shifted by -c
grid <- seq(-6, 3, by = 0.02)
plan(multisession, workers = 10)
ratio_P_incl <- future_map(grid, function(cc) {
  Yc <- Y; Yc[post, treated] <- Yc[post, treated] - cc
  sapply(setdiff(units, treated), function(j) ratio(Yc[, j] - fit_dsc(Yc[, setdiff(units, j), drop = FALSE], Yc[, j])))
}, .options = furrr_options(seed = TRUE, globals = c("Y", "post", "treated", "units", "fit_dsc", "ratio")), .progress = TRUE)
plan(sequential)

res <- tibble(c = grid,
              ratio_T = map_dbl(grid, ~ ratio(gap_T - .x * post)),
              p_A = map2_dbl(ratio_P_incl, ratio_T, ~ (1 + sum(.x >= .y)) / length(units)),
              p_B = map_dbl(ratio_T, ~ (1 + sum(ratio_P_excl >= .x)) / length(units)),
              n_placebos_changed = map_int(ratio_P_incl, ~ sum(abs(.x - ratio_P_incl[[which.min(abs(grid))]]) > 1e-6)))
set_of <- function(keep) { x <- res$c[keep]; if (!length(x)) c(NA, NA) else range(x) }
fit <- dsc_placebo(ge, "country", "year", "gdp", treated, treatment_period = c(year = 1991), conf_level = 1 - gamma, diagnostics = FALSE)
e <- fit$effects[fit$effects$post, ]
summary_tbl <- tibble(
  method = c("A: Firpo-Possebom (ratio test, treated in placebo pools)", "B: ratio test, treated removed from placebo pools",
             "C: dsc_placebo() simultaneous band, constant effects inside it every year"),
  lower = c(set_of(res$p_A > gamma)[1], set_of(res$p_B > gamma)[1], max(e$eff_low_sim)),
  upper = c(set_of(res$p_A > gamma)[2], set_of(res$p_B > gamma)[2], min(e$eff_high_sim)),
  connected = c(all(diff(which(res$p_A > gamma)) == 1), all(diff(which(res$p_B > gamma)) == 1), NA))
extra <- list(mean_effect = mean(e$effect), p_at_zero = c(A = res$p_A[which.min(abs(res$c))], B = res$p_B[which.min(abs(res$c))]),
              max_placebos_changed = max(res$n_placebos_changed),
              per_year_band = tibble(year = e$year, effect = e$effect, lo_sim = e$eff_low_sim, hi_sim = e$eff_high_sim, lo_pw = e$eff_low_pw, hi_pw = e$eff_high_pw))
saveRDS(list(grid = res, summary = summary_tbl, extra = extra), file.path(out_dir, "results.rds"))
write_csv(summary_tbl, file.path(out_dir, "summary.csv"))
print(summary_tbl)
