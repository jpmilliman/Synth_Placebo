# ridge_edge_validation.R
#
# Can a different estimator (ridge augmentation) or a different calibration (in-time placebos) restore
# coverage for a treated unit at or beyond the edge of the donors' range? Uses the data generator of
# Placebo simulations/current/factor_model_simulation.R (copied unchanged): log-scale factor model,
# 2 factors (stationary AR 0.8 or random walk), AR(1) errors; treated unit exchangeable (a check),
# at the edge of the donors' loadings, or 1 SD outside them; 19 or 39 donors; 40 or 80 pre-periods
# (in-time calibration needs a long pre-period); 12 post periods; 100 replicates per cell.
# Methods (ridge_validation_helpers.R): estimator dsc or ridge x calibration space, time or both; 95% and 90%.
# All six are shift-invariant (no calibration uses the treated unit's post period), so coverage is
# measured at no effect and power exactly for effects of -10%, +10% and a ramp to -10%.
#
# Written by Claude (Anthropic), at the direction of Jeffrey Milliman.
# Run from the project root: source("simulations/ridge_edge_validation.R")
# Uses future + furrr::future_pmap (n_workers, 14 by default). Results in simulations/results_ridge_edge/.
# Functions: simulations/ridge_validation_helpers.R (the prototype of dsc_placebo(augment = "ridge");
# ridge_function_check.R confirms they agree).

suppressPackageStartupMessages({ library(tidyverse); library(future); library(furrr) })
sim_dir <- "simulations"                     # "Placebo simulations/current" for the master copy
out_dir <- file.path(sim_dir, "results_ridge_edge"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(sim_dir, "ridge_validation_helpers.R"))
if (!exists("n_workers")) n_workers <- 14
n_reps <- 100
T1 <- 12
gamma_level <- 0.08 * exp(5)
seed <- 20261018
levels_conf <- c(0.95, 0.90)
effects <- tibble(effect = c("-10%", "+10%", "ramp to -10%"), pct = c(-0.10, 0.10, -0.10), shape = c("constant", "constant", "ramp"))

# ---- Data generator: copied from Placebo simulations/current/factor_model_simulation.R --------------
simulate_panel <- function(treated, factors, n_donors, T0, covariate, rep_seed) {
  set.seed(rep_seed)
  TT <- T0 + T1; N <- n_donors + 1; K <- 2
  f <- sapply(seq_len(K), function(k) {
    u <- if (factors == "stationary") rnorm(TT, 0, 0.05) else rnorm(TT, 0.005, 0.03)
    if (factors == "stationary") as.numeric(stats::filter(u, 0.8, method = "recursive")) else cumsum(u)
  })
  lam <- matrix(rnorm(N * K), N, K)
  if (treated == "selected") { top <- which.max(lam[, 1]); lam[c(1, top), ] <- lam[c(top, 1), ] }
  if (treated == "edge")    lam[1, ] <- apply(lam[-1, , drop = FALSE], 2, max)
  if (treated == "outside") lam[1, ] <- apply(lam[-1, , drop = FALSE], 2, max) + 1
  sd_e <- c(if (treated == "noisy") 0.04 else 0.02, rep(0.02, n_donors))
  e <- sapply(seq_len(N), function(i) as.numeric(stats::filter(rnorm(TT, 0, sd_e[i] * sqrt(1 - 0.5^2)), 0.5, method = "recursive")))
  z <- matrix(0, TT, N)
  if (covariate != "none") {
    start <- if (covariate == "exogenous") sample(3:TT, N, replace = TRUE) else
      pmin(pmax(round(TT * (0.6 - 0.25 * lam[, 1] + rnorm(N, 0, 0.1))), 3), TT)
    start[1] <- T0 + 1
    z <- sapply(seq_len(N), function(i) as.numeric(seq_len(TT) >= start[i]))
  }
  alpha <- rnorm(N, 0, 0.3)
  logy <- sapply(seq_len(N), function(i) 5 + alpha[i] + 0.01 * seq_len(TT) + f %*% lam[i, ] + e[, i])
  tibble(unit = rep(sprintf("u%02d", seq_len(N)), each = TT), period = rep(seq_len(TT), N),
         y = exp(as.vector(logy)) + gamma_level * as.vector(z), z = as.vector(z))
}

shift_path <- function(pct, shape, y0) pct * (if (shape == "constant") 1 else seq_along(y0) / length(y0)) * y0

run_rep <- function(treated, factors, n_donors, T0, cell_id, rep) {
  d <- simulate_panel(treated, factors, n_donors, T0, "none", rep_seed = seed + 100000 * cell_id + rep)
  TT <- T0 + T1; N <- n_donors + 1; Y <- matrix(d$y, TT, N); post <- seq_len(TT) > T0; h <- T1
  y0 <- Y[post, 1]
  core <- edge_core(Y, 1, post)
  map_dfr(levels_conf, function(cl) {
    b <- edge_bands(core, cl)
    imap_dfr(b, function(x, nm) {
      base <- tibble(measure = c("per period", "cumulative"),
                     err_pct = 100 * c(x$gap[h] / y0[h], x$gapc[h] / sum(y0)),
                     final_cov_pw = c(abs(x$gap[h]) <= x$hw_pw[h], abs(x$gapc[h]) <= x$hwc_pw[h]),
                     path_cov_sim = c(all(abs(x$gap) <= x$hw_sim), all(abs(x$gapc) <= x$hwc_sim)),
                     hw_pw_pct = 100 * c(x$hw_pw[h] / y0[h], x$hwc_pw[h] / sum(y0)))
      pw <- pmap_dfr(effects, function(effect, pct, shape) { s <- shift_path(pct, shape, y0); g <- x$gap + s; gc <- cumsum(g)
        tibble(measure = c("per period", "cumulative"), effect,
               power_final_pw = c(abs(g[h]) > x$hw_pw[h], abs(gc[h]) > x$hwc_pw[h]),
               power_any_sim = c(any(abs(g) > x$hw_sim), any(abs(gc) > x$hwc_sim))) })
      left_join(base, pw, by = "measure", relationship = "one-to-many") |> mutate(method = nm)
    }) |> mutate(conf = cl)
  }) |> mutate(treated, factors, n_donors, T0, rep, c_ridge = core$c_ridge, n_windows = core$n_windows)
}

grid <- expand_grid(treated = c("exchangeable", "edge", "outside"), factors = c("stationary", "random_walk"),
                    n_donors = c(19L, 39L), T0 = c(40L, 80L)) |> mutate(cell_id = row_number()) |>
  crossing(rep = seq_len(n_reps)) |> arrange(desc(n_donors), desc(T0))

if (!exists("dry_run") || !isTRUE(dry_run)) {
  plan(multisession, workers = n_workers)
  t0 <- Sys.time()
  res <- future_pmap(grid, possibly(run_rep, otherwise = NULL),
                     .options = furrr_options(seed = TRUE, scheduling = Inf, packages = c("tidyverse", "quadprog"))) |> list_rbind()
  plan(sequential)
  run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  saveRDS(res, file.path(out_dir, "results.rds"))
  per_rep <- res |> distinct(treated, factors, n_donors, T0, conf, method, measure, rep, err_pct, final_cov_pw, path_cov_sim, hw_pw_pct)
  coverage <- per_rep |> group_by(treated, factors, n_donors, T0, conf, method, measure) |>
    summarise(n = n(), bias = mean(err_pct), rmse = sqrt(mean(err_pct^2)), final_coverage_pw = mean(final_cov_pw),
              path_coverage_sim = mean(path_cov_sim), finite = mean(is.finite(hw_pw_pct)),
              median_hw_pw_pct = median(hw_pw_pct[is.finite(hw_pw_pct)]), .groups = "drop")
  power <- res |> group_by(treated, factors, n_donors, T0, conf, method, measure, effect) |>
    summarise(power_final_pw = mean(power_final_pw), power_any_sim = mean(power_any_sim), .groups = "drop")
  write_csv(coverage, file.path(out_dir, "summary.csv"))
  write_csv(power, file.path(out_dir, "power.csv"))
  message("Finished in ", run_minutes, " min; ", nrow(grid) - n_distinct(paste(res$treated, res$factors, res$n_donors, res$T0, res$rep)), " failed replicates.")
}
