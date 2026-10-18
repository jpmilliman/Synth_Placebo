# ridge_covariate_validation.R
#
# Does ridge augmentation still work when a time-varying confounder must be adjusted for, and does the
# adjustment still work when the treated unit is outside the donors' range?
# Data: dsc_placebo()'s factor-model generator (copied unchanged), covariate = "confounded": a binary
# covariate that switches on for the treated unit at the treatment date (adding about 8% to its outcome)
# and for each donor at a date that depends on its first factor loading, so the covariate is correlated
# with the factors (unadjusted DSC is biased; TWFE adjustment is biased because the switch date depends on
# the loadings; IFE adjustment is the intended fix).
# Grid: treated exchangeable, selected, edge or outside; stationary or random-walk factors; 14 or 19
# donors; 20 or 40 pre-periods; 12 post periods; 100 replicates per cell.
# Methods: adjustment none / twfe / ife x estimator dsc / ridge; in-space placebo calibration; 90%, 95%.
#
# Written by Claude (Anthropic), at the direction of Jeffrey Milliman.
# Run from the project root: source("simulations/ridge_covariate_validation.R")
# Uses future + furrr::future_pmap (n_workers, 14 by default). Results in simulations/results_ridge_covariate/.
# Functions: simulations/ridge_validation_helpers.R (the prototype of dsc_placebo(augment = "ridge");
# ridge_function_check.R confirms they agree).

suppressPackageStartupMessages({ library(tidyverse); library(future); library(furrr) })
sim_dir <- "simulations"                     # "Placebo simulations/current" for the master copy
out_dir <- file.path(sim_dir, "results_ridge_covariate"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(sim_dir, "ridge_validation_helpers.R"))
if (!exists("n_workers")) n_workers <- 14
if (!exists("n_reps")) n_reps <- 100
T1 <- 12
gamma_level <- 0.08 * exp(5)
seed <- 20261020
levels_conf <- c(0.95, 0.90)
effects <- tibble(effect = c("-10%", "+10%"), pct = c(-0.10, 0.10))

simulate_panel <- function(treated, factors, n_donors, T0, covariate, rep_seed) {    # copied from factor_model_simulation.R
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

run_rep <- function(treated, factors, n_donors, T0, cell_id, rep) {
  d <- simulate_panel(treated, factors, n_donors, T0, "confounded", rep_seed = seed + 100000 * cell_id + rep)
  TT <- T0 + T1; N <- n_donors + 1; post <- seq_len(TT) > T0; h <- T1
  Y <- matrix(d$y, TT, N); Z <- matrix(d$z, TT, N)
  y0 <- Y[post, 1] - gamma_level * Z[post, 1]                 # the treated unit's outcome without the covariate (the truth)
  map_dfr(c("none", "twfe", "ife"), function(adj) {
    core <- cov_core(Y, Z, 1, post, adjust = adj)
    map_dfr(levels_conf, function(cl) imap_dfr(edge_bands(core, cl), function(x, nm) {
      g <- x$gap; gc <- cumsum(g)                            # true effect is zero: the gap is the error (unadjusted gaps include the confounder)
      pw <- map_dfr(seq_len(nrow(effects)), function(k) { s <- effects$pct[k] * y0; ge <- x$gap + s; gce <- cumsum(ge)
        tibble(effect = effects$effect[k], measure = c("per period", "cumulative"),
               power = c(abs(ge[h]) > x$hw_pw[h], abs(gce[h]) > x$hwc_pw[h])) }) |>
        pivot_wider(names_from = effect, values_from = power, names_prefix = "power_")
      tibble(adjust = adj, method = sub("_space", "", nm), conf = cl, measure = c("per period", "cumulative"),
             err_pct = 100 * c(g[h] / y0[h], gc[h] / sum(y0)),
             final_cov_pw = c(abs(g[h]) <= x$hw_pw[h], abs(gc[h]) <= x$hwc_pw[h]),
             hw_pw_pct = 100 * c(x$hw_pw[h] / y0[h], x$hwc_pw[h] / sum(y0))) |> left_join(pw, by = "measure")
    })) |> mutate(gamma_hat = core$gamma, gamma_placebo_sd = core$gamma_placebo_sd, r = core$r, c_ridge = core$c_ridge)
  }) |> mutate(treated, factors, n_donors, T0, rep, gamma_true = gamma_level)
}

grid <- expand_grid(treated = c("exchangeable", "selected", "edge", "outside"), factors = c("stationary", "random_walk"),
                    n_donors = c(14L, 19L), T0 = c(20L, 40L)) |> mutate(cell_id = row_number()) |>
  crossing(rep = seq_len(n_reps)) |> arrange(desc(n_donors), desc(T0))

if (!exists("dry_run") || !isTRUE(dry_run)) {
  plan(multisession, workers = n_workers)
  t0 <- Sys.time()
  res <- future_pmap(grid, possibly(run_rep, otherwise = NULL),
                     .options = furrr_options(seed = TRUE, scheduling = Inf, packages = c("tidyverse", "quadprog", "fect", "withr"))) |> list_rbind()
  plan(sequential)
  run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  saveRDS(list(res = res, run_minutes = run_minutes), file.path(out_dir, "results.rds"))
  summary_tbl <- res |> group_by(treated, factors, n_donors, T0, adjust, method, conf, measure) |>
    summarise(n = n(), bias = mean(err_pct), rmse = sqrt(mean(err_pct^2)), coverage = mean(final_cov_pw),
              finite = mean(is.finite(hw_pw_pct)), median_hw_pct = median(hw_pw_pct[is.finite(hw_pw_pct)]),
              across(starts_with("power_"), mean), gamma_bias_pct = 100 * mean(gamma_hat / gamma_true - 1), .groups = "drop")
  write_csv(summary_tbl, file.path(out_dir, "summary.csv"))
  message("Finished in ", run_minutes, " min; ", nrow(grid) - n_distinct(paste(res$treated, res$factors, res$n_donors, res$T0, res$rep)), " failed replicates.")
}
