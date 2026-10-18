# ridge_function_check.R
#
# Does dsc_placebo(augment = "ridge") reproduce the prototype that the ridge simulations validated
# (ridge_validation_helpers.R)? The two differ only in the weight solver (CVXR in dsc_placebo(),
# quadprog in the prototype), so estimates and intervals should agree up to solver tolerance.
# Compared at the final post period: the effect, the pointwise interval, the simultaneous band and the
# cumulative versions, and the cross-validated penalty.
#   A. Real panels: Kansas (quarterly GDP per capita, 2012 Q2), California (cigarette sales, 1989),
#      West Germany (GDP, 1991).
#   B. Simulated panels (dsc_placebo()'s factor-model generator, copied unchanged): treated unit
#      exchangeable, edge or outside; stationary or random-walk factors; 14 donors; 20 or 40 pre-periods;
#      3 replicates per cell (36 panels).
#   C. With covariate adjustment: the confounded covariate, TWFE and IFE (2 factors, fixed so that both
#      implementations fit the same model); treated unit exchangeable or outside, random-walk factors,
#      14 donors, 40 pre-periods, 3 replicates (12 panels x 2 adjustments).
# Differences are reported in % of the treated unit's untreated outcome (cumulative: of its sum).
#
# Written by Claude (Anthropic), at the direction of Jeffrey Milliman.
# Run from the project root: source("simulations/ridge_function_check.R")
# Uses future + furrr::future_pmap (n_workers, 14 by default). Results in simulations/results_ridge_function_check/.

suppressPackageStartupMessages({ library(tidyverse); library(future); library(furrr) })
sim_dir <- "simulations"                     # "Placebo simulations/current" for the master copy
out_dir <- file.path(sim_dir, "results_ridge_function_check"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
functions_src <- c("synthetic_control_functions/synthetic_control_functions_claude.R", file.path(sim_dir, "ridge_validation_helpers.R"))
for (f in functions_src) source(f)
if (!exists("n_workers")) n_workers <- 14
T1 <- 12
gamma_level <- 0.08 * exp(5)
seed <- 20261021
conf <- 0.90

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

# final-period comparison of a dsc_placebo() fit with a prototype band
compare <- function(fit, b, c_proto, y0) {
  e <- fit$effects[fit$effects$post, ]; h <- nrow(e)
  fn <- c(effect = e$effect[h], low_pw = e$eff_low_pw[h], high_pw = e$eff_high_pw[h], low_sim = e$eff_low_sim[h], high_sim = e$eff_high_sim[h],
          cum_effect = e$cum_effect[h], cum_low_pw = e$cum_eff_low_pw[h], cum_high_pw = e$cum_eff_high_pw[h])
  pr <- c(effect = b$gap[h], low_pw = b$gap[h] - b$hw_pw[h], high_pw = b$gap[h] + b$hw_pw[h], low_sim = b$gap[h] - b$hw_sim[h],
          high_sim = b$gap[h] + b$hw_sim[h], cum_effect = b$gapc[h], cum_low_pw = b$gapc[h] - b$hwc_pw[h], cum_high_pw = b$gapc[h] + b$hwc_pw[h])
  den <- ifelse(grepl("^cum", names(fn)), sum(y0), y0[h])
  tibble(quantity = names(fn), dsc_placebo = unname(fn), prototype = unname(pr), diff_pct = 100 * (unname(fn) - unname(pr)) / den,
         c_function = fit$ridge$c, c_prototype = c_proto)
}

# ---- A. Real panels ----------------------------------------------------------------------------------
kansas <- as.data.frame(augsynth::kansas)
data("smoking", package = "tidysynth"); data("scpi_germany", package = "scpi")
cases <- list(
  kansas = list(d = kansas, unit = "state", time = c("year", "qtr"), y = "gdpcapita", tr = "Kansas", tp = c(year = 2012, qtr = 2)),
  california = list(d = as.data.frame(smoking), unit = "state", time = "year", y = "cigsale", tr = "California", tp = c(year = 1989)),
  germany = list(d = as.data.frame(scpi_germany), unit = "country", time = "year", y = "gdp", tr = "West Germany", tp = c(year = 1991)))
real_check <- function(case) {
  cs <- cases[[case]]
  fit <- suppressWarnings(dsc_placebo(cs$d, cs$unit, cs$time, cs$y, cs$tr, cs$tp, conf_level = conf, augment = "ridge", diagnostics = FALSE))
  Y <- cbind(fit$effects$observed, as.matrix(fit$placebo$gaps) + as.matrix(fit$placebo$synthetic))   # treated, then the donors' outcomes
  post <- fit$effects$post
  core <- edge_core(unname(Y), 1, post, estimators = "ridge"); b <- edge_bands(core, conf)$ridge_space
  compare(fit, b, core$c_ridge, fit$effects$observed[post]) |> mutate(set = "real", panel = case, .before = 1)
}

# ---- B and C. Simulated panels -------------------------------------------------------------------------
sim_check <- function(treated, factors, n_donors, T0, adjust, cell_id, rep) {
  cov <- if (adjust == "none") "none" else "confounded"
  d <- simulate_panel(treated, factors, n_donors, T0, cov, rep_seed = seed + 1000 * cell_id + rep)
  TT <- T0 + T1; post <- seq_len(TT) > T0; Y <- matrix(d$y, TT, n_donors + 1); Z <- matrix(d$z, TT, n_donors + 1)
  args <- list(data = as.data.frame(d), unit_col = "unit", time_cols = "period", outcome_col = "y", treated_unit = "u01",
               treatment_period = c(period = T0 + 1), conf_level = conf, augment = "ridge", diagnostics = FALSE)
  if (adjust != "none") args <- c(args, list(adjust_covariates = "z", adjust_model = adjust, ife_r = 2))
  fit <- suppressWarnings(do.call(dsc_placebo, args))
  core <- if (adjust == "none") edge_core(Y, 1, post, estimators = "ridge") else
    cov_core(Y, Z, 1, post, adjust = adjust, estimators = "ridge", ife_r = 2, ife_seed = 12112025)
  b <- edge_bands(core, conf)$ridge_space
  compare(fit, b, core$c_ridge, Y[post, 1]) |>
    mutate(set = if (adjust == "none") "simulated" else paste("simulated, adjust", adjust), panel = paste(treated, factors, T0, rep),
           gamma_function = if (adjust == "none") NA else fit$adjustment$gamma, gamma_prototype = if (adjust == "none") NA else core$gamma, .before = 1)
}
grid <- bind_rows(
  expand_grid(treated = c("exchangeable", "edge", "outside"), factors = c("stationary", "random_walk"), n_donors = 14L, T0 = c(20L, 40L), adjust = "none"),
  expand_grid(treated = c("exchangeable", "outside"), factors = "random_walk", n_donors = 14L, T0 = 40L, adjust = c("twfe", "ife"))) |>
  mutate(cell_id = row_number()) |> crossing(rep = 1:3)

if (!exists("dry_run") || !isTRUE(dry_run)) {
  t0 <- Sys.time()
  real <- map_dfr(names(cases), real_check)
  plan(multisession, workers = n_workers)
  sim <- future_pmap(grid, possibly(sim_check, otherwise = NULL),
                     .options = furrr_options(seed = TRUE, scheduling = Inf, packages = c("tidyverse", "quadprog", "CVXR", "fect", "withr"))) |> list_rbind()
  plan(sequential)
  run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  res <- bind_rows(real, sim)
  summary_tbl <- res |> group_by(set, quantity) |>
    summarise(n = n(), max_abs_diff_pct = max(abs(diff_pct)), median_abs_diff_pct = median(abs(diff_pct)),
              same_penalty = mean(c_function == c_prototype), .groups = "drop")
  saveRDS(list(results = res, summary = summary_tbl, run_minutes = run_minutes), file.path(out_dir, "results.rds"))
  write_csv(summary_tbl, file.path(out_dir, "summary.csv"))
  print(summary_tbl, n = Inf)
  message("Finished in ", run_minutes, " min; ", nrow(grid) - n_distinct(sim$panel, sim$set), " failed simulated panels.")
}
