# ridge_small_sample_validation.R
#
# Checks of the ridge-augmented demeaned synthetic control (ridge_validation_helpers.R), the basis of
# dsc_placebo(augment = "ridge"):
#   A. Small pools and short pre-periods: dsc_placebo()'s factor-model generator (copied unchanged), treated
#      unit exchangeable, selected, noisy, edge or outside; stationary or random-walk factors; 9, 14 or 19
#      donors; 10 or 20 pre-periods; 12 post periods; 100 replicates per cell. DSC vs ridge, in-space
#      placebo calibration, 90% and 95%.
#   B. Real panels: every unit of Kansas (quarterly GDP), California's tobacco panel (cigarette sales) and
#      West Germany's OECD panel (GDP) taken in turn as a fake treated unit at the real treatment date
#      (the real treated unit dropped); coverage of zero, widths and bias, DSC vs ridge.
#   C. The real treated units (Kansas, California, West Germany): estimates and intervals, DSC vs ridge.
#   D. Agreement with augsynth: augsynth(progfunc = "Ridge", fixedeff = TRUE) on the three real panels, with
#      its own cross-validated lambda and with the prototype's lambda, vs the prototype's ridge estimate.
#
# Written by Claude (Anthropic), at the direction of Jeffrey Milliman.
# Run from the project root: source("simulations/ridge_small_sample_validation.R")
# Uses future + furrr::future_pmap (n_workers, 14 by default). Results in simulations/results_ridge_small_sample/.
# Functions: simulations/ridge_validation_helpers.R (the prototype of dsc_placebo(augment = "ridge");
# ridge_function_check.R confirms they agree).

suppressPackageStartupMessages({ library(tidyverse); library(future); library(furrr) })
sim_dir <- "simulations"                     # "Placebo simulations/current" for the master copy
out_dir <- file.path(sim_dir, "results_ridge_small_sample"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(sim_dir, "ridge_validation_helpers.R"))
if (!exists("n_workers")) n_workers <- 14
n_reps <- 100
T1 <- 12
gamma_level <- 0.08 * exp(5)
seed <- 20261019
levels_conf <- c(0.95, 0.90)
effects <- tibble(effect = c("-10%", "+10%"), pct = c(-0.10, 0.10))

# ---- A. Simulations ----------------------------------------------------------------------------------
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

# coverage, width, bias and power for one panel (treated column I), both estimators, in-space calibration
evaluate_panel <- function(Y, I, post) {
  h <- sum(post); y0 <- Y[post, I]
  core <- edge_core(Y, I, post)
  map_dfr(levels_conf, function(cl) {
    b <- edge_bands(core, cl)[c("dsc_space", "ridge_space")]
    imap_dfr(b, function(x, nm) {
      pw <- map_dfr(seq_len(nrow(effects)), function(k) { s <- effects$pct[k] * y0; g <- x$gap + s; gc <- cumsum(g)
        tibble(effect = effects$effect[k], measure = c("per period", "cumulative"),
               power_final_pw = c(abs(g[h]) > x$hw_pw[h], abs(gc[h]) > x$hwc_pw[h]),
               power_any_sim = c(any(abs(g) > x$hw_sim), any(abs(gc) > x$hwc_sim))) }) |>
        pivot_wider(names_from = effect, values_from = c(power_final_pw, power_any_sim))
      tibble(method = sub("_space", "", nm), conf = cl, measure = c("per period", "cumulative"),
             err_pct = 100 * c(x$gap[h] / y0[h], x$gapc[h] / sum(y0)),
             final_cov_pw = c(abs(x$gap[h]) <= x$hw_pw[h], abs(x$gapc[h]) <= x$hwc_pw[h]),
             path_cov_sim = c(all(abs(x$gap) <= x$hw_sim), all(abs(x$gapc) <= x$hwc_sim)),
             hw_pw_pct = 100 * c(x$hw_pw[h] / y0[h], x$hwc_pw[h] / sum(y0))) |> left_join(pw, by = "measure")
    })
  }) |> mutate(c_ridge = core$c_ridge)
}

run_sim_rep <- function(treated, factors, n_donors, T0, cell_id, rep) {
  d <- simulate_panel(treated, factors, n_donors, T0, "none", rep_seed = seed + 100000 * cell_id + rep)
  Y <- matrix(d$y, T0 + T1, n_donors + 1)
  evaluate_panel(Y, 1, seq_len(T0 + T1) > T0) |> mutate(treated, factors, n_donors, T0, rep)
}
sim_grid <- expand_grid(treated = c("exchangeable", "selected", "noisy", "edge", "outside"), factors = c("stationary", "random_walk"),
                        n_donors = c(9L, 14L, 19L), T0 = c(10L, 20L)) |> mutate(cell_id = row_number()) |> crossing(rep = seq_len(n_reps))

# ---- B and C. Real panels ------------------------------------------------------------------------------
wide <- function(d, unit, time, y) { m <- tapply(d[[y]], list(d[[time]], d[[unit]]), sum); m[order(as.numeric(rownames(m))), , drop = FALSE] }
kansas <- as.data.frame(augsynth::kansas); kansas$t <- kansas$year + (kansas$qtr - 1) / 4
data("smoking", package = "tidysynth"); data("scpi_germany", package = "scpi")
panels <- list(
  kansas     = list(Y = wide(kansas, "state", "t", "gdpcapita"), treated = "Kansas", start = 2012.25),
  california = list(Y = wide(as.data.frame(smoking), "state", "year", "cigsale"), treated = "California", start = 1989),
  germany    = list(Y = wide(as.data.frame(scpi_germany), "country", "year", "gdp"), treated = "West Germany", start = 1991))
real_grid <- imap_dfr(panels, function(p, nm) tibble(panel = nm, unit = setdiff(colnames(p$Y), p$treated)))
run_real_fake <- function(panel, unit) {
  P <- panels[[panel]]; Y <- P$Y[, colnames(P$Y) != P$treated]; post <- as.numeric(rownames(Y)) >= P$start
  evaluate_panel(unname(Y), match(unit, colnames(Y)), post) |> mutate(panel, unit)
}
real_treated <- function(panel) {
  P <- panels[[panel]]; post <- as.numeric(rownames(P$Y)) >= P$start; I <- match(P$treated, colnames(P$Y)); h <- sum(post)
  core <- edge_core(unname(P$Y), I, post, estimators = c("dsc", "ridge"))
  map_dfr(levels_conf, function(cl) imap_dfr(edge_bands(core, cl)[c("dsc_space", "ridge_space")], function(x, nm)
    tibble(panel, method = sub("_space", "", nm), conf = cl, effect = x$gap[h], low = x$gap[h] - x$hw_pw[h], high = x$gap[h] + x$hw_pw[h],
           cum_effect = x$gapc[h], cum_low = x$gapc[h] - x$hwc_pw[h], cum_high = x$gapc[h] + x$hwc_pw[h], c_ridge = core$c_ridge)))
}

# ---- D. Agreement with augsynth -----------------------------------------------------------------------
augsynth_compare <- function(panel) {
  P <- panels[[panel]]; Y <- P$Y; post <- as.numeric(rownames(Y)) >= P$start; I <- match(P$treated, colnames(Y))
  pre_idx <- which(!post); post_idx <- which(post); donors <- setdiff(seq_len(ncol(Y)), I)
  c_r <- choose_c(unname(Y), donors, pre_idx, length(post_idx))
  proto <- fit_target(unname(Y)[, I], unname(Y)[, donors], pre_idx, post_idx, c_r)$gap
  rp <- ridge_parts(unname(Y)[, donors], pre_idx, post_idx); lam <- c_r * mean(rowSums(rp$Mc^2))
  long <- tibble(unit = rep(colnames(Y), each = nrow(Y)), time = rep(seq_len(nrow(Y)), ncol(Y)), y = as.vector(Y)) |>
    mutate(treat = as.integer(unit == P$treated & time >= min(post_idx)))
  att <- function(...) { a <- augsynth::augsynth(y ~ treat, unit, time, long, progfunc = "Ridge", scm = TRUE, fixedeff = TRUE, ...)
    s <- suppressMessages(summary(a, inf = FALSE)); s$att$Estimate[s$att$Time %in% post_idx] }
  tibble(panel, period = seq_along(post_idx), prototype = proto, augsynth_own_lambda = att(), augsynth_proto_lambda = att(lambda = lam),
         c_ridge = c_r)
}

if (!exists("dry_run") || !isTRUE(dry_run)) {
  plan(multisession, workers = n_workers)
  t0 <- Sys.time()
  sim <- future_pmap(sim_grid, possibly(run_sim_rep, otherwise = NULL), .options = furrr_options(seed = TRUE, scheduling = Inf, packages = c("tidyverse", "quadprog"))) |> list_rbind()
  real <- future_pmap(real_grid, possibly(run_real_fake, otherwise = NULL), .options = furrr_options(seed = TRUE, scheduling = Inf, packages = c("tidyverse", "quadprog"))) |> list_rbind()
  plan(sequential)
  treated_results <- map_dfr(names(panels), real_treated)
  aug <- map_dfr(names(panels), possibly(augsynth_compare, otherwise = NULL))
  run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  saveRDS(list(sim = sim, real = real, treated = treated_results, augsynth = aug, run_minutes = run_minutes), file.path(out_dir, "results.rds"))
  sim_summary <- sim |> group_by(treated, factors, n_donors, T0, conf, method, measure) |>
    summarise(n = n(), bias = mean(err_pct), rmse = sqrt(mean(err_pct^2)), final_coverage_pw = mean(final_cov_pw), path_coverage_sim = mean(path_cov_sim),
              finite = mean(is.finite(hw_pw_pct)), median_hw_pw_pct = median(hw_pw_pct[is.finite(hw_pw_pct)]),
              across(starts_with("power_"), mean), .groups = "drop")
  real_summary <- real |> group_by(panel, conf, method, measure) |>
    summarise(n = n(), bias = mean(err_pct), rmse = sqrt(mean(err_pct^2)), final_coverage_pw = mean(final_cov_pw), path_coverage_sim = mean(path_cov_sim),
              median_hw_pw_pct = median(hw_pw_pct[is.finite(hw_pw_pct)]), across(starts_with("power_"), mean), .groups = "drop")
  write_csv(sim_summary, file.path(out_dir, "summary.csv"))
  write_csv(real_summary, file.path(out_dir, "real_panels_summary.csv"))
  message("Finished in ", run_minutes, " min.")
}
