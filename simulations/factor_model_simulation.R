##############################################################################
# factor_model_simulation.R
#
# Fully simulated evaluation of the placebo-based methods (dsc_placebo(),
# dsc_symcal_boot()) and the covariate adjustment, using a linear factor model
# (the model synthetic control is justified by), so that the treated unit can be
# made NON-exchangeable with the donors.
#
# Data generating process (log scale, then exponentiated: outcomes positive,
# effects proportional):
#   Y0_it = exp(5 + alpha_i + 0.01 t + sum_k lambda_ik f_kt + e_it) + gamma z_it
#   alpha_i ~ N(0, 0.3^2); 2 factors; e_it AR(1) (rho 0.5), sd 0.02
#   factors: "stationary"   f_kt = 0.8 f_k,t-1 + u, u ~ N(0, 0.05^2)
#            "random_walk"  f_kt = f_k,t-1 + u, u ~ N(0.005, 0.03^2)
#   donor loadings lambda_ik ~ N(0, 1)
# Treated-unit scenarios:
#   exchangeable  loadings and noise as donors
#   edge          loadings at the donors' maximum (each factor)
#   outside       loadings 1 SD beyond the donors' maximum
#   selected      all units drawn alike; the unit with the largest first loading
#                 becomes the treated unit
#   noisy         exchangeable loadings, twice the noise
# Covariate (covariate block only): z_it switches from 0 to 1 at a unit-specific
#   period; gamma = 0.08 x exp(5) (about 11.9, roughly 8% of a typical level), added
#   on the outcome scale so the effect is common to all units. Treated unit: switches on
#   at treatment.
#     exogenous   random timing for donors
#     confounded  timing related to the first loading (more exposed = earlier)
# Effects: methods are shift-invariant, so each replicate is run with no effect;
#   power is computed exactly for constant effects of +/-5, 10, 20% and for an
#   effect ramping linearly to -10%.
#
# Blocks:
#   core        dsc_placebo() (relative, raw, rmspe scales; classic RMSPE ratio)
#   bootstrap   dsc_symcal_boot() (subset of cells; slower)
#   covariate   dsc_placebo() relative with no adjustment / TWFE / IFE adjustment
#
# Output (results_factor_model/): tasks/*.rds, results.rds, summary.csv,
#   factor_model_table.html (DT, fixed first column and header).
#
# Run from the project root:  source("simulations/factor_model_simulation.R")
##############################################################################

library(tidyverse)
library(future)
library(furrr)


# ---- 1. Settings ----------------------------------------------------------------------
script_file   <- "simulations/factor_model_simulation.R"
functions_src <- "synthetic_control_functions/synthetic_control_functions_claude.R"
out_dir       <- "simulations/results_factor_model"
tasks_dir     <- file.path(out_dir, "tasks")
n_workers     <- 10
reps_per_task <- 10                       # replicates saved together in one task file
conf          <- 0.95
T1            <- 12                       # post-treatment periods
effects_pct   <- c(-0.20, -0.10, -0.05, 0.05, 0.10, 0.20)
seed          <- 20261004
gamma_level   <- 0.08 * exp(5)            # covariate effect, outcome scale

# Simulation grid (edit reps here; see the run plan before launching)
grid_core <- expand_grid(block = "core", treated = c("exchangeable", "edge", "outside", "selected", "noisy"),
                         factors = c("stationary", "random_walk"), n_donors = c(20, 40), T0 = c(10, 20, 40, 80),
                         covariate = "none", reps = 200)
grid_boot <- expand_grid(block = "bootstrap", treated = c("exchangeable", "edge", "noisy"), factors = "random_walk",
                         n_donors = 40, T0 = c(20, 40), covariate = "none", reps = 100)
grid_cov  <- expand_grid(block = "covariate", treated = c("exchangeable", "edge"), factors = "random_walk",
                         n_donors = 40, T0 = c(20, 40, 80), covariate = c("exogenous", "confounded"), reps = 100)
dir.create(tasks_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 2. Data generator --------------------------------------------------------------------
simulate_panel <- function(treated, factors, n_donors, T0, covariate, rep_seed) {
  set.seed(rep_seed)
  TT <- T0 + T1; N <- n_donors + 1; K <- 2
  f <- sapply(seq_len(K), function(k) {
    u <- if (factors == "stationary") rnorm(TT, 0, 0.05) else rnorm(TT, 0.005, 0.03)
    if (factors == "stationary") as.numeric(stats::filter(u, 0.8, method = "recursive")) else cumsum(u)
  })                                                                          # TT x K
  lam <- matrix(rnorm(N * K), N, K)                                           # row 1 = treated
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


# ---- 3. One replicate: every method for the block -----------------------------------------
bands_tbl <- function(fit, method) {
  e <- fit$effects |> filter(post)
  bind_rows(
    tibble(measure = "per period", t = seq_len(nrow(e)), gap0 = e$effect, y0 = e$observed,          # no effect injected: observed = untreated outcome
           hw_pw = (e$eff_high_pw - e$eff_low_pw) / 2, hw_sim = (e$eff_high_sim - e$eff_low_sim) / 2),
    tibble(measure = "cumulative", t = seq_len(nrow(e)), gap0 = e$cum_effect, y0 = cumsum(e$observed),
           hw_pw = (e$cum_eff_high_pw - e$cum_eff_low_pw) / 2, hw_sim = (e$cum_eff_high_sim - e$cum_eff_low_sim) / 2)) |>
    mutate(method = method,
           p_ratio = if ("RMSPE ratio (Abadie et al.)" %in% fit$test$test) fit$test$p_value[fit$test$test == "RMSPE ratio (Abadie et al.)"] else NA_real_,
           gamma = if (!is.null(fit$adjustment)) unname(fit$adjustment$gamma[1]) else NA_real_,
           n_factors = if (!is.null(fit$adjustment)) fit$adjustment$n_factors else NA_integer_,
           failed_donors = nrow(fit$failed_donors))
}

run_rep <- function(cell, rep) {
  d <- simulate_panel(cell$treated, cell$factors, cell$n_donors, cell$T0, cell$covariate,
                      rep_seed = seed + 100000 * cell$cell_id + rep)
  tp <- c(period = cell$T0 + 1)
  call <- function(fun, ...) fun(as.data.frame(d), "unit", "period", "y", "u01", treatment_period = tp, conf_level = conf, ...)
  safe <- function(method, fun, ...) tryCatch(bands_tbl(call(fun, ...), method),
                                              error = function(e) tibble(method = method, error = conditionMessage(e)))
  out <- switch(cell$block,
    core = bind_rows(safe("placebo_relative", dsc_placebo, scale = "relative"),
                     safe("placebo_raw", dsc_placebo, scale = "raw"),
                     safe("placebo_rmspe", dsc_placebo, scale = "rmspe")),
    bootstrap = bind_rows(safe("bootstrap", dsc_symcal_boot),
                          safe("placebo_relative", dsc_placebo, scale = "relative")),
    covariate = bind_rows(safe("no_adjustment", dsc_placebo, scale = "relative"),          # the default scale when this was run
                          safe("adjust_twfe", dsc_placebo, scale = "relative", adjust_covariates = "z", adjust_model = "twfe"),
                          safe("adjust_ife", dsc_placebo, scale = "relative", adjust_covariates = "z", adjust_model = "ife")))
  out |> mutate(rep = rep)
}

# One task = one cell x a chunk of replicates (saved; skipped if done)
run_task <- function(cell_id, chunk) {
  task_file <- file.path(tasks_dir, sprintf("cell%03d_chunk%02d.rds", cell_id, chunk))
  if (file.exists(task_file)) return(readRDS(task_file))
  if (!exists("dsc_placebo", envir = globalenv())) source(functions_src)
  cell <- grid |> filter(cell_id == !!cell_id) |> as.list()
  reps <- ((chunk - 1) * reps_per_task + 1):min(chunk * reps_per_task, cell$reps)
  res <- map_dfr(reps, function(r) tryCatch(run_rep(cell, r),
                 error = function(e) tibble(method = "ERROR", error = conditionMessage(e), rep = r))) |>
    mutate(cell_id = cell_id)
  saveRDS(res, task_file)
  res
}


# ---- 4. Grid of tasks, run in parallel with future_pmap --------------------------------------
grid <- bind_rows(grid_core, grid_boot, grid_cov) |> mutate(cell_id = row_number())
tasks <- grid |> select(cell_id, reps) |> mutate(chunk = map(reps, ~ seq_len(ceiling(.x / reps_per_task)))) |>
  unnest(chunk) |> select(cell_id, chunk)

if (!exists("dry_run") || !isTRUE(dry_run)) {
  plan(multisession, workers = n_workers)
  t0 <- Sys.time()
  results <- future_pmap(tasks, possibly(run_task, otherwise = NULL),
                         .options = furrr_options(seed = TRUE, packages = c("tidyverse", "CVXR", "meboot", "fect")))
  plan(sequential)
  run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  long <- bind_rows(results)
  if (!"error" %in% names(long)) long$error <- NA_character_
  errors <- long |> filter(!is.na(error)) |> distinct(cell_id, method, rep, error)          # failed method x replicate
  long <- long |> filter(is.na(error)) |> left_join(grid, by = "cell_id")
  failures <- grid |> select(cell_id, reps) |>
    inner_join(long |> distinct(cell_id, method), by = "cell_id") |>
    left_join(errors |> count(cell_id, method, name = "failed"), by = c("cell_id", "method")) |>
    mutate(failed = coalesce(failed, 0L), failure_rate = failed / reps)


  # ---- 5. Summaries (tidyverse) -------------------------------------------------------------------
  # Accuracy (true effect 0): error at the end of the post period, % of the untreated outcome
  accuracy <- long |>
    group_by(cell_id, method, measure, rep) |> summarise(err = 100 * last(gap0) / last(y0), .groups = "drop") |>
    group_by(cell_id, method, measure) |>
    summarise(bias = mean(err), rmse = sqrt(mean(err^2)), .groups = "drop")
  # Coverage, false positives and width
  coverage <- long |>
    group_by(cell_id, method, measure, rep) |>
    summarise(path_cov_sim = all(abs(gap0) <= hw_sim), final_cov_pw = abs(last(gap0)) <= last(hw_pw),
              hw_pw_pct = 100 * last(hw_pw) / last(y0), p_ratio = first(p_ratio), .groups = "drop") |>
    group_by(cell_id, method, measure) |>
    summarise(n = n(), path_coverage_sim = mean(path_cov_sim), final_coverage_pw = mean(final_cov_pw),
              median_hw_pw_pct = median(hw_pw_pct[is.finite(hw_pw_pct)]),
              ratio_false_pos = mean(p_ratio <= 1 - conf + 1e-12, na.rm = TRUE), .groups = "drop")
  # Power: constant effects and a ramp to -10% (added exactly)
  shifts <- bind_rows(tibble(effect = sprintf("%+.0f%%", 100 * effects_pct), pct = effects_pct, shape = "constant"),
                      tibble(effect = "ramp to -10%", pct = -0.10, shape = "ramp"))
  power <- long |>
    crossing(shifts) |>
    arrange(cell_id, method, measure, rep, effect, t) |>
    group_by(cell_id, method, measure, rep, effect) |>
    mutate(w = if_else(shape == "constant", 1, t / max(t)),
           shift_t = pct * w * if_else(measure == "per period", y0, y0 - lag(y0, default = 0)),
           shifted = gap0 + if_else(measure == "per period", shift_t, cumsum(shift_t))) |>
    summarise(final_pw = abs(last(shifted)) > last(hw_pw), any_sim = any(abs(shifted) > hw_sim), .groups = "drop") |>
    group_by(cell_id, method, measure, effect) |>
    summarise(power_final_pw = mean(final_pw), power_any_sim = mean(any_sim), .groups = "drop")
  gamma_summary <- long |> filter(block == "covariate", method != "no_adjustment") |>
    distinct(cell_id, method, rep, gamma, n_factors) |>
    group_by(cell_id, method) |>
    summarise(gamma_bias = mean(gamma) - gamma_level, gamma_rmse = sqrt(mean((gamma - gamma_level)^2)),
              median_factors = median(n_factors), .groups = "drop")


  # ---- 6. Save and DT table --------------------------------------------------------------------------
  saveRDS(list(settings = list(grid = grid, conf = conf, T1 = T1, effects_pct = effects_pct, seed = seed,
                               n_workers = n_workers, run_minutes = run_minutes),
               accuracy = accuracy, coverage = coverage, power = power, gamma = gamma_summary, errors = errors,
               failures = failures),
          file.path(out_dir, "results.rds"))
  summary_tbl <- grid |> select(-reps) |>
    inner_join(accuracy, by = "cell_id") |> left_join(coverage, by = c("cell_id", "method", "measure"))
  write_csv(summary_tbl, file.path(out_dir, "summary.csv"))

  dt_tab <- summary_tbl |>
    left_join(power |> filter(effect %in% c("-10%", "+10%", "ramp to -10%")) |>
                mutate(val = sprintf("%.0f / %.0f", 100 * power_final_pw, 100 * power_any_sim)) |>
                select(cell_id, method, measure, effect, val) |>
                pivot_wider(names_from = effect, values_from = val, names_prefix = "power "),
              by = c("cell_id", "method", "measure")) |>
    left_join(gamma_summary |> select(cell_id, method, gamma_bias), by = c("cell_id", "method")) |>
    left_join(failures |> select(cell_id, method, failure_rate), by = c("cell_id", "method")) |>
    transmute(Setting = sprintf("%s | %s | %s | N=%d | T0=%d | cov: %s | %s | %s", block, treated, factors, n_donors, T0,
                                covariate, method, measure),
              `bias (%)` = round(bias, 2), `RMSE (%)` = round(rmse, 2),
              `final coverage (pw)` = round(100 * final_coverage_pw), `path coverage (sim)` = round(100 * path_coverage_sim),
              `RMSPE-ratio false pos.` = round(100 * ratio_false_pos), `final half-width (%)` = round(median_hw_pw_pct, 1),
              across(starts_with("power ")), `gamma bias` = round(gamma_bias, 3),
              `failed replicates (%)` = round(100 * failure_rate, 1))
  widget <- DT::datatable(dt_tab, rownames = FALSE, extensions = "FixedColumns", filter = "top",
    caption = "Factor-model simulation. Values in % except gamma bias. Power: pointwise final period / simultaneous anywhere. Coverage and false positives at no effect.",
    options = list(scrollX = TRUE, scrollY = "70vh", scrollCollapse = TRUE, paging = FALSE, fixedColumns = list(leftColumns = 1)))
  tmp_html <- file.path(tempdir(), "factor_model_table.html")
  htmlwidgets::saveWidget(widget, tmp_html, selfcontained = TRUE, title = "Factor-model simulation")
  file.copy(tmp_html, file.path(out_dir, "factor_model_table.html"), overwrite = TRUE)
  message("Finished in ", run_minutes, " min; ", nrow(errors), " replicate errors.")
}

