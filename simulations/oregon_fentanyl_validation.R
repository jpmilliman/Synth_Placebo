# oregon_fentanyl_validation.R
#
# Does dsc_placebo()'s covariate adjustment remove a confounding shock that
# arrives together with the treatment? Semi-synthetic validation on the panel of
# Zoorob et al. (2024): half-yearly overdose death rates per 100k and the share of
# seized drugs containing fentanyl (NFLIS), 2008 H1 - 2022 H2, 49 states + DC
# minus Washington (also treated) and Oregon (the real treated unit).
#
# Design
#   * Each of the 49 remaining units (48 states + DC) in turn is a fake treated
#     unit, treated from 2021 H1 (as Oregon), with the other 48 as donors. Its
#     real effect is taken to be 0.
#   * Confounding shock (as in Oregon): the unit's fentanyl share is raised after
#     treatment by Oregon's own increase since 2020 H2 (capped at 100%), and its
#     overdose death rate by gamma_true x the added fentanyl share. Scenarios:
#     no shock; gamma_true = 0.155 (deaths per 100k per percentage point, about
#     the gamma estimated from the panel, so the adjustment model is right); and
#     gamma_true = 0.12 (the IFE estimate with Oregon as the treated unit, so
#     gamma is over-estimated by about 0.035).
#   * True effect: 0, or +/-25% of the unit's mean death rate in the last two
#     pre-treatment half-years, added in every post period.
#   * Methods: dsc_placebo() unadjusted, adjusted with interactive fixed effects
#     (default), adjusted with two-way fixed effects. RMSPE scale, 95%.
#
# Shift invariance: the weights, gamma-hat (estimated on untreated cells) and all
# intervals depend only on pre-treatment data and the donors, so the shock and
# the effect shift the treated unit's gap exactly:
#     unadjusted: + gamma_true x added + tau
#     adjusted:   + (gamma_true - gamma_hat) x added + tau
# Each unit and method is therefore fitted once on the real data, and every
# scenario is evaluated from that fit.
#
# Run from the project root:  source("simulations/oregon_fentanyl_validation.R")
# Uses future + furrr::future_pmap (10 workers); resumes from saved tasks.

suppressPackageStartupMessages({
  library(tidyverse)
  library(future)
  library(furrr)
})

functions_src <- "synthetic_control_functions/synthetic_control_functions_claude.R"
source(functions_src)
out_dir   <- file.path("simulations", "results_oregon_fentanyl")
tasks_dir <- file.path(out_dir, "tasks")
dir.create(tasks_dir, recursive = TRUE, showWarnings = FALSE)

# ---- data -------------------------------------------------------------------
panel <- read_csv("data/oregon_fentanyl/data_final.csv", show_col_types = FALSE) |>
  filter(Year <= 2022, ST != "WA") |>
  select(State, Year, half, trend, od_death_rate_half, nflis_fentanyl_percent_total) |>
  as.data.frame()

post_trend  <- 27:30                     # 2021 H1 - 2022 H2
oregon_fent <- panel |> filter(State == "Oregon") |> arrange(trend) |> pull(nflis_fentanyl_percent_total)
oregon_inc  <- oregon_fent[post_trend] - oregon_fent[26]   # Oregon's increase since 2020 H2

sim_panel <- panel |> filter(State != "Oregon")
units     <- sort(unique(sim_panel$State))

grid <- expand_grid(unit = units, method = c("unadjusted", "ife", "twfe"))

# ---- one task: fit a fake treated unit with one method ------------------------
run_task <- function(unit, method) {
  task_file <- file.path(tasks_dir, sprintf("%s_%s.rds", gsub(" ", "_", unit), method))
  if (file.exists(task_file)) return(readRDS(task_file))
  if (!exists("dsc_placebo", envir = globalenv())) source(functions_src)
  adj <- if (method == "unadjusted") NULL else "nflis_fentanyl_percent_total"
  res <- tryCatch({
    fit <- dsc_placebo(sim_panel, unit_col = "State", time_cols = c("Year", "half"),
                       outcome_col = "od_death_rate_half", treated_unit = unit,
                       treatment_period = c(Year = 2021, half = 1),
                       adjust_covariates = adj, adjust_model = if (method == "twfe") "twfe" else "ife",
                       diagnostics = FALSE)
    e <- fit$effects[fit$effects$post, ]
    tibble(unit = unit, method = method, k = seq_len(nrow(e)),
           gap = e$effect, hw_pw = (e$eff_high_pw - e$eff_low_pw) / 2,
           hw_sim = (e$eff_high_sim - e$eff_low_sim) / 2,
           cum_gap = e$cum_effect, cum_hw_pw = (e$cum_eff_high_pw - e$cum_eff_low_pw) / 2,
           cum_hw_sim = (e$cum_eff_high_sim - e$cum_eff_low_sim) / 2,
           gamma_hat = if (is.null(adj)) 0 else unname(fit$adjustment$gamma[1]),
           error = NA_character_)
  }, error = function(err) tibble(unit = unit, method = method, error = conditionMessage(err)))
  saveRDS(res, task_file)
  res
}

plan(multisession, workers = 10)
fits <- future_pmap(grid, run_task, .options = furrr_options(seed = TRUE, packages = "tidyverse"),
                    .progress = TRUE) |> bind_rows()
plan(sequential)

failed <- fits |> filter(!is.na(error))
fits   <- fits |> filter(is.na(error))

# ---- scenarios -----------------------------------------------------------------
unit_info <- sim_panel |>
  group_by(State) |> arrange(trend, .by_group = TRUE) |>
  summarise(base_level = mean(od_death_rate_half[25:26]),
            fent_last_pre = nflis_fentanyl_percent_total[26],
            added = list(pmin(100, nflis_fentanyl_percent_total[post_trend] + oregon_inc) -
                           nflis_fentanyl_percent_total[post_trend]),
            .groups = "drop") |>
  rename(unit = State) |>
  mutate(low_fentanyl = fent_last_pre < 15)   # fentanyl still rare in 2020 H2, as in Oregon (6.7%)

scenarios <- expand_grid(tibble(shock = c("none", "gamma 0.155", "gamma 0.12"), gamma_true = c(0, 0.155, 0.12)),
                         effect = c(0, -0.25, 0.25))

evaluated <- scenarios |>
  pmap(function(shock, gamma_true, effect) {
    fits |>
      left_join(unit_info |> select(unit, base_level, added, low_fentanyl), by = "unit") |>
      mutate(added_k = if (shock == "none") 0 else map2_dbl(added, k, ~ .x[.y]),
             tau     = effect * base_level,
             shift   = (gamma_true - gamma_hat) * added_k,
             est     = gap + shift + tau) |>
      group_by(unit, method) |> arrange(k, .by_group = TRUE) |>
      mutate(cum_tau = cumsum(tau), cum_est = cum_gap + cumsum(shift) + cum_tau) |>
      summarise(low_fentanyl = first(low_fentanyl),
                err_last      = last(est) - last(tau),
                err_cum       = last(cum_est) - last(cum_tau),
                base_level    = first(base_level),
                cover_pw      = abs(err_last) <= last(hw_pw),
                cover_cum_pw  = abs(err_cum) <= last(cum_hw_pw),
                cover_sim     = all(abs(est - tau) <= hw_sim),
                cover_cum_sim = all(abs(cum_est - cum_tau) <= cum_hw_sim),
                reject_pw     = abs(last(est)) > last(hw_pw),
                reject_cum_pw = abs(last(cum_est)) > last(cum_hw_pw),
                reject_sim    = any(abs(est) > hw_sim),
                reject_cum_sim = any(abs(cum_est) > cum_hw_sim),
                .groups = "drop") |>
      mutate(shock = shock, effect = effect)
  }) |> bind_rows()

summarise_cells <- function(x) x |>
  summarise(n = n(),
            bias_cum_pct = 100 * mean(err_cum / (4 * base_level)),
            rmse_cum_pct = 100 * sqrt(mean((err_cum / (4 * base_level))^2)),
            coverage_last = 100 * mean(cover_pw), coverage_cum = 100 * mean(cover_cum_pw),
            coverage_sim  = 100 * mean(cover_sim), coverage_cum_sim = 100 * mean(cover_cum_sim),
            reject_last = 100 * mean(reject_pw), reject_cum = 100 * mean(reject_cum_pw),
            reject_sim  = 100 * mean(reject_sim), reject_cum_sim = 100 * mean(reject_cum_sim),
            .groups = "drop")

summary_tbl <- bind_rows(
  evaluated |> group_by(shock, effect, method) |> summarise_cells() |> mutate(units = "all"),
  evaluated |> filter(low_fentanyl) |> group_by(shock, effect, method) |> summarise_cells() |> mutate(units = "fentanyl < 15% in 2020 H2")
) |>
  mutate(method = factor(method, c("unadjusted", "ife", "twfe")),
         shock = factor(shock, c("none", "gamma 0.155", "gamma 0.12"))) |>
  arrange(units, shock, effect, method) |>
  relocate(units, shock, effect, method)

gamma_tbl <- fits |> filter(k == 1, method != "unadjusted") |>
  group_by(method) |> summarise(gamma_median = median(gamma_hat), gamma_min = min(gamma_hat), gamma_max = max(gamma_hat))

saveRDS(list(fits = fits, evaluated = evaluated, summary = summary_tbl, gamma = gamma_tbl,
             failed = failed, unit_info = unit_info, oregon_increase = oregon_inc),
        file.path(out_dir, "results.rds"))
write_csv(summary_tbl, file.path(out_dir, "summary.csv"))

if (requireNamespace("DT", quietly = TRUE)) {
  DT::saveWidget(DT::datatable(mutate(summary_tbl, across(where(is.numeric), ~ round(.x, 1))),
                               caption = "Oregon fentanyl validation: bias and RMSE (% of base level, cumulative), coverage and rejection rates (%), 95%, RMSPE scale",
                               options = list(pageLength = 36)),
                 normalizePath(file.path(out_dir, "oregon_fentanyl_table.html"), mustWork = FALSE), selfcontained = TRUE)
}

print(summary_tbl |> filter(units == "all") |> select(shock, effect, method, bias_cum_pct, coverage_cum, coverage_sim, reject_cum, reject_sim))
