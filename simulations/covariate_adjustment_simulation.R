##############################################################################
# covariate_adjustment_simulation.R
#
# Does the pooled covariate adjustment ("Option D", adjust_covariates in
# dsc_symcal_boot()) remove confounding from a shock that arrives AT THE SAME
# TIME as the treatment (like fentanyl and Oregon's Measure 110)?
#
# Base data: tidysynth::smoking (California dropped, 38 states, 1970-2000);
# each state takes a turn as the fake treated unit, treatment in 1989.
# Confounding shock z_it: a ramp from 0 to 1 over about 3 years, arriving
#   - for the treated state exactly at treatment (1989);
#   - for donors at a random year in 1976-2000 (about 20% never).
# Scenarios for the shock's effect on the outcome:
#   additive      + 25 packs per capita when fully arrived (common effect; the
#                 method's assumption holds)
#   proportional  + 25% of the state's own level (effect differs across states;
#                 the common-gamma assumption is violated)
# Methods (demeaned SC + symmetric calibrated bootstrap, 95%):
#   none        no covariate
#   regressor   controls = "z", control_units = "treated" (coefficient from the
#               unit's own pre-period; not identified when z is 0 before treatment)
#   adjust      adjust_covariates = "z" (Option D)
# Bands are shift-invariant (gamma never uses the treated post period), so
# coverage is computed at 0% injected effect and power exactly.
#
# Output (results_covariate_adjustment/): tasks/*.rds, results.rds, summary.csv,
#   covariate_adjustment_table.html (DT, fixed first column and header).
##############################################################################

library(tidyverse)
library(future)
library(furrr)


# ---- 1. Settings -----------------------------------------------------------------
script_file   <- "simulations/covariate_adjustment_simulation.R"
functions_src <- "synthetic_control_functions/synthetic_control_functions_claude.R"
out_dir       <- "simulations/results_covariate_adjustment"
tasks_dir     <- file.path(out_dir, "tasks")
n_workers     <- 10
shock_size    <- 25                      # additive: packs per capita; proportional: 0.25 x level
effects_pct   <- seq(-0.25, 0.25, by = 0.05)
seed          <- 12112025
dir.create(tasks_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 2. Data with a simulated confounding shock --------------------------------------
base_data <- tidysynth::smoking |>
  filter(state != "California") |>
  select(state, year, cigsale) |>
  as.data.frame()
all_states <- sort(unique(base_data$state))

ramp <- function(year, arrival) if (is.infinite(arrival)) rep(0, length(year)) else plogis((year - arrival - 1) / 0.6)

make_data <- function(scenario, treated) {
  set.seed(seed + match(treated, all_states))                     # same donor arrivals for every method and scenario
  arrivals <- tibble(state = all_states,
                     arrival = ifelse(runif(length(all_states)) < 0.2, Inf, sample(1976:2000, length(all_states), replace = TRUE))) |>
    mutate(arrival = if_else(state == treated, 1989, arrival))
  base_data |>
    left_join(arrivals, by = "state") |>
    group_by(state) |>
    mutate(z = ramp(year, arrival[1]),
           cigsale = if (scenario == "additive") cigsale + shock_size * z else cigsale * (1 + shock_size / 100 * z)) |>
    ungroup() |>
    select(state, year, cigsale, z) |>
    as.data.frame()
}


# ---- 3. One task: one scenario, treated state and method -------------------------------
sim_one <- function(scenario, treated, method) {
  task_file <- file.path(tasks_dir, sprintf("%s_%s_%s.rds", scenario, method, gsub(" ", "_", treated)))
  if (file.exists(task_file)) return(readRDS(task_file))
  if (!exists("dsc_symcal_boot", envir = globalenv())) source(functions_src)

  d <- make_data(scenario, treated)
  args <- list(data = d, unit_col = "state", time_cols = "year", outcome_col = "cigsale",
               treated_unit = treated, treatment_period = c(year = 1989))
  if (method == "regressor") args <- c(args, list(controls = "z", control_units = "treated"))
  if (method == "adjust")    args <- c(args, list(adjust_covariates = "z"))
  r <- do.call(dsc_symcal_boot, args)

  # Untreated outcome = the treated state's data WITHOUT the injected treatment (none injected here);
  # the true effect is 0, so each gap is the estimation error.
  e <- r$effects |> filter(post)
  result <- bind_rows(
    tibble(measure = "per period", t = seq_len(nrow(e)), gap0 = e$effect,
           hw_pw = (e$eff_high_pw - e$eff_low_pw) / 2, hw_sim = (e$eff_high_sim - e$eff_low_sim) / 2, y0 = e$observed),
    tibble(measure = "cumulative", t = seq_len(nrow(e)), gap0 = e$cum_effect,
           hw_pw = (e$cum_eff_high_pw - e$cum_eff_low_pw) / 2, hw_sim = (e$cum_eff_high_sim - e$cum_eff_low_sim) / 2,
           y0 = cumsum(e$observed))) |>
    mutate(scenario = scenario, treated = treated, method = method,
           gamma = if (is.null(r$adjustment)) NA_real_ else unname(r$adjustment$gamma),
           n_failed_donors = nrow(r$failed_donors))
  saveRDS(result, task_file)
  result
}


# ---- 4. Grid of tasks, run in parallel with future_pmap ---------------------------------
grid <- expand_grid(scenario = c("additive", "proportional"), treated = all_states, method = c("none", "regressor", "adjust"))
plan(multisession, workers = n_workers)
t0 <- Sys.time()
results <- future_pmap(grid, possibly(sim_one, otherwise = NULL),
                       .options = furrr_options(seed = TRUE, packages = c("tidyverse", "CVXR", "meboot")))
plan(sequential)
run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
failed <- grid[map_lgl(results, is.null), ]
long <- bind_rows(results)


# ---- 5. Summaries (tidyverse) --------------------------------------------------------------
# Accuracy: error at the end of the post period, % of the observed outcome (true effect 0)
accuracy <- long |>
  group_by(scenario, method, measure, treated) |>
  summarise(error_pct = 100 * last(gap0) / last(y0), .groups = "drop") |>
  group_by(scenario, method, measure) |>
  summarise(mean_bias = mean(error_pct), bias_ci_lo = t.test(error_pct)$conf.int[1], bias_ci_hi = t.test(error_pct)$conf.int[2],
            rmse = sqrt(mean(error_pct^2)), .groups = "drop")

# Coverage: pointwise at the final period; simultaneous over the whole post period
coverage <- long |>
  group_by(scenario, method, measure, treated) |>
  summarise(final_covered_pw = abs(last(gap0)) <= last(hw_pw), path_covered_sim = all(abs(gap0) <= hw_sim), .groups = "drop") |>
  group_by(scenario, method, measure) |>
  summarise(n_states = n(), final_coverage_pw = mean(final_covered_pw), path_coverage_sim = mean(path_covered_sim), .groups = "drop")

# Power: pointwise at the final period, simultaneous anywhere (effects added exactly)
power <- long |>
  crossing(pct = effects_pct) |>
  mutate(shifted = gap0 + pct * y0) |>
  group_by(scenario, method, measure, treated, pct) |>
  summarise(detect_final_pw = abs(last(shifted)) > last(hw_pw), detect_any_sim = any(abs(shifted) > hw_sim), .groups = "drop") |>
  group_by(scenario, method, measure, pct) |>
  summarise(power_final_pw = mean(detect_final_pw), power_any_sim = mean(detect_any_sim), .groups = "drop")

gamma_summary <- long |> filter(method == "adjust") |> distinct(scenario, treated, gamma) |>
  group_by(scenario) |> summarise(median_gamma = median(gamma), q10 = quantile(gamma, 0.1), q90 = quantile(gamma, 0.9))


# ---- 6. Save and DT table ---------------------------------------------------------------------
saveRDS(list(settings = list(shock_size = shock_size, effects_pct = effects_pct, n_workers = n_workers, run_minutes = run_minutes),
             accuracy = accuracy, coverage = coverage, power = power, gamma = gamma_summary, failed = failed, long = long),
        file.path(out_dir, "results.rds"))
write_csv(left_join(accuracy, coverage, by = c("scenario", "method", "measure")), file.path(out_dir, "summary.csv"))

dt_tab <- accuracy |>
  left_join(coverage, by = c("scenario", "method", "measure")) |>
  left_join(power |> filter(round(abs(pct), 2) %in% c(0.10, 0.25)) |>
              mutate(col = sprintf("power %+.0f%% (final pw / any sim)", 100 * pct),
                     val = sprintf("%.0f / %.0f", 100 * power_final_pw, 100 * power_any_sim)) |>
              select(scenario, method, measure, col, val) |> pivot_wider(names_from = col, values_from = val),
            by = c("scenario", "method", "measure")) |>
  transmute(Method = sprintf("%s | %s | %s", scenario, method, measure),
            `bias [95% CI]` = sprintf("%.2f [%.2f, %.2f]", mean_bias, bias_ci_lo, bias_ci_hi), rmse = round(rmse, 2),
            `final coverage, pointwise` = round(100 * final_coverage_pw), `whole-path coverage, simultaneous` = round(100 * path_coverage_sim),
            across(starts_with("power")))
widget <- DT::datatable(dt_tab, rownames = FALSE, extensions = "FixedColumns", filter = "top",
  caption = "Covariate-adjustment simulation (California panel, confounding shock arriving with the treatment). Values in %. Bias and RMSE: error at the end of the post period, % of the outcome. Power: final period (pointwise) / anywhere (simultaneous).",
  options = list(scrollX = TRUE, scrollY = "70vh", scrollCollapse = TRUE, paging = FALSE, fixedColumns = list(leftColumns = 1)))
tmp_html <- file.path(tempdir(), "covariate_adjustment_table.html")
htmlwidgets::saveWidget(widget, tmp_html, selfcontained = TRUE, title = "Covariate adjustment simulation")
file.copy(tmp_html, file.path(out_dir, "covariate_adjustment_table.html"), overwrite = TRUE)

marker <- "# ---- Results (written automatically when this script is run) ----"
script_lines <- read_lines(script_file)
cut_at <- which(script_lines == marker)
if (length(cut_at) > 0) script_lines <- script_lines[seq_len(cut_at[1] - 1)]
write_lines(c(script_lines, marker,
  sprintf("# Run: %s | %d tasks (%d failed) | %s min on %d workers", format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(grid), nrow(failed), run_minutes, n_workers),
  paste0("# ", knitr::kable(dt_tab |> select(1:5), format = "simple"))), script_file)
message("Finished in ", run_minutes, " min.")

# ---- Results (written automatically when this script is run) ----
# Run: 2026-10-04 00:43 | 228 tasks (0 failed) | 11.3 min on 10 workers
# Method                                  bias [95% CI]             rmse   final coverage, pointwise   whole-path coverage, simultaneous
# --------------------------------------  ---------------------  -------  --------------------------  ----------------------------------
# additive | adjust | cumulative          -0.26 [-3.02, 2.50]       8.29                          95                                  95
# additive | adjust | per period          -0.80 [-4.90, 3.29]      12.32                          97                                  95
# additive | none | cumulative            13.22 [10.05, 16.39]     16.29                          92                                  97
# additive | none | per period            12.20 [7.83, 16.57]      17.92                          95                                  92
# additive | regressor | cumulative       56.05 [8.97, 103.12]    152.03                          95                                 100
# additive | regressor | per period       64.98 [6.87, 123.10]    186.18                          97                                 100
# proportional | adjust | cumulative      -0.74 [-4.05, 2.58]       9.99                          95                                  97
# proportional | adjust | per period      -1.83 [-7.13, 3.48]      16.03                          95                                  97
# proportional | none | cumulative        13.56 [10.54, 16.57]     16.30                          95                                  95
# proportional | none | per period        12.15 [7.52, 16.78]      18.46                          92                                  92
# proportional | regressor | cumulative   54.58 [6.35, 102.81]    154.73                          95                                 100
# proportional | regressor | per period   64.93 [3.08, 126.78]    196.70                         100                                 100
