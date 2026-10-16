##############################################################################
# pre_period_length_simulation.R
#
# How does the LENGTH OF THE PRE-TREATMENT PERIOD affect bias, coverage, interval
# width and power? Two panels, each with its pre-period cut to several lengths
# (earliest periods dropped; post period unchanged):
#   Kansas      augsynth::kansas, quarterly GDP per capita, Kansas dropped, 49 fake
#               treated states, fake treatment 2010 Q1, 16 post quarters.
#               Pre-period: 80 (full), 40, 20, 12, 8 quarters.
#   California  tidysynth::smoking, annual cigarette sales, California dropped,
#               38 fake treated states, fake treatment 1989, 12 post years.
#               Pre-period: 19 (full), 12, 8, 5 years.
#
# Methods (each state is a fake treated unit; its donors are the other states):
#   demeaned SC + symmetric calibrated fixed-weight bootstrap (dsc_symcal_boot()):
#     boot_pointwise          pointwise band
#     boot_simultaneous       smoothed-median simultaneous band
#   demeaned SC + absolute placebo test     placebo_abs   (pointwise, relative gaps)
#   current LASSO + absolute placebo test   placebo_abs   (weights in [0, 1],
#                                           intercept, smallest lambda on path)
# Levels 95% and 90%; measures: per period and cumulative.
#
# All bands are shift-invariant (an injected effect moves the estimate and band
# equally), so coverage is computed at 0% effect and power exactly for injected
# effects of 0, +/-5, ..., +/-25%.
#
# Parallel: future (multisession, 10 workers) + furrr::future_pmap() over a grid
# of panel x pre-period length x treated state. Each task saves its own result,
# so a rerun skips finished tasks.
#
# Output (results_pre_period_length/): tasks/*.rds, results.rds, summary.csv,
#   bias.csv, pre_period_length_table.html (DT, fixed first column and header).
##############################################################################

library(tidyverse)
library(future)
library(furrr)


# ---- 1. Settings ----------------------------------------------------------------
script_file   <- "simulations/pre_period_length_simulation.R"
fn_script     <- "simulations/demeaned_sc_symmetric_calboot_test.R"   # dsc_symcal_boot()
functions_src <- "synthetic_control_functions/synthetic_control_functions_claude.R"   # meboot_multi()
out_dir       <- "simulations/results_pre_period_length"
tasks_dir     <- file.path(out_dir, "tasks")
n_workers     <- 10
conf_levels   <- c(0.95, 0.90)
effects_pct   <- seq(-0.25, 0.25, by = 0.05)
pre_lengths   <- list(kansas = c(80, 40, 20, 12, 8), california = c(19, 12, 8, 5))
dir.create(tasks_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 2. Data: a panel with its pre-period cut to n_pre periods -----------------------
# Returns columns unit, time, y, plus the treatment time.
get_panel <- function(panel, n_pre) {
  if (panel == "kansas") {
    d <- augsynth::kansas |>
      filter(state != "Kansas", year <= 2013) |>
      transmute(unit = state, time = year + (qtr - 1) / 4, y = gdpcapita)
    treat_time <- 2010
    step <- 0.25                                   # quarterly
  } else {
    d <- tidysynth::smoking |>
      filter(state != "California") |>
      transmute(unit = state, time = year, y = cigsale)
    treat_time <- 1989
    step <- 1                                      # annual
  }
  first_time <- treat_time - n_pre * step          # keep exactly n_pre pre-treatment periods
  list(data = d |> filter(time >= first_time - 1e-9) |> as.data.frame(), treat_time = treat_time)
}


# ---- 3. Methods ------------------------------------------------------------------------
# k-th largest value used by the exact tests (Inf when the level is unreachable)
kth_largest <- function(x, conf) {
  m <- floor(round(1 - conf, 10) * (length(x) + 1))
  if (m == 0) Inf else sort(x, decreasing = TRUE)[m]
}

# Absolute placebo band on relative gaps (gap / synthetic), per period and cumulative.
# gap_t, synth_t: treated post-period vectors; donor_gaps, donor_synth: post periods x donors.
placebo_abs_band <- function(gap_t, synth_t, donor_gaps, donor_synth, conf) {
  rel_donor     <- abs(donor_gaps / donor_synth)
  rel_donor_cum <- abs(apply(donor_gaps, 2, cumsum) / apply(donor_synth, 2, cumsum))
  bind_rows(
    tibble(measure = "per period", t = seq_along(gap_t), gap0 = gap_t,
           hw = apply(rel_donor, 1, kth_largest, conf = conf) * synth_t),
    tibble(measure = "cumulative", t = seq_along(gap_t), gap0 = cumsum(gap_t),
           hw = apply(rel_donor_cum, 1, kth_largest, conf = conf) * cumsum(synth_t)))
}

# Current LASSO synthetic control: weights in [0, 1], intercept, smallest lambda on the path.
# Returns the synthetic series for unit `u` from donors `donor_units` (wide matrix Y: time x unit).
lasso_synth <- function(Y, u, donor_units, pre) {
  X <- Y[, donor_units, drop = FALSE]
  fit <- glmnet::glmnet(X[pre, ], Y[pre, u], alpha = 1, standardize = FALSE, intercept = TRUE,
                        lower.limits = 0, upper.limits = 1, nlambda = 1000)
  cf <- as.matrix(coef(fit, s = min(fit$lambda)))
  as.numeric(cf[1, 1] + X %*% cf[-1, 1])
}


# ---- 4. One task: every method for one panel, pre-period length and treated state ------
sim_one <- function(panel, n_pre, unit) {
  task_file <- file.path(tasks_dir, sprintf("%s_%03d_%s.rds", panel, n_pre, gsub(" ", "_", unit)))
  if (file.exists(task_file)) return(readRDS(task_file))

  # Make the methods available on this worker (once per worker)
  if (!exists("dsc_symcal_boot", envir = globalenv())) {
    source(functions_src)
    ex <- parse(fn_script)
    is_fn <- map_lgl(ex, ~ is.call(.x) && identical(.x[[1]], as.name("<-")) && identical(.x[[2]], as.name("dsc_symcal_boot")))
    eval(ex[[which(is_fn)]], globalenv())
  }

  p <- get_panel(panel, n_pre)
  d <- p$data
  wide <- d |> pivot_wider(names_from = unit, values_from = y) |> arrange(time)
  Y <- as.matrix(select(wide, -time))
  pre <- wide$time < p$treat_time - 1e-9
  post <- !pre
  y0 <- Y[post, unit]                                           # untreated post-period outcome
  donors <- setdiff(colnames(Y), unit)

  out <- list()

  # (a) Demeaned SC with the symmetric calibrated bootstrap, at each level
  for (cf in conf_levels) {
    r <- dsc_symcal_boot(d, "unit", "time", "y", unit, treatment_time = p$treat_time, conf = cf)
    e <- as_tibble(r$effects)
    out[[length(out) + 1]] <- bind_rows(
      tibble(method = "boot_pointwise", measure = "per period", gap0 = e$effect, hw = (e$eff_hi_pw - e$eff_lo_pw) / 2),
      tibble(method = "boot_pointwise", measure = "cumulative", gap0 = e$cum_effect, hw = (e$cum_eff_hi_pw - e$cum_eff_lo_pw) / 2),
      tibble(method = "boot_simultaneous", measure = "per period", gap0 = e$effect, hw = (e$eff_hi_sim_ms - e$eff_lo_sim_ms) / 2),
      tibble(method = "boot_simultaneous", measure = "cumulative", gap0 = e$cum_effect, hw = (e$cum_eff_hi_sim_ms - e$cum_eff_lo_sim_ms) / 2)) |>
      group_by(method, measure) |> mutate(t = row_number()) |> ungroup() |>
      mutate(estimator = "demeaned SC", level = cf)

    # (b) Demeaned SC with the absolute placebo test (same fits, returned by dsc_symcal_boot)
    out[[length(out) + 1]] <- placebo_abs_band(e$effect, e$synthetic,
                                               as.matrix(r$donor_gaps)[post, , drop = FALSE],
                                               as.matrix(r$donor_synth)[post, , drop = FALSE], cf) |>
      mutate(method = "placebo_abs", estimator = "demeaned SC", level = cf)
  }

  # (c) Current LASSO with the absolute placebo test
  synth_t <- lasso_synth(Y, unit, donors, pre)
  donor_synth <- sapply(donors, function(dn) lasso_synth(Y, dn, setdiff(donors, dn), pre))
  donor_gaps <- Y[, donors] - donor_synth
  for (cf in conf_levels) {
    out[[length(out) + 1]] <- placebo_abs_band((Y[, unit] - synth_t)[post], synth_t[post],
                                               donor_gaps[post, , drop = FALSE], donor_synth[post, , drop = FALSE], cf) |>
      mutate(method = "placebo_abs", estimator = "current LASSO", level = cf)
  }

  result <- bind_rows(out) |>
    mutate(panel = panel, n_pre = n_pre, unit = unit,
           y0 = if_else(measure == "per period", y0[t], cumsum(y0)[t]))   # untreated outcome (running total for cumulative)
  saveRDS(result, task_file)
  result
}


# ---- 5. Grid of tasks, run in parallel with future_pmap ---------------------------------
units_in <- function(panel) sort(unique(get_panel(panel, max(pre_lengths[[panel]]))$data$unit))
grid <- map_dfr(names(pre_lengths), function(pn)
  expand_grid(panel = pn, n_pre = pre_lengths[[pn]], unit = units_in(pn)))
message(nrow(grid), " tasks (", sum(file.exists(file.path(tasks_dir, sprintf("%s_%03d_%s.rds", grid$panel, grid$n_pre, gsub(" ", "_", grid$unit))))),
        " already done).")

plan(multisession, workers = n_workers)
t0 <- Sys.time()
results <- future_pmap(grid, possibly(sim_one, otherwise = NULL),
                       .options = furrr_options(seed = TRUE, packages = c("tidyverse", "CVXR", "meboot", "glmnet")),
                       .progress = FALSE)
plan(sequential)
run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
failed <- grid[map_lgl(results, is.null), ]
long <- bind_rows(results)


# ---- 6. Summaries (tidyverse) -------------------------------------------------------------
# Coverage and width for each panel x pre-period length x estimator x method x level x measure
per_state <- long |>
  arrange(t) |>
  group_by(panel, n_pre, unit, estimator, method, level, measure) |>
  summarise(infinite      = any(!is.finite(hw)),
            path_covered  = all(abs(gap0) <= hw),          # true effect inside the band in every post period
            final_covered = abs(last(gap0)) <= last(hw),
            mean_hw_pct   = mean(100 * hw / y0),
            .groups = "drop")

coverage <- per_state |>
  group_by(panel, n_pre, estimator, method, level, measure) |>
  summarise(n_states        = n(),
            path_coverage   = mean(path_covered),
            final_coverage  = mean(final_covered),
            share_infinite  = mean(infinite),
            median_hw_pct   = median(mean_hw_pct[is.finite(mean_hw_pct)]),
            .groups = "drop")

# Power: add each injected effect to the 0% gaps (exact, since bands are shift-invariant)
power <- long |>
  crossing(pct = effects_pct) |>
  mutate(shifted = gap0 + pct * y0) |>
  arrange(t) |>
  group_by(panel, n_pre, unit, estimator, method, level, measure, pct) |>
  summarise(detect_any = any(abs(shifted) > hw), detect_final = abs(last(shifted)) > last(hw), .groups = "drop") |>
  group_by(panel, n_pre, estimator, method, level, measure, pct) |>
  summarise(power_any = mean(detect_any), power_final = mean(detect_final), .groups = "drop")

# Bias and RMSE of each estimator (final period and cumulative), from the 0% gaps
bias <- long |>
  filter(method == "placebo_abs", level == 0.95) |>
  arrange(t) |>
  group_by(panel, n_pre, unit, estimator, measure) |>
  summarise(error_pct = 100 * last(gap0) / last(y0), .groups = "drop") |>
  group_by(panel, n_pre, estimator, measure) |>
  summarise(mean_bias = mean(error_pct),
            bias_ci_lo = if (n() > 1) t.test(error_pct)$conf.int[1] else NA_real_,
            bias_ci_hi = if (n() > 1) t.test(error_pct)$conf.int[2] else NA_real_,
            rmse = sqrt(mean(error_pct^2)), .groups = "drop")


# ---- 7. Save and build the DT table ------------------------------------------------------
saveRDS(list(settings = list(pre_lengths = pre_lengths, conf_levels = conf_levels, effects_pct = effects_pct,
                             n_workers = n_workers, run_minutes = run_minutes),
             coverage = coverage, power = power, bias = bias, per_state = per_state, failed = failed),
        file.path(out_dir, "results.rds"))
write_csv(coverage, file.path(out_dir, "summary.csv"))
write_csv(bias, file.path(out_dir, "bias.csv"))

power_wide <- power |>
  filter(pct %in% c(-0.25, -0.10, 0.10, 0.25)) |>
  mutate(col = sprintf("power %+.0f%% (any / final)", 100 * pct),
         val = sprintf("%.0f / %.0f", 100 * power_any, 100 * power_final)) |>
  select(-pct, -power_any, -power_final) |>
  pivot_wider(names_from = col, values_from = val)

dt_tab <- coverage |>
  left_join(power_wide, by = c("panel", "n_pre", "estimator", "method", "level", "measure")) |>
  left_join(bias |> transmute(panel, n_pre, estimator, measure,
                              `bias [95% CI]` = sprintf("%.2f [%.2f, %.2f]", mean_bias, bias_ci_lo, bias_ci_hi),
                              rmse = round(rmse, 2)),
            by = c("panel", "n_pre", "estimator", "measure")) |>
  mutate(Method = sprintf("%s | %d pre | %s | %s | %s | %.0f%%", panel, n_pre, estimator, method, measure, 100 * level),
         `whole-path coverage` = round(100 * path_coverage), `final coverage` = round(100 * final_coverage),
         `infinite` = round(100 * share_infinite), `median half-width (%)` = round(median_hw_pct, 1)) |>
  arrange(panel, desc(n_pre), estimator, method, measure, desc(level)) |>
  select(Method, `whole-path coverage`, `final coverage`, infinite, `median half-width (%)`, starts_with("power"),
         `bias [95% CI]`, rmse)

widget <- DT::datatable(dt_tab, rownames = FALSE, extensions = "FixedColumns", filter = "top",
  caption = "Pre-period length simulation. Values in %. Row label: panel | pre-treatment periods | estimator | method | measure | level. Power: share of states detected anywhere in the post period / in the final period. Bias and RMSE: estimate error in % of the untreated outcome (final period or cumulative total).",
  options = list(scrollX = TRUE, scrollY = "70vh", scrollCollapse = TRUE, paging = FALSE,
                 fixedColumns = list(leftColumns = 1), autoWidth = TRUE))
tmp_html <- file.path(tempdir(), "pre_period_length_table.html")
htmlwidgets::saveWidget(widget, tmp_html, selfcontained = TRUE, title = "Pre-period length simulation")
file.copy(tmp_html, file.path(out_dir, "pre_period_length_table.html"), overwrite = TRUE)


# Compact summary at the end of this script (95%, per period)
compact <- coverage |>
  filter(level == 0.95, measure == "per period") |>
  left_join(power |> filter(level == 0.95, measure == "per period", pct %in% c(-0.25, 0.25)) |>
              mutate(pct = sprintf("pw%+.0f", 100 * pct)) |> select(-power_final, -level, -measure) |>
              pivot_wider(names_from = pct, values_from = power_any), by = c("panel", "n_pre", "estimator", "method")) |>
  transmute(panel, n_pre, estimator, method, path = round(100 * path_coverage), final = round(100 * final_coverage),
            hw = round(median_hw_pct, 1), `pw-25` = round(100 * `pw-25`), `pw+25` = round(100 * `pw+25`)) |>
  arrange(panel, desc(n_pre), estimator, method)
marker <- "# ---- Results (written automatically when this script is run) ----"
script_lines <- read_lines(script_file)
cut_at <- which(script_lines == marker)
if (length(cut_at) > 0) script_lines <- script_lines[seq_len(cut_at[1] - 1)]
write_lines(c(script_lines, marker,
  sprintf("# Run: %s | %d tasks (%d failed) | %s min on %d workers", format(Sys.time(), "%Y-%m-%d %H:%M"),
          nrow(grid), nrow(failed), run_minutes, n_workers),
  "# 95%, per period: path / final = whole-path / final-period coverage (%), hw = median half-width (%), pw-25 / pw+25 = power anywhere (%)",
  paste0("# ", knitr::kable(compact, format = "simple"))), script_file)
message("Finished in ", run_minutes, " min.")

# ---- Results (written automatically when this script is run) ----
# Run: 2026-10-03 20:39 | 397 tasks (0 failed) | 26.6 min on 10 workers
# 95%, per period: path / final = whole-path / final-period coverage (%), hw = median half-width (%), pw-25 / pw+25 = power anywhere (%)
# panel         n_pre  estimator       method               path   final     hw   pw-25   pw+25
# -----------  ------  --------------  ------------------  -----  ------  -----  ------  ------
# california       19  current LASSO   placebo_abs            95      97   40.0     100     100
# california       19  demeaned SC     boot_pointwise         92      97   38.1     100      97
# california       19  demeaned SC     boot_simultaneous      95     100   57.2      37      53
# california       19  demeaned SC     placebo_abs            84      97   26.8     100      97
# california       12  current LASSO   placebo_abs            95      95   53.3     100     100
# california       12  demeaned SC     boot_pointwise         92      97   41.6      92      92
# california       12  demeaned SC     boot_simultaneous      97     100   66.3      45      55
# california       12  demeaned SC     placebo_abs            89      97   33.1     100     100
# california        8  current LASSO   placebo_abs            89     100   28.5     100     100
# california        8  demeaned SC     boot_pointwise         95      97   50.8      92      97
# california        8  demeaned SC     boot_simultaneous      97     100   69.2      58      53
# california        8  demeaned SC     placebo_abs            89      97   29.6     100     100
# california        5  current LASSO   placebo_abs            95     100   53.0      97     100
# california        5  demeaned SC     boot_pointwise         84      97   44.7      97     100
# california        5  demeaned SC     boot_simultaneous      97     100   72.9      47      55
# california        5  demeaned SC     placebo_abs            87      97   27.9     100     100
# kansas           80  current LASSO   placebo_abs            88      94   12.1     100     100
# kansas           80  demeaned SC     boot_pointwise         84      96   10.7     100     100
# kansas           80  demeaned SC     boot_simultaneous      96      98   15.5      98      98
# kansas           80  demeaned SC     placebo_abs            88      96   10.5     100     100
# kansas           40  current LASSO   placebo_abs            80      96   10.0     100     100
# kansas           40  demeaned SC     boot_pointwise         80      96    9.3     100     100
# kansas           40  demeaned SC     boot_simultaneous      96      98   19.8     100      94
# kansas           40  demeaned SC     placebo_abs            90      96    9.4     100     100
# kansas           20  current LASSO   placebo_abs            82      94   10.9     100     100
# kansas           20  demeaned SC     boot_pointwise         82      96    9.5     100     100
# kansas           20  demeaned SC     boot_simultaneous      96     100   16.9      98      96
# kansas           20  demeaned SC     placebo_abs            86      96    9.6     100     100
# kansas           12  current LASSO   placebo_abs            90      96   15.6     100     100
# kansas           12  demeaned SC     boot_pointwise         78      94    9.9     100     100
# kansas           12  demeaned SC     boot_simultaneous      96      98   17.7     100      98
# kansas           12  demeaned SC     placebo_abs            86      96   12.2     100     100
# kansas            8  current LASSO   placebo_abs            86      96   17.3     100     100
# kansas            8  demeaned SC     boot_pointwise         82      98   10.8     100     100
# kansas            8  demeaned SC     boot_simultaneous      96      98   21.2     100      96
# kansas            8  demeaned SC     placebo_abs            80      96   10.9     100     100
