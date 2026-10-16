##############################################################################
# simultaneous_band_ablation.R
#
# ABLATION: what does the maximum entropy bootstrap contribute to the simultaneous
# band? The band's construction is held fixed and only the SCALE each gap is
# divided by is changed:
#   bootstrap   the unit's fixed-weight bootstrap half-width h_{u,t} (current method)
#   prefit      the unit's pre-treatment RMSPE (constant over time)
#   relative    the unit's synthetic value (cumulative: its running total), i.e.
#               relative gaps as in the absolute placebo test
#   raw         nothing (s = 1): raw gaps in the outcome's units, as in the classic
#               placebo plots (added 2026-10-04)
# For every scale s, with standardized gaps k_{u,t} = |gap_{u,t}| / s_{u,t}:
#   pointwise     k*_t = m-th largest donor k_t;            half-width = k*_t x s_{T,t}
#   simultaneous  q~_t = loess-smoothed (log scale, span 0.75) median donor k_t;
#                 lambda = m-th largest donor max_t (k_{p,t} / q~_t);
#                 half-width = lambda x q~_t x s_{T,t}
# m = floor(alpha (N + 1)). Same demeaned SC fits, donors and rank rule throughout.
# Check: the bootstrap-scale bands must reproduce dsc_symcal_boot()'s own bands.
#
# Designs (each untreated state as a fake treated unit; donors = the others):
#   kansas      augsynth::kansas, Kansas dropped, data to 2013, fake treatment 2010 Q1
#               (80 pre-treatment quarters, 16 post), 49 states
#   california  tidysynth::smoking, California dropped, fake treatment 1989
#               (19 pre-treatment years, 12 post), 38 states
# Levels 95% and 90%; per-period and cumulative effects. Bands are shift-invariant,
# so coverage is computed at 0% effect and power exactly for injected effects.
#
# Output (results_simultaneous_band_ablation/): tasks/*.rds, results.rds,
#   summary.csv, ablation_table.html (DT, fixed first column and header).
##############################################################################

library(tidyverse)
library(future)
library(furrr)


# ---- 1. Settings --------------------------------------------------------------------
script_file   <- "simulations/simultaneous_band_ablation.R"
functions_src <- "synthetic_control_functions/synthetic_control_functions_claude.R"
out_dir       <- "simulations/results_simultaneous_band_ablation"
tasks_dir     <- file.path(out_dir, "tasks")
n_workers     <- 10
conf_levels   <- c(0.95, 0.90)
effects_pct   <- seq(-0.25, 0.25, by = 0.05)
span          <- 0.75
dir.create(tasks_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 2. Data -------------------------------------------------------------------------
get_design <- function(design) {
  if (design == "kansas") {
    list(data = augsynth::kansas |> filter(state != "Kansas", year <= 2013) |>
                  select(state, year, qtr, gdpcapita) |> as.data.frame(),
         time_cols = c("year", "qtr"), outcome = "gdpcapita", treatment_period = c(year = 2010, qtr = 1))
  } else {
    list(data = tidysynth::smoking |> filter(state != "California") |> select(state, year, cigsale) |> as.data.frame(),
         time_cols = "year", outcome = "cigsale", treatment_period = c(year = 1989))
  }
}
design_units <- function(design) sort(unique(get_design(design)$data$state))


# ---- 3. Bands from a scale -------------------------------------------------------------
kth_largest <- function(x, m) if (m == 0) Inf else sort(x, decreasing = TRUE)[m]
smooth_pos <- function(v) { tt <- seq_along(v); exp(as.numeric(predict(loess(log(v) ~ tt, span = span, degree = 1)))) }

# gap_T, s_T: treated post-period vectors; G, S: post periods x donors. Returns pointwise and simultaneous half-widths.
bands_from_scale <- function(gap_T, s_T, G, S, conf) {
  m <- floor(round(1 - conf, 10) * (ncol(G) + 1))
  K <- abs(G) / S
  k_pw <- apply(K, 1, kth_largest, m = m)
  q <- smooth_pos(apply(K, 1, median))
  lambda <- kth_largest(apply(K / q, 2, max), m)
  tibble(t = seq_along(gap_T), gap0 = gap_T, hw_pw = k_pw * s_T, hw_sim = lambda * q * s_T)
}


# ---- 4. One task: one design, treated unit and level --------------------------------------
sim_one <- function(design, unit, conf) {
  task_file <- file.path(tasks_dir, sprintf("%s_%02.0f_%s.rds", design, 100 * conf, gsub(" ", "_", unit)))
  if (file.exists(task_file)) return(readRDS(task_file))
  if (!exists("dsc_symcal_boot", envir = globalenv())) source(functions_src)

  ds <- get_design(design)
  r <- dsc_symcal_boot(ds$data, "state", ds$time_cols, ds$outcome, unit, treatment_period = ds$treatment_period, conf_level = conf)
  e <- r$effects; post <- e$post
  gaps_all <- as.matrix(r$placebo$gaps); synth_all <- as.matrix(r$placebo$synthetic)

  # Treated unit and donors: per-period and cumulative gaps, and the three scales
  gT <- e$effect[post]; gTc <- cumsum(gT); sT <- e$synthetic[post]
  G <- gaps_all[post, , drop = FALSE]; Gc <- apply(G, 2, cumsum); Sy <- synth_all[post, , drop = FALSE]
  rmse_T <- sqrt(mean(e$effect[!post]^2)); rmse_P <- sqrt(colMeans(gaps_all[!post, , drop = FALSE]^2))
  scales <- list(
    bootstrap = list(pp = list(e$halfwidth[post], as.matrix(r$placebo$halfwidth)),
                     cu = list(e$cum_halfwidth[post], as.matrix(r$placebo$cum_halfwidth))),
    prefit    = list(pp = list(rep(rmse_T, length(gT)), matrix(rmse_P, nrow(G), ncol(G), byrow = TRUE)),
                     cu = list(rep(rmse_T, length(gT)), matrix(rmse_P, nrow(G), ncol(G), byrow = TRUE))),
    relative  = list(pp = list(sT, Sy), cu = list(cumsum(sT), apply(Sy, 2, cumsum))),
    raw       = list(pp = list(rep(1, length(gT)), matrix(1, nrow(G), ncol(G))),
                     cu = list(rep(1, length(gT)), matrix(1, nrow(G), ncol(G)))))

  y0 <- e$observed[post]
  result <- imap_dfr(scales, function(sc, nm) bind_rows(
    bands_from_scale(gT, sc$pp[[1]], G, sc$pp[[2]], conf) |> mutate(measure = "per period", y0 = y0),
    bands_from_scale(gTc, sc$cu[[1]], Gc, sc$cu[[2]], conf) |> mutate(measure = "cumulative", y0 = cumsum(y0))) |>
      mutate(scale = nm)) |>
    mutate(design = design, unit = unit, level = conf)

  # Check: bootstrap-scale bands reproduce the function's own bands
  b <- filter(result, scale == "bootstrap")
  attr(result, "check") <- max(abs(c(b$hw_pw[b$measure == "per period"] - (e$eff_high_pw[post] - e$eff_low_pw[post]) / 2,
                                    b$hw_sim[b$measure == "per period"] - (e$eff_high_sim[post] - e$eff_low_sim[post]) / 2,
                                    b$hw_sim[b$measure == "cumulative"] - (e$cum_eff_high_sim[post] - e$cum_eff_low_sim[post]) / 2)))
  result$check <- attr(result, "check")
  saveRDS(result, task_file)
  result
}


# ---- 5. Grid of tasks, run in parallel with future_pmap -----------------------------------
grid <- map_dfr(c("kansas", "california"), function(dz) expand_grid(design = dz, unit = design_units(dz), conf = conf_levels))
plan(multisession, workers = n_workers)
t0 <- Sys.time()
results <- future_pmap(grid, possibly(sim_one, otherwise = NULL),
                       .options = furrr_options(seed = TRUE, packages = c("tidyverse", "CVXR", "meboot")))
plan(sequential)
run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
failed <- grid[map_lgl(results, is.null), ]
long <- bind_rows(results)


# ---- 6. Summaries (tidyverse) -----------------------------------------------------------------
per_unit <- long |>
  group_by(design, level, scale, measure, unit) |>
  summarise(path_cov_sim = all(abs(gap0) <= hw_sim), final_cov_pw = abs(last(gap0)) <= last(hw_pw),
            hw_sim_pct = median(100 * hw_sim / y0), hw_pw_final_pct = 100 * last(hw_pw) / last(y0), .groups = "drop")
coverage <- per_unit |>
  group_by(design, level, scale, measure) |>
  summarise(n = n(), path_coverage_sim = mean(path_cov_sim), final_coverage_pw = mean(final_cov_pw),
            median_hw_sim_pct = median(hw_sim_pct), median_hw_pw_final_pct = median(hw_pw_final_pct), .groups = "drop")
power <- long |>
  crossing(pct = effects_pct) |>
  mutate(shifted = gap0 + pct * y0) |>
  group_by(design, level, scale, measure, unit, pct) |>
  summarise(any_sim = any(abs(shifted) > hw_sim), final_pw = abs(last(shifted)) > last(hw_pw), .groups = "drop") |>
  group_by(design, level, scale, measure, pct) |>
  summarise(power_any_sim = mean(any_sim), power_final_pw = mean(final_pw), .groups = "drop")


# ---- 7. Save and DT table ------------------------------------------------------------------------
saveRDS(list(settings = list(conf_levels = conf_levels, effects_pct = effects_pct, span = span, n_workers = n_workers,
                             run_minutes = run_minutes, max_check_diff = max(long$check)),
             coverage = coverage, power = power, per_unit = per_unit, failed = failed), file.path(out_dir, "results.rds"))
write_csv(coverage, file.path(out_dir, "summary.csv"))

dt_tab <- coverage |>
  left_join(power |> filter(round(abs(pct), 2) %in% c(0.10, 0.25)) |>
              mutate(col = sprintf("power %+.0f%% (sim anywhere / pw final)", 100 * pct),
                     val = sprintf("%.0f / %.0f", 100 * power_any_sim, 100 * power_final_pw)) |>
              select(design, level, scale, measure, col, val) |> pivot_wider(names_from = col, values_from = val),
            by = c("design", "level", "scale", "measure")) |>
  mutate(Method = sprintf("%s | %.0f%% | %s | scale: %s", design, 100 * level, measure, scale),
         `whole-path coverage (sim)` = round(100 * path_coverage_sim), `final coverage (pw)` = round(100 * final_coverage_pw),
         `median half-width, sim (%)` = round(median_hw_sim_pct, 1), `final half-width, pw (%)` = round(median_hw_pw_final_pct, 1)) |>
  arrange(design, desc(level), measure, scale) |>
  select(Method, `whole-path coverage (sim)`, `final coverage (pw)`, `median half-width, sim (%)`, `final half-width, pw (%)`, starts_with("power"))
widget <- DT::datatable(dt_tab, rownames = FALSE, extensions = "FixedColumns", filter = "top",
  caption = "Simultaneous-band ablation: the same band construction with three scales (bootstrap, pre-treatment RMSPE, relative level). Values in %. Power: simultaneous band anywhere in the post period / pointwise interval in the final period.",
  options = list(scrollX = TRUE, scrollY = "70vh", scrollCollapse = TRUE, paging = FALSE, fixedColumns = list(leftColumns = 1)))
tmp_html <- file.path(tempdir(), "ablation_table.html")
htmlwidgets::saveWidget(widget, tmp_html, selfcontained = TRUE, title = "Simultaneous band ablation")
file.copy(tmp_html, file.path(out_dir, "ablation_table.html"), overwrite = TRUE)

marker <- "# ---- Results (written automatically when this script is run) ----"
script_lines <- read_lines(script_file)
cut_at <- which(script_lines == marker)
if (length(cut_at) > 0) script_lines <- script_lines[seq_len(cut_at[1] - 1)]
write_lines(c(script_lines, marker,
  sprintf("# Run: %s | %d tasks (%d failed) | %s min on %d workers | bootstrap-scale check max diff: %.3g",
          format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(grid), nrow(failed), run_minutes, n_workers, max(long$check)),
  paste0("# ", knitr::kable(dt_tab |> select(1:5), format = "simple"))), script_file)
message("Finished in ", run_minutes, " min.")

# ---- Results (written automatically when this script is run) ----
# Run: 2026-10-04 14:01 | 174 tasks (0 failed) | 5 min on 10 workers | bootstrap-scale check max diff: 2.91e-11
# Method                                              whole-path coverage (sim)   final coverage (pw)   median half-width, sim (%)   final half-width, pw (%)
# -------------------------------------------------  --------------------------  --------------------  ---------------------------  -------------------------
# california | 95% | cumulative | scale: bootstrap                           97                    97                         38.1                       38.4
# california | 95% | cumulative | scale: prefit                              97                    97                         41.5                       44.8
# california | 95% | cumulative | scale: raw                                 97                    97                         63.6                       26.9
# california | 95% | cumulative | scale: relative                            97                    97                         26.7                       22.8
# california | 95% | per period | scale: bootstrap                           95                    97                         54.0                       39.4
# california | 95% | per period | scale: prefit                              97                    97                         86.0                       72.4
# california | 95% | per period | scale: raw                                 97                    97                         92.5                       40.6
# california | 95% | per period | scale: relative                            97                    97                         46.5                       35.1
# california | 90% | cumulative | scale: bootstrap                           95                    95                         16.4                       13.7
# california | 90% | cumulative | scale: prefit                              92                    92                         20.5                       19.4
# california | 90% | cumulative | scale: raw                                 92                    92                         22.6                       17.2
# california | 90% | cumulative | scale: relative                            92                    95                         16.1                       17.0
# california | 90% | per period | scale: bootstrap                           89                    92                         25.1                       28.9
# california | 90% | per period | scale: prefit                              92                    92                         34.9                       28.4
# california | 90% | per period | scale: raw                                 92                    95                         37.5                       27.3
# california | 90% | per period | scale: relative                            92                    92                         32.5                       26.6
# kansas | 95% | cumulative | scale: bootstrap                               98                    96                          9.5                        9.6
# kansas | 95% | cumulative | scale: prefit                                  96                    96                          7.7                        8.3
# kansas | 95% | cumulative | scale: raw                                     96                    96                         18.4                       13.3
# kansas | 95% | cumulative | scale: relative                                96                    96                         10.4                        9.0
# kansas | 95% | per period | scale: bootstrap                               96                    96                         15.8                       11.1
# kansas | 95% | per period | scale: prefit                                  96                    96                         13.6                       14.4
# kansas | 95% | per period | scale: raw                                     96                    96                         21.9                       13.2
# kansas | 95% | per period | scale: relative                                94                    96                         15.0                       11.3
# kansas | 90% | cumulative | scale: bootstrap                               94                    92                          7.2                        6.2
# kansas | 90% | cumulative | scale: prefit                                  92                    92                          6.2                        6.8
# kansas | 90% | cumulative | scale: raw                                     92                    90                         10.8                        8.9
# kansas | 90% | cumulative | scale: relative                                90                    92                          9.4                        7.3
# kansas | 90% | per period | scale: bootstrap                               92                    92                         13.0                        9.2
# kansas | 90% | per period | scale: prefit                                  92                    92                         12.4                       10.1
# kansas | 90% | per period | scale: raw                                     92                    92                         18.2                       12.9
# kansas | 90% | per period | scale: relative                                94                    94                         13.1                       10.0
