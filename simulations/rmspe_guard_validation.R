# rmspe_guard_validation.R
#
# Validity and power of dsc_placebo() on the RMSPE scale with the guard
# (rmspe_floor: every unit's pre-treatment RMSPE floored at rmspe_floor x the
# median over all units) at floor = 0 (no guard), 0.1 and 0.25, with the relative
# scale as a reference.
#
# Why: with short pre-periods and many donors, the demeaned synthetic control can
# reproduce some placebos exactly before treatment (pre-RMSPE ~ 0). Their
# standardized gaps then become near-infinite and so do the RMSPE-scale intervals
# (found on the Texas prison panel: 8 pre-years, 49 donors).
#
# Designs (true effect 0; power by shifting, since all intervals shift exactly
# with an added effect):
#   Real panels, every untreated unit in turn as a fake treated unit:
#     Kansas GDP per capita (augsynth::kansas), pre-period cut to 80 / 20 / 8 quarters, 16 post
#     California cigarette sales (tidysynth::smoking), pre-period 19 / 8 / 5 years, 12 post
#     Texas prison panel (Cunningham's Mixtape, texas.dta: black male prisoners), 8 pre, 8 post;
#       Texas itself excluded; downloaded from GitHub
#   Factor model (as in augmented_comparison_simulation.R): treated unit exchangeable /
#     edge / outside, stationary or random-walk factors, T0 = 10 / 20 / 40, 40 donors,
#     12 post periods, 100 replicates per cell
#
# Each panel is fitted once with dsc_placebo(); the calibration is then applied
# for each scale / floor by calibrate_floor(), which reproduces dsc_placebo()'s
# calibration exactly (checked at the top of the run).
#
# Run from the project root:  source("simulations/rmspe_guard_validation.R")
# Uses future + furrr::future_pmap (10 workers); resumes from saved tasks.

suppressPackageStartupMessages({
  library(tidyverse)
  library(future)
  library(furrr)
})

functions_src <- "synthetic_control_functions/synthetic_control_functions_claude.R"
source(functions_src)
out_dir   <- file.path("simulations", "results_rmspe_guard")
tasks_dir <- file.path(out_dir, "tasks")
dir.create(tasks_dir, recursive = TRUE, showWarnings = FALSE)
conf <- 0.95; span <- 0.75; n_workers <- 10
methods <- tribble(~method, ~scale, ~floor,
                   "rmspe, no guard", "rmspe", 0,
                   "rmspe, floor 0.1", "rmspe", 0.1,
                   "rmspe, floor 0.25", "rmspe", 0.25,
                   "relative", "relative", 0)

# ---- calibration with a floor (same rules as dsc_placebo) ------------------------------------
calibrate_floor <- function(gT, sT, GP, SP, post, scale, floor) {
  n_post <- sum(post); N <- ncol(GP); m <- floor(round(1 - conf, 10) * (N + 1))
  kth <- function(x) if (m == 0) Inf else sort(x, decreasing = TRUE)[m]
  smooth_pos <- function(v) { tt <- seq_along(v); exp(as.numeric(predict(loess(log(v) ~ tt, span = span, degree = 1)))) }
  rmspe_T <- sqrt(mean(gT[!post]^2)); rmspe_P <- sqrt(colMeans(GP[!post, , drop = FALSE]^2))
  fv <- if (scale == "rmspe" && floor > 0) floor * median(c(rmspe_T, rmspe_P)) else 0
  s_T <- max(rmspe_T, fv); s_P <- pmax(rmspe_P, fv)
  ones <- matrix(1, n_post, N)
  g <- gT[post]; G <- GP[post, , drop = FALSE]; Gc <- apply(G, 2, cumsum)
  sc_T  <- switch(scale, relative = sT[post], rmspe = rep(s_T, n_post))
  sc_Tc <- switch(scale, relative = cumsum(sT[post]), rmspe = rep(s_T, n_post))
  sc_P  <- switch(scale, relative = SP[post, , drop = FALSE], rmspe = sweep(ones, 2, s_P, "*"))
  sc_Pc <- switch(scale, relative = apply(SP[post, , drop = FALSE], 2, cumsum), rmspe = sweep(ones, 2, s_P, "*"))
  if (scale == "relative" && (any(sT[post] <= 0) || any(SP[post, ] <= 0))) stop("relative scale: non-positive synthetic values")
  one <- function(K, est, s) {
    q <- smooth_pos(apply(K, 1, median)); lam <- kth(apply(K / q, 2, max))
    tibble(t = seq_len(n_post), est = est, hw_pw = apply(K, 1, kth) * s, hw_sim = lam * q * s)
  }
  bind_rows(one(abs(G) / sc_P, g, sc_T) |> mutate(measure = "per period"),
            one(abs(Gc) / sc_Pc, cumsum(g), sc_Tc) |> mutate(measure = "cumulative")) |>
    mutate(n_floored = sum(c(rmspe_T, rmspe_P) < fv))
}

calibrate_all <- function(fit, y0) {
  e <- fit$effects; post <- e$post
  GP <- as.matrix(fit$placebo$gaps); SP <- as.matrix(fit$placebo$synthetic)
  pmap_dfr(methods, function(method, scale, floor)
    tryCatch(calibrate_floor(e$effect, e$synthetic, GP, SP, post, scale, floor) |>
               mutate(method = method, y0 = if_else(measure == "per period", y0[t], cumsum(y0)[t])),
             error = function(err) tibble(method = method, error = conditionMessage(err))))
}

# ---- check: calibrate_floor() reproduces dsc_placebo() --------------------------------------
local({
  ks <- as.data.frame(augsynth::kansas[, c("state", "year", "qtr", "gdpcapita")])
  for (fl in c(0, 0.25)) {
    f <- dsc_placebo(ks, "state", c("year", "qtr"), "gdpcapita", "Iowa", treatment_period = c(year = 2012, qtr = 2),
                     diagnostics = FALSE, rmspe_floor = fl)
    e <- f$effects[f$effects$post, ]
    c1 <- calibrate_floor(f$effects$effect, f$effects$synthetic, as.matrix(f$placebo$gaps), as.matrix(f$placebo$synthetic),
                          f$effects$post, "rmspe", fl)
    pp <- filter(c1, measure == "per period"); cu <- filter(c1, measure == "cumulative")
    stopifnot(isTRUE(all.equal(e$effect + pp$hw_pw, e$eff_high_pw)), isTRUE(all.equal(e$effect + pp$hw_sim, e$eff_high_sim)),
              isTRUE(all.equal(e$cum_effect + cu$hw_sim, e$cum_eff_high_sim)))
  }
  message("calibrate_floor() matches dsc_placebo() (floor 0 and 0.25)")
})

# ---- real panels -------------------------------------------------------------------------------
kansas <- as.data.frame(augsynth::kansas[, c("state", "year", "qtr", "gdpcapita")]) |>
  mutate(t = (year - min(year)) * 4 + qtr) |> select(unit = state, t, y = gdpcapita)
kansas_treat <- (2012 - min(augsynth::kansas$year)) * 4 + 2                     # 2012 Q2
california <- as.data.frame(tidysynth::smoking)[, c("state", "year", "cigsale")] |> select(unit = state, t = year, y = cigsale)
texas_file <- file.path(out_dir, "texas.rds")
if (!file.exists(texas_file)) {
  tx <- haven::zap_labels(haven::read_dta("https://github.com/scunning1975/mixtape/raw/master/texas.dta"))
  saveRDS(as.data.frame(tx)[, c("state", "year", "bmprison")] |> mutate(state = as.character(state)), texas_file)
}
texas <- readRDS(texas_file) |> filter(state != "Texas") |> select(unit = state, t = year, y = bmprison)

panels <- list(Kansas = list(d = filter(kansas, unit != "Kansas"), treat = kansas_treat, T0 = c(80, 20, 8), effect = 0.10),
               California = list(d = filter(california, unit != "California"), treat = 1989, T0 = c(19, 8, 5), effect = 0.25),
               Texas = list(d = texas, treat = 1993, T0 = 8, effect = 0.25))
real_grid <- imap_dfr(panels, function(p, nm) expand_grid(panel = nm, T0 = p$T0, unit = sort(unique(p$d$unit))))

run_real <- function(panel, T0, unit) {
  task_file <- file.path(tasks_dir, sprintf("real_%s_%02d_%s.rds", panel, T0, gsub(" ", "_", unit)))
  if (file.exists(task_file)) return(readRDS(task_file))
  if (!exists("dsc_placebo", envir = globalenv())) source(functions_src)
  p <- panels[[panel]]
  d <- p$d |> filter(t >= p$treat - T0) |> as.data.frame()
  res <- tryCatch({
    fit <- dsc_placebo(d, "unit", "t", "y", unit, treatment_period = c(t = p$treat), diagnostics = FALSE)
    y0 <- d |> filter(unit == !!unit, t >= p$treat) |> arrange(t) |> pull(y)
    calibrate_all(fit, y0)
  }, error = function(err) tibble(error = conditionMessage(err)))
  res <- res |> mutate(design = "real panel", panel = panel, T0 = T0, unit = unit, effect_size = p$effect)
  saveRDS(res, task_file)
  res
}

# ---- factor model ------------------------------------------------------------------------------
T1 <- 12
simulate_panel <- function(treated, factors, n_donors, T0, rep_seed) {      # as in augmented_comparison_simulation.R
  set.seed(rep_seed)
  TT <- T0 + T1; N <- n_donors + 1; K <- 2
  f <- sapply(seq_len(K), function(k) {
    u <- if (factors == "stationary") rnorm(TT, 0, 0.05) else rnorm(TT, 0.005, 0.03)
    if (factors == "stationary") as.numeric(stats::filter(u, 0.8, method = "recursive")) else cumsum(u)
  })
  lam <- matrix(rnorm(N * K), N, K)
  if (treated == "edge")    lam[1, ] <- apply(lam[-1, , drop = FALSE], 2, max)
  if (treated == "outside") lam[1, ] <- apply(lam[-1, , drop = FALSE], 2, max) + 1
  e <- sapply(seq_len(N), function(i) as.numeric(stats::filter(rnorm(TT, 0, 0.02 * sqrt(1 - 0.5^2)), 0.5, method = "recursive")))
  alpha <- rnorm(N, 0, 0.3)
  logy <- sapply(seq_len(N), function(i) 5 + alpha[i] + 0.01 * seq_len(TT) + f %*% lam[i, ] + e[, i])
  tibble(unit = rep(sprintf("u%02d", seq_len(N)), each = TT), period = rep(seq_len(TT), N), y = exp(as.vector(logy)))
}
fm_grid <- expand_grid(treated = c("exchangeable", "edge", "outside"), factors = c("stationary", "random_walk"),
                       T0 = c(10, 20, 40)) |> mutate(cell_id = row_number())
reps <- 100; reps_per_task <- 5; seed <- 20261015
fm_tasks <- expand_grid(cell_id = fm_grid$cell_id, chunk = seq_len(reps / reps_per_task))

run_fm <- function(cell_id, chunk) {
  task_file <- file.path(tasks_dir, sprintf("fm_cell%02d_chunk%02d.rds", cell_id, chunk))
  if (file.exists(task_file)) return(readRDS(task_file))
  if (!exists("dsc_placebo", envir = globalenv())) source(functions_src)
  cell <- fm_grid[fm_grid$cell_id == cell_id, ]
  res <- map_dfr(((chunk - 1) * reps_per_task + 1):(chunk * reps_per_task), function(r) {
    d <- simulate_panel(cell$treated, cell$factors, 40, cell$T0, rep_seed = seed + 100000 * cell_id + r)
    tryCatch({
      fit <- dsc_placebo(as.data.frame(d), "unit", "period", "y", "u01", treatment_period = c(period = cell$T0 + 1), diagnostics = FALSE)
      y0 <- d |> filter(unit == "u01", period > cell$T0) |> arrange(period) |> pull(y)
      calibrate_all(fit, y0) |> mutate(rep = r)
    }, error = function(err) tibble(error = conditionMessage(err), rep = r))
  })
  res <- res |> mutate(design = "factor model", treated = cell$treated, factors = cell$factors, T0 = cell$T0, effect_size = 0.10)
  saveRDS(res, task_file)
  res
}

# ---- run -------------------------------------------------------------------------------------------
t_start <- Sys.time()
plan(multisession, workers = n_workers)
opts <- furrr_options(seed = TRUE, packages = c("tidyverse", "CVXR"))
real <- future_pmap(real_grid, run_real, .options = opts, .progress = TRUE) |> bind_rows()
fm   <- future_pmap(fm_tasks, run_fm, .options = opts, .progress = TRUE) |> bind_rows()
plan(sequential)
run_minutes <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

all_rows <- bind_rows(real, fm)
failed <- all_rows |> filter(!is.na(error)) |> count(design, method, error)
rows <- all_rows |> filter(is.na(error)) |>
  mutate(case = if_else(design == "real panel", paste0(panel, ", T0 = ", T0), paste0(treated, ", ", factors, ", T0 = ", T0)),
         id = if_else(design == "real panel", unit, as.character(rep)))

# ---- evaluate: true effect 0 and +/- effect_size x untreated outcome (exact shift) ------------
evaluate <- function(r, delta) r |>
  mutate(tau = delta * y0, est_d = est + tau) |>
  group_by(design, case, method, id, measure) |> arrange(t, .by_group = TRUE) |>
  summarise(cover_last = abs(last(est)) <= last(hw_pw),            # interval at the last period contains the true effect
            cover_sim  = all(abs(est) <= hw_sim),                   # band contains the true effect in every period
            reject_last = abs(last(est_d)) > last(hw_pw),
            reject_sim  = any(abs(est_d) > hw_sim),
            hw_last_pct = 100 * last(hw_pw) / abs(last(y0)),        # half-width, % of the untreated outcome
            n_floored = first(n_floored), .groups = "drop") |>
  mutate(delta = delta)

ev <- rows |> group_split(effect_size) |>
  map_dfr(function(r) map_dfr(c(0, -1, 1) * r$effect_size[1], ~ evaluate(r, .x)))

summary_tbl <- ev |>
  group_by(design, case, measure, method) |>
  summarise(n = sum(delta == 0),
            coverage_last = 100 * mean(cover_last[delta == 0]), coverage_sim = 100 * mean(cover_sim[delta == 0]),
            false_pos_last = 100 * mean(reject_last[delta == 0]), false_pos_sim = 100 * mean(reject_sim[delta == 0]),
            power_last_neg = 100 * mean(reject_last[delta < 0]), power_last_pos = 100 * mean(reject_last[delta > 0]),
            power_sim_neg = 100 * mean(reject_sim[delta < 0]), power_sim_pos = 100 * mean(reject_sim[delta > 0]),
            median_hw_last_pct = median(hw_last_pct[delta == 0]),
            share_hw_over_100pct = 100 * mean(hw_last_pct[delta == 0] > 100),
            mean_units_floored = mean(n_floored[delta == 0]), .groups = "drop") |>
  mutate(method = factor(method, methods$method)) |>
  arrange(design, case, measure, method)

saveRDS(list(settings = list(conf = conf, span = span, methods = methods, panels = map(panels, ~ .x[c("treat", "T0", "effect")]),
                             fm_grid = fm_grid, reps = reps, seed = seed, run_minutes = run_minutes),
             summary = summary_tbl, failed = failed), file.path(out_dir, "results.rds"))
write_csv(summary_tbl, file.path(out_dir, "summary.csv"))
if (requireNamespace("DT", quietly = TRUE)) {
  DT::saveWidget(DT::datatable(mutate(summary_tbl, across(where(is.numeric), ~ round(.x, 1))),
                               caption = "RMSPE guard validation (95%): coverage, false positives and power (%), median half-width (% of untreated outcome)",
                               filter = "top", options = list(pageLength = 40)),
                 normalizePath(file.path(out_dir, "rmspe_guard_table.html"), mustWork = FALSE), selfcontained = TRUE)
}
message(sprintf("Done in %.1f minutes", run_minutes))
