##############################################################################
# augmented_comparison_simulation.R
#
# Compares, on the same simulated panels (linear factor model, as in
# factor_model_simulation.R):
#   dsc              demeaned synthetic control, placebo calibration (RMSPE / relative scale)
#   adsc             AUGMENTED demeaned synthetic control (dsc_placebo(augment = "ife")),
#                    placebo calibration (RMSPE / relative scale)
#   fect_parametric  fect IFEct, fect's own parametric bootstrap
#   fect_placebo     fect IFEct counterfactuals with placebo calibration
#                    (every donor refitted as a pseudo-treated unit; RMSPE / relative)
#   ascm_fe          augsynth ridge augmented SC with fixedeff = TRUE, conformal inference
#   ascm             augsynth standard ridge augmented SC, conformal inference
# Treated unit exchangeable, at the edge, or outside the donors' loading range;
# stationary or random-walk factors; T0 = 20 / 40 / 80; 40 donors; 12 post periods.
#
# Each replicate is run with no effect (Option A). Coverage, bias and false positives
# are computed at no effect for all methods. Power is computed by shifting for the
# placebo-calibrated methods only (dsc, adsc, fect_placebo), whose intervals shift
# exactly with an added effect; fect's parametric and augsynth's conformal intervals
# do not (checked before the run), so their power is not estimated here.
#
# Output (results_augmented_comparison/): tasks/*.rds, results.rds, summary.csv,
#   augmented_comparison_table.html (DT, fixed first column and header).
##############################################################################

library(tidyverse)
library(future)
library(furrr)


# ---- 1. Settings ---------------------------------------------------------------------
functions_src <- "synthetic_control_functions/synthetic_control_functions_claude.R"
out_dir       <- "simulations/results_augmented_comparison"
tasks_dir     <- file.path(out_dir, "tasks")
n_workers     <- 10
reps_per_task <- 5
conf          <- 0.95
T1            <- 12
span          <- 0.75
effects_pct   <- c(-0.20, -0.10, -0.05, 0.05, 0.10, 0.20)
seed          <- 20261014
grid <- expand_grid(treated = c("exchangeable", "edge", "outside"), factors = c("stationary", "random_walk"),
                    n_donors = 40, T0 = c(20, 40, 80), reps = 100) |> mutate(cell_id = row_number())
dir.create(tasks_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 2. Data generator (as in factor_model_simulation.R, no covariate) ----------------------
simulate_panel <- function(treated, factors, n_donors, T0, rep_seed) {
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


# ---- 3. Placebo calibration for any estimator's gaps (same rules as dsc_placebo) -------------
# gT, sT: treated gap and counterfactual (all periods); GP, SP: placebo gaps and counterfactuals (periods x donors)
calibrate <- function(gT, sT, GP, SP, post, scale) {
  n_post <- sum(post); N <- ncol(GP); m <- floor(round(1 - conf, 10) * (N + 1))
  kth <- function(x) if (m == 0) Inf else sort(x, decreasing = TRUE)[m]
  smooth_pos <- function(v) { tt <- seq_along(v); exp(as.numeric(predict(loess(log(v) ~ tt, span = span, degree = 1)))) }
  rmspe_T <- sqrt(mean(gT[!post]^2)); rmspe_P <- sqrt(colMeans(GP[!post, , drop = FALSE]^2))
  ones <- matrix(1, n_post, N)
  g <- gT[post]; G <- GP[post, , drop = FALSE]; Gc <- apply(G, 2, cumsum)
  sc_T  <- switch(scale, relative = sT[post], rmspe = rep(rmspe_T, n_post))
  sc_Tc <- switch(scale, relative = cumsum(sT[post]), rmspe = rep(rmspe_T, n_post))
  sc_P  <- switch(scale, relative = SP[post, , drop = FALSE], rmspe = sweep(ones, 2, rmspe_P, "*"))
  sc_Pc <- switch(scale, relative = apply(SP[post, , drop = FALSE], 2, cumsum), rmspe = sweep(ones, 2, rmspe_P, "*"))
  one <- function(K, kT, est, s) {
    q <- smooth_pos(apply(K, 1, median)); lam <- kth(apply(K / q, 2, max))
    hw_pw <- apply(K, 1, kth) * s; hw_sim <- lam * q * s
    tibble(t = seq_len(n_post), est = est, lo_pw = est - hw_pw, hi_pw = est + hw_pw, lo_sim = est - hw_sim, hi_sim = est + hw_sim)
  }
  bind_rows(one(abs(G) / sc_P, g / sc_T, g, sc_T) |> mutate(measure = "per period"),
            one(abs(Gc) / sc_Pc, cumsum(g) / sc_Tc, cumsum(g), sc_Tc) |> mutate(measure = "cumulative"))
}
from_dsc <- function(fit, method) {                                  # dsc_placebo() output -> both scales
  e <- fit$effects; post <- e$post
  GP <- as.matrix(fit$placebo$gaps); SP <- as.matrix(fit$placebo$synthetic)
  map_dfr(c("rmspe", "relative"), ~ calibrate(e$effect, e$synthetic, GP, SP, post, .x) |> mutate(method = paste0(method, "_", .x)))
}


# ---- 4. fect: IFEct with its own inference, and with placebo calibration ----------------------
fect_fit <- function(dd, treated, r, cv, se) {
  dd$D <- as.integer(dd$unit == treated & dd$post)
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(fect::fect(
    y ~ D, data = as.data.frame(dd), index = c("unit", "period"), method = "ife", r = r, CV = cv, force = "two-way",
    se = se, vartype = "parametric", nboots = 100, time.component.from = "nevertreated", parallel = FALSE))))
  out
}
run_fect <- function(d, T0) {
  dd <- d |> mutate(post = period > T0); post <- sort(unique(dd$period)) > T0
  main <- fect_fit(dd, "u01", r = c(0, 5), cv = TRUE, se = TRUE)
  r_use <- as.integer(main$r.cv)
  ids <- as.character(main$id); Yd <- main$Y.dat; colnames(Yd) <- ids
  Yct <- if (is.matrix(main$Y.ct.full)) main$Y.ct.full else main$Y.ct; colnames(Yct) <- ids
  gT <- Yd[, "u01"] - Yct[, "u01"]; sT <- Yct[, "u01"]
  # fect's own inference: per-period intervals (relative periods 1..T1) and the average effect (x T1 for the cumulative total)
  ea <- as.data.frame(main$est.att); ea$rel <- as.numeric(rownames(ea)); ea <- ea[ea$rel >= 1, ]
  own <- bind_rows(
    tibble(t = ea$rel, est = ea$ATT, lo_pw = ea$CI.lower, hi_pw = ea$CI.upper, measure = "per period"),
    tibble(t = T1, est = T1 * main$est.avg[1, "ATT.avg"], lo_pw = T1 * main$est.avg[1, "CI.lower"], hi_pw = T1 * main$est.avg[1, "CI.upper"],
           measure = "cumulative")) |>
    mutate(lo_sim = NA_real_, hi_sim = NA_real_, method = "fect_parametric", p_whole = main$est.avg[1, "p.value"])
  # placebo calibration: each donor as a pseudo-treated unit (treated unit dropped), same number of factors
  donors <- setdiff(ids, "u01")
  pl <- map(donors, function(p) {
    fp <- fect_fit(dd |> filter(unit != "u01"), p, r = r_use, cv = FALSE, se = FALSE)
    idp <- as.character(fp$id); Ydp <- fp$Y.dat; colnames(Ydp) <- idp
    Yp <- if (is.matrix(fp$Y.ct.full)) fp$Y.ct.full else fp$Y.ct; colnames(Yp) <- idp
    list(gap = Ydp[, p] - Yp[, p], synth = Yp[, p])
  })
  GP <- sapply(pl, `[[`, "gap"); SP <- sapply(pl, `[[`, "synth")
  plac <- map_dfr(c("rmspe", "relative"), ~ calibrate(gT, sT, GP, SP, post, .x) |> mutate(method = paste0("fect_placebo_", .x)))
  list(rows = bind_rows(own, plac), y0 = Yd[post, "u01"], n_factors = r_use)
}


# ---- 5. augsynth: augmented SC with conformal inference -------------------------------------------
run_ascm <- function(d, T0, fixedeff) {
  dd <- d |> mutate(trt = as.integer(unit == "u01" & period > T0)) |> as.data.frame()
  fit <- NULL; sm <- NULL
  utils::capture.output({
    fit <- suppressMessages(suppressWarnings(augsynth::augsynth(y ~ trt, unit = unit, time = period, data = dd,
                                                                 progfunc = "Ridge", scm = TRUE, fixedeff = fixedeff)))
    sm <- suppressMessages(suppressWarnings(summary(fit, inf_type = "conformal")))
  })
  a <- sm$att |> filter(Time > T0)
  tibble(t = seq_len(nrow(a)), est = a$Estimate, lo_pw = a$lower_bound, hi_pw = a$upper_bound, lo_sim = NA_real_, hi_sim = NA_real_,
         measure = "per period", method = if (fixedeff) "ascm_fe" else "ascm", p_whole = sm$average_att$p_val[1])
}


# ---- 6. One replicate ------------------------------------------------------------------------------
run_rep <- function(cell, rep) {
  d <- simulate_panel(cell$treated, cell$factors, cell$n_donors, cell$T0, rep_seed = seed + 100000 * cell$cell_id + rep)
  tp <- c(period = cell$T0 + 1)
  y0 <- d |> filter(unit == "u01", period > cell$T0) |> arrange(period) |> pull(y)
  safe <- function(name, expr) tryCatch(expr, error = function(e) tibble(method = name, error = conditionMessage(e)))
  out <- bind_rows(
    safe("dsc", { f <- dsc_placebo(as.data.frame(d), "unit", "period", "y", "u01", treatment_period = tp, scale = "rmspe")
                  from_dsc(f, "dsc") |> mutate(p_whole = f$test$p_value[f$test$test == "RMSPE ratio (Abadie et al.)"]) }),
    safe("adsc", { f <- dsc_placebo(as.data.frame(d), "unit", "period", "y", "u01", treatment_period = tp, scale = "rmspe", augment = "ife")
                   from_dsc(f, "adsc") |> mutate(p_whole = f$test$p_value[f$test$test == "RMSPE ratio (Abadie et al.)"],
                                                 n_factors = f$adjustment$n_factors) }),
    safe("fect", { rf <- run_fect(d, cell$T0); rf$rows |> mutate(n_factors = rf$n_factors) }),
    safe("ascm_fe", run_ascm(d, cell$T0, fixedeff = TRUE)),
    safe("ascm", run_ascm(d, cell$T0, fixedeff = FALSE)))
  out |>
    mutate(y0 = if_else(measure == "per period", y0[pmin(t, length(y0))], cumsum(y0)[pmin(t, length(y0))]), rep = rep)
}

run_task <- function(cell_id, chunk) {
  task_file <- file.path(tasks_dir, sprintf("cell%02d_chunk%02d.rds", cell_id, chunk))
  if (file.exists(task_file)) return(readRDS(task_file))
  if (!exists("dsc_placebo", envir = globalenv())) source(functions_src)
  cell <- grid |> filter(cell_id == !!cell_id) |> as.list()
  reps <- ((chunk - 1) * reps_per_task + 1):min(chunk * reps_per_task, cell$reps)
  res <- map_dfr(reps, ~ run_rep(cell, .x)) |> mutate(cell_id = cell_id)
  saveRDS(res, task_file)
  res
}


# ---- 7. Run (skipped when dry_run is TRUE) -----------------------------------------------------------
tasks <- grid |> select(cell_id, reps) |> mutate(chunk = map(reps, ~ seq_len(ceiling(.x / reps_per_task)))) |>
  unnest(chunk) |> select(cell_id, chunk)
if (!exists("dry_run") || !isTRUE(dry_run)) {
  plan(multisession, workers = n_workers)
  t0 <- Sys.time()
  results <- future_pmap(tasks, possibly(run_task, otherwise = NULL),
                         .options = furrr_options(seed = TRUE, packages = c("tidyverse", "CVXR", "fect", "augsynth")))
  plan(sequential)
  run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  long <- bind_rows(results)
  if (!"error" %in% names(long)) long$error <- NA_character_
  errors <- long |> filter(!is.na(error)) |> distinct(cell_id, method, rep, error)
  long <- long |> filter(is.na(error)) |> left_join(grid, by = "cell_id")

  # ---- 8. Summaries ---------------------------------------------------------------------------------
  final <- long |> group_by(cell_id, method, measure, rep) |> filter(t == max(t)) |> ungroup()
  accuracy <- final |> mutate(err = 100 * est / y0) |>
    group_by(cell_id, method, measure) |> summarise(bias = mean(err), rmse = sqrt(mean(err^2)), .groups = "drop")
  coverage <- long |> group_by(cell_id, method, measure, rep) |>
    summarise(final_cov = { k <- which.max(t); lo_pw[k] <= 0 & hi_pw[k] >= 0 },
              path_cov = if (all(is.na(lo_sim))) NA else all(lo_sim <= 0 & hi_sim >= 0),
              hw_final_pct = { k <- which.max(t); 100 * (hi_pw[k] - lo_pw[k]) / 2 / y0[k] },
              p_whole = first(p_whole), .groups = "drop") |>
    group_by(cell_id, method, measure) |>
    summarise(n = n(), final_coverage = mean(final_cov, na.rm = TRUE), path_coverage = mean(path_cov),
              median_hw_pct = median(hw_final_pct, na.rm = TRUE),
              whole_false_pos = mean(p_whole <= 1 - conf + 1e-12, na.rm = TRUE), .groups = "drop")
  shifts <- bind_rows(tibble(effect = sprintf("%+.0f%%", 100 * effects_pct), pct = effects_pct, shape = "constant"),
                      tibble(effect = "ramp to -10%", pct = -0.10, shape = "ramp"))
  # Power by shifting is valid only for the shift-invariant (placebo-calibrated) methods; fect's parametric and
  # augsynth's conformal intervals change shape when an effect is added (checked before the run), so they are excluded.
  power <- long |> filter(str_detect(method, "^(dsc|adsc|fect_placebo)_")) |>
    crossing(shifts) |> arrange(cell_id, method, measure, rep, effect, t) |>
    group_by(cell_id, method, measure, rep, effect) |>
    mutate(w = if_else(shape == "constant", 1, t / max(t)),
           step = pct * w * if_else(measure == "per period", y0, y0 - lag(y0, default = 0)),
           s = if_else(measure == "per period", step, cumsum(step))) |>
    summarise(final_det = { k <- which.max(t); lo_pw[k] + s[k] > 0 | hi_pw[k] + s[k] < 0 },
              any_det = if (all(is.na(lo_sim))) NA else any(lo_sim + s > 0 | hi_sim + s < 0), .groups = "drop") |>
    group_by(cell_id, method, measure, effect) |>
    summarise(power_final = mean(final_det, na.rm = TRUE), power_any = mean(any_det), .groups = "drop")
  failures <- errors |> count(cell_id, method, name = "failed")

  saveRDS(list(settings = list(grid = grid, conf = conf, T1 = T1, run_minutes = run_minutes, seed = seed),
               accuracy = accuracy, coverage = coverage, power = power, errors = errors, failures = failures),
          file.path(out_dir, "results.rds"))
  summary_tbl <- grid |> select(-reps) |> inner_join(accuracy, by = "cell_id") |>
    left_join(coverage, by = c("cell_id", "method", "measure"))
  write_csv(summary_tbl, file.path(out_dir, "summary.csv"))
  dt_tab <- summary_tbl |>
    left_join(power |> filter(effect %in% c("-10%", "+10%")) |>
                mutate(val = sprintf("%.0f / %s", 100 * power_final, ifelse(is.na(power_any), "-", sprintf("%.0f", 100 * power_any)))) |>
                select(cell_id, method, measure, effect, val) |> pivot_wider(names_from = effect, values_from = val, names_prefix = "power "),
              by = c("cell_id", "method", "measure")) |>
    transmute(Setting = sprintf("%s | %s | T0=%d | %s | %s", treated, factors, T0, method, measure),
              `bias (%)` = round(bias, 2), `RMSE (%)` = round(rmse, 2), `final coverage` = round(100 * final_coverage),
              `path coverage` = round(100 * path_coverage), `whole-period false pos.` = round(100 * whole_false_pos),
              `final half-width (%)` = round(median_hw_pct, 1), across(starts_with("power ")))
  widget <- DT::datatable(dt_tab, rownames = FALSE, extensions = "FixedColumns", filter = "top",
    caption = "Augmented comparison simulation (factor model, 40 donors). Values in %. Power: final period / anywhere (simultaneous band, placebo methods only).",
    options = list(scrollX = TRUE, scrollY = "70vh", scrollCollapse = TRUE, paging = FALSE, fixedColumns = list(leftColumns = 1)))
  tmp_html <- file.path(tempdir(), "augmented_comparison_table.html")
  htmlwidgets::saveWidget(widget, tmp_html, selfcontained = TRUE, title = "Augmented comparison")
  file.copy(tmp_html, file.path(out_dir, "augmented_comparison_table.html"), overwrite = TRUE)
  message("Finished in ", run_minutes, " min; ", nrow(errors), " method failures.")
}

