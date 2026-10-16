##############################################################################
# diagnostics_validation.R
#
# Validates dsc_placebo()'s $diagnostics on the 1,800 simulated panels of
# augmented_comparison_simulation.R (same seeds, so the panels are identical).
# For each panel: the diagnostics, the treated unit's TRUE position (whether its
# loadings lie outside the donors' convex hull, from the generator), and later a
# link to whether each method's intervals covered in the comparison run.
#
# Output (results_diagnostics_validation/): tasks/*.rds, diagnostics.rds
##############################################################################
library(tidyverse); library(future); library(furrr)
dry_run <- TRUE; source("simulations/augmented_comparison_simulation.R")    # grid, seed, T1, simulate_panel()
dv_dir <- "simulations/results_diagnostics_validation"; dir.create(file.path(dv_dir, "tasks"), recursive = TRUE, showWarnings = FALSE)

# True treated loadings, regenerated exactly as simulate_panel() draws them (same seed, same order of draws)
true_outside <- function(treated, factors, n_donors, T0, rep_seed) {
  set.seed(rep_seed); TT <- T0 + T1; N <- n_donors + 1; K <- 2
  invisible(sapply(seq_len(K), function(k) { u <- if (factors == "stationary") rnorm(TT, 0, 0.05) else rnorm(TT, 0.005, 0.03); u }))
  lam <- matrix(rnorm(N * K), N, K)
  if (treated == "edge")    lam[1, ] <- apply(lam[-1, , drop = FALSE], 2, max)
  if (treated == "outside") lam[1, ] <- apply(lam[-1, , drop = FALSE], 2, max) + 1
  1 %in% grDevices::chull(lam)
}

diag_task <- function(cell_id, chunk) {
  f <- file.path(dv_dir, "tasks", sprintf("cell%02d_chunk%02d.rds", cell_id, chunk))
  if (file.exists(f)) return(readRDS(f))
  if (!exists("dsc_placebo", envir = globalenv())) source(functions_src)
  cell <- grid |> filter(cell_id == !!cell_id) |> as.list()
  reps <- ((chunk - 1) * reps_per_task + 1):min(chunk * reps_per_task, cell$reps)
  res <- map_dfr(reps, function(r) {
    s <- seed + 100000 * cell$cell_id + r
    d <- simulate_panel(cell$treated, cell$factors, cell$n_donors, cell$T0, rep_seed = s)
    out <- tryCatch({
      dg <- dsc_placebo(as.data.frame(d), "unit", "period", "y", "u01", treatment_period = c(period = cell$T0 + 1),
                        scale = "rmspe", warn = FALSE)$diagnostics
      tibble(fit_share_worse = dg$fit_rank$share_placebos_worse, fit_flag = dg$fit_rank$flag,
             n_factors = dg$hull$n_factors, hull_share_beyond = dg$hull$share_donors_beyond %||% NA_real_,
             hull_flag = dg$hull$flag, drift_pct = dg$extrapolation$max_post_drift_pct %||% NA_real_,
             max_ar1 = if (is.null(dg$persistence$factor_ar1)) NA_real_ else max(dg$persistence$factor_ar1),
             trending = dg$persistence$trending, error = NA_character_)
    }, error = function(e) tibble(error = conditionMessage(e)))
    out |> mutate(cell_id = cell$cell_id, rep = r,
                  true_outside_hull = true_outside(cell$treated, cell$factors, cell$n_donors, cell$T0, s))
  })
  saveRDS(res, f); res
}

plan(multisession, workers = 10)
t0 <- Sys.time()
res <- future_pmap(tasks, possibly(diag_task, otherwise = NULL),
                   .options = furrr_options(seed = TRUE, packages = c("tidyverse", "CVXR", "fect")))
plan(sequential)
diag_all <- bind_rows(res) |> left_join(grid, by = "cell_id")
saveRDS(list(diagnostics = diag_all, run_minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)),
        file.path(dv_dir, "diagnostics.rds"))
message("done")

