##############################################################################
#
# NOTE: kept in simulations/ as a dependency of pre_period_length_simulation.R, which uses
# this test version of dsc_symcal_boot() (its interface differs from the version in
# synthetic_control_functions_claude.R).
# demeaned_sc_symmetric_calboot_test.R
#
# TEST: demeaned synthetic control with a symmetric, mean-centred, donor-
# calibrated fixed-weight bootstrap, giving POINTWISE and SIMULTANEOUS effect
# intervals and prediction intervals over the whole post period, per period and
# cumulative, with ggplot figures. Not yet in synthetic_control_functions_claude.R.
#
# Method
#   1. Demeaned SC: intercept + non-negative donor weights summing to one, fitted
#      on the pre-treatment period (CVXR). Intercept and weights are then fixed.
#   2. Fixed-weight bootstrap: by-unit meboot replicates (reps, trim); in each
#      replicate, synthetic = intercept + replicate donor values x fixed weights;
#      replicate effect = observed - replicate synthetic (per period, and its
#      running total for the cumulative effect).
#   3. Symmetric half-width in every post period t: h_t = (1 - alpha) quantile of
#      |replicate effect_t - mean of replicate effects_t|.
#   4. Standardized gap: k_t = |gap_t| / h_t.
#   5. Calibration: every donor is fitted as a placebo (treated unit removed) with
#      steps 1-4. m = floor(alpha x (N + 1)).
#        pointwise:     k*_t = m-th largest donor k_t, separately for each t
#        simultaneous:  K*  = m-th largest donor max_t k_t (whole post period)
#   6. Effect intervals: gap_t +/- factor x h_t. Prediction intervals for the
#      untreated outcome: synthetic_t +/- factor x h_t (cumulative: running totals).
#      Pointwise p_t = (1 + #{donor k_t >= treated k_t}) / (N + 1); simultaneous
#      p = (1 + #{donor max_t k_t >= treated max_t k_t}) / (N + 1).
#   7. Improved simultaneous bands (keep a time profile instead of the bootstrap's):
#        loo:     r_pt = k_pt / k*_t(without donor p); lambda = m-th largest max_t r_pt;
#                 band = gap_t +/- lambda x k*_t x h_t (pointwise shape, uniformly widened)
#        median:  q_t = median donor k_pt; lambda = m-th largest max_t (k_pt / q_t);
#                 band = gap_t +/- lambda x q_t x h_t (robust median profile)
#      (Dividing by the all-donor k*_t would cap every donor ratio at 1 when m = 1.)
#   8. Smoother profiles (loess on the log scale over post periods, span = span [default 0.75], degree 1):
#        median_smooth: q~_t = smoothed median profile; lambda = m-th largest max_t (k_pt / q~_t);
#                       band = gap_t +/- lambda x q~_t x h_t
#        gapsd_smooth:  s~_t = smoothed MAD of the donors' raw gaps in period t (no bootstrap);
#                       lambda = m-th largest max_t |gap_pt| / s~_t; band = gap_t +/- lambda x s~_t
#
# Output (results_dsc_symcal_california/): dsc_symcal_california.rds,
# effects.csv, three PNG figures; summary written at the end of this script.
##############################################################################


# Setup
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
source("synthetic_control_functions/synthetic_control_functions_claude.R")   # meboot_multi()


# ---- Test function ------------------------------------------------------------
dsc_symcal_boot <- function(data, unit_col, time_col, outcome_col, treated_unit, treatment_time,
                            reps = 500, conf = 0.95, trim = list(trim = 0.10, xmin_rel = 0.9), seed = 12112025,
                            span = 0.75) {
  alpha <- round(1 - conf, 10)
  d <- data.frame(unit = as.character(data[[unit_col]]), time = data[[time_col]], y = data[[outcome_col]])
  stopifnot(!anyNA(d), length(unique(table(d$unit))) == 1, treated_unit %in% d$unit)
  times <- sort(unique(d$time)); post <- times >= treatment_time; n_post <- sum(post)
  units <- sort(unique(d$unit)); donors <- setdiff(units, treated_unit)
  Y <- sapply(units, function(u) { du <- d[d$unit == u, ]; du$y[order(du$time)] })

  # 1. Demeaned SC fit (CVXR): intercept + non-negative weights summing to one
  fit_dsc <- function(X, yv) {
    # Rescale before solving: large outcome values (e.g. GDP per capita) make the solver report a false
    # "infeasible". With weights summing to one the weights are scale-free; only the intercept is scaled back.
    sc <- mean(abs(yv[!post]))
    Xp <- X[!post, , drop = FALSE] / sc; yp <- yv[!post] / sc
    w <- CVXR::Variable(ncol(Xp)); b0 <- CVXR::Variable(1)
    prob <- CVXR::Problem(CVXR::Minimize(CVXR::sum_squares(yp - b0 - Xp %*% w)), list(w >= 0, sum(w) == 1))
    CVXR::psolve(prob)
    st <- CVXR::status(prob)
    if (!st %in% c("optimal", "optimal_inaccurate")) stop("CVXR status: ", st)
    wv <- setNames(as.numeric(CVXR::value(w)), colnames(X)); wv[abs(wv) < 1e-8] <- 0
    b <- as.numeric(CVXR::value(b0)) * sc                  # intercept back on the original scale
    list(b0 = b, w = wv, synth = b + as.numeric(X %*% wv))
  }
  # Bootstrap replicates for every unit in a dataset (list of time x unit matrices)
  replicates_for <- function(dd) {
    dd <- dd[order(dd$unit, dd$time), ]
    r <- withCallingHandlers(
      meboot_multi(dd, boot_cols = "y", id_cols = c("unit", "time"), reps = reps, trim = trim, group_col = "unit",
                   output = "wide", names_from = "unit", seed = seed, as_matrix = TRUE),
      warning = function(w) if (grepl("expand.sd", conditionMessage(w))) invokeRestart("muffleWarning"))
    lapply(r, function(m) m[order(m[, "time"]), setdiff(colnames(m), "time"), drop = FALSE])
  }
  # 2-4. Fixed-weight effects, symmetric mean-centred half-widths, standardized gaps (post periods)
  unit_stats <- function(fit, yv, reps_list) {
    gap <- yv - fit$synth
    E <- sapply(reps_list, function(m) yv[post] - as.numeric(fit$b0 + m[post, names(fit$w), drop = FALSE] %*% fit$w))
    Ec <- apply(E, 2, cumsum)
    hw <- function(M) apply(M, 1, function(v) unname(quantile(abs(v - mean(v)), 1 - alpha)))
    h <- hw(E); hc <- hw(Ec); g <- gap[post]; gc <- cumsum(g)
    list(gap = gap, synth = fit$synth, h = h, hc = hc, k = abs(g) / h, kc = abs(gc) / hc, w = fit$w, b0 = fit$b0)
  }

  # Treated unit (replicates of all units) and donors as placebos (treated unit removed)
  reps_all <- replicates_for(d)
  tr <- unit_stats(fit_dsc(Y[, donors, drop = FALSE], Y[, treated_unit]), Y[, treated_unit], reps_all)
  reps_ex <- replicates_for(d[d$unit != treated_unit, ])
  pl <- lapply(donors, function(p) tryCatch(unit_stats(fit_dsc(Y[, setdiff(donors, p), drop = FALSE], Y[, p]), Y[, p], reps_ex),
                                            error = function(e) NULL))
  ok <- !vapply(pl, is.null, logical(1))
  if (any(!ok)) warning(sum(!ok), " placebo fit(s) failed: ", paste(donors[!ok], collapse = ", "))
  pl <- pl[ok]; N <- length(pl); m <- floor(alpha * (N + 1))
  kth <- function(x) if (m == 0) Inf else sort(x, decreasing = TRUE)[m]
  kth_m <- function(x, mm) if (mm == 0) Inf else sort(x, decreasing = TRUE)[mm]

  # 5. Calibration factors
  Kp  <- sapply(pl, `[[`, "k");  Kpc <- sapply(pl, `[[`, "kc")          # post periods x donors
  Kp  <- matrix(Kp, nrow = n_post); Kpc <- matrix(Kpc, nrow = n_post)
  k_pw <- apply(Kp, 1, kth); k_pwc <- apply(Kpc, 1, kth)
  K_sim <- kth(apply(Kp, 2, max)); K_simc <- kth(apply(Kpc, 2, max))

  # 6. Intervals and p-values
  g <- tr$gap[post]; gc <- cumsum(g); s <- tr$synth[post]; sc <- cumsum(s)
  effects <- data.frame(time = times[post], observed = Y[post, treated_unit], synthetic = s, effect = g,
    halfwidth = tr$h, factor_pointwise = k_pw, factor_simultaneous = K_sim,
    eff_lo_pw = g - k_pw * tr$h, eff_hi_pw = g + k_pw * tr$h, eff_lo_sim = g - K_sim * tr$h, eff_hi_sim = g + K_sim * tr$h,
    pred_lo_pw = s - k_pw * tr$h, pred_hi_pw = s + k_pw * tr$h, pred_lo_sim = s - K_sim * tr$h, pred_hi_sim = s + K_sim * tr$h,
    p_pointwise = (1 + rowSums(Kp >= tr$k)) / (N + 1),
    cum_observed = cumsum(Y[post, treated_unit]), cum_synthetic = sc, cum_effect = gc, cum_halfwidth = tr$hc,
    cum_factor_pointwise = k_pwc, cum_factor_simultaneous = K_simc,
    cum_eff_lo_pw = gc - k_pwc * tr$hc, cum_eff_hi_pw = gc + k_pwc * tr$hc,
    cum_eff_lo_sim = gc - K_simc * tr$hc, cum_eff_hi_sim = gc + K_simc * tr$hc,
    cum_pred_lo_pw = sc - k_pwc * tr$hc, cum_pred_hi_pw = sc + k_pwc * tr$hc,
    cum_pred_lo_sim = sc - K_simc * tr$hc, cum_pred_hi_sim = sc + K_simc * tr$hc,
    cum_p_pointwise = (1 + rowSums(Kpc >= tr$kc)) / (N + 1))

  # 7. Improved simultaneous bands
  improved <- function(K, kT, k_pw, h, centre_eff, centre_pred) {
    m_loo <- floor(alpha * N)                                   # rank among the other N - 1 donors (N units)
    R_loo <- vapply(seq_len(ncol(K)), function(j) max(K[, j] / apply(K[, -j, drop = FALSE], 1, kth_m, mm = m_loo)), numeric(1))
    lam_loo <- kth(R_loo); f_loo <- lam_loo * k_pw
    q <- apply(K, 1, median); R_med <- apply(K / q, 2, max)
    lam_med <- kth(R_med); f_med <- lam_med * q
    list(f_loo = f_loo, f_med = f_med, lam_loo = lam_loo, lam_med = lam_med,
         p_loo = (1 + sum(R_loo >= max(kT / k_pw))) / (N + 1), p_med = (1 + sum(R_med >= max(kT / q))) / (N + 1),
         eff = cbind(lo_loo = centre_eff - f_loo * h, hi_loo = centre_eff + f_loo * h, lo_med = centre_eff - f_med * h, hi_med = centre_eff + f_med * h),
         pred = cbind(lo_loo = centre_pred - f_loo * h, hi_loo = centre_pred + f_loo * h, lo_med = centre_pred - f_med * h, hi_med = centre_pred + f_med * h))
  }
  imp  <- improved(Kp,  tr$k,  k_pw,  tr$h,  g,  s)
  impc <- improved(Kpc, tr$kc, k_pwc, tr$hc, gc, sc)

  # 8. Smoother profiles
  smooth_pos <- function(v) { tt <- seq_along(v); exp(as.numeric(predict(loess(log(v) ~ tt, span = span, degree = 1)))) }
  smoothers <- function(K, kT, h, Gd, gT, centre_eff, centre_pred) {
    qs <- smooth_pos(apply(K, 1, median)); R1 <- apply(K / qs, 2, max); lam1 <- kth(R1); f1 <- lam1 * qs * h
    ss <- smooth_pos(apply(Gd, 1, mad)); R2 <- apply(abs(Gd) / ss, 2, max); lam2 <- kth(R2); f2 <- lam2 * ss
    list(hw_ms = f1, hw_gs = f2, lam = c(median_smooth = lam1, gapsd_smooth = lam2),
         p = c(median_smooth = (1 + sum(R1 >= max(kT / qs))) / (N + 1), gapsd_smooth = (1 + sum(R2 >= max(abs(gT) / ss))) / (N + 1)),
         cols = cbind(lo_ms = centre_eff - f1, hi_ms = centre_eff + f1, lo_gs = centre_eff - f2, hi_gs = centre_eff + f2,
                      plo_ms = centre_pred - f1, phi_ms = centre_pred + f1, plo_gs = centre_pred - f2, phi_gs = centre_pred + f2))
  }
  Gd  <- matrix(sapply(pl, function(z) z$gap[post]), nrow = n_post)        # donor raw gaps, post periods x donors
  Gdc <- apply(Gd, 2, cumsum)
  sm  <- smoothers(Kp,  tr$k,  tr$h,  Gd,  g,  g,  s)
  smc <- smoothers(Kpc, tr$kc, tr$hc, Gdc, gc, gc, sc)
  effects <- cbind(effects,
    factor_sim_loo = imp$f_loo, factor_sim_median = imp$f_med,
    eff_lo_sim_loo = imp$eff[, "lo_loo"], eff_hi_sim_loo = imp$eff[, "hi_loo"], eff_lo_sim_med = imp$eff[, "lo_med"], eff_hi_sim_med = imp$eff[, "hi_med"],
    pred_lo_sim_loo = imp$pred[, "lo_loo"], pred_hi_sim_loo = imp$pred[, "hi_loo"], pred_lo_sim_med = imp$pred[, "lo_med"], pred_hi_sim_med = imp$pred[, "hi_med"],
    cum_factor_sim_loo = impc$f_loo, cum_factor_sim_median = impc$f_med,
    cum_eff_lo_sim_loo = impc$eff[, "lo_loo"], cum_eff_hi_sim_loo = impc$eff[, "hi_loo"], cum_eff_lo_sim_med = impc$eff[, "lo_med"], cum_eff_hi_sim_med = impc$eff[, "hi_med"],
    cum_pred_lo_sim_loo = impc$pred[, "lo_loo"], cum_pred_hi_sim_loo = impc$pred[, "hi_loo"], cum_pred_lo_sim_med = impc$pred[, "lo_med"], cum_pred_hi_sim_med = impc$pred[, "hi_med"])
  effects <- cbind(effects,
    eff_lo_sim_ms = sm$cols[, "lo_ms"], eff_hi_sim_ms = sm$cols[, "hi_ms"], eff_lo_sim_gs = sm$cols[, "lo_gs"], eff_hi_sim_gs = sm$cols[, "hi_gs"],
    pred_lo_sim_ms = sm$cols[, "plo_ms"], pred_hi_sim_ms = sm$cols[, "phi_ms"], pred_lo_sim_gs = sm$cols[, "plo_gs"], pred_hi_sim_gs = sm$cols[, "phi_gs"],
    cum_eff_lo_sim_ms = smc$cols[, "lo_ms"], cum_eff_hi_sim_ms = smc$cols[, "hi_ms"], cum_eff_lo_sim_gs = smc$cols[, "lo_gs"], cum_eff_hi_sim_gs = smc$cols[, "hi_gs"],
    cum_pred_lo_sim_ms = smc$cols[, "plo_ms"], cum_pred_hi_sim_ms = smc$cols[, "phi_ms"], cum_pred_lo_sim_gs = smc$cols[, "plo_gs"], cum_pred_hi_sim_gs = smc$cols[, "phi_gs"])
  p_improved <- rbind(per_period = c(loo = imp$p_loo, median = imp$p_med, sm$p), cumulative = c(loo = impc$p_loo, median = impc$p_med, smc$p))
  lambdas <- rbind(per_period = c(loo = imp$lam_loo, median = imp$lam_med, sm$lam), cumulative = c(loo = impc$lam_loo, median = impc$lam_med, smc$lam))
  list(effects = effects,
       pre = data.frame(time = times[!post], observed = Y[!post, treated_unit], synthetic = tr$synth[!post]),
       p_simultaneous = c(per_period = (1 + sum(apply(Kp, 2, max) >= max(tr$k))) / (N + 1),
                          cumulative = (1 + sum(apply(Kpc, 2, max) >= max(tr$kc))) / (N + 1)),
       p_simultaneous_improved = p_improved, lambdas_improved = lambdas,
       weights = sort(tr$w[tr$w > 0], decreasing = TRUE), intercept = tr$b0,
       pre_rmspe = sqrt(mean(tr$gap[!post]^2)),
       donor_factors = list(pointwise = Kp, cumulative = Kpc),
       donor_gaps = setNames(as.data.frame(sapply(pl, `[[`, "gap")), donors[ok]),        # all periods x donors
       donor_synth = setNames(as.data.frame(sapply(pl, `[[`, "synth")), donors[ok]),
       settings = list(treated_unit = treated_unit, treatment_time = treatment_time, conf = conf, reps = reps,
                       trim = trim, seed = seed, n_donors = N, rank_used = m, outcome = outcome_col, span = span))
}


# ---- Plotting function ----------------------------------------------------------
plot_dsc_symcal <- function(res, y_label = res$settings$outcome) {
  e <- res$effects; tt <- res$settings$treatment_time; lvl <- paste0(100 * res$settings$conf, "%")
  band_note <- paste0("Shaded: simultaneous ", lvl, " band (whole post period); dashed: pointwise ", lvl, " band")
  vline <- geom_vline(xintercept = tt - 0.5, linetype = "dotted", colour = "grey40")
  p_pred <- ggplot() +
    geom_ribbon(data = e, aes(time, ymin = pred_lo_sim, ymax = pred_hi_sim), fill = "steelblue", alpha = 0.25) +
    geom_line(data = e, aes(time, pred_lo_pw), linetype = "dashed", colour = "steelblue4") +
    geom_line(data = e, aes(time, pred_hi_pw), linetype = "dashed", colour = "steelblue4") +
    geom_line(data = rbind(res$pre[, c("time", "synthetic")], e[, c("time", "synthetic")]), aes(time, synthetic, colour = "Synthetic")) +
    geom_line(data = rbind(res$pre[, c("time", "observed")], e[, c("time", "observed")]), aes(time, observed, colour = "Observed")) +
    vline + scale_colour_manual(values = c(Observed = "black", Synthetic = "steelblue4"), name = NULL) +
    labs(title = paste0(res$settings$treated_unit, ": observed vs synthetic, with prediction intervals"),
         subtitle = band_note, x = NULL, y = y_label) + theme_minimal() + theme(legend.position = "bottom")
  p_eff <- ggplot(e, aes(time)) +
    geom_ribbon(aes(ymin = eff_lo_sim, ymax = eff_hi_sim), fill = "firebrick", alpha = 0.20) +
    geom_line(aes(y = eff_lo_pw), linetype = "dashed", colour = "firebrick4") +
    geom_line(aes(y = eff_hi_pw), linetype = "dashed", colour = "firebrick4") +
    geom_hline(yintercept = 0, colour = "grey40") + geom_line(aes(y = effect)) + geom_point(aes(y = effect), size = 1) +
    labs(title = paste0(res$settings$treated_unit, ": effect by period (observed - synthetic)"),
         subtitle = sprintf("%s. Simultaneous p = %.3f", band_note, res$p_simultaneous[["per_period"]]),
         x = NULL, y = paste("Effect,", y_label)) + theme_minimal()
  p_cum <- ggplot(e, aes(time)) +
    geom_ribbon(aes(ymin = cum_eff_lo_sim, ymax = cum_eff_hi_sim), fill = "darkgreen", alpha = 0.20) +
    geom_line(aes(y = cum_eff_lo_pw), linetype = "dashed", colour = "darkgreen") +
    geom_line(aes(y = cum_eff_hi_pw), linetype = "dashed", colour = "darkgreen") +
    geom_hline(yintercept = 0, colour = "grey40") + geom_line(aes(y = cum_effect)) + geom_point(aes(y = cum_effect), size = 1) +
    labs(title = paste0(res$settings$treated_unit, ": cumulative effect"),
         subtitle = sprintf("%s. Simultaneous p = %.3f", band_note, res$p_simultaneous[["cumulative"]]),
         x = NULL, y = paste("Cumulative effect,", y_label)) + theme_minimal()
  # Comparison of simultaneous bands
  long <- function(pre, lab) do.call(rbind, list(
    data.frame(time = e$time, band = "Pointwise", lo = e[[paste0(pre, "lo_pw")]], hi = e[[paste0(pre, "hi_pw")]]),
    data.frame(time = e$time, band = "Simultaneous: current", lo = e[[paste0(pre, "lo_sim")]], hi = e[[paste0(pre, "hi_sim")]]),
    data.frame(time = e$time, band = "Simultaneous: leave-one-out", lo = e[[paste0(pre, "lo_sim_loo")]], hi = e[[paste0(pre, "hi_sim_loo")]]),
    data.frame(time = e$time, band = "Simultaneous: median profile", lo = e[[paste0(pre, "lo_sim_med")]], hi = e[[paste0(pre, "hi_sim_med")]]),
    data.frame(time = e$time, band = "Simultaneous: smoothed median", lo = e[[paste0(pre, "lo_sim_ms")]], hi = e[[paste0(pre, "hi_sim_ms")]]),
    data.frame(time = e$time, band = "Simultaneous: donor-gap spread", lo = e[[paste0(pre, "lo_sim_gs")]], hi = e[[paste0(pre, "hi_sim_gs")]]))) |>
    transform(panel = lab)
  cmp <- rbind(long("eff_", "Effect by period"), long("cum_eff_", "Cumulative effect"))
  est <- rbind(data.frame(time = e$time, est = e$effect, panel = "Effect by period"),
               data.frame(time = e$time, est = e$cum_effect, panel = "Cumulative effect"))
  p_cmp <- ggplot(cmp, aes(time)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = band), alpha = 0.25) +
    geom_line(data = est, aes(time, est)) + geom_hline(yintercept = 0, colour = "grey40") +
    facet_grid(panel ~ band, scales = "free_y") + guides(fill = "none") +
    labs(title = paste0(res$settings$treated_unit, ": pointwise vs simultaneous ", lvl, " bands"), x = NULL, y = y_label) +
    theme_minimal() + theme(strip.text = element_text(size = 8))
  list(prediction = p_pred, effect = p_eff, cumulative = p_cum, simultaneous_comparison = p_cmp)
}


# ---- California test (Proposition 99) -------------------------------------------
script_file <- "simulations/demeaned_sc_symmetric_calboot_test.R"
out_dir <- file.path("simulations", "results_dsc_symcal_california")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
smoking <- as.data.frame(tidysynth::smoking)
t0 <- Sys.time()
ca <- dsc_symcal_boot(smoking, "state", "year", "cigsale", "California", treatment_time = 1989)
run_seconds <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
plots <- plot_dsc_symcal(ca, y_label = "cigarette packs per capita")
saveRDS(ca, file.path(out_dir, "dsc_symcal_california.rds"))
write.csv(ca$effects, file.path(out_dir, "effects.csv"), row.names = FALSE)
for (nm in names(plots)) ggsave(file.path(out_dir, paste0("california_", nm, ".png")), plots[[nm]], width = 8, height = 5, dpi = 150)

# Summary written at the end of this script
e <- ca$effects
summ <- data.frame(year = e$time, effect = round(e$effect, 1),
                   pointwise = sprintf("%.1f to %.1f", e$eff_lo_pw, e$eff_hi_pw), simultaneous = sprintf("%.1f to %.1f", e$eff_lo_sim, e$eff_hi_sim),
                   p_pw = round(e$p_pointwise, 3), cum_effect = round(e$cum_effect, 1),
                   cum_simultaneous = sprintf("%.1f to %.1f", e$cum_eff_lo_sim, e$cum_eff_hi_sim), cum_p_pw = round(e$cum_p_pointwise, 3))
marker <- "# ---- Results (written automatically when this script is run) ----"
script_lines <- readLines(script_file)
cut_at <- which(script_lines == marker)
if (length(cut_at) > 0) script_lines <- script_lines[seq_len(cut_at[1] - 1)]
writeLines(c(script_lines, marker,
  sprintf("# Run: %s | California, treatment 1989 | %d donors, rank %d | %d reps | %s s | pre-RMSPE %.2f | simultaneous p: per period %.3f, cumulative %.3f",
          format(Sys.time(), "%Y-%m-%d %H:%M"), ca$settings$n_donors, ca$settings$rank_used, ca$settings$reps, run_seconds,
          ca$pre_rmspe, ca$p_simultaneous[["per_period"]], ca$p_simultaneous[["cumulative"]]),
  sprintf("# Simultaneous factors: per period %.2f, cumulative %.2f", e$factor_simultaneous[1], e$cum_factor_simultaneous[1]),
  paste0("# ", knitr::kable(summ, format = "simple"))), script_file)

# ---- Results (written automatically when this script is run) ----
# Run: 2026-10-01 21:51 | California, treatment 1989 | 38 donors, rank 1 | 500 reps | 25 s | pre-RMSPE 0.96 | simultaneous p: per period 0.333, cumulative 0.462
# Simultaneous factors: per period 14.23, cumulative 10.95
#  year   effect  pointwise       simultaneous        p_pw   cum_effect  cum_simultaneous    cum_p_pw
# -----  -------  --------------  ----------------  ------  -----------  -----------------  ---------
#  1989     -5.8  -22.2 to 10.6   -130.3 to 118.8    0.513         -5.8  -101.6 to 90.1         0.513
#  1990     -4.3  -57.6 to 49.0   -125.6 to 117.0    0.615        -10.1  -193.1 to 173.0        0.590
#  1991     -7.3  -75.9 to 61.3   -130.1 to 115.5    0.615        -17.4  -282.3 to 247.5        0.590
#  1992     -6.1  -64.3 to 52.2   -120.3 to 108.2    0.718        -23.5  -365.6 to 318.7        0.564
#  1993     -8.9  -73.3 to 55.6   -97.0 to 79.3      0.513        -32.3  -428.5 to 363.9        0.564
#  1994    -10.8  -86.5 to 64.9   -90.8 to 69.2      0.462        -43.1  -490.6 to 404.3        0.462
#  1995    -13.0  -71.6 to 45.6   -88.2 to 62.2      0.436        -56.1  -538.2 to 425.9        0.487
#  1996    -12.5  -55.9 to 30.8   -55.9 to 30.8      0.385        -68.7  -588.1 to 450.8        0.513
#  1997    -12.9  -49.2 to 23.3   -62.7 to 36.9      0.359        -81.6  -606.2 to 443.1        0.436
#  1998    -15.7  -49.0 to 17.7   -65.7 to 34.4      0.256        -97.3  -646.3 to 451.7        0.462
#  1999    -18.7  -48.8 to 11.5   -72.4 to 35.1      0.256       -115.9  -696.5 to 464.6        0.462
#  2000    -17.4  -43.1 to 8.3    -67.7 to 33.0      0.154       -133.3  -736.2 to 469.6        0.436
