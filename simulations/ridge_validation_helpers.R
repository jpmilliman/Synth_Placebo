# ridge_validation_helpers.R
#
# Functions used by the ridge-augmentation simulations (ridge_edge_validation.R,
# ridge_small_sample_validation.R, ridge_covariate_validation.R, ridge_function_check.R).
# They are the prototype of dsc_placebo(augment = "ridge"), written for speed (quadprog weights instead
# of CVXR; one call returns every estimator x calibration). ridge_function_check.R confirms that
# dsc_placebo(augment = "ridge") gives the same estimates and intervals.
#
#   dsc_fit_qp()       demeaned synthetic control (intercept + simplex weights) by quadprog
#   smooth_profile()   LOESS profile on the log scale (as in dsc_placebo())
#   ridge_parts(), ridge_correction(), fit_target()
#                      the ridge correction: counterfactual = DSC + (x_T - sum_j w_j x_j)' beta_t,
#                      x_j unit j's demeaned pre-period path, beta_t a ridge regression across donors
#                      of their demeaned, centred period-t outcome on those paths; penalty
#                      c x mean squared centred path norm
#   choose_c()         penalty by cross-validation inside the pre-period (each donor predicts its last
#                      min(T1, T0/4) pre-periods from the earlier ones)
#   thresholds()       rank-rule thresholds, pointwise and simultaneous (as in dsc_placebo())
#   edge_core(), edge_bands()
#                      estimators dsc / ridge x calibrations space (in-space placebos), time (the
#                      treated unit's backdated windows) and both (the wider of the two)
#   gamma_twfe(), gamma_ife(), cov_core()
#                      covariate adjustment as in dsc_placebo(adjust_covariates = ...) followed by DSC
#                      or ridge; each placebo re-estimates gamma without the treated unit and without
#                      its own post period
#
# Written by Claude (Anthropic), at the direction of Jeffrey Milliman.
# Needs quadprog (and fect and withr for gamma_ife()).

# Demeaned synthetic control by quadprog: weights >= 0 summing to 1, intercept = mean pre-period gap
dsc_fit_qp <- function(y, X, pre) {
  X <- as.matrix(X); J <- ncol(X)
  if (J == 1) { w <- 1 } else {
    sc <- mean(abs(y[pre])); sc <- if (is.finite(sc) && sc > 0) sc else 1
    yp <- y[pre] / sc; Xp <- X[pre, , drop = FALSE] / sc
    yc <- yp - mean(yp); Xc <- sweep(Xp, 2, colMeans(Xp))
    D <- crossprod(Xc); D <- D + diag(1e-10 * max(mean(diag(D)), 1e-12), J)
    sol <- tryCatch(quadprog::solve.QP(D, crossprod(Xc, yc), cbind(rep(1, J), diag(J)), c(1, rep(0, J)), meq = 1)$solution,
                    error = function(e) rep(1 / J, J))
    w <- pmax(sol, 0); w <- w / sum(w)
  }
  b0 <- mean(y[pre] - X[pre, , drop = FALSE] %*% w)
  list(w = w, b0 = b0, synthetic = as.numeric(b0 + X %*% w))
}

# LOESS profile on the log scale
smooth_profile <- function(v, span = 0.75) { v <- pmax(v, 1e-8); tt <- seq_along(v)
  if (length(v) < 4) return(v)
  exp(as.numeric(stats::predict(stats::loess(log(v) ~ tt, span = span, degree = 1)))) }

# demeaned pre-period paths of the units in X (rows = units) and their outcomes in out_idx, demeaned
ridge_parts <- function(X, fit_idx, out_idx) {
  mu <- colMeans(X[fit_idx, , drop = FALSE])
  M <- t(sweep(X[fit_idx, , drop = FALSE], 2, mu)); Yo <- t(sweep(X[out_idx, , drop = FALSE], 2, mu))
  list(M = M, Mc = sweep(M, 2, colMeans(M)), Yoc = sweep(Yo, 2, colMeans(Yo)))
}
ridge_correction <- function(d, rp, c_ridge) {          # d: treated minus weighted donor pre-path
  lam <- c_ridge * mean(rowSums(rp$Mc^2))
  A <- t(rp$Mc) %*% solve(tcrossprod(rp$Mc) + diag(lam, nrow(rp$Mc)))
  as.numeric(crossprod(d, A %*% rp$Yoc))
}

# Fit a target series y on donor matrix X using periods fit_idx; gaps on out_idx; c_ridge NULL = plain DSC
fit_target <- function(y, X, fit_idx, out_idx, c_ridge = NULL) {
  f <- dsc_fit_qp(y, X, seq_along(y) %in% fit_idx)            # weights and intercept from fit_idx only
  syn <- f$synthetic; gap_out <- (y - syn)[out_idx]
  if (!is.null(c_ridge)) {
    rp <- ridge_parts(X, fit_idx, out_idx)
    d <- (y[fit_idx] - mean(y[fit_idx])) - as.numeric(crossprod(rp$M, f$w))
    gap_out <- gap_out - ridge_correction(d, rp, c_ridge)
  }
  list(gap = gap_out, s = sqrt(mean((y - syn)[fit_idx]^2)), w = f$w)
}

# Cross-validated penalty: each donor predicted for the last k pre-periods from the earlier ones
choose_c <- function(Y, donors, pre_idx, T1, grid = 10^seq(-3, 2, by = 0.5)) {
  k <- min(T1, floor(length(pre_idx) / 4)); F_idx <- head(pre_idx, -k); V_idx <- tail(pre_idx, k)
  err <- sapply(donors, function(j) {
    others <- setdiff(donors, j); y <- Y[, j]; X <- Y[, others, drop = FALSE]
    f <- dsc_fit_qp(y, X, seq_along(y) %in% F_idx)
    g <- (y - f$synthetic)[V_idx]
    rp <- ridge_parts(X, F_idx, V_idx); d <- (y[F_idx] - mean(y[F_idx])) - as.numeric(crossprod(rp$M, f$w))
    vapply(grid, function(cc) mean((g - ridge_correction(d, rp, cc))^2), numeric(1))
  })
  grid[which.min(rowMeans(err))]
}

# Thresholds (on the standardized scale) from a reference set K (periods x references) by the rank rule
thresholds <- function(K, alpha, span = 0.75) {
  n <- ncol(K); m <- floor(round(alpha, 10) * (n + 1) + 1e-9)
  kth <- function(v) if (m >= 1) sort(v, decreasing = TRUE)[m] else Inf
  q <- smooth_profile(apply(K, 1, stats::median), span); R <- apply(K / q, 2, max)
  list(pw = apply(K, 1, kth), sim = kth(R) * q)
}

# Everything for one panel: Y (T x units, column I treated), post logical. Returns the treated unit's
# gaps (per period and cumulative), its RMSPE, and the standardized reference sets for each calibration
edge_core <- function(Y, I, post, estimators = c("dsc", "ridge"), min_fit = 10, rmspe_floor = 0.25) {
  pre_idx <- which(!post); post_idx <- which(post); T1 <- length(post_idx); T0 <- length(pre_idx)
  donors <- setdiff(seq_len(ncol(Y)), I)
  c_r <- if ("ridge" %in% estimators) choose_c(Y, donors, pre_idx, T1) else NULL
  onsets <- if (T0 - T1 >= min_fit) (min_fit + 1):(T0 - T1 + 1) else integer(0)
  out <- lapply(estimators, function(est) {
    cr <- if (est == "ridge") c_r else NULL
    tr <- fit_target(Y[, I], Y[, donors, drop = FALSE], pre_idx, post_idx, cr)
    pl <- lapply(donors, function(j) fit_target(Y[, j], Y[, setdiff(donors, j), drop = FALSE], pre_idx, post_idx, cr))
    sP <- vapply(pl, `[[`, numeric(1), "s"); flP <- rmspe_floor * stats::median(c(tr$s, sP))
    sT <- max(tr$s, flP); sP <- pmax(sP, flP)
    Gp <- sapply(pl, `[[`, "gap")
    space <- list(K = abs(Gp) / rep(sP, each = T1), Kc = abs(apply(Gp, 2, cumsum)) / rep(sP, each = T1))
    time <- NULL
    if (length(onsets)) {
      win <- lapply(onsets, function(o) fit_target(Y[, I], Y[, donors, drop = FALSE], seq_len(o - 1), o:(o + T1 - 1), cr))
      sW <- vapply(win, `[[`, numeric(1), "s"); sW <- pmax(sW, rmspe_floor * stats::median(sW))
      Gw <- sapply(win, `[[`, "gap")
      time <- list(K = abs(Gw) / rep(sW, each = T1), Kc = abs(apply(Gw, 2, cumsum)) / rep(sW, each = T1))
    }
    list(gap = tr$gap, s = sT, space = space, time = time)
  })
  names(out) <- estimators
  list(est = out, c_ridge = c_r, n_windows = length(onsets))
}

# Half-widths (per period and cumulative; pointwise and simultaneous) for each estimator x calibration
edge_bands <- function(core, conf) {
  a <- 1 - conf
  res <- list()
  for (est in names(core$est)) {
    e <- core$est[[est]]
    th <- list(space = list(per = thresholds(e$space$K, a), cum = thresholds(e$space$Kc, a)))
    if (!is.null(e$time)) {
      th$time <- list(per = thresholds(e$time$K, a), cum = thresholds(e$time$Kc, a))
      th$both <- lapply(c(per = "per", cum = "cum"), function(m) list(pw = pmax(th$space[[m]]$pw, th$time[[m]]$pw),
                                                                       sim = pmax(th$space[[m]]$sim, th$time[[m]]$sim)))
    }
    for (cal in names(th)) res[[paste(est, cal, sep = "_")]] <- list(
      gap = e$gap, gapc = cumsum(e$gap),
      hw_pw = th[[cal]]$per$pw * e$s, hw_sim = th[[cal]]$per$sim * e$s,
      hwc_pw = th[[cal]]$cum$pw * e$s, hwc_sim = th[[cal]]$cum$sim * e$s)
  }
  res
}

# ---- Covariate adjustment + DSC / ridge ----------------------------------------------------------
# long data of the kept units; cells with D = 1 (drop_post unit's post period) are excluded from estimation
cov_long <- function(Y, Z, keep, drop_post, post) {
  tibble(t = rep(seq_len(nrow(Y)), length(keep)), unit = rep(keep, each = nrow(Y)),
         y = as.vector(Y[, keep]), z = as.vector(Z[, keep]),
         D = as.integer(unit %in% drop_post & rep(post, length(keep))))
}
gamma_twfe <- function(Y, Z, keep, drop_post, post) {
  lg <- cov_long(Y, Z, keep, drop_post, post) |> filter(D == 0)
  unname(coef(lm(y ~ z + factor(unit) + factor(t), data = lg))["z"])
}
gamma_ife <- function(Y, Z, keep, drop_post, post, r = 0:5, cv = TRUE, ife_seed = 1) {
  lg <- as.data.frame(cov_long(Y, Z, keep, drop_post, post))
  if (length(r) == 1) cv <- FALSE
  fit <- NULL
  utils::capture.output(fit <- withr::with_seed(ife_seed, suppressMessages(suppressWarnings(fect::fect(
    y ~ D + z, data = lg, index = c("unit", "t"), method = "ife", r = if (cv) range(r) else r, CV = cv,
    force = "two-way", se = FALSE, parallel = FALSE)))))
  list(gamma = as.numeric(fit$beta["z", 1]), r = if (cv) as.integer(fit$r.cv) else as.integer(r))
}

# Y, Z: T x units (outcome, covariate); column I treated; post logical
# adjust: "none", "twfe" or "ife"; estimators: "dsc" and/or "ridge"
# ife_r: candidate numbers of factors (one value skips the cross-validation); ife_seed: fect's seed
cov_core <- function(Y, Z, I, post, adjust = c("none", "twfe", "ife"), estimators = c("dsc", "ridge"), rmspe_floor = 0.25,
                     ife_r = 0:5, ife_seed = 1) {
  adjust <- match.arg(adjust)
  pre_idx <- which(!post); post_idx <- which(post); T1 <- length(post_idx)
  units <- seq_len(ncol(Y)); donors <- setdiff(units, I)
  r_use <- NA_integer_
  g_main <- switch(adjust, none = 0, twfe = gamma_twfe(Y, Z, units, I, post),
                   ife = { r_use <- if (length(ife_r) == 1) ife_r else gamma_ife(Y, Z, units, I, post, r = ife_r, ife_seed = ife_seed)$r
                           gamma_ife(Y, Z, units, I, post, r = r_use, ife_seed = ife_seed)$gamma })
  g_pl <- vapply(donors, function(p) switch(adjust, none = 0, twfe = gamma_twfe(Y, Z, donors, p, post),
                                            ife = gamma_ife(Y, Z, donors, p, post, r = r_use, ife_seed = ife_seed)$gamma), numeric(1))
  Ya <- Y - g_main * Z
  c_r <- if ("ridge" %in% estimators) choose_c(Ya, donors, pre_idx, T1) else NULL
  out <- lapply(estimators, function(est) {
    cr <- if (est == "ridge") c_r else NULL
    tr <- fit_target(Ya[, I], Ya[, donors, drop = FALSE], pre_idx, post_idx, cr)
    pl <- lapply(seq_along(donors), function(k) { p <- donors[k]; Yp <- Y - g_pl[k] * Z
      fit_target(Yp[, p], Yp[, setdiff(donors, p), drop = FALSE], pre_idx, post_idx, cr) })
    sP <- vapply(pl, `[[`, numeric(1), "s"); flP <- rmspe_floor * stats::median(c(tr$s, sP))
    sT <- max(tr$s, flP); sP <- pmax(sP, flP); Gp <- sapply(pl, `[[`, "gap")
    list(gap = tr$gap, s = sT, time = NULL,
         space = list(K = abs(Gp) / rep(sP, each = T1), Kc = abs(apply(Gp, 2, cumsum)) / rep(sP, each = T1)))
  })
  names(out) <- estimators
  list(est = out, c_ridge = c_r, gamma = g_main, gamma_placebo_sd = stats::sd(g_pl), r = r_use, n_windows = 0)
}
