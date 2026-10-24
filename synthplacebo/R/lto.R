# dsc_placebo_lto.R
#
# Leave-two-out (LTO) placebo inference for the demeaned synthetic control:
# an experimental companion to dsc_placebo() for panels with few donors.
#
# Written by Claude (Anthropic), at the direction of Jeffrey Milliman. Experimental:
# not part of the main method, and validated only by the simulations in this folder.
#
# Based on Lei, L. and Sudijono, T. (2025), "Inference for synthetic controls via refined
# placebo tests", arXiv:2401.07152.
#
# The estimator and the statistic are those of dsc_placebo() on its default scale:
#   - demeaned synthetic control: an intercept plus non-negative donor weights that sum to
#     one, fitted to the pre-treatment outcomes;
#   - each unit's post-treatment gap divided by its own pre-treatment RMSPE (with the same
#     RMSPE floor), at the pre-specified test period (pointwise) or cumulated through it.
# Only the calibration changes. dsc_placebo() ranks the treated unit's statistic among N
# placebo statistics, so its p-value lives on {1/(N+1), ..., 1} and a 95% test needs at least
# 19 donors. The LTO test plays a "tournament": for every pair of donors (i, j), the treated
# unit I, i and j are each fitted with every OTHER unit as donors, and I "wins" the match if
# its statistic is larger than both of theirs. The naive LTO p-value is the share of the
# N_d (N_d - 1) / 2 matches that I does not win (N_d = number of donors).
#
# Validity (Lei and Sudijono, Theorems 2.1-2.3): under the sharp null and a uniformly
# assigned treated unit among the N = N_d + 1 units, P(p_naive <= alpha) <= floor(N f(N, alpha)) / N,
# which equals the approximate-placebo bound (floor(N alpha) + 1) / N for the (N, alpha) of
# most applications. That bound can exceed alpha: the guarantee is weaker than the exact
# placebo test's, and the method relies on the actual Type-I error being lower in practice,
# which the simulations in this folder check. The "powered" version rejects when
# p_naive <= alpha + c(N, alpha) - delta, with the same bound.
#
# Intervals: the effect theta in the test period (or the cumulative effect) is in the 1 - alpha
# set when the test of "effect = theta" does not reject. Subtracting theta from the treated
# unit's post-treatment outcomes changes only its own statistic (all fits use pre-treatment
# outcomes only), so in each match I does not win exactly when theta lies in
#   [g_ij - s_ij M_ij, g_ij + s_ij M_ij],
# where g_ij and s_ij are I's gap and pre-treatment RMSPE with i and j removed from its pool
# and M_ij = max(k_i, k_j) is the larger placebo statistic in the match. The confidence set is
# the set of theta covered by more than alpha x (number of matches) of these intervals
# (alpha + c - delta for the powered version), found exactly by sweeping the endpoints.
# The prediction interval for the treated unit's counterfactual is observed - (confidence set).

# ---- the demeaned synthetic control, fitted by quadratic programming ---------------------
# Same problem as in dsc_placebo() (which uses CVXR): minimise sum_pre (y - b0 - X w)^2 over
# b0 and w >= 0, sum(w) = 1. b0 is profiled out by demeaning over the pre-period.
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

# ---- ridge augmentation (as dsc_placebo(augment = "ridge")) ----------------------------------
# For unit y fitted on pool X with DSC weights w: post-period correction d' A Yoc, where d is the pre-
# period imbalance the weights leave (demeaned y path minus the weighted demeaned pool paths), A the ridge
# map from the pool's centred pre-period paths, and Yoc the pool's demeaned, centred post outcomes.
# Returns the correction over all periods (0 before treatment) and, because it is linear in the pool's
# post outcomes, the coefficient on each pool member's post outcome (for exact effect shifts).
lto_ridge <- function(y, X, w, c_ridge, pre) {
  post <- !pre; mu <- colMeans(X[pre, , drop = FALSE])
  M <- t(sweep(X[pre, , drop = FALSE], 2, mu)); Mc <- sweep(M, 2, colMeans(M))
  Yo <- t(sweep(X[post, , drop = FALSE], 2, mu)); Yoc <- sweep(Yo, 2, colMeans(Yo))
  dv <- (y[pre] - mean(y[pre])) - as.numeric(crossprod(M, w))
  lam <- c_ridge * mean(rowSums(Mc^2))
  b <- as.numeric(crossprod(dv, t(Mc) %*% solve(tcrossprod(Mc) + diag(lam, nrow(Mc)))))   # weight on each pool member
  corr <- numeric(length(y)); corr[post] <- as.numeric(crossprod(b, Yoc))
  list(correction = corr, coef = b - mean(b))                      # coef: d corr_t / d (member's post outcome_t)
}
# Penalty by cross-validation inside the pre-period. Each unit (fitted on the others in `units`; all units
# by default, so the choice does not depend on which unit is treated) predicts k = min(T1, T0 / 4)
# held-out pre-periods from the periods before them. Options:
#   windows  1 (default): the last k pre-periods only. More: rolling origins, i.e. also the k periods
#            before those, and so on, each predicted from the periods before it, while at least
#            max(5, T0 / 2) periods remain to fit on; errors are averaged over the windows.
#   rule     "min" (default): the grid value with the lowest mean error. "1se": the largest (most
#            shrinkage) value whose mean error is within one standard error (over units) of the
#            minimum, as in augsynth's min_1se.
# Returns the chosen value; with details = TRUE also the error table and the number of windows used.
lto_choose_ridge <- function(Y, pre, grid = 10^seq(-3, 2, by = 0.5), rule = c("min", "1se"), windows = 1,
                             units = seq_len(ncol(Y)), details = FALSE) {
  rule <- match.arg(rule)
  pre_idx <- which(pre); T0 <- length(pre_idx); T1 <- sum(!pre); k <- min(T1, floor(T0 / 4))
  if (T0 < 8) stop("ridge augmentation with ridge_c = \"cv\" needs at least 8 pre-treatment periods; use a numeric ridge_c (default 0.1).")
  min_fit <- max(5, ceiling(T0 / 2))
  ends <- T0 - k * (seq_len(windows) - 1); ends <- ends[ends - k >= min(min_fit, T0 - k)]   # last index of each validation block
  if (!length(ends)) ends <- T0
  err <- sapply(units, function(u) {
    others <- setdiff(units, u)
    rowMeans(sapply(ends, function(e) {
      Fr <- seq_len(nrow(Y)) %in% pre_idx[seq_len(e - k)]; Vr <- seq_len(nrow(Y)) %in% pre_idx[(e - k + 1):e]; rows <- Fr | Vr
      X <- Y[rows, others, drop = FALSE]; y <- Y[rows, u]; fr <- Fr[rows]
      f <- dsc_fit_qp(y, X, fr); g <- (y - f$synthetic)[!fr]
      vapply(grid, function(cc) mean((g - lto_ridge(y, X, f$w, cc, fr)$correction[!fr])^2), numeric(1))
    }))
  })
  err <- matrix(err, nrow = length(grid))
  m <- rowMeans(err); se <- apply(err, 1, stats::sd) / sqrt(ncol(err)); i <- which.min(m)
  pick <- if (rule == "min") grid[i] else max(grid[m <= m[i] + se[i]])
  if (details) list(c = pick, cv = data.frame(c = grid, mse = m, se = se), windows = length(ends)) else pick
}

# ---- covariate adjustment and IFE augmentation (as in dsc_placebo()) --------------------------------
# One model per unit u, estimated on every cell except u's own post period and the treated unit's post
# period (so the treated unit's post-treatment outcomes never enter any estimate, and lto_shift() stays
# exact). Every fit of u (in any match, the leave-one-out fit and the standard placebo fit) uses the
# outcomes adjusted by u's model:
#   adjust  : y - gamma_u' z, with gamma_u from two-way fixed effects ("twfe") or interactive fixed
#             effects via fect ("ife"; the number of factors chosen by cross-validation once, on the
#             treated unit's model, and reused)
#   augment : y - M_u, with M_u fect's IFE prediction (unit and time effects, factors and covariates)
# dsc_placebo() estimates the placebos' models without the treated unit at all; here the treated unit's
# pre-period stays in every model (its post period is always left out). Per-match models (without the
# other two members) would be exact but cost N x (N-1)(N-2)/2 fits.
lto_with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed); expr
}
lto_long <- function(Y, Zs, excl) {
  d <- data.frame(unit = rep(colnames(Y), each = nrow(Y)), t = rep(seq_len(nrow(Y)), ncol(Y)), y = as.vector(Y), D = as.integer(as.vector(excl)))
  for (v in names(Zs)) d[[v]] <- as.vector(Zs[[v]])
  d
}
lto_gamma_twfe <- function(Y, Zs, excl) {
  d <- lto_long(Y, Zs, excl); d <- d[d$D == 0, ]
  X <- stats::model.matrix(stats::as.formula(paste("~ factor(unit) + factor(t) +", paste(names(Zs), collapse = " + "))), d)
  b <- qr.coef(qr(X), d$y)[names(Zs)]
  if (anyNA(b)) stop("adjust_covariates: gamma not identified (no variation beyond unit and period effects).")
  b
}
lto_fect <- function(Y, Zs, excl, r, cv, seed) {
  d <- lto_long(Y, Zs, excl); f <- stats::as.formula(paste("y ~", paste(c("D", names(Zs)), collapse = " + ")))
  if (length(unique(r)) == 1) { cv <- FALSE; r <- unique(r) }
  fit <- NULL
  utils::capture.output(fit <- lto_with_seed(seed, suppressMessages(suppressWarnings(fect::fect(f, data = d, index = c("unit", "t"), method = "ife",
    r = if (cv) range(r) else r, CV = cv, force = "two-way", se = FALSE, parallel = FALSE)))))
  fit
}
# lapply with an optional PSOCK cluster (n_cores > 1) and an optional progress bar (sequential only).
# FUN's environment is sent to the workers; packages called with :: load there on first use.
lto_par_lapply <- function(X, FUN, n_cores = 1, progress = FALSE) {
  n <- length(X); if (!n) return(list())
  if (n_cores > 1 && n > 1) {
    cl <- parallel::makeCluster(min(n_cores, n)); on.exit(parallel::stopCluster(cl), add = TRUE)
    .sp_cluster_setup(cl)                               # the workers start empty: load the package there
    t0 <- Sys.time(); out <- parallel::parLapplyLB(cl, X, FUN)
    if (progress) message("  done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
    return(out)
  }
  if (progress) { pb <- utils::txtProgressBar(min = 0, max = n, style = 3); on.exit(close(pb), add = TRUE) }
  out <- vector("list", n)
  for (i in seq_len(n)) { out[[i]] <- FUN(X[[i]]); if (progress) utils::setTxtProgressBar(pb, i) }
  out
}
# Returns list(Ylist = per-unit outcome matrices, gamma = units x covariates (adjust) or NULL, r)
lto_unit_models <- function(Y, Zs, treated, post, mode = c("adjust", "augment"), model = c("ife", "twfe"), ife_r = 0:5, ife_seed = 12112025,
                            n_cores = 1, verbose = FALSE) {
  mode <- match.arg(mode); model <- match.arg(model); N <- ncol(Y)
  excl_for <- function(u) { e <- matrix(FALSE, nrow(Y), N); e[post, c(u, treated)] <- TRUE; e }
  use_ife <- mode == "augment" || model == "ife"
  r_use <- NA_integer_
  if (use_ife) {
    if (!requireNamespace("fect", quietly = TRUE)) stop("interactive fixed effects need the fect package.")
    if (verbose && length(unique(ife_r)) > 1) message("dsc_placebo_lto: choosing the number of factors by cross-validation (fect; the slowest single step)")
    f0 <- lto_fect(Y, Zs, excl_for(treated), ife_r, TRUE, ife_seed)
    r_use <- if (length(unique(ife_r)) == 1) unique(ife_r) else as.integer(f0$r.cv)
  }
  gam <- if (mode == "adjust" && length(Zs)) matrix(NA_real_, N, length(Zs), dimnames = list(colnames(Y), names(Zs))) else NULL
  if (verbose && use_ife) message("dsc_placebo_lto: fitting ", N, " unit models (one fect fit each", if (n_cores > 1) paste0(", on ", min(n_cores, N), " cores") else "", ")")
  fits <- lto_par_lapply(seq_len(N), function(u) {
    ex <- excl_for(u)
    if (mode == "adjust") {
      g <- if (model == "twfe") lto_gamma_twfe(Y, Zs, ex) else {
        fit <- lto_fect(Y, Zs, ex, r_use, FALSE, ife_seed); stats::setNames(as.numeric(fit$beta[names(Zs), 1]), names(Zs)) }
      Ya <- Y; for (v in names(Zs)) Ya <- Ya - g[[v]] * Zs[[v]]; list(Y = Ya, g = g)
    } else {
      fit <- lto_fect(Y, Zs, ex, r_use, FALSE, ife_seed)
      M <- if (is.matrix(fit$Y.ct.full)) fit$Y.ct.full else fit$Y.ct; colnames(M) <- as.character(fit$id)
      list(Y = Y - M[, colnames(Y), drop = FALSE], g = NULL)
    }
  }, n_cores = if (use_ife) n_cores else 1, progress = verbose && use_ife)
  Ylist <- lapply(fits, `[[`, "Y")
  if (!is.null(gam)) for (u in seq_len(N)) gam[u, ] <- fits[[u]]$g
  list(Ylist = Ylist, gamma = gam, r = r_use)
}

# ---- diagnostics: is the treated unit beyond what the donors can reproduce? --------------------
# The same checks as dsc_placebo(diagnostics = TRUE), computed from this tournament's own fits, so the LTO
# functions stand alone (no dsc_placebo(), no CVXR):
#   fit_rank      the treated unit's pre-treatment RMSPE relative to its level, against the standard
#                 placebos' (flag: fewer than 10% of placebos fit worse)
#   hull          distance of the treated unit's IFE factor loadings (fect, every cell except the treated
#                 unit's post period) from the donors' convex hull, against each donor's distance from the
#                 other donors' hull (flag: fewer than 10% of donors lie further out); needs fect
#   extrapolation post-treatment drift implied by the loading mismatch (% of the counterfactual)
#   persistence   AR(1) of each estimated factor (trending: largest above 0.95)
#   residual_autocorrelation  persistence of the untreated gaps: pooled AR(1) of the placebos' standardized
#                 post-treatment gaps, with a bootstrap interval and a guide sentence (information only)
# core: lto_core() output (standard placebo gaps in core$standard$gaps); w: the treated unit's weights on the
# donors; Zs: covariate matrices (enter the fect model as in dsc_placebo()); r_known: number of factors when
# already chosen (e.g. by lto_unit_models()), else ife_r is cross-validated.
lto_diagnostics <- function(Y, tr, post, core, w, Zs = list(), ife_r = 0:5, ife_seed = 12112025, r_known = NA,
                            season_length = NA_integer_) {
  pre <- !post; dn <- core$donors; gT <- core$gap_treated; GP <- core$standard$gaps; n_post <- sum(post)
  rel_rmspe <- function(gap, yv) sqrt(mean(gap[pre]^2)) / mean(abs(yv[pre]))
  fit_T <- rel_rmspe(gT, Y[, tr]); fit_P <- vapply(seq_along(dn), function(j) rel_rmspe(GP[, j], Y[, dn[j]]), numeric(1))
  share_worse <- mean(fit_P > fit_T)
  diag <- list(fit_rank = list(treated_rel_rmspe = fit_T, placebo_rel_rmspe_median = stats::median(fit_P),
                               share_placebos_worse = share_worse, flag = share_worse < 0.10))
  # residual autocorrelation (as in dsc_placebo())
  pooled_ar <- function(G, k = 1) {
    if (nrow(G) <= k + 1) return(NA_real_)
    a <- G[-seq_len(k), , drop = FALSE]; b <- G[seq_len(nrow(G) - k), , drop = FALSE]
    if (sum(b^2) == 0) NA_real_ else sum(a * b) / sum(b^2)
  }
  acf_in <- function(g) { g <- g - mean(g); if (sum(g^2) == 0) NA_real_ else sum(g[-1] * g[-length(g)]) / sum(g^2) }
  T0 <- sum(pre); rmspe_P <- sqrt(colMeans(GP[pre, , drop = FALSE]^2))
  Zac <- sweep(GP[post, , drop = FALSE], 2, pmax(rmspe_P, 0.25 * stats::median(rmspe_P)), "/")
  rho_hat <- pooled_ar(Zac)
  rho_int <- if (is.na(rho_hat)) c(NA_real_, NA_real_) else
    stats::quantile(lto_with_seed(ife_seed, replicate(200, pooled_ar(Zac[, sample(ncol(Zac), replace = TRUE), drop = FALSE]))), c(0.05, 0.95), na.rm = TRUE)
  rho_in <- stats::median(c(acf_in(gT[pre]), apply(GP[pre, , drop = FALSE], 2, acf_in)), na.rm = TRUE)
  s_len <- season_length
  rho_seas <- if (!is.na(s_len) && s_len %in% 2:12 && n_post >= s_len + 3) pooled_ar(Zac, s_len) else NA_real_
  rc <- if (is.na(rho_hat)) NA_real_ else min(max(rho_hat, 0), 0.99)
  half_life <- if (is.na(rc)) NA_real_ else if (rho_hat >= 1) Inf else if (rc == 0) 0 else log(0.5) / log(rc)
  backdate <- if (is.na(rc) || rc >= 0.9) NA_real_ else if (rc == 0) 1 else max(1, ceiling(log(0.25) / log(rc)))
  seas_excess <- if (is.na(rho_seas) || is.na(rc)) NA_real_ else rho_seas - rc^s_len
  eff_pre <- if (is.na(rc)) NA_real_ else T0 * (1 - rc) / (1 + rc)
  level <- if (is.na(rho_hat)) NA_character_ else as.character(cut(rho_hat, c(-Inf, 0.3, 0.7, 0.9, Inf), c("low", "moderate", "high", "very high"), right = FALSE))
  guide <- if (is.na(rho_hat)) "not estimated: needs at least 3 post-treatment periods" else paste0(switch(level,
    "low" = sprintf("Low persistence (rho = %.2f): untreated gaps behave close to independent noise. The pre-treatment fit is informative (about %.0f effective pre-periods), a shock just before treatment fades within about a period, and a backdating window of 1-2 periods is enough.", rho_hat, eff_pre),
    "moderate" = sprintf("Moderate persistence (rho = %.2f, half-life %.1f periods): shocks fade over a few periods. The pre-treatment fit is worth about %.0f independent periods; a shock just before treatment would carry into the first post-treatment periods, so backdate by about %d periods to check for anticipation or a pre-treatment shock.", rho_hat, half_life, eff_pre, as.integer(backdate)),
    "high" = sprintf("High persistence (rho = %.2f, half-life %.1f periods): smooth, sustained drifts in the gap occur by chance, so read the gap plot with caution. The pre-treatment fit is worth about %.1f independent periods; a pre-treatment shock carries well into the post period, and backdating needs about %d periods (check that enough pre-period remains to fit on).", rho_hat, half_life, eff_pre, as.integer(backdate)),
    "very high" = sprintf("Very high persistence (rho = %.2f%s): the untreated gaps behave like a random walk, common for outcomes in levels (GDP, sales, rates that trend). The placebo calibration and the growing drift profile of the simultaneous band account for this when the treated unit is like the donors, but a close pre-treatment fit is weak evidence (about %.1f independent periods), gaps fan out with the horizon, and backdating cannot settle whether a pre-treatment shock carries over; see also the factor persistence check.", rho_hat, if (rho_hat >= 1) "; at or above 1, the placebos' gaps drift away rather than revert" else "", eff_pre)),
    if (level != "very high") " The placebo calibration accounts for persistence the treated unit shares with the donors.",
    if (!is.na(rho_in) && rho_hat - rho_in > 0.3) sprintf(" The pre-treatment residuals look much less persistent (%.2f): the weights absorb part of the noise, so the pre-treatment fit overstates how closely the synthetic control will track the treated unit.", rho_in),
    if (!is.na(seas_excess) && seas_excess > 0.2) sprintf(" Seasonal persistence: the lag-%d autocorrelation (%.2f) exceeds what lag-1 persistence implies by %.2f, so a shock in one season tends to recur a year later; backdate in whole years and consider deseasonalizing.", as.integer(s_len), rho_seas, seas_excess))
  diag$residual_autocorrelation <- list(rho = rho_hat, rho_interval_90 = unname(rho_int), level = level, half_life = half_life,
                                        backdate_periods = backdate, effective_pre_periods = eff_pre, rho_in_sample = rho_in,
                                        seasonal_lag = s_len, rho_seasonal = rho_seas, rho_seasonal_excess = seas_excess, n_post_periods = n_post,
                                        n_placebos = ncol(Zac), guide = guide)
  # factor loadings and factors (fect IFE on every cell except the treated unit's post period)
  comps <- NULL; note <- NULL
  if (!requireNamespace("fect", quietly = TRUE)) note <- "fect not installed: the hull and factor-persistence checks were not run" else {
    excl <- matrix(FALSE, nrow(Y), ncol(Y)); excl[post, tr] <- TRUE
    comps <- tryCatch({
      r_d <- if (!is.na(r_known)) r_known else if (length(unique(ife_r)) == 1) unique(ife_r) else as.integer(lto_fect(Y, Zs, excl, ife_r, TRUE, ife_seed)$r.cv)
      fit <- lto_fect(Y, Zs, excl, r_d, FALSE, ife_seed)
      L <- if (!is.null(fit$lambda)) { L <- as.matrix(fit$lambda); rownames(L) <- as.character(fit$id); L[colnames(Y), , drop = FALSE] } else NULL
      list(r = as.integer(r_d), lambda = L, factor = if (!is.null(fit$factor)) as.matrix(fit$factor) else NULL)
    }, error = function(e) { note <<- paste("fect failed:", conditionMessage(e)); NULL })
  }
  if (!is.null(comps) && !is.null(comps$lambda) && comps$r >= 1) {
    L <- comps$lambda; Fm <- comps$factor
    hull_dist <- function(target, pool) {             # min || target - pool' w || over the simplex (quadprog)
      J <- nrow(pool); if (J == 1) return(sqrt(sum((target - pool[1, ])^2)))
      D <- tcrossprod(pool); D <- D + diag(1e-10 * max(mean(diag(D)), 1e-12), J)
      wq <- tryCatch(quadprog::solve.QP(D, as.numeric(pool %*% target), cbind(rep(1, J), diag(J)), c(1, rep(0, J)), meq = 1)$solution,
                     error = function(e) rep(1 / J, J))
      wq <- pmax(wq, 0); wq <- wq / sum(wq); sqrt(sum((target - as.numeric(crossprod(pool, wq)))^2))
    }
    d_T <- hull_dist(L[tr, ], L[dn, , drop = FALSE])
    d_D <- vapply(seq_along(dn), function(j) hull_dist(L[dn[j], ], L[dn[-j], , drop = FALSE]), numeric(1))
    share_beyond <- mean(d_D >= d_T)
    delta <- L[tr, ] - colSums(L[dn, , drop = FALSE] * w)
    corr <- as.numeric(Fm %*% delta); drift <- corr - mean(corr[pre]); synth_T <- Y[, tr] - gT
    rho <- apply(Fm, 2, function(f) unname(stats::coef(stats::lm(f[-1] ~ f[-length(f)]))[2]))
    diag$hull <- list(n_factors = comps$r, treated_distance = d_T, donor_distance_median = stats::median(d_D),
                      share_donors_beyond = share_beyond, flag = d_T > 1e-8 && share_beyond < 0.10)
    diag$extrapolation <- list(max_post_drift_pct = 100 * max(abs(drift[post])) / mean(abs(synth_T[post])),
                               final_drift_pct = 100 * drift[post][n_post] / synth_T[post][n_post])
    diag$persistence <- list(factor_ar1 = rho, trending = max(rho) > 0.95)
  } else {
    diag$hull <- list(n_factors = if (is.null(comps)) NA else comps$r, flag = NA,
                      note = if (!is.null(note)) note else "no loadings (0 factors chosen)")
    diag$persistence <- list(trending = NA)
  }
  diag
}

# ---- covariate balance (as in dsc_placebo()) ----------------------------------------------------
# Z: units x covariates (pre-treatment means; rownames = units); w_T: the treated unit's weights, named by
# donor; W_P: list, one per placebo (named by placebo), of its weights named by its donors. Each gap is
# standardized by the donors' SD; share_placebos_worse is the share of placebos with a larger absolute gap.
lto_covariate_balance <- function(Z, treated, w_T, W_P) {
  syn <- function(w) vapply(colnames(Z), function(v) { z <- Z[names(w), v]; ok <- !is.na(z)
    if (!any(ok)) NA_real_ else sum(w[ok] * z[ok]) / sum(w[ok]) }, numeric(1))
  sdv <- apply(Z[names(w_T), , drop = FALSE], 2, stats::sd, na.rm = TRUE); sdv[!is.finite(sdv) | sdv == 0] <- NA
  sT <- syn(w_T); dT <- (Z[treated, ] - sT) / sdv
  dP <- matrix(vapply(names(W_P), function(p) (Z[p, ] - syn(W_P[[p]])) / sdv, numeric(ncol(Z))), nrow = ncol(Z))
  share <- vapply(seq_len(ncol(Z)), function(k) if (is.na(dT[k])) NA_real_ else mean(abs(dP[k, ]) > abs(dT[k]), na.rm = TRUE), numeric(1))
  data.frame(covariate = colnames(Z), treated = as.numeric(Z[treated, ]), synthetic = as.numeric(sT), std_diff = as.numeric(dT),
             placebo_median_abs_std_diff = apply(abs(dP), 1, stats::median, na.rm = TRUE), share_placebos_worse = share,
             flag = share < 0.10, row.names = NULL)
}

# ---- the Type-I bound and the power adjustment of Lei and Sudijono -----------------------
lto_f <- function(N, alpha) (3 - 3 / N - sqrt(9 * (1 - 1 / N)^2 - 12 * (-4 / (3 * N^2) + 1 / N + alpha * (1 - 1 / N) * (1 - 2 / N)))) / 2
#' Lei and Sudijono's Type-I bound for the LTO test
#'
#' The bound \eqn{\lfloor N f(N, \alpha) \rfloor / N} on the false-positive rate of the naive LTO test with
#' \eqn{N} units at level \eqn{\alpha}, under the sharp null and a uniformly assigned treated unit.
#' @param N Number of units (treated unit plus donors).
#' @param alpha Level.
#' @return The bound.
#' @examples
#' lto_bound(16, 0.05)   # 1/16
#' lto_bound(40, 0.10)
#' @export
lto_bound <- function(N, alpha) floor(N * lto_f(N, alpha)) / N
lto_c <- function(N, alpha) {                       # smallest c with f(N, alpha + c) at the next step
  target <- (floor(N * lto_f(N, alpha)) + 1) / N; hi <- 2 / 3 - alpha - 1e-12
  if (lto_f(N, hi) < target) return(hi)
  lo <- 0; for (i in 1:100) { m <- (lo + hi) / 2; if (lto_f(N, alpha + m) >= target) hi <- m else lo <- m }
  hi
}

# ---- the tournament ------------------------------------------------------------------------
# Y: T x units outcome matrix; treated: column index; post: logical vector of post periods.
# Returns, for every post period, the pieces needed for any theta: the treated unit's gap and
# RMSPE in each match, and the larger placebo statistic in that match (pointwise and cumulative),
# plus the standard placebo pieces (treated fitted on all donors, donors fitted without it).
lto_core <- function(Y, treated, post, rmspe_floor = 0.25, cache = NULL, c_ridge = NULL, Ylist = NULL, tournament = TRUE) {   # tournament = FALSE: standard placebo and leave-one-out fits only   # Ylist: per-unit adjusted outcomes (lto_unit_models())   # c_ridge: ridge penalty (NULL = plain DSC)   # cache: an environment shared across calls on the same Y (fits do not depend on which unit is treated)
  Y <- as.matrix(Y); n <- ncol(Y); pre <- !post; donors <- setdiff(seq_len(n), treated)
  fit_gap <- function(u, pool) { Yu <- if (is.null(Ylist)) Y else Ylist[[u]]
    f <- dsc_fit_qp(Yu[, u], Yu[, pool, drop = FALSE], pre); g <- Yu[, u] - f$synthetic
    wt <- if (treated %in% pool) f$w[match(treated, pool)] else 0
    if (!is.null(c_ridge)) { r <- lto_ridge(Yu[, u], Yu[, pool, drop = FALSE], f$w, c_ridge, pre); g <- g - r$correction
      if (treated %in% pool) wt <- wt + r$coef[match(treated, pool)] }
    attr(g, "w_treated") <- wt; g }              # w_treated: how much the gap moves (negatively) per unit added to the treated unit's post outcome
  rmspe <- function(g) sqrt(mean(g[pre]^2))
  # leave-one-out fits of every unit on all the others. Pre-treatment only, they give the RMSPE
  # floor (rmspe_floor x the median pre-treatment RMSPE), which does not depend on which unit is
  # treated; their post-treatment gaps give the simultaneous band's drift profile (lto_profile()).
  loo <- lapply(seq_len(n), function(u) fit_gap(u, setdiff(seq_len(n), u)))
  s_all <- vapply(loo, rmspe, numeric(1))
  fl <- if (rmspe_floor > 0) rmspe_floor * stats::median(s_all) else 0
  stat <- function(g) { s <- max(rmspe(g), fl); gp <- g[post]
    list(s = s, k = abs(gp) / s, kc = abs(cumsum(gp)) / s, ks = gp / s, kcs = cumsum(gp) / s, g = gp, gc = cumsum(gp), ratio = sqrt(mean(gp^2)) / rmspe(g)) }   # ks, kcs: signed (one-sided tests)
  # standard placebo (as in dsc_placebo): treated on all donors, each donor without the treated unit
  gT <- fit_gap(treated, donors); sT <- stat(gT)
  plac_g <- lapply(donors, function(j) fit_gap(j, setdiff(donors, j))); plac <- lapply(plac_g, stat)
  # the tournament: unit u fitted without the other two members of each match
  pairs <- if (tournament) utils::combn(donors, 2) else matrix(integer(0), 2, 0); P <- ncol(pairs)
  if (is.null(cache)) cache <- new.env()
  member <- function(u, a, b) {                      # u fitted with pool = all \ {u, a, b}
    key <- paste(c(if (!is.null(Ylist)) paste0("T", treated), u, sort(c(a, b))), collapse = "_")   # adjusted fits depend on the treated unit
    if (is.null(cache[[key]])) cache[[key]] <- stat(fit_gap(u, setdiff(seq_len(n), c(u, a, b))))
    cache[[key]]
  }
  nP <- sum(post)
  G <- M <- Gc <- Mc <- Ki <- Kj <- Kci <- Kcj <- Ksi <- Ksj <- Kcsi <- Kcsj <- matrix(NA_real_, nP, P)
  S <- RatioI <- RatioM <- numeric(P)
  for (p in seq_len(P)) {
    i <- pairs[1, p]; j <- pairs[2, p]
    I <- member(treated, i, j); Ai <- member(i, treated, j); Aj <- member(j, treated, i)
    G[, p] <- I$g; Gc[, p] <- I$gc; S[p] <- I$s
    Ki[, p] <- Ai$k; Kj[, p] <- Aj$k; Kci[, p] <- Ai$kc; Kcj[, p] <- Aj$kc
    Ksi[, p] <- Ai$ks; Ksj[, p] <- Aj$ks; Kcsi[, p] <- Ai$kcs; Kcsj[, p] <- Aj$kcs
    M[, p] <- pmax(Ai$k, Aj$k); Mc[, p] <- pmax(Ai$kc, Aj$kc)
    RatioI[p] <- I$ratio; RatioM[p] <- max(Ai$ratio, Aj$ratio)
  }
  list(G = G, Gc = Gc, S = S, M = M, Mc = Mc, Ki = Ki, Kj = Kj, Kci = Kci, Kcj = Kcj, Ksi = Ksi, Ksj = Ksj, Kcsi = Kcsi, Kcsj = Kcsj,
       ratio_I = RatioI, ratio_M = RatioM, pairs = pairs, treated = treated, donors = donors,
       n_units = n, n_pairs = P, floor = fl,
       loo = list(g = sapply(loo, function(g) g[post]), s = pmax(s_all, fl),
                  w_treated = vapply(loo, function(g) attr(g, "w_treated"), numeric(1))),
       standard = list(g = sT$g, gc = sT$gc, s = sT$s, k = sapply(plac, `[[`, "k"), kc = sapply(plac, `[[`, "kc"),
                       gaps = sapply(plac_g, as.numeric),   # each donor's gap in every period, fitted on the other donors
                       gs = matrix(sapply(plac, `[[`, "g"), nrow = nP),   # each donor's post-period gaps (signed)
                       s_donors = sapply(plac, `[[`, "s"), ratio = sT$ratio, ratio_donors = sapply(plac, `[[`, "ratio")),
       synthetic_treated = Y[, treated] - gT, gap_treated = gT, c_ridge = c_ridge)
}

# Add a constant (or per-period) effect tau to the treated unit's post-treatment outcomes of a
# lto_core() object without refitting: every fit uses pre-treatment data only, so the treated
# unit's gaps shift by tau, and every other unit's leave-one-out gap shifts by -(its weight on
# the treated unit) x tau; with ridge augmentation, -(weight + ridge coefficient on the treated
# unit) x tau, since the correction is linear in the pool's post outcomes. Used by the simulations.
lto_shift <- function(core, tau) {
  nP <- nrow(core$G); tt <- rep_len(tau, nP); ct <- cumsum(tt)
  if (ncol(core$G)) { core$G <- core$G + tt; core$Gc <- core$Gc + ct }   # no tournament (lto_core(tournament = FALSE)): nothing to shift
  core$standard$g <- core$standard$g + tt; core$standard$gc <- core$standard$gc + ct
  w <- core$loo$w_treated; w[core$treated] <- -1                           # the treated unit's own gap shifts by +tau
  core$loo$g <- core$loo$g - outer(tt, w)
  core
}

# Drift profiles for the simultaneous band (as in dsc_placebo(): the median standardized gap in
# each period, smoothed by local linear loess on the log scale). For the LTO test the profile is
# computed from the leave-one-out fits of ALL units, so it is the same whichever unit is treated
# (dsc_placebo() uses the donors only).
smooth_profile <- function(v, span = 0.75) { v <- pmax(v, 1e-8); tt <- seq_along(v)
  if (length(v) < 4) return(v)
  exp(as.numeric(stats::predict(stats::loess(log(v) ~ tt, span = span, degree = 1)))) }
lto_profile <- function(core, span = 0.75) {
  K <- abs(core$loo$g) / rep(core$loo$s, each = nrow(core$loo$g))
  Kc <- abs(apply(core$loo$g, 2, cumsum)) / rep(core$loo$s, each = nrow(core$loo$g))
  if (is.null(dim(Kc))) Kc <- matrix(Kc, nrow = 1)
  list(q = smooth_profile(apply(K, 1, stats::median), span), qc = smooth_profile(apply(Kc, 1, stats::median), span))
}

# p-values for a given effect theta (vector over periods; pointwise and cumulative), from lto_core()
lto_pvalues <- function(core, theta = 0, theta_cum = cumsum(rep_len(theta, nrow(core$G)))) {
  th <- rep_len(theta, nrow(core$G)); thc <- rep_len(theta_cum, nrow(core$G))
  kI <- abs(core$G - th) / rep(core$S, each = nrow(core$G)); kIc <- abs(core$Gc - thc) / rep(core$S, each = nrow(core$G))
  st <- core$standard; N_d <- ncol(st$k)
  kT <- abs(st$g - th) / st$s; kTc <- abs(st$gc - thc) / st$s
  data.frame(p_lto = rowMeans(kI <= core$M), p_lto_cum = rowMeans(kIc <= core$Mc),
             p_placebo = (1 + rowSums(st$k >= kT)) / (N_d + 1), p_placebo_cum = (1 + rowSums(st$kc >= kTc)) / (N_d + 1))
}

# One-sided p-values from lto_core(), for the direction stated in advance ("greater": the effect is
# predicted to raise the outcome; "less": to lower it). Every test ranks the signed standardized gap
# (d x gap / RMSPE, d = +1 or -1) instead of its absolute value, so only gaps in the predicted direction
# count against the treated unit: placebo rank tests (p = (1 + donors at least as extreme in that
# direction) / (donors + 1)) and the naive LTO test (share of matches the treated unit does not win,
# with the placebos' signed scores). Simultaneous: each unit's score is max_t d x k_t / q_t with the same
# drift profiles as the two-sided tests. The smallest possible p-values are those of the two-sided tests.
lto_pvalues_onesided <- function(core, direction = c("greater", "less"), span = 0.75, prof = lto_profile(core, span)) {
  direction <- match.arg(direction); d <- if (direction == "greater") 1 else -1
  st <- core$standard; Nd <- ncol(st$k); nP <- length(st$g)
  sd_rep <- rep(st$s_donors, each = nP)
  kT <- d * st$g / st$s; kTc <- d * st$gc / st$s
  kD <- d * st$gs / sd_rep; kDc <- d * matrix(apply(st$gs, 2, cumsum), nrow = nP) / sd_rep
  q0 <- smooth_profile(apply(matrix(st$k, nrow = nP), 1, stats::median), span)
  q0c <- smooth_profile(apply(matrix(st$kc, nrow = nP), 1, stats::median), span)
  out <- list(p_placebo = (1 + rowSums(kD >= kT)) / (Nd + 1), p_placebo_cum = (1 + rowSums(kDc >= kTc)) / (Nd + 1),
              p_placebo_sim = (1 + sum(apply(kD / q0, 2, max) >= max(kT / q0))) / (Nd + 1),
              p_placebo_sim_cum = (1 + sum(apply(kDc / q0c, 2, max) >= max(kTc / q0c))) / (Nd + 1),
              p_lto = rep(NA_real_, nP), p_lto_cum = rep(NA_real_, nP), p_lto_sim = NA_real_, p_lto_sim_cum = NA_real_)
  if (ncol(core$G)) {
    S <- rep(core$S, each = nP)
    out$p_lto <- rowMeans(d * core$G / S <= pmax(d * core$Ksi, d * core$Ksj))
    out$p_lto_cum <- rowMeans(d * core$Gc / S <= pmax(d * core$Kcsi, d * core$Kcsj))
    q <- prof$q; qc <- prof$qc
    out$p_lto_sim <- mean(apply(d * core$G / outer(q, core$S), 2, max) <= pmax(apply(d * core$Ksi / q, 2, max), apply(d * core$Ksj / q, 2, max)))
    out$p_lto_sim_cum <- mean(apply(d * core$Gc / outer(qc, core$S), 2, max) <= pmax(apply(d * core$Kcsi / qc, 2, max), apply(d * core$Kcsj / qc, 2, max)))
  }
  out
}

# confidence set from n intervals [lo, hi]: theta covered by more than `need` intervals
cover_set <- function(lo, hi, need) {
  if (sum(is.finite(lo)) <= need) return(c(lower = Inf, upper = -Inf, connected = NA))
  ev <- rbind(data.frame(x = lo, d = 1), data.frame(x = hi, d = -1)); ev <- ev[order(ev$x, -ev$d), ]
  cnt <- cumsum(ev$d); inside <- cnt > need
  if (!any(inside)) return(c(lower = Inf, upper = -Inf, connected = NA))
  starts <- which(inside & !c(FALSE, head(inside, -1))); ends <- which(inside & !c(tail(inside, -1), FALSE)) + 1
  c(lower = ev$x[starts[1]], upper = ev$x[ends[length(ends)]], connected = length(starts) == 1)
}

# intervals for every post period from lto_core(); version "naive" or "powered"
lto_intervals <- function(core, conf_level = 0.95, version = c("naive", "powered"), delta = 1e-10) {
  version <- match.arg(version); alpha <- round(1 - conf_level, 10); N <- core$n_units
  thr <- if (version == "powered") alpha + lto_c(N, alpha) - delta else alpha
  need <- floor(thr * core$n_pairs + 1e-9)          # theta is rejected when at most `need` matches are not won
  sc <- core$S                                         # the treated unit's RMSPE in each match
  do.call(rbind, lapply(seq_len(nrow(core$G)), function(t) {
    a <- cover_set(core$G[t, ] - sc * core$M[t, ], core$G[t, ] + sc * core$M[t, ], need)
    b <- cover_set(core$Gc[t, ] - sc * core$Mc[t, ], core$Gc[t, ] + sc * core$Mc[t, ], need)
    data.frame(eff_low = a[["lower"]], eff_high = a[["upper"]], eff_connected = as.logical(a[["connected"]]),
               cum_eff_low = b[["lower"]], cum_eff_high = b[["upper"]], cum_connected = as.logical(b[["connected"]]))
  }))
}

# ---- simultaneous test and band ----------------------------------------------------------------
# Each unit's score is its largest standardized gap relative to the drift profile, max_t k_t / q_t
# (cumulative: max_t kc_t / qc_t). In each match the treated unit fails to win for the effect
# path theta exactly when |G_tp - theta_t| <= S_p q_t M_p in every period, with M_p the larger
# placebo score in the match: a box. The 1 - alpha confidence set for the path is the set of paths
# inside more than `need` boxes; the band reported for period t is the set of theta_t covered by
# more than `need` of the boxes' sides in that period, which contains the confidence set's
# projection, so the band covers the whole path whenever the set does.
lto_simultaneous <- function(core, conf_level = 0.95, version = c("naive", "powered"), span = 0.75, prof = lto_profile(core, span)) {
  version <- match.arg(version); alpha <- round(1 - conf_level, 10); N <- core$n_units
  thr <- if (version == "powered") alpha + lto_c(N, alpha) - 1e-10 else alpha
  need <- floor(thr * core$n_pairs + 1e-9)
  q <- prof$q; qc <- prof$qc
  Mi <- pmax(apply(core$Ki / q, 2, max), apply(core$Kj / q, 2, max))
  Mic <- pmax(apply(core$Kci / qc, 2, max), apply(core$Kcj / qc, 2, max))
  RI <- apply(abs(core$G) / outer(q, core$S), 2, max); RIc <- apply(abs(core$Gc) / outer(qc, core$S), 2, max)
  band <- do.call(rbind, lapply(seq_along(q), function(t) {
    a <- cover_set(core$G[t, ] - core$S * q[t] * Mi, core$G[t, ] + core$S * q[t] * Mi, need)
    b <- cover_set(core$Gc[t, ] - core$S * qc[t] * Mic, core$Gc[t, ] + core$S * qc[t] * Mic, need)
    data.frame(eff_low_sim = a[["lower"]], eff_high_sim = a[["upper"]], cum_eff_low_sim = b[["lower"]], cum_eff_high_sim = b[["upper"]])
  }))
  list(p_sim = mean(RI <= Mi), p_sim_cum = mean(RIc <= Mic), reject = mean(RI <= Mi) <= thr + 1e-12,
       reject_cum = mean(RIc <= Mic) <= thr + 1e-12, flagged = which(band$eff_low_sim > 0 | band$eff_high_sim < 0),
       band = band, profile = prof)
}

# the standard simultaneous band of dsc_placebo() from the same fits (donor-only profile, rank rule)
standard_simultaneous <- function(core, conf_level = 0.95, span = 0.75) {
  st <- core$standard; N_d <- ncol(st$k); m <- floor(round(1 - conf_level, 10) * (N_d + 1))
  one <- function(K, gT) { q <- smooth_profile(apply(K, 1, stats::median), span); R <- apply(K / q, 2, max)
    lam <- if (m >= 1) sort(R, decreasing = TRUE)[m] else Inf; RT <- max(abs(gT) / (st$s * q))
    list(low = gT - lam * q * st$s, high = gT + lam * q * st$s, p = (1 + sum(R >= RT)) / (N_d + 1)) }
  a <- one(st$k, st$g); b <- one(st$kc, st$gc)
  list(p_sim = a$p, p_sim_cum = b$p, band = data.frame(eff_low_sim = a$low, eff_high_sim = a$high, cum_eff_low_sim = b$low, cum_eff_high_sim = b$high))
}

# ---- sensitivity to non-uniform assignment (Lei and Sudijono, Section 3.2) ------------------------
# notwin: logical vector over the matches (columns of core$pairs) - TRUE when the treated unit does
# not beat both placebos. With assignment probabilities pi_i in [1/(Gamma N), Gamma/N] summing to
# one, the weighted LTO p-value weights match (j, k) by pi_j pi_k. lto_max_weighted_p() finds its
# largest value over that set (a non-convex problem). The weighted p-value depends only on the donors'
# relative weights and is linear-fractional in each one, so near its maximum each donor sits at a low or
# a high weight: the search runs over which donors are high (exhaustive for up to 14 donors, greedy
# starts plus single-donor flips beyond) and the high/low ratio the constraints allow, then polishes
# the best candidates one weight at a time. In checks against continuous optimisation it is within
# 1e-4 of the best value found.
lto_max_weighted_p <- function(notwin, pairs, donors, N, Gamma, exhaustive_max = 14, n_random = 200, random_seed = 1) {
  Nd <- length(donors); Gm <- matrix(0, Nd, Nd); id <- match(pairs, donors)
  Gm[cbind(id[c(TRUE, FALSE)], id[c(FALSE, TRUE)])] <- notwin; Gm <- Gm + t(Gm)
  a0 <- sum(Gm); base <- a0 / (Nd^2 - Nd); if (Gamma <= 1) return(base)
  l <- 1 / (Gamma * N); u <- Gamma / N
  rmax <- function(k) {                               # largest high/low weight ratio with k high donors
    feas <- function(r) { d <- k * r + Nd - k; max(l, (1 - u) / d) <= min(u / r, (1 - l) / d) + 1e-15 }
    lo <- 1; hi <- Gamma^2; if (feas(hi)) return(hi)
    for (i in 1:50) { m <- (lo + hi) / 2; if (feas(m)) lo <- m else hi <- m }; lo }
  rm <- vapply(0:Nd, function(k) if (k %in% c(0, Nd)) 1 else rmax(k), numeric(1))
  # two-level weights v = 1 + (r - 1) H: the weighted p-value in closed form from three sums
  value <- function(k, a1, a2) { if (k == 0 || k == Nd || rm[k + 1] <= 1) return(c(base, 1))
    r <- seq(1, rm[k + 1], length.out = 41); d <- r - 1
    v <- (a0 + 2 * d * a1 + d^2 * a2) / ((Nd + d * k)^2 - (Nd + (r^2 - 1) * k)); c(max(v), r[which.max(v)]) }
  rs <- rowSums(Gm)
  evalH <- function(H) value(sum(H), sum(rs[H]), sum(Gm[H, H]))
  cand <- list()
  if (Nd <= exhaustive_max) {
    for (code in 1:(2^Nd - 2)) { H <- bitwAnd(code, 2^(0:(Nd - 1))) > 0; cand[[length(cand) + 1]] <- list(H = H, v = evalH(H)) }
  } else {
    ord <- order(rs, decreasing = TRUE)
    for (k in 1:(Nd - 1)) {                          # greedy start: the k donors the treated unit most often fails to beat
      H <- seq_len(Nd) %in% ord[1:k]; cur <- evalH(H)
      repeat { improved <- FALSE                      # single-donor flips
        for (j in seq_len(Nd)) { H2 <- H; H2[j] <- !H2[j]; if (sum(H2) %in% c(0, Nd)) next
          v2 <- evalH(H2); if (v2[1] > cur[1] + 1e-12) { H <- H2; cur <- v2; improved <- TRUE } }
        if (!improved) break }
      cand[[length(cand) + 1]] <- list(H = H, v = cur)
    }
  }
  vals <- vapply(cand, function(x) x$v[1], numeric(1)); best <- max(base, vals)
  # polish the best candidates: the weighted p-value is linear-fractional in each donor's weight, so
  # move one weight at a time to an end of its feasible range (the treated unit's weight absorbs the slack)
  pw <- function(v) as.numeric(crossprod(v, Gm %*% v)) / (sum(v)^2 - sum(v^2))
  for (ci in utils::head(order(vals, decreasing = TRUE), 10)) {
    v <- ifelse(cand[[ci]]$H, cand[[ci]]$v[2], 1); S <- sum(v)
    a <- (max(l / min(v), (1 - u) / S) + min(u / max(v), (1 - l) / S)) / 2; v <- a * v
    for (sweep in 1:50) { improved <- FALSE
      for (j in seq_len(Nd)) { rest <- sum(v) - v[j]
        for (x in c(max(l, 1 - u - rest), min(u, 1 - l - rest))) { w <- v; w[j] <- x
          if (pw(w) > pw(v) + 1e-13) { v <- w; improved <- TRUE } } }
      if (!improved) break }
    best <- max(best, pw(v))
  }
  # random continuous starts, each polished one weight at a time (lto_sensitivity_search_check.R: the high/low
  # search above can fall short, by up to 0.04 at Gamma = 3 and about 0.0003 near the reported Gamma). These do
  # not assume the high/low structure; the session's random-number state is left untouched.
  if (n_random > 0) best <- max(best, lto_with_seed(random_seed, {
    rb <- base
    for (s in seq_len(n_random)) {
      w0 <- if (s %% 2) exp(stats::runif(Nd, 0, 2 * log(Gamma))) else ifelse(stats::runif(Nd) < stats::runif(1), Gamma^2, 1)
      S0 <- sum(w0); lo_a <- max(l / min(w0), (1 - u) / S0); hi_a <- min(u / max(w0), (1 - l) / S0)
      if (lo_a > hi_a) next
      rb <- max(rb, lto_polish_fast(stats::runif(1, lo_a, hi_a) * w0, Gm, l, u))
    }
    rb }))
  best
}
# coordinate ascent for the weighted p-value v'Gv / (sum(v)^2 - sum(v^2)), each weight moved to an end of its
# feasible range in turn (the treated unit's weight absorbs the slack); Gv, sum(v) and sum(v^2) are updated in
# place, so a move costs O(1) to evaluate. Same moves as the polish in lto_max_weighted_p().
lto_polish_fast <- function(v, Gm, l, u) {
  Gv <- as.numeric(Gm %*% v); num <- sum(v * Gv); S <- sum(v); Q <- sum(v^2); val <- num / (S^2 - Q)
  for (sweep in 1:50) { improved <- FALSE
    for (j in seq_along(v)) { rest <- S - v[j]
      for (x in c(max(l, 1 - u - rest), min(u, 1 - l - rest))) { dlt <- x - v[j]; if (abs(dlt) < 1e-15) next
        num2 <- num + 2 * dlt * Gv[j]; S2 <- S + dlt; Q2 <- Q - v[j]^2 + x^2; val2 <- num2 / (S2^2 - Q2)
        if (val2 > val + 1e-13) { Gv <- Gv + dlt * Gm[, j]; num <- num2; S <- S2; Q <- Q2; v[j] <- x; val <- val2; improved <- TRUE } } }
    if (!improved) break }
  val
}

# smallest Gamma at which the largest weighted p-value exceeds alpha (Inf if none up to gamma_max)
lto_sensitivity_gamma <- function(notwin, pairs, donors, N, alpha, gamma_max = 3) {
  f <- function(g) lto_max_weighted_p(notwin, pairs, donors, N, g) > alpha
  if (mean(notwin) > alpha) return(1)                 # not significant to begin with
  if (!f(gamma_max)) return(Inf)
  lo <- 1; hi <- gamma_max; for (i in 1:20) { m <- (lo + hi) / 2; if (f(m)) hi <- m else lo <- m }
  hi
}

#' Leave-two-out placebo inference for the demeaned synthetic control
#'
#' The leave-two-out (LTO) placebo test of Lei and Sudijono (2025) applied to the demeaned synthetic control,
#' for panels with few donors. The estimator, the RMSPE-scaled statistic, the RMSPE floor and the drift profile
#' are those of [dsc_placebo()]; only the calibration changes. The test plays a tournament: for every pair of
#' donors \eqn{(i, j)}, the treated unit and the two donors are each fitted with every *other* unit as donors,
#' and the treated unit "wins" the match if its statistic is larger than both of theirs. The naive p-value is
#' the share of matches it does not win, so it can reject at 5% with far fewer than 19 donors and gives finite
#' 95% intervals.
#'
#' @details
#' **Validity.** Under the sharp null and a uniformly assigned treated unit among the \eqn{N} units,
#' \eqn{P(p \le \alpha) \le \lfloor N f(N, \alpha) \rfloor / N} (Lei and Sudijono, Theorems 2.1-2.3). The bound
#' can exceed \eqn{\alpha} (it is \eqn{1/N} below 10 units); the simulations in the repository find the actual
#' false-positive rate close to or below \eqn{\alpha} on real panels with 10-20 units.
#'
#' **Sensitivity.** `gamma_sensitivity` is the smallest \eqn{\Gamma} at which the worst-case p-value, with
#' assignment probabilities allowed to range over \eqn{[1/(\Gamma N), \Gamma/N]}, exceeds \eqn{\alpha}
#' (Lei and Sudijono, Section 3.2). It measures how unequal assignment would have to be to explain the result
#' away, not robustness to chance.
#'
#' **Tests of a sharp null.** The p-values test no effect (or a constant effect, when inverted) on the treated
#' unit; the intervals can be read as prediction intervals for the realized effect.
#'
#' **Run time.** Seconds for up to about 20 donors (\eqn{3 N_d (N_d - 1)/2} fits); Kansas (48 donors) takes
#' about 40 seconds. With `adjust_model = "ife"` or `augment = "ife"` every unit needs its own `fect` model
#' (Kansas with three covariates: about 11 minutes, 3.5 with `n_cores = 4`).
#'
#' @inheritParams dsc_placebo
#' @param conf_level Confidence level; the tests reject at \eqn{1 -} `conf_level`.
#' @param version `"naive"` (the default: at or below nominal in the exact enumeration) or `"powered"`
#'   (Lei and Sudijono's threshold \eqn{\alpha + c(N, \alpha)}; more liberal, not recommended below 10 units).
#' @param rmspe_floor Floor on each unit's pre-treatment RMSPE, as a share of the median over all units.
#' @param span LOESS span of the simultaneous band's drift profile.
#' @param sensitivity Compute the sensitivity \eqn{\Gamma} for each test.
#' @param gamma_max Largest \eqn{\Gamma} searched (`Inf` is returned when the result survives every value up to
#'   it).
#' @param augment `"none"`, `"ridge"` (every fit in the tournament is ridge-augmented, as in
#'   `dsc_placebo(augment = "ridge")`; recommended by the diagnostics for an extreme treated unit with trending
#'   factors) or `"ife"` (IFE-augmented; needs `fect`).
#' @param ridge_c Ridge penalty (default 0.1); `"cv"` or `NULL` chooses it once by cross-validation over all
#'   units' pre-periods, so it does not depend on which unit is treated.
#' @param ridge_grid,ridge_rule,ridge_windows Grid, rule (`"min"` or `"1se"`) and number of holdout windows of
#'   that cross-validation; used only with `ridge_c = "cv"`.
#' @param diagnostics Run the same checks as [dsc_placebo()]'s diagnostics from the tournament's own fits (the
#'   hull and factor-persistence checks need `fect`), and recommend `augment = "ridge"` only for a treated unit
#'   flagged as extreme with trending factors.
#' @param warn Warn when the diagnostics recommend ridge and it is not used, when `fect` is missing, and with
#'   fewer than 10 units when `version = "powered"` or the level is stricter than 10%.
#' @param adjust_covariates,adjust_model Time-varying covariate adjustment as in [dsc_placebo()]: the effect is
#'   estimated once per unit, leaving out that unit's and the treated unit's post periods, by interactive
#'   (`"ife"`, needs `fect`) or two-way (`"twfe"`) fixed effects. Combines with `augment = "ridge"`.
#' @param n_cores Cores for the per-unit `fect` models (results identical to `n_cores = 1`).
#'
#' @return An object of class `"dsc_placebo_lto"`: a list with
#' * `effects`: one row per period with `observed`, `synthetic`, `effect` and, after treatment, the LTO effect
#'   interval (`eff_low`, `eff_high`), the prediction interval (`pred_low`, `pred_high`), the cumulative
#'   versions, the simultaneous versions (`*_sim`) and the naive LTO and standard placebo p-values;
#' * `test`: the pointwise and cumulative tests at the test period and the two simultaneous tests, with the
#'   LTO p-value, rejection, the standard placebo p-value and `gamma_sensitivity`;
#' * `flagged_periods` (where the simultaneous band excludes zero), `weights`, `diagnostics`, `balance`,
#'   `adjustment`;
#' * `bound`: the Type-I bound, the approximate- and exact-placebo bounds and the rejection threshold;
#' * `settings` and `core` (the tournament, for further tests without refitting).
#' @references Lei, L., and Sudijono, T. (2025). Inference for synthetic controls via refined placebo tests.
#'   arXiv:2401.07152.
#' @seealso [dsc_placebo()], [synth_power()], [summary.dsc_placebo_lto()], [plot.dsc_placebo_lto()].
#' @examples
#' \donttest{
#' if (requireNamespace("scpi", quietly = TRUE)) {
#'   data("scpi_germany", package = "scpi")
#'   germany <- as.data.frame(scpi_germany)[, c("country", "year", "gdp")]
#'   fit <- dsc_placebo_lto(germany, "country", "year", "gdp", treated_unit = "West Germany",
#'                          treatment_period = c(year = 1991), diagnostics = FALSE)
#'   fit
#'   plot(fit, type = "effect")
#' }
#' }
#' @export
dsc_placebo_lto <- function(data, unit_col, time_cols, outcome_col, treated_unit, treatment_period,
                            test_period = NULL, conf_level = 0.95, version = c("naive", "powered"), rmspe_floor = 0.25,
                            span = 0.75, sensitivity = TRUE, gamma_max = 3, augment = c("none", "ridge", "ife"), ridge_c = 0.1,
                            ridge_grid = 10^seq(-3, 2, by = 0.5), ridge_rule = "min", ridge_windows = 1,
                            diagnostics = TRUE, warn = TRUE, adjust_covariates = NULL, adjust_model = c("ife", "twfe"),
                            ife_r = 0:5, ife_seed = 12112025, balance_covariates = adjust_covariates, n_cores = 1, verbose = interactive()) {
  version <- match.arg(version); augment <- match.arg(augment, c("none", "ridge", "ife")); adjust_model <- match.arg(adjust_model)
  d <- as.data.frame(data); d <- d[do.call(order, d[time_cols]), , drop = FALSE]
  key <- do.call(paste, c(d[time_cols], sep = "_"))
  wide <- stats::reshape(data.frame(key = key, unit = d[[unit_col]], y = d[[outcome_col]]),
                         idvar = "key", timevar = "unit", direction = "wide")
  periods <- unique(d[time_cols]); rownames(periods) <- NULL
  Y <- as.matrix(wide[match(unique(key), wide$key), -1, drop = FALSE]); colnames(Y) <- sub("^y\\.", "", colnames(Y))
  if (anyNA(Y)) stop("The panel must be balanced with no missing outcomes.")
  tr <- match(treated_unit, colnames(Y)); if (is.na(tr)) stop("treated_unit not found.")
  after <- function(tp) vapply(seq_len(nrow(periods)), function(r) {   # lexicographic: on or after tp
    for (v in names(tp)) if (periods[[v]][r] != tp[[v]]) return(periods[[v]][r] > tp[[v]])
    TRUE }, logical(1))
  post <- after(treatment_period)
  test_row <- if (is.null(test_period)) sum(post) else {
    hit <- which(Reduce(`&`, lapply(names(test_period), function(v) periods[[v]] == test_period[[v]])))
    if (length(hit) != 1 || !post[hit]) stop("test_period must be one post-treatment period."); hit - sum(!post) }
  if (identical(ridge_c, "cv")) ridge_c <- NULL                  # "cv" (or NULL): choose the penalty by cross-validation
  if (!is.null(ridge_c) && (!is.numeric(ridge_c) || length(ridge_c) != 1 || !is.finite(ridge_c) || ridge_c <= 0)) stop("ridge_c must be one positive number, or \"cv\" to choose it by cross-validation.")
  N_d <- ncol(Y) - 1
  if (N_d < 3) stop("At least 3 donors are needed.")
  c_r <- if (augment == "ridge") { if (is.null(ridge_c)) lto_choose_ridge(Y, !post, ridge_grid, ridge_rule, ridge_windows) else ridge_c } else NULL
  # covariates (wide, same rows and columns as Y) and the per-unit models
  Zs <- setNames(lapply(adjust_covariates, function(v) {
    w <- stats::reshape(data.frame(key = key, unit = d[[unit_col]], z = d[[v]]), idvar = "key", timevar = "unit", direction = "wide")
    m <- as.matrix(w[match(unique(key), w$key), -1, drop = FALSE]); colnames(m) <- sub("^z\\.", "", colnames(m)); m[, colnames(Y), drop = FALSE] }), adjust_covariates)
  if (length(Zs) && any(vapply(Zs, anyNA, logical(1)))) stop("adjust_covariates must have no missing values.")
  if (length(balance_covariates)) {
    if (!all(balance_covariates %in% names(d))) stop("balance_covariates not found in data: ", paste(setdiff(balance_covariates, names(d)), collapse = ", "))
    if (!all(vapply(d[balance_covariates], is.numeric, logical(1)))) stop("balance_covariates must be numeric columns.")
  }
  Zbal <- if (length(balance_covariates)) {             # pre-treatment mean of each balance covariate, units x covariates
    Z <- do.call(cbind, lapply(balance_covariates, function(v) {
      w <- stats::reshape(data.frame(key = key, unit = d[[unit_col]], z = d[[v]]), idvar = "key", timevar = "unit", direction = "wide")
      m <- as.matrix(w[match(unique(key), w$key), -1, drop = FALSE]); colnames(m) <- sub("^z\\.", "", colnames(m))
      colMeans(m[!post, colnames(Y), drop = FALSE], na.rm = TRUE) }))
    dimnames(Z) <- list(colnames(Y), balance_covariates); Z } else NULL
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) stop("n_cores must be a single number >= 1.")
  um <- if (augment == "ife") lto_unit_models(Y, Zs, tr, post, "augment", ife_r = ife_r, ife_seed = ife_seed, n_cores = n_cores, verbose = isTRUE(verbose)) else
        if (length(Zs)) lto_unit_models(Y, Zs, tr, post, "adjust", adjust_model, ife_r, ife_seed, n_cores = n_cores, verbose = isTRUE(verbose)) else NULL
  if (!is.null(um) && augment == "ridge" && is.null(ridge_c)) c_r <- lto_choose_ridge(um$Ylist[[tr]], !post, ridge_grid, ridge_rule, ridge_windows)
  core <- lto_core(Y, tr, post, rmspe_floor, c_ridge = c_r, Ylist = um$Ylist)
  ci <- lto_intervals(core, conf_level, version)
  pv <- lto_pvalues(core, 0)
  sim <- lto_simultaneous(core, conf_level, version, span)
  obs <- Y[post, tr]; col <- function(v) c(rep(NA, sum(!post)), v)
  effects <- cbind(periods, data.frame(post = post, observed = Y[, tr], synthetic = core$synthetic_treated, effect = core$gap_treated,
    eff_low = col(ci$eff_low), eff_high = col(ci$eff_high), pred_low = col(obs - ci$eff_high), pred_high = col(obs - ci$eff_low),
    cum_effect = col(cumsum(core$gap_treated[post])), cum_eff_low = col(ci$cum_eff_low), cum_eff_high = col(ci$cum_eff_high),
    eff_low_sim = col(sim$band$eff_low_sim), eff_high_sim = col(sim$band$eff_high_sim),
    pred_low_sim = col(obs - sim$band$eff_high_sim), pred_high_sim = col(obs - sim$band$eff_low_sim),
    cum_eff_low_sim = col(sim$band$cum_eff_low_sim), cum_eff_high_sim = col(sim$band$cum_eff_high_sim),
    p_lto = col(pv$p_lto), p_lto_cum = col(pv$p_lto_cum), p_placebo = col(pv$p_placebo), p_placebo_cum = col(pv$p_placebo_cum)))
  alpha <- round(1 - conf_level, 10); N <- ncol(Y); thr <- if (version == "powered") alpha + lto_c(N, alpha) - 1e-10 else alpha
  # small pools (exact_size_small_n.R): below 10 units the powered version over-rejects, and the naive test does not attain 5%
  small_pool_note <- if (N >= 10) NULL else c(
    if (version == "powered") paste0("with ", N, " units the powered version over-rejected in validation (16% at 5 units and 10-12% at 6-7 units at the 5% level); use version = \"naive\""),
    if (alpha < 0.10 - 1e-9) paste0("with ", N, " units the naive test's exact false-positive rate at the ", 100 * alpha, "% level was 5-9% in validation ",
                                    "(about 9% at 5 units, where 5% and 10% reject the same units); interpret the result at the 10% level and report the p-value"))
  if (warn && length(small_pool_note)) warning("dsc_placebo_lto: ", paste(small_pool_note, collapse = "; "), ".", call. = FALSE)
  e <- effects[effects$post, ][test_row, ]
  stdsim <- standard_simultaneous(core, conf_level, span)
  pvec <- c(e$p_lto, e$p_lto_cum, sim$p_sim, sim$p_sim_cum)
  test <- data.frame(test = c("pointwise", "pointwise", "simultaneous", "simultaneous"),
                     measure = c("effect", "cumulative effect", "effect path", "cumulative effect path"),
                     estimate = c(e$effect, e$cum_effect, NA, NA), lower = c(e$eff_low, e$cum_eff_low, NA, NA),
                     upper = c(e$eff_high, e$cum_eff_high, NA, NA),
                     p_lto_naive = pvec, reject_lto = pvec <= thr + 1e-12,
                     p_placebo = c(e$p_placebo, e$p_placebo_cum, stdsim$p_sim, stdsim$p_sim_cum))
  if (sensitivity) {
    # which matches the treated unit fails to win, for each test at effect zero
    h <- test_row; kI <- abs(core$G[h, ]) / core$S; kIc <- abs(core$Gc[h, ]) / core$S
    q <- sim$profile$q; qc <- sim$profile$qc
    Mi <- pmax(apply(core$Ki / q, 2, max), apply(core$Kj / q, 2, max)); Mic <- pmax(apply(core$Kci / qc, 2, max), apply(core$Kcj / qc, 2, max))
    notwin <- list(kI <= core$M[h, ], kIc <= core$Mc[h, ], apply(abs(core$G) / outer(q, core$S), 2, max) <= Mi,
                   apply(abs(core$Gc) / outer(qc, core$S), 2, max) <= Mic)
    test$gamma_sensitivity <- vapply(notwin, function(nw) lto_sensitivity_gamma(nw, core$pairs, core$donors, ncol(Y), alpha, gamma_max), numeric(1))
  }
  Yw <- if (is.null(um)) Y else um$Ylist[[tr]]
  w <- dsc_fit_qp(Yw[, tr], Yw[, -tr, drop = FALSE], !post)$w
  # Covariate balance: the outcome-fitted weights applied to the donors' covariates, against the standard placebos'
  balance <- if (!is.null(Zbal)) {
    dn <- core$donors
    W_P <- stats::setNames(lapply(dn, function(j) { Yu <- if (is.null(um)) Y else um$Ylist[[j]]; pool <- setdiff(dn, j)
      stats::setNames(dsc_fit_qp(Yu[, j], Yu[, pool, drop = FALSE], !post)$w, colnames(Y)[pool]) }), colnames(Y)[dn])
    lto_covariate_balance(Zbal, colnames(Y)[tr], stats::setNames(w, colnames(Y)[-tr]), W_P)
  } else NULL
  # Diagnostics: the same checks as dsc_placebo() (lto_diagnostics()); ridge is recommended
  # only for a treated unit flagged as extreme with trending factors, the case where plain LTO failed
  diag <- NULL
  if (diagnostics) {
    {
      if (isTRUE(verbose)) message("dsc_placebo_lto: diagnostics", if (requireNamespace("fect", quietly = TRUE)) " (fect)" else "")
      s_len <- if (length(time_cols) >= 2) length(unique(periods[[time_cols[length(time_cols)]]])) else NA_integer_
      dp <- lto_diagnostics(Y, tr, post, core, w, Zs, ife_r, ife_seed, r_known = if (!is.null(um) && !is.na(um$r)) um$r else NA,
                            season_length = s_len)
      extreme <- isTRUE(dp$fit_rank$flag) || isTRUE(dp$hull$flag); trending <- isTRUE(dp$persistence$trending)
      rec <- NULL
      txt <- if (extreme && trending && augment != "ridge") {
        rec <- list(augment = "ridge")
        if (augment == "ife") "treated unit appears beyond the donors' range and the factors are trending: re-run with augment = \"ridge\" (with adjust_covariates if needed): IFE augmentation over-rejected badly here in the LTO validation (edge unit, random-walk factors: 30% false positives at 10%; outside: 66%), while IFE adjustment with ridge stayed near nominal (6-11%)" else
        "treated unit appears beyond the donors' range and the factors are trending: re-run with augment = \"ridge\". In the small-pool validation, plain LTO over-rejected here (edge unit, random-walk factors: 16-18% false positives at 10%; 1 SD outside: 20-27%), and ridge augmentation brought this to 6% (outside: 5-12%)"
      } else if (extreme && trending) {
        "treated unit appears beyond the donors' range and the factors are trending: augment = \"ridge\" is in use, as recommended for this case; some bias can remain for a unit far outside the range"
      } else if (extreme && augment == "ife") {
        rec <- list(augment = "none")
        "treated unit appears beyond the donors' range (factors not trending): IFE augmentation over-rejected for extreme treated units in the LTO validation (outside the range, stationary factors: 32-36% false positives at 10%); re-run with augment = \"none\" (with adjust_covariates if a confounder must be adjusted for), which stayed valid here"
      } else if (augment == "ridge") {
        rec <- list(augment = "none")
        paste0(if (extreme) "treated unit appears beyond the donors' range, but the factors are not trending" else "no sign that the treated unit lies beyond the donors",
               ": ridge augmentation is not needed. Plain LTO was valid here and kept more power for ordinary treated units in short real panels (Basque subsamples: cumulative power 38.5% against 28.5%); consider re-running with augment = \"none\"")
      } else if (extreme) {
        "treated unit appears beyond the donors' range (factors not trending): plain LTO stayed valid (conservative) in this situation in simulation"
      } else "no sign that the treated unit lies beyond the donors"
      if (!is.null(balance) && any(balance$flag %in% TRUE))
        txt <- paste0(txt, "; covariate imbalance: ", paste(balance$covariate[balance$flag %in% TRUE], collapse = ", "),
                      " (the treated unit is further from its synthetic control on this covariate than 90% or more of the placebos are from theirs; this biases the estimate if the covariate's effect on the outcome changes over time; see $balance)")
      if (is.na(dp$persistence$trending)) txt <- paste0(txt, " (", dp$hull$note, ": the check for trending factors was not run)")
      if (warn && is.na(dp$persistence$trending) && !is.null(dp$hull$note) && dp$hull$note != "no loadings (0 factors chosen)")
        warning("dsc_placebo_lto diagnostics: ", dp$hull$note, ". The loading-hull and factor-persistence checks (and so the ridge recommendation for an extreme unit with trending factors) were not run; install fect to enable them.", call. = FALSE)
      if (warn && extreme && trending && augment != "ridge") warning("dsc_placebo_lto diagnostics: ", txt, ". See $diagnostics.", call. = FALSE)
      diag <- c(dp[c("fit_rank", "hull", "extrapolation", "persistence", "residual_autocorrelation")], list(assessment = txt, recommendation = rec))
    }
  }
  structure(list(effects = effects, test = test, flagged_periods = periods[post, , drop = FALSE][sim$flagged, , drop = FALSE], weights = sort(stats::setNames(w, colnames(Y)[-tr]), decreasing = TRUE),
       diagnostics = diag,
       balance = balance,
       adjustment = if (is.null(um)) NULL else list(model = if (augment == "ife") "ife_augment" else adjust_model, n_factors = um$r,
                                                    gamma_treated = if (is.null(um$gamma)) NULL else um$gamma[tr, ],
                                                    gamma_range = if (is.null(um$gamma)) NULL else apply(um$gamma, 2, range),
                                                    covariate_adjustment = if (is.null(um$gamma)) NULL else Y[, tr] - um$Ylist[[tr]][, tr]),
       bound = c(lto = lto_bound(N, alpha), approximate_placebo = (floor(N * alpha) + 1) / N,
                 exact_placebo = floor(alpha * N) / N, rejection_threshold = thr),
       settings = list(treated_unit = treated_unit, treatment_period = treatment_period, time_cols = time_cols, outcome = outcome_col,
                       conf_level = conf_level, version = version, balance_covariates = balance_covariates, small_pool_note = small_pool_note, n_units = N, n_donors = N_d, n_matches = core$n_pairs,
                       test_row = test_row, rmspe_floor = rmspe_floor, floor_value = core$floor, span = span, gamma_max = gamma_max,
                       augment = augment, ridge_c = if (is.null(c_r)) NA else c_r, ridge_penalty = if (augment != "ridge") NA else if (is.null(ridge_c)) "cross-validation" else "fixed", ridge_rule = ridge_rule, ridge_windows = ridge_windows,
                       adjust_covariates = adjust_covariates, adjust_model = if (length(adjust_covariates) && augment != "ife") adjust_model else NA),
       core = core), class = "dsc_placebo_lto")
}
