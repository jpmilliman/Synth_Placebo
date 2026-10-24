# covariate_trajectory_check.R
#
# covariate_trajectory_check(): does a covariate's path for the treated unit break away after treatment from
# what its donors predict, by more than the donors' own covariate paths break away from each other?
#
# Each covariate is treated as an outcome: the demeaned synthetic control is fitted to the treated unit's
# pre-treatment covariate path, and its post-treatment gap is ranked against the gaps of every donor fitted on
# the other donors (dsc_placebo()'s placebo test: RMSPE-standardized gaps, single period, cumulative and
# simultaneous). It is a screen for the two cases in which time-varying covariate adjustment matters:
#   - the covariate moved because of the treatment (a mediator or another outcome): adjusting for it would
#     remove part of the effect;
#   - something else moved it for the treated unit at about the same time (a coincident confounder, such as
#     fentanyl in Oregon): not adjusting for it would bias the estimate.
# The test cannot tell these apart; that takes knowledge of the case. If a covariate shows no divergence,
# adjusting for it changes the estimate little either way: with linear adjustment the adjusted gap is about the
# unadjusted gap minus gamma x (the covariate's gap), which is reported as implied_change when gamma is given.
#
# NOT A SUBSTITUTE FOR JUDGEMENT. A non-significant result does not show that a covariate is unaffected by the
# treatment or that adjustment is unnecessary: this test has the same limited power as the outcome test, so a
# covariate can move by several percent without the test noticing. Plot the covariates' trajectories
# (plot_covariate_trajectories()) and decide from what is known about the case.
#
# Needs Experimental_inference/LTO/functions/dsc_placebo_lto.R (source it first; it uses dsc_fit_qp(), lto_core(),
# lto_pvalues() and standard_simultaneous()) and quadprog; ggplot2 for the plot. No fect: the run takes well
# under a second per covariate.
# Written by Claude (Anthropic), at the direction of Jeffrey Milliman.

#' Screen covariates for a post-treatment break from their synthetic controls
#'
#' Does a covariate's path for the treated unit break away after treatment from what its donors predict, by
#' more than the donors' own covariate paths break away from each other? Each covariate is treated as an
#' outcome: the demeaned synthetic control is fitted to the treated unit's pre-treatment covariate path, and its
#' post-treatment gap is ranked against every donor's placebo gap (single period, cumulative and simultaneous,
#' as in [dsc_placebo()]).
#'
#' @details
#' Time-varying covariate adjustment matters in two cases that call for opposite decisions: the treatment moved
#' the covariate (a mediator or another outcome: adjusting removes part of the effect), or something else moved
#' it for the treated unit at about the same time (a coincident confounder: not adjusting biases the estimate).
#' In both the covariate breaks from its synthetic control. If it does not, adjustment changes the estimate
#' little either way: the adjusted gap is about the unadjusted gap minus \eqn{\gamma} times the covariate's gap,
#' reported as `implied_change` when `gamma` is given.
#'
#' **This is a screen, not a decision rule.** A covariate that is not flagged has not been shown to be
#' unaffected by the treatment: the test has the same limited power as the outcome test. A flagged covariate may
#' be a mediator or a confounder; the test cannot tell which. Plot the trajectories with
#' [plot_covariate_trajectories()] and decide from what is known about the case.
#'
#' Takes well under a second per covariate; no `fect`.
#'
#' @inheritParams dsc_placebo
#' @param covariates Numeric covariate columns to check. Units with missing values for a covariate are left out
#'   of its check; a covariate missing for the treated unit is skipped.
#' @param test_period Period of the single-period and cumulative tests (default: the last post period).
#' @param alpha Level at which a covariate is flagged (by any of the three tests).
#' @param gamma Optional named vector of the covariates' effects on the outcome (e.g.
#'   `fit$adjustment$gamma_treated` from an adjusted [dsc_placebo_lto()] fit), for the implied change in the
#'   estimate.
#' @param rmspe_floor,span As in [dsc_placebo()].
#' @return An object of class `"covariate_trajectory_check"`: a list with `summary` (one row per covariate:
#'   gaps, `ratio_post_pre`, `p_final`, `p_cumulative`, `p_simultaneous`, `smallest_p`, `flag` and, with
#'   `gamma`, `implied_change_mean` and `implied_change_cumulative`), `paths` (every unit's covariate, synthetic
#'   control and standardized gap by period, for plotting), `skipped`, `note` and `settings`.
#' @seealso [plot_covariate_trajectories()], [dsc_placebo()].
#' @examples
#' \donttest{
#' if (requireNamespace("augsynth", quietly = TRUE)) {
#'   kansas <- as.data.frame(augsynth::kansas)
#'   ck <- covariate_trajectory_check(kansas, "state", c("year", "qtr"),
#'                                    c("estabscapita", "emplvlcapita", "avg_wkly_wage"),
#'                                    "Kansas", c(year = 2012, qtr = 2))
#'   ck
#'   plot(ck)
#' }
#' }
#' @export
covariate_trajectory_check <- function(data, unit_col, time_cols, covariates, treated_unit, treatment_period,
                                       test_period = NULL, alpha = 0.10, gamma = NULL, rmspe_floor = 0.25, span = 0.75) {
  d <- as.data.frame(data)
  need <- c(unit_col, time_cols, covariates)
  if (!all(need %in% names(d))) stop("columns not found: ", paste(setdiff(need, names(d)), collapse = ", "))
  if (!is.null(gamma) && (is.null(names(gamma)) || !all(names(gamma) %in% covariates)))
    stop("gamma must be a named numeric vector with names among the covariates (e.g. fit$adjustment$gamma_treated).")
  d <- d[do.call(order, d[time_cols]), , drop = FALSE]
  key <- do.call(paste, c(d[time_cols], sep = "_"))
  periods <- unique(d[time_cols]); rownames(periods) <- NULL
  wide <- function(v) { w <- stats::reshape(data.frame(key = key, unit = as.character(d[[unit_col]]), z = d[[v]]), idvar = "key", timevar = "unit", direction = "wide")
    m <- as.matrix(w[match(unique(key), w$key), -1, drop = FALSE]); colnames(m) <- sub("^z\\.", "", colnames(m)); m }
  after <- function(tp) vapply(seq_len(nrow(periods)), function(r) {
    for (v in names(tp)) if (periods[[v]][r] != tp[[v]]) return(periods[[v]][r] > tp[[v]]); TRUE }, logical(1))
  post <- after(treatment_period); T1 <- sum(post)
  if (T1 < 1 || sum(!post) < 2) stop("treatment_period leaves no post or too few pre periods.")
  h <- if (is.null(test_period)) T1 else {
    hit <- which(Reduce(`&`, lapply(names(test_period), function(v) periods[[v]] == test_period[[v]])))
    if (length(hit) != 1 || !post[hit]) stop("test_period must be one post-treatment period."); hit - sum(!post) }
  plabel <- do.call(paste, c(periods, sep = "-"))

  rows <- list(); paths <- list(); skipped <- character(0)
  for (v in covariates) {
    Z <- wide(v)
    keep <- colSums(is.na(Z)) == 0
    if (!treated_unit %in% colnames(Z) || !keep[treated_unit]) { skipped <- c(skipped, paste0(v, " (missing for the treated unit)")); next }
    if (sum(keep) < 5) { skipped <- c(skipped, paste0(v, " (fewer than 4 donors without missing values)")); next }
    Z <- Z[, keep, drop = FALSE]; tr <- match(treated_unit, colnames(Z))
    core <- lto_core(Z, tr, post, rmspe_floor = rmspe_floor, tournament = FALSE)
    st <- core$standard; pv <- lto_pvalues(core, 0)
    sim <- suppressWarnings(standard_simultaneous(core, 0.95, span))   # loess warns on flat, near-zero drift profiles
    g <- st$g; gc <- st$gc
    rows[[v]] <- data.frame(
      covariate = v, donors = ncol(Z) - 1, pre_rmspe = st$s,
      mean_post_gap = mean(g), final_gap = g[h], cumulative_gap = gc[h],
      ratio_post_pre = sqrt(mean(g^2)) / st$s,
      p_final = pv$p_placebo[h], p_cumulative = pv$p_placebo_cum[h], p_simultaneous = sim$p_sim,
      smallest_p = 1 / ncol(Z),
      flag = min(pv$p_placebo[h], pv$p_placebo_cum[h], sim$p_sim) <= alpha + 1e-12,
      implied_change_mean = if (!is.null(gamma) && v %in% names(gamma)) -gamma[[v]] * mean(g) else NA_real_,
      implied_change_cumulative = if (!is.null(gamma) && v %in% names(gamma)) -gamma[[v]] * gc[h] else NA_real_)
    # paths: the treated unit and its synthetic control, and every unit's standardized gap (treated and placebos)
    gT <- core$gap_treated; sdon <- st$s_donors
    paths[[v]] <- rbind(
      data.frame(covariate = v, period = plabel, t = seq_along(plabel), post = post, unit = treated_unit, role = "treated",
                 value = Z[, tr], synthetic = core$synthetic_treated, gap = gT, std_gap = gT / st$s),
      do.call(rbind, lapply(seq_along(core$donors), function(j) {
        gj <- st$gaps[, j]
        data.frame(covariate = v, period = plabel, t = seq_along(plabel), post = post, unit = colnames(Z)[core$donors[j]], role = "donor",
                   value = Z[, core$donors[j]], synthetic = Z[, core$donors[j]] - gj, gap = gj, std_gap = gj / sdon[j]) })))
  }
  summary <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(summary)) rownames(summary) <- NULL
  note <- paste("A covariate that is not flagged has not been shown to be unaffected by the treatment: the test has the same",
                "limited power as the outcome test. A flagged covariate may be a mediator (do not adjust) or a coincident",
                "confounder (adjust): the test cannot tell which. Plot the trajectories and use knowledge of the case.")
  structure(list(summary = summary, paths = if (length(paths)) do.call(rbind, paths) else NULL, skipped = skipped, note = note,
                 settings = list(treated_unit = treated_unit, treatment_period = treatment_period, test_period_row = h, alpha = alpha,
                                 n_post = T1, rmspe_floor = rmspe_floor, span = span, gamma = gamma)),
            class = "covariate_trajectory_check")
}

#' @rdname covariate_trajectory_check
#' @export
print.covariate_trajectory_check <- function(x, digits = 3, ...) {
  cat("Covariate trajectory check: ", x$settings$treated_unit, ", ", x$settings$n_post, " post-treatment periods, alpha = ", x$settings$alpha, "\n\n", sep = "")
  if (!is.null(x$summary)) {
    s <- x$summary[, intersect(c("covariate", "mean_post_gap", "ratio_post_pre", "p_final", "p_cumulative",
                                 "p_simultaneous", "flag", "implied_change_mean"), names(x$summary))]
    if (all(is.na(s$implied_change_mean))) s$implied_change_mean <- NULL
    num <- vapply(s, is.numeric, logical(1)); s[num] <- lapply(s[num], signif, digits)
    print(s, row.names = FALSE)
  }
  if (length(x$skipped)) cat("\nSkipped:", paste(x$skipped, collapse = "; "), "\n")
  cat("\n", strwrap(paste("Note:", x$note), width = 100), sep = "\n")
  invisible(x)
}

# Trajectories for a visual check: top, the treated unit's covariate and its synthetic control; bottom, the treated
# unit's RMSPE-standardized gap (coloured) against every donor's placebo gap (grey). One column per covariate.
#' Plot covariate trajectories from covariate_trajectory_check()
#'
#' For each covariate (one column): the treated unit's covariate and its synthetic control (top), and the
#' treated unit's RMSPE-standardized gap against every donor's (bottom, grey). Each panel has its own scale.
#'
#' @param x Output of [covariate_trajectory_check()].
#' @param covariates Covariates to draw (default: all).
#' @param title Plot title.
#' @return A `ggplot` object.
#' @export
plot_covariate_trajectories <- function(x, covariates = NULL, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("plot_covariate_trajectories() needs ggplot2.")
  p <- x$paths; if (!is.null(covariates)) p <- p[p$covariate %in% covariates, ]
  t0 <- min(p$t[p$post]) - 0.5; brks <- unique(round(seq(1, max(p$t), length.out = 6)))
  tr <- p[p$role == "treated", ]
  lev <- rbind(data.frame(covariate = tr$covariate, t = tr$t, value = tr$value, series = x$settings$treated_unit),
               data.frame(covariate = tr$covariate, t = tr$t, value = tr$synthetic, series = paste("Synthetic", x$settings$treated_unit)))
  pl <- c("Covariate", "Standardized gap")
  lev$panel <- factor(pl[1], pl); gaps <- p; gaps$panel <- factor(pl[2], pl)
  ggplot2::ggplot() +
    ggplot2::geom_line(data = transform(gaps[gaps$role == "donor", ], value = std_gap), ggplot2::aes(t, value, group = unit), colour = "grey75", linewidth = 0.3) +
    ggplot2::geom_line(data = transform(gaps[gaps$role == "treated", ], value = std_gap), ggplot2::aes(t, value), colour = "#c8553d", linewidth = 0.9) +
    ggplot2::geom_line(data = lev, ggplot2::aes(t, value, linetype = series), colour = "#c8553d", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = t0, linetype = "dashed", colour = "grey40") +
    ggplot2::facet_wrap(ggplot2::vars(panel, covariate), scales = "free_y", ncol = length(unique(p$covariate)),
                        labeller = ggplot2::label_wrap_gen(width = 40, multi_line = FALSE)) +   # each panel on its own scale
    ggplot2::scale_x_continuous(breaks = brks, labels = tr$period[match(brks, tr$t)]) +
    ggplot2::labs(x = NULL, y = NULL, linetype = NULL,
                  title = if (is.null(title)) paste("Covariate trajectories:", x$settings$treated_unit, "and its donors") else title,
                  subtitle = "Top: covariate and synthetic control. Bottom: standardized gaps (grey: donors). Dashed: treatment.",
                  caption = "Not flagging a covariate does not show it is unaffected by the treatment.") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
