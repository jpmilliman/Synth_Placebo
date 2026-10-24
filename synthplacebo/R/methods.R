# print(), summary() and plot() methods. Printing shows only what is needed to read a result; the full
# results stay in the list (names(x)).

.sp_period_label <- function(tp) if (is.null(tp)) "last post period" else paste(paste0(names(tp), " ", tp), collapse = ", ")
.sp_fmt <- function(df, digits = 3) { num <- vapply(df, is.numeric, logical(1)); df[num] <- lapply(df[num], signif, digits); df }
.sp_rule <- function(width = 78) cat(strrep("-", width), "\n", sep = "")

# ---- dsc_placebo -------------------------------------------------------------------------------------------

#' Print, summarize and plot dsc_placebo() results
#'
#' `print()` shows the setting, the tests, the largest weights and the diagnostics' assessment. `summary()` adds
#' the post-treatment effects with their intervals, the balance check, the covariate adjustment and the
#' diagnostics' flags. `plot()` draws one of the standard figures. The full results remain in the object
#' (`names(x)`).
#'
#' @param x,object Output of [dsc_placebo()].
#' @param n_weights Number of donor weights to show.
#' @param digits Significant digits.
#' @param type `"full_path"` ([plot_full_path()]), `"effect"` or `"cumulative"` ([plot_post_effect()]), or
#'   `"prediction"` ([plot_prediction()]).
#' @param ... Passed to the plotting function (e.g. `title`, `y_label`).
#' @return `print()` returns `x` invisibly; `summary()` an object of class `"summary.dsc_placebo"`; `plot()` a
#'   `ggplot` object.
#' @name dsc_placebo_methods
NULL

#' @rdname dsc_placebo_methods
#' @export
print.dsc_placebo <- function(x, n_weights = 5, digits = 3, ...) {
  s <- x$settings; e <- x$effects
  cat("Demeaned synthetic control with placebo inference\n"); .sp_rule()
  cat("Treated unit:", s$treated_unit, "| donors:", s$n_donors, "| pre/post periods:", sum(!e$post), "/", sum(e$post), "\n")
  cat("Treatment from:", .sp_period_label(s$treatment_period), "| test period:", if (is.null(s$test_period)) "last post period" else paste(s$test_period, collapse = " "),
      "| level:", paste0(100 * s$conf_level, "%"), "| scale:", s$scale, "\n")
  extras <- c(if (length(s$adjust_covariates)) paste0("adjusted for ", paste(s$adjust_covariates, collapse = ", "), " (", s$adjust_model, ")"),
              if (!identical(s$augment, "none")) paste0("augment = \"", s$augment, "\""), if (!identical(s$estimator, "dsc")) paste0("estimator = \"", s$estimator, "\""))
  if (length(extras)) cat(paste(extras, collapse = "; "), "\n")
  if (isTRUE(s$rank_used == 0)) cat("Note: too few donors for this level; intervals are infinite (use dsc_placebo_lto() or a lower conf_level).\n")
  cat("\nTests\n"); print(.sp_fmt(x$test[, intersect(c("test", "measure", "period", "estimate", "lower", "upper", "p_value", "significant"), names(x$test))], digits), row.names = FALSE)
  if (length(x$weights)) { w <- utils::head(x$weights, n_weights)
    cat("\nLargest weights:", paste0(names(w), " ", sprintf("%.3f", w), collapse = ", "), if (length(x$weights) > n_weights) paste0(" (+", length(x$weights) - n_weights, " more)") else "", "\n") }
  if (!is.null(x$diagnostics$assessment)) cat("\n", strwrap(paste("Diagnostics:", x$diagnostics$assessment), width = 100), sep = "\n")
  if (!is.null(x$balance) && any(x$balance$flag %in% TRUE)) cat("Balance: flagged ", paste(x$balance$covariate[x$balance$flag %in% TRUE], collapse = ", "), "\n", sep = "")
  cat("\nsummary() for the effects by period and the checks; plot(x, type = ...) for figures.\n")
  invisible(x)
}

#' @rdname dsc_placebo_methods
#' @export
summary.dsc_placebo <- function(object, ...) {
  e <- object$effects; tc <- object$settings$time_cols
  cols <- intersect(c(tc, "observed", "synthetic", "effect", "eff_low_pw", "eff_high_pw", "eff_low_sim", "eff_high_sim", "p_pointwise", "cum_effect", "cum_p_pointwise"), names(e))
  d <- object$diagnostics
  flags <- if (is.null(d)) NULL else data.frame(
    check = c("fit rank", "loading hull", "trending factors", "residual autocorrelation"),
    result = c(if (is.null(d$fit_rank$flag)) NA else if (isTRUE(d$fit_rank$flag)) "flagged" else "ok",
               if (is.null(d$hull$flag) || is.na(d$hull$flag)) "not run" else if (isTRUE(d$hull$flag)) "flagged" else "ok",
               if (is.null(d$persistence$trending) || is.na(d$persistence$trending)) "not run" else if (isTRUE(d$persistence$trending)) "yes" else "no",
               if (is.null(d$residual_autocorrelation$level)) NA else paste0(d$residual_autocorrelation$level, " (rho = ", signif(d$residual_autocorrelation$rho, 2), ")")))
  structure(list(settings = object$settings, test = object$test, effects = e[e$post, cols, drop = FALSE], weights = object$weights,
                 balance = object$balance, adjustment = object$adjustment[intersect(c("model", "gamma", "n_factors", "unadjusted"), names(object$adjustment))],
                 diagnostics = flags, assessment = d$assessment, recommendation = d$recommendation), class = "summary.dsc_placebo")
}

#' @rdname dsc_placebo_methods
#' @export
print.summary.dsc_placebo <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("Demeaned synthetic control: ", s$treated_unit, ", ", s$n_donors, " donors, ", 100 * s$conf_level, "% level, scale \"", s$scale, "\"\n", sep = ""); .sp_rule()
  cat("Tests\n"); print(.sp_fmt(x$test[, intersect(c("test", "measure", "period", "estimate", "lower", "upper", "p_value", "significant"), names(x$test))], digits), row.names = FALSE)
  cat("\nPost-treatment effects (_pw: pointwise, valid for one period chosen in advance; _sim: simultaneous band)\n")
  print(.sp_fmt(x$effects, digits), row.names = FALSE)
  if (length(x$weights)) cat("\nWeights:", paste0(names(x$weights), " ", sprintf("%.3f", x$weights), collapse = ", "), "\n")
  if (length(x$adjustment)) { cat("\nCovariate adjustment (", x$adjustment$model, "): gamma = ", paste0(names(x$adjustment$gamma), " ", signif(x$adjustment$gamma, digits), collapse = ", "), "\n", sep = "")
    if (!is.null(x$adjustment$unadjusted)) cat("Unadjusted estimates at the test period: effect ", signif(x$adjustment$unadjusted[["effect"]], digits),
                                               ", cumulative ", signif(x$adjustment$unadjusted[["cumulative"]], digits), "\n", sep = "") }
  if (!is.null(x$balance)) { cat("\nCovariate balance (pre-treatment means)\n"); print(.sp_fmt(x$balance, digits), row.names = FALSE) }
  if (!is.null(x$diagnostics)) { cat("\nDiagnostics\n"); print(x$diagnostics, row.names = FALSE)
    if (!is.null(x$assessment)) cat("\n", strwrap(x$assessment, width = 100), sep = "\n")
    if (length(x$recommendation)) cat("Recommended: re-run with ", paste0(names(x$recommendation), " = ", vapply(x$recommendation, deparse, ""), collapse = ", "), "\n", sep = "") }
  invisible(x)
}

#' @rdname dsc_placebo_methods
#' @export
plot.dsc_placebo <- function(x, type = c("full_path", "effect", "cumulative", "prediction"), ...) {
  type <- match.arg(type)
  switch(type, full_path = plot_full_path(x, ...), effect = plot_post_effect(x, measure = "effect", ...),
         cumulative = plot_post_effect(x, measure = "cumulative", ...), prediction = plot_prediction(x, ...))
}

# ---- dsc_placebo_lto ---------------------------------------------------------------------------------------

#' Print, summarize and plot dsc_placebo_lto() results
#'
#' `print()` shows the setting, the tests (LTO p-value, rejection, the standard placebo p-value and the
#' sensitivity \eqn{\Gamma}), the Type-I bound and any small-pool note. `summary()` adds the post-treatment
#' effects with their intervals, the balance check and the diagnostics. `plot()` draws the effect, the
#' cumulative effect or the whole series with the LTO intervals and simultaneous band.
#'
#' @param x,object Output of [dsc_placebo_lto()].
#' @param digits Significant digits.
#' @param type `"effect"`, `"cumulative"` or `"full_path"`.
#' @param title Plot title.
#' @param ... Unused.
#' @return `print()` returns `x` invisibly; `summary()` an object of class `"summary.dsc_placebo_lto"`;
#'   `plot()` a `ggplot` object.
#' @name dsc_placebo_lto_methods
NULL

#' @rdname dsc_placebo_lto_methods
#' @export
print.dsc_placebo_lto <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("Leave-two-out placebo inference (", s$version, ")\n", sep = ""); .sp_rule()
  cat("Treated unit:", s$treated_unit, "| units:", s$n_units, "| matches:", s$n_matches, "| level:", paste0(100 * s$conf_level, "%"), "\n")
  cat("Treatment from:", .sp_period_label(s$treatment_period),
      if (length(s$adjust_covariates)) paste0("| adjusted for ", paste(s$adjust_covariates, collapse = ", "), " (", s$adjust_model, ")") else "",
      if (!identical(s$augment, "none")) paste0("| augment = \"", s$augment, "\"") else "", "\n")
  cat("\nTests\n"); print(.sp_fmt(x$test, digits), row.names = FALSE)
  cat("\nType-I bound at this N and level: ", signif(x$bound[["lto"]], digits), " (rejection threshold ", signif(x$bound[["rejection_threshold"]], digits), ")\n", sep = "")
  if (length(s$small_pool_note)) cat(strwrap(paste("Small pool:", paste(s$small_pool_note, collapse = "; ")), width = 100), sep = "\n")
  if (!is.null(x$diagnostics$assessment)) cat("\n", strwrap(paste("Diagnostics:", x$diagnostics$assessment), width = 100), sep = "\n")
  if (!is.null(x$balance) && any(x$balance$flag %in% TRUE)) cat("Balance: flagged ", paste(x$balance$covariate[x$balance$flag %in% TRUE], collapse = ", "), "\n", sep = "")
  cat("\nsummary() for the effects by period; plot(x, type = ...) for figures.\n")
  invisible(x)
}

#' @rdname dsc_placebo_lto_methods
#' @export
summary.dsc_placebo_lto <- function(object, ...) {
  e <- object$effects; tc <- object$settings$time_cols
  cols <- intersect(c(tc, "observed", "synthetic", "effect", "eff_low", "eff_high", "eff_low_sim", "eff_high_sim", "p_lto", "cum_effect", "cum_eff_low", "cum_eff_high", "p_lto_cum"), names(e))
  structure(list(settings = object$settings, test = object$test, effects = e[e$post, cols, drop = FALSE], weights = object$weights,
                 bound = object$bound, balance = object$balance, assessment = object$diagnostics$assessment,
                 recommendation = object$diagnostics$recommendation), class = "summary.dsc_placebo_lto")
}

#' @rdname dsc_placebo_lto_methods
#' @export
print.summary.dsc_placebo_lto <- function(x, digits = 3, ...) {
  s <- x$settings
  cat("Leave-two-out placebo inference: ", s$treated_unit, ", ", s$n_units, " units, ", 100 * s$conf_level, "% level\n", sep = ""); .sp_rule()
  cat("Tests\n"); print(.sp_fmt(x$test, digits), row.names = FALSE)
  cat("\nPost-treatment effects (eff_: pointwise LTO interval; _sim: simultaneous band)\n"); print(.sp_fmt(x$effects, digits), row.names = FALSE)
  if (length(x$weights)) cat("\nWeights:", paste0(names(x$weights[x$weights > 1e-6]), " ", sprintf("%.3f", x$weights[x$weights > 1e-6]), collapse = ", "), "\n")
  cat("\nBounds:", paste0(names(x$bound), " ", signif(x$bound, digits), collapse = ", "), "\n")
  if (!is.null(x$balance)) { cat("\nCovariate balance (pre-treatment means)\n"); print(.sp_fmt(x$balance, digits), row.names = FALSE) }
  if (!is.null(x$assessment)) cat("\n", strwrap(paste("Diagnostics:", x$assessment), width = 100), sep = "\n")
  invisible(x)
}

#' @rdname dsc_placebo_lto_methods
#' @export
plot.dsc_placebo_lto <- function(x, type = c("effect", "cumulative", "full_path"), title = NULL, ...) {
  type <- match.arg(type)
  e <- x$effects; e$t <- seq_len(nrow(e)); s <- x$settings
  lab <- do.call(paste, c(e[s$time_cols], sep = "-"))
  brks <- unique(round(seq(1, nrow(e), length.out = min(nrow(e), 10))))
  t0 <- min(e$t[e$post]) - 0.5; lvl <- paste0(100 * s$conf_level, "%")
  base <- ggplot2::ggplot() + ggplot2::geom_vline(xintercept = t0, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_x_continuous(breaks = brks, labels = lab[brks]) + ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
  p <- e[e$post, ]
  if (type == "full_path") {
    d <- rbind(data.frame(t = e$t, value = e$observed, series = s$treated_unit), data.frame(t = e$t, value = e$synthetic, series = paste("Synthetic", s$treated_unit)))
    base + ggplot2::geom_ribbon(data = p, ggplot2::aes(.data$t, ymin = .data$pred_low_sim, ymax = .data$pred_high_sim), fill = "grey85") +
      ggplot2::geom_line(data = d, ggplot2::aes(.data$t, .data$value, linetype = .data$series)) +
      ggplot2::labs(x = NULL, y = s$outcome, linetype = NULL,
                    title = if (is.null(title)) paste0(s$treated_unit, " and its synthetic control") else title,
                    subtitle = paste0("Shaded: ", lvl, " LTO simultaneous prediction band"))
  } else {
    cum <- type == "cumulative"
    est <- if (cum) p$cum_effect else p$effect
    lo <- if (cum) p$cum_eff_low_sim else p$eff_low_sim; hi <- if (cum) p$cum_eff_high_sim else p$eff_high_sim
    pw <- p[s$test_row, ]; pwlo <- if (cum) pw$cum_eff_low else pw$eff_low; pwhi <- if (cum) pw$cum_eff_high else pw$eff_high
    d <- data.frame(t = p$t, est = est, lo = lo, hi = hi)
    base + ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
      ggplot2::geom_ribbon(data = d, ggplot2::aes(.data$t, ymin = .data$lo, ymax = .data$hi), fill = if (cum) "#b7d7c4" else "#f2c4b8") +
      ggplot2::geom_line(data = d, ggplot2::aes(.data$t, .data$est), linewidth = 0.8) +
      ggplot2::geom_errorbar(data = data.frame(t = pw$t, lo = pwlo, hi = pwhi), ggplot2::aes(.data$t, ymin = .data$lo, ymax = .data$hi), width = 0.3) +
      ggplot2::labs(x = NULL, y = if (cum) paste("Cumulative effect,", s$outcome) else paste("Effect,", s$outcome),
                    title = if (is.null(title)) paste0(s$treated_unit, ": ", if (cum) "cumulative effect" else "effect by period") else title,
                    subtitle = paste0("Shaded: ", lvl, " LTO simultaneous band; error bar: pointwise interval at the test period"))
  }
}

# ---- synth_power -------------------------------------------------------------------------------------------

#' Print, summarize and plot synth_power() results
#'
#' `print()` shows the false-positive rates and the minimum detectable effects at the first level in the
#' object, one row per shape, test and direction, one column per method and sidedness. `summary()` returns that
#' table for every level. `plot()` calls [plot_synth_power()].
#'
#' @param x,object Output of [synth_power()].
#' @param alpha Level to show (default: the first in the object).
#' @param ... Passed to [plot_synth_power()] (`plot()`); unused otherwise.
#' @return `print()` returns `x` invisibly; `summary()` a data frame; `plot()` a `ggplot` object.
#' @name synth_power_methods
NULL

.sp_mde_table <- function(x, alpha) {
  m <- x$mde[abs(x$mde$alpha - alpha) < 1e-12, ]
  m$value <- ifelse(is.na(m$mde_interpolated), paste0("nr (", round(100 * m$max_power_tested), "%)"),
                    if (x$settings$effect_scale == "percent") paste0(round(100 * m$mde_interpolated, 1), "%") else as.character(signif(m$mde_interpolated, 3)))
  m$key <- paste(m$method, ifelse(m$sided == "two.sided", "two-sided", "one-sided"), sep = ", ")
  w <- stats::reshape(m[, c("shape", "test", "direction", "key", "value")], idvar = c("shape", "test", "direction"), timevar = "key", direction = "wide")
  names(w) <- sub("^value\\.", "", names(w)); rownames(w) <- NULL
  w[order(match(w$shape, c("constant", "ramp", "plateau", "decay", "shock")), match(w$test, c("pointwise", "cumulative", "simultaneous")), w$direction), ]
}

#' @rdname synth_power_methods
#' @export
print.synth_power <- function(x, alpha = x$settings$alpha[1], ...) {
  s <- x$settings
  cat("Power of the placebo and LTO tests\n"); .sp_rule()
  cat(s$n_units, " pseudo-treated units (", s$n_donors, " donors each)", if (!is.null(s$dropped)) paste0(", ", s$dropped, " dropped"), " | ",
      s$n_post, " post periods | effects: ", s$effect_scale, " | level: ", alpha, " | target power: ", s$target_power, "\n", sep = "")
  z <- x$power[x$power$effect == 0 & abs(x$power$alpha - alpha) < 1e-12, ]
  fp <- tapply(z$power, paste(z$method, z$test), function(v) paste(sprintf("%.1f%%", 100 * sort(unique(v))), collapse = "-"))
  cat("\nFalse-positive rate (no effect): ", paste(names(fp), fp, collapse = "; "), "\n", sep = "")
  cat("Smallest possible placebo p-value: ", signif(s$smallest_placebo_p, 3), "\n", sep = "")
  cat("\nMinimum detectable effect (", 100 * s$target_power, "% power; nr: not reached, with the highest power tested)\n", sep = "")
  print(.sp_mde_table(x, alpha), row.names = FALSE)
  invisible(x)
}

#' @rdname synth_power_methods
#' @export
summary.synth_power <- function(object, ...) do.call(rbind, lapply(object$settings$alpha, function(a) cbind(alpha = a, .sp_mde_table(object, a))))

#' @rdname synth_power_methods
#' @export
plot.synth_power <- function(x, ...) plot_synth_power(x, ...)

# ---- covariate_trajectory_check ----------------------------------------------------------------------------

#' @rdname covariate_trajectory_check
#' @param x Output of `covariate_trajectory_check()` (methods).
#' @param digits Significant digits.
#' @param ... Passed to [plot_covariate_trajectories()] (`plot()`); unused otherwise.
#' @export
plot.covariate_trajectory_check <- function(x, ...) plot_covariate_trajectories(x, ...)
