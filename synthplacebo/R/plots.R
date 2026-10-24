#' Plots for dsc_placebo() results
#'
#' `plot_full_path()` draws the whole series: observed and synthetic outcomes, with the pointwise (dark) and
#' simultaneous (light) prediction bands after treatment. `plot_post_effect()` draws the post-treatment effect,
#' per period (`measure = "effect"`) or cumulative, with the simultaneous band shaded and the pointwise interval
#' at the last post-treatment period as an error bar. `plot_prediction()` draws the post-treatment period only:
#' observed and synthetic outcomes with the simultaneous prediction band. `added_effect()` returns the difference
#' between two fits' observed series in each post-treatment period (a known, injected effect), for
#' `plot_post_effect(added = )`.
#'
#' With two or more time columns, the x-axis is labelled at the first sub-period of each value of the first
#' column (e.g. Q1 of each year); with one, at about 12 evenly spaced periods.
#'
#' @param fit Output of [dsc_placebo()].
#' @param title Plot title (`NULL`: automatic).
#' @param unit_label Name of the treated unit in the legend.
#' @param y_label y-axis label (`NULL` for `plot_post_effect()`: automatic).
#' @param measure `"effect"` or `"cumulative"`.
#' @param fill Colour of the band (default: red for the effect, green for the cumulative effect).
#' @param added Optional known effect to draw as a dashed line, one value per post-treatment period (per-period
#'   values; cumulated for `measure = "cumulative"`), typically from `added_effect()`.
#' @param fit_real,fit_constructed Fits on the original data and on the data with an effect added.
#' @return A `ggplot` object (`added_effect()`: a numeric vector).
#' @seealso [plot.dsc_placebo()], which calls these.
#' @name dsc_plots
NULL

## ===========================================================================
## 7. Plotting functions for dsc_symcal_boot() / dsc_placebo() output
## ===========================================================================

##############################################################################
# Generic figures for any output of dsc_symcal_boot() or dsc_placebo(). They use
# only the standard output columns and $settings (time columns, treated unit,
# outcome, confidence level), so they work for any panel and time structure.
#
#   plot_full_path(fit)      whole series: observed vs synthetic, with pointwise
#                            and simultaneous prediction bands after treatment
#   plot_post_effect(fit)    post period: effect (measure = "effect") or
#                            cumulative effect ("cumulative"), simultaneous band
#                            shaded, pointwise interval at the last period
#                            (error bar), optional known added effect (dashed)
#   plot_prediction(fit)     post period: observed vs synthetic, simultaneous
#                            prediction band shaded, pointwise interval at the
#                            last period
#   added_effect(fit_real, fit_constructed)
#                            the post-period difference between two fits'
#                            observed series (e.g. an injected effect), for
#                            plot_post_effect(added = ...)
#
# Common arguments: title (default: automatic), unit_label (default: the treated
# unit), y_label (default: the outcome column). Axis breaks: with two or more
# time columns, the first period of each value of the first column (e.g. Q1 of
# each year); with one time column, about 12 evenly spaced periods.
# Requires ggplot2, dplyr, tidyr and forcats. Value: a ggplot object.
##############################################################################

# Internal: period index, x-axis breaks and labels from the fit's time columns
.dsc_time_axis <- function(e, time_cols) {
  e <- dplyr::mutate(e, x = dplyr::row_number())
  if (length(time_cols) >= 2) {
    first_sub <- min(e[[time_cols[2]]])
    br <- dplyr::filter(e, .data[[time_cols[2]]] == first_sub)
  } else {
    step <- max(1, ceiling(nrow(e) / 12))
    br <- dplyr::slice(e, seq(1, dplyr::n(), by = step))
  }
  list(data = e, breaks = br$x, labels = br[[time_cols[1]]])
}

.dsc_check_fit <- function(fit) {
  for (p in c("ggplot2", "dplyr", "tidyr", "forcats")) if (!requireNamespace(p, quietly = TRUE)) stop("Plotting needs the ", p, " package.")
  need <- c("observed", "synthetic", "effect", "cum_effect", "eff_low_sim", "eff_high_sim", "eff_low_pw", "eff_high_pw",
            "pred_low_sim", "pred_high_sim", "pred_low_pw", "pred_high_pw", "post")
  if (!is.list(fit) || is.null(fit$effects) || !all(need %in% names(fit$effects)) || is.null(fit$settings$time_cols))
    stop("`fit` must be the output of dsc_placebo().")
}

#' @rdname dsc_plots
#' @export
plot_full_path <- function(fit, title = NULL, unit_label = fit$settings$treated_unit, y_label = fit$settings$outcome) {
  .dsc_check_fit(fit)
  level <- paste0(100 * fit$settings$conf_level, "%")
  ax <- .dsc_time_axis(fit$effects, fit$settings$time_cols)
  e <- ax$data
  t0 <- e |> dplyr::filter(post) |> dplyr::slice(1) |> dplyr::pull(x)
  bands <- e |>
    dplyr::filter(post) |>
    dplyr::select(x, pred_low_sim, pred_high_sim, pred_low_pw, pred_high_pw) |>
    tidyr::pivot_longer(-x, names_to = c(".value", "band"), names_pattern = "pred_(low|high)_(sim|pw)") |>
    dplyr::mutate(band = dplyr::if_else(band == "sim", paste(level, "simultaneous band"), paste(level, "pointwise band")),
                  band = forcats::fct_relevel(band, paste(level, "simultaneous band")))
  series <- e |>
    dplyr::select(x, observed, synthetic) |>
    tidyr::pivot_longer(-x, names_to = "series", values_to = "value") |>
    dplyr::mutate(series = dplyr::if_else(series == "observed", paste(unit_label, "(observed)"), paste("Synthetic", unit_label)),
                  series = forcats::fct_relevel(series, paste(unit_label, "(observed)")))
  if (is.null(title)) title <- paste0(unit_label, " vs synthetic ", unit_label, ", with ", level, " prediction bands")
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = bands, ggplot2::aes(x, ymin = low, ymax = high, fill = band), alpha = 0.35) +
    ggplot2::geom_line(data = series, ggplot2::aes(x, value, colour = series, linetype = series), linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = t0 - 0.5, linetype = "dotted", colour = "grey30") +
    ggplot2::scale_fill_manual(values = c("steelblue1", "steelblue4"), name = NULL) +
    ggplot2::scale_colour_manual(values = c("black", "steelblue4"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = ax$breaks, labels = ax$labels, expand = c(0.01, 0)) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = 8),
                   panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::labs(title = title, x = NULL, y = y_label)
}

#' @rdname dsc_plots
#' @export
plot_post_effect <- function(fit, measure = c("effect", "cumulative"), title = NULL, y_label = NULL,
                             fill = NULL, added = NULL, unit_label = fit$settings$treated_unit) {
  .dsc_check_fit(fit)
  measure <- match.arg(measure)
  level <- paste0(100 * fit$settings$conf_level, "%")
  pre <- if (measure == "effect") "eff_" else "cum_eff_"
  if (is.null(fill)) fill <- if (measure == "effect") "firebrick" else "darkgreen"
  if (is.null(y_label)) y_label <- paste(if (measure == "effect") "Effect on" else "Cumulative effect on", fit$settings$outcome)
  if (is.null(title)) title <- paste0(unit_label, ": ", if (measure == "effect") "estimated effect by period" else "cumulative effect")
  ax <- .dsc_time_axis(dplyr::filter(fit$effects, post), fit$settings$time_cols)
  e <- ax$data |>
    dplyr::mutate(estimate = if (measure == "effect") effect else cum_effect,
                  lo_sim = .data[[paste0(pre, "low_sim")]], hi_sim = .data[[paste0(pre, "high_sim")]],
                  lo_pw  = .data[[paste0(pre, "low_pw")]],  hi_pw  = .data[[paste0(pre, "high_pw")]])
  if (!is.null(added)) {
    if (length(added) != nrow(e)) stop("`added` must have one value per post-treatment period (", nrow(e), ").")
    e <- dplyr::mutate(e, added = if (measure == "effect") added else cumsum(added))
  }
  last_period <- dplyr::slice_tail(e, n = 1)
  p <- ggplot2::ggplot(e, ggplot2::aes(x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo_sim, ymax = hi_sim), fill = fill, alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40")
  if (!is.null(added)) p <- p + ggplot2::geom_line(ggplot2::aes(y = added, linetype = "Added effect"), colour = "grey35")
  p +
    ggplot2::geom_line(ggplot2::aes(y = estimate, linetype = "Estimated effect")) +
    ggplot2::geom_point(ggplot2::aes(y = estimate), size = 1.2) +
    ggplot2::geom_errorbar(data = last_period, ggplot2::aes(ymin = lo_pw, ymax = hi_pw), width = 0.35, linewidth = 0.9) +
    ggplot2::geom_point(data = last_period, ggplot2::aes(y = estimate), size = 2.8, shape = 21, fill = "white", stroke = 1.1) +
    ggplot2::scale_linetype_manual(values = c("Estimated effect" = "solid", "Added effect" = "dashed"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = if (is.null(added)) "none" else "bottom",
                   plot.caption = ggplot2::element_text(hjust = 0, size = 9)) +
    ggplot2::labs(title = title, x = NULL, y = y_label,
                  caption = paste0("Shaded: ", level, " simultaneous band (whole post period). Error bar: ", level,
                                   " pointwise interval (last period only)."))
}

#' @rdname dsc_plots
#' @export
plot_prediction <- function(fit, title = NULL, unit_label = fit$settings$treated_unit, y_label = fit$settings$outcome) {
  .dsc_check_fit(fit)
  level <- paste0(100 * fit$settings$conf_level, "%")
  ax <- .dsc_time_axis(dplyr::filter(fit$effects, post), fit$settings$time_cols)
  e <- ax$data
  series <- e |>
    dplyr::select(x, observed, synthetic) |>
    tidyr::pivot_longer(-x, names_to = "series", values_to = "value") |>
    dplyr::mutate(series = dplyr::if_else(series == "observed", paste(unit_label, "(observed)"), paste("Synthetic", unit_label)),
                  series = forcats::fct_relevel(series, paste(unit_label, "(observed)")))
  last_period <- dplyr::slice_tail(e, n = 1)
  if (is.null(title)) title <- paste0(unit_label, ": observed vs synthetic (prediction intervals)")
  ggplot2::ggplot(e, ggplot2::aes(x)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pred_low_sim, ymax = pred_high_sim), fill = "steelblue", alpha = 0.2) +
    ggplot2::geom_line(data = series, ggplot2::aes(x, value, colour = series), linewidth = 0.8) +
    ggplot2::geom_errorbar(data = last_period, ggplot2::aes(ymin = pred_low_pw, ymax = pred_high_pw), width = 0.35,
                           linewidth = 0.9, colour = "steelblue4") +
    ggplot2::geom_point(data = last_period, ggplot2::aes(y = synthetic), size = 2.8, shape = 21, fill = "white",
                        colour = "steelblue4", stroke = 1.1) +
    ggplot2::scale_colour_manual(values = c("black", "steelblue4"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = ax$breaks, labels = ax$labels) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0, size = 9)) +
    ggplot2::labs(title = title, x = NULL, y = y_label,
                  caption = paste0("Shaded: ", level, " simultaneous prediction band. Error bar: ", level,
                                   " pointwise prediction interval (last period only)."))
}

#' @rdname dsc_plots
#' @export
added_effect <- function(fit_real, fit_constructed) {
  .dsc_check_fit(fit_real); .dsc_check_fit(fit_constructed)
  if (!identical(fit_real$effects$post, fit_constructed$effects$post)) stop("The two fits must cover the same periods.")
  (fit_constructed$effects$observed - fit_real$effects$observed)[fit_real$effects$post]
}


