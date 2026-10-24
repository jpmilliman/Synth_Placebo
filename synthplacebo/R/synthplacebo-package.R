#' synthplacebo: placebo-calibrated inference for the demeaned synthetic control
#'
#' Synthetic control for a single treated unit, with inference calibrated against placebo fits of the donors.
#'
#' @section Main functions:
#' * [dsc_placebo()]: the demeaned synthetic control with placebo-calibrated pointwise intervals, simultaneous
#'   bands and prediction intervals; covariate adjustment, ridge augmentation, a balance check and diagnostics.
#'   Plots: [plot_full_path()], [plot_post_effect()], [plot_prediction()], or `plot(fit)`.
#' * [dsc_placebo_lto()]: the leave-two-out placebo test (Lei and Sudijono 2025) for small donor pools, with a
#'   simultaneous band and a sensitivity analysis for unequal assignment probabilities.
#' * [synth_power()]: power and minimum detectable effects of both tests on a given panel, for several effect
#'   shapes, two- or one-sided ([plot_synth_power()]).
#' * [covariate_trajectory_check()]: a screen for covariates whose post-treatment path breaks from their
#'   synthetic control (possible mediators or coincident confounders), with [plot_covariate_trajectories()].
#'
#' Every result has `print()`, `summary()` and `plot()` methods; the full results stay in the returned list.
#'
#' @section Status:
#' Written by Claude (Anthropic), at the direction of Jeffrey Milliman. The methods are validated by the
#' simulations in the repository (<https://github.com/jpmilliman/Synth_Placebo>), not by independent peer
#' review. [dsc_placebo_lto()], [synth_power()] and [covariate_trajectory_check()] are experimental.
#' The leave-two-out test follows Lei and Sudijono (2025), a working paper; its adaptation to the demeaned
#' synthetic control, the simultaneous band, the sensitivity analysis and the ridge option are extensions
#' written for this package.
#'
#' @keywords internal
#' @importFrom rlang .data
#' @importFrom stats as.formula median model.matrix predict setNames
#' @importFrom utils head tail
"_PACKAGE"

# Column names used inside ggplot2::aes() and dplyr verbs (non-standard evaluation).
utils::globalVariables(c("band", "covariate", "cum_effect", "effect", "estimate", "hi_pw", "hi_sim", "high", "line", "lo_pw",
                         "lo_sim", "low", "method", "observed", "panel", "post", "power", "pred_high_pw", "pred_high_sim",
                         "pred_low_pw", "pred_low_sim", "series", "shape", "sided", "std_gap", "synthetic", "unit", "value", "x"))
