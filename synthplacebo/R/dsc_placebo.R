# dsc_placebo(): the main method, and its internal helpers. Moved from
# synthetic_control_functions/synthetic_control_functions_claude.R (frozen copy in archive/).

# Workers started by parallel::makeCluster() are empty R sessions: they need this package loaded to run its
# functions. An installed package is loaded by name; under devtools::load_all() (development) the source
# folder is loaded with pkgload.
.sp_cluster_setup <- function(cl) {
  ns <- asNamespace("synthplacebo")
  dev <- requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("synthplacebo")
  if (dev) {
    path <- getNamespaceInfo(ns, "path")
    parallel::clusterCall(cl, function(p) { suppressMessages(pkgload::load_all(p, quiet = TRUE, export_all = FALSE)); NULL }, path)
  } else parallel::clusterCall(cl, function() { loadNamespace("synthplacebo"); NULL })
  invisible(NULL)
}

# lapply with an optional PSOCK cluster (n_cores > 1) and an optional progress bar (sequential only).
# FUN's environment is sent to the workers; packages called with :: load there on first use.
dsc_par_lapply <- function(X, FUN, n_cores = 1, progress = FALSE) {
  n <- length(X); if (!n) return(list())
  if (n_cores > 1 && n > 1) {
    cl <- parallel::makeCluster(min(n_cores, n)); on.exit(parallel::stopCluster(cl), add = TRUE)
    .sp_cluster_setup(cl)
    t0 <- Sys.time(); out <- parallel::parLapplyLB(cl, X, FUN)
    if (progress) message("  done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
    return(out)
  }
  if (progress) { pb <- utils::txtProgressBar(min = 0, max = n, style = 3); on.exit(close(pb), add = TRUE) }
  out <- vector("list", n)
  for (i in seq_len(n)) { out[[i]] <- FUN(X[[i]]); if (progress) utils::setTxtProgressBar(pb, i) }
  out
}

# Covariate balance of synthetic controls (used by dsc_placebo(); the same function is in dsc_placebo_lto.R).
# Z: units x covariates (pre-treatment means; rownames = units); w_T: the treated unit's weights, named by
# donor; W_P: list, one per placebo (named by placebo), of its weights named by its donors. Each gap is
# standardized by the donors' SD; share_placebos_worse is the share of placebos with a larger absolute gap.
covariate_balance <- function(Z, treated, w_T, W_P) {
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

##############################################################################
# dsc_placebo(): demeaned synthetic control with placebo-calibrated inference
#   (no bootstrap). Full argument and output reference:
#   docs/reference/function_reference.md
#
# Model: intercept + non-negative donor weights summing to one, fitted on the
# pre-treatment period (as in dsc_symcal_boot()).
#
# Inference: the treated unit and every donor (fitted as a placebo with the
# treated unit removed) get standardized gaps k_{u,t} = |gap_{u,t}| / s_{u,t},
# where s depends on `scale` (cumulative effect: cumulative gap, same s for
# rmspe/raw, cumulative synthetic for relative). With N donors and
# m = floor((1 - conf_level) x (N + 1)):
#   pointwise     k*_t = m-th largest donor k_{p,t}; interval gap_t +/- k*_t x s_{T,t};
#                 p = rank of k_{T,t} among the donors (pre-specified test_period)
#   simultaneous  q~_t = loess-smoothed (log scale) median donor k_{p,t};
#                 lambda = m-th largest donor max_t k_{p,t} / q~_t;
#                 band gap_t +/- lambda x q~_t x s_{T,t}
# Prediction intervals: synthetic_t +/- the same half-widths. With m = 0 (fewer
# than 19 donors at 95%) the intervals are infinite.
#
# Scale (`scale`): what each gap is divided by before ranking against placebos.
#   "rmspe" (default)     the unit's pre-treatment RMSPE (floored, see the guard
#                         below): gaps standardized by pre-treatment fit (the logic of
#                         the RMSPE ratio). In simulation it stayed largely valid when
#                         the treated unit was not exchangeable with the donors, where
#                         the other scales failed.
#   "relative"            the unit's synthetic value: percentage gaps. Most powerful
#                         when the treated unit resembles the donors; not robust when
#                         it does not. Needs positive synthetic values.
#   "raw"                 nothing: raw gaps in the outcome's units, as in the classic
#                         placebo plots of Abadie, Diamond and Hainmueller.
#   Whatever the scale, $test also reports the classic post/pre RMSPE-ratio test
#   (whole post period, p-value only).
#
# Arguments: data (long, balanced panel), unit_col, time_cols (coarsest first, e.g.
#   c("year", "qtr")), outcome_col, treated_unit, treatment_period (named vector:
#   first treated period), test_period (same format; default: last post period),
#   conf_level (0.95), span (loess span of the simultaneous profile, 0.75), scale,
#   rmspe_floor, adjust_covariates / adjust_model / ife_r, augment (with ridge_c /
#   ridge_grid), estimator, diagnostics / warn / ife_seed, balance_covariates, n_cores / verbose -
#   described below and in the reference.
#
# Covariate adjustment (optional): adjust_covariates = column name(s) of time-
#   varying covariates. Step 1 estimates their effect gamma on UNTREATED cells
#   only (the treated unit's post period is excluded), then every unit's outcome
#   is adjusted (outcome - gamma x its actual covariate) and the method runs on
#   the adjusted outcomes. Each placebo re-estimates gamma with the treated unit
#   dropped and its own post period excluded.
#     adjust_model = "ife" (default)  interactive fixed effects via fect (unit and
#                                     time effects + latent factors); the number
#                                     of factors is chosen by cross-validation
#                                     over ife_r for the treated unit's analysis
#                                     and reused for every placebo
#     adjust_model = "twfe"           two-way fixed effects (no factors)
#   Assumes the treatment does not affect the covariate and a common linear effect.
#   $adjustment reports gamma, the number of factors, the two-way fixed-effects
#   gamma (with a unit-clustered SE) for comparison, and the unadjusted estimates.
#
# Augmentation (optional): augment = "ife" gives the AUGMENTED demeaned synthetic
#   control. An interactive fixed-effects outcome model (fect: unit and time
#   effects + latent factors, plus adjust_covariates if given) is fitted on
#   UNTREATED cells only; the demeaned synthetic control is then fitted to the
#   residuals (outcome - model prediction), and the counterfactual is the model
#   prediction plus the synthetic control of the residuals. Equivalently:
#   demeaned SC + (lambda_T - sum_j w_j lambda_j) x f_t, i.e. a correction for the
#   factor exposure the weights cannot match (e.g. a treated unit outside the
#   donors' range). The number of factors is chosen by cross-validation over
#   ife_r for the treated unit's analysis and reused for every placebo; with 0
#   factors the estimator reduces exactly to the demeaned synthetic control.
#   Each placebo refits the outcome model with the treated unit dropped and its
#   own post period excluded.
#   augment = "ridge" gives the RIDGE-AUGMENTED demeaned synthetic control (in the
#   spirit of Ben-Michael, Feller and Rothstein's augmented SCM). Each post-period
#   counterfactual is the demeaned SC plus (x_T - sum_j w_j x_j)' beta_t, where x_j is
#   unit j's demeaned pre-period path and beta_t a ridge regression, across donors,
#   of their period-t outcome on their pre-period paths: it extrapolates the pre-period
#   imbalance the weights leave (a treated unit at or beyond the edge of the donors'
#   range) and is about zero when the weights already match. The penalty is ridge_c x
#   the donors' mean squared path norm; ridge_c = 0.1 by default (chosen in
#   Experimental_inference/LTO/simulations/ridge_default_validation.R: as valid as cross-validation in
#   simulation, lower error on real panels, and stable). ridge_c = "cv" (or NULL) chooses it
#   from ridge_grid by cross-validation inside the pre-period instead (each donor predicts its
#   last min(n_post, T0/4) pre-periods from the earlier ones; needs 8+ pre-periods). With
#   adjust_covariates (either adjust_model) the outcomes are adjusted first and the
#   ridge model is fitted to the adjusted outcomes. Pre-treatment gaps (and the RMSPE
#   scale) are the demeaned SC's; each placebo refits its own correction with the same
#   penalty. $ridge reports the penalty, the cross-validation errors and the
#   correction. Validated in simulations/ridge_*_validation.R: valid with 9-39 donors
#   and 10-80 pre-periods; restores coverage for treated units at or beyond the edge
#   of the donors' range with trending factors (89-95% at 90%, vs 49-84% without);
#   lower RMSE than the plain estimator almost everywhere; with IFE covariate
#   adjustment the best combination tested for a confounded treated unit.
#
# Guard (rmspe_floor, RMSPE scale only): each unit's pre-treatment RMSPE is floored at
#   rmspe_floor x the median pre-treatment RMSPE over the treated unit and all placebos
#   (default 0.25; 0 switches the guard off). Without it, a placebo that the demeaned
#   synthetic control reproduces (almost) exactly before treatment - likely with short
#   pre-periods and many donors - gets near-infinite standardized gaps and so do the
#   intervals. $rmspe_guard reports the floor and the units it applied to. The classic
#   RMSPE-ratio test is reported unguarded. Validated in simulations/rmspe_guard_validation.R.
#
# Estimator (estimator): "dsc" (default) is the demeaned synthetic control. "fect"
#   replaces the counterfactual with fect's interactive fixed-effects (IFEct)
#   prediction: unit and time effects + latent factors (+ adjust_covariates), fitted on
#   untreated cells only, with the number of factors chosen by cross-validation over
#   ife_r for the treated unit and reused for every placebo (the treated unit's model is
#   then refitted with that number, exactly as each placebo's is). Each donor is refitted as
#   a pseudo-treated unit (treated unit dropped), and the usual placebo calibration is
#   applied. In simulation, fect with placebo calibration on the RMSPE scale was valid
#   and much more powerful than DSC when the treated unit lay at or beyond the edge of
#   the donors' range with stationary factors; it was not valid with trending factors.
#   Use scale = "rmspe" with it. Cannot be combined with augment = "ife" or "ridge". With
#   estimator = "fect", $weights and $intercept are NULL (there are no donor weights).
#
# Run time and parallel fitting (n_cores, verbose): with adjust_model = "ife", augment = "ife" or
#   estimator = "fect", every placebo needs its own fect model, which dominates the run time (Kansas, 48
#   donors, three covariates: 10.9 minutes, 97% in fect; about 11 s per fit and 35 s for the cross-
#   validation). n_cores > 1 fits the placebo models on that many cores (base R parallel; workers are
#   started and stopped inside the call, costing a few seconds); results are identical to n_cores = 1
#   because every fect fit sets its own seed. n_cores is ignored when there are no fect models (the
#   fits then take seconds). verbose (default interactive()) shows the slow stages and a progress bar.
# Covariate balance (balance_covariates; default: the adjust_covariates): a check, not an adjustment.
#   Each covariate is summarized by its pre-treatment mean per unit. The weights fitted on the outcome
#   path are applied to the donors' values, and the treated unit's gap from its synthetic control is
#   standardized by the donors' SD; each placebo's gap is computed with its own weights, and the treated
#   unit is flagged when fewer than 10% of placebos have a larger gap. Name any numeric column (time-
#   invariant characteristics too); they are not adjusted for. A flagged imbalance biases the estimate
#   only if the covariate's effect on the outcome changes over time. Reported in $balance (estimator =
#   "dsc"; the weights are the synthetic control's, before any ridge correction).
# Diagnostics (diagnostics = TRUE, default): checks whether the treated unit lies
#   beyond what the donors can reproduce, the case where placebo-calibrated
#   intervals were found to under-cover in simulation. $diagnostics reports:
#     fit_rank         treated pre-treatment RMSPE relative to its level, and the
#                      share of placebos that fit worse (flag: worse than 90%)
#     hull             distance of the treated unit's estimated factor loadings
#                      (fect IFE on untreated cells) from the donors' convex hull,
#                      against each donor's distance from the other donors' hull
#                      (flag: beyond 90% of donors)
#     extrapolation    post-treatment drift implied by the loading mismatch,
#                      (lambda_T - sum_j w_j lambda_j) x f_t, as % of the level
#     persistence      first-order autocorrelation of each estimated factor
#                      (trending if above 0.95)
#     residual_autocorrelation
#                      persistence of the untreated gaps: pooled AR(1) of the placebos'
#                      post-treatment gaps (out-of-sample for the untreated donors),
#                      each standardized by its pre-treatment RMSPE; 90% bootstrap
#                      interval over placebos; level (low < 0.3 <= moderate < 0.7 <=
#                      high < 0.9 <= very high), half-life, periods for a shock to fall
#                      to 25% (a backdating guide), effective number of independent
#                      pre-periods T0 (1 - rho) / (1 + rho), the in-sample pre-treatment
#                      value (understates persistence; the difference shows how much the
#                      weights absorb), seasonal lag (finest time column with 2-12
#                      values; reported as lag-s autocorrelation minus rho^s, the seasonal
#                      persistence beyond lag 1) and a guide sentence. backdate_periods is NA
#                      when rho >= 0.9. Information only (no flag).
#   fect's cross-validation is random; ife_seed (default 12112025) fixes it locally (the
#   session's random-number state is restored), so results are reproducible.
#   A warning is issued (warn = TRUE) when the treated unit is flagged as extreme;
#   it is strongest when the factors are trending. Loading-based checks need fect.
#   $diagnostics$recommendation gives the arguments to re-run with when the
#   diagnostics point to a better-suited method: a treated unit flagged as extreme
#   with trending factors -> re-run with augment = "ridge"; with factors NOT trending
#   (and estimator = "dsc", augment not "ridge") -> estimator = "fect", scale =
#   "rmspe". NULL when no change is recommended.
#   Validated on 1,800 simulated panels (diagnostics_validation.R): the extreme flag
#   caught 100% of 'outside' and 81% of 'edge' treated units, with 4-5% false alarms
#   for treated units inside the donors' hull; runs with the combined warning (flag
#   AND trending) had 74% coverage (DSC, RMSPE scale) vs 96% for the others.
#
# Value: a list with effects (one row per period: observed, synthetic, effect, effect
#   and prediction intervals _pw / _sim, p_pointwise, and the cumulative versions cum_*),
#   test (pointwise, simultaneous and RMSPE-ratio tests), weights, intercept, pre_rmspe,
#   rmspe_guard, adjustment (with adjust_covariates), ridge (with augment = "ridge"), placebo (gaps and synthetic
#   values), diagnostics, failed_donors and settings. See the reference for every column.
#   effects also holds std_effect (the gap on the calibration scale) and threshold_ratio_pw /
#   _sim (|std_effect| / the pointwise or simultaneous threshold; above 1 = outside the
#   untreated units' chance range), with cum_ versions: the fit-based benchmark.
##############################################################################
#' Demeaned synthetic control with placebo-calibrated inference
#'
#' Fits the demeaned synthetic control (an intercept plus non-negative donor weights that sum to one) to a
#' single treated unit's pre-treatment outcomes, and calibrates pointwise intervals, simultaneous bands and
#' prediction intervals against placebo fits of every donor (each fitted on the other donors, with the treated
#' unit removed from the pool). Optional time-varying covariate adjustment, ridge or interactive fixed-effects
#' augmentation, an interactive fixed-effects estimator, a covariate balance check and diagnostics for treated
#' units beyond the donors' range.
#'
#' @details
#' Each unit's post-treatment gaps are put on a common scale (by default divided by the unit's own
#' pre-treatment RMSPE) and the treated unit's gaps are ranked among the placebos'. With \eqn{N} donors and
#' \eqn{m = \lfloor (1 - conf\_level)(N + 1) \rfloor}{m = floor((1 - conf_level)(N + 1))}:
#' * **pointwise**: the threshold in each period is the \eqn{m}-th largest placebo value; the interval is the
#'   gap plus or minus the threshold times the treated unit's scale;
#' * **simultaneous**: a smoothed profile of the median placebo value over time and one multiplier, so that
#'   the band covers the whole post-treatment path at once.
#'
#' With fewer than 19 donors a 95% interval is infinite (the smallest p-value is \eqn{1/(N + 1)}); use
#' [dsc_placebo_lto()] for small donor pools.
#'
#' **Run time.** A few seconds for about 50 donors without covariate adjustment. Covariate adjustment with
#' `adjust_model = "ife"`, `augment = "ife"` and `estimator = "fect"` fit one `fect` model per placebo
#' (Kansas, 48 donors, three covariates: about 11 minutes, 3.6 minutes with `n_cores = 4`).
#'
#' @param data A data frame in long format: one row per unit and period. The panel must be balanced, with no
#'   missing values in the outcome or in `adjust_covariates`.
#' @param unit_col Name of the column identifying units.
#' @param time_cols Name(s) of the column(s) identifying periods, from the coarsest to the finest, e.g.
#'   `"year"` or `c("year", "qtr")`.
#' @param outcome_col Name of the outcome column (numeric).
#' @param treated_unit The treated unit: a value of `unit_col`. All other units are donors.
#' @param treatment_period The first treated period, as a named vector with one element per time column,
#'   e.g. `c(year = 1991)` or `c(year = 2012, qtr = 2)`.
#' @param test_period The period of the pointwise tests, in the same format (`NULL`: the last post-treatment
#'   period). The cumulative test sums the effect from the start of treatment through it. Choose it before
#'   looking at the results.
#' @param conf_level Confidence level of the intervals and bands.
#' @param span LOESS span for smoothing the simultaneous band's time profile (log scale).
#' @param scale What each unit's gaps are divided by before ranking: `"rmspe"` (the unit's pre-treatment
#'   RMSPE; the most robust when the treated unit is unlike the donors), `"relative"` (its synthetic value:
#'   percentage gaps; needs positive outcomes) or `"raw"` (nothing).
#' @param adjust_covariates Time-varying covariate columns to adjust for. Their effect \eqn{\gamma} is estimated
#'   on untreated cells only; every unit's outcome is replaced by outcome minus \eqn{\hat\gamma} times the
#'   covariate. Each placebo re-estimates \eqn{\gamma} without its own post-treatment period. Assumes the
#'   treatment does not affect the covariates; see [covariate_trajectory_check()].
#' @param adjust_model Model for \eqn{\gamma}: `"ife"` (interactive fixed effects via `fect`) or `"twfe"`
#'   (two-way fixed effects).
#' @param ife_r Candidate numbers of latent factors for every interactive fixed-effects model, chosen by
#'   cross-validation for the treated unit and reused for every placebo. A single value skips cross-validation.
#' @param augment `"none"`; `"ridge"` (ridge-augmented demeaned synthetic control: corrects the counterfactual
#'   for pre-treatment imbalance the weights leave, for a treated unit at or beyond the edge of the donors'
#'   range); or `"ife"` (augmented with an interactive fixed-effects outcome model; needs `fect`).
#' @param diagnostics Check whether the treated unit lies beyond what the donors can reproduce (fit rank,
#'   distance from the donors' convex hull in factor loadings, implied extrapolation, factor persistence) and
#'   report the residual autocorrelation of the untreated gaps. The loading-based checks need `fect`.
#' @param warn Warn when the diagnostics flag the treated unit as extreme.
#' @param ife_seed Seed for `fect`'s cross-validation, applied locally (the session's random-number state is
#'   restored).
#' @param estimator `"dsc"` (the demeaned synthetic control) or `"fect"` (the counterfactual is `fect`'s
#'   interactive fixed-effects prediction, with the same placebo calibration; use with `scale = "rmspe"`).
#' @param rmspe_floor For `scale = "rmspe"`: each unit's pre-treatment RMSPE is floored at `rmspe_floor` times
#'   the median over the treated unit and the placebos, so a unit fitted almost exactly before treatment
#'   cannot produce near-infinite intervals. `0` switches it off.
#' @param ridge_c For `augment = "ridge"`: the penalty as a multiple of the donors' mean squared pre-treatment
#'   path norm (default 0.1, chosen by simulation). `"cv"` or `NULL` chooses it from `ridge_grid` by
#'   cross-validation inside the pre-treatment period (needs at least 8 pre-treatment periods).
#' @param ridge_grid Candidate values of `ridge_c` for the cross-validation.
#' @param balance_covariates Numeric columns whose pre-treatment balance is checked (default: the
#'   `adjust_covariates`). A check, not an adjustment; time-invariant characteristics can be named too.
#' @param n_cores Cores for the per-placebo `fect` models (results are identical to `n_cores = 1`). Ignored
#'   when there are no `fect` models.
#' @param verbose Report the slow stages and show a progress bar.
#'
#' @return An object of class `"dsc_placebo"`: a list with
#' * `effects`: one row per period, with `observed`, `synthetic`, `effect`, effect intervals (`eff_low_pw`,
#'   `eff_high_pw`, `eff_low_sim`, `eff_high_sim`), prediction intervals (`pred_*`), `p_pointwise`, the
#'   calibration-scale values (`std_effect`, `factor_pointwise`, `factor_simultaneous`, `threshold_ratio_*`)
#'   and the cumulative versions (`cum_*`);
#' * `test`: the pointwise tests of the effect and the cumulative effect at `test_period`, the simultaneous
#'   tests, and the classic post/pre RMSPE-ratio test;
#' * `weights`, `intercept`, `pre_rmspe`, `rmspe_guard`;
#' * `adjustment` (with `adjust_covariates`), `ridge` (with `augment = "ridge"`), `balance`;
#' * `placebo` (each placebo's gaps and synthetic values), `diagnostics` (including an `assessment` sentence and
#'   a `recommendation`), `failed_donors` and `settings`.
#'
#' Printing shows a compact summary; [summary.dsc_placebo()] gives more detail and [plot.dsc_placebo()] draws
#' the standard figures.
#' @seealso [dsc_placebo_lto()] for small donor pools, [synth_power()] for power,
#'   [covariate_trajectory_check()] before adjusting for covariates, [plot_full_path()].
#' @references Abadie, A., Diamond, A., and Hainmueller, J. (2010). Synthetic control methods for comparative
#'   case studies. *Journal of the American Statistical Association*, 105(490), 493-505.
#' @examples
#' \donttest{
#' if (requireNamespace("scpi", quietly = TRUE)) {
#'   data("scpi_germany", package = "scpi")
#'   germany <- as.data.frame(scpi_germany)[, c("country", "year", "gdp")]
#'   fit <- dsc_placebo(germany, "country", "year", "gdp", treated_unit = "West Germany",
#'                      treatment_period = c(year = 1991), conf_level = 0.90, diagnostics = FALSE)
#'   fit
#'   summary(fit)
#'   plot(fit, type = "effect")
#' }
#' }
#' @export
dsc_placebo <- function(data, unit_col, time_cols, outcome_col, treated_unit, treatment_period,
                        test_period = NULL, conf_level = 0.95, span = 0.75,
                        scale = c("rmspe", "relative", "raw"),
                        adjust_covariates = NULL, adjust_model = c("ife", "twfe"), ife_r = 0:5,
                        augment = c("none", "ife", "ridge"), diagnostics = TRUE, warn = TRUE, ife_seed = 12112025,
                        estimator = c("dsc", "fect"), rmspe_floor = 0.25,
                        ridge_c = 0.1, ridge_grid = 10^seq(-3, 2, by = 0.5), balance_covariates = adjust_covariates,
                        n_cores = 1, verbose = interactive()) {

  if (!requireNamespace("CVXR", quietly = TRUE)) stop("dsc_placebo() needs the CVXR package.")
  scale <- match.arg(scale); adjust_model <- match.arg(adjust_model); augment <- match.arg(augment)
  estimator <- match.arg(estimator)
  if (estimator == "fect" && augment != "none") stop("augment = \"", augment, "\" applies to estimator = \"dsc\" only.")
  if (identical(ridge_c, "cv")) ridge_c <- NULL                  # "cv" (or NULL): choose the penalty by cross-validation
  if (!is.null(ridge_c) && (!is.numeric(ridge_c) || length(ridge_c) != 1 || !is.finite(ridge_c) || ridge_c <= 0)) stop("ridge_c must be one positive number, or \"cv\" to choose it by cross-validation.")
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) stop("n_cores must be a single number >= 1.")
  say <- function(...) if (isTRUE(verbose)) message("dsc_placebo: ", ...)
  if (estimator == "fect" && !requireNamespace("fect", quietly = TRUE)) stop("estimator = \"fect\" needs the fect package.")
  if (augment == "ife" && !requireNamespace("fect", quietly = TRUE)) stop("augment = \"ife\" needs the fect package.")
  if (length(adjust_covariates) && adjust_model == "ife" && !requireNamespace("fect", quietly = TRUE)) stop("adjust_model = \"ife\" needs the fect package.")
  alpha <- round(1 - conf_level, 10)

  # ---- Data -------------------------------------------------------------------------
  data <- as.data.frame(data)
  stopifnot(all(c(unit_col, time_cols, outcome_col, adjust_covariates) %in% names(data)), all(names(treatment_period) %in% time_cols))
  if (length(balance_covariates)) {
    if (!all(balance_covariates %in% names(data))) stop("balance_covariates not found in data: ", paste(setdiff(balance_covariates, names(data)), collapse = ", "))
    if (!all(vapply(data[balance_covariates], is.numeric, logical(1)))) stop("balance_covariates must be numeric columns.")
  }
  periods <- unique(data[, time_cols, drop = FALSE])
  periods <- periods[do.call(order, unname(as.list(periods))), , drop = FALSE]; rownames(periods) <- NULL
  find_period <- function(pv, what) {
    hit <- which(Reduce(`&`, lapply(names(pv), function(nm) periods[[nm]] == pv[[nm]])))
    if (length(hit) != 1) stop(what, " (", paste(names(pv), pv, sep = " = ", collapse = ", "), ") must match exactly one period.")
    hit
  }
  t_treat <- find_period(treatment_period, "treatment_period")
  post <- seq_len(nrow(periods)) >= t_treat; n_post <- sum(post)
  t_test <- if (is.null(test_period)) n_post else find_period(test_period, "test_period") - t_treat + 1
  if (t_test < 1) stop("test_period must be a post-treatment period.")
  key <- do.call(paste, unname(as.list(data[, time_cols, drop = FALSE])))
  d <- data.frame(unit = as.character(data[[unit_col]]), .t = match(key, do.call(paste, unname(as.list(periods)))), y = data[[outcome_col]])
  for (v in union(adjust_covariates, balance_covariates)) d[[v]] <- data[[v]]
  if (anyNA(d[, c("y", adjust_covariates)]) || any(table(d$unit) != nrow(periods))) stop("The panel must be balanced with no missing outcomes or covariates.")
  if (!treated_unit %in% d$unit) stop("treated_unit not found in ", unit_col, ".")
  units <- sort(unique(d$unit)); donors <- setdiff(units, treated_unit)
  wide <- function(v) sapply(units, function(u) { du <- d[d$unit == u, ]; du[[v]][order(du$.t)] })
  Y <- wide("y")
  Aobs <- setNames(lapply(adjust_covariates, wide), adjust_covariates)
  Zbal <- if (length(balance_covariates)) {           # pre-treatment mean of each balance covariate, units x covariates
    Z <- do.call(cbind, lapply(balance_covariates, function(v) colMeans(wide(v)[!post, , drop = FALSE], na.rm = TRUE)))
    dimnames(Z) <- list(units, balance_covariates); Z } else NULL

  # ---- Step 1 (optional): covariate effect gamma on untreated cells --------------------------
  # Long data for the units kept; D = 1 marks the excluded cells (drop_post_unit's post period)
  long_for <- function(drop_units, drop_post_unit) {
    keep <- setdiff(units, drop_units)
    out <- expand.grid(.t = seq_len(nrow(periods)), unit = keep, stringsAsFactors = FALSE)
    out$y <- Y[cbind(out$.t, match(out$unit, units))]
    for (v in adjust_covariates) out[[v]] <- Aobs[[v]][cbind(out$.t, match(out$unit, units))]
    out$D <- as.integer(out$unit == drop_post_unit & post[out$.t])
    out
  }
  gamma_twfe <- function(drop_units = NULL, drop_post_unit) {
    lg <- long_for(drop_units, drop_post_unit); lg <- lg[lg$D == 0, ]
    X <- model.matrix(as.formula(paste("~ factor(unit) + factor(.t) +", paste(adjust_covariates, collapse = " + "))), lg)
    qx <- qr(X); b <- qr.coef(qx, lg$y)[adjust_covariates]
    if (anyNA(b)) stop("adjust_covariates: gamma not identified (no variation beyond unit and period effects).")
    e <- lg$y - X[, !is.na(qr.coef(qx, lg$y)), drop = FALSE] %*% qr.coef(qx, lg$y)[!is.na(qr.coef(qx, lg$y))]
    Xk <- X[, !is.na(qr.coef(qx, lg$y)), drop = FALSE]; XtX_inv <- solve(crossprod(Xk))
    meat <- Reduce(`+`, lapply(split(seq_len(nrow(Xk)), lg$unit), function(ix) { sc <- crossprod(Xk[ix, , drop = FALSE], e[ix]); sc %*% t(sc) }))
    G <- length(unique(lg$unit)); V <- G / (G - 1) * XtX_inv %*% meat %*% XtX_inv
    list(gamma = b, se = sqrt(diag(V))[match(adjust_covariates, colnames(Xk))], n_obs = nrow(lg))
  }
  gamma_ife <- function(drop_units = NULL, drop_post_unit, r = ife_r, cv = TRUE) {
    lg <- long_for(drop_units, drop_post_unit)
    if (length(unique(r)) == 1) { cv <- FALSE; r <- unique(r) }   # one candidate: no cross-validation
    fit <- NULL
    utils::capture.output(fit <- with_seed(suppressMessages(suppressWarnings(fect::fect(
      as.formula(paste("y ~ D +", paste(adjust_covariates, collapse = " + "))), data = lg, index = c("unit", ".t"),
      method = "ife", r = if (cv) range(r) else r, CV = cv, force = "two-way", se = FALSE, parallel = FALSE)))))
    b <- setNames(as.numeric(fit$beta[adjust_covariates, 1]), adjust_covariates)
    list(gamma = b, r = if (cv) as.integer(fit$r.cv) else as.integer(r), n_obs = sum(lg$D == 0))
  }
  with_seed <- function(expr) {                       # fect's cross-validation is random: fix it locally, restore the session's RNG
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
    set.seed(ife_seed); expr
  }
  model_ife <- function(drop_units = NULL, drop_post_unit, r = ife_r, cv = TRUE) {
    lg <- long_for(drop_units, drop_post_unit)
    if (length(unique(r)) == 1) { cv <- FALSE; r <- unique(r) }   # one candidate: no cross-validation
    rhs <- paste(c("D", adjust_covariates), collapse = " + ")
    fit <- NULL
    utils::capture.output(fit <- with_seed(suppressMessages(suppressWarnings(fect::fect(
      as.formula(paste("y ~", rhs)), data = lg, index = c("unit", ".t"),
      method = "ife", r = if (cv) range(r) else r, CV = cv, force = "two-way", se = FALSE, parallel = FALSE)))))
    M <- if (is.matrix(fit$Y.ct.full)) fit$Y.ct.full else fit$Y.ct          # Y.ct.full is NULL after cross-validation
    if (!is.matrix(M) || !all(dim(M) == c(nrow(periods), length(unique(lg$unit))))) stop("augment: unexpected fect prediction matrix")
    colnames(M) <- as.character(fit$id)
    lam <- if (!is.null(fit$lambda)) { L <- as.matrix(fit$lambda); rownames(L) <- as.character(fit$id); L } else NULL
    list(M = M, r = if (cv) as.integer(fit$r.cv) else as.integer(r), n_obs = sum(lg$D == 0),
         lambda = lam, factor = if (!is.null(fit$factor)) as.matrix(fit$factor) else NULL,
         gamma = if (length(adjust_covariates)) setNames(as.numeric(fit$beta[adjust_covariates, 1]), adjust_covariates) else NULL)
  }
  adjust_mat <- function(cols, gam) {
    Ya <- Y[, cols, drop = FALSE]
    if (!is.null(gam)) for (v in names(gam)) Ya <- Ya - gam[[v]] * Aobs[[v]][, cols, drop = FALSE]
    Ya
  }
  adj_main <- NULL; r_use <- NA_integer_; tw_main <- NULL; aug_main <- NULL
  if (augment == "ife" || estimator == "fect") {
    say("choosing the number of factors by cross-validation (fect; the slowest single step)")
    r_use <- model_ife(drop_post_unit = treated_unit)$r                  # cross-validation chooses the number of factors
    aug_main <- model_ife(drop_post_unit = treated_unit, r = r_use, cv = FALSE)   # refit with it, exactly as for the placebos
    adj_main <- list(gamma = aug_main$gamma, n_obs = aug_main$n_obs)
    if (length(adjust_covariates)) tw_main <- gamma_twfe(drop_post_unit = treated_unit)
  } else if (length(adjust_covariates)) {
    tw_main <- gamma_twfe(drop_post_unit = treated_unit)
    if (adjust_model == "ife") { say("choosing the number of factors by cross-validation (fect; the slowest single step)"); r_use <- gamma_ife(drop_post_unit = treated_unit)$r
                                  adj_main <- gamma_ife(drop_post_unit = treated_unit, r = r_use, cv = FALSE) } else adj_main <- tw_main
  }
  gamma_for_placebo <- function(p) {
    if (!length(adjust_covariates)) return(NULL)
    if (adjust_model == "ife") gamma_ife(drop_units = treated_unit, drop_post_unit = p, r = r_use, cv = FALSE)$gamma
    else gamma_twfe(drop_units = treated_unit, drop_post_unit = p)$gamma
  }

  # ---- Demeaned SC fit (rescaled before solving) -------------------------------------------
  fit_dsc <- function(X, yv, fit = !post) {                   # fit: the periods the weights are fitted on
    sc <- mean(abs(yv[fit]))
    Xp <- X[fit, , drop = FALSE] / sc; yp <- yv[fit] / sc
    w <- CVXR::Variable(ncol(Xp)); b0 <- CVXR::Variable(1)
    prob <- CVXR::Problem(CVXR::Minimize(CVXR::sum_squares(yp - b0 - Xp %*% w)), list(w >= 0, sum(w) == 1))
    CVXR::psolve(prob)
    st <- CVXR::status(prob)
    if (!st %in% c("optimal", "optimal_inaccurate")) {        # retry once: explicit solver, more iterations
      try(CVXR::psolve(prob, solver = "CLARABEL", max_iter = 10000), silent = TRUE)
      st <- CVXR::status(prob)
    }
    if (!st %in% c("optimal", "optimal_inaccurate")) stop("CVXR status: ", st)
    wv <- setNames(as.numeric(CVXR::value(w)), colnames(X)); wv[abs(wv) < 1e-8] <- 0
    b <- as.numeric(CVXR::value(b0)) * sc
    list(b0 = b, w = wv, synth = b + as.numeric(X %*% wv))
  }

  # ---- Ridge augmentation (augment = "ridge") ------------------------------------------------
  # Counterfactual = demeaned SC + (x_T - sum_j w_j x_j)' beta_t in each post period: x_j is unit j's
  # demeaned pre-period path and beta_t a ridge regression, across the donors, of their (demeaned,
  # centred) period-t outcome on their pre-period paths. The penalty is ridge_c x the donors' mean
  # squared (centred) path norm. The correction extrapolates the pre-period imbalance the weights
  # leave; it is ~0 when the weights already match the treated unit. Pre-treatment gaps (and so the
  # RMSPE scale) are the demeaned SC's. Works on the covariate-adjusted outcomes when adjust_covariates
  # is given. Each placebo refits its own correction with the same penalty.
  ridge_correction <- function(yv, X, w, c_val, fit = !post, out = post) {
    mu <- colMeans(X[fit, , drop = FALSE])
    M <- t(sweep(X[fit, , drop = FALSE], 2, mu)); Mc <- sweep(M, 2, colMeans(M))                 # donors x fit periods
    Yo <- t(sweep(X[out, , drop = FALSE], 2, mu)); Yoc <- sweep(Yo, 2, colMeans(Yo))            # donors x out periods
    dv <- (yv[fit] - mean(yv[fit])) - as.numeric(crossprod(M, w))                                 # imbalance the weights leave
    lam <- c_val * mean(rowSums(Mc^2))
    A <- t(Mc) %*% solve(tcrossprod(Mc) + diag(lam, nrow(Mc)))
    corr <- numeric(length(yv)); corr[out] <- as.numeric(crossprod(dv, A %*% Yoc)); corr
  }
  ridge_info <- NULL; c_ridge <- NULL
  if (augment == "ridge") {
    YaR <- adjust_mat(c(treated_unit, donors), adj_main$gamma)
    pre_idx <- which(!post); T0 <- length(pre_idx)
    if (is.null(ridge_c)) {
      # Cross-validation inside the pre-period: each donor (fitted on the other donors) predicts its
      # last k pre-periods from the earlier ones; the grid value with the lowest mean squared error wins
      k <- min(n_post, floor(T0 / 4))
      if (T0 < 8) stop("augment = \"ridge\" with ridge_c = \"cv\" needs at least 8 pre-treatment periods; use a numeric ridge_c (default 0.1).")
      F_rows <- seq_len(nrow(periods)) %in% head(pre_idx, -k); V_rows <- seq_len(nrow(periods)) %in% tail(pre_idx, k)
      cv_err <- sapply(donors, function(j) tryCatch({
        oth <- setdiff(donors, j); fj <- fit_dsc(YaR[, oth, drop = FALSE], YaR[, j], fit = F_rows)
        g <- (YaR[, j] - fj$synth)[V_rows]
        vapply(ridge_grid, function(cc) mean((g - ridge_correction(YaR[, j], YaR[, oth, drop = FALSE], fj$w, cc, F_rows, V_rows)[V_rows])^2), numeric(1))
      }, error = function(e) rep(NA_real_, length(ridge_grid))))
      cv_mean <- rowMeans(matrix(cv_err, nrow = length(ridge_grid)), na.rm = TRUE)
      c_ridge <- ridge_grid[which.min(cv_mean)]
      cv_table <- data.frame(c = ridge_grid, cv_mse = cv_mean)
    } else { c_ridge <- ridge_c; k <- NA_integer_; cv_table <- NULL }
    at_edge <- is.null(ridge_c) && c_ridge %in% range(ridge_grid)
    ridge_info <- list(c = c_ridge, chosen_by = if (is.null(ridge_c)) "cross-validation" else if (ridge_c == 0.1) "default (0.1)" else "ridge_c", holdout_periods = k,
                       cv = cv_table, at_grid_edge = at_edge,
                       note = if (at_edge) paste0("The cross-validated penalty is at the ", if (c_ridge == max(ridge_grid)) "top" else "bottom",
                         " of ridge_grid (", if (c_ridge == max(ridge_grid)) "the correction is close to switched off" else "the least shrinkage on offer",
                         "); this was more common with about 10 pre-periods and 9 donors in simulation. Compare with augment = \"none\".") else NULL)
  }

  # ---- Treated unit and placebos ------------------------------------------------------------
  YaT <- if (augment == "ife") Y[, c(treated_unit, donors), drop = FALSE] - aug_main$M[, c(treated_unit, donors), drop = FALSE] else
    adjust_mat(c(treated_unit, donors), adj_main$gamma)
  ft <- fit_dsc(YaT[, donors, drop = FALSE], YaT[, treated_unit])
  corr_T <- if (augment == "ridge") ridge_correction(YaT[, treated_unit], YaT[, donors, drop = FALSE], ft$w, c_ridge) else 0
  gT_all <- if (estimator == "fect") Y[, treated_unit] - aug_main$M[, treated_unit] else
    YaT[, treated_unit] - ft$synth - corr_T                                                     # gap (adjusted = original scale)
  if (augment == "ridge") ridge_info$correction <- corr_T
  slow_fits <- (length(adjust_covariates) && adjust_model == "ife") || augment == "ife" || estimator == "fect"   # a fect model per placebo
  if (slow_fits) say("fitting ", length(donors), " placebo models (one fect fit each", if (n_cores > 1) paste0(", on ", min(n_cores, length(donors)), " cores") else "", ")")
  pl <- dsc_par_lapply(donors, function(p) tryCatch({
    if (estimator == "fect")                                                                 # IFEct, donor p as pseudo-treated
      return(list(gap = Y[, p] - model_ife(drop_units = treated_unit, drop_post_unit = p, r = r_use, cv = FALSE)$M[, p]))
    cols_p <- c(p, setdiff(donors, p))
    Yp <- if (augment == "ife") Y[, cols_p, drop = FALSE] - model_ife(drop_units = treated_unit, drop_post_unit = p, r = r_use, cv = FALSE)$M[, cols_p, drop = FALSE] else
      adjust_mat(cols_p, gamma_for_placebo(p))
    fp <- fit_dsc(Yp[, setdiff(donors, p), drop = FALSE], Yp[, p])
    corr_p <- if (augment == "ridge") ridge_correction(Yp[, p], Yp[, setdiff(donors, p), drop = FALSE], fp$w, c_ridge) else 0
    list(gap = Yp[, p] - fp$synth - corr_p, w = stats::setNames(as.numeric(fp$w), setdiff(donors, p)))
  }, error = function(e) conditionMessage(e)), n_cores = if (slow_fits) n_cores else 1, progress = slow_fits && isTRUE(verbose))
  failed <- vapply(pl, is.character, logical(1))
  failed_donors <- if (any(failed)) data.frame(donor = donors[failed], error = unlist(pl[failed])) else
    data.frame(donor = character(0), error = character(0))
  if (any(failed)) warning(sum(failed), " placebo fit(s) failed and were left out: ", paste(donors[failed], collapse = ", "), ". See $failed_donors.")
  pl <- pl[!failed]; dn <- donors[!failed]; N <- length(pl); m <- floor(alpha * (N + 1))
  kth <- function(x) if (m == 0) Inf else sort(x, decreasing = TRUE)[m]
  if (m == 0) warning("With ", N, " donors the ", 100 * conf_level, "% level cannot be reached exactly: intervals are infinite.")

  GP <- sapply(pl, `[[`, "gap"); colnames(GP) <- dn; SP <- Y[, dn, drop = FALSE] - GP                            # gaps and counterfactuals, all periods x donors
  synth_T <- Y[, treated_unit] - gT_all                                                         # counterfactual on the original scale
  balance <- if (!is.null(Zbal) && estimator == "dsc")
    covariate_balance(Zbal, treated_unit, stats::setNames(as.numeric(ft$w), donors), stats::setNames(lapply(pl, `[[`, "w"), dn)) else NULL
  sT <- synth_T[post]; gT <- gT_all[post]
  rmspe_T <- sqrt(mean(gT_all[!post]^2)); rmspe_P <- sqrt(colMeans(GP[!post, , drop = FALSE]^2))
  # Guard (RMSPE scale): floor every unit's pre-treatment RMSPE at rmspe_floor x the median over all
  # units, so a placebo reproduced (almost) exactly before treatment cannot get near-infinite
  # standardized gaps. The same floor applies to the treated unit and every placebo.
  floor_val <- if (scale == "rmspe" && rmspe_floor > 0) rmspe_floor * stats::median(c(rmspe_T, rmspe_P)) else 0
  s_T <- max(rmspe_T, floor_val); s_P <- pmax(rmspe_P, floor_val)
  floored_units <- c(if (rmspe_T < floor_val) treated_unit, dn[rmspe_P < floor_val])
  if (scale == "relative" && (any(sT <= 0) || any(SP[post, ] <= 0)))
    stop("Relative gaps need positive synthetic values; use scale = \"rmspe\" (the default) or \"raw\" for outcomes that can be zero or negative.")
  # Scale for the treated unit (post vector) and placebos (post x donors), per period and cumulative
  ones <- matrix(1, n_post, length(dn))
  sc_T <- switch(scale, relative = sT, raw = rep(1, n_post), rmspe = rep(s_T, n_post))
  sc_Tc <- switch(scale, relative = cumsum(sT), raw = rep(1, n_post), rmspe = rep(s_T, n_post))
  sc_P <- switch(scale, relative = SP[post, , drop = FALSE], raw = ones, rmspe = sweep(ones, 2, s_P, "*"))
  sc_Pc <- switch(scale, relative = apply(SP[post, , drop = FALSE], 2, cumsum), raw = ones, rmspe = sweep(ones, 2, s_P, "*"))
  R  <- matrix(abs(GP[post, , drop = FALSE]) / sc_P, nrow = n_post)                            # standardized |gaps|, post x donors
  Rc <- matrix(abs(apply(GP[post, , drop = FALSE], 2, cumsum)) / sc_Pc, nrow = n_post)
  rT <- gT / sc_T; rTc <- cumsum(gT) / sc_Tc

  smooth_pos <- function(v) { tt <- seq_along(v); exp(as.numeric(predict(stats::loess(log(v) ~ tt, span = span, degree = 1)))) }
  calib <- function(K, kT) {
    q <- smooth_pos(apply(K, 1, median)); Rp <- apply(K / q, 2, max)
    list(pw = apply(K, 1, kth), sim = kth(Rp) * q,
         p_pw = (1 + rowSums(K >= abs(kT))) / (N + 1), p_sim = (1 + sum(Rp >= max(abs(kT) / q))) / (N + 1))
  }
  cp <- calib(R, rT); cc <- calib(Rc, rTc)

  # ---- Diagnostics: is the treated unit beyond what the donors can reproduce? -----------------
  diag <- NULL
  if (diagnostics) {
    rel_rmspe <- function(gap, yv) sqrt(mean(gap[!post]^2)) / mean(abs(yv[!post]))
    fit_T <- rel_rmspe(gT_all, Y[, treated_unit])
    fit_P <- vapply(seq_along(dn), function(j) rel_rmspe(GP[, j], Y[, dn[j]]), numeric(1))
    fit_share_worse <- mean(fit_P > fit_T)
    diag <- list(fit_rank = list(treated_rel_rmspe = fit_T, placebo_rel_rmspe_median = median(fit_P),
                                 share_placebos_worse = fit_share_worse, flag = fit_share_worse < 0.10))
    # Residual autocorrelation: how persistent the untreated gaps are. Estimated from the placebos'
    # POST-treatment gaps (the donors are untreated and those periods lie outside their fit, so these
    # are genuine out-of-sample errors), each standardized by its pre-treatment RMSPE and pooled as an
    # AR(1) through the origin. In-sample pre-treatment residuals understate persistence badly (the
    # weights absorb smooth noise: in simulation a true 0.7 read as about 0.2 with 20 pre-periods),
    # so they are reported only to show how much persistence the fit absorbs. Information only.
    pooled_ar <- function(G, k = 1) {
      if (nrow(G) <= k + 1) return(NA_real_)
      a <- G[-seq_len(k), , drop = FALSE]; b <- G[seq_len(nrow(G) - k), , drop = FALSE]
      if (sum(b^2) == 0) NA_real_ else sum(a * b) / sum(b^2)
    }
    acf_in <- function(g) { g <- g - mean(g); if (sum(g^2) == 0) NA_real_ else sum(g[-1] * g[-length(g)]) / sum(g^2) }
    T0 <- sum(!post)
    Zac <- sweep(GP[post, , drop = FALSE], 2, pmax(rmspe_P, 0.25 * stats::median(rmspe_P)), "/")
    rho_hat <- pooled_ar(Zac)
    rho_int <- if (is.na(rho_hat)) c(NA_real_, NA_real_) else
      stats::quantile(with_seed(replicate(200, pooled_ar(Zac[, sample(ncol(Zac), replace = TRUE), drop = FALSE]))), c(0.05, 0.95), na.rm = TRUE)
    rho_in <- stats::median(c(acf_in(gT_all[!post]), apply(GP[!post, , drop = FALSE], 2, acf_in)), na.rm = TRUE)
    s_len <- if (length(time_cols) >= 2) length(unique(periods[[time_cols[length(time_cols)]]])) else NA_integer_
    rho_seas <- if (!is.na(s_len) && s_len %in% 2:12 && n_post >= s_len + 3) pooled_ar(Zac, s_len) else NA_real_
    rc <- if (is.na(rho_hat)) NA_real_ else min(max(rho_hat, 0), 0.99)
    half_life <- if (is.na(rc)) NA_real_ else if (rho_hat >= 1) Inf else if (rc == 0) 0 else log(0.5) / log(rc)   # Inf: gaps do not revert
    backdate <- if (is.na(rc) || rc >= 0.9) NA_real_ else if (rc == 0) 1 else max(1, ceiling(log(0.25) / log(rc)))   # not useful when shocks barely fade
    seas_excess <- if (is.na(rho_seas) || is.na(rc)) NA_real_ else rho_seas - rc^s_len   # seasonal persistence beyond what lag 1 implies
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
    comps <- if (!is.null(aug_main)) aug_main else if (requireNamespace("fect", quietly = TRUE))
      tryCatch({ say("diagnostics (fect)")
        r_d <- if (!is.na(r_use)) r_use else model_ife(drop_post_unit = treated_unit)$r   # reuse the number of factors already chosen
        model_ife(drop_post_unit = treated_unit, r = r_d, cv = FALSE) },
               error = function(e) NULL) else NULL
    if (!is.null(comps) && !is.null(comps$lambda) && comps$r >= 1) {
      L <- comps$lambda; Fm <- comps$factor
      hull_dist <- function(target, pool) {
        w <- CVXR::Variable(nrow(pool))
        pr <- CVXR::Problem(CVXR::Minimize(CVXR::sum_squares(target - t(pool) %*% w)), list(w >= 0, sum(w) == 1))
        CVXR::psolve(pr); sqrt(max(0, as.numeric(CVXR::value(pr))))
      }
      d_T <- hull_dist(L[treated_unit, ], L[dn, , drop = FALSE])
      d_D <- vapply(dn, function(j) hull_dist(L[j, ], L[setdiff(dn, j), , drop = FALSE]), numeric(1))
      share_donors_beyond <- mean(d_D >= d_T)
      wts <- ft$w[dn]; wts[is.na(wts)] <- 0
      delta <- L[treated_unit, ] - colSums(L[dn, , drop = FALSE] * wts)
      corr <- as.numeric(Fm %*% delta); drift <- corr - mean(corr[!post])
      rho <- apply(Fm, 2, function(f) unname(stats::coef(stats::lm(f[-1] ~ f[-length(f)]))[2]))
      diag$hull <- list(n_factors = comps$r, treated_distance = d_T, donor_distance_median = median(d_D),
                        share_donors_beyond = share_donors_beyond, flag = d_T > 0 && share_donors_beyond < 0.10)
      diag$extrapolation <- list(max_post_drift_pct = 100 * max(abs(drift[post])) / mean(abs(synth_T[post])),
                                 final_drift_pct = 100 * drift[post][sum(post)] / synth_T[post][sum(post)])
      diag$persistence <- list(factor_ar1 = rho, trending = max(rho) > 0.95)
    } else {
      diag$hull <- list(n_factors = if (is.null(comps)) NA else comps$r, flag = NA,
                        note = "no loadings (0 factors chosen, or fect unavailable)")
      diag$persistence <- list(trending = NA)
    }
    extreme <- isTRUE(diag$fit_rank$flag) || isTRUE(diag$hull$flag)
    trending <- isTRUE(diag$persistence$trending)
    diag$recommendation <- NULL
    diag$assessment <- if (!extreme) {
      if (estimator == "fect") "no sign that the treated unit lies beyond the donors; the simpler estimator = \"dsc\" with scale = \"rmspe\" was as valid and at least as powerful in this situation in simulation" else
        "no sign that the treated unit lies beyond the donors"
    } else if (trending && augment == "ridge") {
      "treated unit appears beyond the donors' range and the factors are trending: augment = \"ridge\" is in use, which restored coverage in this situation in simulation (89-95% at 90%) where the plain demeaned synthetic control, augment = \"ife\" and fect under-covered; some bias can remain, so read the estimate with caution"
    } else if (trending) {
      diag$recommendation <- list(estimator = "dsc", augment = "ridge")
      "treated unit appears beyond the donors' range and the factors are trending: re-run with augment = \"ridge\" (ridge-augmented demeaned synthetic control with placebo calibration), which restored coverage in this situation in simulation (89-95% at 90%); the plain demeaned synthetic control, augment = \"ife\" and fect under-covered here"
    } else if (augment == "ridge") {
      "treated unit appears beyond the donors' range (factors not trending): augment = \"ridge\" was valid (conservative) in this situation in simulation; estimator = \"fect\" with scale = \"rmspe\" was also valid and more powerful than the plain demeaned synthetic control"
    } else if (estimator == "dsc") {
      diag$recommendation <- list(estimator = "fect", scale = "rmspe")
      "treated unit appears beyond the donors' range (factors not trending): re-run with estimator = \"fect\" and scale = \"rmspe\" (fect with placebo calibration), which was valid and much more powerful than DSC in this situation in simulation"
    } else if (scale != "rmspe") {
      diag$recommendation <- list(scale = "rmspe")
      "treated unit appears beyond the donors' range (factors not trending): fect suits this situation, but re-run with scale = \"rmspe\"; the relative and raw scales under-covered here in simulation"
    } else "treated unit appears beyond the donors' range (factors not trending): fect with placebo calibration on the RMSPE scale is the recommended method for this situation"
    if (!is.null(balance) && any(balance$flag %in% TRUE))
      diag$assessment <- paste0(diag$assessment, "; covariate imbalance: ", paste(balance$covariate[balance$flag %in% TRUE], collapse = ", "),
                                " (the treated unit is further from its synthetic control on this covariate than 90% or more of the placebos are from theirs; this biases the estimate if the covariate's effect on the outcome changes over time; see $balance)")
    if (warn && extreme) warning("dsc_placebo diagnostics: ", diag$assessment, ". See $diagnostics.", call. = FALSE)
  }

  # ---- Effects table (same column names as dsc_symcal_boot) ----------------------------------
  gc <- cumsum(gT); sc <- cumsum(sT)
  col <- function(v) c(rep(NA_real_, sum(!post)), v)
  effects <- cbind(periods, data.frame(
    post = post, observed = Y[, treated_unit], synthetic = synth_T, effect = gT_all,
    covariate_adjustment = Y[, treated_unit] - YaT[, treated_unit],
    rel_effect = col(gT / sT), std_effect = col(rT), factor_pointwise = col(cp$pw), factor_simultaneous = col(cp$sim),
    threshold_ratio_pw = col(abs(rT) / cp$pw), threshold_ratio_sim = col(abs(rT) / cp$sim),
    eff_low_pw = col(gT - cp$pw * sc_T), eff_high_pw = col(gT + cp$pw * sc_T),
    eff_low_sim = col(gT - cp$sim * sc_T), eff_high_sim = col(gT + cp$sim * sc_T),
    pred_low_pw = col(sT - cp$pw * sc_T), pred_high_pw = col(sT + cp$pw * sc_T),
    pred_low_sim = col(sT - cp$sim * sc_T), pred_high_sim = col(sT + cp$sim * sc_T),
    p_pointwise = col(cp$p_pw),
    cum_effect = col(gc), cum_rel_effect = col(gc / sc), cum_std_effect = col(rTc),
    cum_factor_pointwise = col(cc$pw), cum_factor_simultaneous = col(cc$sim),
    cum_threshold_ratio_pw = col(abs(rTc) / cc$pw), cum_threshold_ratio_sim = col(abs(rTc) / cc$sim),
    cum_eff_low_pw = col(gc - cc$pw * sc_Tc), cum_eff_high_pw = col(gc + cc$pw * sc_Tc),
    cum_eff_low_sim = col(gc - cc$sim * sc_Tc), cum_eff_high_sim = col(gc + cc$sim * sc_Tc),
    cum_pred_low_pw = col(sc - cc$pw * sc_Tc), cum_pred_high_pw = col(sc + cc$pw * sc_Tc),
    cum_pred_low_sim = col(sc - cc$sim * sc_Tc), cum_pred_high_sim = col(sc + cc$sim * sc_Tc),
    cum_p_pointwise = col(cc$p_pw)))

  # ---- Tests -------------------------------------------------------------------------------------
  e_post <- effects[post, ]
  test_label <- paste(names(periods), unlist(periods[t_treat + t_test - 1, ]), sep = " = ", collapse = ", ")
  flagged <- function(lo, hi) {
    hit <- which(lo > 0 | hi < 0)
    if (length(hit)) paste(apply(periods[post, , drop = FALSE][hit, , drop = FALSE], 1, paste, collapse = "-"), collapse = "; ") else "none"
  }
  ratio_T <- sqrt(mean(gT^2)) / rmspe_T; ratio_P <- sqrt(colMeans(GP[post, , drop = FALSE]^2)) / rmspe_P
  p_ratio <- (1 + sum(ratio_P >= ratio_T)) / (N + 1)
  test <- data.frame(
    test = c("pointwise", "pointwise", "simultaneous", "simultaneous"),
    measure = c("effect", "cumulative effect", "effect", "cumulative effect"),
    period = c(test_label, paste("through", test_label), "whole post period", "whole post period"),
    estimate = c(e_post$effect[t_test], e_post$cum_effect[t_test], NA, NA),
    lower = c(e_post$eff_low_pw[t_test], e_post$cum_eff_low_pw[t_test], NA, NA),
    upper = c(e_post$eff_high_pw[t_test], e_post$cum_eff_high_pw[t_test], NA, NA),
    p_value = c(e_post$p_pointwise[t_test], e_post$cum_p_pointwise[t_test], cp$p_sim, cc$p_sim),
    significant = c(e_post$eff_low_pw[t_test] > 0 | e_post$eff_high_pw[t_test] < 0,
                    e_post$cum_eff_low_pw[t_test] > 0 | e_post$cum_eff_high_pw[t_test] < 0, cp$p_sim <= alpha, cc$p_sim <= alpha),
    periods_flagged = c(NA, NA, flagged(e_post$eff_low_sim, e_post$eff_high_sim), flagged(e_post$cum_eff_low_sim, e_post$cum_eff_high_sim)))
  test <- rbind(test, data.frame(test = "RMSPE ratio (Abadie et al.)", measure = "post/pre RMSPE ratio", period = "whole post period",
                                 estimate = ratio_T, lower = NA, upper = NA, p_value = p_ratio, significant = p_ratio <= alpha,
                                 periods_flagged = NA))

  structure(list(effects = effects, test = test,
       weights = if (estimator == "fect") NULL else sort(ft$w[ft$w > 0], decreasing = TRUE),
       intercept = if (estimator == "fect") NULL else ft$b0,
       pre_rmspe = sqrt(mean(gT_all[!post]^2)), rmspe_guard = list(floor = floor_val, units_floored = floored_units),
       adjustment = if (is.null(adj_main)) NULL else {
         gu <- (Y[, treated_unit] - fit_dsc(Y[, donors, drop = FALSE], Y[, treated_unit])$synth)[post]
         list(model = if (estimator == "fect") "fect_ife" else if (augment == "ife") "ife_augment" else adjust_model,
              gamma = adj_main$gamma, n_factors = r_use,
              gamma_twfe = tw_main$gamma, se_twfe = tw_main$se, n_obs = adj_main$n_obs,
              unadjusted = c(effect = gu[t_test], cumulative = sum(gu[seq_len(t_test)])))
       },
       ridge = ridge_info,
       placebo = list(gaps = as.data.frame(GP), synthetic = as.data.frame(SP)),
       diagnostics = diag,
       balance = balance,
       failed_donors = failed_donors,
       settings = list(treated_unit = treated_unit, treatment_period = treatment_period, test_period = test_label,
                       time_cols = time_cols, outcome = outcome_col, conf_level = conf_level, span = span,
                       scale = scale, rmspe_floor = if (scale == "rmspe") rmspe_floor else NA, adjust_covariates = adjust_covariates,
                       balance_covariates = balance_covariates, adjust_model = if (length(adjust_covariates)) adjust_model else NA, augment = augment, estimator = estimator,
                       ridge_c = if (augment == "ridge") c_ridge else NA, n_donors = N, rank_used = m)), class = "dsc_placebo")
}
