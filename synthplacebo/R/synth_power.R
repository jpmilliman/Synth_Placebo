# synth_power.R
#
# synth_power(): how large an effect can the placebo test and the LTO test detect on a given panel?
#
# The real treated unit is dropped. Each remaining unit is made the (pseudo-)treated unit in turn, a known
# effect is added to its post-treatment outcomes, and the tests are run with the other units as donors. Power
# is the share of pseudo-treated units for which a test rejects; with no effect added, that share is the
# exact false-positive rate of the design (as in exact_size_validation.R). Every fit uses pre-treatment
# data only, so the effect is added without refitting (lto_shift()): each unit is fitted once, and any
# number of effect sizes and shapes cost almost nothing more.
#
# Effects (effects x effect_scale x shapes):
#   effect_scale  "percent" (share of the unit's own untreated outcome in each post period: -0.1 = 10% lower),
#                 "level" (outcome units), or "rmspe" (multiples of the unit's pre-treatment fit error: the
#                 effect relative to the noise, which is what power depends on)
#   shapes        "constant" (a step at treatment that lingers), "ramp" (rising linearly to the full size in the
#                 last post period; ramp_jump > 0 starts it with an immediate jump of that share of the full size),
#                 "plateau" (phase-in then plateau: rising linearly to the full size over plateau_periods post
#                 periods, then flat; ramp_jump applies too), "decay" (the full size in the first post period,
#                 then fading geometrically: decay_rate is the share kept each period, so 0.5 halves it each
#                 period), "shock" (one period only, at shock_period; default the first post period)
# Tests (tests x methods):
#   "pointwise" (the test period; for "shock" the shock period and for "decay" the first post period, where
#   the effect is largest), "cumulative" (through the test period),
#   "simultaneous" (the whole post period) - for methods "placebo" (dsc_placebo()'s rank tests, with the
#   treated unit's gap on the RMSPE scale; the smallest possible p-value is 1 / (number of donors + 1)) and
#   "lto" (the naive LTO test of dsc_placebo_lto()).
# Sidedness (alternative): "two.sided" (the default, as in the two functions) and/or "one.sided". A one-sided
#   test ranks the signed gap in the direction of the effect (a test of "less" for effects below zero, of
#   "greater" for effects above), as if that direction had been stated in advance (lto_pvalues_onesided()).
#   With no effect both directions are tested, so each has its own false-positive rate.
# Covariates (optional): adjust_covariates as in the two functions. Each pseudo-treated unit's covariate effect
#   gamma is estimated once, on every cell except that unit's post period (two-way or interactive fixed
#   effects, fect), and all units' outcomes are adjusted by it. This approximates the functions, which
#   re-estimate gamma for every placebo (dsc_placebo()) or every unit (dsc_placebo_lto()), at 1/N of the cost.
#   confounder = list(covariate, shift, gamma): a confounding shock added to each pseudo-treated unit after
#   treatment, as in simulations/oregon_fentanyl_validation.R: its covariate rises by `shift` (a number, or one
#   value per post period) and its outcome by gamma x shift. Because gamma-hat is estimated without the unit's
#   post period, the shock moves the unit's adjusted gap by exactly (gamma - gamma-hat) x shift (unadjusted:
#   gamma x shift), so it is added like the effect. With effects = 0 this gives the false-positive rate under
#   confounding.
# Run time: the fits are those of one dsc_placebo_lto() run per unit, shared across units (seconds to a few
#   minutes without covariates; with adjust_model = "ife", one fect fit per unit, about 10 s each for Kansas).
#   methods = "placebo" skips the tournament. n_cores > 1 spreads the units over that many cores.
#
# Returns: power (one row per method, test, shape, effect and alpha: power and the number of units), mde (the
#   smallest effect in the grid with power >= target_power, by direction, and a linear interpolation between
#   grid points), p_values (one row per unit and scenario), settings. plot_synth_power() draws the power curves.
#
# Needs dsc_placebo_lto.R (Experimental_inference/LTO/functions/; source it first) and quadprog; fect for
# adjust_model = "ife"; ggplot2 for plot_synth_power(). Documentation: Experimental_inference/synth_power/README.md.
# Written by Claude (Anthropic), at the direction of Jeffrey Milliman.

#' Power and minimum detectable effects of the placebo and LTO tests on a panel
#'
#' How large an effect could [dsc_placebo()]'s placebo tests and the leave-two-out test of [dsc_placebo_lto()]
#' detect on the panel at hand? The real treated unit is dropped and every remaining unit is made the treated
#' unit in turn, with the others as donors. A known effect is added to its post-treatment outcomes and the tests
#' are run. Power is the share of units for which a test rejects; with no effect added, that share is the
#' design's exact false-positive rate. Every fit uses pre-treatment data only, so effects are added without
#' refitting and a fine grid of sizes and shapes costs little.
#'
#' @details
#' Power on these panels is limited less by the pre-treatment fit than by how far placebo gaps drift after
#' treatment; see the `power_across_datasets` document in the repository. A null result needs a power
#' statement: running this on the study's own panel shows which effects the data could have detected.
#'
#' **Covariates.** With `adjust_covariates`, \eqn{\gamma} is estimated once per pseudo-treated unit (an
#' approximation of the functions, which re-estimate it per placebo or per unit). `confounder` adds a
#' confounding shock: after treatment the pseudo-treated unit's covariate rises by `shift` and its outcome by
#' `gamma` times that, which moves its adjusted gap by exactly \eqn{(\gamma - \hat\gamma)} times the shift.
#'
#' @inheritParams dsc_placebo
#' @param treated_unit The real treated unit, dropped before the analysis (`NULL` keeps every unit).
#' @param effects Effect sizes; 0 is always added (the false-positive rate).
#' @param effect_scale `"percent"` (share of the unit's own untreated outcome: -0.1 is 10% lower), `"level"`
#'   (outcome units) or `"rmspe"` (multiples of the unit's pre-treatment RMSPE).
#' @param shapes Any of `"constant"` (a step that lingers), `"ramp"` (linear build-up to the full size in the
#'   last period), `"plateau"` (build-up over `plateau_periods`, then flat), `"decay"` (full size at once, then
#'   fading by `decay_rate` per period) and `"shock"` (one period only).
#' @param ramp_jump Immediate share of the full effect for `"ramp"` and `"plateau"` (0: no jump).
#' @param plateau_periods Post periods over which `"plateau"` reaches its full size.
#' @param decay_rate Share of the effect kept each period for `"decay"`.
#' @param methods `"placebo"` and/or `"lto"`; `"placebo"` alone skips the tournament.
#' @param tests `"pointwise"` (single period: `test_period`, or the peak for decay and shock), `"cumulative"`
#'   and/or `"simultaneous"`.
#' @param alternative `"two.sided"` and/or `"one.sided"` (in the effect's direction, as if stated in advance).
#' @param alpha Level(s) of the tests.
#' @param test_period Period of the single-period and cumulative tests (default: the last post period).
#' @param shock_period Period of the `"shock"` (default: the first post period).
#' @param target_power Power at which the minimum detectable effect is reported.
#' @param augment,ridge_c `"none"` or `"ridge"`, and the ridge penalty, as in [dsc_placebo_lto()].
#' @param adjust_covariates,adjust_model,ife_r,ife_seed Covariate adjustment as in [dsc_placebo_lto()].
#' @param confounder Optional `list(covariate, shift, gamma)`: a confounding shock after treatment (`shift`: one
#'   number or one per post period).
#' @param n_cores Spread the units over this many cores.
#' @return An object of class `"synth_power"`: a list with `power` (by `alpha`, method, test, shape,
#'   sidedness, direction and effect), `mde` (the smallest effect in the grid with power at or above
#'   `target_power`, interpolated linearly; `NA` when not reached, with `max_power_tested`), `p_values` (one row
#'   per unit and scenario) and `settings`.
#' @seealso [plot_synth_power()], [print.synth_power()].
#' @examples
#' \donttest{
#' if (requireNamespace("scpi", quietly = TRUE)) {
#'   data("scpi_germany", package = "scpi")
#'   germany <- as.data.frame(scpi_germany)[, c("country", "year", "gdp")]
#'   pw <- synth_power(germany, "country", "year", "gdp", c(year = 1991), treated_unit = "West Germany",
#'                     shapes = c("constant", "ramp"), verbose = FALSE)
#'   pw
#'   plot(pw)
#' }
#' }
#' @export
synth_power <- function(data, unit_col, time_cols, outcome_col, treatment_period, treated_unit = NULL,
                        effects = c(-0.5, -0.4, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5),
                        effect_scale = c("percent", "level", "rmspe"), shapes = c("constant", "ramp", "plateau", "decay", "shock"),
                        ramp_jump = 0, plateau_periods = 3, decay_rate = 0.5, methods = c("placebo", "lto"), tests = c("pointwise", "cumulative", "simultaneous"),
                        alternative = "two.sided", alpha = 0.10, test_period = NULL, shock_period = NULL, target_power = 0.8,
                        rmspe_floor = 0.25, span = 0.75, augment = c("none", "ridge"), ridge_c = 0.1,
                        adjust_covariates = NULL, adjust_model = c("ife", "twfe"), ife_r = 0:5, ife_seed = 12112025,
                        confounder = NULL, n_cores = 1, verbose = interactive()) {
  effect_scale <- match.arg(effect_scale); augment <- match.arg(augment); adjust_model <- match.arg(adjust_model)
  shapes <- match.arg(shapes, c("constant", "ramp", "plateau", "decay", "shock"), several.ok = TRUE)
  if (!is.numeric(ramp_jump) || length(ramp_jump) != 1 || ramp_jump < 0 || ramp_jump > 1) stop("ramp_jump must be a number in [0, 1].")
  if (!is.numeric(decay_rate) || length(decay_rate) != 1 || decay_rate < 0 || decay_rate >= 1) stop("decay_rate must be a number in [0, 1).")
  if (!is.numeric(plateau_periods) || length(plateau_periods) != 1 || plateau_periods < 1) stop("plateau_periods must be a number >= 1.")
  methods <- match.arg(methods, c("placebo", "lto"), several.ok = TRUE)
  tests <- match.arg(tests, c("pointwise", "cumulative", "simultaneous"), several.ok = TRUE)
  alternative <- match.arg(alternative, c("two.sided", "one.sided"), several.ok = TRUE)
  if (!is.numeric(effects) || !length(effects)) stop("effects must be a numeric vector.")
  if (!is.numeric(alpha) || any(alpha <= 0 | alpha >= 1)) stop("alpha must be in (0, 1).")
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) stop("n_cores must be a single number >= 1.")
  effects <- sort(unique(c(0, effects)))                 # 0: the false-positive rate
  say <- function(...) if (isTRUE(verbose)) message("synth_power: ", ...)

  # ---- data: units x periods, the real treated unit dropped ---------------------------------------
  d <- as.data.frame(data)
  need <- c(unit_col, time_cols, outcome_col, adjust_covariates, if (!is.null(confounder)) confounder$covariate)
  if (!all(need %in% names(d))) stop("columns not found: ", paste(setdiff(need, names(d)), collapse = ", "))
  if (!is.null(treated_unit)) d <- d[as.character(d[[unit_col]]) != treated_unit, , drop = FALSE]
  d <- d[do.call(order, d[time_cols]), , drop = FALSE]
  key <- do.call(paste, c(d[time_cols], sep = "_"))
  wide <- function(v) { w <- stats::reshape(data.frame(key = key, unit = as.character(d[[unit_col]]), z = d[[v]]), idvar = "key", timevar = "unit", direction = "wide")
    m <- as.matrix(w[match(unique(key), w$key), -1, drop = FALSE]); colnames(m) <- sub("^z\\.", "", colnames(m)); m }
  Y <- wide(outcome_col); if (anyNA(Y)) stop("The panel must be balanced with no missing outcomes.")
  N <- ncol(Y); if (N < 4) stop("At least 4 units are needed (3 donors per pseudo-treated unit).")
  periods <- unique(d[time_cols]); rownames(periods) <- NULL
  after <- function(tp) vapply(seq_len(nrow(periods)), function(r) {
    for (v in names(tp)) if (periods[[v]][r] != tp[[v]]) return(periods[[v]][r] > tp[[v]]); TRUE }, logical(1))
  post <- after(treatment_period); T1 <- sum(post); if (T1 < 1 || sum(!post) < 2) stop("treatment_period leaves no post or too few pre periods.")
  row_of <- function(tp, what) { if (is.null(tp)) return(NULL)
    hit <- which(Reduce(`&`, lapply(names(tp), function(v) periods[[v]] == tp[[v]])))
    if (length(hit) != 1 || !post[hit]) stop(what, " must be one post-treatment period."); hit - sum(!post) }
  h_test <- if (is.null(test_period)) T1 else row_of(test_period, "test_period")
  h_shock <- if (is.null(shock_period)) 1L else row_of(shock_period, "shock_period")
  Zs <- stats::setNames(lapply(adjust_covariates, function(v) wide(v)[, colnames(Y), drop = FALSE]), adjust_covariates)
  if (length(Zs) && any(vapply(Zs, anyNA, logical(1)))) stop("adjust_covariates must have no missing values.")
  adjust <- length(Zs) > 0

  # ---- the confounding shock (a post-period path) ----------------------------------------------------
  conf_path <- NULL
  if (!is.null(confounder)) {
    if (!all(c("covariate", "shift", "gamma") %in% names(confounder))) stop("confounder needs covariate, shift and gamma.")
    conf_path <- rep_len(as.numeric(confounder$shift), T1)
    if (length(confounder$shift) != 1 && length(confounder$shift) != T1) stop("confounder$shift must be one number or one per post period.")
  }

  # ---- ridge penalty, number of factors ---------------------------------------------------------------
  c_r <- if (augment == "ridge") { if (identical(ridge_c, "cv") || is.null(ridge_c)) lto_choose_ridge(Y, !post) else ridge_c } else NULL
  excl_for <- function(u) { e <- matrix(FALSE, nrow(Y), N); e[post, u] <- TRUE; e }
  r_use <- NA_integer_
  if (adjust && adjust_model == "ife") {
    if (!requireNamespace("fect", quietly = TRUE)) stop("adjust_model = \"ife\" needs the fect package.")
    r_use <- if (length(unique(ife_r)) == 1) unique(ife_r) else {
      say("choosing the number of factors by cross-validation (fect, once)")
      as.integer(lto_fect(Y, Zs, excl_for(1), ife_r, TRUE, ife_seed)$r.cv) }
  }
  gamma_for <- function(u) {                             # gamma-hat for pseudo-treated unit u (its post period left out)
    if (adjust_model == "twfe") lto_gamma_twfe(Y, Zs, excl_for(u)) else {
      fit <- lto_fect(Y, Zs, excl_for(u), r_use, FALSE, ife_seed); stats::setNames(as.numeric(fit$beta[names(Zs), 1]), names(Zs)) }
  }
  tt <- seq_len(T1)
  shape_vec <- function(sh) switch(sh, constant = rep(1, T1),
    ramp = ramp_jump + (1 - ramp_jump) * tt / T1,
    plateau = ramp_jump + (1 - ramp_jump) * pmin(tt / plateau_periods, 1),
    decay = decay_rate^(tt - 1),
    shock = as.numeric(tt == h_shock))
  lto_sim_p <- function(core) {                          # the naive LTO simultaneous p-value (as lto_simultaneous())
    q <- lto_profile(core, span)$q
    Mi <- pmax(apply(core$Ki / q, 2, max), apply(core$Kj / q, 2, max))
    mean(apply(abs(core$G) / outer(q, core$S), 2, max) <= Mi)
  }

  # ---- one pseudo-treated unit: fit once, then every scenario ----------------------------------------
  unit_task <- function(u, cache) {
    g <- if (adjust) gamma_for(u) else NULL
    Ya <- Y; if (adjust) for (v in names(Zs)) Ya <- Ya - g[[v]] * Zs[[v]]
    core0 <- lto_core(Ya, u, post, rmspe_floor, cache = if (adjust) NULL else cache, c_ridge = c_r, tournament = "lto" %in% methods)
    y0 <- Y[post, u]
    base <- switch(effect_scale, percent = y0, level = rep(1, T1), rmspe = rep(core0$standard$s, T1))
    conf <- if (is.null(conf_path)) 0 else
      (confounder$gamma - if (adjust && confounder$covariate %in% names(Zs)) g[[confounder$covariate]] else 0) * conf_path
    out <- list()
    for (sh in shapes) {
      sv <- shape_vec(sh); hp <- switch(sh, shock = h_shock, decay = 1L, h_test)
      for (e in effects) {
        core <- lto_shift(core0, e * base * sv + conf); pv <- lto_pvalues(core, 0)
        add <- function(method, test, p, sided = "two.sided", direction = "both") out[[length(out) + 1]] <<-
          data.frame(unit = colnames(Y)[u], shape = sh, effect = e, method = method, test = test, sided = sided, direction = direction, p = p)
        if ("one.sided" %in% alternative) for (dr in if (e < 0) "less" else if (e > 0) "greater" else c("less", "greater")) {
          os <- lto_pvalues_onesided(core, dr, span)
          for (m in methods) {
            if ("pointwise" %in% tests) add(m, "pointwise", os[[paste0("p_", m)]][hp], "one.sided", dr)
            if ("cumulative" %in% tests) add(m, "cumulative", os[[paste0("p_", m, "_cum")]][h_test], "one.sided", dr)
            if ("simultaneous" %in% tests) add(m, "simultaneous", os[[paste0("p_", m, "_sim")]], "one.sided", dr)
          }
        }
        if (!"two.sided" %in% alternative) next
        if ("placebo" %in% methods) {
          if ("pointwise" %in% tests) add("placebo", "pointwise", pv$p_placebo[hp])
          if ("cumulative" %in% tests) add("placebo", "cumulative", pv$p_placebo_cum[h_test])
          if ("simultaneous" %in% tests) add("placebo", "simultaneous", standard_simultaneous(core, 0.9, span)$p_sim)
        }
        if ("lto" %in% methods) {
          if ("pointwise" %in% tests) add("lto", "pointwise", pv$p_lto[hp])
          if ("cumulative" %in% tests) add("lto", "cumulative", pv$p_lto_cum[h_test])
          if ("simultaneous" %in% tests) add("lto", "simultaneous", lto_sim_p(core))
        }
      }
    }
    out <- do.call(rbind, out)
    if (adjust) out$gamma_hat <- paste(signif(g, 4), collapse = ", ")
    out
  }

  say(N, " pseudo-treated units, ", length(effects) - 1, " effect sizes x ", length(shapes), " shapes", if (adjust) paste0("; ", N, " covariate models") else "",
      if (n_cores > 1) paste0(", on ", n_cores, " cores") else "")
  t0 <- Sys.time()
  shared <- new.env()
  chunks <- if (n_cores > 1) split(seq_len(N), cut(seq_len(N), min(N, 2 * n_cores), labels = FALSE)) else as.list(seq_len(N))
  res <- lto_par_lapply(chunks, function(us) { cache <- if (n_cores > 1) new.env() else shared
    do.call(rbind, lapply(us, function(u) unit_task(u, cache))) }, n_cores = n_cores, progress = isTRUE(verbose))
  p_values <- do.call(rbind, res); rownames(p_values) <- NULL

  # ---- power and minimum detectable effects -------------------------------------------------------------
  power <- do.call(rbind, lapply(alpha, function(a) {
    agg <- stats::aggregate(list(power = p_values$p <= a + 1e-12), p_values[c("method", "test", "shape", "sided", "direction", "effect")], mean)
    agg$alpha <- a; agg$n_units <- N; agg }))
  power <- power[order(power$alpha, power$method, power$test, power$shape, power$sided, power$effect, power$direction), ]; rownames(power) <- NULL
  mde <- do.call(rbind, lapply(split(power, power[c("alpha", "method", "test", "shape", "sided")], drop = TRUE), function(pw) {
    do.call(rbind, lapply(c(-1, 1), function(sg) {
      s <- pw[sign(pw$effect) == sg, ]; s <- s[order(abs(s$effect)), ]
      z <- pw[pw$effect == 0, ]; fp <- if (pw$sided[1] == "two.sided") z$power[1] else z$power[z$direction == if (sg < 0) "less" else "greater"]
      hit <- which(s$power >= target_power)[1]
      interp <- if (is.na(hit)) NA_real_ else {
        x0 <- if (hit == 1) 0 else abs(s$effect[hit - 1]); y0 <- if (hit == 1) fp else s$power[hit - 1]
        sg * (x0 + (target_power - y0) / (s$power[hit] - y0) * (abs(s$effect[hit]) - x0)) }
      data.frame(alpha = pw$alpha[1], method = pw$method[1], test = pw$test[1], shape = pw$shape[1], sided = pw$sided[1], direction = if (sg < 0) "decrease" else "increase",
                 false_positive = fp, mde_grid = if (is.na(hit)) NA_real_ else s$effect[hit], mde_interpolated = interp,
                 max_power_tested = if (nrow(s)) max(s$power) else NA_real_)
    }))
  }))
  rownames(mde) <- NULL
  say("done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " minutes")
  structure(list(power = power, mde = mde, p_values = p_values,
       settings = list(units = colnames(Y), dropped = treated_unit, n_units = N, n_donors = N - 1, effect_scale = effect_scale,
                       shapes = shapes, ramp_jump = ramp_jump, plateau_periods = plateau_periods, decay_rate = decay_rate, methods = methods, tests = tests, alternative = alternative, alpha = alpha, target_power = target_power,
                       test_period_row = h_test, shock_period_row = h_shock, n_post = T1, augment = augment, ridge_c = c_r,
                       adjust_covariates = adjust_covariates, adjust_model = if (adjust) adjust_model else NA, n_factors = r_use,
                       confounder = confounder, smallest_placebo_p = 1 / N)), class = "synth_power")
}

# Power curves: power against effect size; rows are effect shapes, columns tests, colour the method.
# shapes / tests / sided: subsets to draw (default all in res; sided "two.sided" or "one.sided").
#' Power curves from synth_power()
#'
#' Power against effect size: one row per effect shape, one column per test, coloured by method (line type by
#' sidedness when both are present).
#'
#' @param res Output of [synth_power()].
#' @param alpha The level to draw (default: the first in `res`).
#' @param title Plot title.
#' @param shapes,tests,sided Subsets to draw (`sided`: `"two.sided"` or `"one.sided"`).
#' @return A `ggplot` object.
#' @export
plot_synth_power <- function(res, alpha = res$settings$alpha[1], title = NULL, shapes = NULL, tests = NULL, sided = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("plot_synth_power() needs ggplot2.")
  pw <- res$power[abs(res$power$alpha - alpha) < 1e-12, ]
  if (!is.null(shapes)) pw <- pw[pw$shape %in% shapes, ]
  if (!is.null(tests)) pw <- pw[pw$test %in% tests, ]
  if (!is.null(sided)) pw <- pw[pw$sided %in% sided, ]
  pw$test <- factor(pw$test, c("pointwise", "cumulative", "simultaneous"), c("Single period", "Cumulative", "Simultaneous"))
  pw$shape <- factor(pw$shape, intersect(c("constant", "ramp", "plateau", "decay", "shock"), pw$shape))
  pw$line <- interaction(pw$method, pw$sided, pw$direction, drop = TRUE)   # one-sided: a curve for each direction, each starting at its own rate at 0
  xl <- switch(res$settings$effect_scale, percent = "Effect (share of the unit's own outcome)", level = "Effect (outcome units)",
               rmspe = "Effect (multiples of the unit's pre-treatment RMSPE)")
  two <- length(unique(pw$sided)) > 1
  map <- if (two) ggplot2::aes(effect, power, colour = method, group = line, linetype = sided) else ggplot2::aes(effect, power, colour = method, group = line)
  ggplot2::ggplot(pw, map) +
    ggplot2::geom_hline(yintercept = alpha, colour = "grey60", linetype = "dotted") +
    ggplot2::geom_hline(yintercept = res$settings$target_power, colour = "grey60", linetype = "dashed") +
    ggplot2::geom_line() + ggplot2::geom_point(size = 1) +
    ggplot2::facet_grid(shape ~ test) + ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    ggplot2::labs(x = xl, y = paste0("Power (share of ", res$settings$n_units, " units rejecting at ", 100 * alpha, "%)"),
                  title = if (is.null(title)) "Power by effect size" else title,
                  subtitle = paste("Single period: the test period (the peak for decay and shock). Dotted: the level.",
                                   "Dashed: the target power.", sep = "\n")) +
    ggplot2::theme_minimal()
}
