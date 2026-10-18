##############################################################################
# synthetic_control_functions_claude.R
#
# Synthetic control for a single treated unit with placebo-calibrated inference.
# Written by Claude (Anthropic) at the direction of the repository owner.
# Full argument and output reference: docs/reference/function_reference.md
#
# Main method:
#   6. dsc_placebo()       - demeaned synthetic control (intercept + non-negative donor
#                            weights summing to one) with placebo-calibrated pointwise
#                            intervals, simultaneous bands and prediction intervals;
#                            covariate adjustment, an optional fect estimator, and
#                            diagnostics
#   7. plot_full_path(), plot_post_effect(), plot_prediction(), added_effect()
#                          - ggplot figures for dsc_placebo() / dsc_symcal_boot() output
#
# Earlier versions, kept for the record of the method's development:
#   1. meboot_multi()      - maximum entropy bootstrap of multiple columns (used inside
#                            synth_lasso_boot())
#   2. synth_lasso_boot()  - LASSO synthetic control with bootstrap inference
#   3. synth_lasso_placebo() - the same model with exact placebo (permutation) tests
#   4. dsc_symcal_boot()   - demeaned synthetic control with a symmetric, donor-
#                            calibrated fixed-weight bootstrap
#   5. plot_dsc_symcal()   - ggplot figures for dsc_symcal_boot() results
#
# Usage (from the project root):
#   source("synthetic_control_functions/synthetic_control_functions_claude.R")
##############################################################################

## ===========================================================================
## 1. meboot_multi()
## ===========================================================================

##############################################################################
# meboot_multi(): Maximum entropy bootstrap for multiple columns
#
# meboot::meboot() only bootstraps one series at a time. This function loops
# meboot over several columns, lets you set trim / limits (xmin, xmax) per
# column, and returns a list of `reps` bootstrapped data frames in either
# long or wide format.
#
# Arguments
#   data        data frame (long format, e.g. the `kansas` panel)
#   boot_cols   character vector of numeric columns to bootstrap
#   id_cols     character vector of columns carried along unchanged
#               (e.g. c("state", "year", "qtr")). Default: all columns not
#               in boot_cols.
#   reps        number of bootstrap replicates (passed to meboot)
#   trim        EITHER a single list applied to every column, e.g.
#                 list(trim = 0.10, xmin = 0, xmax = NULL)
#               OR a named list with one list per column, e.g.
#                 list(gdpcapita   = list(trim = 0.10, xmin = 0),
#                      lngdpcapita = list(trim = 0.10, xmin = NULL))
#               Settings within each list:
#                 trim      trimmed-mean proportion (default 0.10)
#                 xmin/xmax absolute lower / upper tail limits. One number
#                           for every series, or - with group_col - a
#                           vector with at most one unnamed default and
#                           values named by group, e.g.
#                           xmin = c(15000, Alaska = 30000)
#                 xmin_rel  lower limit as a fraction of EACH series' own
#                           minimum, e.g. 0.9 (positive data only)
#                 xmax_rel  upper limit as a multiple of each series' own
#                           maximum, e.g. 1.1 (positive data only)
#               Give xmin or xmin_rel (xmax or xmax_rel), not both. Unset
#               limits use the meboot defaults; columns not listed get the
#               defaults. A limit inside the observed range is an error.
#   group_col   optional column name. NULL (default) bootstraps each column
#               as one stacked series, exactly like the Kansas example.
#               If set (e.g. "state"), each group's series is bootstrapped
#               separately. Rows are used in their current order within each
#               group, so arrange the data by time first.
#   output      "long" (default) or "wide"
#   names_from  column to pivot into new columns when output = "wide"
#               (e.g. "state"). Required for wide output.
#   names_sep   separator for wide column names when more than one column
#               is bootstrapped (e.g. "gdpcapita_Kansas")
#   seed        optional seed, set once before any bootstrapping
#   shared_draws  with group_col only: TRUE gives every group the same
#               random numbers within each column and replicate, so groups
#               are nudged up or down together (by rank). FALSE (default)
#               draws independently for every group.
#   as_matrix   if TRUE, each list element is converted with as.matrix()
#               (useful for glmnet). All columns must be numeric for a
#               numeric matrix.
#   ...         further arguments passed to meboot::meboot()
#               (e.g. reachbnd, expand.sd, force.clt, sym)
#
# Value
#   A list of length `reps` (named rep_1, rep_2, ...). Each element is a
#   data frame (or matrix) holding id_cols plus the bootstrapped columns.
##############################################################################

meboot_multi <- function(data,
                         boot_cols,
                         id_cols = NULL,
                         reps = 999,
                         trim = list(trim = 0.10, xmin = NULL, xmax = NULL),
                         group_col = NULL,
                         output = c("long", "wide"),
                         names_from = NULL,
                         names_sep = "_",
                         seed = NULL,
                         as_matrix = FALSE,
                         shared_draws = FALSE,
                         ...) {
  
  output <- match.arg(output)
  data <- as.data.frame(data)
  
  ## ---- Checks ------------------------------------------------------------
  missing_cols <- setdiff(c(boot_cols, id_cols, group_col, names_from), names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }
  
  not_numeric <- boot_cols[!vapply(data[boot_cols], is.numeric, logical(1))]
  if (length(not_numeric) > 0) {
    stop("boot_cols must be numeric: ", paste(not_numeric, collapse = ", "))
  }
  
  has_na <- boot_cols[vapply(data[boot_cols], anyNA, logical(1))]
  if (length(has_na) > 0) {
    stop("meboot cannot handle missing values. NAs found in: ",
         paste(has_na, collapse = ", "))
  }
  
  if (output == "wide" && is.null(names_from)) {
    stop("output = 'wide' requires names_from (e.g. names_from = 'state').")
  }
  
  # Default id columns: everything that is not being bootstrapped
  if (is.null(id_cols)) {
    id_cols <- setdiff(names(data), boot_cols)
  }
  
  # Grouping and pivoting columns must be kept with the output
  id_cols <- unique(c(id_cols, group_col, names_from))
  
  ## ---- Build a trim list for each column ---------------------------------
  # Keys: trim, xmin, xmax (absolute limits) and xmin_rel, xmax_rel (limits
  # relative to each series' own minimum / maximum).
  trim_default <- list(trim = 0.10, xmin = NULL, xmax = NULL,
                       xmin_rel = NULL, xmax_rel = NULL)
  trim_keys <- names(trim_default)
  
  # One list for all columns if every name is a trim key
  is_global_trim <- !is.null(names(trim)) && all(names(trim) %in% trim_keys)
  
  if (is_global_trim) {
    trim_list <- rep(list(utils::modifyList(trim_default, trim, keep.null = TRUE)),
                     length(boot_cols))
    names(trim_list) <- boot_cols
  } else {
    bad_names <- setdiff(names(trim), boot_cols)
    if (is.null(names(trim)) || length(bad_names) > 0) {
      stop("trim must be a single list (", paste(trim_keys, collapse = "/"),
           ") or a list named by boot_cols. Unrecognised name(s): ",
           paste(bad_names, collapse = ", "))
    }
    trim_list <- lapply(boot_cols, function(col) {
      if (col %in% names(trim)) {
        bad_keys <- setdiff(names(trim[[col]]), trim_keys)
        if (length(bad_keys) > 0) {
          stop("trim for ", col, ": unknown setting(s) ", paste(bad_keys, collapse = ", "),
               ". Use ", paste(trim_keys, collapse = ", "), ".")
        }
        utils::modifyList(trim_default, trim[[col]], keep.null = TRUE)
      } else {
        trim_default
      }
    })
    names(trim_list) <- boot_cols
  }
  
  groups <- if (is.null(group_col)) NULL else unique(as.character(data[[group_col]]))
  
  # Check limits named by unit before any bootstrapping
  for (col in boot_cols) {
    for (key in c("xmin", "xmax")) {
      val <- trim_list[[col]][[key]]
      if (!is.null(trim_list[[col]][[paste0(key, "_rel")]]) && !is.null(val)) {
        stop("trim for ", col, ": give ", key, " or ", key, "_rel, not both.")
      }
      nm <- names(val)
      if (!is.null(nm) && any(nm != "")) {
        if (is.null(group_col)) {
          stop("trim for ", col, ": limits named by unit (", key, ") need group_col ",
               "(boot_by_unit = TRUE in synth_lasso_boot()).")
        }
        if (sum(nm == "") > 1) {
          stop("trim for ", col, ": ", key, " can have at most one unnamed default value.")
        }
        unknown <- setdiff(nm[nm != ""], groups)
        if (length(unknown) > 0) {
          stop("trim for ", col, ": ", key, " names unit(s) not in the data: ",
               paste(unknown, collapse = ", "))
        }
      }
    }
  }
  
  # The meboot trim list for one series (one column, and one group if grouped)
  series_trim <- function(col, x, group = NULL) {
    spec <- trim_list[[col]]
    where <- if (is.null(group)) col else paste0(col, " (", group, ")")
    
    pick <- function(val) {
      if (is.null(val)) return(NULL)
      nm <- names(val)
      if (is.null(nm) || all(nm == "")) return(unname(val[1]))
      if (!is.null(group) && group %in% nm) return(unname(val[group]))
      if (any(nm == "")) return(unname(val[nm == ""][1]))
      NULL
    }
    
    xmin <- pick(spec$xmin)
    xmax <- pick(spec$xmax)
    
    if (!is.null(spec$xmin_rel)) {
      if (any(x <= 0)) stop("xmin_rel needs strictly positive values; ", where, " has values <= 0. Use xmin instead.")
      if (spec$xmin_rel > 1) stop("xmin_rel must be <= 1 (a fraction of the series minimum).")
      xmin <- spec$xmin_rel * min(x)
    }
    if (!is.null(spec$xmax_rel)) {
      if (any(x <= 0)) stop("xmax_rel needs strictly positive values; ", where, " has values <= 0. Use xmax instead.")
      if (spec$xmax_rel < 1) stop("xmax_rel must be >= 1 (a multiple of the series maximum).")
      xmax <- spec$xmax_rel * max(x)
    }
    
    if (!is.null(xmin) && xmin > min(x)) {
      stop("xmin for ", where, " (", signif(xmin, 6), ") is above its lowest value (",
           signif(min(x), 6), ").")
    }
    if (!is.null(xmax) && xmax < max(x)) {
      stop("xmax for ", where, " (", signif(xmax, 6), ") is below its highest value (",
           signif(max(x), 6), ").")
    }
    list(trim = spec$trim, xmin = xmin, xmax = xmax)
  }
  
  ## ---- Run meboot for each column (and group) ----------------------------
  if (shared_draws && is.null(group_col)) {
    warning("shared_draws only applies when bootstrapping by group (group_col); ignored.")
    shared_draws <- FALSE
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # With shared draws, every group of a column restarts from the same seed, so
  # meboot draws the same random numbers for each group in each replicate.
  col_seeds <- if (shared_draws) sample.int(.Machine$integer.max, length(boot_cols)) else NULL
  
  run_meboot <- function(x, col_trim) {
    ens <- meboot::meboot(x, reps = reps, trim = col_trim, ...)$ensemble
    as.matrix(ens)  # rows = observations, columns = replicates
  }
  
  ensembles <- lapply(seq_along(boot_cols), function(j) {
    col <- boot_cols[j]
    x <- data[[col]]
    
    if (is.null(group_col)) {
      return(run_meboot(x, series_trim(col, x)))
    }
    
    # Bootstrap each group separately and put results back in original rows
    ens <- matrix(NA_real_, nrow = length(x), ncol = reps)
    row_index <- split(seq_along(x), as.character(data[[group_col]]))
    for (g in names(row_index)) {
      rows <- row_index[[g]]
      if (shared_draws) set.seed(col_seeds[j])
      ens[rows, ] <- run_meboot(x[rows], series_trim(col, x[rows], g))
    }
    ens
  })
  names(ensembles) <- boot_cols
  
  ## ---- Assemble one data frame per replicate -----------------------------
  base_df <- data[id_cols]
  
  if (output == "long") {
    boot_list <- lapply(seq_len(reps), function(r) {
      df <- base_df
      for (col in boot_cols) {
        df[[col]] <- ensembles[[col]][, r]
      }
      if (as_matrix) df <- as.matrix(df)
      df
    })
  }
  
  if (output == "wide") {
    wide_id_cols <- setdiff(id_cols, names_from)
    
    # Each id / names_from combination must identify a single row
    if (anyDuplicated(data[c(wide_id_cols, names_from)]) > 0) {
      stop("id_cols and names_from do not uniquely identify rows, so the ",
           "data cannot be pivoted wider. Check id_cols.")
    }
    
    # Pivot ONCE: a template whose cells hold the original row numbers.
    # Every replicate is then filled in by indexing, instead of calling
    # pivot_wider() reps times.
    template <- base_df[c(wide_id_cols, names_from)]
    template$.row_index <- seq_len(nrow(template))
    template <- tidyr::pivot_wider(template,
                                   id_cols = tidyselect::all_of(wide_id_cols),
                                   names_from = tidyselect::all_of(names_from),
                                   values_from = ".row_index")
    template <- as.data.frame(template)
    
    unit_names <- setdiff(names(template), wide_id_cols)
    row_index <- as.matrix(template[unit_names])  # NA where a unit has no row
    id_part <- template[wide_id_cols]
    
    # Column names follow pivot_wider(): unit names for one column,
    # <column><names_sep><unit> when several columns are bootstrapped
    wide_names <- if (length(boot_cols) == 1) {
      unit_names
    } else {
      unlist(lapply(boot_cols, function(col) paste(col, unit_names, sep = names_sep)))
    }
    
    boot_list <- lapply(seq_len(reps), function(r) {
      vals <- do.call(cbind, lapply(boot_cols, function(col) {
        v <- ensembles[[col]][, r]
        matrix(v[row_index], nrow = nrow(row_index))
      }))
      colnames(vals) <- wide_names
      
      if (as_matrix) {
        return(cbind(as.matrix(id_part), vals))
      }
      df <- cbind(id_part, as.data.frame(vals, stringsAsFactors = FALSE))
      rownames(df) <- NULL
      df
    })
  }
  
  names(boot_list) <- paste0("rep_", seq_len(reps))
  
  boot_list
}


## ===========================================================================
## 2. synth_lasso_boot()
## ===========================================================================

##############################################################################
# synth_lasso_boot(): LASSO synthetic control with maximum entropy bootstrap
#                     effect and cumulative-effect intervals
#
# Combines three steps in one call:
#   1. synth_lasso_function()        - LASSO synthetic control on the observed
#                                      data (weights, predictions, ATT,
#                                      sensitivity bounds, optional jackknife
#                                      conformal intervals)
#   2. meboot_multi()                - maximum entropy bootstrap of the outcome
#                                      (and any control variables) for every
#                                      unit, split into pre / post
#   3. bootstrap_cum_effects_fixed() - refit the LASSO on each bootstrap
#                                      replicate at the fixed lambda, and
#                                      summarise the effects and cumulative
#                                      effects
#
# Uses meboot_multi(), defined earlier in this file.
#
# Data: LONG format panel, one row per unit and time period
#   (e.g. the augsynth `kansas` data).
#
# ---- Data arguments --------------------------------------------------------
#   data              long data frame
#   unit_col          column holding unit names (e.g. "state")
#   time_cols         column(s) that order time (e.g. c("year", "qtr")).
#                     Periods are numbered 1, 2, ... in sorted order; this
#                     number is the `rowid` column in the output.
#   outcome_col       outcome column (e.g. "gdpcapita")
#   treated_unit      name of the treated unit (e.g. "Kansas")
#   treatment_period  the FIRST post-treatment period, given as the value(s)
#                     of time_cols in that period:
#                       - one time column:  c(year = 2013), or just 2013
#                       - several columns:  c(year = 2012, qtr = 2), or
#                                           unnamed in time_cols order,
#                                           c(2012, 2)
#                     Values are matched exactly against the data (a list
#                     can hold mixed types, e.g. a Date). Every period
#                     before it is pre-treatment; it and every later period
#                     are post-treatment. At least 2 pre-treatment periods
#                     are required.
#
# ---- Control variables ------------------------------------------------------
#   controls          optional character vector of extra numeric columns in
#                     `data` (e.g. c("popestimate", "avg_wkly_wage")). Each
#                     is pivoted wide and added to the LASSO as predictors
#                     named <control>_<unit> (e.g. "popestimate_Texas"), and
#                     is bootstrapped by meboot along with the outcome.
#   control_units     which units' control series enter the LASSO:
#                       "all" (default), "donors", or "treated".
#                     NOTE: the treated unit's own controls in the post
#                     period may themselves be affected by the treatment.
#                     "donors" avoids that.
#
# ---- LASSO arguments --------------------------------------------------------
#   intercept       include an intercept (default TRUE)
#   lower_constr    lower limit(s) for the coefficients. Either
#                     - one number: applies to every donor column, or
#                     - a vector with at most one unnamed value (the donor
#                       default) and named values for specific columns,
#                       e.g. c(0, Texas = -0.5, popestimate_Kansas = -Inf)
#                   glmnet requires every lower limit to be <= 0.
#                   Default 0.
#   upper_constr    upper limit(s), same format as lower_constr. glmnet
#                   requires every upper limit to be >= 0. Default 1.
#   control_lower   default lower limit for control columns (default -Inf)
#   control_upper   default upper limit for control columns (default Inf)
#                   Named entries in lower_constr / upper_constr override
#                   these for individual control columns.
#   unpenalized     optional character vector of predictor columns that are
#                   NOT penalized (penalty.factor = 0), i.e. always kept in
#                   the model: donor names (e.g. "Texas") and/or control
#                   columns (e.g. "popestimate_Kansas"). Their limits still
#                   apply.
#   standardize     passed to glmnet (default FALSE). TRUE standardizes
#                   EVERY predictor, donors included, which changes the
#                   donor weights even without controls. Coefficients and
#                   limits stay on the original scale.
#   scale_controls  TRUE rescales only the control columns so their
#                   pre-period standard deviation matches the average
#                   pre-period SD of the donor columns (default FALSE).
#                   Donors are untouched; control coefficients and limits
#                   are reported / applied in the original units. The same
#                   scale factors are used for every bootstrap refit.
#                   Ignored (with a warning) when standardize = TRUE.
#   nlambda         lambda path length for cv.glmnet (default 1000)
#   lambda_choice   which lambda is used for predictions, weights and the
#                   bootstrap refits:
#                     "min_path"   smallest lambda on the path,
#                                  min(cv.fit$lambda). This is what the
#                                  existing scripts use, and it does not
#                                  depend on the CV folds.
#                     "lambda.min" lambda with the lowest CV error
#                     "lambda.1se" largest lambda within 1 SE of the minimum
#   seed            seed for the cv.glmnet folds
#
# ---- Inference on the observed fit ------------------------------------------
#   M               sensitivity-bound multiplier on the largest absolute
#                   pre-period residual (default 1)
#   jack_conform    TRUE adds jackknife conformal prediction intervals for the
#                   post period (needs conformalInference)
#                   The intervals come from conformalInference::conformal.pred.jack().
#                   That function mishandles a single prediction column (it then
#                   uses only the first leave-one-out residual, whatever the level),
#                   so the prediction function here returns two identical columns
#                   and the first is used.
#   conform_level   coverage of the conformal intervals (default 0.95)
#
# ---- Bootstrap arguments ------------------------------------------------------
#   boot            TRUE (default) runs the bootstrap
#   reps            number of bootstrap replicates (default 1000)
#   trim            meboot limits, passed to meboot_multi(). One list for
#                   every bootstrapped column, or a list named by column.
#                   Settings: trim, xmin, xmax, xmin_rel, xmax_rel, e.g.
#                     list(gdpcapita   = list(xmin = 15000),
#                          popestimate = list(xmin = 450000))
#                   With boot_by_unit = TRUE, limits can differ by unit:
#                     xmin = c(15000, Alaska = 30000)   named by unit, or
#                     xmin_rel = 0.9   90% of EACH unit's own minimum
#                     xmax_rel = 1.1   110% of each unit's own maximum
#                   (relative limits need strictly positive data)
#   boot_by_unit    FALSE (default) bootstraps each variable as one stacked
#                   series in the row order of `data` (reproduces the Kansas
#                   script). TRUE bootstraps each unit's series separately.
#   shared_draws    with boot_by_unit = TRUE: TRUE uses the same random
#                   numbers for every unit (within each variable and
#                   replicate), so units are nudged up or down together by
#                   rank and keep their co-movement. FALSE (default) draws
#                   independently for each unit. Ignored, with a warning,
#                   when boot_by_unit = FALSE.
#   boot_nlambda    lambda path length for the bootstrap glmnet refits
#                   (default 100, as in bootstrap_cum_effects_fixed)
#   conf_int        interval level for bootstrap effects (default 0.95)
#   boot_seed       seed for meboot (default = seed)
#   keep_boot       TRUE also returns the bootstrap replicate matrices
#   ...             further arguments passed to meboot::meboot()
#
# ---- Value --------------------------------------------------------------------
# A list:
#   effects    one row per period: rowid, time columns, period label,
#              post (TRUE from the first treated period on),
#              treated_observed, treated_predicted, att, att_post, cum_att,
#              bounds_lo / bounds_hi, (pred_jack_lo / pred_jack_hi),
#              eff_low / eff_high / average_eff (bootstrap effect intervals),
#              cum_eff_low / cum_eff_high / average_cum_effect (bootstrap
#              cumulative-effect intervals). Pre-period inference columns
#              are NA.
#   fit_stats  CV MSE, lambda used, pre-period MSE of the fit, M bound
#   weights    non-zero coefficients at the lambda used, with their type
#              (intercept / donor / control) and whether they were
#              unpenalized
#   model      the cv.glmnet object
#   boot       list: effects and cum_effects (post periods x reps matrices),
#              and the replicates if keep_boot = TRUE (NULL if boot = FALSE)
#   settings   the main settings used, including the limits, penalty
#              factors and scale factors for every predictor
##############################################################################

synth_lasso_boot <- function(data,
                             unit_col,
                             time_cols,
                             outcome_col,
                             treated_unit,
                             treatment_period,
                             controls = NULL,
                             control_units = c("all", "donors", "treated"),
                             intercept = TRUE,
                             lower_constr = 0,
                             upper_constr = 1,
                             control_lower = -Inf,
                             control_upper = Inf,
                             unpenalized = NULL,
                             standardize = FALSE,
                             scale_controls = FALSE,
                             nlambda = 1000,
                             lambda_choice = c("min_path", "lambda.min", "lambda.1se"),
                             seed,
                             M = 1,
                             jack_conform = FALSE,
                             conform_level = 0.95,
                             boot = TRUE,
                             reps = 1000,
                             trim = list(trim = 0.10, xmin = NULL, xmax = NULL),
                             boot_by_unit = FALSE,
                             shared_draws = FALSE,
                             boot_nlambda = 100,
                             conf_int = 0.95,
                             boot_seed = seed,
                             keep_boot = FALSE,
                             ...) {
  
  lambda_choice <- match.arg(lambda_choice)
  control_units <- match.arg(control_units)
  data <- as.data.frame(data)
  vars <- c(outcome_col, controls)
  
  ## ---- Checks ------------------------------------------------------------
  if (missing(seed)) stop("Please supply a seed for reproducibility.")
  
  missing_cols <- setdiff(c(unit_col, time_cols, vars), names(data))
  if (length(missing_cols) > 0) {
    stop("Column(s) not found in data: ", paste(missing_cols, collapse = ", "))
  }
  if (outcome_col %in% controls) stop("outcome_col cannot also be a control.")
  if (anyDuplicated(controls)) stop("controls contains duplicates.")
  
  not_numeric <- vars[!vapply(data[vars], is.numeric, logical(1))]
  if (length(not_numeric) > 0) {
    stop("Outcome and controls must be numeric: ", paste(not_numeric, collapse = ", "))
  }
  has_na <- vars[vapply(data[vars], anyNA, logical(1))]
  if (length(has_na) > 0) {
    stop("Missing values in: ", paste(has_na, collapse = ", "))
  }
  if (!treated_unit %in% data[[unit_col]]) {
    stop("treated_unit '", treated_unit, "' not found in ", unit_col, ".")
  }
  if (boot && !exists("meboot_multi", mode = "function")) {
    stop("meboot_multi() is not loaded. Run: ",
         "source(\"Synthetic_control_functions/synthetic_control_functions_claude.R\")")
  }
  if (jack_conform && !requireNamespace("conformalInference", quietly = TRUE)) {
    stop("jack_conform = TRUE needs the conformalInference package.")
  }
  
  ## ---- Number the time periods (rowid) -----------------------------------
  time_key <- do.call(paste, c(data[time_cols], sep = "\r"))
  periods <- unique(data[time_cols])
  periods <- periods[do.call(order, unname(as.list(periods))), , drop = FALSE]
  rownames(periods) <- NULL
  period_key <- do.call(paste, c(periods, sep = "\r"))
  data$.rowid <- match(time_key, period_key)
  n_periods <- nrow(periods)
  
  ## ---- Locate the first treated period -----------------------------------
  period_labels <- do.call(paste, c(periods, sep = "_"))
  tp <- as.list(treatment_period)
  
  if (length(tp) != length(time_cols)) {
    stop("treatment_period needs one value for each time column (",
         paste(time_cols, collapse = ", "), "), e.g. ",
         if (length(time_cols) == 1) paste0("c(", time_cols, " = 2013)")
         else paste0("c(", paste(time_cols, "= ...", collapse = ", "), ")"), ".")
  }
  tp_names <- names(tp)
  if (is.null(tp_names) || all(tp_names == "")) {
    names(tp) <- time_cols
  } else if (any(tp_names == "") || !setequal(tp_names, time_cols)) {
    stop("treatment_period names must match time_cols (",
         paste(time_cols, collapse = ", "), "), or leave them all unnamed.")
  } else {
    tp <- tp[time_cols]
  }
  
  is_tp <- Reduce(`&`, lapply(time_cols, function(tc) {
    as.character(periods[[tc]]) == as.character(tp[[tc]])
  }))
  if (!any(is_tp)) {
    stop("treatment_period (",
         paste(names(tp), vapply(tp, as.character, character(1)), sep = " = ", collapse = ", "),
         ") does not match any period in the data. Periods run from ",
         period_labels[1], " to ", period_labels[n_periods], ".")
  }
  treatment_rowid <- which(is_tp)
  if (treatment_rowid < 3) {
    stop("At least 2 pre-treatment periods are needed; the first treated period is ",
         period_labels[treatment_rowid], ".")
  }
  pre  <- seq_len(n_periods) < treatment_rowid
  post <- !pre
  
  ## ---- Units and predictor layout ----------------------------------------
  units <- unique(as.character(data[[unit_col]]))
  if (anyDuplicated(data[c(".rowid", unit_col)]) > 0) {
    stop("More than one row per unit and period. Check unit_col / time_cols.")
  }
  donors <- setdiff(units, treated_unit)
  
  ctrl_units <- switch(control_units,
                       all     = units,
                       donors  = donors,
                       treated = treated_unit)
  control_names <- unlist(lapply(controls, function(v) paste(v, ctrl_units, sep = "_")))
  x_names <- c(donors, control_names)
  
  if (anyDuplicated(x_names)) {
    stop("Predictor names clash (a unit name matches a <control>_<unit> name).")
  }
  
  # Build the X matrix from a list of period x unit matrices, one per variable
  build_X <- function(var_mats, scaled = TRUE) {
    X <- var_mats[[outcome_col]][, donors, drop = FALSE]
    for (v in controls) {
      Xc <- var_mats[[v]][, ctrl_units, drop = FALSE]
      colnames(Xc) <- paste(v, ctrl_units, sep = "_")
      X <- cbind(X, Xc)
    }
    if (scaled) X <- sweep(X, 2, x_scale, `*`)
    X
  }
  
  ## ---- Coefficient limits and penalty factors ----------------------------
  limit_vector <- function(spec, donor_fallback, control_default, arg) {
    nm <- names(spec)
    if (is.null(nm)) nm <- rep("", length(spec))
    unnamed <- spec[nm == ""]
    if (length(unnamed) > 1) {
      stop(arg, ": give at most one unnamed value (the default for donor columns).")
    }
    donor_default <- if (length(unnamed) == 1) unname(unnamed) else donor_fallback
    
    out <- c(stats::setNames(rep(donor_default, length(donors)), donors),
             stats::setNames(rep(control_default, length(control_names)), control_names))
    named <- spec[nm != ""]
    bad <- setdiff(names(named), x_names)
    if (length(bad) > 0) {
      stop(arg, ": unknown predictor name(s): ", paste(bad, collapse = ", "),
           ". Use donor names or <control>_<unit> names.")
    }
    out[names(named)] <- named
    out[x_names]
  }
  
  lower_vec <- limit_vector(lower_constr, 0, control_lower, "lower_constr")
  upper_vec <- limit_vector(upper_constr, 1, control_upper, "upper_constr")
  
  if (any(lower_vec > 0)) {
    stop("glmnet requires every lower limit to be <= 0. Check: ",
         paste(names(lower_vec)[lower_vec > 0], collapse = ", "))
  }
  if (any(upper_vec < 0)) {
    stop("glmnet requires every upper limit to be >= 0. Check: ",
         paste(names(upper_vec)[upper_vec < 0], collapse = ", "))
  }
  
  bad_unpen <- setdiff(unpenalized, x_names)
  if (length(bad_unpen) > 0) {
    stop("unpenalized: unknown predictor name(s): ", paste(bad_unpen, collapse = ", "),
         ". Use donor names or <control>_<unit> names.")
  }
  penalty_vec <- stats::setNames(ifelse(x_names %in% unpenalized, 0, 1), x_names)
  if (all(penalty_vec == 0)) stop("At least one predictor must be penalized.")
  
  # glmnet ignores the vectors' names; order matches x_names.
  # lower_fit / upper_fit are the limits on the (possibly rescaled) columns.
  fit_glmnet <- function(x, y, nl) {
    glmnet::glmnet(x, y, family = "gaussian", alpha = 1, standardize = standardize,
                   intercept = intercept, nlambda = nl,
                   lower.limits = unname(lower_fit), upper.limits = unname(upper_fit),
                   penalty.factor = unname(penalty_vec))
  }
  
  ## ---- Observed wide matrices: periods x units, one per variable ---------
  unit_index <- match(as.character(data[[unit_col]]), units)
  obs_mats <- lapply(vars, function(v) {
    m <- matrix(NA_real_, n_periods, length(units), dimnames = list(NULL, units))
    m[cbind(data$.rowid, unit_index)] <- data[[v]]
    m
  })
  names(obs_mats) <- vars
  if (anyNA(obs_mats[[outcome_col]])) {
    stop("The panel is unbalanced: some units are missing periods.")
  }
  
  y_obs <- obs_mats[[outcome_col]][, treated_unit]
  
  ## ---- Optional rescaling of the control columns -------------------------
  if (standardize && scale_controls) {
    warning("scale_controls is ignored when standardize = TRUE ",
            "(glmnet already standardizes every column).")
    scale_controls <- FALSE
  }
  
  x_scale <- stats::setNames(rep(1, length(x_names)), x_names)
  if (scale_controls && length(controls) > 0) {
    X_raw <- build_X(obs_mats, scaled = FALSE)
    pre_sd <- apply(X_raw[pre, , drop = FALSE], 2, stats::sd)
    if (any(pre_sd[control_names] == 0)) {
      stop("scale_controls: these control columns are constant in the pre-period: ",
           paste(control_names[pre_sd[control_names] == 0], collapse = ", "))
    }
    donor_sd <- mean(pre_sd[donors])
    x_scale[control_names] <- donor_sd / pre_sd[control_names]
  }
  
  # Limits on the scale glmnet sees: coefficient on (x * k) is b / k
  lower_fit <- lower_vec / x_scale
  upper_fit <- upper_vec / x_scale
  
  X_obs <- build_X(obs_mats)
  
  ## ---- 1. LASSO synthetic control on the observed data -------------------
  set.seed(seed)
  cv.fit <- glmnet::cv.glmnet(x = X_obs[pre, , drop = FALSE], y = y_obs[pre],
                              family = "gaussian", alpha = 1, standardize = standardize,
                              intercept = intercept, nlambda = nlambda,
                              lower.limits = unname(lower_fit),
                              upper.limits = unname(upper_fit),
                              penalty.factor = unname(penalty_vec))
  
  lambda_used <- switch(lambda_choice,
                        min_path   = min(cv.fit$lambda),
                        lambda.min = cv.fit$lambda.min,
                        lambda.1se = cv.fit$lambda.1se)
  
  # Coefficients at the lambda actually used for the predictions
  w <- as.matrix(stats::coef(cv.fit, s = lambda_used))
  # Convert coefficients on rescaled control columns back to original units
  w[x_names, 1] <- w[x_names, 1] * x_scale
  weights <- data.frame(Unit = rownames(w), Weights = w[, 1], row.names = NULL)
  weights$Type <- ifelse(weights$Unit == "(Intercept)", "intercept",
                         ifelse(weights$Unit %in% donors, "donor", "control"))
  weights$Unpenalized <- weights$Unit %in% unpenalized
  weights <- weights[weights$Weights != 0 | weights$Unpenalized, ]
  weights <- weights[order(-weights$Weights), ]
  rownames(weights) <- NULL
  
  treated_predicted <- as.numeric(stats::predict(cv.fit, newx = X_obs, s = lambda_used))
  att <- y_obs - treated_predicted
  
  effects <- data.frame(rowid = seq_len(n_periods), periods,
                        period = period_labels,
                        post = post,
                        treated_observed = y_obs,
                        treated_predicted = treated_predicted,
                        att = att,
                        att_post = ifelse(post, att, NA_real_),
                        cum_att = NA_real_)
  effects$cum_att[post] <- cumsum(att[post])
  
  # Sensitivity bounds: M x largest absolute pre-period residual
  bound <- M * max(abs(att[pre]))
  effects$bounds_lo <- ifelse(post, y_obs - (treated_predicted + bound), NA_real_)
  effects$bounds_hi <- ifelse(post, y_obs - (treated_predicted - bound), NA_real_)
  
  fit_stats <- data.frame(
    Name  = c("CV MSE (min)", "Lambda used", "Pre-period MSE", "M bound"),
    Stats = c(min(cv.fit$cvm), lambda_used, mean(att[pre]^2), bound)
  )
  
  ## ---- Optional jackknife conformal intervals ----------------------------
  if (jack_conform) {
    train_fun <- function(x, y, out = NULL) fit_glmnet(x, y, nlambda)
    # conformal.pred.jack() mishandles a single prediction column (t() of a vector makes res 1 x n, so only
    # the first leave-one-out residual is used). Returning two identical columns avoids this; column 1 is used.
    pred_fun <- function(out, newx) { p <- stats::predict(out, newx, s = lambda_used); cbind(p, p) }
    
    limits_jack <- conformalInference::conformal.pred.jack(
      x = X_obs[pre, , drop = FALSE], y = y_obs[pre],
      x0 = X_obs[post, , drop = FALSE],
      train.fun = train_fun, predict.fun = pred_fun,
      alpha = 1 - conform_level, verbose = FALSE, plus = FALSE)
    
    effects$pred_jack_lo <- NA_real_
    effects$pred_jack_hi <- NA_real_
    effects$pred_jack_lo[post] <- as.numeric(as.matrix(limits_jack$lo)[, 1])
    effects$pred_jack_hi[post] <- as.numeric(as.matrix(limits_jack$up)[, 1])
  }
  
  ## ---- 2 + 3. Bootstrap and refit ----------------------------------------
  boot_out <- NULL
  
  if (boot) {
    boot_df <- data[c(unit_col, ".rowid", vars)]
    boot_df[[unit_col]] <- as.character(boot_df[[unit_col]])
    if (boot_by_unit) {
      boot_df <- boot_df[order(boot_df[[unit_col]], boot_df$.rowid), ]
    }
    
    replicates <- meboot_multi(boot_df,
                               boot_cols  = vars,
                               id_cols    = c(unit_col, ".rowid"),
                               reps       = reps,
                               trim       = trim,
                               group_col  = if (boot_by_unit) unit_col else NULL,
                               shared_draws = shared_draws,
                               output     = "wide",
                               names_from = unit_col,
                               seed       = boot_seed,
                               as_matrix  = TRUE,
                               ...)
    
    # Split each replicate back into period-ordered matrices, one per variable
    replicates <- lapply(replicates, function(m) {
      m <- m[order(m[, ".rowid"]), , drop = FALSE]
      mats <- lapply(vars, function(v) {
        cols <- if (length(vars) == 1) units else paste(v, units, sep = "_")
        out <- m[, cols, drop = FALSE]
        colnames(out) <- units
        out
      })
      names(mats) <- vars
      mats
    })
    
    # Refit the LASSO on each bootstrapped pre period, predict the
    # bootstrapped post period at the fixed lambda, and compare with the
    # OBSERVED treated outcome (as in bootstrap_cum_effects_fixed)
    eff_sim <- vapply(replicates, function(mats) {
      X_b <- build_X(mats)
      y_b <- mats[[outcome_col]][, treated_unit]
      fit <- fit_glmnet(X_b[pre, , drop = FALSE], y_b[pre], boot_nlambda)
      pred <- stats::predict(fit, newx = X_b[post, , drop = FALSE], s = lambda_used)
      y_obs[post] - as.numeric(pred)
    }, numeric(sum(post)))
    eff_sim <- matrix(eff_sim, nrow = sum(post))  # post periods x reps
    
    cum_eff_sim <- apply(eff_sim, 2, cumsum)
    cum_eff_sim <- matrix(cum_eff_sim, nrow = sum(post))
    
    a <- (1 - conf_int) / 2
    row_q <- function(mat, p) apply(mat, 1, stats::quantile, probs = p, names = FALSE)
    
    for (nm in c("eff_low", "eff_high", "average_eff",
                 "cum_eff_low", "cum_eff_high", "average_cum_effect")) {
      effects[[nm]] <- NA_real_
    }
    
    effects$eff_low[post]            <- row_q(eff_sim, a)
    effects$eff_high[post]           <- row_q(eff_sim, 1 - a)
    effects$average_eff[post]        <- rowMeans(eff_sim)
    effects$cum_eff_low[post]       <- row_q(cum_eff_sim, a)
    effects$cum_eff_high[post]      <- row_q(cum_eff_sim, 1 - a)
    effects$average_cum_effect[post] <- rowMeans(cum_eff_sim)
    
    boot_out <- list(effects = eff_sim, cum_effects = cum_eff_sim)
    if (keep_boot) boot_out$replicates <- replicates
  }
  
  effects$treated_unit <- treated_unit
  
  settings <- list(treated_unit = treated_unit, treatment_period = tp,
                   treatment_rowid = treatment_rowid,
                   treatment_label = period_labels[treatment_rowid],
                   controls = controls, control_units = control_units,
                   intercept = intercept,
                   standardize = standardize, scale_controls = scale_controls,
                   predictors = data.frame(predictor = x_names,
                                           lower = unname(lower_vec),
                                           upper = unname(upper_vec),
                                           penalty_factor = unname(penalty_vec),
                                           scale_factor = unname(x_scale)),
                   lambda_choice = lambda_choice, lambda_used = lambda_used,
                   n_pre = sum(pre), n_post = sum(post),
                   boot = boot, reps = if (boot) reps else NA,
                   boot_by_unit = boot_by_unit, shared_draws = shared_draws,
                   trim = trim, conf_int = conf_int,
                   seed = seed, boot_seed = boot_seed)
  
  list(effects = effects, fit_stats = fit_stats, weights = weights,
       model = cv.fit, boot = boot_out, settings = settings)
}


## ===========================================================================
## 3. synth_lasso_placebo()
## ===========================================================================

##############################################################################
# synth_lasso_placebo(): LASSO synthetic control with exact placebo
#                        (permutation) tests and confidence intervals
#
# Same model and arguments as synth_lasso_boot(), but the bootstrap intervals
# are replaced by exact placebo tests:
#   1. Fit the treated unit (identical to synth_lasso_boot(), boot = FALSE).
#   2. Drop the treated unit. Fit every donor in turn as a placebo "treated"
#      unit, from the remaining donors, with the same settings.
#   3. Optionally keep only placebos with a good pre-treatment fit (Abadie,
#      Diamond and Hainmueller 2010): pre-period MSPE <= abadie_cutoff x the
#      treated unit's pre-period MSPE (typically 5 or 2).
#   4. For every post-treatment period, test the treated unit's gap against the
#      placebo gaps, for the per-period effect and for the cumulative effect
#      through that period. Gaps are relative: (observed - synthetic) /
#      synthetic; cumulative gaps are sum(gap) / sum(synthetic).
#
# Exact tests with N placebos (N + 1 units), alpha = 1 - conf_level:
#   absolute  p = (1 + #{|placebo gap| >= |treated gap|}) / (N + 1)
#             CI = treated gap +/- the k-th largest |placebo gap|,
#             k = floor(alpha x (N + 1))
#   signed    p = min(1, 2 x min(lower-tail p, upper-tail p))
#             CI = [treated gap - k2-th largest placebo gap,
#                   treated gap - k2-th smallest placebo gap],
#             k2 = floor(alpha / 2 x (N + 1))
#   The CIs are the effect sizes the test does not reject (test inversion);
#   a CI excludes zero exactly when p <= alpha. If k (or k2) is 0 there are
#   too few placebos for conf_level and the CI is infinite (with a warning):
#   the absolute test needs N >= 1/alpha - 1 placebos, the signed test
#   N >= 2/alpha - 1 (19 and 39 at 95%).
#   CIs are computed on the relative scale and converted to outcome units by
#   multiplying by the treated unit's synthetic value (or its running total).
#
# Arguments: as synth_lasso_boot() for the data, controls, LASSO fit, M bounds
# and jackknife (see its documentation). Bootstrap arguments are replaced by:
#   conf_level     confidence level of the placebo CIs (default 0.95)
#   abadie_cutoff  NULL (default) keeps every placebo; a number (e.g. 5 or 2)
#                  keeps placebos whose pre-period MSPE is at most that many
#                  times the treated unit's
#   n_cores        number of cores for the placebo fits (default 1; > 1 runs
#                  them in parallel)
#
# Named settings (lower_constr / upper_constr / unpenalized entries naming a
# donor or <control>_<unit> column) are applied in a placebo run only when that
# column exists in it; e.g. unpenalized = "Colorado" is dropped in the run
# where Colorado is the placebo.
#
# Value: a list
#   effects    one row per period, as in synth_lasso_boot() (without bootstrap
#              columns), plus for post-treatment periods:
#              att_rel, cum_att_rel         relative gap and cumulative gap
#              p_abs, eff_low_abs, eff_high_abs              per-period, absolute
#              p_signed, eff_low_signed, eff_high_signed     per-period, signed
#              cum_p_abs, cum_eff_low_abs, cum_eff_high_abs  cumulative, absolute
#              cum_p_signed, cum_eff_low_signed, cum_eff_high_signed
#                                                            cumulative, signed
#              (intervals in outcome units)
#   fit_stats, weights, model   as in synth_lasso_boot()
#   placebo    list: gaps and cum_gaps (post periods x placebos, relative),
#              fit (each placebo's pre-period MSPE, its ratio to the treated
#              unit's, and whether it was kept), n_kept, n_failed
#   settings   the settings used
##############################################################################

synth_lasso_placebo <- function(data,
                                unit_col,
                                time_cols,
                                outcome_col,
                                treated_unit,
                                treatment_period,
                                controls = NULL,
                                control_units = c("all", "donors", "treated"),
                                intercept = TRUE,
                                lower_constr = 0,
                                upper_constr = 1,
                                control_lower = -Inf,
                                control_upper = Inf,
                                unpenalized = NULL,
                                standardize = FALSE,
                                scale_controls = FALSE,
                                nlambda = 1000,
                                lambda_choice = c("min_path", "lambda.min", "lambda.1se"),
                                seed,
                                M = 1,
                                jack_conform = FALSE,
                                conform_level = 0.95,
                                conf_level = 0.95,
                                abadie_cutoff = NULL,
                                n_cores = 1) {
  
  control_units <- match.arg(control_units)
  lambda_choice <- match.arg(lambda_choice)
  data <- as.data.frame(data)
  
  ## ---- Checks ------------------------------------------------------------
  if (missing(seed)) stop("Please supply a seed for reproducibility.")
  if (!is.numeric(conf_level) || length(conf_level) != 1 || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be a single number between 0 and 1, e.g. 0.95.")
  }
  if (!is.null(abadie_cutoff) &&
      (!is.numeric(abadie_cutoff) || length(abadie_cutoff) != 1 || abadie_cutoff <= 0)) {
    stop("abadie_cutoff must be NULL or a single positive number, e.g. 5 or 2.")
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) {
    stop("n_cores must be a single number >= 1.")
  }
  alpha <- round(1 - conf_level, 10)   # avoid floating-point error in the rank (e.g. 1 - 0.9)
  
  # Arguments shared by every fit (treated and placebo)
  common <- list(unit_col = unit_col, time_cols = time_cols, outcome_col = outcome_col,
                 treatment_period = treatment_period, controls = controls,
                 control_units = control_units, intercept = intercept,
                 control_lower = control_lower, control_upper = control_upper,
                 standardize = standardize, scale_controls = scale_controls,
                 nlambda = nlambda, lambda_choice = lambda_choice, seed = seed, boot = FALSE)
  
  ## ---- 1. Treated unit (same as synth_lasso_boot, no bootstrap) ----------
  main <- do.call(synth_lasso_boot, c(common, list(
    data = data, treated_unit = treated_unit, lower_constr = lower_constr,
    upper_constr = upper_constr, unpenalized = unpenalized, M = M,
    jack_conform = jack_conform, conform_level = conform_level)))
  eff <- main$effects
  post <- eff$post
  treated_pre_mspe <- mean(eff$att[!post]^2)
  
  ## ---- 2. Placebo fits: every donor, treated unit removed ----------------
  pl_data <- data[as.character(data[[unit_col]]) != treated_unit, ]
  pl_units <- unique(as.character(pl_data[[unit_col]]))
  if (length(pl_units) < 3) stop("At least 3 donor units are needed for placebo runs.")
  
  # Keep only named settings that exist as predictors in a given placebo run
  predictor_names <- function(p) {
    donors_p <- setdiff(pl_units, p)
    ctrl_p <- switch(control_units, all = pl_units, donors = donors_p, treated = p)
    c(donors_p, unlist(lapply(controls, function(v) paste(v, ctrl_p, sep = "_"))))
  }
  keep_named <- function(spec, ok) {
    nm <- names(spec)
    if (is.null(nm)) return(spec)
    spec[nm == "" | nm %in% ok]
  }
  placebo_fit <- function(p) {
    ok <- predictor_names(p)
    unpen_p <- intersect(unpenalized, ok)
    m <- tryCatch(
      do.call(synth_lasso_boot, c(common, list(
        data = pl_data, treated_unit = p,
        lower_constr = keep_named(lower_constr, ok), upper_constr = keep_named(upper_constr, ok),
        unpenalized = if (length(unpen_p) > 0) unpen_p else NULL))),
      error = function(e) conditionMessage(e))
    if (is.character(m)) return(list(error = m))
    e <- m$effects
    list(gap_rel = e$att[e$post] / e$treated_predicted[e$post],
         cum_gap_rel = cumsum(e$att[e$post]) / cumsum(e$treated_predicted[e$post]),
         pre_mspe = mean(e$att[!e$post]^2))
  }
  
  if (n_cores > 1) {
    cl <- parallel::makeCluster(min(n_cores, length(pl_units)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    fn_env <- environment(sys.function())
    parallel::clusterExport(cl, intersect(c("synth_lasso_boot", "meboot_multi"), ls(fn_env)), envir = fn_env)
    fits <- parallel::parLapplyLB(cl, pl_units, placebo_fit)
  } else {
    fits <- lapply(pl_units, placebo_fit)
  }
  names(fits) <- pl_units
  
  failed <- vapply(fits, function(f) !is.null(f$error), logical(1))
  if (any(failed)) {
    warning(sum(failed), " placebo fit(s) failed and were left out: ",
            paste(names(fits)[failed], collapse = ", "), ". First error: ", fits[[which(failed)[1]]]$error)
  }
  fits <- fits[!failed]
  if (length(fits) == 0) stop("All placebo fits failed.")
  
  G  <- sapply(fits, `[[`, "gap_rel");     G  <- matrix(G,  nrow = sum(post), dimnames = list(NULL, names(fits)))
  Gc <- sapply(fits, `[[`, "cum_gap_rel"); Gc <- matrix(Gc, nrow = sum(post), dimnames = list(NULL, names(fits)))
  pl_mspe <- vapply(fits, `[[`, numeric(1), "pre_mspe")
  
  ## ---- 3. Abadie filter on pre-treatment fit -------------------------------
  kept <- if (is.null(abadie_cutoff)) rep(TRUE, length(fits)) else pl_mspe <= abadie_cutoff * treated_pre_mspe
  n_kept <- sum(kept)
  fit_table <- data.frame(unit = names(fits), pre_mspe = pl_mspe,
                          mspe_ratio_to_treated = pl_mspe / treated_pre_mspe, kept = kept, row.names = NULL)
  
  if (floor(alpha * (n_kept + 1)) == 0) {
    warning(n_kept, " placebos kept: too few for a ", 100 * conf_level,
            "% absolute CI (needs ", ceiling(1 / alpha - 1), "); absolute CIs are infinite.")
  }
  if (floor(alpha / 2 * (n_kept + 1)) == 0) {
    warning(n_kept, " placebos kept: too few for a ", 100 * conf_level,
            "% signed CI (needs ", ceiling(2 / alpha - 1), "); signed CIs are infinite.")
  }
  
  ## ---- 4. Exact tests and inverted CIs, every post period ----------------
  exact_tests <- function(t, g) {
    n1 <- length(g) + 1
    p_abs <- (1 + sum(abs(g) >= abs(t))) / n1
    k <- floor(alpha * n1)
    a_k <- if (k == 0) Inf else sort(abs(g), decreasing = TRUE)[k]
    p_sgn <- min(1, 2 * min((1 + sum(g <= t)) / n1, (1 + sum(g >= t)) / n1))
    k2 <- floor(alpha / 2 * n1)
    gs <- sort(g)
    c(p_abs = p_abs, abs_lo = t - a_k, abs_hi = t + a_k, p_signed = p_sgn,
      signed_lo = if (k2 == 0) -Inf else t - gs[length(g) + 1 - k2],
      signed_hi = if (k2 == 0)  Inf else t - gs[k2])
  }
  
  pred_post <- eff$treated_predicted[post]
  att_rel <- eff$att[post] / pred_post
  cum_att_rel <- cumsum(eff$att[post]) / cumsum(pred_post)
  test_names <- c("p_abs", "abs_lo", "abs_hi", "p_signed", "signed_lo", "signed_hi")
  per <- t(vapply(seq_along(att_rel), function(h) exact_tests(att_rel[h], G[h, kept]),
                  setNames(numeric(6), test_names)))
  cum <- t(vapply(seq_along(cum_att_rel), function(h) exact_tests(cum_att_rel[h], Gc[h, kept]),
                  setNames(numeric(6), test_names)))
  
  add_col <- function(name, values) {
    eff[[name]] <<- NA_real_
    eff[[name]][post] <<- values
  }
  add_col("att_rel", att_rel)
  add_col("cum_att_rel", cum_att_rel)
  add_col("p_abs", per[, "p_abs"])
  add_col("eff_low_abs", per[, "abs_lo"] * pred_post)
  add_col("eff_high_abs", per[, "abs_hi"] * pred_post)
  add_col("p_signed", per[, "p_signed"])
  add_col("eff_low_signed", per[, "signed_lo"] * pred_post)
  add_col("eff_high_signed", per[, "signed_hi"] * pred_post)
  add_col("cum_p_abs", cum[, "p_abs"])
  add_col("cum_eff_low_abs", cum[, "abs_lo"] * cumsum(pred_post))
  add_col("cum_eff_high_abs", cum[, "abs_hi"] * cumsum(pred_post))
  add_col("cum_p_signed", cum[, "p_signed"])
  add_col("cum_eff_low_signed", cum[, "signed_lo"] * cumsum(pred_post))
  add_col("cum_eff_high_signed", cum[, "signed_hi"] * cumsum(pred_post))
  eff$treated_unit <- NULL
  eff$treated_unit <- treated_unit
  
  settings <- main$settings
  settings$boot <- NULL; settings$reps <- NULL; settings$boot_by_unit <- NULL
  settings$shared_draws <- NULL; settings$trim <- NULL; settings$conf_int <- NULL
  settings$boot_seed <- NULL
  settings <- c(settings, list(inference = "placebo", conf_level = conf_level,
                               abadie_cutoff = abadie_cutoff, treated_pre_mspe = treated_pre_mspe,
                               n_placebos = length(fits), n_placebos_kept = n_kept))
  
  list(effects = eff, fit_stats = main$fit_stats, weights = main$weights, model = main$model,
       placebo = list(gaps = G, cum_gaps = Gc, fit = fit_table, n_kept = n_kept,
                      n_failed = sum(failed)),
       settings = settings)
}


## ===========================================================================
## 4. dsc_symcal_boot()
## ===========================================================================

##############################################################################
# dsc_symcal_boot(): demeaned synthetic control with a symmetric, donor-
# calibrated fixed-weight bootstrap
#
# Model
#   Demeaned synthetic control: an intercept plus non-negative donor weights
#   that sum to one, fitted on the pre-treatment period (constrained least
#   squares, CVXR). The intercept absorbs level differences; the sum-to-one
#   constraint keeps the synthetic control following the donors' full movement.
#   Optional control series (e.g. population) enter with FREE, unpenalized
#   coefficients.
#
# Inference
#   1. Fixed-weight bootstrap: by-unit meboot replicates of every series (the
#      outcome and any controls); in each replicate, synthetic = intercept +
#      replicate donor values x the FIXED weights + replicate controls x the FIXED
#      control coefficients; replicate effect = observed - replicate synthetic
#      (per period, and running total for the cumulative effect).
#   2. Symmetric, mean-centred half-width in every post period t:
#      h_t = (conf_level) quantile of |replicate effect_t - mean replicate effect_t|.
#   3. Standardized gap k_t = |effect_t| / h_t.
#   4. Calibration: every donor is fitted as a placebo (treated unit removed)
#      with steps 1-3, using the same model (including controls). With N donors,
#      m = floor((1 - conf_level) x (N + 1)):
#        pointwise factor    k*_t = m-th largest donor k_t in period t
#        simultaneous band   smoothed-median profile q~_t (loess on the log of the
#                            median donor k_t over post periods); lambda = m-th
#                            largest donor max_t (k_t / q~_t); factor = lambda x q~_t
#   5. Intervals: effect_t +/- factor_t x h_t; prediction intervals for the
#      untreated outcome: synthetic_t +/- factor_t x h_t (cumulative: running
#      totals).
#
# Which interval answers which question
#   pointwise     the effect in ONE pre-specified period (test_period; default
#                 the last post period). Valid for that period only: scanning
#                 many periods for a striking one inflates false positives.
#   simultaneous  whether there was an effect at ANY point in the post period,
#                 which periods, and the shape of the effect path (covers every
#                 post period at once).
#
# In simulations on the Kansas and California panels (Placebo simulations/),
# the simultaneous band kept whole-path coverage at 95-97%, and pointwise
# intervals kept 94-100% coverage at a pre-specified period, for pre-periods of
# 5-80 periods (no controls). For a final-period effect, the absolute placebo
# test with the same model was more powerful.
#
# With fewer than 19 donors, m = 0 at 95% and the intervals are infinite (the
# exact calibration cannot reach that level); use conf_level = 0.90 or more
# donors.
#
# Arguments
#   data, unit_col, time_cols, outcome_col, treated_unit
#                     long panel, as in synth_lasso_boot(); balanced, no NAs
#   treatment_period  named vector of time column values for the first treated
#                     period, e.g. c(year = 2012, qtr = 2) or c(year = 1989)
#   test_period       named vector (same form) of the post period at which the
#                     pointwise tests are made, for the per-period effect and the
#                     cumulative effect through that period. NULL (default) = the
#                     last post period.
#   controls          optional character vector of numeric columns in `data`
#                     (e.g. "popestimate"), added as predictors named
#                     <control>_<unit> with free, unpenalized coefficients and
#                     bootstrapped alongside the outcome (same trim proportion,
#                     no limits). Not validated by simulation;
#                     with the LASSO model, population controls added bias.
#   control_units     whose control series enter: "treated" (default: the unit
#                     being fitted, one column per control), "donors", or "all".
#                     "donors" / "all" add one free coefficient per unit; the
#                     function stops if the free coefficients would outnumber the
#                     pre-treatment periods. NOTE: with "treated" or "all", the
#                     treated unit's own post-period controls may themselves be
#                     affected by the treatment.
#   adjust_covariates optional character vector of time-varying covariates for a
#                     POOLED adjustment ("Option D"): their effect gamma is
#                     estimated by a two-way fixed-effects regression (unit and
#                     period effects) on untreated observations only (all donor
#                     periods + the treated unit's pre-period); every unit's
#                     outcome is adjusted (outcome - gamma x covariate, using each
#                     unit's ACTUAL covariate, post-treatment included) and the
#                     method runs on the adjusted outcomes. Each placebo
#                     re-estimates gamma with its own post period left out (and
#                     the treated unit excluded). Suited to external shocks that
#                     may start at the same time as the treatment (e.g. the spread
#                     of fentanyl); assumes the treatment does not affect the
#                     covariate and a common, linear covariate effect. Cannot be
#                     combined with `controls`.
#   conf_level        confidence level (default 0.95)
#   reps              bootstrap replicates (default 500)
#   trim              meboot trim / limits (default list(trim = 0.10, xmin_rel = 0.9))
#   span              loess span for the simultaneous band's profile (default 0.75)
#   seed              random seed for the bootstrap (default 12112025)
#
# Value: a list with
#   effects   one row per period: time columns, post, observed, synthetic
#             (the counterfactual; with adjust_covariates it includes
#             gamma x the treated unit's covariate, shown in covariate_adjustment),
#             effect; for post periods also halfwidth, factor_pointwise,
#             factor_simultaneous, eff_low_pw / eff_high_pw, eff_low_sim /
#             eff_high_sim, pred_low_pw / pred_high_pw, pred_low_sim /
#             pred_high_sim, p_pointwise, and the same for the cumulative
#             effect (cum_ prefix)
#   test      the pointwise tests at test_period and the simultaneous test
#   adjustment (with adjust_covariates: gamma, clustered SE, n_obs, the
#             unadjusted estimates, and the estimates at gamma's 95% CI limits)
#   weights, intercept, control_coefs, pre_rmspe, donor_factors,
#   placebo   each donor-as-placebo's gaps and synthetic values (all periods) and
#             bootstrap half-widths (post periods, per period and cumulative)
#   failed_donors, settings
##############################################################################
dsc_symcal_boot <- function(data, unit_col, time_cols, outcome_col, treated_unit, treatment_period,
                            test_period = NULL, controls = NULL, control_units = c("treated", "donors", "all"),
                            adjust_covariates = NULL,
                            conf_level = 0.95, reps = 500,
                            trim = list(trim = 0.10, xmin_rel = 0.9), span = 0.75, seed = 12112025) {

  for (p in c("CVXR", "meboot")) if (!requireNamespace(p, quietly = TRUE)) stop("dsc_symcal_boot() needs the ", p, " package.")
  control_units <- match.arg(control_units)
  alpha <- round(1 - conf_level, 10)

  # ---- Data: ordered periods, post indicator, wide matrices --------------------
  data <- as.data.frame(data)
  if (length(controls) && length(adjust_covariates)) stop("Use either `controls` or `adjust_covariates`, not both.")
  stopifnot(all(c(unit_col, time_cols, outcome_col, controls, adjust_covariates) %in% names(data)),
            all(names(treatment_period) %in% time_cols))
  periods <- unique(data[, time_cols, drop = FALSE])
  periods <- periods[do.call(order, unname(as.list(periods))), , drop = FALSE]
  rownames(periods) <- NULL
  find_period <- function(pv, what) {
    hit <- which(Reduce(`&`, lapply(names(pv), function(nm) periods[[nm]] == pv[[nm]])))
    if (length(hit) != 1) stop(what, " (", paste(names(pv), pv, sep = " = ", collapse = ", "), ") must match exactly one period.")
    hit
  }
  t_treat <- find_period(treatment_period, "treatment_period")
  post <- seq_len(nrow(periods)) >= t_treat
  n_post <- sum(post)
  t_test <- if (is.null(test_period)) n_post else find_period(test_period, "test_period") - t_treat + 1
  if (t_test < 1) stop("test_period must be a post-treatment period.")

  d <- data.frame(unit = as.character(data[[unit_col]]),
                  .t = match(do.call(paste, unname(as.list(data[, time_cols, drop = FALSE]))),
                             do.call(paste, unname(as.list(periods)))),
                  y = data[[outcome_col]])
  for (v in c(controls, adjust_covariates)) d[[v]] <- data[[v]]
  if (anyNA(d[, c("y", controls, adjust_covariates)]) || any(table(d$unit) != nrow(periods)))
    stop("The panel must be balanced with no missing outcomes or controls.")
  if (!treated_unit %in% d$unit) stop("treated_unit not found in ", unit_col, ".")
  units <- sort(unique(d$unit)); donors <- setdiff(units, treated_unit)
  wide <- function(v) sapply(units, function(u) { du <- d[d$unit == u, ]; du[[v]][order(du$.t)] })
  Y <- wide("y")
  Zobs <- setNames(lapply(controls, wide), controls)            # one period x unit matrix per control
  Aobs <- setNames(lapply(adjust_covariates, wide), adjust_covariates)   # covariates for the pooled adjustment

  # Option D: gamma from a two-way FE regression on untreated observations, excluding the units in
  # `drop_units` entirely and the post period of `drop_post_unit`; standard errors clustered by unit
  twfe_gamma <- function(drop_units = NULL, drop_post_unit) {
    keep_u <- setdiff(units, drop_units)
    long <- expand.grid(t = seq_len(nrow(periods)), unit = keep_u, stringsAsFactors = FALSE)
    long$y <- Y[cbind(long$t, match(long$unit, units))]
    for (v in adjust_covariates) long[[v]] <- Aobs[[v]][cbind(long$t, match(long$unit, units))]
    long <- long[!(long$unit == drop_post_unit & post[long$t]), ]
    X <- model.matrix(as.formula(paste("~ factor(unit) + factor(t) +", paste(adjust_covariates, collapse = " + "))), long)
    fit <- lm.fit(X, long$y)
    b <- fit$coefficients[adjust_covariates]
    if (anyNA(b)) stop("adjust_covariates: gamma not identified (no variation beyond unit and period effects).")
    XtX_inv <- chol2inv(qr.R(qr(X)))[order(qr(X)$pivot), order(qr(X)$pivot)]
    e <- fit$residuals
    meat <- Reduce(`+`, lapply(split(seq_len(nrow(X)), long$unit), function(ix) { s <- crossprod(X[ix, , drop = FALSE], e[ix]); s %*% t(s) }))
    G <- length(unique(long$unit)); V <- G / (G - 1) * XtX_inv %*% meat %*% XtX_inv
    idx <- match(adjust_covariates, colnames(X))
    list(gamma = b, se = sqrt(diag(V)[idx]), n_obs = nrow(X))
  }
  adjust_mat <- function(Ymat, Alist, gam) {
    if (is.null(gam)) return(Ymat)
    for (v in names(gam)) Ymat <- Ymat - gam[[v]] * Alist[[v]][, colnames(Ymat), drop = FALSE]
    Ymat
  }

  # Control columns for fitting unit u from donor set D (observed or replicate matrices)
  control_cols <- function(Zlist, u, D) {
    if (!length(controls)) return(NULL)
    who <- switch(control_units, treated = u, donors = D, all = c(u, D))
    do.call(cbind, lapply(controls, function(v) { M <- Zlist[[v]][, who, drop = FALSE]; colnames(M) <- paste0(v, "_", who); M }))
  }
  n_ctrl <- length(controls) * switch(control_units, treated = 1, donors = length(donors) - 1, all = length(donors))
  if (n_ctrl + 1 >= sum(!post))
    stop("Too many free coefficients (", n_ctrl, " control coefficients + intercept) for ", sum(!post),
         " pre-treatment periods. Use control_units = \"treated\" or fewer controls.")

  # ---- Demeaned SC fit (rescaled before solving: large values can make the -----
  # ---- solver report a false "infeasible"; weights are scale-free) --------------
  fit_dsc <- function(X, yv, C = NULL) {
    sc <- mean(abs(yv[!post]))
    Xp <- X[!post, , drop = FALSE] / sc; yp <- yv[!post] / sc
    w <- CVXR::Variable(ncol(Xp)); b0 <- CVXR::Variable(1)
    if (is.null(C)) {
      resid <- yp - b0 - Xp %*% w
    } else {
      cs <- colMeans(abs(C[!post, , drop = FALSE])); cs[cs == 0] <- 1
      Cp <- sweep(C[!post, , drop = FALSE], 2, cs, "/")
      g <- CVXR::Variable(ncol(Cp))
      resid <- yp - b0 - Xp %*% w - Cp %*% g
    }
    prob <- CVXR::Problem(CVXR::Minimize(CVXR::sum_squares(resid)), list(w >= 0, sum(w) == 1))
    CVXR::psolve(prob)
    st <- CVXR::status(prob)
    if (!st %in% c("optimal", "optimal_inaccurate")) {        # retry once: explicit solver, more iterations
      try(CVXR::psolve(prob, solver = "CLARABEL", max_iter = 10000), silent = TRUE)
      st <- CVXR::status(prob)
    }
    if (!st %in% c("optimal", "optimal_inaccurate")) stop("CVXR status: ", st)
    wv <- setNames(as.numeric(CVXR::value(w)), colnames(X)); wv[abs(wv) < 1e-8] <- 0
    b <- as.numeric(CVXR::value(b0)) * sc
    gv <- if (is.null(C)) NULL else setNames(as.numeric(CVXR::value(g)) * sc / cs, colnames(C))   # original units
    synth <- b + as.numeric(X %*% wv) + if (is.null(C)) 0 else as.numeric(C %*% gv)
    list(b0 = b, w = wv, g = gv, synth = synth)
  }

  # ---- Bootstrap replicates: list of list(y = period x unit, ctrl = list) -------
  # The outcome uses `trim` (including any limits); covariates are bootstrapped separately with the
  # same trim proportion but NO limits (relative limits need strictly positive series).
  replicates_for <- function(dd) {
    dd <- dd[order(dd$unit, dd$.t), ]
    us <- sort(unique(dd$unit))
    boot_wide <- function(vars, tr, sd) {
      r <- withCallingHandlers(
        meboot_multi(dd[, c("unit", ".t", vars)], boot_cols = vars, id_cols = c("unit", ".t"), reps = reps, trim = tr,
                     group_col = "unit", output = "wide", names_from = "unit", seed = sd, as_matrix = TRUE),
        warning = function(w) if (grepl("expand.sd", conditionMessage(w))) invokeRestart("muffleWarning"))
      lapply(r, function(m) {
        m <- m[order(m[, ".t"]), , drop = FALSE]
        setNames(lapply(vars, function(v) { cn <- if (length(vars) == 1) us else paste0(v, "_", us)
          M <- m[, cn, drop = FALSE]; colnames(M) <- us; M }), vars)
      })
    }
    ry <- boot_wide("y", trim, seed)
    cov_vars <- c(controls, adjust_covariates)
    rc <- if (length(cov_vars)) boot_wide(cov_vars, list(trim = if (is.null(trim$trim)) 0.10 else trim$trim), seed + 1) else NULL
    if (length(cov_vars)) {                       # constant series (e.g. never-exposed units) cannot be resampled:
      obs <- setNames(lapply(cov_vars, function(v) sapply(us, function(u) dd[[v]][dd$unit == u][order(dd$.t[dd$unit == u])])), cov_vars)
      rc <- lapply(rc, function(rr) setNames(lapply(cov_vars, function(v) { M <- rr[[v]]
        bad <- colSums(!is.finite(M)) > 0; M[, bad] <- obs[[v]][, bad]; M }), cov_vars))   # keep them as observed
    }
    lapply(seq_along(ry), function(i) list(y = ry[[i]]$y,
      ctrl = if (length(controls)) rc[[i]][controls] else list(),
      adj  = if (length(adjust_covariates)) rc[[i]][adjust_covariates] else list()))
  }

  # ---- Fixed-weight effects, symmetric half-widths, standardized gaps -----------
  unit_stats <- function(u, D, reps_list, gam = NULL) {
    Ya <- adjust_mat(Y[, c(u, D), drop = FALSE], Aobs, gam)
    yv <- Ya[, u]
    fit <- fit_dsc(Ya[, D, drop = FALSE], yv, control_cols(Zobs, u, D))
    gap <- yv - fit$synth
    E <- sapply(reps_list, function(m) {
      ya_r <- adjust_mat(m$y[, D, drop = FALSE], m$adj, gam)
      s_r <- fit$b0 + as.numeric(ya_r[post, , drop = FALSE] %*% fit$w[D])
      if (!is.null(fit$g)) s_r <- s_r + as.numeric(control_cols(m$ctrl, u, D)[post, , drop = FALSE] %*% fit$g)
      yv[post] - s_r
    })
    E <- matrix(E, nrow = n_post)
    Ec <- matrix(apply(E, 2, cumsum), nrow = n_post)
    hw <- function(M) apply(M, 1, function(v) unname(quantile(abs(v - mean(v)), 1 - alpha)))
    g <- gap[post]; gc <- cumsum(g); h <- hw(E); hc <- hw(Ec)
    list(gap = gap, synth = fit$synth, h = h, hc = hc, k = abs(g) / h, kc = abs(gc) / hc, w = fit$w, b0 = fit$b0, g = fit$g)
  }

  # ---- Treated unit, and every donor as a placebo (treated unit removed) --------
  adj_main <- if (length(adjust_covariates)) twfe_gamma(drop_post_unit = treated_unit) else NULL
  tr <- unit_stats(treated_unit, donors, replicates_for(d), adj_main$gamma)
  reps_ex <- replicates_for(d[d$unit != treated_unit, ])
  pl <- lapply(donors, function(p) tryCatch({
    gp <- if (length(adjust_covariates)) twfe_gamma(drop_units = treated_unit, drop_post_unit = p)$gamma else NULL
    unit_stats(p, setdiff(donors, p), reps_ex, gp)
  }, error = function(e) conditionMessage(e)))
  failed <- vapply(pl, is.character, logical(1))
  failed_donors <- if (any(failed)) data.frame(donor = donors[failed], error = unlist(pl[failed])) else
    data.frame(donor = character(0), error = character(0))
  if (any(failed)) warning(sum(failed), " placebo fit(s) failed and were left out of the calibration: ",
                           paste(donors[failed], collapse = ", "), ". See $failed_donors.")
  pl <- pl[!failed]; N <- length(pl); m <- floor(alpha * (N + 1))
  kth <- function(x) if (m == 0) Inf else sort(x, decreasing = TRUE)[m]
  if (m == 0) warning("With ", N, " donors the ", 100 * conf_level, "% level cannot be reached exactly: intervals are infinite.")

  # ---- Calibration factors ------------------------------------------------------
  Kp  <- matrix(sapply(pl, `[[`, "k"),  nrow = n_post)                  # post periods x donors
  Kpc <- matrix(sapply(pl, `[[`, "kc"), nrow = n_post)
  smooth_pos <- function(v) { tt <- seq_along(v); exp(as.numeric(predict(stats::loess(log(v) ~ tt, span = span, degree = 1)))) }
  sim_factor <- function(K, kT) {
    q <- smooth_pos(apply(K, 1, median)); R <- apply(K / q, 2, max)
    list(factor = kth(R) * q, p = (1 + sum(R >= max(kT / q))) / (N + 1))
  }
  k_pw <- apply(Kp, 1, kth); k_pwc <- apply(Kpc, 1, kth)
  sim <- sim_factor(Kp, tr$k); simc <- sim_factor(Kpc, tr$kc)

  # ---- Effects table --------------------------------------------------------------
  g <- tr$gap[post]; gc <- cumsum(g); s <- tr$synth[post]; sc <- cumsum(s)
  na_pre <- rep(NA_real_, sum(!post))
  col <- function(v) c(na_pre, v)
  cov_adj <- if (is.null(adj_main)) rep(0, nrow(periods)) else
    as.numeric(Reduce(`+`, lapply(adjust_covariates, function(v) adj_main$gamma[[v]] * Aobs[[v]][, treated_unit])))
  s <- s + cov_adj[post]; sc <- cumsum(s)                        # counterfactual on the original scale
  effects <- cbind(periods, data.frame(
    post = post, observed = Y[, treated_unit], synthetic = Y[, treated_unit] - tr$gap, effect = tr$gap,
    covariate_adjustment = cov_adj,
    halfwidth = col(tr$h), factor_pointwise = col(k_pw), factor_simultaneous = col(sim$factor),
    eff_low_pw = col(g - k_pw * tr$h), eff_high_pw = col(g + k_pw * tr$h),
    eff_low_sim = col(g - sim$factor * tr$h), eff_high_sim = col(g + sim$factor * tr$h),
    pred_low_pw = col(s - k_pw * tr$h), pred_high_pw = col(s + k_pw * tr$h),
    pred_low_sim = col(s - sim$factor * tr$h), pred_high_sim = col(s + sim$factor * tr$h),
    p_pointwise = col((1 + rowSums(Kp >= tr$k)) / (N + 1)),
    cum_effect = col(gc), cum_halfwidth = col(tr$hc),
    cum_factor_pointwise = col(k_pwc), cum_factor_simultaneous = col(simc$factor),
    cum_eff_low_pw = col(gc - k_pwc * tr$hc), cum_eff_high_pw = col(gc + k_pwc * tr$hc),
    cum_eff_low_sim = col(gc - simc$factor * tr$hc), cum_eff_high_sim = col(gc + simc$factor * tr$hc),
    cum_pred_low_pw = col(sc - k_pwc * tr$hc), cum_pred_high_pw = col(sc + k_pwc * tr$hc),
    cum_pred_low_sim = col(sc - simc$factor * tr$hc), cum_pred_high_sim = col(sc + simc$factor * tr$hc),
    cum_p_pointwise = col((1 + rowSums(Kpc >= tr$kc)) / (N + 1))))

  # ---- Tests: pointwise at test_period; simultaneous over the whole post period -----
  e_post <- effects[post, ]
  test_label <- paste(names(periods), unlist(periods[t_treat + t_test - 1, ]), sep = " = ", collapse = ", ")
  flagged <- function(lo, hi) {
    hit <- which(lo > 0 | hi < 0)
    if (length(hit)) paste(apply(periods[post, , drop = FALSE][hit, , drop = FALSE], 1, paste, collapse = "-"), collapse = "; ") else "none"
  }
  test <- data.frame(
    test = c("pointwise", "pointwise", "simultaneous", "simultaneous"),
    measure = c("effect", "cumulative effect", "effect", "cumulative effect"),
    period = c(test_label, paste("through", test_label), "whole post period", "whole post period"),
    estimate = c(e_post$effect[t_test], e_post$cum_effect[t_test], NA, NA),
    lower = c(e_post$eff_low_pw[t_test], e_post$cum_eff_low_pw[t_test], NA, NA),
    upper = c(e_post$eff_high_pw[t_test], e_post$cum_eff_high_pw[t_test], NA, NA),
    p_value = c(e_post$p_pointwise[t_test], e_post$cum_p_pointwise[t_test], sim$p, simc$p),
    significant = c(e_post$eff_low_pw[t_test] > 0 | e_post$eff_high_pw[t_test] < 0,
                    e_post$cum_eff_low_pw[t_test] > 0 | e_post$cum_eff_high_pw[t_test] < 0,
                    sim$p <= alpha, simc$p <= alpha),
    periods_flagged = c(NA, NA, flagged(e_post$eff_low_sim, e_post$eff_high_sim),
                        flagged(e_post$cum_eff_low_sim, e_post$cum_eff_high_sim)))

  adjustment <- NULL
  if (!is.null(adj_main)) {
    point_at <- function(gam) {                                 # point estimates only (no bootstrap)
      Ya <- adjust_mat(Y[, c(treated_unit, donors), drop = FALSE], Aobs, gam)
      gp <- (Ya[, treated_unit] - fit_dsc(Ya[, donors, drop = FALSE], Ya[, treated_unit])$synth)[post]
      c(effect = gp[t_test], cumulative = sum(gp[seq_len(t_test)]))
    }
    unadj <- { gp <- (Y[, treated_unit] - fit_dsc(Y[, donors, drop = FALSE], Y[, treated_unit])$synth)[post]
               c(effect = gp[t_test], cumulative = sum(gp[seq_len(t_test)])) }
    adjustment <- list(gamma = adj_main$gamma, se = adj_main$se, n_obs = adj_main$n_obs,
                       estimates = rbind(unadjusted = unadj, adjusted = point_at(adj_main$gamma),
                                         gamma_lower = point_at(adj_main$gamma - 1.96 * adj_main$se),
                                         gamma_upper = point_at(adj_main$gamma + 1.96 * adj_main$se)))
  }

  list(effects = effects, test = test, adjustment = adjustment,
       weights = sort(tr$w[tr$w > 0], decreasing = TRUE), intercept = tr$b0, control_coefs = tr$g,
       pre_rmspe = sqrt(mean(tr$gap[!post]^2)),
       donor_factors = list(pointwise = Kp, cumulative = Kpc),
       placebo = list(gaps = setNames(as.data.frame(sapply(pl, `[[`, "gap")), donors[!failed]),        # all periods x donors
                      synthetic = setNames(as.data.frame(sapply(pl, `[[`, "synth")), donors[!failed]),
                      halfwidth = setNames(as.data.frame(matrix(sapply(pl, `[[`, "h"), nrow = n_post)), donors[!failed]),
                      cum_halfwidth = setNames(as.data.frame(matrix(sapply(pl, `[[`, "hc"), nrow = n_post)), donors[!failed])),
       failed_donors = failed_donors,
       settings = list(treated_unit = treated_unit, treatment_period = treatment_period, test_period = test_label,
                       time_cols = time_cols, outcome = outcome_col, controls = controls,
                       adjust_covariates = adjust_covariates,
                       control_units = if (length(controls)) control_units else NA,
                       conf_level = conf_level, reps = reps, trim = trim, span = span, seed = seed,
                       n_donors = N, rank_used = m))
}

## ===========================================================================
## 5. plot_dsc_symcal()
## ===========================================================================

##############################################################################
# plot_dsc_symcal(): ggplot figures for a dsc_symcal_boot() result
#   prediction   observed vs synthetic, with prediction bands
#   effect       effect by period, with effect bands
#   cumulative   cumulative effect, with effect bands
# Shaded: simultaneous band; dashed: pointwise band; the pre-specified test
# period is marked with a point and a vertical line.
# Arguments: res (dsc_symcal_boot() result), y_label (axis label).
# Value: a list of three ggplot objects.
##############################################################################
plot_dsc_symcal <- function(res, y_label = res$settings$outcome) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("plot_dsc_symcal() needs the ggplot2 package.")
  e <- res$effects
  e$x <- seq_len(nrow(e))
  labs_x <- apply(e[, res$settings$time_cols, drop = FALSE], 1, paste, collapse = "-")
  brks <- unique(round(seq(1, nrow(e), length.out = min(nrow(e), 8))))
  t0 <- which(e$post)[1]
  ep <- e[e$post, ]
  tt <- which(e$post)[1] - 1 + which(apply(ep[, res$settings$time_cols, drop = FALSE], 1, function(r)
    paste(names(r), r, sep = " = ", collapse = ", ")) == res$settings$test_period)
  lvl <- paste0(100 * res$settings$conf_level, "%")
  note <- paste0("Shaded: simultaneous ", lvl, " band; dashed: pointwise ", lvl, " band; marked: test period")
  common <- list(ggplot2::geom_vline(xintercept = t0 - 0.5, linetype = "dotted", colour = "grey40"),
                 ggplot2::geom_vline(xintercept = tt, linetype = "dashed", colour = "grey60"),
                 ggplot2::scale_x_continuous(breaks = brks, labels = labs_x[brks]),
                 ggplot2::theme_minimal())
  band_plot <- function(est, lo_s, hi_s, lo_p, hi_p, fill, title, ylab) {
    ggplot2::ggplot(ep, ggplot2::aes(x)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lo_s]], ymax = .data[[hi_s]]), fill = fill, alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(y = .data[[lo_p]]), linetype = "dashed", colour = fill) +
      ggplot2::geom_line(ggplot2::aes(y = .data[[hi_p]]), linetype = "dashed", colour = fill) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
      ggplot2::geom_line(ggplot2::aes(y = .data[[est]])) +
      ggplot2::geom_point(data = e[tt, ], ggplot2::aes(x, .data[[est]]), size = 2.5) +
      common + ggplot2::labs(title = title, subtitle = note, x = NULL, y = ylab)
  }
  p_pred <- ggplot2::ggplot(e, ggplot2::aes(x)) +
    ggplot2::geom_ribbon(data = ep, ggplot2::aes(ymin = pred_low_sim, ymax = pred_high_sim), fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(data = ep, ggplot2::aes(y = pred_low_pw), linetype = "dashed", colour = "steelblue4") +
    ggplot2::geom_line(data = ep, ggplot2::aes(y = pred_high_pw), linetype = "dashed", colour = "steelblue4") +
    ggplot2::geom_line(ggplot2::aes(y = synthetic, colour = "Synthetic")) +
    ggplot2::geom_line(ggplot2::aes(y = observed, colour = "Observed")) +
    ggplot2::scale_colour_manual(values = c(Observed = "black", Synthetic = "steelblue4"), name = NULL) +
    common + ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(title = paste0(res$settings$treated_unit, ": observed vs synthetic, with prediction intervals"),
                  subtitle = note, x = NULL, y = y_label)
  list(prediction = p_pred,
       effect = band_plot("effect", "eff_low_sim", "eff_high_sim", "eff_low_pw", "eff_high_pw", "firebrick",
                          paste0(res$settings$treated_unit, ": effect by period"), paste("Effect,", y_label)),
       cumulative = band_plot("cum_effect", "cum_eff_low_sim", "cum_eff_high_sim", "cum_eff_low_pw", "cum_eff_high_pw", "darkgreen",
                              paste0(res$settings$treated_unit, ": cumulative effect"), paste("Cumulative effect,", y_label)))
}



## ===========================================================================
## 6. dsc_placebo()
## ===========================================================================

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
#   ridge_grid), estimator, diagnostics / warn / ife_seed - described below and in
#   the reference.
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
#   Experimental_inference/ridge_default_validation.R: as valid as cross-validation in
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
dsc_placebo <- function(data, unit_col, time_cols, outcome_col, treated_unit, treatment_period,
                        test_period = NULL, conf_level = 0.95, span = 0.75,
                        scale = c("rmspe", "relative", "raw"),
                        adjust_covariates = NULL, adjust_model = c("ife", "twfe"), ife_r = 0:5,
                        augment = c("none", "ife", "ridge"), diagnostics = TRUE, warn = TRUE, ife_seed = 12112025,
                        estimator = c("dsc", "fect"), rmspe_floor = 0.25,
                        ridge_c = 0.1, ridge_grid = 10^seq(-3, 2, by = 0.5)) {

  if (!requireNamespace("CVXR", quietly = TRUE)) stop("dsc_placebo() needs the CVXR package.")
  scale <- match.arg(scale); adjust_model <- match.arg(adjust_model); augment <- match.arg(augment)
  estimator <- match.arg(estimator)
  if (estimator == "fect" && augment != "none") stop("augment = \"", augment, "\" applies to estimator = \"dsc\" only.")
  if (identical(ridge_c, "cv")) ridge_c <- NULL                  # "cv" (or NULL): choose the penalty by cross-validation
  if (!is.null(ridge_c) && (!is.numeric(ridge_c) || length(ridge_c) != 1 || !is.finite(ridge_c) || ridge_c <= 0)) stop("ridge_c must be one positive number, or \"cv\" to choose it by cross-validation.")
  if (estimator == "fect" && !requireNamespace("fect", quietly = TRUE)) stop("estimator = \"fect\" needs the fect package.")
  if (augment == "ife" && !requireNamespace("fect", quietly = TRUE)) stop("augment = \"ife\" needs the fect package.")
  if (length(adjust_covariates) && adjust_model == "ife" && !requireNamespace("fect", quietly = TRUE)) stop("adjust_model = \"ife\" needs the fect package.")
  alpha <- round(1 - conf_level, 10)

  # ---- Data -------------------------------------------------------------------------
  data <- as.data.frame(data)
  stopifnot(all(c(unit_col, time_cols, outcome_col, adjust_covariates) %in% names(data)), all(names(treatment_period) %in% time_cols))
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
  for (v in adjust_covariates) d[[v]] <- data[[v]]
  if (anyNA(d[, c("y", adjust_covariates)]) || any(table(d$unit) != nrow(periods))) stop("The panel must be balanced with no missing outcomes or covariates.")
  if (!treated_unit %in% d$unit) stop("treated_unit not found in ", unit_col, ".")
  units <- sort(unique(d$unit)); donors <- setdiff(units, treated_unit)
  wide <- function(v) sapply(units, function(u) { du <- d[d$unit == u, ]; du[[v]][order(du$.t)] })
  Y <- wide("y")
  Aobs <- setNames(lapply(adjust_covariates, wide), adjust_covariates)

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
    r_use <- model_ife(drop_post_unit = treated_unit)$r                  # cross-validation chooses the number of factors
    aug_main <- model_ife(drop_post_unit = treated_unit, r = r_use, cv = FALSE)   # refit with it, exactly as for the placebos
    adj_main <- list(gamma = aug_main$gamma, n_obs = aug_main$n_obs)
    if (length(adjust_covariates)) tw_main <- gamma_twfe(drop_post_unit = treated_unit)
  } else if (length(adjust_covariates)) {
    tw_main <- gamma_twfe(drop_post_unit = treated_unit)
    if (adjust_model == "ife") { r_use <- gamma_ife(drop_post_unit = treated_unit)$r
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
  pl <- lapply(donors, function(p) tryCatch({
    if (estimator == "fect")                                                                 # IFEct, donor p as pseudo-treated
      return(list(gap = Y[, p] - model_ife(drop_units = treated_unit, drop_post_unit = p, r = r_use, cv = FALSE)$M[, p]))
    cols_p <- c(p, setdiff(donors, p))
    Yp <- if (augment == "ife") Y[, cols_p, drop = FALSE] - model_ife(drop_units = treated_unit, drop_post_unit = p, r = r_use, cv = FALSE)$M[, cols_p, drop = FALSE] else
      adjust_mat(cols_p, gamma_for_placebo(p))
    fp <- fit_dsc(Yp[, setdiff(donors, p), drop = FALSE], Yp[, p])
    corr_p <- if (augment == "ridge") ridge_correction(Yp[, p], Yp[, setdiff(donors, p), drop = FALSE], fp$w, c_ridge) else 0
    list(gap = Yp[, p] - fp$synth - corr_p)
  }, error = function(e) conditionMessage(e)))
  failed <- vapply(pl, is.character, logical(1))
  failed_donors <- if (any(failed)) data.frame(donor = donors[failed], error = unlist(pl[failed])) else
    data.frame(donor = character(0), error = character(0))
  if (any(failed)) warning(sum(failed), " placebo fit(s) failed and were left out: ", paste(donors[failed], collapse = ", "), ". See $failed_donors.")
  pl <- pl[!failed]; dn <- donors[!failed]; N <- length(pl); m <- floor(alpha * (N + 1))
  kth <- function(x) if (m == 0) Inf else sort(x, decreasing = TRUE)[m]
  if (m == 0) warning("With ", N, " donors the ", 100 * conf_level, "% level cannot be reached exactly: intervals are infinite.")

  GP <- sapply(pl, `[[`, "gap"); colnames(GP) <- dn; SP <- Y[, dn, drop = FALSE] - GP                            # gaps and counterfactuals, all periods x donors
  synth_T <- Y[, treated_unit] - gT_all                                                         # counterfactual on the original scale
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
      tryCatch({ r_d <- model_ife(drop_post_unit = treated_unit)$r; model_ife(drop_post_unit = treated_unit, r = r_d, cv = FALSE) },
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

  list(effects = effects, test = test,
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
       failed_donors = failed_donors,
       settings = list(treated_unit = treated_unit, treatment_period = treatment_period, test_period = test_label,
                       time_cols = time_cols, outcome = outcome_col, conf_level = conf_level, span = span,
                       scale = scale, rmspe_floor = if (scale == "rmspe") rmspe_floor else NA, adjust_covariates = adjust_covariates,
                       adjust_model = if (length(adjust_covariates)) adjust_model else NA, augment = augment, estimator = estimator,
                       ridge_c = if (augment == "ridge") c_ridge else NA, n_donors = N, rank_used = m))
}



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
    stop("`fit` must be the output of dsc_symcal_boot() or dsc_placebo().")
}

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

added_effect <- function(fit_real, fit_constructed) {
  .dsc_check_fit(fit_real); .dsc_check_fit(fit_constructed)
  if (!identical(fit_real$effects$post, fit_constructed$effects$post)) stop("The two fits must cover the same periods.")
  (fit_constructed$effects$observed - fit_real$effects$observed)[fit_real$effects$post]
}

