##############################################################################
# california_estimator_accuracy_test.R
#
# TEST (step 1): accuracy of synthetic control point estimates on the California
# smoking panel (tidysynth::smoking; California dropped; each of the 38 other
# states as a fake treated unit in 1989, with its 37 donors). No effect is
# injected, so each state's gap is its estimation error (= bias at any effect
# size, since the estimate shifts exactly with an injected effect).
#
# Variants (all with an intercept unless noted):
#   lasso_01            LASSO, weights in [0, 1]                      (current)
#   sum1_nonneg         weights >= 0 summing to 1 ("demeaned" SC; no penalty
#                       needed because the L1 norm is fixed at 1; CVXR)
#   lasso_nonneg        LASSO, weights >= 0, no upper bound
#   lasso_free          LASSO, no weight constraints (SCUL-style)
#   ridge_01            ridge, weights in [0, 1]
#   ridge_free          ridge, no weight constraints
#   enet_free           elastic net (alpha = 0.5), no weight constraints
#   lasso_01_noint      LASSO, weights in [0, 1], NO intercept
# Penalty rules (glmnet variants): min_path (smallest lambda on the path),
#   cv_min (glmnet 10-fold CV lambda.min), rolling (rolling-origin CV, SCUL
#   median rule; validation window min(post, T0/3), first window T0/2).
#
# Diagnostic: for the current model, relation between the sum of the weights
# and each state's error.
#
# Output: results_california_accuracy/accuracy_results.rds; tables at the end
# of this script.
##############################################################################


# Setup
suppressPackageStartupMessages({ library(dplyr); library(CVXR) })
source("synthetic_control_functions/synthetic_control_functions_claude.R")
script_file <- "simulations/california_estimator_accuracy_test.R"
out_dir <- file.path("simulations", "results_california_accuracy")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
treat_year <- 1989; seed <- 12112025

sm <- as.data.frame(tidysynth::smoking)[, c("state", "year", "cigsale")]
stopifnot(!anyNA(sm), length(unique(table(sm$state))) == 1)
pl <- sm[sm$state != "California", ]
states <- sort(unique(pl$state)); years <- sort(unique(pl$year)); post <- years >= treat_year
Y <- sapply(states, function(u) { du <- pl[pl$state == u, ]; du$cigsale[order(du$year)] })

variants <- list(
  lasso_01       = list(alpha = 1,   lower = 0,    upper = 1,   intercept = TRUE),
  lasso_nonneg   = list(alpha = 1,   lower = 0,    upper = Inf, intercept = TRUE),
  lasso_free     = list(alpha = 1,   lower = -Inf, upper = Inf, intercept = TRUE),
  ridge_01       = list(alpha = 0,   lower = 0,    upper = 1,   intercept = TRUE),
  ridge_free     = list(alpha = 0,   lower = -Inf, upper = Inf, intercept = TRUE),
  enet_free      = list(alpha = 0.5, lower = -Inf, upper = Inf, intercept = TRUE),
  lasso_01_noint = list(alpha = 1,   lower = 0,    upper = 1,   intercept = FALSE))

gfit <- function(X, y, v, nl) glmnet::glmnet(X, y, alpha = v$alpha, standardize = FALSE, intercept = v$intercept,
                                             lower.limits = v$lower, upper.limits = v$upper, nlambda = nl)
rolling_lambda <- function(X, y, v) {
  T0 <- nrow(X); h <- min(sum(post), floor(T0 / 3)); init <- floor(T0 / 2)
  median(vapply(init:(T0 - h), function(i) {
    f <- gfit(X[1:i, , drop = FALSE], y[1:i], v, 100)
    pr <- predict(f, newx = X[(i + 1):(i + h), , drop = FALSE])
    f$lambda[which.min(colMeans((pr - y[(i + 1):(i + h)])^2))]
  }, numeric(1)))
}
summ_fit <- function(unit, variant, rule, b0, w, X, y, lambda = NA) {
  gap <- y - as.numeric(b0 + X %*% w)
  data.frame(unit = unit, variant = variant, rule = rule, lambda = lambda,
             err_final = 100 * tail(gap, 1) / tail(y, 1), err_cum = 100 * sum(gap[post]) / sum(y[post]),
             pre_rmspe = sqrt(mean(gap[!post]^2)), weight_sum = sum(w), n_donors = sum(abs(w) > 1e-8))
}


# Fit every variant and rule for every state
t0 <- Sys.time()
res <- do.call(rbind, lapply(states, function(j) {
  X <- Y[, setdiff(states, j)]; y <- Y[, j]; Xp <- X[!post, ]; yp <- y[!post]
  out <- list()
  for (vn in names(variants)) {
    v <- variants[[vn]]
    f <- gfit(Xp, yp, v, 1000)
    cvl <- withCallingHandlers({ set.seed(seed)
      glmnet::cv.glmnet(Xp, yp, alpha = v$alpha, standardize = FALSE, intercept = v$intercept,
                        lower.limits = v$lower, upper.limits = v$upper, nlambda = 1000)$lambda.min },
      warning = function(w) if (grepl("grouped=FALSE", conditionMessage(w))) invokeRestart("muffleWarning"))
    lams <- c(min_path = min(f$lambda), cv_min = cvl, rolling = rolling_lambda(Xp, yp, v))
    for (r in names(lams)) {
      cf <- as.matrix(coef(f, s = lams[[r]]))
      out[[length(out) + 1]] <- summ_fit(j, vn, r, cf[1, 1], cf[-1, 1], X, y, lams[[r]])
    }
  }
  # Demeaned SC: intercept + non-negative weights summing to one (CVXR)
  w <- Variable(ncol(Xp)); b0 <- Variable(1)
  sol <- tryCatch(solve(Problem(Minimize(sum_squares(yp - b0 - Xp %*% w)), list(w >= 0, sum(w) == 1))),
                  error = function(e) NULL)
  if (!is.null(sol) && sol$status %in% c("optimal", "optimal_inaccurate")) {
    wv <- as.numeric(sol$getValue(w)); wv[abs(wv) < 1e-8] <- 0
    out[[length(out) + 1]] <- summ_fit(j, "sum1_nonneg", "none", as.numeric(sol$getValue(b0)), wv, X, y)
  }
  do.call(rbind, out)
}))
run_minutes <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)


# Check: baseline (lasso_01, min_path) reproduces synth_lasso_boot() for the first state
chk_state <- states[1]
m <- suppressWarnings(synth_lasso_boot(pl, "state", "year", "cigsale", chk_state, treatment_period = c(year = treat_year),
                                       seed = seed, boot = FALSE))
base_chk <- res[res$unit == chk_state & res$variant == "lasso_01" & res$rule == "min_path", ]
check_diff <- abs(base_chk$err_final - 100 * tail(m$effects$att, 1) / tail(Y[, chk_state], 1))


# Accuracy summary (one error per state)
acc <- res |> group_by(variant, rule) |>
  summarise(n_states = n(),
            bias_final = mean(err_final), bias_final_lo = t.test(err_final)$conf.int[1], bias_final_hi = t.test(err_final)$conf.int[2],
            rmse_final = sqrt(mean(err_final^2)), mae_final = mean(abs(err_final)),
            bias_cum = mean(err_cum), bias_cum_lo = t.test(err_cum)$conf.int[1], bias_cum_hi = t.test(err_cum)$conf.int[2],
            rmse_cum = sqrt(mean(err_cum^2)), mae_cum = mean(abs(err_cum)),
            median_pre_rmspe = median(pre_rmspe), median_weight_sum = median(weight_sum),
            median_donors = median(n_donors), .groups = "drop") |>
  arrange(rmse_cum) |> as.data.frame()

# Diagnostic: weight sum vs error for the current model
base <- res[res$variant == "lasso_01" & res$rule == "min_path", ]
diag <- list(
  weight_sum_summary = summary(base$weight_sum),
  cor_weight_sum_err_final = cor.test(base$weight_sum, base$err_final)[c("estimate", "p.value")],
  cor_weight_sum_err_cum = cor.test(base$weight_sum, base$err_cum)[c("estimate", "p.value")],
  by_group = base |> mutate(group = ifelse(weight_sum < 0.9, "sum < 0.9", "sum >= 0.9")) |> group_by(group) |>
    summarise(n = n(), mean_err_final = mean(err_final), mean_err_cum = mean(err_cum), .groups = "drop") |> as.data.frame())


# Save and write the tables at the end of this script
saveRDS(list(accuracy = acc, per_state = res, diagnostic = diag, check_diff = check_diff, run_minutes = run_minutes),
        file.path(out_dir, "accuracy_results.rds"))
disp <- acc |> transmute(variant, rule,
  bias_final = sprintf("%.2f [%.2f, %.2f]", bias_final, bias_final_lo, bias_final_hi), rmse_final = round(rmse_final, 2),
  bias_cum = sprintf("%.2f [%.2f, %.2f]", bias_cum, bias_cum_lo, bias_cum_hi), rmse_cum = round(rmse_cum, 2),
  pre_rmspe = round(median_pre_rmspe, 2), weight_sum = round(median_weight_sum, 2), donors = median_donors)
marker <- "# ---- Results (written automatically when this script is run) ----"
script_lines <- readLines(script_file)
cut_at <- which(script_lines == marker)
if (length(cut_at) > 0) script_lines <- script_lines[seq_len(cut_at[1] - 1)]
writeLines(c(script_lines, marker,
  sprintf("# Run: %s | %d states | %s min | baseline check vs synth_lasso_boot(): %.3g",
          format(Sys.time(), "%Y-%m-%d %H:%M"), length(states), run_minutes, check_diff),
  "# Errors in % of the untreated outcome; sorted by cumulative RMSE.",
  paste0("# ", knitr::kable(disp, format = "simple"))), script_file)

# ---- Results (written automatically when this script is run) ----
# Run: 2026-10-01 15:30 | 38 states | 2.1 min | baseline check vs synth_lasso_boot(): 0.00326
# Errors in % of the untreated outcome; sorted by cumulative RMSE.
# variant          rule       bias_final                rmse_final  bias_cum                 rmse_cum   pre_rmspe   weight_sum   donors
# ---------------  ---------  -----------------------  -----------  ----------------------  ---------  ----------  -----------  -------
# lasso_01_noint   min_path   -1.37 [-6.15, 3.40]            14.41  1.95 [-1.10, 5.00]           9.36        8.37         0.60      2.0
# lasso_01_noint   cv_min     -1.36 [-6.14, 3.42]            14.41  1.96 [-1.09, 5.01]           9.37        8.37         0.60      2.0
# lasso_01_noint   rolling    -0.28 [-4.97, 4.42]            14.10  2.88 [-0.16, 5.93]           9.59        8.45         0.59      2.0
# sum1_nonneg      none       -0.89 [-6.11, 4.33]            15.70  0.14 [-3.14, 3.42]           9.85        1.94         1.00      6.0
# ridge_01         min_path   -6.49 [-12.25, -0.74]          18.46  -3.41 [-6.96, 0.14]         11.19        2.80         0.75     31.0
# ridge_01         cv_min     -6.52 [-12.28, -0.75]          18.48  -3.43 [-6.98, 0.12]         11.20        2.80         0.75     31.0
# ridge_01         rolling    -8.80 [-14.42, -3.17]          19.03  -4.64 [-8.03, -1.24]        11.20        2.92         0.70     32.5
# lasso_01         cv_min     -7.83 [-13.82, -1.83]          19.63  -3.40 [-7.18, 0.38]         11.85        2.10         0.73      4.0
# lasso_01         min_path   -4.95 [-11.29, 1.39]           19.67  -1.57 [-5.51, 2.36]         11.92        1.78         0.79      5.0
# lasso_nonneg     cv_min     -7.99 [-13.96, -2.01]          19.65  -3.32 [-7.15, 0.51]         11.97        2.10         0.73      4.0
# ridge_free       min_path   -9.55 [-15.33, -3.77]          19.81  -5.50 [-9.06, -1.94]        12.02        2.37         0.68     37.0
# ridge_free       cv_min     -9.57 [-15.36, -3.79]          19.83  -5.52 [-9.08, -1.96]        12.03        2.37         0.68     37.0
# lasso_nonneg     min_path   -5.11 [-11.44, 1.22]           19.69  -1.50 [-5.47, 2.48]         12.04        1.78         0.79      5.0
# enet_free        cv_min     -11.08 [-17.02, -5.14]         20.99  -5.53 [-9.37, -1.69]        12.79        1.55         0.64      6.0
# lasso_01         rolling    -9.54 [-15.87, -3.21]          21.27  -4.31 [-8.33, -0.29]        12.81        2.24         0.71      4.0
# lasso_nonneg     rolling    -9.73 [-16.01, -3.45]          21.22  -4.27 [-8.31, -0.24]        12.85        2.24         0.71      4.0
# lasso_free       cv_min     -11.29 [-17.24, -5.34]         21.13  -5.71 [-9.57, -1.85]        12.92        1.55         0.64      5.5
# lasso_free       min_path   -10.08 [-16.26, -3.90]         21.12  -5.09 [-9.20, -0.98]        13.35        1.45         0.71      7.0
# enet_free        min_path   -10.13 [-16.35, -3.92]         21.24  -5.12 [-9.24, -1.00]        13.39        1.45         0.71      7.0
# ridge_free       rolling    -12.06 [-18.04, -6.09]         21.62  -7.14 [-10.93, -3.35]       13.43        2.86         0.58     37.0
# lasso_free       rolling    -15.67 [-22.09, -9.25]         24.84  -8.36 [-12.46, -4.26]       14.88        2.29         0.60      4.0
# enet_free        rolling    -15.72 [-22.17, -9.27]         24.94  -8.45 [-12.56, -4.35]       14.95        2.31         0.59      4.0
