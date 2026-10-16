# Function reference

*Written by Claude, an AI model made by Anthropic, at the direction of Jeffrey Milliman.*

All functions are in [`synthetic_control_functions/synthetic_control_functions_claude.R`](../../synthetic_control_functions/synthetic_control_functions_claude.R). Load them with

```r
source("synthetic_control_functions/synthetic_control_functions_claude.R")
```

- [`dsc_placebo()`](#dsc_placebo): the main method
  - [Arguments](#arguments)
  - [Value](#value)
  - [Errors, warnings and limits](#errors-warnings-and-limits)
- [Plotting functions](#plotting-functions): `plot_full_path()`, `plot_post_effect()`, `plot_prediction()`, `added_effect()`
- [Earlier functions](#earlier-functions)

---

## dsc_placebo()

Demeaned synthetic control for a single treated unit, with placebo-calibrated pointwise intervals, simultaneous bands, prediction intervals and diagnostics.

```r
dsc_placebo(data, unit_col, time_cols, outcome_col, treated_unit, treatment_period,
            test_period = NULL, conf_level = 0.95, span = 0.75,
            scale = c("rmspe", "relative", "raw"),
            adjust_covariates = NULL, adjust_model = c("ife", "twfe"), ife_r = 0:5,
            augment = c("none", "ife"), diagnostics = TRUE, warn = TRUE, ife_seed = 12112025,
            estimator = c("dsc", "fect"), rmspe_floor = 0.25)
```

**What it does.** The counterfactual is an intercept plus non-negative donor weights that sum to one, fitted to the treated unit's pre-treatment outcomes. Every donor is then fitted the same way as if it had been treated (a *placebo*), with the treated unit removed from the pool. Each unit's post-treatment gaps are put on a common scale (by default, divided by the unit's own pre-treatment RMSPE), and the treated unit's gaps are ranked among the placebos'. With $N$ donors and $m = \lfloor (1 - \text{conf\_level})(N + 1) \rfloor$:

- **pointwise**: the threshold in each period is the $m$-th largest placebo value; the interval is the gap $\pm$ threshold $\times$ the treated unit's scale;
- **simultaneous**: a smoothed profile of the median placebo value over time, and one multiplier (the $m$-th largest placebo maximum relative to that profile), so that the band covers the whole post-treatment path at once.

The construction is explained, with a worked example, in [Placebo-calibrated intervals and simultaneous bands](../placebo/placebo_calibrated_bands.md).

### Arguments

**Data**

| Argument | Description |
|---|---|
| `data` | A data frame in long format: one row per unit and period. The panel must be balanced (every unit observed in every period) with no missing values in the outcome or in `adjust_covariates`. |
| `unit_col` | Name of the column identifying units (character). |
| `time_cols` | Name(s) of the column(s) identifying periods, from the coarsest to the finest, e.g. `"year"` or `c("year", "qtr")`. Periods are ordered by these columns in turn. |
| `outcome_col` | Name of the outcome column (numeric). |
| `treated_unit` | The treated unit: a value of `unit_col`. All other units are donors. |
| `treatment_period` | The first treated period, as a named vector with one element per time column, e.g. `c(year = 1991)` or `c(year = 2012, qtr = 2)`. Every earlier period is pre-treatment. |

**Inference**

| Argument | Default | Description |
|---|---|---|
| `test_period` | `NULL` | The period for the pointwise tests, in the same format as `treatment_period`. `NULL` uses the last post-treatment period. The per-period test is for that period, and the cumulative test for the effect summed from the start of treatment through that period. Choose it before looking at the results. |
| `conf_level` | `0.95` | Confidence level for the intervals and bands. With $N$ donors the smallest possible p-value is $1/(N + 1)$: a 95% interval needs at least 19 donors and a 90% interval at least 9; with fewer, the intervals are infinite. |
| `span` | `0.75` | Loess span for smoothing the simultaneous band's time profile (on the log scale). Larger values give a smoother profile. |
| `scale` | `"rmspe"` | What each unit's gaps are divided by before ranking. `"rmspe"`: the unit's pre-treatment RMSPE (gaps measured in units of the unit's own pre-treatment fit; the most robust when the treated unit is unlike the donors). `"relative"`: the unit's synthetic value, i.e. percentage gaps (the most powerful when the treated unit resembles the donors; needs positive outcomes). `"raw"`: nothing (gaps in the outcome's units). |
| `rmspe_floor` | `0.25` | Guard for `scale = "rmspe"`: each unit's pre-treatment RMSPE is floored at `rmspe_floor` $\times$ the median pre-treatment RMSPE over the treated unit and all placebos, so that a unit reproduced (almost) exactly before treatment cannot produce near-infinite intervals. Matters mainly with short pre-treatment periods and many donors. `0` switches it off. Ignored for the other scales. |

**Covariate adjustment** (for measured time-varying confounders)

| Argument | Default | Description |
|---|---|---|
| `adjust_covariates` | `NULL` | Name(s) of time-varying covariate columns to adjust for. Their effect $\gamma$ is estimated on untreated cells only (all donor periods and the treated unit's pre-treatment periods); every unit's outcome is then replaced by outcome $- \hat\gamma \times$ covariate, and the method runs on the adjusted outcomes. Each placebo re-estimates $\gamma$ without its own post-treatment period. Assumes the treatment does not affect the covariates and that their effect is common and linear. See the [Oregon case study](../oregon_fentanyl/oregon_fentanyl_case_study.md). |
| `adjust_model` | `"ife"` | Model used to estimate $\gamma$: `"ife"`, interactive fixed effects (unit and time effects plus latent factors, via `fect`), or `"twfe"`, two-way fixed effects. |
| `ife_r` | `0:5` | Candidate numbers of latent factors for every interactive fixed-effects model (covariate adjustment, `augment = "ife"`, `estimator = "fect"` and the diagnostics). The number is chosen by cross-validation for the treated unit and reused for every placebo. A single value skips cross-validation. |

**Estimator**

| Argument | Default | Description |
|---|---|---|
| `estimator` | `"dsc"` | `"dsc"`: the demeaned synthetic control. `"fect"`: the counterfactual is `fect`'s interactive fixed-effects prediction, fitted on untreated cells, with each donor refitted as a pseudo-treated unit and the same placebo calibration (use with `scale = "rmspe"`). Recommended by the diagnostics only for a treated unit beyond the donors' range with stationary factors. |
| `augment` | `"none"` | `"ife"`: augmented demeaned synthetic control. An interactive fixed-effects outcome model is fitted on untreated cells, the synthetic control is fitted to its residuals, and the counterfactual is the model prediction plus the synthetic control of the residuals. Reduces bias for a treated unit outside the donors' range, but raised variance and made intervals overconfident in simulation; not combined with `estimator = "fect"`. |

**Diagnostics and reproducibility**

| Argument | Default | Description |
|---|---|---|
| `diagnostics` | `TRUE` | Check whether the treated unit lies beyond what the donors can reproduce (pre-treatment fit rank, distance from the donors' convex hull in estimated factor loadings, implied extrapolation drift, factor persistence). Loading-based checks need `fect` and take most of the run time; `FALSE` skips them. See [When placebo inference fails](../diagnostics/diagnostics_toy_cases.md). |
| `warn` | `TRUE` | Issue a warning when the diagnostics flag the treated unit as extreme. |
| `ife_seed` | `12112025` | Seed for `fect`'s cross-validation, applied locally (the session's random-number state is restored afterwards), so results with interactive fixed effects are reproducible. |

### Value

A list:

| Element | Contents |
|---|---|
| `effects` | Data frame, one row per period (all periods; post-treatment rows have `post = TRUE`). Columns below. |
| `test` | Data frame with five rows: the pointwise test of the effect at `test_period`, the pointwise test of the cumulative effect through `test_period`, the simultaneous tests of the effect and of the cumulative effect over the whole post-treatment period, and the classic post/pre RMSPE-ratio test (Abadie et al.). Columns: `test`, `measure`, `period`, `estimate`, `lower`, `upper` (the interval; pointwise tests only), `p_value`, `significant` (the interval excludes zero, or for the simultaneous and ratio tests `p_value` $\le 1 -$ `conf_level`), `periods_flagged` (simultaneous tests: the periods where the band excludes zero). |
| `weights` | Named vector of non-zero donor weights, largest first (`NULL` with `estimator = "fect"`). |
| `intercept` | The intercept of the synthetic control (`NULL` with `estimator = "fect"`). |
| `pre_rmspe` | The treated unit's pre-treatment RMSPE (unfloored). |
| `rmspe_guard` | `floor`: the floor applied to the pre-treatment RMSPEs (0 if none); `units_floored`: the units it applied to. |
| `adjustment` | With `adjust_covariates`: `model`, `gamma` (the estimates), `n_factors`, `gamma_twfe` and `se_twfe` (two-way fixed-effects estimate with a unit-clustered standard error, for comparison), `n_obs`, and `unadjusted` (the effect and cumulative effect at `test_period` without the adjustment). `NULL` otherwise. |
| `placebo` | `gaps` and `synthetic`: data frames (periods $\times$ donors) of each placebo's gaps and counterfactuals. |
| `diagnostics` | With `diagnostics = TRUE`: `fit_rank` (treated relative pre-treatment RMSPE, placebo median, share of placebos that fit worse, flag), `hull` (number of factors, treated and median donor distance from the convex hull, share of donors further out, flag), `extrapolation` (largest and final implied drift, % of the counterfactual), `persistence` (first-order autocorrelation of each estimated factor, and whether any exceeds 0.95), `assessment` (a sentence summarising the checks) and `recommendation` (arguments to re-run with, or `NULL`). |
| `failed_donors` | Data frame of placebo fits that failed (donor and error message); these donors are left out of the calibration. |
| `settings` | The settings used: `treated_unit`, `treatment_period`, `test_period`, `time_cols`, `outcome`, `conf_level`, `span`, `scale`, `rmspe_floor`, `adjust_covariates`, `adjust_model`, `augment`, `estimator`, `n_donors` and `rank_used` ($m$). |

**Columns of `effects`** (`_pw` = pointwise, `_sim` = simultaneous; interval columns are `NA` before treatment):

| Column(s) | Contents |
|---|---|
| *time columns*, `post` | The period, and whether it is post-treatment. |
| `observed`, `synthetic`, `effect` | The treated unit's outcome, its counterfactual, and the gap (observed $-$ synthetic). |
| `covariate_adjustment` | The amount removed from the treated unit's outcome by `adjust_covariates` (0 without it). |
| `rel_effect`, `std_effect` | The effect as a share of the counterfactual, and on the calibration scale (e.g. in pre-treatment RMSPEs). |
| `factor_pointwise`, `factor_simultaneous` | The calibrated thresholds on the calibration scale (pointwise threshold; simultaneous multiplier $\times$ profile). |
| `eff_low_pw`, `eff_high_pw`, `eff_low_sim`, `eff_high_sim` | Effect intervals. The pointwise columns are filled for every post-treatment period, but each is valid only for a single period chosen in advance. |
| `pred_low_pw`, `pred_high_pw`, `pred_low_sim`, `pred_high_sim` | Prediction intervals for the treated unit's untreated outcome. |
| `threshold_ratio_pw`, `threshold_ratio_sim` | The absolute value of `std_effect` divided by `factor_pointwise` or `factor_simultaneous`: above 1 when the interval or band excludes zero; see [Reading the fit-based benchmark](#reading-the-fit-based-benchmark). |
| `p_pointwise` | Pointwise p-value in each post-treatment period. |
| `cum_effect`, `cum_rel_effect`, `cum_std_effect` | The cumulative effect since the start of treatment, as a share of the cumulative counterfactual, and on the calibration scale. |
| `cum_factor_pointwise`, `cum_factor_simultaneous` | The calibrated thresholds for the cumulative effect. |
| `cum_threshold_ratio_pw`, `cum_threshold_ratio_sim` | The absolute value of `cum_std_effect` divided by `cum_factor_pointwise` or `cum_factor_simultaneous`. |
| `cum_eff_low_pw`, `cum_eff_high_pw`, `cum_eff_low_sim`, `cum_eff_high_sim` | Pointwise intervals and simultaneous bands for the cumulative effect. |
| `cum_pred_low_pw`, `cum_pred_high_pw`, `cum_pred_low_sim`, `cum_pred_high_sim` | Prediction intervals for the treated unit's cumulative untreated outcome. |
| `cum_p_pointwise` | Pointwise p-value for the cumulative effect through each post-treatment period. |

### Reading the fit-based benchmark

On the default RMSPE scale, three columns of `effects` show how large the effect is compared with the synthetic control's pre-treatment fit, benchmarked against the untreated units:

| Column | Meaning |
|---|---|
| `std_effect` | The treated unit's gap in its own pre-treatment RMSPEs: gap $/ s_T$ (with the `rmspe_floor` guard applied to $s_T$). |
| `factor_pointwise` | The threshold: how many of their own pre-treatment RMSPEs the untreated units' synthetic controls miss by in the same period (the $m$-th largest placebo value). `factor_simultaneous` is the threshold for the simultaneous band. |
| `threshold_ratio_pw`, `threshold_ratio_sim` | The absolute value of `std_effect` divided by the pointwise threshold (`factor_pointwise`) or the simultaneous one (`factor_simultaneous`). |

A ratio above 1 means the effect is larger, relative to the pre-treatment fit, than untreated units' gaps get by chance, and the corresponding interval excludes zero. Its size is the margin: 1.3 means the effect would still be significant if untreated drift were 30% larger than the placebo values that set the threshold. This is the calibrated counterpart of the breakdown value $M$ of a sensitivity bound (effect $\pm M \times$ the largest pre-treatment miss, as in the earlier `synth_lasso_boot()`): instead of a multiplier chosen by the analyst, the benchmark is how far untreated synthetic controls actually drift, which grows with time since treatment.

For West Germany at 90%, the 2003 effect is 52.0 pre-treatment RMSPEs and the threshold 40.8 (set by Italy), so `threshold_ratio_pw` = 1.28. The cumulative ratio is 1.49; the simultaneous ratio is above 1 in 1999-2001 (peaking at 1.12 in 2000) and 0.99 in 2003, which is why the band flags 1999-2001 but not 2003.

Read the pointwise ratio only at the pre-specified `test_period`, and `threshold_ratio_sim` for statements about the whole path. With few donors the threshold rests on one or two placebos, so the ratio is coarse. On the relative and raw scales the columns compare percentage or raw gaps instead of multiples of the pre-treatment fit.

### Errors, warnings and limits

- **Too few donors for the level.** If $\lfloor (1 - \text{conf\_level})(N + 1) \rfloor = 0$, the intervals are infinite (`Inf`) and the smallest p-value is above $1 -$ `conf_level`. Lower `conf_level`.
- **Relative scale with non-positive values.** `scale = "relative"` stops if any synthetic value is zero or negative; use `"rmspe"` or `"raw"`.
- **Failed placebo fits** are left out with a warning and listed in `failed_donors`.
- **Constant treated outcome.** A treated unit whose pre-treatment outcome is constant (e.g. all zeros) cannot be fitted.
- **Diagnostics warning.** With `warn = TRUE`, a treated unit flagged as extreme produces a warning that repeats `$diagnostics$assessment`.
- **Run time.** A few seconds for about 50 donors without diagnostics. The diagnostics, covariate adjustment with `adjust_model = "ife"`, `augment = "ife"` and `estimator = "fect"` fit `fect` models and take longer (roughly 10-50 seconds for the examples in this repository).
- **Requirements.** `CVXR`; `fect` for interactive fixed effects (diagnostics, covariate adjustment, augmentation, the fect estimator).

### Example

```r
data("scpi_germany", package = "scpi")
germany <- as.data.frame(scpi_germany)[, c("country", "year", "gdp")]

fit <- dsc_placebo(germany, unit_col = "country", time_cols = "year", outcome_col = "gdp",
                   treated_unit = "West Germany", treatment_period = c(year = 1991),
                   conf_level = 0.90)
fit$test
fit$diagnostics$assessment
```

---

## Plotting functions

Each takes the output of `dsc_placebo()` (or `dsc_symcal_boot()`) and returns a `ggplot` object that can be modified or saved with `ggplot2::ggsave()`. They need `ggplot2`, `dplyr`, `tidyr` and `forcats`. With two or more time columns, the x-axis is labelled at the first sub-period of each value of the first column (e.g. Q1 of each year); with one, at about 12 evenly spaced periods.

### plot_full_path()

```r
plot_full_path(fit, title = NULL, unit_label = fit$settings$treated_unit, y_label = fit$settings$outcome)
```

The whole series: observed and synthetic outcomes, with the pointwise (dark) and simultaneous (light) prediction bands after treatment.

| Argument | Description |
|---|---|
| `fit` | Output of `dsc_placebo()`. |
| `title` | Plot title (`NULL`: automatic). |
| `unit_label` | Name of the treated unit in the legend (default: the treated unit). |
| `y_label` | y-axis label (default: the outcome column name). |

### plot_post_effect()

```r
plot_post_effect(fit, measure = c("effect", "cumulative"), title = NULL, y_label = NULL,
                 fill = NULL, added = NULL, unit_label = fit$settings$treated_unit)
```

The post-treatment effect, per period (`measure = "effect"`) or cumulative (`"cumulative"`): the estimate, the simultaneous band (shaded), and the pointwise interval at the last post-treatment period (error bar).

| Argument | Description |
|---|---|
| `fit` | Output of `dsc_placebo()`. |
| `measure` | `"effect"` or `"cumulative"`. |
| `title`, `y_label` | Plot title and y-axis label (`NULL`: automatic). |
| `fill` | Colour of the band (default: red for the effect, green for the cumulative effect). |
| `added` | Optional known effect to draw as a dashed line, one value per post-treatment period (per-period values; they are cumulated for `measure = "cumulative"`). Typically from `added_effect()`. |
| `unit_label` | Name of the treated unit (default: the treated unit). |

The error bar is always drawn at the last post-treatment period, which is the default `test_period`; with another `test_period`, read the pointwise interval from `$test`.

### plot_prediction()

```r
plot_prediction(fit, title = NULL, unit_label = fit$settings$treated_unit, y_label = fit$settings$outcome)
```

The post-treatment period only: observed and synthetic outcomes, with the simultaneous prediction band (shaded) and the pointwise prediction interval at the last period. Arguments as for `plot_full_path()`.

### added_effect()

```r
added_effect(fit_real, fit_constructed)
```

For constructed examples in which a known effect is added to the treated unit's data: returns the difference between the two fits' observed series in each post-treatment period (the added effect), for `plot_post_effect(..., added = )`. Both fits must cover the same periods.

| Argument | Description |
|---|---|
| `fit_real` | Fit on the original data. |
| `fit_constructed` | Fit on the data with the effect added. |

---

## Earlier functions

Kept for the record of how the method was developed (see [`experimental_methods_tracker.md`](../../experimental_methods_tracker.md)); `dsc_placebo()` replaces them.

| Function | Description |
|---|---|
| `dsc_symcal_boot()`, `plot_dsc_symcal()` | Demeaned synthetic control with a symmetric, donor-calibrated fixed-weight maximum entropy bootstrap; see [the calibrated bootstrap](../bootstrap/fixed_weight_calibrated_bootstrap.md). |
| `synth_lasso_boot()`, `synth_lasso_placebo()`, `meboot_multi()` | The original LASSO synthetic control with bootstrap or exact placebo inference, and a helper that bootstraps several series; documented in [`archive/README_2026-10-15.md`](../../archive/README_2026-10-15.md). They need `glmnet`, `meboot` and, for conformal intervals, `conformalInference`. |
