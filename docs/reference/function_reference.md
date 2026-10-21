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
            augment = c("none", "ife", "ridge"), diagnostics = TRUE, warn = TRUE, ife_seed = 12112025,
            estimator = c("dsc", "fect"), rmspe_floor = 0.25,
            ridge_c = 0.1, ridge_grid = 10^seq(-3, 2, by = 0.5),
            balance_covariates = adjust_covariates, n_cores = 1, verbose = interactive())
```

**What it does.** The counterfactual is an intercept plus non-negative donor weights that sum to one, fitted to the treated unit's pre-treatment outcomes. Every donor is then fitted the same way as if it had been treated (a *placebo*), with the treated unit removed from the pool. Each unit's post-treatment gaps are put on a common scale (by default, divided by the unit's own pre-treatment RMSPE), and the treated unit's gaps are ranked among the placebos'. With $`N`$ donors and $`m = \lfloor (1 - \mathtt{conf\_level})(N + 1) \rfloor`$:

- **pointwise**: the threshold in each period is the $`m`$-th largest placebo value; the interval is the gap $`\pm`$ threshold $`\times`$ the treated unit's scale;
- **simultaneous**: a smoothed profile of the median placebo value over time, and one multiplier (the $`m`$-th largest placebo maximum relative to that profile), so that the band covers the whole post-treatment path at once.

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
| `conf_level` | `0.95` | Confidence level for the intervals and bands. With $`N`$ donors the smallest possible p-value is $`1/(N + 1)`$: a 95% interval needs at least 19 donors and a 90% interval at least 9; with fewer, the intervals are infinite. |
| `span` | `0.75` | Loess span for smoothing the simultaneous band's time profile (on the log scale). Larger values give a smoother profile. |
| `scale` | `"rmspe"` | What each unit's gaps are divided by before ranking. `"rmspe"`: the unit's pre-treatment RMSPE (gaps measured in units of the unit's own pre-treatment fit; the most robust when the treated unit is unlike the donors). `"relative"`: the unit's synthetic value, i.e. percentage gaps (the most powerful when the treated unit resembles the donors; needs positive outcomes). `"raw"`: nothing (gaps in the outcome's units). |
| `rmspe_floor` | `0.25` | Guard for `scale = "rmspe"`: each unit's pre-treatment RMSPE is floored at `rmspe_floor` $`\times`$ the median pre-treatment RMSPE over the treated unit and all placebos, so that a unit reproduced (almost) exactly before treatment cannot produce near-infinite intervals. Matters mainly with short pre-treatment periods and many donors. `0` switches it off. Ignored for the other scales. |

**Covariate adjustment** (for measured time-varying confounders)

| Argument | Default | Description |
|---|---|---|
| `adjust_covariates` | `NULL` | Name(s) of time-varying covariate columns to adjust for. Their effect $`\gamma`$ is estimated on untreated cells only (all donor periods and the treated unit's pre-treatment periods); every unit's outcome is then replaced by outcome $`- \hat\gamma \times`$ covariate, and the method runs on the adjusted outcomes. Each placebo re-estimates $`\gamma`$ without its own post-treatment period. Assumes the treatment does not affect the covariates and that their effect is common and linear. See the [Oregon case study](../oregon_fentanyl/oregon_fentanyl_case_study.md). |
| `adjust_model` | `"ife"` | Model used to estimate $`\gamma`$: `"ife"`, interactive fixed effects (unit and time effects plus latent factors, via `fect`), or `"twfe"`, two-way fixed effects. Ignored with `augment = "ife"` or `estimator = "fect"`, where the covariates enter the interactive fixed-effects outcome model directly. Used with `augment = "ridge"` (the outcomes are adjusted first). |
| `ife_r` | `0:5` | Candidate numbers of latent factors for every interactive fixed-effects model (covariate adjustment, `augment = "ife"`, `estimator = "fect"` and the diagnostics). The number is chosen by cross-validation for the treated unit and reused for every placebo. A single value skips cross-validation. |

**Estimator**

| Argument | Default | Description |
|---|---|---|
| `estimator` | `"dsc"` | `"dsc"`: the demeaned synthetic control. `"fect"`: the counterfactual is `fect`'s interactive fixed-effects prediction, fitted on untreated cells, with each donor refitted as a pseudo-treated unit and the same placebo calibration (use with `scale = "rmspe"`). Recommended by the diagnostics only for a treated unit beyond the donors' range with stationary factors. |
| `augment` | `"none"` | Corrects the counterfactual for a treated unit the weights can't match (at or beyond the edge of the donors' range).<br>• `"ridge"`: **ridge-augmented demeaned synthetic control**, in the spirit of Ben-Michael, Feller and Rothstein's augmented SCM. In each post-treatment period the counterfactual is the demeaned synthetic control plus a correction: a ridge regression, across the donors, of their outcome in that period on their demeaned pre-treatment paths, applied to the part of the treated unit's pre-treatment path the weights leave unmatched. The correction is about zero when the weights already match. Pre-treatment gaps, and so the RMSPE scale, are the demeaned synthetic control's; each placebo refits its own correction. Validated in simulation (see [Ridge augmentation](#ridge-augmentation)): it restored coverage for a treated unit at or beyond the edge of the donors' range with trending factors, and lowered the error almost everywhere. The diagnostics recommend it in that situation.<br>• `"ife"`: augmented demeaned synthetic control with an interactive fixed-effects outcome model fitted on untreated cells; the synthetic control is fitted to its residuals, and the counterfactual is the model prediction plus the synthetic control of the residuals. Reduces bias for a treated unit outside the donors' range, but raised variance and made intervals overconfident in simulation.<br>**Combinations:**<br>• `"ridge"` with `adjust_covariates`: the outcomes are adjusted first (with `adjust_model`), and the ridge model is fitted to the adjusted outcomes. With `adjust_model = "ife"` this was the best combination tested for a confounded treated unit.<br>• `"ife"` with `adjust_covariates`: the covariates enter the same outcome model, and γ is estimated jointly with the factors (`adjust_model` is then ignored). This combination has not been validated in simulation.<br>• Either with `estimator = "fect"`: not allowed; the function stops with an error. |
| `ridge_c` | `0.1` | `augment = "ridge"` only. The ridge penalty, as a multiple of the donors' mean squared (centred) pre-treatment path norm. The default 0.1 was chosen by simulation (see "Choosing the penalty" below). `"cv"` (or `NULL`) chooses it from `ridge_grid` by cross-validation inside the pre-treatment period instead: each donor, fitted on the other donors, predicts its last min(post-treatment periods, T₀/4) pre-treatment periods from the earlier ones. Cross-validation needs at least 8 pre-treatment periods. Larger values shrink the correction towards zero (the plain demeaned synthetic control). |
| `ridge_grid` | `10^seq(-3, 2, by = 0.5)` | `augment = "ridge"` only. The candidate values of `ridge_c` for the cross-validation (`ridge_c = "cv"`). |
| `n_cores` | `1` | Number of cores for the per-placebo `fect` models (covariate adjustment with `adjust_model = "ife"`, `augment = "ife"`, `estimator = "fect"`). Uses base R's `parallel`; the workers are started and stopped inside the call, which costs a few seconds, so it pays off for long runs only. Results are identical to `n_cores = 1`. Ignored when there are no `fect` models. |
| `verbose` | `interactive()` | Report the slow stages (cross-validation, placebo models, diagnostics) and show a progress bar when the placebo models are fitted sequentially. |
| `balance_covariates` | `adjust_covariates` | Numeric columns whose balance is checked; a check, not an adjustment. Each is summarized by its pre-treatment mean per unit. The weights fitted on the outcome path are applied to the donors' values, and the treated unit's gap from its synthetic control (standardized by the donors' SD) is ranked against each placebo's gap, computed with that placebo's own weights. Time-invariant characteristics can be named too. `NULL` (with no `adjust_covariates`) skips it. See "Covariate balance" below. |

**Diagnostics and reproducibility**

| Argument | Default | Description |
|---|---|---|
| `diagnostics` | `TRUE` | Check whether the treated unit lies beyond what the donors can reproduce (pre-treatment fit rank, distance from the donors' convex hull in estimated factor loadings, implied extrapolation drift, factor persistence), and report how persistent the untreated gaps are (residual autocorrelation). Loading-based checks need `fect` and take most of the run time; `FALSE` skips them. See [When placebo inference fails](../diagnostics/diagnostics_toy_cases.md). |
| `warn` | `TRUE` | Issue a warning when the diagnostics flag the treated unit as extreme. |
| `ife_seed` | `12112025` | Seed for `fect`'s cross-validation, applied locally (the session's random-number state is restored afterwards), so results with interactive fixed effects are reproducible. |

### Value

A list:

| Element | Contents |
|---|---|
| `effects` | Data frame, one row per period (all periods; post-treatment rows have `post = TRUE`). Columns below. |
| `test` | Data frame with five rows: the pointwise test of the effect at `test_period`, the pointwise test of the cumulative effect through `test_period`, the simultaneous tests of the effect and of the cumulative effect over the whole post-treatment period, and the classic post/pre RMSPE-ratio test (Abadie et al.). Columns: `test`, `measure`, `period`, `estimate`, `lower`, `upper` (the interval; pointwise tests only), `p_value`, `significant` (the interval excludes zero, or for the simultaneous and ratio tests `p_value` $`\le 1 -`$ `conf_level`), `periods_flagged` (simultaneous tests: the periods where the band excludes zero). |
| `weights` | Named vector of non-zero donor weights, largest first (`NULL` with `estimator = "fect"`). |
| `intercept` | The intercept of the synthetic control (`NULL` with `estimator = "fect"`). |
| `pre_rmspe` | The treated unit's pre-treatment RMSPE (unfloored). |
| `rmspe_guard` | `floor`: the floor applied to the pre-treatment RMSPEs (0 if none); `units_floored`: the units it applied to. |
| `ridge` | With `augment = "ridge"`: `c` (the penalty used), `chosen_by` (`"default (0.1)"`, `"ridge_c"` or `"cross-validation"`), `holdout_periods` (the cross-validation's held-out pre-treatment periods), `cv` (with cross-validation: data frame of `c` and cross-validated mean squared error), `at_grid_edge` (whether the chosen value is at either end of `ridge_grid`), `note` (a caution when it is), and `correction` (the correction subtracted from the treated unit's gap in each period; 0 before treatment). `NULL` otherwise. |
| `balance` | With balance covariates and `estimator = "dsc"`: one row per covariate with `treated` and `synthetic` (pre-treatment means), `std_diff` (their difference in donor SDs), `placebo_median_abs_std_diff`, `share_placebos_worse` (share of placebos with a larger absolute gap) and `flag` (`share_placebos_worse` below 0.10). A flag is also noted in `$diagnostics$assessment`. `NULL` otherwise. |
| `adjustment` | With `adjust_covariates`: `model`, `gamma` (the estimates), `n_factors`, `gamma_twfe` and `se_twfe` (two-way fixed-effects estimate with a unit-clustered standard error, for comparison), `n_obs`, and `unadjusted` (the effect and cumulative effect at `test_period` without the adjustment). `NULL` otherwise. |
| `placebo` | `gaps` and `synthetic`: data frames (periods $`\times`$ donors) of each placebo's gaps and counterfactuals. |
| `diagnostics` | With `diagnostics = TRUE`: `fit_rank` (treated relative pre-treatment RMSPE, placebo median, share of placebos that fit worse, flag), `hull` (number of factors, treated and median donor distance from the convex hull, share of donors further out, flag), `extrapolation` (largest and final implied drift, % of the counterfactual), `persistence` (first-order autocorrelation of each estimated factor, and whether any exceeds 0.95), `residual_autocorrelation` (persistence of the untreated gaps: `rho`, `rho_interval_90`, `level`, `half_life`, `backdate_periods`, `effective_pre_periods`, `rho_in_sample`, `seasonal_lag`, `rho_seasonal`, `rho_seasonal_excess`, `n_post_periods`, `n_placebos`, `guide`; see [Reading the residual autocorrelation](#reading-the-residual-autocorrelation)), `assessment` (a sentence summarising the checks) and `recommendation` (arguments to re-run with, or `NULL`: `augment = "ridge"` for a treated unit flagged as extreme with trending factors; `estimator = "fect"`, `scale = "rmspe"` for one flagged with stationary factors). |
| `failed_donors` | Data frame of placebo fits that failed (donor and error message); these donors are left out of the calibration. |
| `settings` | The settings used: `treated_unit`, `treatment_period`, `test_period`, `time_cols`, `outcome`, `conf_level`, `span`, `scale`, `rmspe_floor`, `adjust_covariates`, `adjust_model`, `augment`, `estimator`, `ridge_c` (the penalty used with `augment = "ridge"`, else `NA`), `n_donors` and `rank_used` ($`m`$). |

**Columns of `effects`** (`_pw` = pointwise, `_sim` = simultaneous; interval columns are `NA` before treatment):

| Column(s) | Contents |
|---|---|
| *time columns*, `post` | The period, and whether it is post-treatment. |
| `observed`, `synthetic`, `effect` | The treated unit's outcome, its counterfactual (including the ridge correction with `augment = "ridge"`), and the gap (observed $`-`$ synthetic). |
| `covariate_adjustment` | The amount removed from the treated unit's outcome by `adjust_covariates` (0 without it). |
| `rel_effect`, `std_effect` | The effect as a share of the counterfactual, and on the calibration scale (e.g. in pre-treatment RMSPEs). |
| `factor_pointwise`, `factor_simultaneous` | The calibrated thresholds on the calibration scale (pointwise threshold; simultaneous multiplier $`\times`$ profile). |
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
| `std_effect` | The treated unit's gap in its own pre-treatment RMSPEs: gap $`/ s_T`$ (with the `rmspe_floor` guard applied to $`s_T`$). |
| `factor_pointwise` | The threshold: how many of their own pre-treatment RMSPEs the untreated units' synthetic controls miss by in the same period (the $`m`$-th largest placebo value). `factor_simultaneous` is the threshold for the simultaneous band. |
| `threshold_ratio_pw`, `threshold_ratio_sim` | The absolute value of `std_effect` divided by the pointwise threshold (`factor_pointwise`) or the simultaneous one (`factor_simultaneous`). |

A ratio above 1 means the effect is larger, relative to the pre-treatment fit, than untreated units' gaps get by chance, and the corresponding interval excludes zero. Its size is the margin: 1.3 means the effect would still be significant if untreated drift were 30% larger than the placebo values that set the threshold. This is the calibrated counterpart of the breakdown value $`M`$ of a sensitivity bound (effect $`\pm M \times`$ the largest pre-treatment miss, as in the earlier `synth_lasso_boot()`): instead of a multiplier chosen by the analyst, the benchmark is how far untreated synthetic controls actually drift, which grows with time since treatment.

For West Germany at 90%, the 2003 effect is 52.0 pre-treatment RMSPEs and the threshold 40.8 (set by Italy), so `threshold_ratio_pw` = 1.28. The cumulative ratio is 1.49; the simultaneous ratio is above 1 in 1999-2001 (peaking at 1.12 in 2000) and 0.99 in 2003, which is why the band flags 1999-2001 but not 2003.

Read the pointwise ratio only at the pre-specified `test_period`, and `threshold_ratio_sim` for statements about the whole path. With few donors the threshold rests on one or two placebos, so the ratio is coarse. On the relative and raw scales the columns compare percentage or raw gaps instead of multiples of the pre-treatment fit.

### Reading the residual autocorrelation

`$diagnostics$residual_autocorrelation` measures how persistent the untreated gaps are: how much of this period's gap is still there next period. It is information, not a flag. It doesn't change any result or the assessment; it tells you how to read them.

**How it is estimated.**
- **Data:** the placebos' post-treatment gaps. The donors are untreated, and those periods lie outside their fit, so these gaps are genuine out-of-sample errors: exactly the kind of error the treated unit's counterfactual makes.
- **Estimator:** each placebo's gaps are divided by its pre-treatment RMSPE, then pooled as a first-order autoregression through the origin, giving `rho`.
- **Uncertainty:** `rho_interval_90` is a 90% bootstrap interval over placebos.
- **Requirement:** at least 3 post-treatment periods.

**Why not the pre-treatment residuals?** They understate persistence badly, because the weights absorb smooth noise. In simulation (two factors, AR(1) errors, 29 donors), the median in-sample residual autocorrelation compared with the truth was:

| True persistence | 10 pre-periods | 20 pre-periods | 40 pre-periods |
|---|---|---|---|
| 0.7 | −0.12 | 0.19 | 0.38 |
| 0.9 | −0.07 | 0.23 | 0.50 |

The out-of-sample estimate recovered 0.70–0.71 for 0.7, and 0.83–0.86 for 0.9, at all three lengths. The in-sample value is still reported as `rho_in_sample`: a large gap between the two means the weights absorb much of the noise, so the pre-treatment fit overstates how closely the synthetic control will track the treated unit. The guide says so when the gap exceeds 0.3.

**What the level means:**

| `level` | `rho` | What it means |
|---|---|---|
| low | < 0.3 | Gaps behave close to independent noise. A close pre-treatment fit is informative; a shock just before treatment fades within a period; backdating by 1–2 periods is enough. |
| moderate | 0.3–0.7 | Shocks fade over a few periods. The pre-treatment fit is worth fewer independent periods than it seems; a pre-treatment shock carries into the first post-treatment periods, so backdate by `backdate_periods` to check for anticipation or such a shock. |
| high | 0.7–0.9 | Smooth, sustained drifts in the gap happen by chance, so a convincing-looking gap plot is weaker evidence than it appears. Backdating needs more periods; check that enough pre-period remains to fit on. |
| very high | ≥ 0.9 | The gaps behave like a random walk, which is common for outcomes in levels: GDP, sales, trending rates. A close pre-treatment fit is weak evidence, gaps fan out with the horizon, and backdating can't settle whether a pre-treatment shock carries over (`backdate_periods` is `NA`). A value at or above 1 means the placebos' gaps drift away rather than revert. |

**The supporting numbers:**
- **`half_life`:** periods for a shock to halve, log(0.5)/log(ρ); `Inf` when ρ ≥ 1.
- **`backdate_periods`:** periods for a shock to fall to 25%, a guide to how far to backdate `treatment_period`. Choose the window from substance first (announcement and debate dates), then check it against this.
- **`effective_pre_periods`:** roughly how many independent periods the pre-period is worth, T₀(1 − ρ)/(1 + ρ). At ρ = 0.9, 20 pre-periods are worth about one.

**What it doesn't change.** If the treated unit is like the donors, the placebo calibration already accounts for persistence. The placebos drift as much as the treated unit would, and the simultaneous band's drift profile grows over the post-period accordingly. So a high value is not a reason to distrust the p-values. It is a reason to:
- give less weight to a close pre-treatment fit and to the shape of the gap plot;
- think harder about whether a pre-treatment shock or anticipation is carrying over;
- prefer later post-treatment periods if that is a concern.

**Seasonal persistence.** When the finest time column has 2–12 values (quarters or months), `rho_seasonal` is the pooled autocorrelation at that lag (4 or 12). It needs at least lag + 3 post-treatment periods. Because a persistent series is also correlated at the seasonal lag, the guide uses `rho_seasonal_excess` = `rho_seasonal` − ρ^lag: the seasonal persistence beyond what lag 1 implies. When this exceeds 0.2, a shock in one season tends to recur a year later. Backdate in whole years and consider deseasonalizing each unit before fitting.

**Examples:**

| Panel | `rho` | `level` | `rho_in_sample` | Seasonal excess |
|---|---|---|---|---|
| West Germany (annual GDP, 1991) | 1.05 | very high | 0.56 | — |
| California (annual cigarette sales, 1989) | 0.96 | very high | −0.01 | — |
| Kansas (quarterly GDP, 2012 Q2) | 0.99 | very high | 0.74 | 0.03 (lag 4) |
| Simulated (stationary factors, AR(1) errors 0.3) | 0.52 | moderate | 0.11 | — |
| Simulated (AR(1) errors 0.7) | 0.78 | high | 0.21 | — |

All three real panels have random-walk-like untreated gaps, which matches how much their simultaneous bands widen over the post-period. The simulated panels come out higher than their error persistence: the gap also contains the factor-exposure mismatch, which has persistence of its own (AR 0.5 here).

**Limits:**
- **It describes the donors' gaps, not the treated unit's.** A treated unit more persistent than its donors isn't detected.
- **It needs post-treatment periods,** so it isn't available for a design check before treatment.
- **It is wide with few placebos or few post-treatment periods;** check `rho_interval_90`.
- **Disturbed donors distort it.** Spillovers or shocks to donors in the post-period affect it.
- **Steady divergence also reads as persistence.** A systematic drift of the placebos' gaps (a deterministic divergence) produces a value near 1, just as a random walk does.

### Ridge augmentation

`augment = "ridge"` targets the case where placebo inference otherwise fails: a treated unit at or beyond the edge of the donors' range. There, non-negative weights summing to one can't reproduce its pre-treatment path, and the leftover imbalance becomes bias. The ridge correction learns from the donors how pre-treatment differences carry into each post-treatment period, and applies that to the treated unit's leftover imbalance.

**How it differs from `augsynth`.** It is the same estimator as `augsynth(progfunc = "Ridge", fixedeff = TRUE)`: given the same penalty, the two agree to within 0.00003 on Kansas, California and West Germany. They differ in how the penalty is chosen and in the inference.
- **Penalty:** `augsynth` fits the treated unit's own pre-period, leaving out one period at a time; `dsc_placebo()` predicts donors' late pre-treatment periods from their earlier ones.
- **Inference:** `augsynth` uses conformal inference, which over-rejected in simulation; `dsc_placebo()` uses the same placebo calibration as without augmentation, each placebo refitting its own correction.

**Evidence** (simulations on the factor-model generator; coverage of 90% intervals):

| Situation | Without augmentation | `augment = "ridge"` |
|---|---|---|
| Treated unit at the edge of the donors' range, trending factors (19-39 donors, 40-80 pre-periods) | 62% (bias +16.5%) | 91% (bias +3.3%) |
| Treated unit outside the donors' range, trending factors (19-39 donors, 40-80 pre-periods) | 49% (bias +38%) | 91% (bias +9%) |
| Treated unit outside the donors' range, trending factors, 9-19 donors, 20 pre-periods | 76% | 90% |
| Treated unit like the donors (no covariate; 9-39 donors, 10-80 pre-periods) | 89-93% | 88-94% |
| Confounded covariate, outside the range, trending factors, with `adjust_covariates` (IFE) | 73% | 89% (cumulative 93%) |
| Fake treated units in Kansas, California and West Germany | 88-92% | 88-92% |

Point estimates had lower RMSE than without augmentation almost everywhere, by half or more for extreme treated units. Power was about the same when the treated unit is like the donors. On the three real cases the estimates barely change (Kansas −2,006 vs −1,891 in the final quarter; California −17.3 vs −17.4; West Germany −3.51 vs −3.49).

**Choosing the penalty.** The default is a fixed `ridge_c = 0.1`; `ridge_c = "cv"` cross-validates it. The evidence above used cross-validation. `Experimental_inference/LTO/simulations/ridge_default_validation.R` compared the two, along with rolling-window cross-validation, for `dsc_placebo()` and the LTO test. It used 9-14 donors, 10-40 pre-periods and the confounded covariate, plus subsamples of four real panels.
- **Extreme treated units with trending factors:** both remove the excess false positives (24-28% → 6-8% at 10%).
- **Real panels:** the fixed value had lower error (RMSE 15.2% against 18.8% for flagged units) and more power.
- **10 pre-periods:** the fixed value had lower error (19.1% against 22.8%).
- **40 pre-periods:** the fixed value had slightly more false positives (8.8% against 7.1%).
- **Stability:** a fixed value does not change when one donor is dropped; the cross-validated one changed about half the time with 9 donors.

Report cross-validation as a sensitivity check when the correction matters. On the Basque Country, for example, cross-validation corrected more (final-year effect −1.00 against −1.39 with 0.1), and these data can't tell which is right.

**Cautions:**
- **Very short pre-periods.** With about 10 pre-treatment periods, the cross-validated penalty (`ridge_c = "cv"`) was less stable: at an end of `ridge_grid` in about 5-10% of replicates. For a noisy treated unit with trending factors the error was larger than without augmentation (RMSE 20% vs 16%), although coverage held. `$ridge$note` flags a penalty at the grid's edge; compare with `augment = "none"`.
- **Bias can remain.** For a treated unit well outside the donors' range with trending factors, bias fell from about +38% to +9% but did not vanish.
- **With stationary factors it is conservative** for extreme treated units (99-100% coverage at 90%).
- **Covariates.** Adjust with `adjust_model = "ife"` when the covariate's timing may be related to the factors: the two-way fixed-effects estimate of γ was badly biased in that design.

Simulations: `simulations/ridge_edge_validation.R`, `simulations/ridge_small_sample_validation.R`, `simulations/ridge_covariate_validation.R`; agreement of `dsc_placebo()` with the validated prototype: `simulations/ridge_function_check.R`.

### Covariate balance

`balance_covariates` checks whether the synthetic control matched units that are structurally similar to the treated unit, not just its outcome path. The weights come only from the fit on pre-treatment outcomes; the covariates play no part in the fit, so they give an independent check.

1. **Summarize** each covariate by its pre-treatment mean per unit (missing values are skipped).
2. **Compare** the treated unit's value with its synthetic control's: the outcome-fitted weights applied to the donors' values. The gap is standardized by the donors' standard deviation (`std_diff`).
3. **Rank** that gap against the placebos' gaps, each computed with that placebo's own weights. The covariate is flagged when fewer than 10% of placebos have a larger gap, as in the fit-rank check.

It reports only: no estimate, interval or p-value changes. By default it checks the `adjust_covariates`; any other numeric column can be named, including characteristics that don't change over time, and these are not adjusted for. A flagged imbalance matters only if the covariate's effect on the outcome changes over time; a constant effect is removed by the demeaning. Pickett, Hill and Cowan (2026) argue that a close pre-treatment fit alone can mislead with short pre-periods, because it can fit noise; this is one way to check it. The weights used are the synthetic control's, before any ridge correction.

### Errors, warnings and limits

- **Too few donors for the level.** If $`\lfloor (1 - \mathtt{conf\_level})(N + 1) \rfloor = 0`$, the intervals are infinite (`Inf`) and the smallest p-value is above $`1 -`$ `conf_level`. Lower `conf_level`.
- **Relative scale with non-positive values.** `scale = "relative"` stops if any synthetic value is zero or negative; use `"rmspe"` or `"raw"`.
- **Failed placebo fits** are left out with a warning and listed in `failed_donors`.
- **Constant treated outcome.** A treated unit whose pre-treatment outcome is constant (e.g. all zeros) cannot be fitted.
- **Diagnostics warning.** With `warn = TRUE`, a treated unit flagged as extreme produces a warning that repeats `$diagnostics$assessment`.
- **Ridge augmentation with short pre-periods.** `augment = "ridge"` with `ridge_c = "cv"` stops if there are fewer than 8 pre-treatment periods; the default fixed penalty has no such limit.
- **Run time.** A few seconds for about 50 donors without covariate adjustment (`augment = "ridge"` with `ridge_c = "cv"` adds one cross-validation fit per donor; the diagnostics add one or two `fect` fits). Covariate adjustment with `adjust_model = "ife"`, `augment = "ife"` and `estimator = "fect"` fit one `fect` model per placebo. On Kansas (48 donors, 105 quarters, three covariates) that took 10.9 minutes: about 11 seconds per placebo model and 35 seconds for the cross-validation, 97% of the run. With `n_cores = 4` it took 3.6 minutes, with identical results. The time grows with the number of donors, periods and covariates.
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
