# Simulations

*Written by Claude, an AI model made by Anthropic, at the direction of Jeffrey Milliman.*

The simulations behind the main design choices in `dsc_placebo()`, with their summary results. Each script writes its results to the `results_*` folder of the same name. The per-task files a run produces are not kept here (see `.gitignore`), so the folders hold only the summaries: `results.rds`, and for most scripts a `summary.csv` and an interactive DT table (`.html`).

Run scripts from the **project root**, e.g. `source("simulations/factor_model_simulation.R")`. They use `synthetic_control_functions/synthetic_control_functions_claude.R`, run in parallel with `future` and `furrr::future_pmap()` (10 workers by default), and resume from saved tasks if interrupted. The ridge scripts are the exception: they use the faster prototype in `ridge_validation_helpers.R` (14 workers by default, set `n_workers` to change), don't resume, and their run times are for 14 workers; `ridge_function_check.R` confirms the prototype matches `dsc_placebo()`.

| Script | Question | Data | Run time (10 workers) | Results |
|---|---|---|---|---|
| `california_estimator_accuracy_test.R` | Which point estimator? | California panel, each state as fake treated | ~2 min | `results_california_accuracy/` |
| `pre_period_length_simulation.R` | How short can the pre-treatment period be? | Kansas and California, pre-periods cut down | ~27 min | `results_pre_period_length/` |
| `covariate_adjustment_simulation.R` | Can a confounder that arrives with the treatment be adjusted for? | California panel with a simulated shock | ~11 min | `results_covariate_adjustment/` |
| `simultaneous_band_ablation.R` | Does the bootstrap add anything once placebos calibrate the intervals? Which scale? | Kansas and California panels | ~5 min | `results_simultaneous_band_ablation/` |
| `factor_model_simulation.R` | Where do placebo methods fail? (treated units unlike the donors) | Simulated factor model | ~12 hours | `results_factor_model/` |
| `augmented_comparison_simulation.R` | Augmented demeaned SC vs fect vs augsynth | Simulated factor model | ~4.5 hours | `results_augmented_comparison/` |
| `diagnostics_validation.R` | Do the diagnostics flag the failing cases? | The 1,800 panels of the augmented comparison | ~1 hour | `results_diagnostics_validation/` |
| `oregon_fentanyl_validation.R` | Does the covariate adjustment remove a confounding shock that arrives with the treatment? | Oregon overdose panel (Zoorob et al. 2024), each state as fake treated with an Oregon-like fentanyl shock | ~4 min | `results_oregon_fentanyl/` |
| `rmspe_guard_validation.R` | Does the guard on the RMSPE scale (`rmspe_floor`) keep the intervals valid and usable with short pre-periods? | Kansas, California and Texas prison panels with short pre-periods; simulated factor model | ~9 min | `results_rmspe_guard/` |
| `firpo_possebom_comparison.R` | How do `dsc_placebo()`'s intervals relate to Firpo and Possebom's (2018) confidence sets? | West Germany (`scpi::scpi_germany`) | ~2 min | `results_firpo_possebom/` |
| `ridge_edge_validation.R` | Does ridge augmentation (or in-time placebo calibration) restore coverage for a treated unit at or beyond the edge of the donors' range? | Simulated factor model (19-39 donors, 40-80 pre-periods) | ~2 min | `results_ridge_edge/` |
| `ridge_small_sample_validation.R` | Does ridge augmentation stay valid with small pools and short pre-periods, on real panels, and does it match `augsynth`? | Simulated factor model (9-19 donors, 10-20 pre-periods); Kansas, California and West Germany | ~4 min | `results_ridge_small_sample/` |
| `ridge_covariate_validation.R` | Ridge augmentation after covariate adjustment, with a confounded covariate | Simulated factor model with the confounded covariate | ~36 min | `results_ridge_covariate/` |
| `ridge_function_check.R` | Does `dsc_placebo(augment = "ridge")` reproduce the prototype the three simulations above validated? | Kansas, California, West Germany; 48 simulated panels | ~3 min | `results_ridge_function_check/` |
| `ridge_validation_helpers.R` | Helper functions for the four ridge scripts (the prototype of `augment = "ridge"`) | | | |
| `diagnostics_toy_cases.R` | Helper functions and toy panels for the diagnostics write-up | Simulated | seconds | figures in `docs/diagnostics/` |
| `demeaned_sc_symmetric_calboot_test.R` | Dependency of `pre_period_length_simulation.R` (test version of `dsc_symcal_boot()`) | | | |

## Which simulation supports which choice

The method was built in stages, and each stage changed one part of it. The table maps each choice in the current `dsc_placebo()` to the simulation that motivated it. The sections below give the details. Row numbers refer to `../experimental_methods_tracker.md`.

| Choice in `dsc_placebo()` | Alternative it replaced or beat | Evidence | Tracker rows |
|---|---|---|---|
| **Demeaned synthetic control** (intercept + non-negative weights summing to one) | LASSO synthetic control with weights in [0, 1]; relaxed, ridge and elastic-net variants | `california_estimator_accuracy_test.R`: essentially unbiased, 17-20% lower RMSE than the LASSO | 17, 18 |
| **No minimum pre-period beyond what the fit needs** (the method is used with short pre-periods) | LASSO, which lost almost all power with short pre-periods | `pre_period_length_simulation.R`: valid and unbiased down to 8 quarters / 5 years | 22 |
| **Placebo calibration with no bootstrap** | The maximum entropy bootstrap (`dsc_symcal_boot()`) | `simultaneous_band_ablation.R`: the same construction without the bootstrap was valid and at least as powerful almost everywhere | 25, 26 |
| **Raw gaps not the default** (kept as an option) | Raw gaps, the conventional placebo-plot scale | `simultaneous_band_ablation.R`: valid but by far the least powerful | 27 |
| **RMSPE scale as the default** (relative kept as an option) | Relative scale, the original default | `factor_model_simulation.R`: RMSPE stayed largely valid for treated units unlike the donors; relative and raw failed | 29, 33 |
| **Guard on the RMSPE scale** (`rmspe_floor = 0.25`: each unit's pre-treatment RMSPE floored at 0.25 x the median) | No guard | `rmspe_guard_validation.R`: without it, intervals are effectively infinite with short pre-periods; with it, coverage 92-98% and usable intervals | 35 |
| **Pooled covariate adjustment estimated on untreated cells** (`adjust_covariates`) | No adjustment; unit-specific regressors (`controls`) | `covariate_adjustment_simulation.R`: removed the confounding bias; unit-specific regressors were badly biased | 24 |
| **Covariate adjustment for a confounder that arrives with the treatment, in a real panel** | Unadjusted synthetic control | `oregon_fentanyl_validation.R`: with an Oregon-like fentanyl shock, unadjusted false positives 33-39%; adjusted 0-4%, coverage 96-100% | 34 |
| **Interactive fixed effects as the default adjustment model** (`adjust_model = "ife"`) | Two-way fixed effects | `factor_model_simulation.R`: IFE recovered the covariate effect; TWFE was badly biased with trending factors | 28, 29 |
| **No augmentation by default** (`augment = "none"`) | Augmented demeaned SC (`augment = "ife"`) | `augmented_comparison_simulation.R`: augmentation reduced bias outside the donors' range but raised variance and made intervals overconfident | 30 |
| **fect offered only as an option, with placebo calibration** (`estimator = "fect"`) | fect's own parametric inference; augsynth's conformal inference | `augmented_comparison_simulation.R`: placebo-calibrated fect valid with stationary factors; fect's parametric and augsynth's conformal inference over-rejected | 30, 32 |
| **Diagnostics and a warning that needs an extreme unit AND trending factors** (persistence threshold 0.95) | No diagnostics; an extreme-unit flag on its own | `diagnostics_validation.R`: the combined warning separated failing runs (74% coverage when warned vs 96% when not) | 31 |
| **Recommendation to switch to fect only for stationary factors** | Always DSC; always fect | `diagnostics_toy_cases.R` (checked) and `augmented_comparison_simulation.R` (coverage) | 32 |
| **Ridge augmentation as an option** (`augment = "ridge"`), recommended by the diagnostics for an extreme treated unit with trending factors | No remedy (results flagged as unreliable); `augment = "ife"`; in-time placebo calibration | `ridge_edge_validation.R`: coverage 91% at 90% for edge and outside units with random-walk factors, against 49-62% without; in-time calibration fell short. `ridge_small_sample_validation.R`: valid with 9-19 donors and 10-20 pre-periods, and on fake treated units in three real panels | 38 |
| **Ridge augmentation combinable with covariate adjustment** (adjust first, then ridge) | Either one alone | `ridge_covariate_validation.R`: with a confounded covariate, IFE adjustment + ridge covered 89-99% in every scenario; adjustment alone under-covered extreme units and ridge alone didn't remove the confounding | 38 |

Some earlier choices come from simulations that are **not kept in this folder** (they remain in the git-ignored `Placebo simulations/` folder and are summarized in `../experimental_methods_tracker.md`):

- the **smoothed-median time profile** for the simultaneous band, and its insensitivity to the loess span: `california_band_coverage_simulation.R` (tracker row 20);
- the **exact rank rule** rather than a parametric multiplier for small donor pools: `california_donor_subsample_simulation.R` (row 21);
- the validity of the **absolute (two-sided) placebo test**: the `Placebo_exact_*.R` simulations (row 7).

## The simulations in order

### 1. Estimator accuracy (`california_estimator_accuracy_test.R`)

**Question.** Which point estimator should the method use? The earlier LASSO synthetic control was biased towards zero.

**Design.** Every state in the California tobacco panel (38 states, 19 pre-treatment years) is treated in turn as a fake treated unit with no effect. Each estimator's error in the cumulative and final-year effect is recorded. The estimators compared are the LASSO (weights in [0, 1] with an intercept, under several ways of choosing the penalty), relaxed and unconstrained LASSOs, ridge, elastic net, a no-intercept LASSO, and the demeaned synthetic control.

**Results.**

- The **demeaned synthetic control** was essentially unbiased (cumulative +0.1%, final year -0.9%), with 17-20% lower RMSE than the LASSO.
- The LASSO's bias came from **weights summing to less than one** (about 0.6-0.8): the penalty shrinks the weights while the intercept absorbs the level, which damps the donors' trend.
- Relaxing the weight bounds, ridge and elastic net **did not remove the bias**, because the weights still summed to less than one. The no-intercept LASSO had little bias but tracked the pre-period poorly (median pre-RMSPE 8.4 vs 1.8).

**Supports.** The demeaned synthetic control as the estimator, with no penalty to tune (under the sum-to-one constraint the LASSO penalty is constant). The cost is a slightly worse pre-treatment fit (median RMSPE 1.94 vs 1.78).

### 2. Pre-period length (`pre_period_length_simulation.R`)

**Question.** Does the method stay valid and useful when the pre-treatment period is short?

**Design.** The pre-treatment period is cut to 80, 40, 20, 12 and 8 quarters (Kansas, 49 states) and 19, 12, 8 and 5 years (California, 38 states), with the post period fixed. Each state takes a turn as the treated unit, with effects injected. Three methods are compared at 95% and 90%: the demeaned SC with the absolute placebo test on relative gaps; the demeaned SC with the symmetric calibrated bootstrap (pointwise and smoothed-median simultaneous bands); and the LASSO with the absolute placebo test.

**Results.**

- **Validity held at every length.** Final-period coverage was 94-100% for all methods, and the simultaneous bands kept whole-path coverage at 95-97%.
- **The demeaned SC stayed unbiased at every length**, and its cumulative RMSE barely grew (Kansas 4.4-4.7%; California 9.9-11.7%).
- **The demeaned SC lost far less power than the LASSO.** At a -10% cumulative effect in Kansas with 8 quarters, the demeaned placebo test had 51% power against the LASSO's 4%.
- Asked the same final-period question, the simultaneous band was less powerful than the pointwise test (Kansas, 80 quarters: 49% vs 65%).

**Supports.** The demeaned SC over the LASSO; using the method with short pre-periods; and the division of labour between the two kinds of interval: the pointwise test for a pre-specified period, the simultaneous band for statements about the whole path.

**Caveats.** The simulation predates `dsc_placebo()` and the RMSPE default: it used the relative scale and the test version of `dsc_symcal_boot()`. In California, the simultaneous band's power to detect an effect anywhere rose as the pre-period shortened, while coverage was unchanged; this is not yet explained.

### 3. Covariate adjustment (`covariate_adjustment_simulation.R`)

**Question.** Can the method adjust for a time-varying confounder, including a shock that starts at the same time as the treatment?

**Design.** California panel (38 states). A covariate shock with a known effect (gamma = 25) arrives in some donors at different times and in the treated unit together with the treatment, so an unadjusted analysis attributes the shock to the treatment. The pooled adjustment estimates gamma by two-way fixed effects on untreated observations only (all donor periods plus the treated unit's pre-period), subtracts gamma × covariate from every unit, and then runs the usual method. Each placebo re-estimates gamma without its own post period. The comparisons are no adjustment and unit-specific regressors (`controls`).

**Results.**

- **Bias removed:** cumulative bias -0.3% with adjustment vs +13.2% without (RMSE 8.3 vs 16.3). Gamma was recovered (median 25.7 vs true 25).
- **Still valid:** coverage 95-97%. Power was balanced (±25%: 26% / 34% final-period pointwise; 61% / 71% anywhere, simultaneous), whereas the unadjusted analysis had inflated power for increases and hid decreases.
- **Robust** to a moderately heterogeneous (proportional) shock effect: bias -0.7%, RMSE 10.0.
- **Unit-specific regressors were badly biased** (+55-65%, RMSE 150-200) when the covariate is zero before treatment.

**Supports.** The pooled adjustment estimated on untreated cells (`adjust_covariates`), rather than covariates as regressors in the synthetic control. The adjustment model was later changed from TWFE to interactive fixed effects; see simulation 5.

**Caveats.** Run with the bootstrap version and TWFE; it assumes the treatment does not affect the covariate, and the uncertainty in gamma is only partly reflected in the intervals.

### 4. Simultaneous-band ablation (`simultaneous_band_ablation.R`)

**Question.** The method began as a bootstrap: the maximum entropy bootstrap gave each period's interval its width, and placebos calibrated it. Does the bootstrap contribute anything once the placebos do the calibrating? And what should each gap be measured against?

**Design.** Kansas (49 states, 80 pre-treatment quarters) and California (38 states, 19 pre-treatment years), each state taking a turn as the treated unit, with effects injected (Kansas ±10%, California ±25%), at 95% and 90%. The same demeaned SC fits, donors, band construction (smoothed-median profile, donor maximum ratio, rank rule) and pointwise test are used throughout; only the **scale** changes: the bootstrap half-width, the pre-treatment RMSPE, the relative level (the synthetic value) and raw gaps. The bootstrap-scale bands reproduce `dsc_symcal_boot()` to 3e-11, so the comparison isolates what the bootstrap contributes.

**Results.**

- **All four scales were valid:** whole-path and final-period coverage 94-98% at 95%.
- **The bootstrap added nothing.** Without it, the relative scale matched or beat it almost everywhere. Simultaneous power at 95% (cumulative, -/+ effect) was 100% / 89% vs 39% / 39% in California and 88% / 84% vs 76% / 80% in Kansas. The bootstrap kept an edge only for small effects per period at 90% in California.
- **Relative was the most powerful scale and raw the least** (Kansas cumulative, 95%: 88% / 84% vs 22% / 20%; California: 100% / 89% vs 5% / 5%). With raw gaps, high-level states such as Alaska and North Dakota set the thresholds for everyone.
- **RMSPE was valid and close to relative in Kansas** (78% / 80% cumulative), but much weaker in California (42% / 47%).

**Supports.** Dropping the bootstrap: `dsc_placebo()` replaced `dsc_symcal_boot()`, running in seconds instead of about a minute. Not using raw gaps as the default. The relative scale was chosen as the original default on this evidence, before simulation 5 showed that its power advantage came with a robustness problem.

**Caveats.** In real panels the treated unit is always one of the states, so it is **exchangeable** with the donors by construction. This simulation cannot show what happens when the treated unit is unlike them.

### 5. Factor model (`factor_model_simulation.R`)

**Question.** Where do placebo methods fail? The real-panel simulations only test treated units that look like the donors.

**Design.** Fully simulated panels from a log-scale linear factor model: 2 common factors, either stationary (AR(1)) or random walks (trending), with unit-specific loadings and AR(1) errors. The treated unit is exchangeable with the donors, noisy, selected on its loadings, at the **edge** of the donors' loading range, or **outside** it. The grid also varies donors (20 / 40), pre-period length (10-80) and effects, with 12 post periods. A covariate block adds a covariate with a common additive effect, whose timing is either exogenous or related to the loadings. 98 cells, 100-200 replicates each; per-method failures are recorded. The methods compared are `dsc_placebo()` on each scale, `dsc_symcal_boot()`, and covariate adjustment by TWFE and by IFE.

**Results.**

- **Exchangeable treated units: all scales valid.**
- **Unlike the donors: the relative and raw scales failed, and the RMSPE scale stayed largely valid.** Final-period cumulative coverage at 95%:

  | Treated unit | RMSPE | Relative / raw |
  |---|---|---|
  | Selected on loadings | 95% | 81% |
  | Edge of the donors' range | 91% | 68% |
  | Outside the donors' range | 89% | 37% |

  The reason is that a treated unit the donors cannot reproduce also fits badly **before** treatment, so its pre-treatment RMSPE, and with it its intervals, grow. Percentage and dollar thresholds are set by the donors and are the same for every treated unit.
- **Power reversed relative to the real panels:** here RMSPE was the most powerful (exchangeable, cumulative -10%: 62% vs 45% for relative).
- **The bootstrap** was valid when exchangeable, better than relative at the edge, but worse than RMSPE.
- **Covariate adjustment:** IFE estimated gamma well (bias within about 1.7 of 11.9 when exchangeable), while TWFE failed with random-walk factors and long panels (bias -37 to -64 at T0 = 80; 10-14% of runs failed).

**Supports.** Changing the default scale from relative to RMSPE (tracker row 33): losing some power when the treated unit resembles the donors is a smaller cost than invalid intervals when it does not. Making IFE the default covariate-adjustment model (`adjust_model = "ife"`).

**Caveats.** Even the RMSPE scale under-covered with random-walk factors and long pre-periods for edge and outside units (56-82%). Outside the donors' range the estimator itself was biased (+10.9% cumulative), which no inference method can fix; this motivated simulations 6 and 7. The covariate block used only the relative scale, so the adjustment could not be judged at the edge.

### 6. Augmented comparison (`augmented_comparison_simulation.R`)

**Question.** Simulation 5 showed the estimator is biased when the treated unit lies outside the donors' range. Can an outcome model fix this, and do other established methods do better?

**Design.** The factor model of simulation 5, with the treated unit exchangeable, at the edge, or outside the range; stationary or random-walk factors; T0 = 20, 40, 80; 40 donors; 100 replicates (1,800 panels). The methods compared are:

- the demeaned SC (DSC);
- the **augmented demeaned SC** (ADSC: an IFE outcome model fitted on untreated cells, with the demeaned SC applied to its residuals; `augment = "ife"`);
- **fect** IFEct, with its own parametric inference and with placebo calibration;
- **augsynth** (ridge ASCM, standard and with fixed effects), with conformal inference.

Power was estimated only for methods whose intervals shift exactly with an added effect (the placebo-calibrated ones); for the others only bias, RMSE and size were assessed.

**Results.**

- **Bias:** augmentation reduced DSC's bias outside the range (random walk: +37% to +22%) but did not remove it. Standard augsynth had the best point estimates (lowest RMSE in 5 of 6 scenario × factor combinations).
- **Validity:** DSC with the RMSPE scale was the most reliable, valid with stationary factors even outside the range (98%). fect with placebo calibration on the RMSPE scale was valid with stationary factors (93-96%) and more powerful there.
- **Augmentation hurt inference:** it raised variance when the treated unit was ordinary (random-walk RMSE 25.7 vs 8.0), made RMSPE-calibrated intervals overconfident when bias remained (outside, random walk: 26% coverage), and its RMSPE-ratio test over-rejected (19-78%).
- **Other methods' own inference failed:** fect's parametric intervals covered only 25-65% off exchangeability, and augsynth's conformal average-effect test over-rejected even when exchangeable (13-29% at 5%).
- **No method was valid** for a treated unit outside the range with random-walk factors.

**Supports.** Keeping DSC + RMSPE as the default and augmentation as an option only (`augment = "none"`). Offering fect only with placebo calibration (`estimator = "fect"`), for stationary settings. The need for diagnostics that flag the remaining failure case.

**Caveats.** Run before the fix to how the treated unit's fect model is refitted (see *Notes on reproducibility*). Real-effect power for fect's and augsynth's own inference was not estimated.

### 7. Diagnostics validation (`diagnostics_validation.R`)

**Question.** Can the data tell the user when the method is likely to fail? Simulations 5 and 6 located the failure case (a treated unit beyond the donors' range with trending factors), but in practice the true loadings are unknown.

**Design.** The diagnostics in `dsc_placebo()` are computed on each of the 1,800 panels of simulation 6 and linked to that replicate's coverage. They are:

- **pre-treatment fit rank:** the treated unit's relative RMSPE against the placebos';
- **loading hull distance:** the treated unit's distance from the donors' convex hull in fect's estimated loadings, against each donor's distance from the others';
- **implied extrapolation drift**;
- **factor persistence:** the AR(1) coefficient of the estimated factors.

**Results.**

- **The extreme-unit flag** caught 100% of outside and 81% of edge treated units, with 4-5% false alarms for exchangeable units truly inside the hull (38-43% for those outside it by chance).
- **The flag alone does not predict failure:** flagged units with stationary factors still had 98% coverage.
- **The combined warning** (extreme unit AND trending factors, persistence above 0.95) separated failing runs: 74% coverage when warned vs 96% when not (DSC, RMSPE scale).
- **Kansas and California are not flagged** (inside the donors' range; implied drift 1.4% and 3.5%), although both have trending factors.

**Supports.** The `$diagnostics` output, the rule that the warning requires both conditions, and the 0.95 persistence threshold.

**Caveats.** The persistence check misclassifies some stationary factors (15% at 0.95, more with long pre-periods), and the thresholds were calibrated on one factor-model design.

### 8. Toy cases (`diagnostics_toy_cases.R`)

**Question.** Does the diagnostics-based recommendation send the user to the right estimator?

**Design.** The script holds the helpers (`run_toy_case()` and the plotting and table functions) that the diagnostics write-up, `docs/diagnostics/diagnostics_toy_cases.qmd`, uses to simulate six toy panels from the factor model. The panels cover exchangeable, edge and outside treated units, each with stationary or random-walk factors. Each panel is fitted with `dsc_placebo()` (RMSPE scale, with diagnostics) and then re-run as `$diagnostics$recommendation` advises.

**Results.**

- **The recommendation was correct in all six cases:** it suggests fect only for stationary edge and outside units.
- **Following it helped.** Stationary edge: cumulative error -8.6% to +1.1%, half-width 32.1% to 5.2%. Stationary outside: error -20.6% to -18.4%, half-width 58.2% to 16.0%.
- **The fect refit matches a direct fect call** (difference 0).

**Supports.** `estimator = "fect"` and the recommendation rule (fect + RMSPE scale only when the treated unit is extreme and the factors do not trend). Coverage for these cases comes from simulation 6 (91% for the stationary outside case).

### 9. Oregon fentanyl validation (`oregon_fentanyl_validation.R`)

**Question.** Does the covariate adjustment work on real data, for the kind of confounder it was designed for? In Oregon, fentanyl arrived in the drug supply just as drug possession was decriminalized (2021), and an unadjusted synthetic control attributes the resulting rise in overdose deaths to the policy (see `docs/oregon_fentanyl/`).

**Design.** The panel of Zoorob et al. (2024) (`data/oregon_fentanyl/`): half-yearly overdose deaths per 100k and the fentanyl share of seized drugs, 2008-2022, without Washington and Oregon. Each of the 49 remaining units (48 states and DC) is a fake treated unit from 2021 H1, with the other 48 as donors. An Oregon-like shock is added after treatment: the unit's fentanyl share rises by Oregon's actual increase (up to +40 points), and its death rate by gamma_true × the added share. Scenarios: no shock; gamma_true = 0.155 (about the gamma the adjustment estimates, so the adjustment model is right); gamma_true = 0.12 (gamma over-estimated by 0.035); true effects 0 and ±25%. Methods: `dsc_placebo()` unadjusted and adjusted (IFE and TWFE; identical here, because cross-validation chose no factors for 48 of the 49 units). Because weights, gamma-hat and interval widths depend only on pre-treatment data and donors, each scenario shifts the gap exactly, so each unit is fitted once and every scenario is evaluated from that fit.

**Results** (95%, RMSPE scale, true effect zero):

| Scenario | Method | Cumulative bias | Coverage (cumulative / simultaneous) | False positives (cumulative / simultaneous) |
|---|---|---|---|---|
| No shock | Unadjusted | -1.9% | 96% / 96% | 4% / 4% |
| No shock | Adjusted | -1.8% | 96% / 96% | 4% / 4% |
| Shock, gamma correct | Unadjusted | +26.3% | 61% / 61% | 39% / 39% |
| Shock, gamma correct | Adjusted | -1.7% | 96% / 96% | 4% / 4% |
| Shock, gamma over-estimated | Unadjusted | +19.9% | 67% / 65% | 33% / 35% |
| Shock, gamma over-estimated | Adjusted | -8.1% | 100% / 96% | 0% / 4% |

- Without a shock, adjusting costs nothing in validity; with an Oregon-like shock, the unadjusted method finds a non-existent effect in a third or more of the units, and the adjusted method stays valid, even when gamma is mis-estimated (then with some bias).
- The price is power: without a shock, cumulative power at -25% / +25% falls from 55% / 41% (unadjusted) to 29% / 27% (adjusted). With the shock, the unadjusted method's power is distorted by its bias (8% / 80%); the adjusted method's is unchanged (29% / 29%).

**Supports.** The covariate adjustment (`adjust_covariates`) as the remedy for a measured confounder that arrives with the treatment, on real data rather than the simulated shock of simulation 3; and the Oregon case study's conclusion that the unadjusted effect is largely the arrival of fentanyl.

**Caveats.** The shock acts through the covariate with a common linear effect, as the adjustment assumes; heterogeneous effects, a covariate affected by the treatment, and unmeasured confounders are not tested. Ranks within one panel are not independent across fake treated units.

### 10. RMSPE guard (`rmspe_guard_validation.R`)

**Question.** With a short pre-period and many donors, the demeaned synthetic control can reproduce some placebos exactly before treatment. Their pre-treatment RMSPE is then about zero, their standardized gaps near-infinite, and so are the RMSPE-scale intervals (found on the Texas prison panel of Cunningham (2021): 8 pre-treatment years, 49 donors, 95% interval about ±9 billion). Does flooring each unit's pre-treatment RMSPE at a fraction of the median (`rmspe_floor`) fix this without breaking validity?

**Design.** Floors of 0 (no guard), 0.1 and 0.25, with the relative scale as a reference, on: Kansas (pre-period cut to 80 / 20 / 8 quarters) and California (19 / 8 / 5 years), each untreated state in turn as the fake treated unit; the Texas prison panel of Cunningham (2021) (8 pre-years; Texas excluded; downloaded from the book's GitHub repository); and the factor model of simulation 6 with T0 = 10 / 20 / 40 (exchangeable, edge and outside treated units, stationary or random-walk factors, 40 donors, 100 replicates). Each panel is fitted once; the calibration is applied for each floor by a function checked to reproduce `dsc_placebo()` exactly. Power by exact shifting (effects of ±10% for Kansas and the factor model, ±25% for California and Texas).

**Results** (cumulative effect, 95%):

| Design | Coverage, no guard / floor 0.25 | Median half-width (% of outcome), no guard / floor 0.25 | Simultaneous power at a negative effect, no guard / floor 0.25 |
|---|---|---|---|
| Kansas, 80 quarters | 96% / 96% | 6% / 6% | 80% / 80% |
| Kansas, 8 quarters | 94% / 96% | 2,090% / 11% | 6% / 82% |
| California, 8 years | 97% / 97% | 626,550% / 53% | 3% / 50% |
| California, 5 years | 100% / 97% | ~4 million % / 86% | 0% / 47% |
| Texas, 8 years | 94% / 96% | ~4 million % / 157% | 6% / 25% |
| Factor model, exchangeable, T0 = 10 | 94-96% / 92-96% | 30-35% / 16-17% | 19-22% / 43-45% |

- With long enough pre-periods (Kansas 20+ quarters, California 19 years, factor model T0 = 20 or 40) the results are **identical** with and without the guard.
- With short pre-periods, the unguarded RMSPE scale was valid only because its intervals were effectively infinite. With the floor at 0.25, coverage stayed at 92-98% in the real panels and for exchangeable units, and power rose sharply. A floor of 0.1 was in between.
- For edge and outside treated units, the guard left coverage as it was.

**Supports.** The guard, with `rmspe_floor = 0.25` as the default.

**Caveats.** The lowest coverage (92%) was for exchangeable units with random-walk factors and T0 = 10 (100 replicates; standard error about 2 points). With very short pre-periods the intervals stay wide (Texas 157%). One fake treated unit could not be fitted at all (Vermont in the Texas panel, whose pre-treatment outcome is zero throughout).

### 11. Ridge augmentation (`ridge_edge_validation.R`, `ridge_small_sample_validation.R`, `ridge_covariate_validation.R`, `ridge_function_check.R`)

**Question.** Simulations 5 and 6 found one case where no placebo calibration rescued the demeaned synthetic control: a treated unit at or beyond the edge of the donors' range with trending factors. The estimate itself is biased there, because non-negative weights summing to one can't reproduce the treated unit's path. Does a ridge-augmented estimator fix it, and does it hold up elsewhere?

**The estimator** (`augment = "ridge"`). In each post-treatment period, the counterfactual is the demeaned synthetic control plus (x_T − Σ w_j x_j)′β_t. Here x_j is unit j's demeaned pre-treatment path, and β_t a ridge regression, across the donors, of their period-t outcome on those paths. The correction extrapolates the pre-treatment imbalance the weights leave, and is about zero when they match. The penalty is chosen by cross-validation inside the pre-period: each donor predicts its last pre-treatment periods from the earlier ones. Inference is the usual in-space placebo calibration, each placebo refitting its own correction. With `adjust_covariates`, the outcomes are adjusted first.

**Designs** (all on the factor-model generator of simulation 5, 100 replicates per cell, 12 post periods):
- `ridge_edge_validation.R`: treated unit exchangeable, at the edge or outside; stationary or random-walk factors; 19 or 39 donors; 40 or 80 pre-periods. Also tests in-time placebo calibration (the treated unit's own backdated windows) and the wider of the two calibrations. 2,400 replicates.
- `ridge_small_sample_validation.R`: five kinds of treated unit (adding selected and noisy); 9, 14 or 19 donors; 10 or 20 pre-periods (6,000 replicates). Also every untreated unit of Kansas, California and West Germany as a fake treated unit, the three real treated units, and agreement with `augsynth`.
- `ridge_covariate_validation.R`: the confounded covariate of simulation 5 (it switches on for the treated unit at treatment, adding about 8%, and for donors at dates that depend on their loadings); no, TWFE or IFE adjustment, each with or without ridge; 14 or 19 donors; 20 or 40 pre-periods (3,200 replicates).

**Results** (coverage at 90%):

| Design | Without ridge | With ridge |
|---|---|---|
| Edge, random walk, 19-39 donors, 40-80 pre-periods (cumulative) | 62% | 91% |
| Outside, random walk, 19-39 donors, 40-80 pre-periods (cumulative) | 49% | 91% |
| Outside, random walk, 9-19 donors, 20 pre-periods (final period) | 76% | 90% |
| Exchangeable, no covariate (all of the above designs) | 89-93% | 88-94% |
| Fake treated units: Kansas / California / West Germany | 90 / 92 / 88% | 90 / 92 / 88% |
| Confounded covariate, IFE adjustment, outside, random walk (final / cumulative) | 73 / 73% | 89 / 93% |

- **Bias and accuracy.** Bias fell by three-quarters or more for edge and outside units (outside, random walk: +38% to +9%), and RMSE by half or more in every scenario. Intervals were narrower. Power for exchangeable units was the same as without ridge.
- **Small samples.** Ridge stayed valid with 9 donors and 10 pre-periods (89-99% at 90% averaged over pool sizes). With about 10 pre-periods the cross-validated penalty was less stable, and for a noisy treated unit with random-walk factors the RMSE was higher than without ridge (20% vs 16%).
- **Stationary factors.** Ridge was conservative for edge and outside units (99-100%).
- **In-time calibration** on its own under-covered (81-83%); taking the wider of the in-space and in-time thresholds reached 84-87% for the plain estimator but doubled the widths. It was not adopted.
- **Covariates.** IFE adjustment followed by ridge was the best combination in every scenario (89-99% final period, 86-99% cumulative). Without adjustment, cumulative coverage fell to 69-82%; ridge alone did not remove the confounding. TWFE's estimate of the covariate effect was badly biased (RMSE 79-317% of the true value), because the covariate's timing depends on the factors.
- **Real cases.** Estimates barely change: Kansas −2,006 vs −1,891 in the final quarter, California −17.3 vs −17.4, West Germany −3.51 vs −3.49.
- **augsynth.** With the same penalty, `augsynth(progfunc = "Ridge", fixedeff = TRUE)` reproduces the estimates to within 0.00003; with its own penalty choice they differ (Kansas −2,852).
- **The function.** `ridge_function_check.R` compares `dsc_placebo(augment = "ridge")` with the prototype on the three real panels and 48 simulated panels (including both covariate adjustments). They chose the same penalty in every panel, and the effects, pointwise intervals, simultaneous bands and cumulative intervals agreed to within 0.0008% of the outcome (median 10⁻⁸%), the difference being the weight solver's tolerance. The covariate effect estimates were identical. So the three simulations above apply to `dsc_placebo()` itself.

**Supports.** `augment = "ridge"` as an option, combinable with `adjust_covariates`, and the diagnostics' recommendation of it for an extreme treated unit with trending factors. The default stays `augment = "none"`: for treated units like the donors the two are equally valid, and ridge adds little.

**Caveats.** The simulations use the prototype in `ridge_validation_helpers.R` (quadprog weights rather than CVXR), not `dsc_placebo()` itself; the function check covers that. Some bias remains for a treated unit well outside the donors' range with random-walk factors. Exogenous covariate timing and `augment = "ife"` with covariates were not tested.

## Notes on reproducibility

- `augmented_comparison_simulation.R` and `diagnostics_validation.R` were run before a fix to how `dsc_placebo()` refits fect models (2026-10-15: the treated unit's model is now refitted with the cross-validated number of factors, as each placebo's is). Rerunning them with the current functions will give slightly different numbers.
- `factor_model_simulation.R`'s covariate block uses `scale = "relative"` explicitly: it was run when that was `dsc_placebo()`'s default.
- `pre_period_length_simulation.R`, `covariate_adjustment_simulation.R` and `simultaneous_band_ablation.R` predate the RMSPE default; the scales they used are stated above.
- All simulations before `rmspe_guard_validation.R` were run without the guard (`rmspe_floor`, added 2026-10-15). It changes nothing where no unit is fitted almost exactly before treatment, which held for every design with 20 or more pre-periods in simulation 10; the T0 = 10 cells of `factor_model_simulation.R` would change.

Full details, including the simulations not kept here, are in `../experimental_methods_tracker.md`.

## References

- Cunningham, S. (2021). *Causal Inference: The Mixtape*. New Haven: Yale University Press. Data: https://github.com/scunning1975/mixtape (`texas.dta`).
- Zoorob, M. J., Park, J. N., Kral, A. H., Lambdin, B. H., and del Pozo, B. (2024). Drug decriminalization, fentanyl, and fatal overdoses in Oregon. *JAMA Network Open*, 7(9), e2431612.
- Firpo, S., and Possebom, V. (2018). Synthetic control method: Inference, sensitivity analysis and confidence sets. *Journal of Causal Inference*, 6(2), 1-26.
