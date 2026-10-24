# Write-ups

*Written by Claude, an AI model made by Anthropic, at the direction of Jeffrey Milliman.*

Each folder holds one Quarto document and the figures it uses (`figures/`), so a folder works wherever it is placed (`reference/` holds a plain Markdown file, and `readme_figures/` the figures of the main README).
Each renders to a GitHub-flavoured Markdown file (`.md`), which is what GitHub displays and what the links point to, and to a self-contained HTML file (`embed-resources: true`), which is not kept in the repository.

| Folder | Document | Topic |
|---|---|---|
| `reference/` | `function_reference.md` | Every argument of `dsc_placebo()` and the plotting functions, and every element of the output (plain Markdown; no rendering needed) |
| `case_study/` | `kansas_case_study.qmd` | **Start here:** the Kansas tax cuts analysed with `dsc_placebo()`, showing its features, diagnostics and sensitivity checks |
| `oregon_fentanyl/` | `oregon_fentanyl_case_study.qmd` | Oregon's drug decriminalization and the arrival of fentanyl: covariate adjustment for a confounder that arrives with the treatment, validated on the same panel (reads `data/oregon_fentanyl/`) |
| `diagnostics/` | `diagnostics_toy_cases.qmd` | When placebo inference fails (extreme treated units, trending factors), and the diagnostics in `dsc_placebo()` |
| `placebo/` | `placebo_calibrated_bands.qmd` | The placebo-calibrated method (pointwise intervals and simultaneous bands), explained on the default RMSPE scale, with validation and the alternative relative and raw scales |
| `validation/` | `validation_summary.md` | How the experimental functions (`dsc_placebo_lto()`, `synth_power()`, `covariate_trajectory_check()`) were checked by simulation (plain Markdown; no rendering needed) |
| `bootstrap/` | `fixed_weight_calibrated_bootstrap.qmd` | The earlier approach: fixed-weight maximum entropy bootstrap with donor calibration (kept as part of the method's development; the placebo method is preferred) |

To render both formats, run `quarto render docs/<folder>/<file>.qmd` from the project root (add `--to gfm` for the Markdown only). RStudio's Render button produces only the HTML version.

The code chunks are not run when rendering (`eval: false`); they show how each figure and result was produced.
Their paths (e.g. `synthetic_control_functions/synthetic_control_functions_claude.R`, `simulations/diagnostics_toy_cases.R`) are relative to the
project root, so run copied code from there.
