# Variety Trial Spatial Modeling Workflow
Created by: AJ Brown, Ag Data Scientist, Agricultural Water Quality Program  
Updated: 23 April 2026

## Overview

This repository contains an R implementation of a two stage variety trial analysis workflow designed to mirror and improve upon an existing SAS based process used for crop variety testing. The pipeline fits within environment models first, then summarizes variety performance across environments using weighted least squares means and a one stage BLUP model.

The current code supports both simulated and real input files, writes a reproducible set of tabular outputs, and can generate rendered HTML summaries from the frontend R Markdown workflow.

## What is in this repository

```text
code/
  r_equivalent.R                  Main R workflow
  variety_trial_frontend.Rmd      Rendered reporting frontend
  Allsites-variety_trial_frontend.html
  Millet2025_trial_frontend.html
  validation.Rmd
  validation.html
  sim_data.r
  aj_stats.sas                    Legacy or comparison SAS code

sim_data/
  Example simulated input files

sim_output/
  Example simulated outputs produced by the R workflow

real_data/
  Real input files, intended to stay out of the public repo

real_output/
  Real outputs, intended to stay out of the public repo

docs/
  Background documents and collaborator materials
```

## Workflow summary

The implemented workflow has three main analytical components.

### Stage 1, within environment analysis

Each environment, defined as a site year combination in `env`, is analyzed separately after dropping rows with missing yield values.

The code attempts spatial mixed models first using `nlme::lme` with:

- fixed effect: `yield ~ entry`
- random effect: `~ 1 | rep`
- spatial correlation over field coordinates `row` and `col`
- REML estimation
- model comparison by AICc

Candidate covariance labels are:

- `expa`
- `exp`
- `sph`
- `gau`

Important implementation note: in the current R code, `expa` and `exp` are both mapped to `nlme::corExp(...)`. That means the present implementation does **not** yet fit a distinct anisotropic exponential model, even though that is part of the intended conceptual design.

For each environment, the workflow extracts:

- LS means by entry
- standard errors and confidence intervals
- the selected covariance label
- plot level fitted values and residual diagnostics
- a model selection table covering both successful and failed fits
- environment level RMSE and CV summaries

### Fallback ladder in Stage 1

If all spatial models fail for an environment, the code falls back in this order:

1. RCBD mixed model with `lme4::lmer(yield ~ entry + (1 | rep))`
2. fixed effects model with `lm(yield ~ entry)` when the RCBD fit is not feasible or also fails

This is an important part of the actual workflow and should be understood as part of the intended robustness of the pipeline.

### Stage 2B, across environment LS means

If two or more environments contribute usable Stage 1 LS means with positive variances, the code combines them using inverse variance weighting.

- with at least two environments, it uses `nlme::lme` on Stage 1 LS means with a random environment intercept and `varFixed(~ var_lsmean)`
- with only one usable environment, it skips meta analysis and instead computes within environment LS means and Tukey summaries directly from the raw data

The Stage 2B outputs include:

- overall adjusted LS means across environments
- Tukey adjusted pairwise comparisons
- compact letter displays
- approximate `LSD_0.05` and `LSD_0.30` columns derived from unadjusted pairwise standard errors and degrees of freedom

Those LSD columns are still being reported for continuity with historical practice, but the inferential emphasis in this workflow is on the Tukey adjusted results and compact letter displays.

### Stage 2A, one stage BLUP model

The workflow also fits a one stage mixed model for entry BLUPs.

When multiple environments are available, the intended model is:

```r
yield ~ 1 +
  (1 | env) +
  (1 | env:rep) +
  (1 | entry) +
  (1 | entry:env)
```

When only one environment is available, the BLUP model is simplified accordingly. If the mixed model fails, the code falls back to `lm(yield ~ entry)` and reports `emmeans` summaries instead.

## Input requirements

The workflow expects a CSV file with at least these columns:

| Column | Description |
|---|---|
| `site` | Site name or code |
| `year` | Trial year |
| `env` | Environment identifier, usually site plus year |
| `rep` | Replication or block |
| `row` | Plot row coordinate |
| `col` | Plot column coordinate |
| `entry` | Variety, genotype, or treatment identifier |
| `yield` | Response variable used for modeling |

Rows with missing `yield` are excluded from model fitting.

## Current default file paths in the script

The main script currently defines these defaults:

```r
DATA_MODE    <- "sim"
SIM_IN_CSV   <- "sim_data/trial_sim.csv"
SIM_OUT_DIR  <- "sim_output"
REAL_IN_CSV  <- "real_data/MilletVT2025.csv"
REAL_OUT_DIR <- "real_output"
```

The main function is:

```r
analyze_trial(in_csv = SIM_IN_CSV, out_dir = SIM_OUT_DIR, alpha = 0.05)
```

## Primary outputs actually produced by the code

The current `r_equivalent.R` workflow writes all of the following:

| File | Description |
|---|---|
| `lsm_stage1.csv` | Per environment LS means and variances from Stage 1 |
| `lsm_across.csv` | Across environment LS means, confidence intervals, and LSD helper columns |
| `diffs_across.csv` | Pairwise contrasts from Stage 2B |
| `cld_across.csv` | Compact letter display table |
| `entry_blups.csv` | Entry BLUPs or fallback summaries |
| `stage1_model_selection.csv` | Full Stage 1 model ranking table, including failures |
| `model_specs_stage1.csv` | Explicit record of the fitted Stage 1 model specification per environment |
| `diagnostics_stage1.csv` | Plot level fitted values, residuals, and outlier flags |
| `cv_by_env.csv` | Environment level mean yield, RMSE, and CV |
| `run_report.txt` | Run audit trail and decision summary |

## Frontend reporting layer

The repository also includes an R Markdown reporting frontend and rendered HTML reports. The HTML output indicates the report currently summarizes:

- workflow overview
- backend run and configuration summary
- data characteristics
- Stage 1 model results
- Stage 2B LS means and CLD
- Stage 2A BLUPs
- CV summary
- field layout and heatmap style plots
- backend run report
- output descriptions

This means the repository is not only an analysis engine, it also includes a report generation layer for communicating results to collaborators and reviewers.

## Key differences from the older SAS style workflow

Compared with a traditional per site GLM or MIXED workflow that emphasizes Fisher style LSD reporting, this implementation adds several improvements.

- Spatial model comparison is automated by AICc within environment.
- If spatial fitting fails, the code degrades gracefully rather than stopping.
- Across environment summaries are weighted by Stage 1 uncertainty.
- Tukey adjusted inference and compact letter displays are used for reporting.
- A one stage BLUP model is available for breeder style interpretation.
- Diagnostics, model specifications, and run audit files are exported explicitly.

## Important current caveats

These points are especially important for anyone reusing or publishing from this workflow.

- `expa` is currently only a label and is not distinct from `exp` in the present implementation.
- The code assumes a single random intercept for `rep` in Stage 1, with no explicit nested or crossed structure beyond the environment specific fits.
- Confidence intervals and degrees of freedom come from standard approximations in `nlme`, `lme4`, and `emmeans`, and may be unstable in sparse or highly unbalanced settings.
- The approximate LSD columns are convenience summaries, not the main inferential basis.
- Real collaborator documents in `docs/from sally/` should be reviewed before making the repository public.

## Typical usage

Example run on simulated data:

```r
source("code/r_equivalent.R")
analyze_trial(
  in_csv = "sim_data/trial_sim.csv",
  out_dir = "sim_output",
  alpha = 0.05
)
```

Example run on real data:

```r
source("code/r_equivalent.R")
analyze_trial(
  in_csv = "real_data/MilletVT2025.csv",
  out_dir = "real_output",
  alpha = 0.05
)
```

## Suggested next development steps

The code and repository suggest the following near term priorities.

- Implement a truly anisotropic spatial option if `expa` is meant to remain in the candidate set.
- Decide whether `LSD_0.30` should remain for legacy compatibility or move to a historical appendix only.
- Add a short public facing note explaining simulated versus real example content.
- Review `docs/from sally/` before publication, since those files may still contain private collaborator material.
- Consider adding a minimal reproducible example dataset and a one command render instruction for the frontend report.

## Package dependencies

The script loads:

- `dplyr`
- `tidyr`
- `purrr`
- `nlme`
- `lme4`
- `emmeans`
- `ggplot2`
- `multcompView`
- `broom`
- `tibble`

## Repository publication note

Before making the repository public, the real data and real outputs should remain ignored by Git. At minimum, `real_data/` and `real_output/` should be excluded. Depending on your privacy needs, you may also want to review whether `docs/from sally/` belongs in the public version.