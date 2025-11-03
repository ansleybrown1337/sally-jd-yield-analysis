# Variety Trial Spatial Modeling Workflow
Created by: AJ Brown, Ag Data Scientist, Agricultural Water Quality Program
Date: 3 November 2025

## Overview

This project implements an updated spatial mixed-model workflow for analyzing multi-site and multi-year crop variety trial data. The goal is to improve statistical rigor, precision, and interpretability of entry (variety) means while preserving compatibility with existing SAS-based workflows used by CSU Crops Testing Program (via Sally Jones-Diamond).

The workflow addresses limitations of the previous approach by explicitly modeling spatial autocorrelation within trial fields and by harmonizing analyses across sites and years through mixed-model integration or meta-analysis.

---

## Previous Approach ("Sally’s Original Code")

### Summary

Sally’s code analyzed each site-year ("location") independently using a combination of `PROC GLM` and `PROC MIXED` in SAS. It relied on Fisher’s LSD (α = 0.30) for entry comparison and did not consistently propagate spatial variance structures across environments.

### Key Features

* **PROC GLM:** Fit fixed-effects model (entry + replication) ignoring spatial correlation.
* **PROC MIXED:** Added spatial covariance structure using the `repeated / type=sp(exp)` statement to reduce residual variance.
* **Model Selection:** Compared spatial covariance forms (e.g., exponential, spherical, Gaussian) by AIC/AICc and log-likelihood improvement.
* **Reporting:** Exported pairwise `diffs` and averaged individual LSD values to report a single mean LSD per trial.

### Limitations

| Limitation                                    | Impact                                                                |
| --------------------------------------------- | --------------------------------------------------------------------- |
| Averaged LSDs lack inferential meaning        | Cannot control type I error or ensure consistent thresholds           |
| Single environment focus                      | Cannot estimate entry × environment (G×E) effects or cross-site means |
| Uniform α = 0.30                              | Inflated false positives, limited reproducibility                     |
| No weighting by precision across sites        | Overweights noisy trials                                              |
| Identical spatial parameters assumed per site | Ignores real variation in field heterogeneity                         |

---

## New Approach (2025 Workflow)

### Goals

* Increase statistical power through appropriate spatial modeling.
* Generate comparable and reproducible least-squares means or BLUPs across sites.
* Provide flexibility for both **fixed-effects reporting (LS-means)** and **random-effects inference (BLUPs)**.
* Eliminate use of averaged LSDs; replace with Tukey-adjusted contrasts or compact letter displays.
* Extend analyses to multi-site and multi-year trials with robust weighting.

### Stage 1 — Within-Environment Spatial Modeling

Each site-year ("environment") is analyzed separately using candidate spatial covariance structures.

**Models tested:**

* `sp(exp)` — isotropic exponential
* `sp(expa)` — anisotropic exponential (different range per axis)
* `sp(sph)` — spherical
* `sp(gau)` — Gaussian

The model with the **lowest AICc** is selected per environment.

**Output per environment:**

* LS-means for each entry
* Standard errors and variance of LS-means
* Chosen covariance type

These results are stored in a tidy dataset for downstream meta-analysis.

### Stage 2B — Fixed-Effects Meta-Analysis Across Environments

Stage 1 LS-means are combined across sites (and optionally years) using inverse-variance weighting.

```sas
proc mixed method=reml;
  class env entry;
  model estimate = entry / s;
  random env;
  weight 1/var_lsmean;
  lsmeans entry / pdiff=all adjust=tukey cl;
run;
```

**Key Outputs:**

* Environment-weighted entry LS-means
* Tukey-adjusted pairwise comparisons (controlled FWER)
* Optional compact letter display (CLD) for report tables

### Stage 2A — Random-Effects MET Model (Optional)

A unified one-stage model can also estimate entry BLUPs across all environments:

```sas
proc mixed method=reml;
  class env rep entry;
  model yield = ;
  random env rep(env) entry entry*env;
  repeated / subject=env type=sp(expa) (row col) group=env;
run;
```

**Advantages:**

* Accounts for unbalanced entries across sites/years.
* Produces BLUPs and variance components (σ²_entry, σ²_entry×env).
* Allows calculation of heritability and stability indices.

---

## Improvements Over Previous Workflow

| Improvement                                 | Description                                             | Benefit                                                   |
| ------------------------------------------- | ------------------------------------------------------- | --------------------------------------------------------- |
| **Model-based LSD replacement**             | Tukey-adjusted pairwise contrasts replace averaged LSDs | Controls experiment-wise error; interpretable letters     |
| **AICc-based model selection**              | Automated per-site evaluation of spatial covariance     | Objectively identifies best-fitting correlation structure |
| **Weighted meta-analysis**                  | Combines LS-means using inverse variance                | Properly balances high- vs low-precision sites            |
| **Optional BLUP framework**                 | Random entry model spans unbalanced trials              | Increases stability and predictive reliability            |
| **Environment-specific spatial parameters** | `group=env` option in `repeated`                        | Captures heterogeneous field variability                  |
| **Modular SAS macros**                      | Stage 1/Stage 2 separation                              | Easy updates, reproducible code, transparent output       |

---

## Desired Goals (per Sally’s notes) and Achievements

| Goal                                | Addressed By                                             |
| ----------------------------------- | -------------------------------------------------------- |
| Reduce uncontrolled field variation | Spatial covariance (`sp(expa)` and AICc selection)       |
| Improve power and precision         | Mixed-model REML estimation with spatial residuals       |
| Facilitate clear reporting          | LS-means or BLUPs + Tukey or CLD tables                  |
| Maintain simplicity across sites    | Unified Stage 2 analysis with inverse-variance weighting |
| Enable cross-site inference         | Entry × environment modeling and weighting               |

---

## Deliverables

* **`spatial_trial_analysis.sas`** — Main script containing Stage 1 and Stage 2 macros.
* **`lsm_stage1.csv`** — Per-environment LS-means and variances.
* **`lsm_across.csv`** — Weighted LS-means and Tukey comparisons across environments.
* **`entry_blups.csv`** — Entry BLUPs from MET model.
* **`readme.md` (this file)** — Documentation of methodology and improvements.

---

## Next Steps

* Validate covariance model selection on several representative sites.
* Integrate optional R-based variogram diagnostics for visual verification.
* Add heritability and stability metrics for BLUP outputs.
* Prepare publication-ready summary tables of entry means and grouping letters for extension reports.

