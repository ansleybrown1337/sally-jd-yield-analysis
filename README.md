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
This new workflow reproduces the logic of Sally’s SAS variety-trial analysis but extends it to a more rigorous and transparent multi-environment framework. It first analyzes each trial site (or environment) individually using spatial mixed models that account for field position (row and column) and block structure, automatically testing multiple spatial covariance types (exponential, spherical, Gaussian, anisotropic exponential) and selecting the best by AICc. From each best-fitting model, it extracts spatially adjusted least-squares means (LS-means) and their variances—these represent each variety’s yield corrected for within-field variability. The second stage then combines these site-level LS-means across all environments using an inverse-variance meta-analysis, so sites with more precise estimates contribute proportionally more weight. Tukey-adjusted pairwise comparisons and compact letter displays (CLDs) summarize which varieties differ significantly overall. A complementary mixed model using BLUPs (Best Linear Unbiased Predictions) estimates variety effects directly across all environments in a single joint model, accounting for random environment and genotype×environment interaction. Together, this framework increases statistical power, handles spatial autocorrelation explicitly, and provides both adjusted means (for reporting) and BLUPs (for breeding insight), offering a more modern, reproducible alternative to traditional per-site ANOVA and LSD methods.

### Goals

* Increase statistical power through appropriate spatial modeling.
* Generate comparable and reproducible least-squares means or BLUPs across sites.
* Provide flexibility for both **fixed-effects reporting (LS-means)** and **random-effects inference (BLUPs)**.
* Eliminate use of averaged LSDs; replace with Tukey-adjusted contrasts or compact letter displays.
* Extend analyses to multi-site and multi-year trials with robust weighting.


### Key Terms
| Term                                       | Definition                                                                                                                                                               |
| ------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Trial site**                             | A physical field location where a variety trial was planted and measured (e.g., “Fort Collins 2024”). Each site has its own soil type, management, and local conditions. |
| **Environment (env)**                      | The combination of a trial site and year (e.g., “FortCollins_2024”). It represents a unique set of environmental conditions under which genotypes were tested.           |
| **Replication (rep)**                      | A complete experimental block or repetition of all varieties within a trial. Replicates help estimate field variability and improve statistical power.                   |
| **Row and Column (row, col)**              | The spatial coordinates of each plot in the field grid. Used in the spatial model to correct for field trends such as fertility or moisture gradients.                   |
| **Entry (entry)**                          | The genotype, cultivar, or variety being evaluated (e.g., “Millet_14”).                                                                                                  |
| **Yield (yield)**                          | The observed response variable (e.g., grain yield, biomass, protein) measured on each plot.                                                                              |
| **Covariance structure (cov)**             | A model describing how spatial correlation declines with distance (e.g., exponential, spherical). Determines how residual variation is modeled across the grid.          |
| **LS-mean (Least Squares Mean)**           | The adjusted mean yield per variety after accounting for random and spatial effects within a site. Equivalent to the SAS “adjusted mean.”                                |
| **BLUP (Best Linear Unbiased Prediction)** | A prediction of a random effect (e.g., variety performance) that accounts for shrinkage toward the mean. BLUPs are often used to predict true genetic potential.         |
| **Compact Letter Display (CLD)**           | A letter-based summary of pairwise comparisons, where varieties sharing a letter are not significantly different.                                                        |
| **Inverse-variance weighting**             | A method that gives more weight to site means with higher precision (smaller variance) when combining results across environments.                                       |

### Input Data Structure
| Column  | Type                 | Description                                                                                |
| ------- | -------------------- | ------------------------------------------------------------------------------------------ |
| `site`  | character            | Site name or code (e.g., “FortCollins”).                                                   |
| `year`  | integer              | Year of trial (e.g., 2024).                                                                |
| `env`   | character            | Concatenation of site and year (e.g., “FortCollins_2024”). Used as the environment factor. |
| `rep`   | integer              | Replicate number (e.g., 1–4).                                                              |
| `row`   | integer              | Field grid row number.                                                                     |
| `col`   | integer              | Field grid column number.                                                                  |
| `entry` | character or integer | Variety identifier (e.g., “M1”, “14”, “HybridA”).                                          |
| `yield` | numeric              | Measured trait (e.g., grain yield in kg/ha or bu/ac).                                      |

*Example Data Header*
| site        | year | env              | rep | row | col | entry | yield  |
| ----------- | ---- | ---------------- | --- | --- | --- | ----- | ------ |
| FortCollins | 2024 | FortCollins_2024 | 1   | 1   | 3   | 14    | 5023.4 |
| FortCollins | 2024 | FortCollins_2024 | 1   | 1   | 4   | 9     | 4872.1 |
| Fruita      | 2024 | Fruita_2024      | 2   | 3   | 7   | 18    | 4630.9 |
| Fruita      | 2024 | Fruita_2024      | 2   | 3   | 8   | 5     | 4459.7 |


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

**Key Outputs:**

* Environment-weighted entry LS-means
* Tukey-adjusted pairwise comparisons (controlled FWER)
* Optional compact letter display (CLD) for report tables

### Stage 2A — Random-Effects MET Model (Optional)

A unified one-stage model can also estimate entry BLUPs across all environments:


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

* `lsm_stage1.csv`— per-environment LS-means (+ SE, chosen covariance)
* `lsm_across.csv` — inverse-variance weighted across-env LS-means
* `diffs_across.csv` — Tukey-adjusted pairwise contrasts across environments
* `cld_across.csv` — compact letter display (groups) across environments
* `entry_blups.csv` — entry BLUPs (and SE) from the MET model
* **`readme.md` (this file)** — Documentation of methodology and improvements.

---

## Next Steps

* Validate covariance model selection on several representative sites.
* Integrate optional R-based variogram diagnostics for visual verification.
* Add heritability and stability metrics for BLUP outputs.
* Prepare publication-ready summary tables of entry means and grouping letters for extension reports.

