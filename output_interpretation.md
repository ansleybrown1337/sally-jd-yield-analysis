# Output File Descriptions and Interpretation Guidelines

This document describes the files written by the current backend workflow in `r_equivalent.R`. The code writes these files into the selected output directory, typically `sim_output/` for simulated runs or `real_output/` for real datasets. The file set below reflects the current implementation, not the older generic `output/` examples. fileciteturn4file0turn4file1

## Output file list

The current backend writes the following outputs:

- `lsm_stage1.csv`
- `lsm_across.csv`
- `diffs_across.csv`
- `cld_across.csv`
- `entry_blups.csv`
- `stage1_model_selection.csv`
- `model_specs_stage1.csv`
- `diagnostics_stage1.csv`
- `cv_by_env.csv`
- `run_report.txt` fileciteturn4file0turn4file1

---

## 1. `lsm_stage1.csv`

**Purpose:** Contains Stage 1, within-environment least-squares means for each entry after the selected per-environment model is fit. Depending on the environment, this may come from a spatial `nlme::lme` model, an RCBD fallback via `lme4::lmer`, or an `lm` fallback when mixed modeling is not feasible. fileciteturn4file0turn4file1

**Typical columns:**

| Column | Description |
| --- | --- |
| `env` | Environment identifier, usually site-year. |
| `cov` | Selected Stage 1 model label for that environment, such as `expa`, `exp`, `sph`, `gau`, `none_rcbd`, or `none_lm`. |
| `entry` | Variety or genotype identifier. |
| `estimate` | Adjusted mean for that entry within that environment. |
| `stderr` | Standard error of the adjusted mean. |
| `df` | Degrees of freedom used by the mean comparison machinery. |
| `lower.CL`, `upper.CL` | Confidence limits for the adjusted mean. |
| `var_lsmean` | Variance of the LS-mean, calculated as `stderr^2`, used for Stage 2 inverse-variance weighting. |

**Interpretation:**
This is the core Stage 1 output and the foundation for Stage 2 summarization. Each row is an environment-specific adjusted mean, not a statewide or across-environment mean. Smaller `stderr` and `var_lsmean` values indicate more precise estimation for that entry in that environment. The `cov` column also tells you whether the environment used a spatial model or a fallback model. fileciteturn4file0turn4file5

**Important note:**
The code currently labels both `expa` and `exp` as separate candidate covariance types, but both are presently mapped to `nlme::corExp(...)` in the helper function. So `expa` should currently be interpreted as a placeholder label rather than a fully distinct anisotropic implementation. fileciteturn4file0

---

## 2. `lsm_across.csv`

**Purpose:** Contains the Stage 2 adjusted means for each entry. If multiple environments are available, these come from the across-environment inverse-variance meta-analysis. If only one environment is available, the backend skips meta-analysis and computes the table directly within that environment. fileciteturn4file0turn4file1

**Typical columns:**

| Column | Description |
| --- | --- |
| `entry` | Variety identifier. |
| `estimate` | Across-environment adjusted mean, or within-environment adjusted mean when only one environment exists. |
| `stderr` | Standard error of the adjusted mean. |
| `df` | Degrees of freedom used in the summary. |
| `lower.CL`, `upper.CL` | Confidence limits for the adjusted mean. |
| `LSD_0.05` | Single model-based least significant difference threshold at nominal α = 0.05. |
| `LSD_0.30` | Single model-based least significant difference threshold at nominal α = 0.30. |

**Interpretation:**
This is the main summary table for comparing overall entry performance. In multi-environment mode, entries with more precise Stage 1 LS-means contribute more weight because the model uses inverse-variance weighting. In single-environment mode, the values come from a direct within-environment model instead. fileciteturn4file1turn4file5

`LSD_0.05` and `LSD_0.30` are included for continuity with traditional reporting, but they are not Tukey-adjusted grouping thresholds. The frontend describes them as single, model-based least significant difference thresholds derived from the average standard error of unadjusted pairwise contrasts under an equal-precision approximation. They should be treated as rough reporting aids, while formal inference should come from `diffs_across.csv` and `cld_across.csv`. fileciteturn4file5turn4file7turn4file10

---

## 3. `diffs_across.csv`

**Purpose:** Contains pairwise entry comparisons from Stage 2 using Tukey adjustment whenever those contrasts can be computed. fileciteturn4file3turn4file5

**Typical columns:**

| Column | Description |
| --- | --- |
| `contrast` | Pairwise comparison, such as `A - B`. |
| `estimate` | Estimated difference in adjusted means. |
| `SE` | Standard error of the contrast. |
| `df` | Degrees of freedom used in the comparison, when available. |
| `t.ratio` | Test statistic for the contrast. |
| `p.value` | Tukey-adjusted p-value. |

**Interpretation:**
Use this table for formal pairwise inference. A small Tukey-adjusted `p.value` indicates evidence that two entries differ after accounting for multiple comparisons. This table is the detailed inferential companion to the more compact `cld_across.csv`. fileciteturn4file3turn4file5

---

## 4. `cld_across.csv`

**Purpose:** Contains the compact letter display derived from the Tukey-adjusted pairwise comparisons. fileciteturn4file3

**Typical columns:**

| Column | Description |
| --- | --- |
| `entry` | Variety identifier. |
| `group` | Compact letter assignment. |

**Interpretation:**
Entries that share at least one letter are not significantly different from one another at the chosen alpha level used for the Tukey comparisons, which in the current workflow is ordinarily α = 0.05. Entries with no letters in common are significantly different under that rule. This file is convenient for extension tables and summary figures, but it is less detailed than `diffs_across.csv`. fileciteturn4file7

---

## 5. `entry_blups.csv`

**Purpose:** Contains entry BLUP-style summaries from the Stage 2A one-stage mixed-effects model, or fallback entry means if the mixed model fails. fileciteturn4file0turn4file2

**Typical columns:**

| Column | Description |
| --- | --- |
| `entry` | Variety identifier. |
| `BLUP` | Entry random-effect estimate when mixed modeling succeeds, or fallback entry mean when the code falls back to `lm`. |
| `SE` | Conditional standard error for the BLUP, or standard error of the fallback entry summary. |

**Interpretation:**
When the mixed model succeeds, positive `BLUP` values indicate entries predicted above the population average random effect, and negative values indicate entries predicted below it. In single-environment or failed-model scenarios, the file may instead contain entry summaries from a reduced mixed model or an `lm` fallback, so the label `BLUP` should be interpreted with that implementation detail in mind. The code itself documents these fallback rules explicitly. fileciteturn4file2

---

## 6. `stage1_model_selection.csv`

**Purpose:** Provides the full Stage 1 model-selection record for each environment, including successful candidate fits and failed candidates. fileciteturn4file0turn4file1

**Typical columns:**

| Column | Description |
| --- | --- |
| `env` | Environment identifier. |
| `cov` | Candidate covariance or fallback label. |
| `engine` | Fitting engine, such as `nlme::lme`, `lme4::lmer`, or `stats::lm`. |
| `converged` | Whether that candidate fit succeeded. |
| `n` | Number of observations used by the fit. |
| `k` | Effective number of estimated parameters tracked for model-selection reporting. |
| `logLik` | Log-likelihood, when available. |
| `AIC` | Akaike information criterion. |
| `AICc` | Small-sample corrected AIC, when available. |
| `BIC` | Bayesian information criterion. |
| `best` | Logical flag indicating the selected model for that environment. |

**Interpretation:**
This is the audit trail for Stage 1 selection. It allows you to verify which spatial structure won in each environment, whether fallbacks were needed, and whether the selected model won by a wide or narrow AICc margin. This file is especially useful for debugging and for assessing whether spatial modeling materially improved fit. fileciteturn4file0turn4file11

---

## 7. `model_specs_stage1.csv`

**Purpose:** Stores a cleaner, human-readable record of the explicit Stage 1 model specification used for each environment. fileciteturn4file1

**Typical columns:**

| Column | Description |
| --- | --- |
| `env` | Environment identifier. |
| `stage` | Analysis stage, typically `stage1`. |
| `engine` | Fitting engine used. |
| `fixed` | Fixed-effects formula description. |
| `random` | Random-effects structure description. |
| `residual` | Residual or correlation structure description. |
| `method` | Estimation method, such as REML or OLS. |
| `cov_selected` | Selected covariance or fallback label. |

**Interpretation:**
Where `stage1_model_selection.csv` is mainly for selection statistics, `model_specs_stage1.csv` is better for documentation. It tells you, in plain model terms, what was actually fit for each environment, including whether the residual structure was spatial or whether the environment fell back to RCBD or `lm`. fileciteturn4file0

---

## 8. `diagnostics_stage1.csv`

**Purpose:** Contains plot-level fitted values and residual diagnostics from the selected Stage 1 model for each environment. fileciteturn4file0turn4file1

**Typical columns:**

| Column | Description |
| --- | --- |
| `env`, `site`, `year`, `rep`, `row`, `col`, `entry`, `yield` | Original identifying and observed-data fields carried into the diagnostics file. |
| `fitted` | Model-fitted value for that plot. |
| `resid` | Raw residual. |
| `resid_norm` | Normalized, standardized, or studentized residual depending on the fitted engine. |
| `resid_kind` | Text label describing what `resid_norm` means for that model. |
| `flag_outlier` | Logical flag for large residuals, currently based on absolute residual metric greater than 3. |

**Interpretation:**
This file supports quality control at the plot level. It is the source for residual diagnostics and residual heatmaps in the frontend. `flag_outlier` should be treated as a screening aid, not an automatic deletion rule. The residual standardization method depends on the model engine, so `resid_kind` matters when comparing environments. fileciteturn4file0turn4file9turn4file11

---

## 9. `cv_by_env.csv`

**Purpose:** Summarizes Stage 1 residual error by environment using RMSE and coefficient of variation. fileciteturn4file0turn4file1

**Typical columns:**

| Column | Description |
| --- | --- |
| `env` | Environment identifier. |
| `cov_selected` | Selected covariance or fallback label used in that environment. |
| `n_used` | Number of non-missing plots used in modeling. |
| `mean_y` | Mean observed yield for the environment. |
| `rmse` | Root mean squared residual error. |
| `cv_pct` | Percent CV, computed as `100 * rmse / mean_y` when possible. |

**Interpretation:**
This is a compact environment-level diagnostics table. Lower RMSE and CV generally indicate tighter fit relative to the scale of the response, though interpretation should always be contextualized by trial type, trait scale, and design quality. This file is useful for identifying unusually noisy environments and for comparing residual performance across selected covariance structures. fileciteturn4file0turn4file11

---

## 10. `run_report.txt`

**Purpose:** A plain-text audit trail summarizing what the analysis did, what data were used, how many environments were modeled, what Stage 1 and Stage 2 modes were used, and what assumptions/rules governed the run. fileciteturn4file1turn4file11

**Typical contents:**

- Timestamp
- Input file path
- Output directory
- Data summary, including missing yields and used rows
- Environment-level summary counts
- Stage 1 candidate covariance list
- Paths written for model-spec, selection, diagnostics, and CV files
- Count of environments using spatial versus fallback models
- Best covariance by environment
- Stage 2B mode description
- Stage 2A BLUP mode description
- Key modeling assumptions and fallback rules fileciteturn4file11

**Interpretation:**
This is the first file to inspect when you want a quick narrative summary of the run. It is especially valuable for reproducibility, troubleshooting, and reporting exactly what happened without needing to infer workflow details from multiple CSV outputs. fileciteturn4file11

---

## Practical reading order

For most users, the most efficient order is:

1. `run_report.txt` to understand what happened overall.
2. `lsm_across.csv` for the main adjusted means table.
3. `cld_across.csv` and `diffs_across.csv` for significance interpretation.
4. `entry_blups.csv` for one-stage entry ranking or breeding-oriented summaries.
5. `stage1_model_selection.csv`, `model_specs_stage1.csv`, `diagnostics_stage1.csv`, and `cv_by_env.csv` when you need technical validation or troubleshooting. fileciteturn4file11turn4file3

## Important cautions

- The current workflow uses α = 0.05 as the main inferential threshold in the code unless you intentionally change it. The frontend also states that Tukey CLD letters reflect α = 0.05. fileciteturn4file0turn4file7
- `LSD_0.30` is retained for continuity with older reporting conventions, but it should not be confused with the Tukey-based multiple-comparison results. fileciteturn4file5turn4file7
- Confidence intervals and p-values come from standard approximations used by `nlme`, `lme4`, and `emmeans`, and may be less stable when designs are highly unbalanced or degrees of freedom are small. fileciteturn4file11
- Because `entry_blups.csv` can fall back to `lm` summaries if mixed models fail, the file name is more stable than the exact statistical interpretation of its contents. Always cross-check the run report if that distinction matters. fileciteturn4file2turn4file11
