# Output File Descriptions and Interpretation Guidelines

### 1. `output/lsm_stage1.csv`

**Purpose:** Contains the within-environment least-squares means (LS-means) for each variety (entry) after accounting for spatial correlation.

**Columns:**

| Column                 | Description                                                       | Interpretation                                                    |
| ---------------------- | ----------------------------------------------------------------- | ----------------------------------------------------------------- |
| `env`                  | Environment identifier (site-year combination).                   | Each unique test environment.                                     |
| `cov`                  | Chosen spatial covariance model (`expa`, `exp`, `sph`, or `gau`). | The model with the lowest AICc was selected for that environment. |
| `entry`                | Variety or genotype identifier.                                   | Factor variable representing a tested entry.                      |
| `estimate`             | Adjusted mean yield (or trait value).                             | The LS-mean corrected for spatial trends and block effects.       |
| `stderr`               | Standard error of the LS-mean.                                    | Smaller SE indicates higher precision within that environment.    |
| `df`                   | Residual degrees of freedom.                                      | Used in confidence interval and test computation.                 |
| `lower.CL`, `upper.CL` | 95% confidence interval for the LS-mean.                          | Non-overlapping intervals suggest differences among entries.      |
| `var_lsmean`           | Variance used for inverse-variance weighting in Stage 2.          | Drives precision weighting across environments.                   |

**Interpretation:**
This file represents “cleaned” per-site means analogous to what SAS produces with spatial `sp(exp)` or similar models. These estimates account for spatial heterogeneity and random block structure, providing the most accurate representation of variety performance within each trial location.

---

### 2. `output/lsm_across.csv`

**Purpose:** Weighted meta-analysis of LS-means across environments.

**Columns:**

| Column                 | Description                                                                     |
| ---------------------- | ------------------------------------------------------------------------------- |
| `entry`                | Variety identifier.                                                             |
| `estimate`             | Overall mean performance across all environments, weighted by inverse variance. |
| `stderr`               | Standard error of the across-environment LS-mean.                               |
| `df`                   | Degrees of freedom from the meta-model.                                         |
| `lower.CL`, `upper.CL` | Confidence bounds for the across-environment mean.                              |

**Interpretation:**
Represents each variety’s average performance, adjusted for within-site uncertainty. Entries with smaller SEs are more stable across environments. Differences between entries can be formally tested using the Tukey-adjusted results in `diffs_across.csv`. This table is what you would report as the “adjusted mean yield” across locations in a statewide or multi-site variety report.

---

### 3. `output/diffs_across.csv`

**Purpose:** Pairwise Tukey-adjusted comparisons among entries across environments.

**Columns:**

| Column     | Description                                 |
| ---------- | ------------------------------------------- |
| `contrast` | Pair of compared entries (e.g., “E1 - E2”). |
| `estimate` | Difference in LS-means between entries.     |
| `SE`       | Standard error of the difference.           |
| `df`       | Degrees of freedom used in the comparison.  |
| `t.ratio`  | t-statistic for the difference.             |
| `p.value`  | Adjusted p-value (Tukey HSD correction).    |

**Interpretation:**
This file allows you to identify which entries differ significantly in mean performance. A small `p.value` (e.g., < 0.05) indicates a statistically significant difference between those two varieties’ adjusted means.

---

### 4. `output/cld_across.csv`

**Purpose:** Compact Letter Display summarizing the Tukey test results.

**Columns:**

| Column  | Description                                |
| ------- | ------------------------------------------ |
| `entry` | Variety identifier.                        |
| `group` | Letter group from multicomparison results. |

**Interpretation:**
Entries sharing at least one letter are **not significantly different** from each other at α = 0.05. Entries with no letters in common differ significantly.
Example:

| Entry                                   | LS-mean | Group |
| --------------------------------------- | ------- | ----- |
| A                                       | 5200    | a     |
| B                                       | 5050    | ab    |
| C                                       | 4900    | b     |
| Here, A = B (not different), but A ≠ C. |         |       |

---

### 5. `output/entry_blups.csv`

**Purpose:** Empirical Best Linear Unbiased Predictions (BLUPs) of entry effects from a one-stage mixed model combining all environments.

**Columns:**

| Column  | Description                                  |
| ------- | -------------------------------------------- |
| `entry` | Variety identifier.                          |
| `BLUP`  | Random-effect deviation from the grand mean. |
| `SE`    | Conditional standard error of the BLUP.      |

**Interpretation:**
BLUPs reflect predicted performance of each entry after accounting for all random sources of variation (environments, blocks, and G×E).

* **Positive BLUP** → Entry performed above the population mean.
* **Negative BLUP** → Entry performed below the population mean.
  You can add the overall intercept (from `fixef(fit2)["(Intercept)"]`) to each BLUP to obtain predicted overall means.
  Because BLUPs shrink extreme values toward the mean, they provide more reliable estimates of true genetic performance, particularly for unbalanced or partially replicated trials.

---

### 6. Optional diagnostic plots (if added later)

* **Residual heatmaps** per environment to visualize spatial correction effectiveness.
  Patterns should be random (no large patches of residual structure).
* **Across-environment mean plots** with CLD letters, showing adjusted means and confidence intervals.

---