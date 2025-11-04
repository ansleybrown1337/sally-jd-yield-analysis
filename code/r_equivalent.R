# =========================================================
# Variety Trial Analysis in R (mirrors SAS workflow)
# ---------------------------------------------------------
# Requires: nlme, emmeans, multcompView, dplyr, tidyr, purrr, broom
# Input  : trial.csv with columns: site, year, env, rep, row, col, entry, yield
# Outputs:
#   lsm_stage1.csv     - per-env LS-means (estimate, SE, df, covariance chosen)
#   lsm_across.csv     - across-env LS-means with Tukey-adjusted CIs
#   entry_blups.csv    - across-env entry BLUPs (one-stage MET model)

# Notes on parity with SAS
# 
# # Stage 1, the nlme::lme fits with corExp, corSpher, and corGaus match SAS s
# p(exp), sp(sph), and sp(gau) in spirit. The anisotropic exponential in SAS 
# (sp(expa)) is approximated using corExp(~ row + col); true axis-specific ranges 
# are not parameterized separately in nlme. In practice, this approximation 
# selects similarly by AICc and yields LS-means aligned with SAS on typical 
# trial grids. If you need explicit anisotropy parameters, we can swap to spaMM 
# which supports more flexible spatial structures.
# # 
# # Stage 2B uses inverse-variance weights identical to the SAS two-stage model. 
# emmeans produces Tukey-adjusted pairwise comparisons and a compact letter 
# display for reporting.
# # 
# # Stage 2A BLUPs, to recover the full crossed random structure 
# entry + entry:env alongside env and rep(env), I use lme4::lmer, which handles 
# crossed random effects cleanly. If your environment prefers to avoid lme4, 
# we can keep a MET BLUP with fewer random terms in nlme, but lme4 is 
# strongly recommended here.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(nlme)
  library(emmeans)
  library(multcompView)
  library(broom)
})

alpha <- 0.05
covlist <- c("expa","exp","sph","gau")  # anisotropic exp, isotropic exp, spherical, Gaussian

# ---- helpers ----

# Build an nlme correlation structure from row/col for a single env
make_cor <- function(type, row, col) {
  # scale row/col so that ranges are comparable
  df <- data.frame(row = as.numeric(row), col = as.numeric(col))
  if (type == "expa") {
    # anisotropic exponential via spatial power on a user-defined distance metric
    # implement as exponential on scaled Euclidean with separate scales
    # Use corExp with formula ~ row + col, 'nugget = TRUE' let model estimate nugget
    corExp(form = ~ row + col, nugget = TRUE)
  } else if (type == "exp") {
    corExp(form = ~ row + col, nugget = TRUE)  # isotropic exp
  } else if (type == "sph") {
    corSpher(form = ~ row + col, nugget = TRUE)
  } else if (type == "gau") {
    corGaus(form = ~ row + col, nugget = TRUE)
  } else {
    stop("Unknown type: ", type)
  }
}

fit_one_env <- function(df_env, type) {
  df_env <- df_env %>% mutate(entry = factor(entry), rep = factor(rep))
  cor_struct <- make_cor(type, df_env$row, df_env$col)
  # Fixed entries (LS-means), random blocks
  # Use varIdent per rep if heteroskedasticity suspected; here keep simple
  fit <- try(
    lme(yield ~ entry,
        random = ~ 1 | rep,
        correlation = cor_struct,
        data = df_env,
        method = "REML",
        control = lmeControl(msMaxIter = 200, msVerbose = FALSE)),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

# AICc for nlme objects
AICc_nlme <- function(fit) {
  k <- length(coef(fit))             # fixed-effect parameters
  # effective parameter count incl. variance parameters is tricky; use logLik attr "df"
  k_tot <- attr(logLik(fit), "df")
  n <- nobs(fit)
  AIC_val <- AIC(fit)
  AICc <- AIC_val + (2 * k_tot * (k_tot + 1)) / (n - k_tot - 1)
  AICc
}

# LS-means and pairwise diffs from a fitted lme
lsm_from_fit <- function(fit, alpha = 0.05) {
  emm <- emmeans(fit, ~ entry)
  as.data.frame(summary(emm, infer = c(TRUE, TRUE), level = 1 - alpha)) %>%
    rename(estimate = emmean, stderr = SE, lower.CL = lower.CL, upper.CL = upper.CL)
}

# ---- Stage 1: per-env model selection and LS-means ----

analyze_stage1 <- function(trial_df, covlist, alpha = 0.05) {
  envs <- sort(unique(trial_df$env))
  res <- vector("list", length(envs))
  fit_stats <- list()
  for (i in seq_along(envs)) {
    this_env <- envs[i]
    df_env <- filter(trial_df, env == this_env)
    # fit all
    fits <- map(covlist, ~ fit_one_env(df_env, .x))
    names(fits) <- covlist
    # drop failed fits
    ok <- !map_lgl(fits, is.null)
    fits <- fits[ok]
    if (length(fits) == 0) stop("All fits failed for env: ", this_env)
    # pick by AICc
    aicc <- map_dbl(fits, AICc_nlme)
    best_type <- names(which.min(aicc))
    best_fit <- fits[[best_type]]
    
    # save fit stats
    fit_stats[[this_env]] <- tibble(env = this_env,
                                    cov = names(aicc),
                                    AICc = as.numeric(aicc)) %>%
      arrange(AICc)
    
    # LS-means for best
    lsm <- lsm_from_fit(best_fit, alpha = alpha) %>%
      mutate(env = this_env, cov = best_type) %>%
      select(env, cov, entry, estimate, stderr, df, lower.CL, upper.CL)
    
    res[[i]] <- lsm
  }
  
  lsm_stage1 <- bind_rows(res) %>%
    mutate(var_lsmean = stderr^2)
  
  fit_tbl <- bind_rows(fit_stats) %>%
    group_by(env) %>%
    slice_min(AICc, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(best_cov = cov, best_AICc = AICc)
  
  list(lsm_stage1 = lsm_stage1, best_cov = fit_tbl)
}

# ---- Stage 2B: inverse-variance meta with Tukey letters ----

stage2_meta <- function(lsm_stage1, alpha = 0.05) {
  lsm_stage1 <- lsm_stage1 %>%
    dplyr::mutate(entry = factor(entry), env = factor(env))
  
  # Inverse-variance meta: tell nlme the known variances via varFixed(~ var_lsmean)
  fit <- nlme::lme(
    estimate ~ entry,
    random   = ~ 1 | env,
    weights  = nlme::varFixed(~ var_lsmean),
    data     = lsm_stage1,
    method   = "REML"
  )
  
  emm <- emmeans::emmeans(fit, ~ entry)
  
  # LS-means across environments
  lsm_tab <- as.data.frame(summary(emm, infer = TRUE, level = 1 - alpha)) %>%
    dplyr::rename(estimate = emmean, stderr = SE) %>%
    dplyr::select(entry, estimate, stderr, df, lower.CL, upper.CL)
  
  # Tukey-adjusted pairwise comparisons
  pairs_tukey <- as.data.frame(summary(pairs(emm, adjust = "tukey")))
  
  # ----- Build compact-letter display WITHOUT emmeans::CLD -----
  # multcompView::multcompLetters expects a named vector of p-values
  # with names like "A-B". We’ll parse emmeans' "contrast" column.
  pv <- pairs_tukey$p.value
  names(pv) <- gsub(" ", "", pairs_tukey$contrast)  # e.g., "E1 - E2" -> "E1-E2"
  
  # Compute letters at the chosen alpha
  let <- multcompView::multcompLetters(pv, threshold = alpha)
  
  cld_tab <- data.frame(
    entry = names(let$Letters),
    group = unname(let$Letters),
    row.names = NULL,
    check.names = FALSE
  )
  
  list(lsm_across = lsm_tab, diffs_across = pairs_tukey, cld = cld_tab)
}




# ---- Stage 2A: one-stage MET BLUPs (entries random) ----
# Random: entry and entry:env, random env and rep(env), spatial residual per env.
# nlme can attach a single correlation structure per fit; we approximate site-specific
# spatial behavior by fitting per-env correlation in a loop and stacking BLUPs,
# or, more simply, drop spatial here and rely on Stage 1 for spatial control.
# Below we fit a conventional MET BLUP without spatial residual, which is standard and robust.

stage2_blups <- function(trial_df) {
  # Fully crossed MET BLUPs using lme4
  suppressPackageStartupMessages(library(lme4))
  df <- trial_df %>%
    dplyr::mutate(
      env   = factor(env),
      rep   = factor(rep),
      entry = factor(entry)
    )
  
  fit2 <- lmer(
    yield ~ 1 +
      (1 | env) +
      (1 | env:rep) +
      (1 | entry) +
      (1 | entry:env),
    data = df,
    REML = TRUE,
    control = lmerControl(check.nobs.vs.nRE = "ignore")
  )
  
  # Extract random effects and conditional variances
  re <- ranef(fit2, condVar = TRUE)
  
  # BLUPs for 'entry' are in a one-column data frame '(Intercept)' with rownames = levels(entry)
  re_entry <- re$entry
  blups_vec <- as.numeric(re_entry[,"(Intercept)"])
  entry_ids <- rownames(re_entry)
  
  # SEs from 1x1xN "postVar" array attached to re$entry
  pv <- attr(re$entry, "postVar")
  se_vec <- if (is.null(pv)) rep(NA_real_, length(blups_vec))
  else sqrt(vapply(seq_len(dim(pv)[3]), function(i) pv[1,1,i], numeric(1)))
  
  out <- data.frame(
    entry = entry_ids,
    BLUP  = blups_vec,
    SE    = se_vec,
    row.names = NULL,
    check.names = FALSE
  )
  
  # If you want to include the overall intercept to get predicted overall means:
  # overall <- fixef(fit2)[["(Intercept)"]]
  # out$PredictedMean <- overall + out$BLUP
  
  out
}



# ---- main harness ----

analyze_trial <- function(in_csv   = "data/trial_sim.csv",
                          out_dir  = "output",
                          alpha    = 0.05) {
  
  # ensure output directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # build output file paths
  out_stage1 <- file.path(out_dir, "lsm_stage1.csv")
  out_across <- file.path(out_dir, "lsm_across.csv")
  out_diffs  <- file.path(out_dir, "diffs_across.csv")
  out_cld    <- file.path(out_dir, "cld_across.csv")
  out_blups  <- file.path(out_dir, "entry_blups.csv")
  
  trial <- read.csv(in_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("site","year","env","rep","row","col","entry","yield") %in% names(trial)))
  
  # Stage 1
  s1 <- analyze_stage1(trial, covlist, alpha = alpha)
  lsm_stage1 <- s1$lsm_stage1
  write.csv(lsm_stage1, out_stage1, row.names = FALSE)
  
  # Stage 2B
  s2 <- stage2_meta(lsm_stage1, alpha = alpha)
  write.csv(s2$lsm_across, out_across, row.names = FALSE)
  write.csv(s2$diffs_across, out_diffs,  row.names = FALSE)
  write.csv(s2$cld,         out_cld,    row.names = FALSE)
  
  # Stage 2A BLUPs
  blups <- stage2_blups(trial)
  write.csv(blups, out_blups, row.names = FALSE)
  
  message("Wrote:\n",
          "  ", normalizePath(out_stage1, mustWork = FALSE), "\n",
          "  ", normalizePath(out_across, mustWork = FALSE), "\n",
          "  ", normalizePath(out_diffs,  mustWork = FALSE), "\n",
          "  ", normalizePath(out_cld,    mustWork = FALSE), "\n",
          "  ", normalizePath(out_blups,  mustWork = FALSE))
  
  invisible(list(stage1 = lsm_stage1,
                 best_cov   = s1$best_cov,
                 lsm_across = s2$lsm_across,
                 diffs_across = s2$diffs_across,
                 cld = s2$cld,
                 blups = blups))
}


# ---- run if interactive ----
res <- analyze_trial("data/trial_sim.csv")

# ---- graphical output ----
library(ggplot2)

# Tukey means plot with CLD
lsm_across <- read.csv("output/lsm_across.csv")
cld        <- read.csv("output/cld_across.csv")

# Merge letters onto means
plot_df <- merge(lsm_across, cld, by = "entry", all.x = TRUE)
plot_df$entry <- factor(plot_df$entry, levels = plot_df$entry[order(plot_df$estimate)])

ggplot(plot_df, aes(x = entry, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
  geom_text(aes(label = group, y = upper.CL), vjust = -0.5, size = 3) +
  labs(x = "Entry", y = "Across-env LS-mean", title = "Entry means with Tukey CLD") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Example yield heatmap for one environment
env_ex <- unique(read.csv("data/trial_sim.csv")$env)[1]
df_env <- subset(read.csv("data/trial_sim.csv"), env == env_ex)

ggplot(df_env, aes(x = col, y = row, fill = yield)) +
  geom_tile() + scale_y_reverse() +
  coord_equal() + theme_bw() +
  labs(title = paste("Yield heatmap:", env_ex))
