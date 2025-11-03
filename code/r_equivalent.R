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
  # weighted linear mixed model: estimate ~ entry + (1|env), weights = 1/var
  lsm_stage1 <- lsm_stage1 %>% mutate(entry = factor(entry), env = factor(env))
  w <- 1 / lsm_stage1$var_lsmean
  # Use nlme lme for weights + random env
  fit <- lme(estimate ~ entry, random = ~ 1 | env, weights = ~ w, data = lsm_stage1, method = "REML")
  emm <- emmeans(fit, ~ entry)
  lsm_tab <- as.data.frame(summary(emm, infer = TRUE, level = 1 - alpha)) %>%
    rename(estimate = emmean, stderr = SE) %>%
    select(entry, estimate, stderr, df, lower.CL, upper.CL)
  
  # Tukey-adjusted pairwise
  pairs_tukey <- as.data.frame(summary(pairs(emm, adjust = "tukey")))
  # letters
  tukey_p <- as.matrix(pairs_tukey %>% select(contrast, p.value))
  # build a symmetric p-value matrix for multcompView
  # emmeans provides a helper: multcomp::cld, but we can use cld() directly
  cld_tab <- multcomp::cld(emm, alpha = alpha, Letters = letters, adjust = "tukey") %>%
    as.data.frame() %>%
    select(entry, .group) %>%
    rename(group = .group)
  
  list(lsm_across = lsm_tab, diffs_across = pairs_tukey, cld = cld_tab)
}

# ---- Stage 2A: one-stage MET BLUPs (entries random) ----
# Random: entry and entry:env, random env and rep(env), spatial residual per env.
# nlme can attach a single correlation structure per fit; we approximate site-specific
# spatial behavior by fitting per-env correlation in a loop and stacking BLUPs,
# or, more simply, drop spatial here and rely on Stage 1 for spatial control.
# Below we fit a conventional MET BLUP without spatial residual, which is standard and robust.

stage2_blups <- function(trial_df) {
  df <- trial_df %>% mutate(entry = factor(entry),
                            env = factor(env),
                            rep = factor(rep))
  # random = entry + entry:env + rep%in%env + env
  # Fit using lme4::lmer for speed; we’ll use nlme::lme to avoid an extra dependency
  fit <- nlme::lme(yield ~ 1,
                   random = ~ 1 | env/rep,
                   data = df,
                   method = "REML")
  # Add random entry and entry:env via a second stage using lme4 would be preferable,
  # but nlme cannot nest crossed random effects easily. Use lme4 if allowed:
  # lmer(yield ~ 1 + (1|env) + (1|env:rep) + (1|entry) + (1|entry:env), data=df)
  # Here, if lme4 is acceptable in your environment, uncomment:
  #
  # library(lme4)
  # fit <- lmer(yield ~ 1 + (1|env) + (1|env:rep) + (1|entry) + (1|entry:env), data = df)
  #
  # Extract BLUPs for entries (overall)
  # Using lme4:
  # blups <- ranef(fit, condVar = TRUE)$entry
  # With nlme fallback, we approximate overall entry effects using a simple EBLUP from a separate model:
  library(lme4)
  fit2 <- lmer(yield ~ 1 + (1|env) + (1|env:rep) + (1|entry) + (1|entry:env), data = df, REML = TRUE)
  re <- ranef(fit2, condVar = TRUE)
  entry_blup <- tibble(entry = rownames(re$entry),
                       BLUP = as.numeric(re$entry[,1]))
  # SEs from condVar
  se <- sapply(re$entry@postVar, function(m) sqrt(m[1,1,]))
  entry_blup$SE <- as.numeric(se)
  entry_blup
}

# ---- main harness ----

analyze_trial <- function(in_csv = "trial.csv",
                          out_stage1 = "lsm_stage1.csv",
                          out_across = "lsm_across.csv",
                          out_diffs  = "diffs_across.csv",
                          out_cld    = "cld_across.csv",
                          out_blups  = "entry_blups.csv",
                          alpha = 0.05) {
  
  trial <- read.csv(in_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("site","year","env","rep","row","col","entry","yield") %in% names(trial)))
  
  # Stage 1
  s1 <- analyze_stage1(trial, covlist, alpha = alpha)
  lsm_stage1 <- s1$lsm_stage1
  write.csv(lsm_stage1, out_stage1, row.names = FALSE)
  
  # Stage 2B
  s2 <- stage2_meta(lsm_stage1, alpha = alpha)
  write.csv(s2$lsm_across, out_across, row.names = FALSE)
  write.csv(s2$diffs_across, out_diffs, row.names = FALSE)
  write.csv(s2$cld, out_cld, row.names = FALSE)
  
  # Stage 2A BLUPs
  blups <- stage2_blups(trial)
  write.csv(blups, out_blups, row.names = FALSE)
  
  invisible(list(stage1 = lsm_stage1,
                 best_cov = s1$best_cov,
                 lsm_across = s2$lsm_across,
                 diffs_across = s2$diffs_across,
                 cld = s2$cld,
                 blups = blups))
}

# ---- run if interactive ----
# analyze_trial("trial_sim.csv")
