# =========================================================
# Variety Trial Analysis in R (mirrors SAS workflow)
# ---------------------------------------------------------
# Requires: nlme, lme4, emmeans, multcompView, dplyr, tidyr, purrr, broom, ggplot2
# Input  : trial.csv with columns: site, year, env, rep, row, col, entry, yield
# Outputs:
#   lsm_stage1.csv     - per-env LS-means (estimate, SE, df, covariance chosen)
#   lsm_across.csv     - across-env LS-means with Tukey-adjusted CIs
#   entry_blups.csv    - across-env entry BLUPs (one-stage MET model)
#
# Stage 1:
#   • Try nlme spatial models (corExp, corSpher, corGaus) per environment.
#   • If ALL spatial fits fail, fall back to:
#       - RCBD via lme4::lmer(yield ~ entry + (1 | rep)) when ≥ 2 usable reps.
#       - Fixed-effects ANOVA via lm(yield ~ entry) when < 2 usable reps.
#   • In fallback cases, we explicitly report that only block-adjusted (or
#     non-block-adjusted) entry means are used, with no spatial correction.
#
# Stage 2B:
#   • Inverse-variance meta analysis across environments (nlme::lme with varFixed).
#   • Tukey-adjusted pairwise comparisons and CLD using emmeans + multcompView.
#
# Stage 2A:
#   • MET BLUPs via lme4::lmer with env, env:rep, entry, and entry:env as random.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(nlme)
  library(lme4)
  library(emmeans)
  library(ggplot2)
  library(multcompView)
  library(broom)
})

alpha   <- 0.05
covlist <- c("expa","exp","sph","gau")  # anisotropic exp, isotropic exp, spherical, Gaussian

# =========================================================
# Helper: correlation structure for nlme spatial models
# =========================================================

make_cor <- function(type, row, col) {
  df <- data.frame(row = as.numeric(row), col = as.numeric(col))
  if (type == "expa") {
    corExp(form = ~ row + col, nugget = TRUE)
  } else if (type == "exp") {
    corExp(form = ~ row + col, nugget = TRUE)
  } else if (type == "sph") {
    corSpher(form = ~ row + col, nugget = TRUE)
  } else if (type == "gau") {
    corGaus(form = ~ row + col, nugget = TRUE)
  } else {
    stop("Unknown covariance type: ", type)
  }
}

# =========================================================
# Helper: AICc for nlme fits
# =========================================================

AICc_nlme <- function(fit) {
  k_tot <- attr(logLik(fit), "df")
  n     <- nobs(fit)
  AIC_v <- AIC(fit)
  AIC_v + (2 * k_tot * (k_tot + 1)) / (n - k_tot - 1)
}

# =========================================================
# Helper: LS-means extractor (works for lme and lmer, lm)
# =========================================================

lsm_from_fit <- function(fit, alpha = 0.05) {
  emm <- emmeans(fit, ~ entry)
  as.data.frame(summary(emm, infer = c(TRUE, TRUE), level = 1 - alpha)) %>%
    dplyr::rename(estimate = emmean, stderr = SE) %>%
    dplyr::select(entry, estimate, stderr, df, lower.CL, upper.CL)
}

# =========================================================
# Stage 1 helpers: nlme spatial fits per environment
# =========================================================

fit_one_env_spatial <- function(df_env, type) {
  # Drop missing yields so nlme only sees usable plots
  df_env <- df_env %>%
    dplyr::filter(!is.na(yield)) %>%
    dplyr::mutate(
      entry = factor(entry),
      rep   = factor(rep)
    )
  
  if (nrow(df_env) == 0L) return(NULL)
  
  cor_struct <- make_cor(type, df_env$row, df_env$col)
  
  fit <- try(
    nlme::lme(
      yield ~ entry,
      random      = ~ 1 | rep,
      correlation = cor_struct,
      data        = df_env,
      method      = "REML",
      control     = nlme::lmeControl(msMaxIter = 200, msVerbose = FALSE)
    ),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) return(NULL)
  fit
}

# =========================================================
# Stage 1: per-env model selection with robust fallbacks
# =========================================================

analyze_stage1 <- function(trial_df, covlist, alpha = 0.05) {
  envs      <- sort(unique(trial_df$env))
  res       <- vector("list", length(envs))
  fit_stats <- vector("list", length(envs))
  
  for (i in seq_along(envs)) {
    this_env <- envs[i]
    
    df_env <- trial_df %>%
      dplyr::filter(env == this_env) %>%
      dplyr::filter(!is.na(yield))
    
    if (nrow(df_env) == 0L) {
      warning("No non-missing yields for env ", this_env,
              "; skipping Stage 1 for this env.")
      next
    }
    
    # -----------------------------------------------------
    # 1) Try spatial nlme models for this environment
    # -----------------------------------------------------
    fits <- purrr::map(covlist, ~ fit_one_env_spatial(df_env, .x))
    names(fits) <- covlist
    
    ok_idx  <- !purrr::map_lgl(fits, is.null)
    fits_ok <- fits[ok_idx]
    
    if (length(fits_ok) > 0L) {
      # At least one spatial model converged -> pick by AICc
      aicc      <- purrr::map_dbl(fits_ok, AICc_nlme)
      best_type <- names(which.min(aicc))
      best_fit  <- fits_ok[[best_type]]
      
      # LS-means for best spatial model
      lsm <- lsm_from_fit(best_fit, alpha = alpha) %>%
        dplyr::mutate(env = this_env, cov = best_type) %>%
        dplyr::select(env, cov, entry, estimate, stderr, df, lower.CL, upper.CL)
      
      res[[i]] <- lsm
      
      # Store only the best covariance choice per env
      fit_stats[[i]] <- tibble::tibble(
        env       = this_env,
        best_cov  = best_type,
        best_AICc = min(aicc)
      )
      
      next
    }
    
    # -----------------------------------------------------
    # 2) All spatial fits failed: fall back to RCBD or lm
    # -----------------------------------------------------
    warning(
      "All spatial fits failed for env ", this_env,
      ". We could not reliably estimate spatial covariance for this trial, ",
      "so we report entry means without spatial correction."
    )
    
    df_env <- df_env %>%
      dplyr::mutate(
        entry = factor(entry),
        rep   = factor(rep)
      )
    
    n_reps <- dplyr::n_distinct(df_env$rep)
    
    # Decide between RCBD (random rep) and one-way ANOVA (no random rep)
    if (n_reps >= 2L) {
      # Attempt RCBD via lme4::lmer
      rcbd_fit <- try(
        lme4::lmer(
          yield ~ entry + (1 | rep),
          data    = df_env,
          REML    = TRUE,
          control = lme4::lmerControl(check.nobs.vs.nRE = "ignore")
        ),
        silent = TRUE
      )
      
      if (!inherits(rcbd_fit, "try-error")) {
        # Successful RCBD fallback
        lsm <- lsm_from_fit(rcbd_fit, alpha = alpha) %>%
          dplyr::mutate(env = this_env, cov = "none_rcbd") %>%
          dplyr::select(env, cov, entry, estimate, stderr, df, lower.CL, upper.CL)
        
        res[[i]] <- lsm
        
        fit_stats[[i]] <- tibble::tibble(
          env       = this_env,
          best_cov  = "none_rcbd",
          best_AICc = NA_real_
        )
        
        next
      } else {
        warning(
          "RCBD fallback (lmer) also failed for env ", this_env,
          ". Falling back further to fixed-effects ANOVA (lm, no random rep)."
        )
      }
    } else {
      warning(
        "Environment ", this_env,
        " has fewer than 2 usable reps; RCBD cannot be fit. ",
        "Using fixed-effects ANOVA (lm) with entry as the only factor."
      )
    }
    
    # Final fallback: fixed-effects ANOVA via lm(yield ~ entry)
    lm_fit <- lm(yield ~ entry, data = df_env)
    
    lsm <- lsm_from_fit(lm_fit, alpha = alpha) %>%
      dplyr::mutate(env = this_env, cov = "none_lm") %>%
      dplyr::select(env, cov, entry, estimate, stderr, df, lower.CL, upper.CL)
    
    res[[i]] <- lsm
    
    fit_stats[[i]] <- tibble::tibble(
      env       = this_env,
      best_cov  = "none_lm",
      best_AICc = NA_real_
    )
  }
  
  # Combine all environments
  lsm_stage1 <- dplyr::bind_rows(res) %>%
    dplyr::mutate(var_lsmean = stderr^2)
  
  best_cov_tbl <- dplyr::bind_rows(fit_stats)
  
  list(
    lsm_stage1 = lsm_stage1,
    best_cov   = best_cov_tbl
  )
}

# =========================================================
# Stage 2B: inverse-variance meta with Tukey letters
# =========================================================

stage2_meta <- function(lsm_stage1, alpha = 0.05) {
  # helper to normalize entry labels like "entry12" -> "12"
  clean_entry_labels <- function(x) {
    sub("^entry", "", as.character(x))
  }
  
  # Ensure factors
  lsm_stage1 <- lsm_stage1 %>%
    dplyr::mutate(
      entry = factor(entry),
      env   = factor(env)
    )
  
  # Drop rows with NA or nonpositive variances; nlme::varFixed requires > 0
  lsm2 <- lsm_stage1 %>%
    dplyr::filter(!is.na(var_lsmean), var_lsmean > 0)
  
  if (nrow(lsm2) == 0L) {
    stop(
      "Stage 2 meta-analysis: no rows with positive var_lsmean. ",
      "Check that Stage 1 produced finite standard errors."
    )
  }
  
  n_env <- dplyr::n_distinct(lsm2$env)
  n_ent <- dplyr::n_distinct(lsm2$entry)
  
  # --------------------------------------------------------
  # Very small case: fewer than 2 entries with usable LS-means
  # --------------------------------------------------------
  if (n_ent < 2L) {
    warning(
      "Stage 2 meta-analysis: fewer than 2 entries with usable LS-means. ",
      "Reporting marginal means without pairwise comparisons or CLD."
    )
    
    fit <- lm(estimate ~ 1, data = lsm2, weights = 1 / var_lsmean)
    emm <- emmeans::emmeans(fit, ~ 1)
    base_mean <- as.data.frame(summary(emm, infer = TRUE, level = 1 - alpha))
    
    # replicate same mean for each entry
    lsm_tab <- data.frame(
      entry     = clean_entry_labels(levels(lsm2$entry)),
      estimate  = base_mean$emmean,
      stderr    = base_mean$SE,
      df        = base_mean$df,
      lower.CL  = base_mean$lower.CL,
      upper.CL  = base_mean$upper.CL,
      stringsAsFactors = FALSE
    )
    
    cld_tab <- data.frame(
      entry = lsm_tab$entry,
      group = "a",
      row.names   = NULL,
      check.names = FALSE
    )
    
    return(list(
      lsm_across   = lsm_tab,
      diffs_across = data.frame(),
      cld          = cld_tab
    ))
  }
  
  # --------------------------------------------------------
  # Fit model for across-env means
  # --------------------------------------------------------
  if (n_env >= 2L) {
    # Preferred case: random env effect with known inverse-variance weights
    fit <- nlme::lme(
      estimate ~ entry,
      random   = ~ 1 | env,
      weights  = nlme::varFixed(~ var_lsmean),
      data     = lsm2,
      method   = "REML"
    )
  } else {
    # Fallback: fewer than 2 environments with usable weights
    warning(
      "Stage 2 meta-analysis: fewer than 2 environments with positive var_lsmean. ",
      "Fitting a weighted fixed-effects model (lm) without random env."
    )
    
    fit <- lm(
      estimate ~ entry,
      data    = lsm2,
      weights = 1 / var_lsmean
    )
  }
  
  # emmeans works for both lme and lm
  emm <- emmeans::emmeans(fit, ~ entry)
  
  # LS-means across environments (or weighted means in the fallback)
  lsm_tab <- as.data.frame(summary(emm, infer = TRUE, level = 1 - alpha)) %>%
    dplyr::rename(estimate = emmean, stderr = SE) %>%
    dplyr::select(entry, estimate, stderr, df, lower.CL, upper.CL) %>%
    dplyr::mutate(entry = clean_entry_labels(entry))
  
  # --------------------------------------------------------
  # Tukey-adjusted pairwise comparisons (if possible)
  # --------------------------------------------------------
  pairs_tukey <- tryCatch(
    as.data.frame(summary(pairs(emm, adjust = "tukey"))),
    error = function(e) NULL
  )
  
  # If pairwise contrasts failed or there are no contrasts, give everyone "a"
  if (is.null(pairs_tukey) || nrow(pairs_tukey) == 0L) {
    warning(
      "Stage 2 meta-analysis: no valid pairwise contrasts for CLD. ",
      "Assigning all entries to a single group 'a'."
    )
    
    cld_tab <- data.frame(
      entry = lsm_tab$entry,
      group = "a",
      row.names   = NULL,
      check.names = FALSE
    )
    
    return(list(
      lsm_across   = lsm_tab,
      diffs_across = if (is.null(pairs_tukey)) data.frame() else pairs_tukey,
      cld          = cld_tab
    ))
  }
  
  # Remove NA p-values before calling multcompLetters
  pv       <- pairs_tukey$p.value
  name_vec <- gsub(" ", "", pairs_tukey$contrast)  # "entry1 - entry2" -> "entry1-entry2"
  
  valid    <- !is.na(pv)
  pv       <- pv[valid]
  name_vec <- name_vec[valid]
  
  if (length(pv) == 0L) {
    warning(
      "Stage 2 meta-analysis: all pairwise p-values are NA. ",
      "Assigning all entries to a single group 'a'."
    )
    
    cld_tab <- data.frame(
      entry = lsm_tab$entry,
      group = "a",
      row.names   = NULL,
      check.names = FALSE
    )
    
    return(list(
      lsm_across   = lsm_tab,
      diffs_across = pairs_tukey,
      cld          = cld_tab
    ))
  }
  
  names(pv) <- name_vec
  
  # multcompLetters can still occasionally choke; protect with tryCatch
  let <- tryCatch(
    multcompView::multcompLetters(pv, threshold = alpha),
    error = function(e) NULL
  )
  
  if (is.null(let)) {
    warning(
      "Stage 2 meta-analysis: multcompLetters failed. ",
      "Assigning all entries to a single group 'a'."
    )
    
    cld_tab <- data.frame(
      entry = lsm_tab$entry,
      group = "a",
      row.names   = NULL,
      check.names = FALSE
    )
  } else {
    # Clean off any "entry" prefix here
    entry_names_clean <- clean_entry_labels(names(let$Letters))
    cld_tab <- data.frame(
      entry = entry_names_clean,
      group = unname(let$Letters),
      row.names   = NULL,
      check.names = FALSE
    )
  }
  
  list(
    lsm_across   = lsm_tab,
    diffs_across = pairs_tukey,
    cld          = cld_tab
  )
}

# =========================================================
# Stage 2B (single environment): within-env LS-means + Tukey CLD
# =========================================================

stage2_single_env <- function(trial_df, env_id, alpha = 0.05) {
  # Work only on the chosen environment, drop missing yields
  df_env <- trial_df %>%
    dplyr::filter(env == env_id, !is.na(yield)) %>%
    dplyr::mutate(
      entry = factor(entry),
      rep   = factor(rep)
    )
  
  n_rep <- dplyr::n_distinct(df_env$rep)
  
  # Prefer a mixed model with random reps if possible
  if (n_rep >= 2L) {
    fit <- try(
      lme4::lmer(
        yield ~ entry + (1 | rep),
        data    = df_env,
        REML    = TRUE,
        control = lme4::lmerControl(check.nobs.vs.nRE = "ignore")
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      warning(
        "Stage 2 single-env: lmer(yield ~ entry + (1 | rep)) failed for env ",
        env_id, "; falling back to lm(yield ~ entry)."
      )
      fit <- lm(yield ~ entry, data = df_env)
    }
  } else {
    warning(
      "Stage 2 single-env: env ", env_id,
      " has fewer than 2 reps; using fixed-effects ANOVA lm(yield ~ entry)."
    )
    fit <- lm(yield ~ entry, data = df_env)
  }
  
  # LS-means and CIs within this environment
  emm <- emmeans::emmeans(fit, ~ entry)
  lsm_tab <- as.data.frame(summary(emm, infer = TRUE, level = 1 - alpha)) %>%
    dplyr::rename(estimate = emmean, stderr = SE) %>%
    dplyr::select(entry, estimate, stderr, df, lower.CL, upper.CL)
  
  # Tukey-adjusted pairwise comparisons
  pairs_tukey <- as.data.frame(summary(pairs(emm, adjust = "tukey")))
  
  # Compact letter display (CLD)
  pv <- pairs_tukey$p.value
  names(pv) <- gsub(" ", "", pairs_tukey$contrast)  # "E1 - E2" -> "E1-E2"
  
  let <- multcompView::multcompLetters(pv, threshold = alpha)
  
  cld_tab <- data.frame(
    entry = names(let$Letters),
    group = unname(let$Letters),
    row.names   = NULL,
    check.names = FALSE
  )
  
  list(
    lsm_across   = lsm_tab,
    diffs_across = pairs_tukey,
    cld          = cld_tab
  )
}

# =========================================================
# Stage 2A: one-stage MET BLUPs (entries random)
# =========================================================

stage2_blups <- function(trial_df) {
  
  df <- trial_df %>%
    dplyr::filter(!is.na(yield)) %>%
    dplyr::mutate(
      env   = factor(env),
      rep   = factor(rep),
      entry = factor(entry)
    )
  
  n_env    <- dplyr::n_distinct(df$env)
  n_rep    <- dplyr::n_distinct(df$rep)
  n_entry  <- dplyr::n_distinct(df$entry)
  
  # ---------------------------------------------------------
  # Case 0: Not enough entries for any random-effects model
  # ---------------------------------------------------------
  if (n_entry < 2L) {
    warning("BLUPs: fewer than 2 entries. Returning raw mean.")
    m <- df %>% summarise(
      BLUP = mean(yield),
      SE   = sd(yield) / sqrt(n())
    )
    return(data.frame(
      entry = levels(df$entry),
      BLUP  = m$BLUP,
      SE    = m$SE
    ))
  }
  
  # ---------------------------------------------------------
  # Case 1: Only one environment
  # -> Fit single-environment BLUP: entry + rep
  # ---------------------------------------------------------
  if (n_env < 2L) {
    warning(
      "BLUPs: only one environment detected. ",
      "Fitting single-environment BLUP model yield ~ 1 + (1 | rep) + (1 | entry)."
    )
    
    # If only one rep level, remove (1|rep)
    if (n_rep < 2L) {
      warning("BLUPs: only one rep. Model reduces to (1 | entry).")
      form <- yield ~ 1 + (1 | entry)
    } else {
      form <- yield ~ 1 + (1 | rep) + (1 | entry)
    }
    
    fit <- try(lme4::lmer(form, data=df, REML=TRUE), silent=TRUE)
    
    if (!inherits(fit,"try-error")) {
      re <- ranef(fit, condVar = TRUE)$entry
      blups_vec <- as.numeric(re[,"(Intercept)"])
      entry_ids <- rownames(re)
      pv <- attr(ranef(fit, condVar=TRUE)$entry, "postVar")
      se_vec <- if(is.null(pv)) rep(NA_real_, length(blups_vec)) else 
        sqrt(vapply(seq_len(dim(pv)[3]), function(i) pv[1,1,i], numeric(1)))
      
      return(data.frame(entry=entry_ids, BLUP=blups_vec, SE=se_vec))
    }
    
    warning("BLUPs: single-environment lmer failed. Falling back to lm().")
    
    # fallback: fixed-effect entry model
    lm_fit <- lm(yield ~ entry, data=df)
    emm <- emmeans(lm_fit, ~ entry)
    tbl <- summary(emm)
    return(data.frame(
      entry   = tbl$entry,
      BLUP    = tbl$emmean,
      SE      = tbl$SE
    ))
  }
  
  # ---------------------------------------------------------
  # Case 2: Normal multi-environment model
  # ---------------------------------------------------------
  form_full <- yield ~ 1 +
    (1 | env) +
    (1 | env:rep) +
    (1 | entry) +
    (1 | entry:env)
  
  fit2 <- try(
    lme4::lmer(
      form_full,
      data=df,
      REML=TRUE,
      control=lme4::lmerControl(check.nobs.vs.nRE="ignore")
    ),
    silent=TRUE
  )
  
  if (!inherits(fit2,"try-error")) {
    re <- ranef(fit2, condVar=TRUE)$entry
    blups_vec <- as.numeric(re[,"(Intercept)"])
    entry_ids <- rownames(re)
    pv <- attr(ranef(fit2, condVar=TRUE)$entry,"postVar")
    se_vec <- if(is.null(pv)) rep(NA_real_, length(blups_vec)) else
      sqrt(vapply(seq_len(dim(pv)[3]), function(i) pv[1,1,i], numeric(1)))
    
    return(data.frame(entry=entry_ids, BLUP=blups_vec, SE=se_vec))
  }
  
  warning("BLUPs: multi-environment lmer failed. Falling back to lm().")
  
  lm_fit <- lm(yield ~ entry, data=df)
  emm <- emmeans(lm_fit, ~ entry)
  tbl <- summary(emm)
  
  data.frame(
    entry   = tbl$entry,
    BLUP    = tbl$emmean,
    SE      = tbl$SE
  )
}

# =========================================================
# Main harness configuration (sim vs real)
# =========================================================

DATA_MODE  <- "sim"   # "sim" or "real"

SIM_IN_CSV  <- "sim_data/trial_sim.csv"
SIM_OUT_DIR <- "sim_output"

REAL_IN_CSV  <- "real_data/MilletVT2025.csv"
REAL_OUT_DIR <- "real_output"

# =========================================================
# Main harness
# =========================================================

analyze_trial <- function(in_csv   = SIM_IN_CSV,
                          out_dir  = SIM_OUT_DIR,
                          alpha    = 0.05) {
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_stage1 <- file.path(out_dir, "lsm_stage1.csv")
  out_across <- file.path(out_dir, "lsm_across.csv")
  out_diffs  <- file.path(out_dir, "diffs_across.csv")
  out_cld    <- file.path(out_dir, "cld_across.csv")
  out_blups  <- file.path(out_dir, "entry_blups.csv")
  out_report <- file.path(out_dir, "run_report.txt")
  
  trial <- read.csv(in_csv, stringsAsFactors = FALSE)
  
  stopifnot(all(c("site","year","env","rep","row","col","entry","yield") %in% names(trial)))
  
  # -------------------------
  # Stage 1: per-env models
  # -------------------------
  s1 <- analyze_stage1(trial, covlist, alpha = alpha)
  lsm_stage1   <- s1$lsm_stage1
  best_cov_tbl <- s1$best_cov
  
  write.csv(lsm_stage1, out_stage1, row.names = FALSE)
  
  # Distinct environments actually used in Stage 1
  env_ids <- sort(unique(lsm_stage1$env))
  n_env   <- length(env_ids)
  
  # -------------------------
  # Stage 2B: across-env LS means
  #   • If > 1 env: meta analysis as before
  #   • If 1 env : skip meta, do within-env LS means + CLD
  # -------------------------
  if (n_env > 1L) {
    stage2_mode <- "multi_env_meta"
    s2 <- stage2_meta(lsm_stage1, alpha = alpha)
  } else {
    stage2_mode <- paste0("single_env_within_trial (env = ", env_ids[1], ")")
    message(
      "Only one environment detected (", env_ids[1],
      "); skipping Stage 2 meta-analysis and computing LS-means + Tukey CLD ",
      "directly from raw data for this environment."
    )
    s2 <- stage2_single_env(trial, env_id = env_ids[1], alpha = alpha)
  }
  
  write.csv(s2$lsm_across,   out_across, row.names = FALSE)
  write.csv(s2$diffs_across, out_diffs,  row.names = FALSE)
  write.csv(s2$cld,          out_cld,    row.names = FALSE)
  
  # -------------------------
  # Stage 2A: BLUPs
  # -------------------------
  blups <- stage2_blups(trial)
  write.csv(blups, out_blups, row.names = FALSE)
  
  # -------------------------
  # Run report: summarize what happened
  # -------------------------
  
  # Basic data summary
  total_rows    <- nrow(trial)
  total_missing <- sum(is.na(trial$yield))
  used_rows     <- total_rows - total_missing
  
  env_stats <- trial %>%
    dplyr::group_by(env) %>%
    dplyr::summarise(
      n_plots   = dplyr::n(),
      n_missing = sum(is.na(yield)),
      n_used    = n_plots - n_missing,
      n_rep     = dplyr::n_distinct(rep),
      n_entry   = dplyr::n_distinct(entry),
      .groups   = "drop"
    )
  
  # Stage 1 covariance usage counts
  spatial_types <- covlist
  n_spatial <- sum(best_cov_tbl$best_cov %in% spatial_types, na.rm = TRUE)
  n_rcbd    <- sum(best_cov_tbl$best_cov == "none_rcbd", na.rm = TRUE)
  n_lm      <- sum(best_cov_tbl$best_cov == "none_lm",   na.rm = TRUE)
  
  # Guess whether this run is "sim" or "real" based on path
  data_type_guess <- if (grepl("sim_data", in_csv)) "simulated" else "real_or_external"
  
  # BLUP mode (single vs multi-env) inferred from trial
  n_env_trial <- dplyr::n_distinct(trial$env)
  n_rep_trial <- dplyr::n_distinct(trial$rep)
  blup_mode <- if (n_env_trial < 2L) {
    if (n_rep_trial < 2L) {
      "single_env_single_rep (random entry only; may fall back to lm)"
    } else {
      "single_env_multi_rep (random entry and random rep; may fall back to lm)"
    }
  } else {
    "multi_env (random env, env:rep, entry, entry:env; may fall back to lm)"
  }
  
  report_lines <- c(
    "Variety Trial Analysis run report",
    "--------------------------------",
    paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Input file:", normalizePath(in_csv, mustWork = FALSE)),
    paste("Output directory:", normalizePath(out_dir, mustWork = FALSE)),
    paste("Data type guess:", data_type_guess),
    paste("Alpha level:", alpha),
    "",
    "Data summary:",
    paste("  Total rows:", total_rows),
    paste("  Total missing yields:", total_missing),
    paste("  Rows used in modeling (non-missing yield):", used_rows),
    paste("  Number of environments with LS-means:", n_env),
    "",
    "Environment-level summary:"
  )
  
  if (nrow(env_stats) > 0L) {
    env_lines <- apply(env_stats, 1, function(row) {
      paste0(
        "  Env ", row[["env"]], ": ",
        row[["n_plots"]], " plots (", row[["n_used"]], " used, ",
        row[["n_missing"]], " missing); reps = ", row[["n_rep"]],
        ", entries = ", row[["n_entry"]]
      )
    })
    report_lines <- c(report_lines, env_lines)
  } else {
    report_lines <- c(report_lines, "  <no environments found in trial data>")
  }
  
  report_lines <- c(
    report_lines,
    "",
    "Stage 1 (per-environment models):",
    paste("  Spatial covariance candidates:", paste(covlist, collapse = ", ")),
    paste("  Environments using spatial covariance (expa/exp/sph/gau):", n_spatial),
    paste("  Environments using RCBD fallback (none_rcbd):", n_rcbd),
    paste("  Environments using fixed-effects ANOVA fallback (none_lm):", n_lm),
    "  Best covariance per environment:"
  )
  
  if (nrow(best_cov_tbl) > 0L) {
    cov_lines <- apply(best_cov_tbl, 1, function(row) {
      aicc_str <- if (is.na(row[["best_AICc"]])) "NA" else sprintf("%.2f", as.numeric(row[["best_AICc"]]))
      paste0("    Env ", row[["env"]], ": best_cov = ", row[["best_cov"]],
             ", best_AICc = ", aicc_str)
    })
    report_lines <- c(report_lines, cov_lines)
  } else {
    report_lines <- c(report_lines, "    <no Stage 1 covariance information available>")
  }
  
  report_lines <- c(
    report_lines,
    "",
    "Stage 2B (across-environment LS-means):",
    paste("  Mode:", stage2_mode),
    if (n_env > 1L) {
      "  Description: inverse-variance meta-analysis across environments with random env (nlme::lme) when possible; Tukey-adjusted CLD via emmeans + multcompView."
    } else {
      "  Description: only one environment present; meta-analysis skipped. Within-environment LS-means and Tukey CLD computed directly from raw data for that environment."
    },
    "",
    "Stage 2A (BLUPs):",
    paste("  Mode:", blup_mode),
    "  Description: entry treated as random; env and env:rep random effects included when multiple environments and reps are available.",
    "",
    "Key modeling assumptions and rules:",
    "  • Rows with NA yield are excluded from all model fitting.",
    "  • Entry is always treated as a factor.",
    "  • Stage 1 prefers spatial nlme models; if all spatial fits fail for an environment, the model falls back to RCBD (lmer) or, if that fails or reps are insufficient, to fixed-effects ANOVA (lm).",
    "  • Stage 2 meta-analysis is only meaningful when ≥ 2 environments are present; with a single environment, across-env meta is skipped and within-env LS-means + Tukey CLD are used instead.",
    "  • BLUPs treat entry as random and, when possible, include env and env:rep random effects; if mixed models fail to converge, the code falls back to lm(yield ~ entry) and uses emmeans to summarize entry effects.",
    "  • Confidence intervals and p-values follow standard approximations from nlme/lme4/emmeans and may be unstable when degrees of freedom are small or designs are highly unbalanced."
  )
  
  writeLines(report_lines, con = out_report)
  
  # -------------------------
  # Console summary
  # -------------------------
  message("Wrote:\n",
          "  ", normalizePath(out_stage1, mustWork = FALSE), "\n",
          "  ", normalizePath(out_across, mustWork = FALSE), "\n",
          "  ", normalizePath(out_diffs,  mustWork = FALSE), "\n",
          "  ", normalizePath(out_cld,    mustWork = FALSE), "\n",
          "  ", normalizePath(out_blups,  mustWork = FALSE), "\n",
          "  ", normalizePath(out_report, mustWork = FALSE))
  
  invisible(list(
    stage1       = lsm_stage1,
    best_cov     = best_cov_tbl,
    lsm_across   = s2$lsm_across,
    diffs_across = s2$diffs_across,
    cld          = s2$cld,
    blups        = blups,
    report_path  = out_report
  ))
}


# =========================================================
# Plotting helpers (unchanged logic, now follow DATA_MODE)
# =========================================================

make_summary_plots <- function(in_csv, out_dir, mode_label = "sim") {
  # helper
  clean_entry_labels <- function(x) sub("^entry", "", as.character(x))
  
  # --------- 1. Across environment LS means with Tukey CLD ----------
  lsm_ac <- read.csv(file.path(out_dir, "lsm_across.csv"),
                     stringsAsFactors = FALSE)
  cld    <- read.csv(file.path(out_dir, "cld_across.csv"),
                     stringsAsFactors = FALSE)
  
  # Ensure CLD has a 'group' column (sometimes called .group)
  if ("group" %in% names(cld)) {
    cld_fix <- cld
  } else if (".group" %in% names(cld)) {
    cld_fix <- dplyr::rename(cld, group = .group)
  } else {
    stop("CLD table must contain a 'group' or '.group' column.")
  }
  
  # Standardize CI column names
  fix_ci_names <- function(df){
    nm <- names(df)
    nm <- sub("^lower\\.CL\\.$","lower.CL", nm)
    nm <- sub("^upper\\.CL\\.$","upper.CL", nm)
    names(df) <- nm
    df
  }
  lsm_ac <- fix_ci_names(lsm_ac)
  
  # Coerce entry and numeric columns to stable types, and clean labels
  lsm_ac <- lsm_ac %>%
    dplyr::mutate(
      entry    = clean_entry_labels(entry),
      estimate = as.numeric(estimate),
      lower.CL = as.numeric(lower.CL),
      upper.CL = as.numeric(upper.CL)
    )
  
  cld_fix <- cld_fix %>%
    dplyr::mutate(entry = clean_entry_labels(entry))
  
  # Merge LS means with grouping letters and order by mean
  plot_df <- lsm_ac %>%
    select(entry, estimate, lower.CL, upper.CL) %>%
    inner_join(cld_fix, by = "entry") %>%
    mutate(entry = as.character(entry)) %>%
    arrange(estimate) %>%                   # <-- THIS restores yield ranking
    mutate(entry = factor(entry, levels = entry))
  
  p_tukey <- ggplot2::ggplot(plot_df,
                             ggplot2::aes(x = entry, y = estimate)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.2
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = group, y = upper.CL),
      vjust = -0.6,
      size = 3
    ) +
    ggplot2::labs(
      title = "Across environment LS means with Tukey CLD",
      x = "Entry",
      y = "Adjusted mean yield"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  
  
  print(p_tukey)
  
  # --------- 2. Heatmap (unchanged) ----------
  trial <- read.csv(in_csv, stringsAsFactors = FALSE)
  
  trial <- trial %>%
    dplyr::mutate(
      row_num = as.numeric(row),
      col_num = as.numeric(col),
      env     = as.factor(env)
    )
  
  env_levels <- levels(trial$env)
  if (length(env_levels) == 0L) {
    warning("Heatmap: no environments found in trial data, skipping heatmap.")
    return(invisible(list(tukey_plot = p_tukey, heatmap = NULL)))
  }
  env_ex <- env_levels[1]
  
  df_env <- trial %>%
    dplyr::filter(env == env_ex, !is.na(yield)) %>%
    dplyr::filter(!is.na(row_num), !is.na(col_num))
  
  if (nrow(df_env) == 0L) {
    warning("Heatmap: no rows with non missing yield and numeric row or column for env ",
            env_ex, ", skipping heatmap.")
    return(invisible(list(tukey_plot = p_tukey, heatmap = NULL)))
  }
  
  p_heat <- ggplot2::ggplot(df_env,
                            ggplot2::aes(x = col_num, y = row_num, fill = yield)) +
    ggplot2::geom_tile() +
    ggplot2::scale_y_reverse() +
    ggplot2::coord_equal() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Yield heatmap:", env_ex, "(", mode_label, "data)"),
      x = "Column",
      y = "Row"
    )
  
  print(p_heat)
  
  invisible(list(tukey_plot = p_tukey, heatmap = p_heat))
}

# =========================================================
# Run block using DATA_MODE switch (interactive only)
# =========================================================

if (interactive()) {
  if (DATA_MODE == "sim") {
    in_csv  <- SIM_IN_CSV
    out_dir <- SIM_OUT_DIR
  } else if (DATA_MODE == "real") {
    in_csv  <- REAL_IN_CSV
    out_dir <- REAL_OUT_DIR
  } else {
    stop("DATA_MODE must be either 'sim' or 'real'")
  }
  
  res <- analyze_trial(in_csv = in_csv, out_dir = out_dir)
  make_summary_plots(in_csv = in_csv, out_dir = out_dir, mode_label = DATA_MODE)
}
