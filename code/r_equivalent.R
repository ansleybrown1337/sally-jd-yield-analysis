# =========================================================
# Variety Trial Analysis in R (mirrors SAS workflow)
# ---------------------------------------------------------
# Requires: nlme, lme4, emmeans, multcompView, dplyr, tidyr, purrr, broom, ggplot2
# Input  : trial.csv with columns: site, year, env, rep, row, col, entry, yield
# Outputs:
#   lsm_stage1.csv              - per-env LS-means (estimate, SE, df, covariance chosen)
#   lsm_across.csv              - across-env LS-means with Tukey-adjusted CIs
#   diffs_across.csv            - Tukey pairwise comparisons (across-env model)
#   cld_across.csv              - Compact letter display for entries
#   entry_blups.csv             - across-env entry BLUPs (one-stage MET model)
#   stage1_model_selection.csv  - full Stage 1 model-selection table (incl. failures)
#   model_specs_stage1.csv      - explicit per-env model specification record
#   diagnostics_stage1.csv      - plot-level diagnostics (fitted + residuals)
#   cv_by_env.csv               - Stage 1 RMSE and CV by environment
#   run_report.txt              - run audit trail
#
# Notes:
#   - Stage 1 attempts spatial nlme::lme fits for each environment and chooses by AICc.
#   - If all spatial fits fail in an environment, falls back to RCBD (lmer) or lm.
#   - Stage 2B uses inverse-variance meta-analysis across environments when possible.
#   - If only one environment contributes usable Stage 1 LS-means, Stage 2B is computed
#     within that environment directly from raw data (RCBD or lm).
#   - Stage 2A computes entry BLUPs from a one-stage MET model when feasible.
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
  library(tibble)
})

alpha   <- 0.05
covlist <- c("expa","exp","sph","gau")  # exp(a) placeholder, exp, spherical, Gaussian

# =========================================================
# Helper: correlation structure for nlme spatial models
# =========================================================
make_cor <- function(type, row, col) {
  # row/col are used in the calling scope; keep signature stable
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
  envs <- sort(unique(trial_df$env))

  res_lsm        <- vector("list", length(envs))
  best_cov_list  <- vector("list", length(envs))
  sel_all_list   <- vector("list", length(envs))
  spec_list      <- vector("list", length(envs))
  diag_list      <- vector("list", length(envs))
  cv_list        <- vector("list", length(envs))

  for (i in seq_along(envs)) {
    this_env <- envs[i]

    df_env_raw <- trial_df %>%
      dplyr::filter(env == this_env) %>%
      dplyr::mutate(
        entry = factor(entry),
        rep   = factor(rep)
      )

    df_env <- df_env_raw %>%
      dplyr::filter(!is.na(yield))

    if (nrow(df_env) == 0L) {
      warning("No non-missing yields for env ", this_env, "; skipping Stage 1 for this env.")
      next
    }

    # -----------------------------------------------------
    # 1) Try spatial nlme models for this environment
    # -----------------------------------------------------
    fits <- purrr::map(covlist, ~ fit_one_env_spatial(df_env, .x))
    names(fits) <- covlist

    # Collect full model-selection diagnostics (including failures)
    sel_rows <- purrr::imap_dfr(fits, function(fit, cov_name) {
      if (is.null(fit)) {
        return(tibble::tibble(
          env = this_env, cov = cov_name, engine = "nlme::lme",
          converged = FALSE, n = NA_integer_, k = NA_integer_,
          logLik = NA_real_, AIC = NA_real_, AICc = NA_real_, BIC = NA_real_,
          best = FALSE
        ))
      }
      k_tot <- attr(logLik(fit), "df")
      n_obs <- nobs(fit)
      ll    <- as.numeric(logLik(fit))
      aic   <- AIC(fit)
      bic   <- BIC(fit)
      aicc  <- AICc_nlme(fit)
      tibble::tibble(
        env = this_env, cov = cov_name, engine = "nlme::lme",
        converged = TRUE, n = n_obs, k = k_tot,
        logLik = ll, AIC = aic, AICc = aicc, BIC = bic,
        best = FALSE
      )
    })

    ok_idx  <- !purrr::map_lgl(fits, is.null)
    fits_ok <- fits[ok_idx]

    if (length(fits_ok) > 0L) {
      aicc <- purrr::map_dbl(fits_ok, AICc_nlme)
      best_type <- names(which.min(aicc))
      best_fit  <- fits_ok[[best_type]]

      sel_rows <- sel_rows %>%
        dplyr::mutate(best = (cov == best_type) & converged)

      lsm <- lsm_from_fit(best_fit, alpha = alpha) %>%
        dplyr::mutate(env = this_env, cov = best_type) %>%
        dplyr::select(env, cov, entry, estimate, stderr, df, lower.CL, upper.CL)

      res_lsm[[i]] <- lsm

      best_cov_list[[i]] <- tibble::tibble(
        env       = this_env,
        best_cov  = best_type,
        best_AICc = min(aicc)
      )

      spec_list[[i]] <- tibble::tibble(
        env        = this_env,
        stage      = "stage1",
        engine     = "nlme::lme",
        fixed      = "yield ~ entry",
        random     = "~ 1 | rep",
        residual   = paste0("correlation = ", best_type, " (nlme cor*; nugget=TRUE)"),
        method     = "REML",
        cov_selected = best_type
      )

      diag_df <- df_env %>%
        dplyr::mutate(
          fitted     = as.numeric(fitted(best_fit)),
          resid      = as.numeric(residuals(best_fit, type = "response")),
          resid_norm = as.numeric(residuals(best_fit, type = "normalized")),
          resid_kind = "normalized (nlme)",
          flag_outlier = abs(resid_norm) > 3
        ) %>%
        dplyr::select(env, site, year, rep, row, col, entry, yield,
                      fitted, resid, resid_norm, resid_kind, flag_outlier)

      diag_list[[i]] <- diag_df

      # CV summary (names match frontend expectation: mean_y, rmse, cv_pct)
      rmse   <- sqrt(mean(diag_df$resid^2, na.rm = TRUE))
      mean_y <- mean(diag_df$yield, na.rm = TRUE)
      cv_pct <- if (is.finite(mean_y) && mean_y != 0) 100 * rmse / mean_y else NA_real_

      cv_list[[i]] <- tibble::tibble(
        env          = this_env,
        cov_selected = best_type,
        n_used       = sum(!is.na(diag_df$yield)),
        mean_y       = mean_y,
        rmse         = rmse,
        cv_pct       = cv_pct
      )

      sel_all_list[[i]] <- sel_rows
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

    n_reps <- dplyr::n_distinct(df_env$rep)

    if (n_reps >= 2L) {
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
        lsm <- lsm_from_fit(rcbd_fit, alpha = alpha) %>%
          dplyr::mutate(env = this_env, cov = "none_rcbd") %>%
          dplyr::select(env, cov, entry, estimate, stderr, df, lower.CL, upper.CL)

        res_lsm[[i]] <- lsm

        best_cov_list[[i]] <- tibble::tibble(
          env       = this_env,
          best_cov  = "none_rcbd",
          best_AICc = NA_real_
        )

        spec_list[[i]] <- tibble::tibble(
          env        = this_env,
          stage      = "stage1",
          engine     = "lme4::lmer",
          fixed      = "yield ~ entry",
          random     = "(1 | rep)",
          residual   = "none (no spatial correlation)",
          method     = "REML",
          cov_selected = "none_rcbd"
        )

        sigma_hat <- tryCatch(sigma(rcbd_fit), error = function(e) NA_real_)
        resid_raw <- as.numeric(residuals(rcbd_fit))
        resid_std <- if (is.finite(sigma_hat) && sigma_hat > 0) resid_raw / sigma_hat else NA_real_

        diag_df <- df_env %>%
          dplyr::mutate(
            fitted     = as.numeric(fitted(rcbd_fit)),
            resid      = resid_raw,
            resid_norm = resid_std,
            resid_kind = "standardized (lmer resid / sigma)",
            flag_outlier = abs(resid_norm) > 3
          ) %>%
          dplyr::select(env, site, year, rep, row, col, entry, yield,
                        fitted, resid, resid_norm, resid_kind, flag_outlier)

        diag_list[[i]] <- diag_df

        rmse   <- sqrt(mean(diag_df$resid^2, na.rm = TRUE))
        mean_y <- mean(diag_df$yield, na.rm = TRUE)
        cv_pct <- if (is.finite(mean_y) && mean_y != 0) 100 * rmse / mean_y else NA_real_

        cv_list[[i]] <- tibble::tibble(
          env          = this_env,
          cov_selected = "none_rcbd",
          n_used       = sum(!is.na(diag_df$yield)),
          mean_y       = mean_y,
          rmse         = rmse,
          cv_pct       = cv_pct
        )

        sel_rows <- sel_rows %>% dplyr::mutate(best = FALSE)
        sel_rows <- dplyr::bind_rows(
          sel_rows,
          tibble::tibble(
            env = this_env, cov = "none_rcbd", engine = "lme4::lmer",
            converged = TRUE, n = nobs(rcbd_fit), k = attr(logLik(rcbd_fit), "df"),
            logLik = as.numeric(logLik(rcbd_fit)),
            AIC = AIC(rcbd_fit), AICc = NA_real_, BIC = BIC(rcbd_fit),
            best = TRUE
          )
        )

        sel_all_list[[i]] <- sel_rows
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

    lm_fit <- lm(yield ~ entry, data = df_env)

    lsm <- lsm_from_fit(lm_fit, alpha = alpha) %>%
      dplyr::mutate(env = this_env, cov = "none_lm") %>%
      dplyr::select(env, cov, entry, estimate, stderr, df, lower.CL, upper.CL)

    res_lsm[[i]] <- lsm

    best_cov_list[[i]] <- tibble::tibble(
      env       = this_env,
      best_cov  = "none_lm",
      best_AICc = NA_real_
    )

    spec_list[[i]] <- tibble::tibble(
      env        = this_env,
      stage      = "stage1",
      engine     = "stats::lm",
      fixed      = "yield ~ entry",
      random     = "none",
      residual   = "none (no spatial correlation)",
      method     = "OLS",
      cov_selected = "none_lm"
    )

    stud <- tryCatch(as.numeric(rstudent(lm_fit)), error = function(e) rep(NA_real_, nrow(df_env)))
    diag_df <- df_env %>%
      dplyr::mutate(
        fitted     = as.numeric(fitted(lm_fit)),
        resid      = as.numeric(residuals(lm_fit)),
        resid_norm = stud,
        resid_kind = "studentized (lm rstudent)",
        flag_outlier = abs(resid_norm) > 3
      ) %>%
      dplyr::select(env, site, year, rep, row, col, entry, yield,
                    fitted, resid, resid_norm, resid_kind, flag_outlier)

    diag_list[[i]] <- diag_df

    rmse   <- sqrt(mean(diag_df$resid^2, na.rm = TRUE))
    mean_y <- mean(diag_df$yield, na.rm = TRUE)
    cv_pct <- if (is.finite(mean_y) && mean_y != 0) 100 * rmse / mean_y else NA_real_

    cv_list[[i]] <- tibble::tibble(
      env          = this_env,
      cov_selected = "none_lm",
      n_used       = sum(!is.na(diag_df$yield)),
      mean_y       = mean_y,
      rmse         = rmse,
      cv_pct       = cv_pct
    )

    sel_rows <- sel_rows %>% dplyr::mutate(best = FALSE)
    sel_rows <- dplyr::bind_rows(
      sel_rows,
      tibble::tibble(
        env = this_env, cov = "none_lm", engine = "stats::lm",
        converged = TRUE, n = nobs(lm_fit), k = length(coef(lm_fit)),
        logLik = as.numeric(logLik(lm_fit)),
        AIC = AIC(lm_fit), AICc = NA_real_, BIC = BIC(lm_fit),
        best = TRUE
      )
    )
    sel_all_list[[i]] <- sel_rows
  }

  lsm_stage1 <- dplyr::bind_rows(res_lsm) %>%
    dplyr::mutate(var_lsmean = stderr^2)

  list(
    lsm_stage1       = lsm_stage1,
    best_cov         = dplyr::bind_rows(best_cov_list),
    model_selection  = dplyr::bind_rows(sel_all_list),
    model_specs      = dplyr::bind_rows(spec_list),
    diagnostics      = dplyr::bind_rows(diag_list),
    cv_by_env        = dplyr::bind_rows(cv_list)
  )
}

# =========================================================
# Stage 2B: inverse-variance meta with Tukey letters + LSD
# =========================================================
stage2_meta <- function(lsm_stage1, alpha = 0.05) {
  lsm_stage1 <- lsm_stage1 %>%
    dplyr::mutate(entry = factor(entry), env = factor(env))

  lsm2 <- lsm_stage1 %>%
    dplyr::filter(!is.na(var_lsmean), var_lsmean > 0)

  if (nrow(lsm2) == 0L) {
    stop("Stage 2 meta-analysis: no rows with positive var_lsmean. Check Stage 1 SEs.")
  }

  n_env <- dplyr::n_distinct(lsm2$env)
  n_ent <- dplyr::n_distinct(lsm2$entry)

  if (n_ent < 2L) {
    warning("Stage 2 meta-analysis: < 2 entries. Reporting marginal mean only.")
    fit <- lm(estimate ~ 1, data = lsm2, weights = 1 / var_lsmean)
    emm <- emmeans::emmeans(fit, ~ 1)
    base_mean <- as.data.frame(summary(emm, infer = TRUE, level = 1 - alpha))

    lsm_tab <- data.frame(
      entry     = levels(lsm2$entry),
      estimate  = base_mean$emmean,
      stderr    = base_mean$SE,
      df        = base_mean$df,
      lower.CL  = base_mean$lower.CL,
      upper.CL  = base_mean$upper.CL,
      LSD_0.05  = NA_real_,
      LSD_0.30  = NA_real_,
      stringsAsFactors = FALSE
    )

    cld_tab <- data.frame(entry = lsm_tab$entry, group = "a", row.names = NULL, check.names = FALSE)

    return(list(lsm_across = lsm_tab, diffs_across = data.frame(), cld = cld_tab,
                lsd_0.05 = NA_real_, lsd_0.30 = NA_real_))
  }

  if (n_env >= 2L) {
    fit <- nlme::lme(
      estimate ~ entry,
      random   = ~ 1 | env,
      weights  = nlme::varFixed(~ var_lsmean),
      data     = lsm2,
      method   = "REML"
    )
  } else {
    warning("Stage 2 meta-analysis: < 2 environments with positive var_lsmean; using weighted lm().")
    fit <- lm(estimate ~ entry, data = lsm2, weights = 1 / var_lsmean)
  }

  emm <- emmeans::emmeans(fit, ~ entry)

  lsm_tab <- as.data.frame(summary(emm, infer = TRUE, level = 1 - alpha)) %>%
    dplyr::rename(estimate = emmean, stderr = SE) %>%
    dplyr::select(entry, estimate, stderr, df, lower.CL, upper.CL)

  pairs_tukey <- tryCatch(as.data.frame(summary(pairs(emm, adjust = "tukey"))), error = function(e) NULL)

  lsd_0.05 <- NA_real_
  lsd_0.30 <- NA_real_

  pairs_none <- tryCatch(as.data.frame(summary(pairs(emm, adjust = "none"))), error = function(e) NULL)
  if (!is.null(pairs_none) && all(c("SE", "df") %in% names(pairs_none))) {
    valid <- with(pairs_none, !is.na(SE) & !is.na(df) & df > 0)
    if (any(valid)) {
      se_eff <- mean(pairs_none$SE[valid])
      df_eff <- mean(pairs_none$df[valid])
      lsd_0.05 <- stats::qt(1 - 0.05/2, df_eff) * se_eff
      lsd_0.30 <- stats::qt(1 - 0.30/2, df_eff) * se_eff
    }
  }

  lsm_tab <- lsm_tab %>% dplyr::mutate(LSD_0.05 = lsd_0.05, LSD_0.30 = lsd_0.30)

  if (is.null(pairs_tukey) || nrow(pairs_tukey) == 0L) {
    warning("Stage 2 meta-analysis: no valid Tukey contrasts for CLD; assigning group 'a' to all.")
    cld_tab <- data.frame(entry = lsm_tab$entry, group = "a", row.names = NULL, check.names = FALSE)
    return(list(lsm_across = lsm_tab, diffs_across = if (is.null(pairs_tukey)) data.frame() else pairs_tukey,
                cld = cld_tab, lsd_0.05 = lsd_0.05, lsd_0.30 = lsd_0.30))
  }

  pv <- pairs_tukey$p.value
  name_vec <- gsub(" ", "", pairs_tukey$contrast)
  valid_p <- !is.na(pv)
  pv <- pv[valid_p]
  name_vec <- name_vec[valid_p]

  if (length(pv) == 0L) {
    warning("Stage 2 meta-analysis: all Tukey p-values are NA; assigning group 'a' to all.")
    cld_tab <- data.frame(entry = lsm_tab$entry, group = "a", row.names = NULL, check.names = FALSE)
    return(list(lsm_across = lsm_tab, diffs_across = pairs_tukey, cld = cld_tab,
                lsd_0.05 = lsd_0.05, lsd_0.30 = lsd_0.30))
  }

  names(pv) <- name_vec
  let <- tryCatch(multcompView::multcompLetters(pv, threshold = alpha), error = function(e) NULL)

  if (is.null(let)) {
    warning("Stage 2 meta-analysis: multcompLetters failed; assigning group 'a' to all.")
    cld_tab <- data.frame(entry = lsm_tab$entry, group = "a", row.names = NULL, check.names = FALSE)
  } else {
    cld_tab <- data.frame(entry = names(let$Letters), group = unname(let$Letters), row.names = NULL, check.names = FALSE)
  }

  list(lsm_across = lsm_tab, diffs_across = pairs_tukey, cld = cld_tab, lsd_0.05 = lsd_0.05, lsd_0.30 = lsd_0.30)
}

# =========================================================
# Stage 2B helper: single-environment LS-means + CLD + LSD
# =========================================================
stage2_single_env <- function(trial_df, env_id, alpha = 0.05) {
  df_env <- trial_df %>%
    dplyr::filter(env == env_id, !is.na(yield)) %>%
    dplyr::mutate(entry = factor(entry), rep = factor(rep))

  if (nrow(df_env) == 0L) stop("stage2_single_env: no non-missing yields for env = ", env_id)

  n_rep <- dplyr::n_distinct(df_env$rep)
  if (n_rep >= 2L) {
    fit <- lme4::lmer(yield ~ entry + (1 | rep), data = df_env, REML = TRUE)
  } else {
    warning("stage2_single_env: env ", env_id, " has < 2 reps; fitting lm(yield ~ entry).")
    fit <- lm(yield ~ entry, data = df_env)
  }

  emm <- emmeans::emmeans(fit, ~ entry)

  lsm_tab <- as.data.frame(summary(emm, infer = TRUE, level = 1 - alpha)) %>%
    dplyr::rename(estimate = emmean, stderr = SE) %>%
    dplyr::select(entry, estimate, stderr, df, lower.CL, upper.CL)

  pairs_tukey <- tryCatch(as.data.frame(summary(pairs(emm, adjust = "tukey"))), error = function(e) NULL)

  lsd_0.05 <- NA_real_
  lsd_0.30 <- NA_real_
  pairs_none <- tryCatch(as.data.frame(summary(pairs(emm, adjust = "none"))), error = function(e) NULL)

  if (!is.null(pairs_none) && all(c("SE", "df") %in% names(pairs_none))) {
    valid <- with(pairs_none, !is.na(SE) & !is.na(df) & df > 0)
    if (any(valid)) {
      se_eff <- mean(pairs_none$SE[valid])
      df_eff <- mean(pairs_none$df[valid])
      lsd_0.05 <- stats::qt(1 - 0.05/2, df_eff) * se_eff
      lsd_0.30 <- stats::qt(1 - 0.30/2, df_eff) * se_eff
    }
  }

  lsm_tab <- lsm_tab %>% dplyr::mutate(LSD_0.05 = lsd_0.05, LSD_0.30 = lsd_0.30)

  if (is.null(pairs_tukey) || nrow(pairs_tukey) == 0L) {
    warning("stage2_single_env: no valid Tukey contrasts for CLD; assigning group 'a' to all.")
    cld_tab <- data.frame(entry = lsm_tab$entry, group = "a", row.names = NULL, check.names = FALSE)
    return(list(lsm_across = lsm_tab, diffs_across = if (is.null(pairs_tukey)) data.frame() else pairs_tukey,
                cld = cld_tab, lsd_0.05 = lsd_0.05, lsd_0.30 = lsd_0.30))
  }

  pv <- pairs_tukey$p.value
  name_vec <- gsub(" ", "", pairs_tukey$contrast)
  valid_p <- !is.na(pv)
  pv <- pv[valid_p]
  name_vec <- name_vec[valid_p]

  if (length(pv) == 0L) {
    warning("stage2_single_env: all Tukey p-values are NA; assigning group 'a' to all.")
    cld_tab <- data.frame(entry = lsm_tab$entry, group = "a", row.names = NULL, check.names = FALSE)
    return(list(lsm_across = lsm_tab, diffs_across = pairs_tukey, cld = cld_tab,
                lsd_0.05 = lsd_0.05, lsd_0.30 = lsd_0.30))
  }

  names(pv) <- name_vec
  let <- tryCatch(multcompView::multcompLetters(pv, threshold = alpha), error = function(e) NULL)

  if (is.null(let)) {
    warning("stage2_single_env: multcompLetters failed; assigning group 'a' to all.")
    cld_tab <- data.frame(entry = lsm_tab$entry, group = "a", row.names = NULL, check.names = FALSE)
  } else {
    cld_tab <- data.frame(entry = names(let$Letters), group = unname(let$Letters), row.names = NULL, check.names = FALSE)
  }

  list(lsm_across = lsm_tab, diffs_across = pairs_tukey, cld = cld_tab, lsd_0.05 = lsd_0.05, lsd_0.30 = lsd_0.30)
}

# =========================================================
# Stage 2A: one-stage MET BLUPs (entries random)
# =========================================================
stage2_blups <- function(trial_df) {
  df <- trial_df %>%
    dplyr::filter(!is.na(yield)) %>%
    dplyr::mutate(env = factor(env), rep = factor(rep), entry = factor(entry))

  n_env   <- dplyr::n_distinct(df$env)
  n_rep   <- dplyr::n_distinct(df$rep)
  n_entry <- dplyr::n_distinct(df$entry)

  if (n_entry < 2L) {
    warning("BLUPs: fewer than 2 entries. Returning raw mean.")
    m <- df %>% summarise(BLUP = mean(yield), SE = sd(yield) / sqrt(n()))
    return(data.frame(entry = levels(df$entry), BLUP = m$BLUP, SE = m$SE))
  }

  if (n_env < 2L) {
    warning("BLUPs: only one environment detected. Fitting yield ~ 1 + (1|entry) (+ rep if possible).")
    form <- if (n_rep < 2L) (yield ~ 1 + (1 | entry)) else (yield ~ 1 + (1 | rep) + (1 | entry))
    fit <- try(lme4::lmer(form, data = df, REML = TRUE), silent = TRUE)

    if (!inherits(fit, "try-error")) {
      re <- ranef(fit, condVar = TRUE)$entry
      blups_vec <- as.numeric(re[, "(Intercept)"])
      entry_ids <- rownames(re)
      pv <- attr(ranef(fit, condVar = TRUE)$entry, "postVar")
      se_vec <- if (is.null(pv)) rep(NA_real_, length(blups_vec)) else
        sqrt(vapply(seq_len(dim(pv)[3]), function(i) pv[1, 1, i], numeric(1)))
      return(data.frame(entry = entry_ids, BLUP = blups_vec, SE = se_vec))
    }

    warning("BLUPs: single-environment lmer failed. Falling back to lm().")
    lm_fit <- lm(yield ~ entry, data = df)
    emm <- emmeans(lm_fit, ~ entry)
    tbl <- summary(emm)
    return(data.frame(entry = tbl$entry, BLUP = tbl$emmean, SE = tbl$SE))
  }

  form_full <- yield ~ 1 +
    (1 | env) +
    (1 | env:rep) +
    (1 | entry) +
    (1 | entry:env)

  fit2 <- try(
    lme4::lmer(form_full, data = df, REML = TRUE,
               control = lme4::lmerControl(check.nobs.vs.nRE = "ignore")),
    silent = TRUE
  )

  if (!inherits(fit2, "try-error")) {
    re <- ranef(fit2, condVar = TRUE)$entry
    blups_vec <- as.numeric(re[, "(Intercept)"])
    entry_ids <- rownames(re)
    pv <- attr(ranef(fit2, condVar = TRUE)$entry, "postVar")
    se_vec <- if (is.null(pv)) rep(NA_real_, length(blups_vec)) else
      sqrt(vapply(seq_len(dim(pv)[3]), function(i) pv[1, 1, i], numeric(1)))
    return(data.frame(entry = entry_ids, BLUP = blups_vec, SE = se_vec))
  }

  warning("BLUPs: multi-environment lmer failed. Falling back to lm().")
  lm_fit <- lm(yield ~ entry, data = df)
  emm <- emmeans(lm_fit, ~ entry)
  tbl <- summary(emm)
  data.frame(entry = tbl$entry, BLUP = tbl$emmean, SE = tbl$SE)
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
analyze_trial <- function(in_csv = SIM_IN_CSV, out_dir = SIM_OUT_DIR, alpha = 0.05) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_stage1   <- file.path(out_dir, "lsm_stage1.csv")
  out_across   <- file.path(out_dir, "lsm_across.csv")
  out_diffs    <- file.path(out_dir, "diffs_across.csv")
  out_cld      <- file.path(out_dir, "cld_across.csv")
  out_blups    <- file.path(out_dir, "entry_blups.csv")
  out_report   <- file.path(out_dir, "run_report.txt")
  out_modelsel <- file.path(out_dir, "stage1_model_selection.csv")
  out_specs    <- file.path(out_dir, "model_specs_stage1.csv")
  out_diag     <- file.path(out_dir, "diagnostics_stage1.csv")
  out_cv       <- file.path(out_dir, "cv_by_env.csv")

  trial <- read.csv(in_csv, stringsAsFactors = FALSE)
  stopifnot(all(c("site","year","env","rep","row","col","entry","yield") %in% names(trial)))

  s1 <- analyze_stage1(trial, covlist, alpha = alpha)
  lsm_stage1   <- s1$lsm_stage1
  best_cov_tbl <- s1$best_cov

  write.csv(lsm_stage1, out_stage1, row.names = FALSE)
  write.csv(s1$model_selection, out_modelsel, row.names = FALSE)
  write.csv(s1$model_specs,     out_specs,    row.names = FALSE)
  write.csv(s1$diagnostics,     out_diag,     row.names = FALSE)
  write.csv(s1$cv_by_env,       out_cv,       row.names = FALSE)

  env_ids <- sort(unique(lsm_stage1$env))
  n_env   <- length(env_ids)

  if (n_env > 1L) {
    stage2_mode <- "multi_env_meta"
    s2 <- stage2_meta(lsm_stage1, alpha = alpha)
  } else {
    stage2_mode <- paste0("single_env_within_trial (env = ", env_ids[1], ")")
    message("Only one environment detected (", env_ids[1], "); Stage 2 meta skipped; computing within-env LS-means + CLD.")
    s2 <- stage2_single_env(trial, env_id = env_ids[1], alpha = alpha)
  }

  write.csv(s2$lsm_across,   out_across, row.names = FALSE)
  write.csv(s2$diffs_across, out_diffs,  row.names = FALSE)
  write.csv(s2$cld,          out_cld,    row.names = FALSE)

  blups <- stage2_blups(trial)
  write.csv(blups, out_blups, row.names = FALSE)

  # -------------------------
  # Run report: summarize what happened
  # -------------------------
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

  spatial_types <- covlist
  n_spatial <- sum(best_cov_tbl$best_cov %in% spatial_types, na.rm = TRUE)
  n_rcbd    <- sum(best_cov_tbl$best_cov == "none_rcbd", na.rm = TRUE)
  n_lm      <- sum(best_cov_tbl$best_cov == "none_lm",   na.rm = TRUE)

  data_type_guess <- if (grepl("sim_data", in_csv)) "simulated" else "real_or_external"

  n_env_trial <- dplyr::n_distinct(trial$env)
  n_rep_trial <- dplyr::n_distinct(trial$rep)
  blup_mode <- if (n_env_trial < 2L) {
    if (n_rep_trial < 2L) "single_env_single_rep (random entry only; may fall back to lm)"
    else "single_env_multi_rep (random entry and random rep; may fall back to lm)"
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
      paste0("  Env ", row[["env"]], ": ",
             row[["n_plots"]], " plots (", row[["n_used"]], " used, ",
             row[["n_missing"]], " missing); reps = ", row[["n_rep"]],
             ", entries = ", row[["n_entry"]])
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
    "  Explicit model specs written to: model_specs_stage1.csv",
    "  Full model-selection table written to: stage1_model_selection.csv",
    "  Plot-level diagnostics written to: diagnostics_stage1.csv",
    "  CV summary written to: cv_by_env.csv",
    paste("  Environments using spatial covariance (expa/exp/sph/gau):", n_spatial),
    paste("  Environments using RCBD fallback (none_rcbd):", n_rcbd),
    paste("  Environments using fixed-effects ANOVA fallback (none_lm):", n_lm),
    "  Best covariance per environment:"
  )

  if (nrow(best_cov_tbl) > 0L) {
    cov_lines <- apply(best_cov_tbl, 1, function(row) {
      aicc_str <- if (is.na(row[["best_AICc"]])) "NA" else sprintf("%.2f", as.numeric(row[["best_AICc"]]))
      paste0("    Env ", row[["env"]], ": best_cov = ", row[["best_cov"]], ", best_AICc = ", aicc_str)
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

  message("Wrote:\n",
          "  ", normalizePath(out_stage1, mustWork = FALSE), "\n",
          "  ", normalizePath(out_across, mustWork = FALSE), "\n",
          "  ", normalizePath(out_diffs,  mustWork = FALSE), "\n",
          "  ", normalizePath(out_cld,    mustWork = FALSE), "\n",
          "  ", normalizePath(out_blups,  mustWork = FALSE), "\n",
          "  ", normalizePath(out_modelsel, mustWork = FALSE), "\n",
          "  ", normalizePath(out_specs,    mustWork = FALSE), "\n",
          "  ", normalizePath(out_diag,     mustWork = FALSE), "\n",
          "  ", normalizePath(out_cv,       mustWork = FALSE), "\n",
          "  ", normalizePath(out_report, mustWork = FALSE))

  invisible(list(
    stage1       = lsm_stage1,
    best_cov     = best_cov_tbl,
    lsm_across   = s2$lsm_across,
    diffs_across = s2$diffs_across,
    cld          = s2$cld,
    blups        = blups,
    report_path  = out_report,
    stage1_model_selection = s1$model_selection,
    model_specs_stage1     = s1$model_specs,
    diagnostics_stage1     = s1$diagnostics,
    cv_by_env              = s1$cv_by_env
  ))
}

# =========================================================
# Plotting helpers
# =========================================================
make_summary_plots <- function(in_csv,
                               out_dir,
                               mode_label = "sim",
                               show_tukey = TRUE,
                               show_heatmap = NULL,
                               show_residual_heatmap = TRUE,
                               show_residual_diagnostics = TRUE,
                               show_yield_heatmap = FALSE) {
  if (!is.null(show_heatmap)) {
    show_residual_heatmap <- isTRUE(show_heatmap)
  }

  clean_entry_labels <- function(x) sub("^entry", "", as.character(x))

  lsm_ac <- read.csv(file.path(out_dir, "lsm_across.csv"), stringsAsFactors = FALSE)
  cld    <- read.csv(file.path(out_dir, "cld_across.csv"), stringsAsFactors = FALSE)

  if ("group" %in% names(cld)) {
    cld_fix <- cld
  } else if (".group" %in% names(cld)) {
    cld_fix <- dplyr::rename(cld, group = .group)
  } else {
    stop("CLD table must contain a 'group' or '.group' column.")
  }

  fix_ci_names <- function(df){
    nm <- names(df)
    nm <- sub("^lower\\.CL\\.$","lower.CL", nm)
    nm <- sub("^upper\\.CL\\.$","upper.CL", nm)
    names(df) <- nm
    df
  }
  lsm_ac <- fix_ci_names(lsm_ac)

  lsm_ac <- lsm_ac %>%
    dplyr::mutate(
      entry    = clean_entry_labels(entry),
      estimate = as.numeric(estimate),
      lower.CL = as.numeric(lower.CL),
      upper.CL = as.numeric(upper.CL)
    )

  cld_fix <- cld_fix %>% dplyr::mutate(entry = clean_entry_labels(entry))

  plot_df <- lsm_ac %>%
    dplyr::select(entry, estimate, stderr, df, lower.CL, upper.CL) %>%
    dplyr::inner_join(cld_fix, by = "entry") %>%
    dplyr::mutate(entry = as.character(entry)) %>%
    dplyr::arrange(estimate) %>%
    dplyr::mutate(entry = factor(entry, levels = entry))

  if (show_tukey) {
    
    alpha_wide <- 0.05  # 95% CI
    alpha_nar  <- 0.30  # 70% CI
    
    plot_df2 <- plot_df %>%
      dplyr::mutate(
        df_eff = suppressWarnings(as.numeric(.data$df)),
        se_eff = suppressWarnings(as.numeric(.data$stderr)),
        
        t_wide = dplyr::if_else(
          is.na(df_eff),
          stats::qnorm(1 - alpha_wide/2),
          stats::qt(1 - alpha_wide/2, df = df_eff)
        ),
        t_nar = dplyr::if_else(
          is.na(df_eff),
          stats::qnorm(1 - alpha_nar/2),
          stats::qt(1 - alpha_nar/2, df = df_eff)
        ),
        
        lower_wide = estimate - t_wide * se_eff,
        upper_wide = estimate + t_wide * se_eff,
        lower_nar  = estimate - t_nar  * se_eff,
        upper_nar  = estimate + t_nar  * se_eff
      )
    
    p_tukey <- ggplot2::ggplot(plot_df2, ggplot2::aes(y = entry, x = estimate)) +
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = lower_wide, xmax = upper_wide),
        linewidth = 5, alpha = 0.20, color = "darkred"
      ) +
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = lower_nar, xmax = upper_nar),
        linewidth = 5, alpha = 0.40, color = "darkred"
      ) +
      ggplot2::geom_point(color = "darkred", size = 2) +
      ggplot2::geom_text(
        ggplot2::aes(x = upper_wide, label = group),
        hjust = -0.2, size = 3
      ) +
      ggplot2::labs(
        title = "Across-environment LS-means with Tukey CLD",
        subtitle = "Light band: 95% CI (α = 0.05). Dark band: 70% CI (α = 0.30).",
        x = "Adjusted mean yield",
        y = "Entry (ordered by estimate)"
      ) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(size = 18, face = "bold"),
        plot.subtitle   = ggplot2::element_text(size = 14),
        axis.title.x    = ggplot2::element_text(size = 14),
        axis.title.y    = ggplot2::element_text(size = 14),
        axis.text.x     = ggplot2::element_text(size = 12),
        axis.text.y     = ggplot2::element_text(size = 12),
        panel.grid.major.y = ggplot2::element_blank(),  # reduces clutter
        panel.grid.minor   = ggplot2::element_blank()
      ) +
      ggplot2::theme_bw()
    
    print(p_tukey)
    
  } else {
    p_tukey <- NULL
  }
  

  diag_path <- file.path(out_dir, "diagnostics_stage1.csv")
  if (!file.exists(diag_path)) {
    warning("No diagnostics_stage1.csv found in output dir; skipping residual plots.")
    return(invisible(list(tukey_plot = p_tukey, residual_heatmap = NULL, residual_diag = NULL)))
  }

  diag <- read.csv(diag_path, stringsAsFactors = FALSE) %>%
    dplyr::mutate(row_num = as.numeric(row), col_num = as.numeric(col))

  env_levels <- sort(unique(diag$env))
  if (length(env_levels) == 0L) {
    warning("Diagnostics file has no environments; skipping residual plots.")
    return(invisible(list(tukey_plot = p_tukey, residual_heatmap = NULL, residual_diag = NULL)))
  }
  env_ex <- env_levels[1]

  df_env <- diag %>%
    dplyr::filter(env == env_ex) %>%
    dplyr::filter(!is.na(yield)) %>%
    dplyr::filter(!is.na(row_num), !is.na(col_num))

  if (nrow(df_env) == 0L) {
    warning("Residual plots: no rows with usable row/col/yield for env ", env_ex, ".")
    return(invisible(list(tukey_plot = p_tukey, residual_heatmap = NULL, residual_diag = NULL)))
  }

  # Residual heatmap (uses resid_norm which is always written by Stage 1)
  if (show_residual_heatmap) {
    p_resid <- ggplot2::ggplot(df_env, ggplot2::aes(x = col_num, y = row_num, fill = resid_norm)) +
      ggplot2::geom_tile() +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_equal() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0("Residual heatmap (", env_ex, "; ", mode_label, ")"),
        subtitle = unique(df_env$resid_kind)[1],
        x = "Column",
        y = "Row",
        fill = "Residual"
      )
    print(p_resid)
  } else {
    p_resid <- NULL
  }

  if (show_yield_heatmap) {
    p_yield <- ggplot2::ggplot(df_env, ggplot2::aes(x = col_num, y = row_num, fill = yield)) +
      ggplot2::geom_tile() +
      ggplot2::scale_y_reverse() +
      ggplot2::coord_equal() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0("Yield heatmap (", env_ex, "; ", mode_label, ")"),
        x = "Column",
        y = "Row",
        fill = "Yield"
      )
    print(p_yield)
  }

  if (show_residual_diagnostics) {
    p_hist <- ggplot2::ggplot(df_env, ggplot2::aes(x = resid_norm)) +
      ggplot2::geom_histogram(bins = 30) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0("Residual distribution (", env_ex, ")"),
        subtitle = unique(df_env$resid_kind)[1],
        x = "Residual (normalized/studentized/standardized)",
        y = "Count"
      )
    print(p_hist)

    p_qq <- ggplot2::ggplot(df_env, ggplot2::aes(sample = resid_norm)) +
      ggplot2::stat_qq() +
      ggplot2::stat_qq_line() +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0("Residual Q-Q plot (", env_ex, ")"),
        subtitle = unique(df_env$resid_kind)[1],
        x = "Theoretical quantiles",
        y = "Sample quantiles"
      )
    print(p_qq)

    p_rvf <- ggplot2::ggplot(df_env, ggplot2::aes(x = fitted, y = resid_norm)) +
      ggplot2::geom_point() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0("Residuals vs fitted (", env_ex, ")"),
        subtitle = unique(df_env$resid_kind)[1],
        x = "Fitted",
        y = "Residual"
      )
    print(p_rvf)
  } else {
    p_hist <- p_qq <- p_rvf <- NULL
  }

  invisible(list(
    tukey_plot = p_tukey,
    residual_heatmap = p_resid,
    residual_diag = list(hist = p_hist, qq = p_qq, rvf = p_rvf)
  ))
}

# =========================================================
# Run block using DATA_MODE switch (interactive only)
# Guarded so the file can be safely sourced by Shiny apps.
# =========================================================
if (interactive() && exists("DATA_MODE", inherits = FALSE)) {
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
  make_summary_plots(in_csv = in_csv, out_dir = out_dir, mode_label = DATA_MODE, show_tukey = TRUE)
}
