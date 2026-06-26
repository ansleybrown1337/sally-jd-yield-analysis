# SAS-style spatial anisotropic exponential covariance for Stage 1 fits.
# Target covariance: sigma_e^2 * prod_k exp(-theta_k * |d_k| ^ power_k),
# plus an independent random intercept variance for rep.

expa_power_upper <- 10

expa_power_from_unconstrained <- function(x) {
  expa_power_upper * stats::plogis(x)
}

expa_make_cor <- function(row, col, theta_row, theta_col, power_row, power_col) {
  d_row <- abs(outer(row, row, "-"))
  d_col <- abs(outer(col, col, "-"))
  exp(-theta_row * (d_row ^ power_row) - theta_col * (d_col ^ power_col))
}

expa_solve_chol <- function(chol_v, b) {
  backsolve(chol_v, forwardsolve(t(chol_v), b))
}

expa_reml_components <- function(par, y, x_mat, z_mat, row, col) {
  sigma_e <- exp(par[[1]])
  sigma_rep <- exp(par[[2]])
  theta_row <- exp(par[[3]])
  theta_col <- exp(par[[4]])
  power_row <- expa_power_from_unconstrained(par[[5]])
  power_col <- expa_power_from_unconstrained(par[[6]])

  cor_mat <- expa_make_cor(row, col, theta_row, theta_col, power_row, power_col)
  v_mat <- sigma_e^2 * cor_mat + sigma_rep^2 * tcrossprod(z_mat)
  diag(v_mat) <- diag(v_mat) + 1e-8

  chol_v <- tryCatch(chol(v_mat), error = function(e) NULL)
  if (is.null(chol_v)) return(NULL)

  v_inv_y <- expa_solve_chol(chol_v, y)
  v_inv_x <- expa_solve_chol(chol_v, x_mat)
  xt_v_inv_x <- crossprod(x_mat, v_inv_x)
  chol_xt <- tryCatch(chol(xt_v_inv_x), error = function(e) NULL)
  if (is.null(chol_xt)) return(NULL)

  beta <- expa_solve_chol(chol_xt, crossprod(x_mat, v_inv_y))
  marginal_resid <- as.numeric(y - x_mat %*% beta)
  quad <- as.numeric(crossprod(marginal_resid, expa_solve_chol(chol_v, marginal_resid)))
  logdet_v <- 2 * sum(log(diag(chol_v)))
  logdet_xt <- 2 * sum(log(diag(chol_xt)))
  n <- length(y)
  p <- qr(x_mat)$rank
  reml_loglik <- -0.5 * ((n - p) * log(2 * pi) + logdet_v + logdet_xt + quad)

  rep_blup <- sigma_rep^2 * crossprod(z_mat, expa_solve_chol(chol_v, marginal_resid))
  fitted_marginal <- as.numeric(x_mat %*% beta)
  fitted_conditional <- fitted_marginal + as.numeric(z_mat %*% rep_blup)
  conditional_resid <- as.numeric(y - fitted_conditional)

  list(
    sigma_e = sigma_e,
    sigma_rep = sigma_rep,
    theta_row = theta_row,
    theta_col = theta_col,
    power_row = power_row,
    power_col = power_col,
    cor_mat = cor_mat,
    v_mat = v_mat,
    chol_v = chol_v,
    xt_v_inv_x = xt_v_inv_x,
    beta = as.numeric(beta),
    vcov_beta = chol2inv(chol_xt),
    reml_loglik = reml_loglik,
    marginal_resid = marginal_resid,
    conditional_resid = conditional_resid,
    fitted_marginal = fitted_marginal,
    fitted_conditional = fitted_conditional,
    rep_blup = as.numeric(rep_blup),
    quad = quad
  )
}

expa_neg_reml <- function(par, y, x_mat, z_mat, row, col) {
  comp <- expa_reml_components(par, y, x_mat, z_mat, row, col)
  if (is.null(comp) || !is.finite(comp$reml_loglik)) return(Inf)
  -comp$reml_loglik
}

fit_expa_reml <- function(df_env, alpha = 0.05) {
  df_env <- df_env[!is.na(df_env$yield), , drop = FALSE]
  df_env$entry <- factor(df_env$entry)
  df_env$rep <- factor(df_env$rep)

  y <- as.numeric(df_env$yield)
  x_mat <- stats::model.matrix(~ entry, data = df_env)
  z_mat <- stats::model.matrix(~ 0 + rep, data = df_env)
  row <- as.numeric(df_env$row)
  col <- as.numeric(df_env$col)
  n <- length(y)
  p <- qr(x_mat)$rank

  if (n <= p + 6L) {
    stop("EXPA REML fit has too few observations for fixed effects plus covariance parameters.")
  }

  y_sd <- stats::sd(y)
  if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1

  power_to_unconstrained <- function(pwr) stats::qlogis(pwr / expa_power_upper)
  starts <- expand.grid(
    theta_row = c(0.03, 0.08, 0.15),
    theta_col = c(0.03, 0.08, 0.15),
    power = c(1.0, 2.0)
  )

  lower <- c(log(y_sd * 1e-5), log(y_sd * 1e-5), log(1e-5), log(1e-5), power_to_unconstrained(0.05), power_to_unconstrained(0.05))
  upper <- c(log(y_sd * 20), log(y_sd * 20), log(20), log(20), power_to_unconstrained(expa_power_upper - 1e-4), power_to_unconstrained(expa_power_upper - 1e-4))

  best <- NULL
  for (i in seq_len(nrow(starts))) {
    start <- c(
      log(y_sd * 0.75),
      log(y_sd * 0.25),
      log(starts$theta_row[[i]]),
      log(starts$theta_col[[i]]),
      power_to_unconstrained(starts$power[[i]]),
      power_to_unconstrained(starts$power[[i]])
    )

    opt <- tryCatch(
      stats::optim(
        par = start,
        fn = expa_neg_reml,
        y = y,
        x_mat = x_mat,
        z_mat = z_mat,
        row = row,
        col = col,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = 500)
      ),
      error = function(e) NULL
    )

    if (!is.null(opt) && is.finite(opt$value) && (is.null(best) || opt$value < best$value)) {
      best <- opt
    }
  }

  if (is.null(best)) stop("EXPA REML optimization failed for all starting values.")

  comp <- expa_reml_components(best$par, y, x_mat, z_mat, row, col)
  if (is.null(comp)) stop("EXPA REML fit could not be reconstructed from optimized parameters.")

  k <- p + 6L
  aic <- -2 * comp$reml_loglik + 2 * k
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
  bic <- -2 * comp$reml_loglik + log(n) * k

  entry_levels <- levels(df_env$entry)
  l_mat <- stats::model.matrix(~ entry, data = data.frame(entry = factor(entry_levels, levels = entry_levels)))
  estimates <- as.numeric(l_mat %*% comp$beta)
  ses <- sqrt(pmax(0, diag(l_mat %*% comp$vcov_beta %*% t(l_mat))))
  df_resid <- max(1, n - p)
  tcrit <- stats::qt(1 - alpha / 2, df = df_resid)

  lsm <- data.frame(
    entry = entry_levels,
    estimate = estimates,
    stderr = ses,
    df = df_resid,
    lower.CL = estimates - tcrit * ses,
    upper.CL = estimates + tcrit * ses,
    stringsAsFactors = FALSE
  )

  norm_resid <- as.numeric(expa_solve_chol(comp$chol_v, comp$marginal_resid))

  structure(
    list(
      engine = "custom::sas_sp_expa_reml",
      converged = isTRUE(best$convergence == 0),
      opt_message = if (is.null(best$message)) "" else best$message,
      n = n,
      p = p,
      k = k,
      logLik = comp$reml_loglik,
      AIC = aic,
      AICc = aicc,
      BIC = bic,
      alpha = alpha,
      data = df_env,
      x_mat = x_mat,
      z_mat = z_mat,
      params = c(
        sigma_e = comp$sigma_e,
        sigma_rep = comp$sigma_rep,
        theta_row = comp$theta_row,
        theta_col = comp$theta_col,
        power_row = comp$power_row,
        power_col = comp$power_col
      ),
      beta = comp$beta,
      vcov_beta = comp$vcov_beta,
      fitted = comp$fitted_conditional,
      resid = comp$conditional_resid,
      resid_norm = norm_resid,
      lsm = lsm
    ),
    class = "expa_reml_fit"
  )
}
