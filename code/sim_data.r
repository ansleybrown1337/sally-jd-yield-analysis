# ---------------------------------------------
# Simulate multi-environment variety trial data
# RCBD + anisotropic exponential spatial errors
# Entry fixed "truth" + random GxE
# ---------------------------------------------

# Minimal deps, base R + stats only
set.seed(20251103)

# ---- Configuration ----
cfg <- list(
  sites       = c("FortCollins","Akron","RockyFord"), # site names
  years       = 2024:2025,                            # years
  n_entry     = 20,                                   # total entries in the program
  reps        = 4,                                    # reps per environment
  nrow        = 16,                                   # rows per environment (grid)
  ncol        = 10,                                   # cols per environment (grid)
  # spatial covariance (anisotropic exponential)
  sigma2_eps  = 500^2,        # spatial field variance
  nugget2     = 250^2,        # independent (plot-level) variance
  range_row   = 3.5,          # range along rows (in grid units)
  range_col   = 6.0,          # range along cols
  # entry main effects (fixed "truth")
  mu_global   = 4500,         # grand mean yield
  entry_sd    = 250,          # SD for entry main-effect deviations
  # GxE (entry × environment) random effects
  gxe_sd      = 200,          # SD for GxE deviations
  # missingness: each env has this fraction of plot yields missing at random
  frac_missing_per_env = 0.20,
  # harvesting only central sub-plot? (kept here for future; grid still full)
  # output control
  make_csv    = TRUE,
  out_csv     = "sim_data/trial_sim.csv"
)

# ---- Helpers ----

# Build anisotropic exponential correlation given distances in row & col
corr_expa <- function(dr, dc, range_row, range_col) {
  # distance metric using axis-specific ranges
  # "product-sum" style: use scaled Euclidean with different ranges
  d_scaled <- sqrt( (dr / range_row)^2 + (dc / range_col)^2 )
  exp(-d_scaled)
}

# Fast generator for Gaussian field via Cholesky (works for moderate grid sizes)
simulate_spatial_field <- function(rows, cols, range_row, range_col, sigma2, nugget2 = 0) {
  # coordinates
  idx <- expand.grid(row = seq_len(rows), col = seq_len(cols))
  n <- nrow(idx)

  # compute pairwise distances efficiently
  # vectorized differences
  dr <- outer(idx$row, idx$row, FUN = function(a,b) abs(a - b))
  dc <- outer(idx$col, idx$col, FUN = function(a,b) abs(a - b))

  R <- corr_expa(dr, dc, range_row, range_col)
  # covariance = sigma^2 * R + nugget^2 * I
  Sigma <- sigma2 * R
  diag(Sigma) <- diag(Sigma) + nugget2

  # Cholesky
  L <- chol(Sigma + 1e-8 * diag(n))  # small jitter for numerics
  z <- rnorm(n)
  y <- as.vector(t(L) %*% z)
  list(field = y, coords = idx)
}

# Assign entries to grid in RCBD: randomization per rep
assign_entries_rcbd <- function(nrow, ncol, reps, entries) {
  plots_per_rep <- (nrow * ncol) / reps
  if ((nrow * ncol) %% reps != 0) stop("Grid must be divisible by reps.")
  if (plots_per_rep %% length(entries) != 0)
    warning("Plots per rep not a multiple of n_entry; some entries will repeat unevenly.")

  layout <- data.frame()
  # partition rows into contiguous stripes per rep for convenience
  rows_per_rep <- nrow / reps
  if (nrow %% reps != 0) stop("nrow must be divisible by reps for striped RCBD.")
  for (r in seq_len(reps)) {
    r_lo <- (r - 1) * rows_per_rep + 1
    r_hi <- r * rows_per_rep
    block <- expand.grid(row = r_lo:r_hi, col = 1:ncol)
    # allocate entries randomly across block
    k <- nrow(block)
    ent_block <- sample(rep(entries, length.out = k))
    block$rep <- r
    block$entry <- ent_block
    layout <- rbind(layout, block)
  }
  layout[order(layout$row, layout$col), ]
}

# Introduce missing harvested plots per environment while preserving layout rows
mark_missing_yields_by_env <- function(df_env, frac_missing) {
  if (frac_missing <= 0) return(df_env)
  n_missing <- floor(nrow(df_env) * frac_missing)
  missing_rows <- if (n_missing > 0) sort(sample(seq_len(nrow(df_env)), n_missing)) else integer(0)
  attr(df_env, "missing_yield_rows") <- missing_rows
  df_env
}

# ---- Simulation driver ----

simulate_trial <- function(cfg) {
  envs <- expand.grid(site = cfg$sites, year = cfg$years, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  envs$env <- paste(envs$site, envs$year, sep = "_")

  # entry main effects (fixed truth)
  entry_ids <- seq_len(cfg$n_entry)
  entry_effect <- rnorm(cfg$n_entry, mean = 0, sd = cfg$entry_sd)
  names(entry_effect) <- entry_ids

  # GxE deviations per env×entry
  gxe <- expand.grid(env = envs$env, entry = entry_ids, KEEP.OUT.ATTRS = FALSE)
  gxe$dev <- rnorm(nrow(gxe), 0, cfg$gxe_sd)

  # container
  out <- vector("list", length = nrow(envs))
  meta_missing <- list()

  for (i in seq_len(nrow(envs))) {
    env <- envs$env[i]

    # spatial field for this env
    sp <- simulate_spatial_field(
      rows = cfg$nrow, cols = cfg$ncol,
      range_row = cfg$range_row, range_col = cfg$range_col,
      sigma2 = cfg$sigma2_eps, nugget2 = cfg$nugget2
    )

    # RCBD layout and entry assignment
    lay <- assign_entries_rcbd(cfg$nrow, cfg$ncol, cfg$reps, entry_ids)

    # potentially mark some harvested plot yields as missing
    lay_env <- mark_missing_yields_by_env(lay, cfg$frac_missing_per_env)
    missing_rows <- attr(lay_env, "missing_yield_rows")
    meta_missing[[env]] <- missing_rows

    # merge spatial field onto layout
    lay_env$plot_id <- seq_len(nrow(lay_env))
    # map spatial value by matching row/col
    idx <- match(paste(lay_env$row, lay_env$col),
                 paste(sp$coords$row, sp$coords$col))
    eps <- sp$field[idx]

    # assemble mean structure
    ee  <- entry_effect[as.character(lay_env$entry)]
    dev <- merge(
      data.frame(env = env, entry = lay_env$entry),
      gxe, by = c("env","entry"), all.x = TRUE, sort = FALSE
    )$dev

    mu <- cfg$mu_global + ee + dev

    # observed yield
    y <- as.numeric(mu + eps)
    if (length(missing_rows) > 0L) {
      y[missing_rows] <- NA_real_
    }

    df <- data.frame(
      site  = envs$site[i],
      year  = envs$year[i],
      env   = env,
      rep   = lay_env$rep,
      row   = lay_env$row,
      col   = lay_env$col,
      entry = lay_env$entry,
      yield = y
    )
    out[[i]] <- df
  }

  trial <- do.call(rbind, out)
  attr(trial, "entry_truth") <- data.frame(entry = entry_ids,
                                           effect = entry_effect)
  attr(trial, "missing_rows_by_env") <- meta_missing
  trial
}

# ---- Run simulation ----
trial <- simulate_trial(cfg)

# Quick sanity checks
cat("Rows:", nrow(trial), " Environments:", length(unique(trial$env)),
    " Entries total:", length(unique(trial$entry)), "\n")

# Inspect a slice
head(trial, 12)

# Optionally save
if (isTRUE(cfg$make_csv)) {
  write.csv(trial, cfg$out_csv, row.names = FALSE)
  message("Wrote: ", cfg$out_csv)
}

# Save truth tables for validation (optional)
truth_entry <- attr(trial, "entry_truth")
write.csv(truth_entry, "sim_data/entry_truth.csv", row.names = FALSE)

drops <- attr(trial, "missing_rows_by_env")
# write a simple log
con <- file("sim_data/env_missing_yields.log", open = "wt")
on.exit(close(con), add = TRUE)
for (nm in names(drops)) {
  cat(nm, ": missing yield row indices -> ",
      if (length(drops[[nm]]) == 0) "none" else paste(drops[[nm]], collapse = ","),
      "\n", file = con)
}

