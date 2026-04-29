suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(DT)
  library(dplyr)
  library(readr)
  library(tibble)
  library(htmltools)
})

app_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
backend_env <- new.env(parent = globalenv())
source(file.path(app_root, "code", "r_equivalent.R"), local = backend_env)

required_cols <- c("site", "year", "env", "rep", "row", "col", "entry", "yield")
spatial_cov_types <- c("expa", "exp", "sph", "gau")
fallback_output_root <- file.path(tempdir(), "sally-yield-analysis-output")

path_label <- function(path) {
  clean_path <- gsub("\\\\", "/", path)
  clean_root <- gsub("\\\\", "/", app_root)
  sub(paste0("^", clean_root, "/?"), "", clean_path)
}

discover_example_files <- function() {
  roots <- c(file.path(app_root, "real_data"), file.path(app_root, "sim_data"))
  files <- unlist(
    lapply(roots, function(root) {
      if (!dir.exists(root)) return(character(0))
      list.files(root, pattern = "\\.csv$", full.names = TRUE)
    }),
    use.names = FALSE
  )
  sort(unique(files))
}

make_run_output_dir <- function() {
  dir.create(fallback_output_root, recursive = TRUE, showWarnings = FALSE)
  path <- tempfile(pattern = "run-", tmpdir = fallback_output_root)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

fix_ci_cols <- function(df) {
  nm <- names(df)
  nm <- sub("^lower\\.CL\\.$", "lower.CL", nm)
  nm <- sub("^upper\\.CL\\.$", "upper.CL", nm)
  names(df) <- nm
  df
}

fix_group_col <- function(df) {
  if ("group" %in% names(df)) return(df)
  if (".group" %in% names(df)) return(dplyr::rename(df, group = .group))
  dplyr::mutate(df, group = "")
}

clean_entry_labels <- function(x) {
  sub("^entry", "", as.character(x))
}

build_env_summary <- function(trial_df) {
  trial_df %>%
    dplyr::group_by(env) %>%
    dplyr::summarise(
      n_plots = dplyr::n(),
      n_missing = sum(is.na(yield)),
      n_used = n_plots - n_missing,
      n_rep = dplyr::n_distinct(rep),
      n_entry = dplyr::n_distinct(entry),
      .groups = "drop"
    )
}

build_stage2_plot_data <- function(result) {
  lsm <- result$lsm_across %>%
    tibble::as_tibble() %>%
    fix_ci_cols() %>%
    dplyr::mutate(
      entry = clean_entry_labels(entry),
      estimate = as.numeric(estimate),
      stderr = suppressWarnings(as.numeric(stderr)),
      df = suppressWarnings(as.numeric(df)),
      lower.CL = as.numeric(lower.CL),
      upper.CL = as.numeric(upper.CL)
    )

  cld <- result$cld %>%
    tibble::as_tibble() %>%
    fix_group_col() %>%
    dplyr::mutate(
      entry = clean_entry_labels(entry),
      group = trimws(as.character(group))
    ) %>%
    dplyr::select(entry, group)

  lsm %>%
    dplyr::left_join(cld, by = "entry") %>%
    dplyr::mutate(
      group = dplyr::coalesce(group, ""),
      hover = paste0(
        "Entry: ", entry,
        "<br>Adjusted mean: ", round(estimate, 2),
        "<br>95% CI: [", round(lower.CL, 2), ", ", round(upper.CL, 2), "]",
        ifelse(is.na(stderr), "", paste0("<br>SE: ", round(stderr, 3))),
        ifelse(is.na(df), "", paste0("<br>df: ", round(df, 1))),
        ifelse(group == "", "", paste0("<br>CLD: ", group))
      )
    ) %>%
    dplyr::arrange(estimate) %>%
    dplyr::mutate(entry = factor(entry, levels = entry))
}

build_blup_plot_data <- function(result) {
  result$blups %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      entry = clean_entry_labels(entry),
      BLUP = as.numeric(BLUP),
      SE = as.numeric(SE),
      lower = BLUP - 1.96 * SE,
      upper = BLUP + 1.96 * SE,
      hover = paste0(
        "Entry: ", entry,
        "<br>BLUP: ", round(BLUP, 3),
        ifelse(is.na(SE), "", paste0("<br>SE: ", round(SE, 3))),
        ifelse(is.na(SE), "", paste0("<br>Approx. 95% interval: [", round(lower, 3), ", ", round(upper, 3), "]"))
      )
    ) %>%
    dplyr::arrange(BLUP) %>%
    dplyr::mutate(entry = factor(entry, levels = entry))
}

make_metric_card <- function(title, value, subtitle) {
  bslib::card(
    class = "metric-card",
    bslib::card_header(title),
    div(class = "metric-value", value),
    div(class = "metric-subtitle", subtitle)
  )
}

make_sidebar_stat <- function(title, value, subtitle) {
  tags$div(
    class = "sidebar-stat",
    tags$div(class = "sidebar-stat-title", title),
    tags$div(class = "sidebar-stat-value", value),
    tags$div(class = "sidebar-stat-subtitle", subtitle)
  )
}

adjusted_plot_height <- function(n_entries) {
  paste0(max(700, min(1500, 28 * n_entries + 220)), "px")
}

diagnostic_plot_height <- function(n_rows) {
  paste0(max(650, min(1100, 24 * n_rows + 180)), "px")
}

build_status_summary <- function(result) {
  n_env <- dplyr::n_distinct(result$stage1$env)
  best_cov <- result$best_cov
  n_spatial <- sum(best_cov$best_cov %in% spatial_cov_types, na.rm = TRUE)
  n_fallback <- sum(!best_cov$best_cov %in% spatial_cov_types, na.rm = TRUE)
  outlier_count <- sum(result$diagnostics_stage1$flag_outlier %in% TRUE, na.rm = TRUE)
  stage2_label <- if (n_env > 1L) {
    paste0("Across-environment meta-analysis was run because ", n_env, " environments contributed Stage 1 means.")
  } else {
    "Only one environment contributed Stage 1 means, so the app used the single-environment Stage 2 fallback."
  }

  env_bits <- if (nrow(best_cov)) {
    apply(best_cov, 1, function(row) {
      paste0(row[["env"]], " -> ", row[["best_cov"]])
    })
  } else {
    "<no environment fits recorded>"
  }

  tags$div(
    class = "status-summary",
    tags$div(
      class = "status-summary-title",
      "Detected model workflow"
    ),
    tags$ul(
      class = "status-summary-list",
      tags$li(paste0("Detected ", n_env, " environment(s).")),
      tags$li(paste0(n_spatial, " environment(s) used a spatial covariance model; ", n_fallback, " used a fallback model.")),
      tags$li(stage2_label),
      tags$li(paste0("Flagged ", outlier_count, " plot-level outlier(s) using |normalized residual| > 3.")),
      tags$li(paste0("Best Stage 1 fit by environment: ", paste(env_bits, collapse = "; "), "."))
    )
  )
}

make_across_plot <- function(plot_df) {
  if (!nrow(plot_df)) {
    return(plotly::plotly_empty(type = "scatter", mode = "markers"))
  }

  p <- plotly::plot_ly(
    data = plot_df,
    x = ~estimate,
    y = ~entry,
    type = "scatter",
    mode = "markers",
    marker = list(size = 10, color = "#9b5d2e"),
    error_x = list(
      type = "data",
      symmetric = FALSE,
      array = plot_df$upper.CL - plot_df$estimate,
      arrayminus = plot_df$estimate - plot_df$lower.CL,
      color = "#9b5d2e",
      thickness = 1.4,
      width = 0
    ),
    text = ~hover,
    hovertemplate = "%{text}<extra></extra>",
    showlegend = FALSE
  )

  if (any(nzchar(plot_df$group))) {
    label_df <- dplyr::filter(plot_df, nzchar(group))
    p <- p %>%
      plotly::add_trace(
        data = label_df,
        x = label_df$upper.CL,
        y = label_df$entry,
        text = label_df$group,
        type = "scatter",
        mode = "text",
        textposition = "middle right",
        textfont = list(color = "#1f5d50", size = 12),
        hoverinfo = "skip",
        showlegend = FALSE,
        inherit = FALSE
      )
  }

  p %>%
    plotly::layout(
      title = list(text = "Across-environment adjusted means"),
      xaxis = list(title = "Adjusted mean yield", zeroline = FALSE),
      yaxis = list(title = "Entry", automargin = TRUE),
      margin = list(l = 110, r = 80, t = 60, b = 60),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)"
    )
}

make_blup_plot <- function(plot_df) {
  if (!nrow(plot_df)) {
    return(plotly::plotly_empty(type = "scatter", mode = "markers"))
  }

  plotly::plot_ly(
    data = plot_df,
    x = ~BLUP,
    y = ~entry,
    type = "scatter",
    mode = "markers",
    marker = list(size = 10, color = "#1f5d50"),
    error_x = list(
      type = "data",
      symmetric = FALSE,
      array = plot_df$upper - plot_df$BLUP,
      arrayminus = plot_df$BLUP - plot_df$lower,
      color = "#1f5d50",
      thickness = 1.4,
      width = 0
    ),
    text = ~hover,
    hovertemplate = "%{text}<extra></extra>",
    showlegend = FALSE
  ) %>%
    plotly::layout(
      title = list(text = "Entry BLUPs"),
      xaxis = list(title = "BLUP", zeroline = TRUE, zerolinecolor = "#7d8590"),
      yaxis = list(title = "Entry", automargin = TRUE),
      margin = list(l = 110, r = 50, t = 60, b = 60),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)"
    )
}

make_heatmap_plot <- function(df_env, fill_var) {
  if (!nrow(df_env)) {
    return(plotly::plotly_empty(type = "heatmap"))
  }

  label_map <- c(
    resid_norm = "Normalized residual",
    resid = "Raw residual",
    yield = "Yield"
  )
  fill_label <- unname(label_map[[fill_var]])
  colorscale <- if (identical(fill_var, "yield")) {
    "Viridis"
  } else {
    list(
      list(0.00, "#8c1d18"),
      list(0.50, "#f7f2e8"),
      list(1.00, "#0f766e")
    )
  }

  df_env <- df_env %>%
    dplyr::mutate(
      plot_value = .data[[fill_var]],
      hover = paste0(
        "Env: ", env,
        "<br>Site: ", site,
        "<br>Year: ", year,
        "<br>Rep: ", rep,
        "<br>Entry: ", entry,
        "<br>Row: ", row,
        "<br>Col: ", col,
        "<br>", fill_label, ": ", round(plot_value, 3)
      )
    )

  row_vals <- sort(unique(df_env$row))
  col_vals <- sort(unique(df_env$col))
  z_mat <- matrix(NA_real_, nrow = length(row_vals), ncol = length(col_vals))
  text_mat <- matrix("", nrow = length(row_vals), ncol = length(col_vals))

  for (i in seq_len(nrow(df_env))) {
    r_idx <- match(df_env$row[i], row_vals)
    c_idx <- match(df_env$col[i], col_vals)
    z_mat[r_idx, c_idx] <- df_env$plot_value[i]
    text_mat[r_idx, c_idx] <- df_env$hover[i]
  }

  plotly::plot_ly(
    x = col_vals,
    y = row_vals,
    z = z_mat,
    text = text_mat,
    type = "heatmap",
    colorscale = colorscale,
    hovertemplate = "%{text}<extra></extra>",
    showscale = TRUE
  ) %>%
    plotly::layout(
      title = list(text = paste0(fill_label, " by field position")),
      xaxis = list(title = "Column", dtick = 1),
      yaxis = list(title = "Row", autorange = "reversed", dtick = 1),
      margin = list(l = 70, r = 20, t = 60, b = 60),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)"
    )
}

make_resid_fit_plot <- function(df_env) {
  if (!nrow(df_env)) {
    return(plotly::plotly_empty(type = "scatter", mode = "markers"))
  }

  x_rng <- range(df_env$fitted, na.rm = TRUE)

  plotly::plot_ly(
    data = df_env,
    x = ~fitted,
    y = ~resid_norm,
    type = "scatter",
    mode = "markers",
    color = ~flag_outlier,
    colors = c("FALSE" = "#1f5d50", "TRUE" = "#b42318"),
    text = ~paste0(
      "Entry: ", entry,
      "<br>Rep: ", rep,
      "<br>Fitted: ", round(fitted, 2),
      "<br>Residual: ", round(resid_norm, 3)
    ),
    hovertemplate = "%{text}<extra></extra>"
  ) %>%
    plotly::layout(
      title = list(text = "Residuals vs fitted"),
      xaxis = list(title = "Fitted value"),
      yaxis = list(title = "Normalized residual"),
      shapes = list(
        list(
          type = "line",
          x0 = x_rng[1],
          x1 = x_rng[2],
          y0 = 0,
          y1 = 0,
          line = list(color = "#7d8590", dash = "dash")
        )
      ),
      legend = list(title = list(text = "Flagged outlier")),
      margin = list(l = 60, r = 20, t = 60, b = 55),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)"
    )
}

make_hist_plot <- function(df_env) {
  if (!nrow(df_env)) {
    return(plotly::plotly_empty(type = "histogram"))
  }

  plotly::plot_ly(
    x = df_env$resid_norm,
    type = "histogram",
    nbinsx = 25,
    marker = list(color = "#2f6d62")
  ) %>%
    plotly::layout(
      title = list(text = "Residual distribution"),
      xaxis = list(title = "Normalized residual"),
      yaxis = list(title = "Count"),
      margin = list(l = 60, r = 20, t = 60, b = 55),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)"
    )
}

make_qq_plot <- function(df_env) {
  vals <- df_env$resid_norm
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(plotly::plotly_empty(type = "scatter", mode = "markers"))
  }

  qq <- stats::qqnorm(vals, plot.it = FALSE)
  q_sample <- stats::quantile(vals, probs = c(0.25, 0.75), na.rm = TRUE)
  q_theory <- stats::qnorm(c(0.25, 0.75))
  slope <- diff(q_sample) / diff(q_theory)
  intercept <- q_sample[1] - slope * q_theory[1]
  x_line <- range(qq$x, na.rm = TRUE)

  plotly::plot_ly(
    x = qq$x,
    y = qq$y,
    type = "scatter",
    mode = "markers",
    marker = list(color = "#9b5d2e", size = 8),
    hovertemplate = paste0(
      "Theoretical quantile: %{x:.3f}",
      "<br>Sample quantile: %{y:.3f}<extra></extra>"
    ),
    showlegend = FALSE
  ) %>%
    plotly::add_lines(
      x = x_line,
      y = intercept + slope * x_line,
      line = list(color = "#1f5d50", width = 2),
      hoverinfo = "skip",
      showlegend = FALSE
    ) %>%
    plotly::layout(
      title = list(text = "Residual Q-Q plot"),
      xaxis = list(title = "Theoretical quantile"),
      yaxis = list(title = "Sample quantile"),
      margin = list(l = 60, r = 20, t = 60, b = 55),
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)"
    )
}

datatable_opts <- list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE)
example_files <- discover_example_files()
example_choices <- if (length(example_files)) {
  stats::setNames(example_files, vapply(example_files, path_label, character(1)))
} else {
  c("No bundled CSV files found" = "")
}
example_selected <- if (length(example_files)) example_files[[1]] else ""

theme <- bslib::bs_theme(
  version = 5,
  bg = "#f5f1e8",
  fg = "#1f2933",
  primary = "#1f5d50",
  secondary = "#9b5d2e",
  base_font = bslib::font_google("IBM Plex Sans"),
  heading_font = bslib::font_google("Newsreader"),
  code_font = bslib::font_google("IBM Plex Mono")
)

ui <- bslib::page_sidebar(
  title = div(
    class = "app-title",
    span("Variety Trial Explorer")
  ),
  theme = theme,
  sidebar = bslib::sidebar(
    width = 340,
    open = "desktop",
    radioButtons(
      "data_source",
      "Data source",
      choices = c("Bundled example" = "bundled", "Upload CSV" = "upload"),
      selected = "upload"
    ),
    conditionalPanel(
      condition = "input.data_source === 'bundled'",
      selectInput("example_file", "Bundled file", choices = example_choices, selected = example_selected)
    ),
    conditionalPanel(
      condition = "input.data_source === 'upload'",
      fileInput("trial_file", "Upload trial CSV", accept = ".csv")
    ),
    numericInput("alpha", "Alpha", value = 0.05, min = 0.001, max = 0.50, step = 0.01),
    actionButton("run_analysis", "Run analysis", class = "btn-primary"),
    downloadButton("download_outputs", "Download outputs (.zip)", class = "btn-primary"),
    tags$div(style = "height: 0.5rem;"),
    downloadButton("download_report", "Download run report", class = "btn-outline-secondary"),
    helpText("On Posit Connect, outputs are prepared in a temporary server-side folder and downloaded through your browser."),
    tags$hr(),
    uiOutput("sidebar_metrics"),
    tags$hr(),
    tags$div(
      class = "sidebar-stat",
      tags$div(class = "sidebar-stat-title", "Author"),
      tags$div(class = "sidebar-stat-subtitle", "A.J. Brown"),
      tags$div(class = "sidebar-stat-subtitle", "Agricultural Data Scientist"),
      tags$a(
        href = "https://sites.google.com/view/ansleyjbrown",
        target = "_blank",
        rel = "noopener noreferrer",
        "Website"
      )
    )
  ),
  tags$head(
    tags$style(HTML("
      body {
        background:
          radial-gradient(circle at top left, rgba(31, 93, 80, 0.10), transparent 28%),
          linear-gradient(180deg, #f5f1e8 0%, #eef4ef 100%);
      }
      .app-title {
        font-size: 1.35rem;
        font-weight: 700;
        letter-spacing: 0.02em;
      }
      .hero-banner {
        padding: 1.2rem 1.35rem;
        margin-bottom: 1rem;
        border-radius: 18px;
        background: rgba(255, 255, 255, 0.72);
        border: 1px solid rgba(114, 98, 75, 0.15);
        box-shadow: 0 12px 28px rgba(38, 52, 61, 0.08);
      }
      .hero-banner h2 {
        margin: 0 0 0.3rem 0;
        font-size: 1.65rem;
      }
      .hero-banner p {
        margin: 0;
        color: #52606d;
        max-width: 60rem;
      }
      .metric-card,
      .bslib-card {
        border-radius: 18px;
        border: 1px solid rgba(114, 98, 75, 0.14);
        box-shadow: 0 14px 30px rgba(38, 52, 61, 0.07);
      }
      .metric-value {
        font-size: 2rem;
        font-weight: 700;
        color: #1f5d50;
        line-height: 1.1;
      }
      .metric-subtitle {
        margin-top: 0.3rem;
        color: #5b646d;
      }
      .sidebar-stat {
        padding: 0.8rem 0.9rem;
        margin-bottom: 0.75rem;
        border-radius: 14px;
        background: rgba(255, 255, 255, 0.72);
        border: 1px solid rgba(114, 98, 75, 0.14);
        box-shadow: 0 8px 18px rgba(38, 52, 61, 0.06);
      }
      .sidebar-stat-title {
        font-size: 0.9rem;
        font-weight: 700;
        color: #415161;
        margin-bottom: 0.2rem;
      }
      .sidebar-stat-value {
        font-size: 1.7rem;
        font-weight: 700;
        color: #1f5d50;
        line-height: 1.05;
      }
      .sidebar-stat-subtitle {
        margin-top: 0.2rem;
        color: #5b646d;
        font-size: 0.9rem;
      }
      .status-summary {
        padding: 1rem 1.2rem;
        margin: 0.25rem 0 1rem 0;
        border-radius: 18px;
        background: rgba(255, 255, 255, 0.74);
        border: 1px solid rgba(114, 98, 75, 0.16);
        box-shadow: 0 10px 24px rgba(38, 52, 61, 0.06);
      }
      .status-summary-title {
        font-size: 1.1rem;
        font-weight: 700;
        margin-bottom: 0.55rem;
        color: #1f5d50;
      }
      .status-summary-list {
        margin-bottom: 0;
        padding-left: 1.25rem;
        color: #415161;
      }
      .status-summary-list li {
        margin-bottom: 0.35rem;
      }
      .btn-primary {
        width: 100%;
        font-weight: 700;
        background-color: #1f5d50;
        border-color: #1f5d50;
      }
      .nav-link.active {
        font-weight: 700;
      }
      .form-label, .control-label {
        font-weight: 700;
      }
    "))
  ),
  div(
    class = "hero-banner",
    h2("Interactive front end for the existing spatial workflow"),
    p("Run the current Stage 1 and Stage 2 models, inspect adjusted means with Plotly, and switch site-years in diagnostics without re-rendering a static HTML report.")
  ),
  bslib::navset_card_tab(
    height = "auto",
    bslib::nav_panel(
      "Adjusted means",
      tags$div(
        class = "status-summary",
        tags$p("This tab reports adjusted LS-means and significance groupings from the LS-means workflow. It is separate from the one-stage BLUP analysis shown later.")
      ),
      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Adjusted means"),
          uiOutput("across_plot_ui")
        ),
        bslib::card(
          bslib::card_header("Adjusted means results"),
          selectInput(
            "adjusted_results_view",
            "Results panel",
            choices = c(
              "LS-means table" = "lsm",
              "Significance groups" = "cld"
            ),
            selected = "lsm"
          ),
          uiOutput("adjusted_results_ui")
        )
      )
    ),
    bslib::nav_panel(
      "Diagnostics",
      tags$div(
        class = "status-summary",
        tags$p("These diagnostics are for the Stage 1 per-environment spatial fits used to build the LS-means workflow. They are not BLUP-specific diagnostics.")
      ),
      bslib::layout_sidebar(
        sidebar = bslib::sidebar(
          width = 280,
          title = "Diagnostic controls",
          uiOutput("env_control"),
          selectInput(
            "heatmap_fill",
            "Heatmap fill",
            choices = c(
              "Normalized residual" = "resid_norm",
              "Raw residual" = "resid",
              "Yield" = "yield"
            ),
            selected = "resid_norm"
          ),
          selectInput(
            "diagnostic_plot_type",
            "Diagnostic plot",
            choices = c(
              "Residuals vs fitted" = "resid_fit",
              "Q-Q plot" = "qq",
              "Residual distribution" = "hist"
            ),
            selected = "resid_fit"
          ),
          helpText("These controls only affect the diagnostics view.")
        ),
        bslib::layout_column_wrap(
          width = 1/2,
          bslib::card(
            full_screen = TRUE,
            bslib::card_header("Field map"),
            uiOutput("heatmap_plot_ui")
          ),
          bslib::card(
            full_screen = TRUE,
            bslib::card_header(uiOutput("diagnostic_plot_title")),
            uiOutput("diagnostic_detail_plot_ui")
          )
        ),
        bslib::card(
          bslib::card_header("Diagnostic rows"),
          DT::DTOutput("diagnostic_rows_table")
        )
      )
    ),
    bslib::nav_panel(
      "Outliers",
      tags$div(
        class = "status-summary",
        tags$p("These flagged outliers come from the Stage 1 residual diagnostics used in the LS-means workflow. They are diagnostic flags, not automatic exclusions, and they are not computed from the BLUP fit.")
      ),
      bslib::card(
        bslib::card_header("Flagged outliers across all environments"),
        DT::DTOutput("outlier_table")
      ),
      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          bslib::card_header("Selected outlier row"),
          DT::DTOutput("outlier_detail_table")
        ),
        bslib::card(
          bslib::card_header("Full input row for the selected outlier"),
          DT::DTOutput("outlier_raw_row_table")
        )
      )
    ),
    bslib::nav_panel(
      "Model fit",
      tags$div(
        class = "status-summary",
        tags$p("This tab summarizes the Stage 1 spatial covariance selection and related LS-means workflow fit summaries by environment. It does not describe the one-stage BLUP model.")
      ),
      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          bslib::card_header("Model specifications by environment"),
          DT::DTOutput("model_specs_table")
        ),
        bslib::card(
          bslib::card_header("Stage 1 covariance ranking"),
          DT::DTOutput("model_selection_table")
        )
      ),
      bslib::card(
        bslib::card_header("CV summary by environment"),
        DT::DTOutput("cv_table")
      )
    ),
    bslib::nav_panel(
      "BLUPs",
      tags$div(
        class = "status-summary",
        tags$p("This tab reports entry BLUPs from the one-stage mixed model. BLUPs are shrunken random-effect estimates and should be interpreted separately from the adjusted LS-means in the earlier tabs."),
        tags$p("A separate BLUP diagnostics tab is not shown here because the current backend exports BLUP estimates and standard errors, but not a dedicated BLUP diagnostic dataset.")
      ),
      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          full_screen = TRUE,
          bslib::card_header("Entry BLUPs"),
          uiOutput("blup_plot_ui")
        ),
        bslib::card(
          bslib::card_header("BLUP results"),
          DT::DTOutput("blup_table")
        )
      )
    ),
    bslib::nav_panel(
      "Run summary",
      bslib::card(
        bslib::card_header("Detected model workflow"),
        uiOutput("status_summary")
      ),
      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          bslib::card_header("Environment summary"),
          DT::DTOutput("env_summary_table")
        ),
        bslib::card(
          bslib::card_header("Input preview"),
          DT::DTOutput("trial_head_table")
        )
      ),
      bslib::card(
        bslib::card_header("Backend run report"),
        verbatimTextOutput("run_report", placeholder = TRUE)
      )
    ),
    bslib::nav_panel(
      "About this tool",
      bslib::card(
        bslib::card_header("Overview"),
        tags$p("This Shiny app is an interactive front end for a spatial mixed-model workflow for crop variety trial analysis. It combines per-environment spatial fitting, across-environment adjusted means, and one-stage BLUP estimation in one place.")
      ),
      bslib::card(
        bslib::card_header("How to use"),
        tags$ol(
          tags$li("Upload a CSV or choose a bundled example."),
          tags$li("Click ", tags$code("Run analysis"), " to fit the workflow."),
          tags$li("Use the ", tags$strong("Adjusted means"), ", ", tags$strong("BLUPs"), ", ", tags$strong("Diagnostics"), ", and ", tags$strong("Outliers"), " tabs to inspect results."),
          tags$li("Download the generated outputs as a ZIP file from the sidebar.")
        )
      ),
      bslib::card(
        bslib::card_header("Important assumptions"),
        tags$ul(
          tags$li("Input data must include columns: ", tags$code("site, year, env, rep, row, col, entry, yield"), "."),
          tags$li("Rows with missing ", tags$code("yield"), " values are excluded from model fitting."),
          tags$li("Stage 1 fits spatial covariance models where possible and falls back to RCBD or fixed-effects ANOVA if needed."),
          tags$li("Across-environment LS-means and BLUPs answer different questions and should be interpreted separately."),
          tags$li("Flagged outliers are diagnostic cues based on large normalized residuals, not automatic deletions.")
        )
      ),
      bslib::card(
        bslib::card_header("Authorship"),
        tags$p(
          "A.J. Brown, Agricultural Data Scientist. More information: ",
          tags$a(
            href = "https://sites.google.com/view/ansleyjbrown",
            target = "_blank",
            rel = "noopener noreferrer",
            "https://sites.google.com/view/ansleyjbrown"
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  analysis_data <- reactiveVal(NULL)

  chosen_input_path <- reactive({
    if (identical(input$data_source, "upload")) {
      req(input$trial_file$datapath)
      return(input$trial_file$datapath)
    }
    req(input$example_file)
    input$example_file
  })

  observeEvent(input$run_analysis, {
    tryCatch({
      in_csv <- chosen_input_path()
      trial <- read.csv(in_csv, stringsAsFactors = FALSE)
      missing_cols <- setdiff(required_cols, names(trial))
      if (length(missing_cols)) {
        stop(
          "Input file is missing required columns: ",
          paste(missing_cols, collapse = ", ")
        )
      }

      out_dir <- make_run_output_dir()

      shiny::withProgress(message = "Running spatial workflow", value = 0.2, {
        result <- backend_env$analyze_trial(
          in_csv = in_csv,
          out_dir = out_dir,
          alpha = input$alpha
        )
        shiny::incProgress(0.8)
        analysis_data(list(
          result = result,
          trial = trial,
          input_path = in_csv,
          out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE)
        ))
      })

      showNotification("Analysis complete.", type = "message", duration = 4)
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })

  output$download_outputs <- downloadHandler(
    filename = function() {
      paste0("sally-yield-analysis-outputs-", Sys.Date(), ".zip")
    },
    content = function(file) {
      state <- analysis_data()
      req(state)
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(state$out_dir)
      files <- list.files(".", all.files = FALSE, no.. = TRUE)
      utils::zip(zipfile = file, files = files)
    }
  )

  output$download_report <- downloadHandler(
    filename = function() {
      paste0("sally-yield-analysis-run-report-", Sys.Date(), ".txt")
    },
    content = function(file) {
      state <- analysis_data()
      req(state)
      file.copy(state$result$report_path, file, overwrite = TRUE)
    }
  )

  current_env_data <- reactive({
    state <- analysis_data()
    req(state)
    diag_df <- state$result$diagnostics_stage1
    req(nrow(diag_df) > 0)
    req(input$env_focus)

    diag_df %>%
      dplyr::mutate(
        row = as.numeric(row),
        col = as.numeric(col),
        fitted = as.numeric(fitted),
        resid = as.numeric(resid),
        resid_norm = as.numeric(resid_norm),
        flag_outlier = as.logical(flag_outlier)
      ) %>%
      dplyr::filter(env == input$env_focus) %>%
      dplyr::filter(!is.na(row), !is.na(col))
  })

  current_outliers <- reactive({
    state <- analysis_data()
    req(state)
    state$result$diagnostics_stage1 %>%
      dplyr::mutate(
        row = as.numeric(row),
        col = as.numeric(col),
        fitted = as.numeric(fitted),
        resid = as.numeric(resid),
        resid_norm = as.numeric(resid_norm),
        flag_outlier = as.logical(flag_outlier)
      ) %>%
      dplyr::filter(flag_outlier %in% TRUE) %>%
      dplyr::arrange(env, desc(abs(resid_norm)), row, col)
  })

  selected_outlier <- reactive({
    rows <- input$outlier_table_rows_selected
    outliers <- current_outliers()
    req(length(rows) == 1, nrow(outliers) >= rows[1])
    outliers[rows[1], , drop = FALSE]
  })

  output$env_control <- renderUI({
    state <- analysis_data()
    if (is.null(state)) {
      return(helpText("Run an analysis to activate environment-level diagnostics."))
    }

    envs <- sort(unique(state$result$diagnostics_stage1$env))
    if (!length(envs)) {
      return(helpText("No environment diagnostics were produced for this run."))
    }

    selectInput("env_focus", "Diagnostic environment", choices = envs, selected = envs[[1]])
  })

  output$status_summary <- renderUI({
    state <- analysis_data()
    if (is.null(state)) {
      return(
        tags$p("Run an analysis to see a human-readable summary of the environment count, Stage 1 spatial fits, Stage 2 mode, and flagged outliers.")
      )
    }

    build_status_summary(state$result)
  })

  output$sidebar_metrics <- renderUI({
    state <- analysis_data()
    if (is.null(state)) {
      return(
        tagList(
          make_sidebar_stat("Environments", "0", "Run the workflow to populate results"),
          make_sidebar_stat("Spatial Fits", "0", "No Stage 1 model selection yet"),
          make_sidebar_stat("Flagged Outliers", "0", "Diagnostics will appear after the first run")
        )
      )
    }

    result <- state$result
    best_cov <- result$best_cov
    n_env <- dplyr::n_distinct(result$stage1$env)
    n_spatial <- sum(best_cov$best_cov %in% spatial_cov_types, na.rm = TRUE)
    n_fallback <- nrow(best_cov) - n_spatial
    outlier_count <- sum(result$diagnostics_stage1$flag_outlier %in% TRUE, na.rm = TRUE)
    stage2_mode <- if (n_env > 1L) "Meta-analysis across environments" else "Single-environment fallback"

    tagList(
      make_sidebar_stat("Environments", n_env, "Site-years with usable Stage 1 means"),
      make_sidebar_stat("Spatial Fits", n_spatial, paste0(n_fallback, " environments used fallback models")),
      make_sidebar_stat("Flagged Outliers", outlier_count, stage2_mode)
    )
  })

  output$across_plot <- plotly::renderPlotly({
    state <- analysis_data()
    req(state)
    make_across_plot(build_stage2_plot_data(state$result))
  })

  output$adjusted_results_ui <- renderUI({
    req(input$adjusted_results_view)
    switch(
      input$adjusted_results_view,
      lsm = DT::DTOutput("lsm_table"),
      cld = DT::DTOutput("cld_table")
    )
  })

  output$across_plot_ui <- renderUI({
    state <- analysis_data()
    plot_df <- if (is.null(state)) tibble::tibble() else build_stage2_plot_data(state$result)
    plotly::plotlyOutput("across_plot", height = adjusted_plot_height(nrow(plot_df)))
  })

  output$heatmap_plot <- plotly::renderPlotly({
    make_heatmap_plot(current_env_data(), input$heatmap_fill)
  })

  output$heatmap_plot_ui <- renderUI({
    plotly::plotlyOutput("heatmap_plot", height = diagnostic_plot_height(nrow(current_env_data())))
  })

  output$diagnostic_plot_title <- renderUI({
    title <- switch(
      input$diagnostic_plot_type,
      resid_fit = "Residuals vs fitted",
      qq = "Q-Q plot",
      hist = "Residual distribution",
      "Diagnostic plot"
    )
    HTML(title)
  })

  output$diagnostic_detail_plot_ui <- renderUI({
    plotly::plotlyOutput("diagnostic_detail_plot", height = "620px")
  })

  output$diagnostic_detail_plot <- plotly::renderPlotly({
    req(input$diagnostic_plot_type)
    switch(
      input$diagnostic_plot_type,
      resid_fit = make_resid_fit_plot(current_env_data()),
      qq = make_qq_plot(current_env_data()),
      hist = make_hist_plot(current_env_data())
    )
  })

  output$blup_plot_ui <- renderUI({
    state <- analysis_data()
    plot_df <- if (is.null(state)) tibble::tibble() else build_blup_plot_data(state$result)
    plotly::plotlyOutput("blup_plot", height = adjusted_plot_height(nrow(plot_df)))
  })

  output$blup_plot <- plotly::renderPlotly({
    state <- analysis_data()
    req(state)
    make_blup_plot(build_blup_plot_data(state$result))
  })

  output$lsm_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(
      fix_ci_cols(state$result$lsm_across),
      options = datatable_opts,
      rownames = FALSE
    )
  })

  output$cld_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(
      fix_group_col(state$result$cld),
      options = datatable_opts,
      rownames = FALSE
    )
  })

  output$blup_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(
      state$result$blups,
      options = c(datatable_opts, list(pageLength = 15, scrollY = "620px")),
      rownames = FALSE
    )
  })

  output$outlier_table <- DT::renderDT({
    df <- current_outliers() %>%
      dplyr::select(env, site, year, rep, row, col, entry, yield, fitted, resid, resid_norm, resid_kind)

    DT::datatable(
      df,
      options = c(datatable_opts, list(pageLength = 15, scrollY = "420px")),
      rownames = FALSE,
      selection = "single",
      filter = "top"
    )
  })

  output$outlier_detail_table <- DT::renderDT({
    outlier <- selected_outlier() %>%
      dplyr::select(env, site, year, rep, row, col, entry, yield, fitted, resid, resid_norm, resid_kind)

    DT::datatable(outlier, options = list(dom = "t"), rownames = FALSE)
  })

  output$outlier_raw_row_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    outlier <- selected_outlier()

    raw_row <- state$trial %>%
      dplyr::mutate(
        env = as.character(env),
        rep = as.character(rep),
        row = as.numeric(row),
        col = as.numeric(col),
        year = as.numeric(year),
        entry = as.character(entry)
      ) %>%
      dplyr::filter(
        env == as.character(outlier$env[[1]]),
        rep == as.character(outlier$rep[[1]]),
        row == outlier$row[[1]],
        col == outlier$col[[1]],
        entry == as.character(outlier$entry[[1]])
      ) %>%
      dplyr::slice_head(n = 1)

    DT::datatable(raw_row, options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  output$diagnostic_rows_table <- DT::renderDT({
    df <- current_env_data() %>%
      dplyr::select(env, site, year, rep, row, col, entry, yield, fitted, resid, resid_norm, resid_kind, flag_outlier)

    DT::datatable(
      df,
      options = c(datatable_opts, list(pageLength = 20, scrollY = "500px")),
      rownames = FALSE,
      filter = "top"
    )
  })

  output$model_specs_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(state$result$model_specs_stage1, options = datatable_opts, rownames = FALSE)
  })

  output$model_selection_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(state$result$stage1_model_selection, options = datatable_opts, rownames = FALSE)
  })

  output$cv_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(state$result$cv_by_env, options = datatable_opts, rownames = FALSE)
  })

  output$env_summary_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(build_env_summary(state$trial), options = datatable_opts, rownames = FALSE)
  })

  output$trial_head_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(utils::head(state$trial, 10), options = datatable_opts, rownames = FALSE)
  })

  output$run_report <- renderText({
    state <- analysis_data()
    req(state)
    report_path <- state$result$report_path
    paste(readLines(report_path), collapse = "\n")
  })
}

shinyApp(ui, server)
