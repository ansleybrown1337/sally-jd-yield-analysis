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
interpretation_md_path <- file.path(app_root, "output_interpretation.md")
bundled_example_path <- file.path(app_root, "sim_data", "trial_sim.csv")

required_cols <- c("site", "year", "env", "rep", "row", "col", "entry", "yield")
spatial_cov_types <- c("expa", "exp", "sph", "gau")
fallback_output_root <- file.path(tempdir(), "sally-yield-analysis-output")

render_markdown_file_html <- function(path) {
  if (!file.exists(path)) {
    return(tags$p("The requested Markdown file was not found."))
  }

  md_text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  HTML(commonmark::markdown_html(md_text))
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

package_version_or_na <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return("not installed")
  as.character(utils::packageVersion(pkg))
}

build_error_log_text <- function(stage, error, context = list()) {
  context_lines <- if (length(context)) {
    paste0("  ", names(context), ": ", unlist(context, use.names = FALSE))
  } else {
    "  <none>"
  }

  call_text <- conditionCall(error)
  call_text <- if (is.null(call_text)) "<not available>" else paste(deparse(call_text), collapse = "\n")

  paste(
    "Sally JD Yield Analysis - error log",
    paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Stage: ", stage),
    paste0("Error class: ", paste(class(error), collapse = ", ")),
    paste0("Message: ", conditionMessage(error)),
    "Call:",
    call_text,
    "Context:",
    paste(context_lines, collapse = "\n"),
    "Package versions:",
    paste0("  R: ", as.character(getRversion())),
    paste0("  shiny: ", package_version_or_na("shiny")),
    paste0("  bslib: ", package_version_or_na("bslib")),
    paste0("  lme4: ", package_version_or_na("lme4")),
    paste0("  reformulas: ", package_version_or_na("reformulas")),
    paste0("  nlme: ", package_version_or_na("nlme")),
    paste0("  emmeans: ", package_version_or_na("emmeans")),
    sep = "\n"
  )
}

make_check_metric <- function(title, value, subtitle) {
  tags$div(
    class = "check-metric",
    tags$div(class = "check-metric-title", title),
    tags$div(class = "check-metric-value", value),
    tags$div(class = "check-metric-subtitle", subtitle)
  )
}

build_data_check_summary <- function(check) {
  if (is.null(check)) {
    return(tags$p("Upload a CSV to run the data check before analysis."))
  }

  report <- check$report
  extra_note <- if (report$extra_cols_ignored > 0L) {
    paste0(
      report$extra_cols_ignored,
      " extra column(s) ignored; ",
      report$extra_blank_cols_ignored,
      " were fully blank."
    )
  } else {
    "No extra columns detected."
  }

  tagList(
    tags$div(
      class = "workflow-steps",
      tags$div(class = "workflow-step workflow-step-complete", "1. Upload"),
      tags$div(class = "workflow-step workflow-step-active", "2. Check data"),
      tags$div(class = "workflow-step", "3. Run analysis")
    ),
    tags$div(
      class = "next-action",
      tags$strong("Next step: "),
      "review the cleaning notes below, then click ",
      tags$code("Run analysis"),
      " in the sidebar."
    ),
    tags$div(
      class = "check-grid",
      make_check_metric(
        "Rows",
        paste0(report$raw_n_rows, " -> ", report$cleaned_n_rows),
        paste0(report$blank_rows_dropped, " fully blank row(s) dropped")
      ),
      make_check_metric(
        "Columns",
        paste0(report$raw_n_cols, " -> ", report$cleaned_n_cols),
        extra_note
      ),
      make_check_metric(
        "Model rows",
        report$cleaned_n_rows - report$rows_excluded_from_model,
        paste0(report$rows_excluded_from_model, " row(s) excluded because yield is missing")
      ),
      make_check_metric(
        "Environments",
        report$n_env,
        paste0(report$n_rep, " rep(s), ", report$n_entry, " entry level(s)")
      )
    )
  )
}

build_data_check_actions <- function(check) {
  if (is.null(check)) {
    return(tags$p("No checked data are available yet."))
  }

  report <- check$report
  action_rows <- tibble::tibble(
    check = c(
      "Required columns",
      "Blank rows",
      "Extra columns",
      "Missing yield",
      "Non-numeric row values",
      "Non-numeric column values",
      "Non-numeric yield values"
    ),
    result = c(
      "Passed",
      ifelse(report$blank_rows_dropped > 0L, "Cleaned", "Passed"),
      ifelse(report$extra_cols_ignored > 0L, "Ignored", "Passed"),
      ifelse(report$yield_missing > 0L, "Review", "Passed"),
      ifelse(report$row_conversion_na > 0L, "Review", "Passed"),
      ifelse(report$col_conversion_na > 0L, "Review", "Passed"),
      ifelse(report$yield_conversion_na > 0L, "Review", "Passed")
    ),
    detail = c(
      paste(required_cols, collapse = ", "),
      paste0(report$blank_rows_dropped, " fully blank row(s) dropped."),
      if (report$extra_cols_ignored > 0L) {
        paste0(
          report$extra_cols_ignored,
          " extra column(s) ignored",
          if (length(report$extra_nonblank_col_names)) {
            paste0("; nonblank ignored: ", paste(report$extra_nonblank_col_names, collapse = ", "))
          } else {
            "."
          }
        )
      } else {
        "No extra columns found."
      },
      paste0(report$yield_missing, " cleaned row(s) have missing yield and will not be modeled."),
      paste0(report$row_conversion_na, " value(s) became NA when row was converted to numeric."),
      paste0(report$col_conversion_na, " value(s) became NA when col was converted to numeric."),
      paste0(report$yield_conversion_na, " value(s) became NA when yield was converted to numeric.")
    )
  )

  action_rows
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
required_cols_tbl <- tibble::tibble(
  column = required_cols,
  description = c(
    "Site name or code.",
    "Trial year.",
    "Environment identifier, typically site-year.",
    "Replication identifier.",
    "Field row coordinate.",
    "Field column coordinate.",
    "Entry, genotype, or variety identifier.",
    "Observed response value to analyze."
  )
)

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
      tags$div(
        class = "sidebar-stat",
        tags$div(class = "sidebar-stat-title", "Bundled example"),
        tags$div(class = "sidebar-stat-value", "trial_sim.csv"),
        tags$div(
          class = "sidebar-stat-subtitle",
          "Simulated yield trial with 6 site-years, 20 entries, 4 reps per environment, spatial field variation, and about 20% of entries missing within each environment."
        )
      )
    ),
    conditionalPanel(
      condition = "input.data_source === 'upload'",
      fileInput("trial_file", "Upload trial CSV", accept = ".csv")
    ),
    numericInput("alpha", "Alpha", value = 0.05, min = 0.001, max = 0.50, step = 0.01),
    actionButton("run_analysis", "Run analysis", class = "btn-primary"),
    downloadButton("download_clean_data", "Download cleaned CSV", class = "btn-outline-secondary"),
    helpText("Upload a file, review the Data check tab, then run the analysis. The analysis uses the cleaned data shown in the app."),
    uiOutput("error_log_controls"),
    tags$div(style = "height: 0.5rem;"),
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
      tags$div(style = "height: 0.35rem;"),
      tags$a(
        href = "https://sites.google.com/view/ansleyjbrown",
        target = "_blank",
        rel = "noopener noreferrer",
        "Website"
      ),
      tags$br(),
      tags$a(
        href = "https://github.com/ansleybrown1337/sally-jd-yield-analysis",
        target = "_blank",
        rel = "noopener noreferrer",
        "GitHub repository"
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
      .error-log-panel {
        margin-top: 0.85rem;
        padding: 0.85rem;
        border-radius: 10px;
        background: rgba(180, 35, 24, 0.07);
        border: 1px solid rgba(180, 35, 24, 0.24);
        border-left: 5px solid #b42318;
      }
      .error-log-panel .help-block {
        margin: 0.35rem 0 0;
        color: #5f1f18;
      }
      .workflow-steps {
        display: flex;
        flex-wrap: wrap;
        gap: 0.5rem;
        margin-bottom: 1rem;
      }
      .workflow-step {
        padding: 0.45rem 0.7rem;
        border-radius: 999px;
        background: rgba(255, 255, 255, 0.68);
        border: 1px solid rgba(114, 98, 75, 0.18);
        color: #52606d;
        font-weight: 700;
        font-size: 0.9rem;
      }
      .workflow-step-active {
        background: #1f5d50;
        border-color: #1f5d50;
        color: #fff;
      }
      .workflow-step-complete {
        background: rgba(31, 93, 80, 0.12);
        border-color: rgba(31, 93, 80, 0.32);
        color: #1f5d50;
      }
      .next-action {
        padding: 0.85rem 1rem;
        margin-bottom: 1rem;
        border-left: 5px solid #1f5d50;
        border-radius: 10px;
        background: rgba(31, 93, 80, 0.10);
        color: #1f2933;
        font-size: 1.03rem;
      }
      .next-action code {
        font-weight: 700;
        color: #1f5d50;
      }
      .check-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
        gap: 0.75rem;
      }
      .check-metric {
        padding: 0.85rem 0.95rem;
        border-radius: 14px;
        background: rgba(255, 255, 255, 0.74);
        border: 1px solid rgba(114, 98, 75, 0.14);
      }
      .check-metric-title {
        font-size: 0.85rem;
        font-weight: 700;
        color: #52606d;
      }
      .check-metric-value {
        margin-top: 0.2rem;
        font-size: 1.55rem;
        font-weight: 700;
        color: #1f5d50;
        line-height: 1.1;
      }
      .check-metric-subtitle {
        margin-top: 0.25rem;
        color: #5b646d;
        font-size: 0.9rem;
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
    h2("Explore spatial variety trial results"),
    p("Use this tool to analyze crop variety trials laid out in reps, rows, and columns, whether you have one site-year or many. It helps you compare entries, account for field-position effects, and review adjusted means, BLUPs, and diagnostics in one place.")
  ),
  bslib::navset_card_tab(
    id = "main_tabs",
    selected = "About this tool",
    height = "auto",
    bslib::nav_panel(
      "About this tool",
      bslib::card(
        bslib::card_header("How to use"),
        tags$div(
          class = "workflow-steps",
          tags$div(class = "workflow-step workflow-step-active", "1. Prepare data"),
          tags$div(class = "workflow-step", "2. Run analysis"),
          tags$div(class = "workflow-step", "3. Review results"),
          tags$div(class = "workflow-step", "4. Check diagnostics")
        ),
        tags$ol(
          tags$li("Upload a CSV or choose the bundled example."),
          tags$li("Review the ", tags$strong("Data check"), " tab to see blank rows, extra columns, missing yields, and the cleaned data that will be analyzed."),
          tags$li("Set the alpha level in the sidebar, then click ", tags$code("Run analysis"), "."),
          tags$li("Review adjusted means, BLUPs, diagnostics, model fit, and the run summary."),
          tags$li("Download the cleaned CSV, run report, or output ZIP from the sidebar if needed.")
        )
      ),
      bslib::card(
        bslib::card_header("Purpose"),
        tags$p("This app is for crop variety trial data collected from field layouts with replications, rows, and columns. You can use it for a single site-year trial or for multiple site-years combined in one file, as long as the data follow the required CSV structure."),
        tags$p("After you run the analysis, the app helps you compare variety performance with adjusted means, review one-stage BLUP results, and check diagnostics that account for field-position patterns and possible outliers."),
        tags$p("Stage 1 fits each site-year separately so the app can account for field-position effects within that environment and produce adjusted entry means plus diagnostics. Stage 2 then uses those Stage 1 results in one of two ways: if you uploaded multiple site-years, it combines the adjusted means across environments so entries can be compared across site-years; if you uploaded only one site-year, it skips the across-environment step and reports that trial on its own."),
        tags$p("BLUPs are also reported as a separate one-stage mixed-model analysis. They are useful for ranking entries, but they answer a different question than the adjusted means shown from the Stage 1 to Stage 2 workflow.")
      ),
      bslib::card(
        bslib::card_header("Important assumptions"),
        tags$ul(
          tags$li("Input data must include columns: ", tags$code("site, year, env, rep, row, col, entry, yield"), "."),
          tags$li("Rows with missing ", tags$code("yield"), " values are excluded from model fitting."),
          tags$li("Stage 1 is the per-environment fitting step. It uses spatial covariance models where possible and falls back to RCBD or fixed-effects ANOVA if needed."),
          tags$li("Stage 2 is the adjusted-means summary step. With multiple environments, it combines Stage 1 adjusted means across site-years; with one environment, it reports that environment directly."),
          tags$li("Across-environment LS-means and BLUPs answer different questions and should be interpreted separately."),
          tags$li("Flagged outliers are diagnostic cues based on large normalized residuals, not automatic deletions.")
        )
      )
    ),
    bslib::nav_panel(
      "Data format requirements",
      bslib::card(
        bslib::card_header("Required columns"),
        tags$p("Input files must be CSVs containing the following columns. The app uses these fields to define environments, field position, entries, and the response variable."),
        DT::DTOutput("required_columns_table")
      ),
      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          bslib::card_header("Single-site example"),
          tags$p("Example rows from one environment in the bundled simulated dataset."),
          DT::DTOutput("single_site_example_table")
        ),
        bslib::card(
          bslib::card_header("Multi-site example"),
          tags$p("Example rows from the bundled simulated dataset with multiple environments included."),
          DT::DTOutput("multi_site_example_table")
        )
      )
    ),
    bslib::nav_panel(
      "Data check",
      bslib::card(
        bslib::card_header("Upload check"),
        uiOutput("data_check_summary")
      ),
      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          bslib::card_header("Cleaning actions"),
          DT::DTOutput("data_check_actions_table")
        ),
        bslib::card(
          bslib::card_header("Environment summary before modeling"),
          DT::DTOutput("data_check_env_table")
        )
      )
    ),
    bslib::nav_panel(
      "Input data",
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Input data"),
        tags$p("This table shows the cleaned data used for the analysis after blank rows, extra columns, and missing-value cleanup are applied."),
        DT::DTOutput("trial_head_table")
      )
    ),
    bslib::nav_panel(
      "Adjusted means",
      tags$div(
        class = "status-summary",
        tags$p("This tab shows adjusted LS-means from the Stage 1 to Stage 2 workflow. The error bars are confidence intervals for each adjusted mean at the user-selected alpha level, and the letters are Tukey-adjusted significance groupings: entries that share a letter are not significantly different at that alpha level. This is separate from the one-stage BLUP analysis shown later.")
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
          full_screen = TRUE,
          bslib::card_header("Backend run report"),
          verbatimTextOutput("run_report", placeholder = TRUE)
        )
      )
    ),
    bslib::nav_panel(
      "Interpreting results",
      bslib::card(
        bslib::card_header("Interpretation guide"),
        uiOutput("interpretation_guide_html")
      )
    )
  )
)

server <- function(input, output, session) {
  analysis_data <- reactiveVal(NULL)
  input_check <- reactiveVal(NULL)
  error_logs <- reactiveVal(character())

  append_error_log <- function(stage, error, context = list()) {
    log_text <- build_error_log_text(stage, error, context)
    error_logs(c(error_logs(), log_text))
  }

  current_error_context <- function(extra = list()) {
    check <- input_check()
    uploaded <- input$trial_file
    base_context <- list(
      data_source = if (is.null(input$data_source)) "<not selected>" else input$data_source,
      uploaded_file = if (is.null(uploaded) || is.null(uploaded$name)) "<none>" else uploaded$name,
      alpha = if (is.null(input$alpha)) "<not set>" else input$alpha,
      checked_source = if (is.null(check) || is.null(check$source_name)) "<none>" else check$source_name,
      checked_rows = if (is.null(check)) "<none>" else nrow(check$trial),
      checked_envs = if (is.null(check)) "<none>" else check$report$n_env
    )
    c(base_context, extra)
  }

  chosen_input_path <- reactive({
    if (identical(input$data_source, "upload")) {
      req(input$trial_file$datapath)
      return(input$trial_file$datapath)
    }
    req(file.exists(bundled_example_path))
    bundled_example_path
  })

  prepare_input_check <- function(path, source_name) {
    raw <- read.csv(path, stringsAsFactors = FALSE)
    checked <- backend_env$clean_trial_data_with_report(raw)
    clean_path <- tempfile(pattern = "cleaned-trial-", fileext = ".csv")
    write.csv(checked$data, clean_path, row.names = FALSE, na = "")

    list(
      raw = raw,
      trial = checked$data,
      report = checked$report,
      source_name = source_name,
      input_path = path,
      clean_path = clean_path
    )
  }

  observeEvent(input$data_source, {
    analysis_data(NULL)
    if (identical(input$data_source, "bundled")) {
      tryCatch({
        input_check(prepare_input_check(bundled_example_path, basename(bundled_example_path)))
        updateTabsetPanel(session, "main_tabs", selected = "Data check")
      }, error = function(e) {
        input_check(NULL)
        append_error_log(
          "Preparing bundled example data check",
          e,
          current_error_context(list(source_path = bundled_example_path))
        )
        showNotification(conditionMessage(e), type = "error", duration = NULL)
      })
    } else if (is.null(input$trial_file)) {
      input_check(NULL)
    }
  }, ignoreInit = TRUE)

  observeEvent(input$trial_file, {
    analysis_data(NULL)
    req(input$trial_file$datapath)

    tryCatch({
      input_check(prepare_input_check(input$trial_file$datapath, input$trial_file$name))
      updateTabsetPanel(session, "main_tabs", selected = "Data check")
      showNotification("Data check complete. Review the Data check tab, then run the analysis.", type = "message", duration = 5)
    }, error = function(e) {
      input_check(NULL)
      append_error_log(
        "Preparing uploaded data check",
        e,
        current_error_context(list(source_path = input$trial_file$datapath))
      )
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  }, ignoreInit = TRUE)

  observeEvent(input$run_analysis, {
    tryCatch({
      check <- input_check()
      if (is.null(check)) {
        check <- prepare_input_check(chosen_input_path(), basename(chosen_input_path()))
        input_check(check)
      }

      in_csv <- check$clean_path
      trial <- check$trial

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
          input_path = check$input_path,
          clean_input_path = in_csv,
          input_check = check,
          out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE)
        ))
      })

      updateTabsetPanel(session, "main_tabs", selected = "Adjusted means")
      showNotification("Analysis complete.", type = "message", duration = 4)
    }, error = function(e) {
      append_error_log("Running analysis", e, current_error_context())
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })

  output$error_log_controls <- renderUI({
    logs <- error_logs()
    if (!length(logs)) return(NULL)

    tags$div(
      class = "error-log-panel",
      downloadButton("download_error_log", "Download error log", class = "btn-outline-danger"),
      helpText("An error was recorded in this session. Download this log and send it with the CSV if the problem persists.")
    )
  })

  output$download_error_log <- downloadHandler(
    filename = function() {
      paste0("sally-yield-analysis-error-log-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".txt")
    },
    content = function(file) {
      logs <- error_logs()
      req(length(logs) > 0)
      writeLines(paste(logs, collapse = "\n\n---\n\n"), con = file, useBytes = TRUE)
    }
  )

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

  output$download_clean_data <- downloadHandler(
    filename = function() {
      check <- input_check()
      base_name <- if (is.null(check) || is.null(check$source_name)) {
        "trial"
      } else {
        tools::file_path_sans_ext(basename(check$source_name))
      }
      paste0(base_name, "-cleaned.csv")
    },
    content = function(file) {
      check <- input_check()
      req(check)
      write.csv(check$trial, file, row.names = FALSE, na = "")
    }
  )

  output$data_check_summary <- renderUI({
    build_data_check_summary(input_check())
  })

  output$data_check_actions_table <- DT::renderDT({
    check <- input_check()
    if (is.null(check)) {
      return(DT::datatable(
        tibble::tibble(check = character(), result = character(), detail = character()),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }

    DT::datatable(
      build_data_check_actions(check),
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
      rownames = FALSE
    )
  })

  output$data_check_env_table <- DT::renderDT({
    check <- input_check()
    if (is.null(check)) {
      return(DT::datatable(
        tibble::tibble(env = character(), n_rows = integer(), n_model_rows = integer()),
        options = list(dom = "t"),
        rownames = FALSE
      ))
    }

    DT::datatable(check$report$env_summary, options = datatable_opts, rownames = FALSE)
  })

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

  output$required_columns_table <- DT::renderDT({
    DT::datatable(required_cols_tbl, options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })

  output$interpretation_guide_html <- renderUI({
    render_markdown_file_html(interpretation_md_path)
  })

  output$single_site_example_table <- DT::renderDT({
    df <- read.csv(bundled_example_path, stringsAsFactors = FALSE)
    env_id <- sort(unique(df$env))[1]
    df_one_env <- df[df$env == env_id, , drop = FALSE]
    DT::datatable(utils::head(df_one_env, 12), options = c(datatable_opts, list(scrollY = "320px")), rownames = FALSE)
  })

  output$multi_site_example_table <- DT::renderDT({
    df <- read.csv(bundled_example_path, stringsAsFactors = FALSE)
    env_ids <- sort(unique(df$env))
    df_multi_env <- do.call(
      rbind,
      lapply(env_ids[seq_len(min(4, length(env_ids)))], function(env_id) {
        utils::head(df[df$env == env_id, , drop = FALSE], 3)
      })
    )
    DT::datatable(df_multi_env, options = c(datatable_opts, list(scrollY = "320px")), rownames = FALSE)
  })

  output$trial_head_table <- DT::renderDT({
    state <- analysis_data()
    req(state)
    DT::datatable(
      state$trial,
      options = c(datatable_opts, list(pageLength = 25, scrollX = TRUE)),
      rownames = FALSE,
      filter = "top"
    )
  })

  output$run_report <- renderText({
    state <- analysis_data()
    req(state)
    report_path <- state$result$report_path
    paste(readLines(report_path), collapse = "\n")
  })
}

shinyApp(ui, server)
