suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(plotly)
  library(DT)
  library(dplyr)
  library(readr)
  library(tibble)
})

app_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
backend_env <- new.env(parent = globalenv())
source(file.path(app_root, "code", "r_equivalent.R"), local = backend_env)

required_cols <- c("site", "year", "env", "rep", "row", "col", "entry", "yield")
spatial_cov_types <- c("expa", "exp", "sph", "gau")

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

make_metric_card <- function(title, value, subtitle) {
  bslib::card(
    class = "metric-card",
    bslib::card_header(title),
    div(class = "metric-value", value),
    div(class = "metric-subtitle", subtitle)
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
    p <- p %>%
      plotly::add_text(
        data = dplyr::filter(plot_df, nzchar(group)),
        x = ~upper.CL,
        y = ~entry,
        text = ~group,
        textposition = "middle right",
        textfont = list(color = "#1f5d50", size = 12),
        hoverinfo = "skip",
        showlegend = FALSE
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
      selected = "bundled"
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
    textInput("out_dir", "Output folder", value = "shiny_output"),
    actionButton("run_analysis", "Run analysis", class = "btn-primary"),
    tags$hr(),
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
    helpText("Diagnostics stay focused on one site-year at a time so multi-environment runs remain readable.")
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
  uiOutput("metric_cards"),
  bslib::navset_card_tab(
    height = "auto",
    bslib::nav_panel(
      "Across environments",
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Adjusted means"),
        plotly::plotlyOutput("across_plot", height = "560px")
      ),
      bslib::layout_column_wrap(
        width = 1/2,
        bslib::card(
          bslib::card_header("Across-environment LS-means"),
          DT::DTOutput("lsm_table")
        ),
        bslib::card(
          bslib::card_header("Compact letter display"),
          DT::DTOutput("cld_table")
        )
      )
    ),
    bslib::nav_panel(
      "Diagnostics",
      bslib::card(
        full_screen = TRUE,
        bslib::card_header("Field map"),
        plotly::plotlyOutput("heatmap_plot", height = "560px")
      ),
      bslib::layout_column_wrap(
        width = 1/3,
        bslib::card(
          bslib::card_header("Residuals vs fitted"),
          plotly::plotlyOutput("resid_fit_plot", height = "360px")
        ),
        bslib::card(
          bslib::card_header("Q-Q plot"),
          plotly::plotlyOutput("qq_plot", height = "360px")
        ),
        bslib::card(
          bslib::card_header("Residual distribution"),
          plotly::plotlyOutput("hist_plot", height = "360px")
        )
      ),
      bslib::card(
        bslib::card_header("Flagged outliers for the selected environment"),
        DT::DTOutput("outlier_table")
      )
    ),
    bslib::nav_panel(
      "Model selection",
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

      out_dir <- input$out_dir
      if (!grepl("^[A-Za-z]:|^/", out_dir)) {
        out_dir <- file.path(app_root, out_dir)
      }
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

  output$metric_cards <- renderUI({
    state <- analysis_data()
    if (is.null(state)) {
      return(
        bslib::layout_column_wrap(
          width = 1/3,
          make_metric_card("Environments", "0", "Run the workflow to populate results"),
          make_metric_card("Spatial Fits", "0", "No Stage 1 model selection yet"),
          make_metric_card("Flagged Outliers", "0", "Diagnostics will appear after the first run")
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

    bslib::layout_column_wrap(
      width = 1/3,
      make_metric_card("Environments", n_env, "Site-years with usable Stage 1 means"),
      make_metric_card("Spatial Fits", n_spatial, paste0(n_fallback, " environments used fallback models")),
      make_metric_card("Flagged Outliers", outlier_count, stage2_mode)
    )
  })

  output$across_plot <- plotly::renderPlotly({
    state <- analysis_data()
    req(state)
    make_across_plot(build_stage2_plot_data(state$result))
  })

  output$heatmap_plot <- plotly::renderPlotly({
    make_heatmap_plot(current_env_data(), input$heatmap_fill)
  })

  output$resid_fit_plot <- plotly::renderPlotly({
    make_resid_fit_plot(current_env_data())
  })

  output$qq_plot <- plotly::renderPlotly({
    make_qq_plot(current_env_data())
  })

  output$hist_plot <- plotly::renderPlotly({
    make_hist_plot(current_env_data())
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

  output$outlier_table <- DT::renderDT({
    df <- current_env_data() %>%
      dplyr::filter(flag_outlier %in% TRUE) %>%
      dplyr::select(env, site, year, rep, row, col, entry, yield, fitted, resid, resid_norm)

    DT::datatable(df, options = datatable_opts, rownames = FALSE)
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
