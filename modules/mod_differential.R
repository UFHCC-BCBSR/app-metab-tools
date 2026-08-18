# ===========================
# Differential Analysis — Shiny module
# ===========================
# UI and server for the Differential Analysis tab. All statistics, plots and
# report content live in modules/de_analysis.R; this file only wires them to
# inputs and outputs.
#
# The module is handed one reactive, `de_inputs`, which returns NULL until
# preprocessing has produced a log2-normalized matrix, and otherwise a list of:
#   normalized       log2-normalized matrix (features x samples)
#   batch_corrected  ComBat-corrected matrix, or NULL before ComBat has run
#   sample_data      matched metadata (sample_name, batch, biological_var)
#   full_sample_data all metadata columns, rows aligned to matrix columns
#   bio_var_name     name of the biological variable column, or NULL
# It returns a reactive holding the saved comparisons, which the download
# handlers fold into the Excel package and the HTML report.

differentialUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = "Set Up Comparison", status = "warning", solidHeader = TRUE, width = 12,
          conditionalPanel(
            condition = "!output.de_ready", ns = ns,
            p("Run preprocessing on the Preprocess Data tab first. Differential analysis needs the log2-normalized matrix.")
          ),
          conditionalPanel(
            condition = "output.de_ready", ns = ns,
            fluidRow(
              column(4,
                     div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; min-height: 220px;",
                         p(strong("What to compare:")),
                         selectInput(ns("group_var"), "Grouping variable:", choices = NULL),
                         radioButtons(ns("mode"), "Comparison type:",
                                      choices = list(
                                        "Two groups (moderated t-test)" = "pairwise",
                                        "All groups (moderated F-test)" = "global"
                                      ), selected = "pairwise")
                     )
              ),
              column(4,
                     div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; min-height: 220px;",
                         conditionalPanel(
                           condition = "input.mode == 'pairwise'", ns = ns,
                           p(strong("Groups:")),
                           selectInput(ns("group_a"), "Group A:", choices = NULL),
                           selectInput(ns("group_b"), "Group B (reference):", choices = NULL),
                           p(style = "color: #666; font-size: 12px;",
                             "Log2 fold change is A minus B: positive values are higher in group A.")
                         ),
                         conditionalPanel(
                           condition = "input.mode == 'global'", ns = ns,
                           p(strong("Groups:")),
                           p(style = "color: #666; font-size: 12px;",
                             "Every level of the grouping variable is tested at once. The F-test asks whether a feature differs across any of them, without picking a direction.")
                         )
                     )
              ),
              column(4,
                     div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; min-height: 220px;",
                         p(strong("Data and thresholds:")),
                         uiOutput(ns("matrix_choice")),
                         numericInput(ns("fdr_cutoff"), "FDR cutoff:",
                                      value = 0.05, min = 0.001, max = 1, step = 0.01),
                         conditionalPanel(
                           condition = "input.mode == 'pairwise'", ns = ns,
                           numericInput(ns("lfc_cutoff"), "Minimum |log2 fold change|:",
                                        value = 1, min = 0, max = 10, step = 0.5)
                         )
                     )
              )
            ),
            br(),
            verbatimTextOutput(ns("model_summary")),
            br(),
            actionButton(ns("run_de"), "Run Differential Analysis",
                         class = "btn-warning btn-lg", icon = icon("flask"))
          )
      )
    ),

    conditionalPanel(
      condition = "output.has_results", ns = ns,
      fluidRow(
        box(title = "Results", status = "success", solidHeader = TRUE, width = 12,
            verbatimTextOutput(ns("result_summary")),
            br(),
            tabsetPanel(
              tabPanel("Results table",
                       br(),
                       DT::dataTableOutput(ns("results_table")),
                       br(),
                       downloadButton(ns("download_csv"), "Download this table (CSV)",
                                      class = "btn-primary")
              ),
              tabPanel("Volcano plot",
                       br(),
                       uiOutput(ns("main_plot_note")),
                       plotlyOutput(ns("main_plot"), height = "500px")
              ),
              tabPanel("Heatmap",
                       br(),
                       uiOutput(ns("heatmap_note")),
                       plotlyOutput(ns("heatmap"), height = "650px")
              ),
              tabPanel("Feature detail",
                       br(),
                       fluidRow(
                         column(6, selectizeInput(ns("feature"), "Feature:", choices = NULL,
                                                  width = "100%"))
                       ),
                       p(style = "color: #666; font-size: 12px;",
                         "Selecting a row in the results table also loads that feature here."),
                       plotlyOutput(ns("boxplot"), height = "450px"),
                       br(),
                       verbatimTextOutput(ns("feature_stats"))
              )
            )
        )
      ),
      fluidRow(
        box(title = "Saved Comparisons", status = "info", solidHeader = TRUE, width = 12,
            p("Every comparison you run is saved and included in the Excel package and HTML report. Re-running the same comparison replaces its saved copy."),
            DT::dataTableOutput(ns("saved_table")),
            br(),
            actionButton(ns("clear_runs"), "Clear saved comparisons", class = "btn-danger")
        )
      )
    )
  )
}

differentialServer <- function(id, de_inputs) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(current = NULL, runs = list())

    output$de_ready <- reactive({ !is.null(de_inputs()) })
    output$has_results <- reactive({ !is.null(rv$current) })
    outputOptions(output, "de_ready", suspendWhenHidden = FALSE)
    outputOptions(output, "has_results", suspendWhenHidden = FALSE)

    # A saved comparison belongs to the matrix it was fitted on. Re-running
    # preprocessing or ComBat would otherwise leave results from the old data
    # sitting in the same report as results from the new data. Only comparisons
    # whose own matrix changed are dropped: running ComBat after a comparison on
    # the normalized matrix leaves that comparison untouched.
    current_fingerprints <- reactive({
      dat <- de_inputs()
      list(
        normalized = de_mat_fingerprint(dat$normalized),
        batch_corrected = de_mat_fingerprint(dat$batch_corrected)
      )
    })

    observe({
      dat <- de_inputs()
      if (is.null(dat)) {
        if (length(rv$runs) > 0 || !is.null(rv$current)) {
          rv$runs <- list()
          rv$current <- NULL
        }
        return()
      }
      if (length(rv$runs) == 0 && is.null(rv$current)) return()
      fps <- current_fingerprints()
      is_stale <- function(run) {
        !identical(run$source_fingerprint, fps[[run$matrix_source]])
      }
      stale <- vapply(rv$runs, is_stale, logical(1))
      if (any(stale)) {
        rv$runs <- rv$runs[!stale]
        showNotification(
          paste0("Data changed — ", sum(stale),
                 " saved comparison(s) no longer match the current data and were removed."),
          type = "warning", duration = 10)
      }
      if (!is.null(rv$current) && is_stale(rv$current)) rv$current <- NULL
    })

    # ---- input choices -------------------------------------------------
    observe({
      dat <- de_inputs()
      req(dat)
      choices <- de_group_candidates(dat$full_sample_data)
      if (length(choices) == 0) return()
      selected <- isolate(input$group_var)
      if (is.null(selected) || !selected %in% choices) {
        selected <- if (!is.null(dat$bio_var_name) && dat$bio_var_name %in% choices) {
          dat$bio_var_name
        } else {
          choices[1]
        }
      }
      updateSelectInput(session, "group_var", choices = choices, selected = selected)
    })

    group_vector <- reactive({
      dat <- de_inputs()
      req(dat, input$group_var)
      req(input$group_var %in% colnames(dat$full_sample_data))
      droplevels(as.factor(dat$full_sample_data[[input$group_var]]))
    })

    observe({
      g <- group_vector()
      lv <- levels(g)
      req(length(lv) >= 2)
      sel_a <- isolate(input$group_a); sel_b <- isolate(input$group_b)
      if (is.null(sel_a) || !sel_a %in% lv) sel_a <- lv[2]
      if (is.null(sel_b) || !sel_b %in% lv || identical(sel_b, sel_a)) sel_b <- lv[1]
      updateSelectInput(session, "group_a", choices = lv, selected = sel_a)
      updateSelectInput(session, "group_b", choices = lv, selected = sel_b)
    })

    output$matrix_choice <- renderUI({
      dat <- de_inputs()
      req(dat)
      combat_ready <- !is.null(dat$batch_corrected)
      choices <- list("Log2-normalized, batch as covariate" = "normalized")
      if (combat_ready) {
        choices[["ComBat-corrected, no batch term"]] <- "batch_corrected"
      }
      selected <- isolate(input$matrix_source)
      if (is.null(selected) || !selected %in% unlist(choices)) selected <- "normalized"
      tagList(
        radioButtons(ns("matrix_source"), "Data to model:",
                     choices = choices, selected = selected),
        if (!combat_ready) {
          p(style = "color: #666; font-size: 12px;",
            "Run batch correction to also model the ComBat-corrected matrix.")
        }
      )
    })

    # ---- model description shown before running ------------------------
    output$model_summary <- renderText({
      dat <- de_inputs()
      req(dat, input$group_var, input$mode)
      g <- group_vector()
      batch <- droplevels(as.factor(dat$sample_data$batch))
      use_batch <- identical(input$matrix_source, "normalized") && nlevels(batch) > 1

      counts <- table(g)
      counts_text <- paste(paste0(names(counts), ": ", as.integer(counts)), collapse = ", ")
      design_text <- paste0("~ 0 + ", input$group_var, if (use_batch) " + batch" else "")
      data_text <- if (identical(input$matrix_source, "batch_corrected")) {
        "ComBat batch-corrected log2 data (batch already removed from the matrix)"
      } else if (use_batch) {
        "log2-normalized data (batch effects controlled for in the model)"
      } else {
        "log2-normalized data (only one batch present, so no batch term)"
      }
      test_text <- if (identical(input$mode, "pairwise")) {
        paste0("moderated t-test: ", input$group_a, " - ", input$group_b)
      } else {
        paste0("moderated F-test across all ", nlevels(g), " levels of ", input$group_var)
      }

      paste0(
        "Model to be fitted: ", design_text, "\n",
        "Test: ", test_text, "\n",
        "Data: ", data_text, "\n",
        "Samples per group: ", counts_text
      )
    })

    # ---- run -----------------------------------------------------------
    observeEvent(input$run_de, {
      dat <- de_inputs()
      req(dat, input$group_var, input$mode)
      shinyjs::disable("run_de")
      showModal(modalDialog("Fitting models...", footer = NULL))
      tryCatch({
        mat <- if (identical(input$matrix_source, "batch_corrected")) {
          dat$batch_corrected
        } else {
          dat$normalized
        }
        if (is.null(mat)) stop("The selected data matrix is not available yet.")

        batch <- droplevels(as.factor(dat$sample_data$batch))
        use_batch <- identical(input$matrix_source, "normalized") && nlevels(batch) > 1

        run <- de_run(
          mat = mat,
          groups = group_vector(),
          batch = batch,
          mode = input$mode,
          group_a = input$group_a,
          group_b = input$group_b,
          matrix_source = if (identical(input$matrix_source, "batch_corrected")) {
            "batch_corrected"
          } else {
            "normalized"
          },
          use_batch = use_batch,
          group_var = input$group_var,
          fdr_cutoff = input$fdr_cutoff,
          lfc_cutoff = if (identical(input$mode, "pairwise")) input$lfc_cutoff else 0
        )

        run$source_fingerprint <- de_mat_fingerprint(mat)
        rv$current <- run
        key <- paste0(run$label, " [", run$matrix_source, "]")
        rv$runs[[key]] <- run

        top_features <- run$results$feature
        updateSelectizeInput(session, "feature", choices = top_features,
                             selected = top_features[1], server = TRUE)

        removeModal()
        shinyjs::enable("run_de")
        showNotification(paste0("Differential analysis complete: ", run$n_sig,
                                " significant features."), type = "message")
      }, error = function(e) {
        removeModal()
        shinyjs::enable("run_de")
        showNotification(paste("Differential analysis failed:", e$message),
                         type = "error", duration = 12)
      })
    })

    observeEvent(input$clear_runs, {
      rv$runs <- list()
      rv$current <- NULL
      showNotification("Saved comparisons cleared.", type = "message")
    })

    # ---- results -------------------------------------------------------
    output$result_summary <- renderText({
      req(rv$current)
      de_design_summary(rv$current)
    })

    results_display <- reactive({
      req(rv$current)
      res <- rv$current$results
      res$significant <- ifelse(res$significant, "Yes", "No")
      res
    })

    output$results_table <- DT::renderDataTable({
      res <- results_display()
      num_cols <- setdiff(colnames(res)[vapply(res, is.numeric, logical(1))],
                          c("p_value", "fdr"))
      DT::datatable(
        res,
        selection = "single",
        rownames = FALSE,
        options = list(scrollX = TRUE, pageLength = 15,
                       order = list())
      ) %>%
        DT::formatSignif(columns = num_cols, digits = 4) %>%
        DT::formatSignif(columns = c("p_value", "fdr"), digits = 3)
    })

    observeEvent(input$results_table_rows_selected, {
      req(rv$current)
      i <- input$results_table_rows_selected
      req(length(i) == 1)
      updateSelectizeInput(session, "feature",
                           choices = rv$current$results$feature,
                           selected = rv$current$results$feature[i],
                           server = TRUE)
    })

    output$main_plot_note <- renderUI({
      req(rv$current)
      if (rv$current$mode == "pairwise") {
        p(style = "color: #666;",
          paste0("Coloured points clear both the FDR and fold change thresholds. ",
                 "Positive log2 fold changes are higher in ", rv$current$group_a, "."))
      } else {
        p(style = "color: #666;",
          "An F-test has no single fold change, so features are plotted against mean intensity. Coloured points clear the FDR threshold.")
      }
    })

    output$main_plot <- renderPlotly({
      req(rv$current)
      if (rv$current$mode == "pairwise") {
        de_volcano_plot(rv$current)
      } else {
        de_significance_plot(rv$current)
      }
    })

    output$heatmap_note <- renderUI({
      req(rv$current)
      sel <- de_display_features(rv$current)
      if (sel$used_fallback) {
        p(style = "color: #666;",
          "No feature cleared the thresholds, so the strongest features by p-value are shown. Values are per-feature z-scores.")
      } else {
        p(style = "color: #666;",
          paste0("Showing ", min(length(sel$features), DE_MAX_HEATMAP_FEATURES), " of ",
                 sel$n_available,
                 " significant features as per-feature z-scores, with samples grouped by ",
                 rv$current$group_var, "."))
      }
    })

    output$heatmap <- renderPlotly({
      req(rv$current)
      p <- de_heatmap_plot(rv$current)
      validate(need(!is.null(p), "At least two features are needed to draw a heatmap."))
      p
    })

    output$boxplot <- renderPlotly({
      req(rv$current, input$feature)
      p <- de_feature_boxplot(rv$current, input$feature)
      validate(need(!is.null(p), "Select a feature to plot."))
      p
    })

    output$feature_stats <- renderText({
      req(rv$current, input$feature)
      res <- rv$current$results
      row <- res[res$feature == input$feature, , drop = FALSE]
      req(nrow(row) == 1)
      lines <- c(paste0("Feature: ", row$feature))
      if (rv$current$mode == "pairwise") {
        lines <- c(lines,
                   paste0("log2 fold change: ", signif(row$log2FC, 4),
                          "  (fold change ", signif(row$fold_change, 4), ")"),
                   paste0("t: ", signif(row$t, 4)))
      } else {
        lines <- c(lines,
                   paste0("F: ", signif(row$F, 4)),
                   paste0("Largest |log2FC|: ", signif(row$max_abs_log2FC, 4)))
      }
      c(lines,
        paste0("p-value: ", signif(row$p_value, 3), "   FDR: ", signif(row$fdr, 3)),
        paste0("Mean log2 intensity: ", signif(row$mean_intensity, 4)),
        paste0("Significant at current thresholds: ",
               ifelse(row$significant, "yes", "no"))) %>%
        paste(collapse = "\n")
    })

    output$saved_table <- DT::renderDataTable({
      req(length(rv$runs) > 0)
      df <- do.call(rbind, lapply(seq_along(rv$runs), function(i) {
        run <- rv$runs[[i]]
        data.frame(
          Sheet = de_sheet_name(run, i),
          Comparison = run$label,
          Test = if (run$mode == "pairwise") "moderated t-test" else "moderated F-test",
          Data = de_matrix_label(run),
          Thresholds = paste0("FDR <= ", run$fdr_cutoff,
                              if (run$mode == "pairwise") {
                                paste0(", |log2FC| >= ", run$lfc_cutoff)
                              } else ""),
          Significant = run$n_sig,
          check.names = FALSE, stringsAsFactors = FALSE
        )
      }))
      DT::datatable(df, rownames = FALSE, selection = "none",
                    options = list(scrollX = TRUE, dom = "t", paging = FALSE))
    })

    output$download_csv <- downloadHandler(
      filename = function() {
        req(rv$current)
        paste0(gsub("[^A-Za-z0-9]+", "_", rv$current$label),
               "_differential_results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(rv$current)
        write.csv(rv$current$results, file, row.names = FALSE)
      }
    )

    reactive(rv$runs)
  })
}
