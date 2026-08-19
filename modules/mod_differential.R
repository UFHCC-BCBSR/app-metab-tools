# ===========================
# Differential Analysis — Shiny module
# ===========================
# UI and server for the Differential Analysis tab. All statistics, plots and
# report content live in modules/de_analysis.R; this file only wires them to
# inputs and outputs.
#
# The module is handed one reactive, `modes`, returning the list of processed
# ionization modes (empty until at least one mode has been preprocessed). Each
# mode carries its own matrices and metadata; the module fits them separately
# and corrects once across all of them.
#
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
            uiOutput(ns("mode_note")),
            fluidRow(
              column(4,
                     div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; min-height: 260px;",
                         p(strong("What to compare:")),
                         selectInput(ns("group_var"), "Grouping variable:", choices = NULL),
                         radioButtons(ns("comparison"), "Comparison type:",
                                      choices = list(
                                        "Two groups (moderated t-test)" = "pairwise",
                                        "Several groups (moderated F-test)" = "global"
                                      ), selected = "pairwise")
                     )
              ),
              column(4,
                     div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; min-height: 260px;",
                         conditionalPanel(
                           condition = "input.comparison == 'pairwise'", ns = ns,
                           p(strong("Groups:")),
                           selectInput(ns("group_a"), "Group A:", choices = NULL),
                           selectInput(ns("group_b"), "Group B (reference):", choices = NULL),
                           p(style = "color: #666; font-size: 12px;",
                             "Log2 fold change is A minus B: positive values are higher in group A.")
                         ),
                         conditionalPanel(
                           condition = "input.comparison == 'global'", ns = ns,
                           p(strong("Groups to include:")),
                           selectInput(ns("groups_included"), NULL, choices = NULL, multiple = TRUE),
                           p(style = "color: #666; font-size: 12px;",
                             "The F-test asks whether a feature differs across any of the selected groups. Unselected groups are left out of the model entirely, so they do not affect the variance estimate.")
                         )
                     )
              ),
              column(4,
                     div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; min-height: 260px;",
                         p(strong("Data and thresholds:")),
                         uiOutput(ns("matrix_choice")),
                         numericInput(ns("fdr_cutoff"), "FDR cutoff:",
                                      value = 0.05, min = 0.001, max = 1, step = 0.01),
                         conditionalPanel(
                           condition = "input.comparison == 'pairwise'", ns = ns,
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
                       uiOutput(ns("heatmap_mode_choice")),
                       uiOutput(ns("heatmap_note")),
                       plotlyOutput(ns("heatmap"), height = "650px")
              ),
              tabPanel("Feature detail",
                       br(),
                       fluidRow(
                         column(6, selectizeInput(ns("feature"), "Feature:", choices = NULL, width = "100%"))
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

differentialServer <- function(id, modes) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    rv <- reactiveValues(current = NULL, runs = list())

    active_modes <- reactive({
      ms <- modes()
      if (is.null(ms)) list() else ms
    })

    output$de_ready <- reactive({ length(active_modes()) > 0 })
    output$has_results <- reactive({ !is.null(rv$current) })
    outputOptions(output, "de_ready", suspendWhenHidden = FALSE)
    outputOptions(output, "has_results", suspendWhenHidden = FALSE)

    output$mode_note <- renderUI({
      ms <- active_modes()
      req(length(ms) > 0)
      if (length(ms) == 1) {
        p(style = "color: #666;",
          paste0("One ionization mode loaded (", ms[[1]]$label, ")."))
      } else {
        p(style = "color: #666;",
          paste0(length(ms), " ionization modes loaded (",
                 paste(vapply(ms, function(m) m$label, character(1)), collapse = ", "),
                 "). Each is modelled separately and the false discovery rate is corrected once across both, so the FDR applies to the experiment rather than to one acquisition."))
      }
    })

    # A saved comparison belongs to the matrices it was fitted on. Re-running
    # preprocessing or ComBat would otherwise leave results from the old data
    # sitting in the same report as results from the new data.
    current_fingerprints <- reactive({
      ms <- active_modes()
      out <- list()
      for (m in ms) {
        out[[m$label]] <- list(
          covariate = de_mat_fingerprint(m$normalized),
          none = de_mat_fingerprint(m$normalized),
          batch_corrected = de_mat_fingerprint(m$batch_corrected)
        )
      }
      out
    })

    run_fingerprint <- function(run, fps) {
      labs <- run$mode_labels
      if (!all(labs %in% names(fps))) return(NULL)
      vapply(labs, function(l) {
        f <- fps[[l]][[run$matrix_source]]
        if (is.null(f)) NA_character_ else f
      }, character(1))
    }

    observe({
      ms <- active_modes()
      if (length(ms) == 0) {
        if (length(rv$runs) > 0 || !is.null(rv$current)) {
          rv$runs <- list()
          rv$current <- NULL
        }
        return()
      }
      if (length(rv$runs) == 0 && is.null(rv$current)) return()
      fps <- current_fingerprints()
      is_stale <- function(run) !identical(run$source_fingerprint, run_fingerprint(run, fps))
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
      ms <- active_modes()
      req(length(ms) > 0)
      choices <- de_shared_group_candidates(ms)
      if (length(choices) == 0) return()
      selected <- isolate(input$group_var)
      if (is.null(selected) || !selected %in% choices) {
        bio <- ms[[1]]$params$bio_var
        selected <- if (!is.null(bio) && bio %in% choices) bio else choices[1]
      }
      updateSelectInput(session, "group_var", choices = choices, selected = selected)
    })

    group_levels <- reactive({
      ms <- active_modes()
      req(length(ms) > 0, input$group_var)
      req(input$group_var %in% colnames(ms[[1]]$matched_data$full_sample_data))
      de_group_levels(ms, input$group_var)
    })

    observe({
      lv <- group_levels()
      req(length(lv) >= 2)
      sel_a <- isolate(input$group_a); sel_b <- isolate(input$group_b)
      if (is.null(sel_a) || !sel_a %in% lv) sel_a <- lv[2]
      if (is.null(sel_b) || !sel_b %in% lv || identical(sel_b, sel_a)) sel_b <- lv[1]
      updateSelectInput(session, "group_a", choices = lv, selected = sel_a)
      updateSelectInput(session, "group_b", choices = lv, selected = sel_b)

      sel_multi <- isolate(input$groups_included)
      sel_multi <- if (is.null(sel_multi) || !all(sel_multi %in% lv)) lv else sel_multi
      updateSelectInput(session, "groups_included", choices = lv, selected = sel_multi)
    })

    # Which matrices can be modelled depends on what each mode actually has.
    output$matrix_choice <- renderUI({
      ms <- active_modes()
      req(length(ms) > 0)
      any_batch <- any(vapply(ms, function(m) isTRUE(m$has_batch), logical(1)))
      all_combat <- all(vapply(ms, function(m) !is.null(m$batch_corrected), logical(1)))

      choices <- list()
      if (any_batch) choices[["Log2-normalized, batch as covariate"]] <- "covariate"
      choices[["Log2-normalized, no batch adjustment"]] <- "none"
      if (all_combat) choices[["ComBat-corrected, no batch term"]] <- "batch_corrected"

      selected <- isolate(input$matrix_source)
      if (is.null(selected) || !selected %in% unlist(choices)) {
        selected <- if (any_batch) "covariate" else "none"
      }
      tagList(
        radioButtons(ns("matrix_source"), "Data to model:", choices = choices, selected = selected),
        if (!any_batch) {
          p(style = "color: #666; font-size: 12px;",
            "No batch variable was supplied, so there is no batch to adjust for.")
        } else if (!all_combat) {
          p(style = "color: #666; font-size: 12px;",
            "Run the scaling step with ComBat enabled for every mode to also model the corrected matrix.")
        }
      )
    })

    selected_groups <- reactive({
      if (identical(input$comparison, "pairwise")) {
        c(input$group_a, input$group_b)
      } else {
        input$groups_included
      }
    })

    # ---- model description shown before running ------------------------
    output$model_summary <- renderText({
      ms <- active_modes()
      req(length(ms) > 0, input$group_var, input$comparison, input$matrix_source)
      grp_sel <- selected_groups()
      req(length(grp_sel) >= 2)

      counts <- vapply(ms, function(m) {
        g <- m$matched_data$full_sample_data[[input$group_var]]
        g <- g[g %in% grp_sel]
        paste0(m$label, " — ",
               paste(vapply(grp_sel, function(lv) paste0(lv, ": ", sum(g == lv)), character(1)),
                     collapse = ", "))
      }, character(1))

      use_batch_modes <- vapply(ms, function(m) {
        identical(input$matrix_source, "covariate") && isTRUE(m$has_batch)
      }, logical(1))
      design_text <- paste0("~ 0 + ", input$group_var,
                            if (any(use_batch_modes)) " + batch" else "")
      data_text <- switch(input$matrix_source,
        "batch_corrected" = "ComBat batch-corrected log2 data (batch already removed from the matrix)",
        "none" = "log2-normalized data, no batch term in the model",
        "log2-normalized data (batch effects controlled for in the model)")
      test_text <- if (identical(input$comparison, "pairwise")) {
        paste0("moderated t-test: ", input$group_a, " - ", input$group_b)
      } else {
        paste0("moderated F-test across ", length(grp_sel), " groups: ", paste(grp_sel, collapse = ", "))
      }

      paste0(
        "Model to be fitted: ", design_text, "\n",
        "Test: ", test_text, "\n",
        "Data: ", data_text, "\n",
        "Modes fitted separately: ", paste(vapply(ms, function(m) m$label, character(1)), collapse = ", "),
        if (length(ms) > 1) "  (one FDR correction across both)" else "", "\n",
        "Samples per group:\n  ", paste(counts, collapse = "\n  ")
      )
    })

    # ---- run -----------------------------------------------------------
    observeEvent(input$run_de, {
      ms <- active_modes()
      req(length(ms) > 0, input$group_var, input$comparison)
      shinyjs::disable("run_de")
      showModal(modalDialog("Fitting models...", footer = NULL))
      tryCatch({
        run <- de_run(
          modes = ms,
          group_var = input$group_var,
          comparison = input$comparison,
          group_a = input$group_a,
          group_b = input$group_b,
          groups_included = if (identical(input$comparison, "global")) input$groups_included else NULL,
          matrix_source = input$matrix_source,
          fdr_cutoff = input$fdr_cutoff,
          lfc_cutoff = if (identical(input$comparison, "pairwise")) input$lfc_cutoff else 0
        )
        run$source_fingerprint <- run_fingerprint(run, current_fingerprints())

        rv$current <- run
        key <- paste0(run$label, " [", run$matrix_source, "]")
        rv$runs[[key]] <- run

        updateSelectizeInput(session, "feature", choices = run$results$feature,
                             selected = run$results$feature[1], server = TRUE)

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

    output$results_table <- DT::renderDataTable({
      req(rv$current)
      res <- rv$current$results
      res$significant <- ifelse(res$significant, "Yes", "No")
      num_cols <- setdiff(colnames(res)[vapply(res, is.numeric, logical(1))], c("p_value", "fdr"))
      dt <- DT::datatable(res, selection = "single", rownames = FALSE,
                          options = list(scrollX = TRUE, pageLength = 15, order = list()))
      if (length(num_cols) > 0) dt <- DT::formatSignif(dt, columns = num_cols, digits = 4)
      DT::formatSignif(dt, columns = c("p_value", "fdr"), digits = 3)
    })

    observeEvent(input$results_table_rows_selected, {
      req(rv$current)
      i <- input$results_table_rows_selected
      req(length(i) == 1)
      updateSelectizeInput(session, "feature", choices = rv$current$results$feature,
                           selected = rv$current$results$feature[i], server = TRUE)
    })

    output$main_plot_note <- renderUI({
      req(rv$current)
      multi <- length(rv$current$mode_labels) > 1
      base <- if (rv$current$mode == "pairwise") {
        paste0("Coloured points clear both the FDR and fold change thresholds. ",
               "Positive log2 fold changes are higher in ", rv$current$group_a, ".")
      } else {
        "An F-test has no single fold change, so features are plotted against mean intensity. Coloured points clear the FDR threshold."
      }
      p(style = "color: #666;",
        paste0(base, if (multi) " Marker shape distinguishes the two ionization modes." else ""))
    })

    output$main_plot <- renderPlotly({
      req(rv$current)
      if (rv$current$mode == "pairwise") de_volcano_plot(rv$current) else de_significance_plot(rv$current)
    })

    output$heatmap_mode_choice <- renderUI({
      req(rv$current)
      labs <- rv$current$mode_labels
      if (length(labs) < 2) return(NULL)
      selectInput(ns("heatmap_mode"), "Ionization mode:", choices = labs, selected = labs[1])
    })

    heatmap_mode <- reactive({
      req(rv$current)
      labs <- rv$current$mode_labels
      if (length(labs) < 2) return(labs[1])
      if (is.null(input$heatmap_mode) || !input$heatmap_mode %in% labs) labs[1] else input$heatmap_mode
    })

    output$heatmap_note <- renderUI({
      req(rv$current)
      ml <- heatmap_mode()
      sel <- de_display_features(rv$current, mode_label = ml)
      multi <- length(rv$current$mode_labels) > 1
      prefix <- if (multi) "The modes have different samples, so each gets its own heatmap. " else ""
      if (sel$used_fallback) {
        p(style = "color: #666;",
          paste0(prefix, "No feature cleared the thresholds in this mode, so the strongest features by p-value are shown. Values are per-feature z-scores."))
      } else {
        p(style = "color: #666;",
          paste0(prefix, "Showing ", min(length(sel$features), DE_MAX_HEATMAP_FEATURES), " of ",
                 sel$n_available, " significant features as per-feature z-scores, with samples grouped by ",
                 rv$current$group_var, "."))
      }
    })

    output$heatmap <- renderPlotly({
      req(rv$current)
      p <- de_heatmap_plot(rv$current, mode_label = heatmap_mode())
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
      req(nrow(row) >= 1)
      row <- row[1, , drop = FALSE]
      lines <- c(paste0("Feature: ", row$feature),
                 paste0("Ionization mode: ", row$mode))
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
      paste(c(lines,
              paste0("p-value: ", signif(row$p_value, 3), "   FDR: ", signif(row$fdr, 3)),
              paste0("Mean log2 intensity: ", signif(row$mean_intensity, 4)),
              paste0("Significant at current thresholds: ", ifelse(row$significant, "yes", "no"))),
            collapse = "\n")
    })

    output$saved_table <- DT::renderDataTable({
      req(length(rv$runs) > 0)
      df <- do.call(rbind, lapply(seq_along(rv$runs), function(i) {
        run <- rv$runs[[i]]
        data.frame(
          Sheet = de_sheet_name(run, i),
          Comparison = run$label,
          Test = if (run$mode == "pairwise") "moderated t-test" else "moderated F-test",
          Modes = paste(run$mode_labels, collapse = ", "),
          Data = de_matrix_label(run),
          Thresholds = paste0("FDR <= ", run$fdr_cutoff,
                              if (run$mode == "pairwise") paste0(", |log2FC| >= ", run$lfc_cutoff) else ""),
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
