# ===========================
# Per-mode processing pipeline — Shiny module
# ===========================
# One instance per ionization mode. Everything from upload through scaling is
# per mode: intensities are not comparable across modes, so sum/median
# normalization, missing-value filtering and batch correction all have to be
# estimated within a mode. Only the differential analysis combines them.
#
# The UI comes in three fragments so the app keeps its stage-based tabs
# (Upload / Preprocess / Batch Correction) with one section per mode inside
# each, rather than a tab per mode.
#
# Batch handling is optional throughout: a mode may have no batch column, or
# have one and skip ComBat. The single source of truth is has_usable_batch().

pipelineUploadUI <- function(id, label) {
  ns <- NS(id)
  fluidRow(
    box(title = paste0("Upload Feature Data — ", label), status = "primary",
        solidHeader = TRUE, width = 6,
        fileInput(ns("count_file"), "Choose Feature Data File (.csv or .xlsx)",
                  accept = c(".csv", ".xlsx")),
        conditionalPanel(
          condition = "output.count_uploaded", ns = ns,
          h4("Preview:"),
          DT::dataTableOutput(ns("count_preview"))
        )
    ),
    box(title = paste0("Upload Sample Metadata — ", label), status = "primary",
        solidHeader = TRUE, width = 6,
        fileInput(ns("sample_file"), "Choose Sample Metadata File (.csv or .xlsx)",
                  accept = c(".csv", ".xlsx")),
        conditionalPanel(
          condition = "output.sample_uploaded", ns = ns,
          h4("Preview:"),
          DT::dataTableOutput(ns("sample_preview"))
        )
    )
  )
}

pipelineConfigureUI <- function(id, label) {
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = paste0("Configure Feature Data — ", label), status = "warning",
          solidHeader = TRUE, width = 6,
          conditionalPanel(
            condition = "output.count_uploaded", ns = ns,
            selectInput(ns("feature_col"), "Select Feature Name Column:", choices = NULL),
            selectInput(ns("drop_cols"), "Select Columns to Drop:", choices = NULL, multiple = TRUE),
            br(),
            verbatimTextOutput(ns("count_config_summary"))
          )
      ),
      box(title = paste0("Configure Sample Data — ", label), status = "warning",
          solidHeader = TRUE, width = 6,
          conditionalPanel(
            condition = "output.sample_uploaded", ns = ns,
            selectInput(ns("sample_name_col"), "Select Sample Name Column:", choices = NULL),
            checkboxInput(ns("has_batch_var"), "My data has a batch variable", value = TRUE),
            conditionalPanel(
              condition = "input.has_batch_var == true", ns = ns,
              selectInput(ns("batch_col"), "Select Batch Variable Column:", choices = NULL)
            ),
            checkboxInput(ns("has_bio_var"), "My data has a biological variable of interest", value = TRUE),
            conditionalPanel(
              condition = "input.has_bio_var == true", ns = ns,
              selectInput(ns("bio_col"), "Select Biological Variable Column:", choices = NULL)
            ),
            br(),
            verbatimTextOutput(ns("sample_config_summary"))
          )
      )
    ),
    fluidRow(
      box(title = paste0("Step 1: Preview Sample Matching — ", label), status = "info",
          solidHeader = TRUE, width = 12,
          conditionalPanel(
            condition = "output.ready_to_preview", ns = ns,
            actionButton(ns("preview_matching"), "Preview Sample Matching",
                         class = "btn-info btn-lg")
          ),
          br(), br(),
          conditionalPanel(
            condition = "output.matching_previewed", ns = ns,
            verbatimTextOutput(ns("matching_summary")),
            br(),
            DT::dataTableOutput(ns("matching_table"))
          )
      )
    ),
    fluidRow(
      box(title = paste0("Step 2: Preprocessing — ", label), status = "warning",
          solidHeader = TRUE, width = 12,
          conditionalPanel(
            condition = "output.matching_previewed", ns = ns,
            fluidRow(
              column(6,
                     div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                         p(strong("Filtering and imputation:")),
                         checkboxInput(ns("do_missing_filter"),
                                       "Filter features with >50% missing values", value = TRUE),
                         checkboxInput(ns("do_imputation"),
                                       "Impute remaining missing values (KNN)", value = TRUE),
                         checkboxInput(ns("do_iqr_filter"),
                                       "Filter low-variance features (IQR)", value = FALSE),
                         conditionalPanel(
                           condition = "input.do_iqr_filter == true", ns = ns,
                           numericInput(ns("iqr_threshold"),
                                        "Remove bottom X% lowest variance features:",
                                        value = 10, min = 1, max = 50, step = 5)
                         )
                     )
              ),
              column(6,
                     div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px;",
                         p(strong("Sample normalization:")),
                         selectInput(ns("row_norm"), "Method:",
                                     choices = list(
                                       "Sum Normalization" = "SumNorm",
                                       "Median Normalization" = "MedianNorm",
                                       "None" = "None"
                                     ), selected = "SumNorm"),
                         p(style = "color: #666; font-size: 12px;",
                           "Log2 transformation will always be applied after sample normalization.")
                     )
              )
            ),
            br(),
            actionButton(ns("run_preprocessing"), "Run Preprocessing", class = "btn-warning btn-lg"),
            br(), br(),
            conditionalPanel(
              condition = "output.preprocessing_complete", ns = ns,
              verbatimTextOutput(ns("preprocessing_summary"))
            )
          )
      )
    ),
    fluidRow(
      box(title = paste0("Step 3: PCA Before Batch Correction — ", label), status = "primary",
          solidHeader = TRUE, width = 12,
          conditionalPanel(
            condition = "output.preprocessing_complete", ns = ns,
            actionButton(ns("run_pca_before"), "Run PCA", class = "btn-primary btn-lg"),
            br(), br()
          ),
          conditionalPanel(
            condition = "output.pca_before_complete", ns = ns,
            uiOutput(ns("pca_before_note")),
            plotlyOutput(ns("pca_before"), height = "450px")
          )
      )
    )
  )
}

pipelineCorrectionUI <- function(id, label) {
  ns <- NS(id)
  fluidRow(
    box(title = paste0("Batch Correction & Scaling — ", label), status = "success",
        solidHeader = TRUE, width = 12,
        conditionalPanel(
          condition = "!output.preprocessing_complete", ns = ns,
          p("Run preprocessing for this mode first.")
        ),
        conditionalPanel(
          condition = "output.preprocessing_complete", ns = ns,
          div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
              uiOutput(ns("combat_control")),
              selectInput(ns("scale_norm"), "Scaling method:",
                          choices = list(
                            "Pareto Scaling" = "ParetoNorm",
                            "Auto Scaling" = "AutoNorm",
                            "None" = "None"
                          ), selected = "ParetoNorm"),
              p(style = "color: #666; font-size: 12px;",
                "Scaling is applied to the normalized data, and to the batch-corrected data as well when ComBat is run. Scaling affects PCA, clustering and heatmaps; differential analysis is always modeled on the unscaled log2 data.")
          ),
          actionButton(ns("run_correction"), "Apply", class = "btn-success btn-lg"),
          br(), br(),
          conditionalPanel(
            condition = "output.correction_complete", ns = ns,
            verbatimTextOutput(ns("correction_summary")),
            br(),
            conditionalPanel(
              condition = "output.combat_ran", ns = ns,
              h4("PCA after batch correction and scaling:"),
              plotlyOutput(ns("pca_after"), height = "450px")
            )
          )
        )
    )
  )
}

# `mode_label` is a reactive holding the display name for this mode ("Positive",
# "Negative", ...); it names the mode in notifications, plots, results and the
# report. `preset` is a reactive carrying a demo dataset to load (or NULL) — it
# is how the app's demo buttons drive a mode without the module knowing about
# them.
pipelineServer <- function(id, mode_label, preset = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    values <- reactiveValues(
      count_data = NULL, sample_data = NULL, matched_data = NULL,
      preprocessed_data = NULL, normalized_data = NULL, normalized_scaled_data = NULL,
      batch_corrected_data = NULL, batch_corrected_scaled_data = NULL,
      pca_before = NULL, pca_after = NULL,
      matching_previewed = FALSE, preprocessing_complete = FALSE,
      pca_before_complete = FALSE, correction_complete = FALSE, combat_ran = FALSE,
      preprocessing_log = NULL, correction_log = NULL, processing_params = NULL
    )

    # Invalidate everything downstream of the stage that just changed. Note
    # "preprocess" means "preprocessing has just been redone", so it clears the
    # correction outputs only — never the matrices preprocessing just wrote.
    reset_downstream <- function(from = c("match", "preprocess")) {
      from <- match.arg(from)
      if (from == "match") {
        values$matching_previewed <- FALSE
        values$preprocessed_data <- NULL
        values$normalized_data <- NULL
        values$pca_before <- NULL
        values$preprocessing_complete <- FALSE
        values$pca_before_complete <- FALSE
      }
      values$normalized_scaled_data <- NULL
      values$batch_corrected_data <- NULL
      values$batch_corrected_scaled_data <- NULL
      values$pca_after <- NULL
      values$correction_complete <- FALSE
      values$combat_ran <- FALSE
    }

    load_file <- function(file_path) {
      ext <- tools::file_ext(file_path)
      if (ext == "csv") return(read.csv(file_path, stringsAsFactors = FALSE))
      if (ext %in% c("xlsx", "xls")) return(as.data.frame(read_excel(file_path)))
      NULL
    }

    observeEvent(input$count_file, {
      req(input$count_file)
      values$count_data <- load_file(input$count_file$datapath)
      updateSelectInput(session, "feature_col", choices = colnames(values$count_data))
      updateSelectInput(session, "drop_cols", choices = colnames(values$count_data))
      reset_downstream("match")
    })

    observeEvent(input$sample_file, {
      req(input$sample_file)
      values$sample_data <- load_file(input$sample_file$datapath)
      updateSelectInput(session, "batch_col", choices = colnames(values$sample_data))
      updateSelectInput(session, "bio_col", choices = colnames(values$sample_data))
      updateSelectInput(session, "sample_name_col", choices = colnames(values$sample_data))
      reset_downstream("match")
    })

    observeEvent(preset(), {
      p <- preset()
      req(p)
      if (!file.exists(p$feature_file) || !file.exists(p$sample_file)) {
        showNotification(paste0("Test data files not found for ", mode_label(), "."), type = "error")
        return()
      }
      values$count_data  <- read.csv(p$feature_file, stringsAsFactors = FALSE)
      values$sample_data <- read.csv(p$sample_file, stringsAsFactors = FALSE)
      feature_col <- if (is.null(p$feature_col)) colnames(values$count_data)[1] else p$feature_col
      updateSelectInput(session, "feature_col", choices = colnames(values$count_data), selected = feature_col)
      updateSelectInput(session, "drop_cols", choices = colnames(values$count_data), selected = p$drop_cols)
      updateSelectInput(session, "sample_name_col", choices = colnames(values$sample_data), selected = p$sample_name_col)
      updateSelectInput(session, "batch_col", choices = colnames(values$sample_data), selected = p$batch_col)
      updateSelectInput(session, "bio_col", choices = colnames(values$sample_data), selected = p$bio_col)
      updateCheckboxInput(session, "has_batch_var", value = !is.null(p$batch_col))
      updateCheckboxInput(session, "has_bio_var", value = !is.null(p$bio_col))
      reset_downstream("match")
    }, ignoreNULL = TRUE)

    output$count_uploaded <- reactive({ !is.null(values$count_data) })
    output$sample_uploaded <- reactive({ !is.null(values$sample_data) })
    output$matching_previewed <- reactive({ values$matching_previewed })
    output$preprocessing_complete <- reactive({ values$preprocessing_complete })
    output$pca_before_complete <- reactive({ values$pca_before_complete })
    output$correction_complete <- reactive({ values$correction_complete })
    output$combat_ran <- reactive({ values$combat_ran })
    output$ready_to_preview <- reactive({
      basic <- !is.null(values$count_data) && !is.null(values$sample_data) &&
        isTruthy(input$feature_col) && isTruthy(input$sample_name_col)
      if (!basic) return(FALSE)
      if (isTRUE(input$has_batch_var) && !isTruthy(input$batch_col)) return(FALSE)
      if (isTRUE(input$has_bio_var) && !isTruthy(input$bio_col)) return(FALSE)
      TRUE
    })
    for (nm in c("count_uploaded", "sample_uploaded", "matching_previewed",
                 "preprocessing_complete", "pca_before_complete",
                 "correction_complete", "combat_ran", "ready_to_preview")) {
      outputOptions(output, nm, suspendWhenHidden = FALSE)
    }

    observeEvent(input$preview_matching, {
      req(values$count_data, values$sample_data, input$feature_col, input$sample_name_col)
      tryCatch({
        matched <- match_samples(
          count_data = values$count_data,
          sample_data = values$sample_data,
          feature_col = input$feature_col,
          sample_name_col = input$sample_name_col,
          batch_col = if (isTRUE(input$has_batch_var)) input$batch_col else NULL,
          bio_col = if (isTRUE(input$has_bio_var)) input$bio_col else NULL,
          drop_cols = input$drop_cols
        )
        if (is.null(matched)) {
          values$matched_data <- NULL
          reset_downstream("match")
          showNotification(paste0(mode_label(), ": no sample matches found. Check your sample names."),
                           type = "error")
          return()
        }
        values$matched_data <- matched
        reset_downstream("match")
        values$matching_previewed <- TRUE
        showNotification(paste0(mode_label(), ": matched ", length(matched$matched_samples), " samples."),
                         type = "message")
      }, error = function(e) {
        values$matched_data <- NULL
        reset_downstream("match")
        showNotification(paste("Error processing data:", e$message), type = "error")
      })
    })

    observeEvent(input$run_preprocessing, {
      req(values$matched_data)
      tryCatch({
        mat <- values$matched_data$count_matrix
        log_lines <- c()
        n_original <- nrow(mat)
        n_imputed <- 0
        n_after_missing <- n_original

        mat <- apply(mat, 2, as.numeric)
        rownames(mat) <- rownames(values$matched_data$count_matrix)

        if (isTRUE(input$do_missing_filter)) {
          before <- nrow(mat)
          mat <- filter_missing(mat, threshold = 0.5)
          n_after_missing <- nrow(mat)
          log_lines <- c(log_lines, paste0("Missing value filter: ", before, " -> ", nrow(mat),
                                           " features (removed ", before - nrow(mat), ")"))
        }
        if (isTRUE(input$do_imputation)) {
          n_imputed <- sum(is.na(mat) | mat == 0)
          mat[mat == 0] <- NA
          mat <- impute_knn(mat)
          log_lines <- c(log_lines, paste0("KNN imputation: ", n_imputed, " values imputed"))
        }
        if (isTRUE(input$do_iqr_filter)) {
          before <- nrow(mat)
          mat <- filter_iqr(mat, top_percent = input$iqr_threshold)
          log_lines <- c(log_lines, paste0("IQR filter: ", before, " -> ", nrow(mat),
                                           " features (removed ", before - nrow(mat), ")"))
        }
        n_after_iqr <- nrow(mat)

        values$preprocessed_data <- mat

        mat_norm <- sample_normalize(mat, input$row_norm)
        if (any(mat_norm <= 0, na.rm = TRUE)) mat_norm[mat_norm <= 0] <- 0.01
        mat_norm <- log2(mat_norm)
        values$normalized_data <- mat_norm

        log_lines <- c(paste0("Input features: ", n_original), log_lines,
                       paste0("Sample normalization: ", input$row_norm),
                       "Log2 transformation applied",
                       paste0("Final features: ", nrow(mat_norm)))

        values$processing_params <- list(
          mode_label = mode_label(),
          do_missing_filter = isTRUE(input$do_missing_filter),
          do_imputation = isTRUE(input$do_imputation),
          do_iqr_filter = isTRUE(input$do_iqr_filter),
          iqr_threshold = if (isTRUE(input$do_iqr_filter)) input$iqr_threshold else NA,
          row_norm = input$row_norm,
          scale_norm = "None",
          n_original = n_original,
          n_after_missing = n_after_missing,
          n_after_iqr = n_after_iqr,
          n_imputed = n_imputed,
          n_final = nrow(mat_norm),
          n_samples = ncol(mat_norm),
          bio_var = if (isTRUE(values$matched_data$has_bio_var)) input$bio_col else "None",
          batch_var = if (has_usable_batch(values$matched_data)) input$batch_col else "None",
          has_batch = has_usable_batch(values$matched_data),
          combat_applied = FALSE
        )

        values$preprocessing_log <- paste(log_lines, collapse = "\n")
        reset_downstream("preprocess")
        values$preprocessing_complete <- TRUE
        showNotification(paste0(mode_label(), ": preprocessing complete."), type = "message")
      }, error = function(e) {
        showNotification(paste("Error in preprocessing:", e$message), type = "error")
      })
    })

    observeEvent(input$run_pca_before, {
      req(values$normalized_data, values$matched_data)
      tryCatch({
        values$pca_before <- perform_pca(
          values$normalized_data,
          values$matched_data$sample_data,
          has_bio_var = isTRUE(values$matched_data$has_bio_var),
          bio_var_name = if (isTRUE(values$matched_data$has_bio_var)) input$bio_col else "Biological Variable",
          has_batch = has_usable_batch(values$matched_data)
        )
        values$pca_before_complete <- TRUE
        showNotification(paste0(mode_label(), ": PCA complete."), type = "message")
      }, error = function(e) {
        showNotification(paste("Error in PCA:", e$message), type = "error")
      })
    })

    # ComBat is only offered when there is a batch variable with more than one
    # level, and even then it can be declined: a dataset with no batch effect
    # is better left uncorrected than pushed through an unnecessary model.
    output$combat_control <- renderUI({
      ns <- session$ns
      if (has_usable_batch(values$matched_data)) {
        tagList(
          checkboxInput(ns("do_combat"),
                        "Remove batch effects with ComBat", value = TRUE),
          p(style = "color: #666; font-size: 12px;",
            "Leave this unchecked if the PCA shows no batch effect. Differential analysis can still include batch as a covariate.")
        )
      } else {
        tagList(
          p(strong("No batch variable for this mode."),
            " Batch correction does not apply; scaling will still be applied."),
          p(style = "color: #666; font-size: 12px;",
            if (is.null(values$matched_data) || !isTRUE(values$matched_data$has_batch)) {
              "No batch column was selected on the Preprocess Data tab."
            } else {
              "The batch column has only one level among the matched samples."
            })
        )
      }
    })

    observeEvent(input$run_correction, {
      req(values$normalized_data, values$matched_data, input$scale_norm)
      run_combat <- has_usable_batch(values$matched_data) && isTRUE(input$do_combat)
      shinyjs::disable("run_correction")
      if (run_combat) showModal(modalDialog("Running ComBat batch correction...", footer = NULL))
      tryCatch({
        mat_norm <- values$normalized_data
        sample_data <- values$matched_data$sample_data
        log_lines <- character(0)

        values$normalized_scaled_data <- scale_matrix(mat_norm, input$scale_norm)

        if (run_combat) {
          batch_factor <- as.factor(sample_data$batch)
          mod <- if (isTRUE(values$matched_data$has_bio_var) &&
                     !"no_bio_var" %in% sample_data$biological_var) {
            model.matrix(~as.factor(sample_data$biological_var))
          } else {
            NULL
          }
          mat_corrected <- ComBat(dat = mat_norm, batch = batch_factor, mod = mod, par.prior = TRUE)
          values$batch_corrected_data <- mat_corrected
          values$batch_corrected_scaled_data <- scale_matrix(mat_corrected, input$scale_norm)
          values$pca_after <- perform_pca(
            values$batch_corrected_scaled_data, sample_data,
            has_bio_var = isTRUE(values$matched_data$has_bio_var),
            bio_var_name = if (isTRUE(values$matched_data$has_bio_var)) isolate(input$bio_col) else "Biological Variable",
            has_batch = TRUE
          )
          log_lines <- c(log_lines, "ComBat batch correction applied")
        } else {
          values$batch_corrected_data <- NULL
          values$batch_corrected_scaled_data <- NULL
          values$pca_after <- NULL
          log_lines <- c(log_lines, if (has_usable_batch(values$matched_data)) {
            "ComBat batch correction skipped by request"
          } else {
            "ComBat batch correction not applicable (no batch variable)"
          })
        }
        log_lines <- c(log_lines, paste0("Scaling: ", input$scale_norm))

        values$processing_params$scale_norm <- input$scale_norm
        values$processing_params$combat_applied <- run_combat
        values$correction_log <- paste(log_lines, collapse = "\n")
        values$combat_ran <- run_combat
        values$correction_complete <- TRUE

        removeModal()
        shinyjs::enable("run_correction")
        showNotification(paste0(mode_label(), ": ",
                                if (run_combat) "batch correction and scaling complete." else "scaling complete."),
                         type = "message")
      }, error = function(e) {
        removeModal()
        shinyjs::enable("run_correction")
        showNotification(paste("Error in batch correction:", e$message), type = "error")
      })
    }, ignoreInit = TRUE)

    # ---- outputs -------------------------------------------------------
    output$count_preview <- DT::renderDataTable({
      req(values$count_data)
      DT::datatable(values$count_data[seq_len(min(10, nrow(values$count_data))),
                                      seq_len(min(10, ncol(values$count_data)))],
                    options = list(scrollX = TRUE))
    })

    output$sample_preview <- DT::renderDataTable({
      req(values$sample_data)
      DT::datatable(values$sample_data, options = list(scrollX = TRUE))
    })

    output$count_config_summary <- renderText({
      req(values$count_data, input$feature_col)
      dropped <- if (is.null(input$drop_cols) || length(input$drop_cols) == 0) "None" else paste(input$drop_cols, collapse = ", ")
      paste0("Feature column: ", input$feature_col, "\n", "Dropped columns: ", dropped)
    })

    output$sample_config_summary <- renderText({
      req(values$sample_data, input$sample_name_col)
      batch_text <- if (isTRUE(input$has_batch_var) && isTruthy(input$batch_col)) {
        paste("Batch column:", input$batch_col)
      } else {
        "No batch variable selected"
      }
      bio_text <- if (isTRUE(input$has_bio_var) && isTruthy(input$bio_col)) {
        paste("Biological variable column:", input$bio_col)
      } else {
        "No biological variable selected"
      }
      paste0("Sample name column: ", input$sample_name_col, "\n", batch_text, "\n", bio_text)
    })

    output$matching_summary <- renderText({
      req(values$matched_data)
      md <- values$matched_data
      unmatched <- md$unmatched_samples
      unmatched_text <- if (length(unmatched) > 0) {
        paste0("\nUnmatched samples (", length(unmatched), "): ", paste(unmatched, collapse = ", "))
      } else {
        "\nAll samples matched successfully."
      }
      paste0(
        "Matched samples: ", length(md$matched_samples), "\n",
        "Features: ", nrow(md$count_matrix), "\n",
        "Batches: ", if (has_usable_batch(md)) length(unique(md$sample_data$batch)) else "None", "\n",
        "Biological groups: ", if (isTRUE(md$has_bio_var)) length(unique(md$sample_data$biological_var)) else "None",
        if (isTRUE(md$n_duplicate_features > 0)) paste0(
          "\nDuplicate feature names: ", md$n_duplicate_features,
          " (suffixed to keep them distinct)") else "",
        unmatched_text
      )
    })

    output$matching_table <- DT::renderDataTable({
      req(values$matched_data)
      md <- values$matched_data
      df <- data.frame(`Count Column` = md$matched_samples,
                       `Metadata Sample` = md$sample_data$sample_name,
                       check.names = FALSE)
      if (has_usable_batch(md)) df$Batch <- md$sample_data$batch
      if (isTRUE(md$has_bio_var)) df$`Biological Group` <- md$sample_data$biological_var
      DT::datatable(df, options = list(scrollX = TRUE, pageLength = 10))
    })

    output$preprocessing_summary <- renderText({ req(values$preprocessing_log); values$preprocessing_log })
    output$correction_summary <- renderText({ req(values$correction_log); values$correction_log })

    output$pca_before_note <- renderUI({
      req(values$pca_before)
      if (isTRUE(values$pca_before$has_batch)) {
        p(style = "color: #666;",
          "If the batch panel shows no separation, there is no batch effect to remove and ComBat can be skipped on the next tab.")
      } else {
        p(style = "color: #666;", "No batch variable for this mode, so samples are colored by biological group only.")
      }
    })

    output$pca_before <- renderPlotly({
      req(values$pca_before)
      create_pca_plots(values$pca_before, paste0(mode_label(), " — before batch correction"))
    })

    output$pca_after <- renderPlotly({
      req(values$pca_after)
      create_pca_plots(values$pca_after, paste0(mode_label(), " — after batch correction and scaling"))
    })

    # ---- exported state ------------------------------------------------
    reactive({
      if (is.null(values$normalized_data) || is.null(values$matched_data)) return(NULL)
      list(
        label = mode_label(),
        matched_data = values$matched_data,
        preprocessed = values$preprocessed_data,
        normalized = values$normalized_data,
        normalized_scaled = values$normalized_scaled_data,
        batch_corrected = values$batch_corrected_data,
        batch_corrected_scaled = values$batch_corrected_scaled_data,
        pca_before = values$pca_before,
        pca_after = values$pca_after,
        params = values$processing_params,
        preprocessing_log = values$preprocessing_log,
        correction_log = values$correction_log,
        has_batch = has_usable_batch(values$matched_data),
        combat_ran = values$combat_ran,
        correction_complete = values$correction_complete
      )
    })
  })
}
