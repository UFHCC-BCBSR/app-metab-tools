# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(readxl)
library(sva)
library(dplyr)
library(tibble)
library(shinyjs)
library(impute)
options(shiny.maxRequestSize = 30*1024^3)

# ===========================
# Helper Functions
# ===========================

filter_missing <- function(count_matrix, threshold = 0.5) {
  keep <- rowMeans(is.na(count_matrix) | count_matrix == 0) <= threshold
  count_matrix[keep, ]
}

impute_knn <- function(count_matrix, k = 10) {
  result <- impute.knn(as.matrix(count_matrix), k = k)
  result$data
}

filter_iqr <- function(count_matrix, top_percent = 10) {
  iqrs <- apply(count_matrix, 1, IQR, na.rm = TRUE)
  threshold <- quantile(iqrs, top_percent / 100)  # remove bottom X%, not keep top X%
  count_matrix[iqrs >= threshold, ]
}

normalize_matrix <- function(count_matrix, row_norm = "SumNorm", 
                             trans_norm = "LogNorm", scale_norm = "ParetoNorm") {
  # Row normalization
  if (row_norm == "SumNorm") {
    col_sums <- colSums(count_matrix, na.rm = TRUE)
    count_matrix <- sweep(count_matrix, 2, col_sums, "/") * median(col_sums)
  } else if (row_norm == "MedianNorm") {
    col_medians <- apply(count_matrix, 2, median, na.rm = TRUE)
    count_matrix <- sweep(count_matrix, 2, col_medians, "/") * median(col_medians)
  }
  
  # Transformation
  if (trans_norm == "LogNorm") {
    count_matrix <- log2(count_matrix + 1)
  } else if (trans_norm == "Log10Norm") {
    count_matrix <- log10(count_matrix + 1)
  }
  
  # Scaling
  if (scale_norm == "ParetoNorm") {
    row_means <- rowMeans(count_matrix, na.rm = TRUE)
    row_sds <- apply(count_matrix, 1, sd, na.rm = TRUE)
    count_matrix <- sweep(count_matrix, 1, row_means, "-")
    count_matrix <- sweep(count_matrix, 1, sqrt(row_sds), "/")
  } else if (scale_norm == "AutoNorm") {
    row_means <- rowMeans(count_matrix, na.rm = TRUE)
    row_sds <- apply(count_matrix, 1, sd, na.rm = TRUE)
    count_matrix <- sweep(count_matrix, 1, row_means, "-")
    count_matrix <- sweep(count_matrix, 1, row_sds, "/")
  }
  
  count_matrix
}

# ===========================
# UI
# ===========================
ui <- dashboardPage(
  dashboardHeader(title = "Metabo Tools"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Preprocess Data", tabName = "configure", icon = icon("cogs")),
      menuItem("Batch Correction", tabName = "correction", icon = icon("magic")),
      menuItem("Download Results", tabName = "download", icon = icon("download"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tabItems(
      # Upload tab
      tabItem(tabName = "upload",
              fluidRow(
                box(title = "Try the Demo", status = "success", solidHeader = TRUE, width = 12,
                    p("Want to try the app with example data first?"),
                    actionButton("load_demo", "Load Test Data",
                                 class = "btn-success btn-lg",
                                 icon = icon("play")),
                    br(), br(),
                    div(id = "demo_info", style = "display: none;",
                        div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; border: 1px solid #c3e6cb;",
                            p(strong("Demo data loaded!"),
                              "This includes simulated RNA-seq count data with 2000 features, 48 samples across 3 batches and 2 biological conditions.
                Go to the Configure Data tab to proceed."))
                    ),
                    hr()
                )
              ),
              fluidRow(
                box(title = "Upload Data", status = "primary", solidHeader = TRUE, width = 6,
                    fileInput("count_file", "Choose Count Data File (.csv or .xlsx)",
                              accept = c(".csv", ".xlsx")),
                    conditionalPanel(
                      condition = "output.count_uploaded",
                      h4("Count Data Preview:"),
                      DT::dataTableOutput("count_preview")
                    )
                ),
                box(title = "Upload Sample Metadata", status = "primary", solidHeader = TRUE, width = 6,
                    fileInput("sample_file", "Choose Sample Metadata File (.csv or .xlsx)",
                              accept = c(".csv", ".xlsx")),
                    conditionalPanel(
                      condition = "output.sample_uploaded",
                      h4("Sample Data Preview:"),
                      DT::dataTableOutput("sample_preview")
                    )
                )
              )
      ),
      
      # Configure tab
      tabItem(tabName = "configure",
              fluidRow(
                box(title = "Configure Count Data", status = "warning", solidHeader = TRUE, width = 6,
                    conditionalPanel(
                      condition = "output.count_uploaded",
                      selectInput("feature_col", "Select Feature Name Column:", choices = NULL),
                      selectInput("drop_cols", "Select Columns to Drop:", choices = NULL, multiple = TRUE),
                      br(),
                      verbatimTextOutput("count_config_summary")
                    )
                ),
                box(title = "Configure Sample Data", status = "warning", solidHeader = TRUE, width = 6,
                    conditionalPanel(
                      condition = "output.sample_uploaded",
                      selectInput("batch_col", "Select Batch Variable Column:", choices = NULL),
                      selectInput("sample_name_col", "Select Sample Name Column:", choices = NULL),
                      br(),
                      checkboxInput("has_bio_var", "My data has biological variation of interest", value = TRUE),
                      conditionalPanel(
                        condition = "input.has_bio_var == true",
                        selectInput("bio_col", "Select Biological Variable Column:", choices = NULL)
                      ),
                      conditionalPanel(
                        condition = "input.has_bio_var == false",
                        div(style = "background-color: #fff3cd; padding: 10px; border-radius: 5px; border: 1px solid #ffeaa7;",
                            p(strong("Note:"), "Batch correction will be performed without preserving biological variation."))
                      ),
                      br(),
                      verbatimTextOutput("sample_config_summary")
                    )
                )
              ),
              
              # Step 1: Preview Sample Matching
              fluidRow(
                box(title = "Step 1: Preview Sample Matching", status = "info", 
                    solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.ready_to_preview",
                      actionButton("preview_matching", "Preview Sample Matching",
                                   class = "btn-info btn-lg")
                    ),
                    br(), br(),
                    conditionalPanel(
                      condition = "output.matching_previewed",
                      verbatimTextOutput("matching_summary"),
                      br(),
                      DT::dataTableOutput("matching_table")
                    )
                )
              ),
              
              # Step 2: Optional Preprocessing
              fluidRow(
                box(title = "Step 2: Optional Preprocessing", status = "warning",
                    solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.matching_previewed",
                      div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
                          h5("Select preprocessing steps to apply before batch correction.
               Uncheck any steps already performed on your data."),
                          br(),
                          div(style = "background-color: #fff3cd; padding: 10px; border-radius: 5px; border-left: 4px solid #ffc107; margin-bottom: 15px;",
                              p(strong("Important:"), "For ComBat (metabolomics/intensity data), do ", 
                                strong("not"), " apply normalization or scaling here — ComBat applies 
                  log2 transformation internally before correction. Normalization and scaling 
                  should be applied ", strong("after"), " batch correction. 
                  For ComBat-Seq (RNA-seq), do ", strong("not"), " log transform or 
                  normalize at all — it requires raw integer counts.")
                          ),
                          fluidRow(
                            column(6,
                                   div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                                       p(strong("Recommended for most datasets:")),
                                       checkboxInput("do_missing_filter",
                                                     "Filter features with >50% missing values",
                                                     value = TRUE),
                                       p(style = "color: #666; font-size: 12px; margin-left: 20px; margin-top: -10px;",
                                         "Removes features that are absent in most samples and would 
                            be unreliable for batch correction. Safe to leave checked 
                            unless you have already filtered your data."),
                                       checkboxInput("do_imputation",
                                                     "Impute remaining missing values (KNN)",
                                                     value = TRUE),
                                       p(style = "color: #666; font-size: 12px; margin-left: 20px; margin-top: -10px;",
                                         "Required if any missing values remain — ComBat and ComBat-Seq 
                            cannot handle NA values. KNN imputation estimates missing values 
                            from similar features. Leave checked unless your data is already 
                            complete.")
                                   ),
                                   br(),
                                   div(style = "background-color: #fff3cd; padding: 10px; border-radius: 5px;",
                                       p(strong("Optional — use with caution:")),
                                       checkboxInput("do_iqr_filter",
                                                     "Filter low-variance features (IQR)",
                                                     value = FALSE),
                                       p(style = "color: #666; font-size: 12px; margin-left: 20px; margin-top: -10px;",
                                         "Removes the lowest-variance features which carry little 
                            information. Not required for batch correction and is better 
                            suited as a downstream analysis step. Only check if your 
                            dataset is very large and you want to reduce noise."),
                                       conditionalPanel(
                                         condition = "input.do_iqr_filter == true",
                                         numericInput("iqr_threshold",
                                                      "Remove bottom X% lowest variance features:",
                                                      value = 10, min = 1, max = 50, step = 5),
                                         p(style = "color: #666; font-size: 12px;",
                                           "Default 10% is conservative. Values above 25% will 
                              remove a substantial number of features.")
                                       )
                                   )
                            ),
                            column(6,
                                   div(style = "background-color: #f8d7da; padding: 10px; border-radius: 5px;",
                                       p(strong("Not recommended before batch correction:")),
                                       checkboxInput("do_normalization",
                                                     "Normalize data",
                                                     value = FALSE),
                                       p(style = "color: #666; font-size: 12px; margin-left: 20px; margin-top: -10px;",
                                         "Normalization and scaling should generally be applied ", 
                                         strong("after"), " batch correction, not before. 
                            Pareto scaling produces negative values which will break 
                            ComBat-Seq, and scaling before ComBat can obscure the 
                            batch signal. Only check this if you have a specific reason 
                            to normalize first."),
                                       conditionalPanel(
                                         condition = "input.do_normalization == true",
                                         div(style = "background-color: #f8d7da; padding: 8px; border-radius: 5px; margin-bottom: 8px;",
                                             p(style = "color: #721c24; font-size: 12px;",
                                               icon("exclamation-triangle"),
                                               strong("Warning:"), " If using ComBat, do not apply 
                                  log transformation here as it will be applied 
                                  automatically. If using ComBat-Seq, do not normalize 
                                  at all — it requires raw counts.")
                                         ),
                                         selectInput("row_norm", "Row Normalization:",
                                                     choices = list(
                                                       "Sum Normalization" = "SumNorm",
                                                       "Median Normalization" = "MedianNorm",
                                                       "None" = "None"
                                                     ), selected = "None"),
                                         p(style = "color: #666; font-size: 12px;",
                                           "Sum/Median normalization corrects for differences in 
                              total sample intensity. Useful if samples have very 
                              different total abundances."),
                                         selectInput("trans_norm", "Transformation:",
                                                     choices = list(
                                                       "Log2" = "LogNorm",
                                                       "Log10" = "Log10Norm",
                                                       "None" = "None"
                                                     ), selected = "None"),
                                         p(style = "color: #666; font-size: 12px;",
                                           "ComBat applies log2 internally — do not apply here 
                              if using ComBat. ComBat-Seq does not use log transformation."),
                                         selectInput("scale_norm", "Scaling:",
                                                     choices = list(
                                                       "Pareto Scaling" = "ParetoNorm",
                                                       "Auto Scaling" = "AutoNorm",
                                                       "None" = "None"
                                                     ), selected = "None"),
                                         p(style = "color: #666; font-size: 12px;",
                                           "Pareto and Auto scaling produce negative values and 
                              will break ComBat-Seq. For ComBat, scaling before 
                              correction is not recommended — apply after.")
                                       )
                                   )
                            )
                          )
                      ),
                      actionButton("run_preprocessing", "Run Preprocessing",
                                   class = "btn-warning btn-lg"),
                      br(), br(),
                      conditionalPanel(
                        condition = "output.preprocessing_complete",
                        verbatimTextOutput("preprocessing_summary")
                      )
                    )
                )
              ),
              
              # Step 3: PCA
              fluidRow(
                box(title = "Step 3: PCA Before Batch Correction", status = "primary",
                    solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.preprocessing_complete",
                      actionButton("run_pca_before", "Run PCA",
                                   class = "btn-primary btn-lg"),
                      br(), br()
                    ),
                    conditionalPanel(
                      condition = "output.pca_before_complete",
                      plotlyOutput("pca_before", height = "400px")
                    )
                )
              )
      ),
      
      # Correction tab
      tabItem(tabName = "correction",
              fluidRow(
                box(title = "Select Batch Correction Method", status = "info", 
                    solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.pca_before_complete",
                      radioButtons("batch_method", "Choose batch correction method:",
                                   choices = list(
                                     "ComBat-Seq (recommended for RNA-seq and other count-based data)" = "combat_seq",
                                     "ComBat (recommended for metabolomics and intensity data)" = "combat_original"
                                   ),
                                   selected = "combat_seq"),
                      br(),
                      div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #007bff;",
                          h5("Method Guide:"),
                          tags$ul(
                            tags$li(strong("ComBat-Seq:"), "Uses negative binomial model, works directly with count data. Best for RNA-seq or other data with count-like properties and overdispersion."),
                            tags$li(strong("ComBat:"), "Uses log transformation + linear model. Best for metabolomics peak intensities, microarray data, or other continuous measurements.")
                          ),
                          p(strong("Note:"), "PCA plots always use log-transformed data for visualization regardless of correction method.")
                      )
                    )
                ),
                box(title = "Batch Correction", status = "success", solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.pca_before_complete",
                      actionButton("run_combat", "Run Batch Correction",
                                   class = "btn-success btn-lg"),
                      br(), br()
                    ),
                    conditionalPanel(
                      condition = "output.correction_complete",
                      h4("Data after batch correction:"),
                      plotlyOutput("pca_after", height = "400px")
                    )
                )
              )
      ),
      
      # Download tab
      tabItem(tabName = "download",
              fluidRow(
                box(title = "Download Results", status = "info", solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.correction_complete",
                      h4("Download batch-corrected data:"),
                      downloadButton("download_corrected", "Download Batch-Corrected Matrix",
                                     class = "btn-primary"),
                      br(), br(),
                      downloadButton("download_metadata", "Download Updated Sample Metadata",
                                     class = "btn-primary"),
                      br(), br(),
                      h4("Usage Notes:"),
                      tags$ul(
                        tags$li("The batch-corrected matrix should be used for visualization and clustering, but NOT for differential expression analysis"),
                        tags$li("For differential expression (DESeq2, limma-voom, edgeR), use the original raw counts with the batch factor from the metadata"),
                        tags$li("Include the batch correction factor as a covariate in your model design")
                      )
                    )
                )
              )
      )
    )
  )
)

# ===========================
# Server
# ===========================
server <- function(input, output, session) {
  values <- reactiveValues(
    count_data = NULL,
    sample_data = NULL,
    matched_data = NULL,
    preprocessed_data = NULL,
    corrected_data = NULL,
    pca_before = NULL,
    pca_after = NULL,
    matching_previewed = FALSE,
    preprocessing_complete = FALSE,
    pca_before_complete = FALSE,
    preprocessing_log = NULL
  )
  
  load_file <- function(file_path) {
    ext <- tools::file_ext(file_path)
    if (ext == "csv") return(read.csv(file_path, stringsAsFactors = FALSE))
    if (ext %in% c("xlsx", "xls")) return(read_excel(file_path))
    NULL
  }
  
  observeEvent(input$count_file, {
    req(input$count_file)
    values$count_data <- load_file(input$count_file$datapath)
    updateSelectInput(session, "feature_col", choices = colnames(values$count_data))
    updateSelectInput(session, "drop_cols", choices = colnames(values$count_data))
    values$matching_previewed <- FALSE
    values$preprocessing_complete <- FALSE
    values$pca_before_complete <- FALSE
  })
  
  observeEvent(input$sample_file, {
    req(input$sample_file)
    values$sample_data <- load_file(input$sample_file$datapath)
    updateSelectInput(session, "batch_col", choices = colnames(values$sample_data))
    updateSelectInput(session, "bio_col", choices = colnames(values$sample_data))
    updateSelectInput(session, "sample_name_col", choices = colnames(values$sample_data))
    values$matching_previewed <- FALSE
    values$preprocessing_complete <- FALSE
    values$pca_before_complete <- FALSE
  })
  
  observeEvent(input$load_demo, {
    tryCatch({
      if (file.exists("test-data/counts-data.csv") && file.exists("test-data/sample-data.csv")) {
        values$count_data <- read.csv("test-data/counts-data.csv", stringsAsFactors = FALSE)
        values$sample_data <- read.csv("test-data/sample-data.csv", stringsAsFactors = FALSE)
        updateSelectInput(session, "feature_col", choices = colnames(values$count_data), selected = "gene_id")
        updateSelectInput(session, "drop_cols", choices = colnames(values$count_data),
                          selected = c("description", "notes", "quality_flag"))
        updateSelectInput(session, "batch_col", choices = colnames(values$sample_data), selected = "batch")
        updateSelectInput(session, "bio_col", choices = colnames(values$sample_data), selected = "biological_var")
        updateSelectInput(session, "sample_name_col", choices = colnames(values$sample_data), selected = "sample_name")
        values$matching_previewed <- FALSE
        values$preprocessing_complete <- FALSE
        values$pca_before_complete <- FALSE
        updateRadioButtons(session, "batch_method", selected = "combat_seq")
        updateCheckboxInput(session, "has_bio_var", value = TRUE)
        shinyjs::show("demo_info")
        showNotification("Demo data loaded! Go to Preprocess Data tab.", type = "message")
      } else {
        showNotification("Test data files not found.", type = "error")
      }
    }, error = function(e) {
      showNotification(paste("Error loading demo data:", e$message), type = "error")
    })
  })
  
  observeEvent({
    list(input$feature_col, input$drop_cols, input$batch_col, input$bio_col,
         input$sample_name_col, input$has_bio_var)
  }, {
    values$matching_previewed <- FALSE
    values$preprocessing_complete <- FALSE
    values$pca_before_complete <- FALSE
  }, ignoreNULL = FALSE)
  
  output$count_uploaded <- reactive({ !is.null(values$count_data) })
  output$sample_uploaded <- reactive({ !is.null(values$sample_data) })
  
  output$ready_to_preview <- reactive({
    basic_ready <- !is.null(values$count_data) && !is.null(values$sample_data) &&
      !is.null(input$feature_col) && !is.null(input$batch_col) &&
      !is.null(input$sample_name_col) &&
      input$feature_col != "" && input$batch_col != "" && input$sample_name_col != ""
    if (basic_ready && !is.null(input$has_bio_var) && input$has_bio_var) {
      return(basic_ready && !is.null(input$bio_col) && input$bio_col != "")
    } else if (basic_ready && !is.null(input$has_bio_var)) {
      return(basic_ready)
    } else {
      return(FALSE)
    }
  })
  
  # ===========================
  # Step 1: Preview Sample Matching
  # ===========================
  observeEvent(input$preview_matching, {
    req(values$count_data, values$sample_data,
        input$feature_col, input$batch_col, input$sample_name_col)
    if (!is.null(input$has_bio_var) && input$has_bio_var) req(input$bio_col)
    
    tryCatch({
      count_processed <- values$count_data
      if (!is.null(input$drop_cols) && length(input$drop_cols) > 0) {
        count_processed <- count_processed[, !colnames(count_processed) %in% input$drop_cols, drop = FALSE]
      }
      
      feature_col <- input$feature_col
      potential_sample_cols <- setdiff(colnames(count_processed), feature_col)
      
      if (!is.null(input$has_bio_var) && input$has_bio_var && !is.null(input$bio_col) && input$bio_col != "") {
        sample_processed <- values$sample_data[, c(input$sample_name_col, input$batch_col, input$bio_col)]
        colnames(sample_processed) <- c("sample_name", "batch", "biological_var")
        has_bio_var <- TRUE
      } else {
        sample_processed <- values$sample_data[, c(input$sample_name_col, input$batch_col)]
        colnames(sample_processed) <- c("sample_name", "batch")
        sample_processed$biological_var <- "no_bio_var"
        has_bio_var <- FALSE
      }
      
      matched_samples <- c()
      sample_mapping <- c()
      unmatched_samples <- c()
      
      for (sample_name in sample_processed$sample_name) {
        normalized_sample <- gsub("[-.]", ".", sample_name)
        normalized_cols <- gsub("[-.]", ".", potential_sample_cols)
        match_indices <- which(normalized_cols == normalized_sample)
        if (length(match_indices) == 0) {
          match_indices <- grep(normalized_sample, normalized_cols, fixed = TRUE)
        }
        if (length(match_indices) > 0) {
          matched_samples <- c(matched_samples, potential_sample_cols[match_indices[1]])
          sample_mapping <- c(sample_mapping, sample_name)
        } else {
          unmatched_samples <- c(unmatched_samples, sample_name)
        }
      }
      
      if (length(matched_samples) > 0) {
        count_matrix <- as.matrix(count_processed[, matched_samples, drop = FALSE])
        rownames(count_matrix) <- count_processed[[feature_col]]
        sample_data_matched <- sample_processed[sample_processed$sample_name %in% sample_mapping, ]
        sample_order <- match(sample_mapping, sample_data_matched$sample_name)
        sample_data_matched <- sample_data_matched[sample_order, ]
        
        values$matched_data <- list(
          count_matrix = count_matrix,
          sample_data = sample_data_matched,
          matched_samples = matched_samples,
          unmatched_samples = unmatched_samples,
          has_bio_var = has_bio_var
        )
        
        # Reset downstream steps
        values$preprocessed_data <- NULL
        values$pca_before <- NULL
        values$preprocessing_complete <- FALSE
        values$pca_before_complete <- FALSE
        values$matching_previewed <- TRUE
        
        showNotification(paste("Matched", length(matched_samples), "samples."), type = "message")
      } else {
        values$matched_data <- NULL
        values$matching_previewed <- FALSE
        showNotification("No sample matches found! Check your sample names.", type = "error")
      }
    }, error = function(e) {
      values$matched_data <- NULL
      values$matching_previewed <- FALSE
      showNotification(paste("Error processing data:", e$message), type = "error")
      print(paste("Debug error:", e$message))
    })
  })
  
  output$matching_summary <- renderText({
    req(values$matched_data)
    unmatched <- values$matched_data$unmatched_samples
    unmatched_text <- if (length(unmatched) > 0) {
      paste0("\nUnmatched samples (", length(unmatched), "): ",
             paste(unmatched, collapse = ", "))
    } else {
      "\nAll samples matched successfully."
    }
    paste0(
      "Matched samples: ", length(values$matched_data$matched_samples), "\n",
      "Features: ", nrow(values$matched_data$count_matrix), "\n",
      "Batches: ", length(unique(values$matched_data$sample_data$batch)), "\n",
      "Biological groups: ", ifelse(
        values$matched_data$has_bio_var,
        length(unique(values$matched_data$sample_data$biological_var)),
        "None"
      ),
      unmatched_text
    )
  })
  
  output$matching_table <- DT::renderDataTable({
    req(values$matched_data)
    df <- data.frame(
      `Count Column` = values$matched_data$matched_samples,
      `Metadata Sample` = values$matched_data$sample_data$sample_name,
      `Batch` = values$matched_data$sample_data$batch,
      check.names = FALSE
    )
    if (values$matched_data$has_bio_var) {
      df$`Biological Group` <- values$matched_data$sample_data$biological_var
    }
    DT::datatable(df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # ===========================
  # Step 2: Run Preprocessing
  # ===========================
  observeEvent(input$run_preprocessing, {
    req(values$matched_data)
    
    tryCatch({
      mat <- values$matched_data$count_matrix
      log_lines <- c()
      features_start <- nrow(mat)
      
      # Convert to numeric
      mat <- apply(mat, 2, as.numeric)
      rownames(mat) <- rownames(values$matched_data$count_matrix)
      
      # Missing value filter
      if (input$do_missing_filter) {
        before <- nrow(mat)
        mat <- filter_missing(mat, threshold = 0.5)
        after <- nrow(mat)
        log_lines <- c(log_lines,
                       paste0("Missing value filter (>50%): ", before, " -> ", after,
                              " features (removed ", before - after, ")"))
      }
      
      # KNN imputation
      if (input$do_imputation) {
        na_count <- sum(is.na(mat) | mat == 0)
        mat[mat == 0] <- NA
        mat <- impute_knn(mat)
        log_lines <- c(log_lines,
                       paste0("KNN imputation: ", na_count, " missing/zero values imputed"))
      }
      
      # IQR filter
      if (input$do_iqr_filter) {
        before <- nrow(mat)
        mat <- filter_iqr(mat, top_percent = input$iqr_threshold)
        after <- nrow(mat)
        log_lines <- c(log_lines,
                       paste0("IQR filter (bottom ", input$iqr_threshold, "% removed): ", 
                              before, " -> ", after,
                              " features (removed ", before - after, ")"))
      }
      
      # Normalization
      if (input$do_normalization) {
        mat <- normalize_matrix(mat,
                                row_norm = input$row_norm,
                                trans_norm = input$trans_norm,
                                scale_norm = input$scale_norm)
        log_lines <- c(log_lines,
                       paste0("Normalization: ", input$row_norm, " / ",
                              input$trans_norm, " / ", input$scale_norm))
      }
      
      log_lines <- c(
        paste0("Features before preprocessing: ", features_start),
        log_lines,
        paste0("Features after preprocessing:  ", nrow(mat))
      )
      
      values$preprocessed_data <- mat
      values$preprocessing_log <- paste(log_lines, collapse = "\n")
      values$preprocessing_complete <- TRUE
      values$pca_before_complete <- FALSE
      
      showNotification("Preprocessing complete!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error in preprocessing:", e$message), type = "error")
      print(paste("Preprocessing error:", e$message))
    })
  })
  
  output$preprocessing_summary <- renderText({
    req(values$preprocessing_log)
    values$preprocessing_log
  })
  
  # ===========================
  # Step 3: Run PCA Before Correction
  # ===========================
  observeEvent(input$run_pca_before, {
    req(values$preprocessed_data, values$matched_data)
    tryCatch({
      values$pca_before <- perform_pca(
        values$preprocessed_data,
        values$matched_data$sample_data,
        values$matched_data$has_bio_var
      )
      values$pca_before_complete <- TRUE
      showNotification("PCA complete!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error in PCA:", e$message), type = "error")
      print(paste("PCA error:", e$message))
    })
  })
  
  output$matching_previewed <- reactive({ values$matching_previewed })
  output$preprocessing_complete <- reactive({ values$preprocessing_complete })
  output$pca_before_complete <- reactive({ values$pca_before_complete })
  output$ready_for_correction <- reactive({ values$pca_before_complete })
  output$correction_complete <- reactive({ !is.null(values$corrected_data) })
  
  outputOptions(output, "count_uploaded", suspendWhenHidden = FALSE)
  outputOptions(output, "sample_uploaded", suspendWhenHidden = FALSE)
  outputOptions(output, "ready_to_preview", suspendWhenHidden = FALSE)
  outputOptions(output, "matching_previewed", suspendWhenHidden = FALSE)
  outputOptions(output, "preprocessing_complete", suspendWhenHidden = FALSE)
  outputOptions(output, "pca_before_complete", suspendWhenHidden = FALSE)
  outputOptions(output, "ready_for_correction", suspendWhenHidden = FALSE)
  outputOptions(output, "correction_complete", suspendWhenHidden = FALSE)
  
  output$count_preview <- DT::renderDataTable({
    req(values$count_data)
    DT::datatable(values$count_data[1:min(10, nrow(values$count_data)),
                                    1:min(10, ncol(values$count_data))],
                  options = list(scrollX = TRUE))
  })
  
  output$sample_preview <- DT::renderDataTable({
    req(values$sample_data)
    DT::datatable(values$sample_data, options = list(scrollX = TRUE))
  })
  
  # ===========================
  # PCA Function
  # ===========================
  perform_pca <- function(count_matrix, sample_data, has_bio_var = TRUE) {
    count_matrix <- log2(abs(count_matrix) + 1)
    
    # Remove non-finite values
    count_matrix <- count_matrix[apply(count_matrix, 1, function(x) all(is.finite(x))), ]
    
    if (nrow(count_matrix) > 1000) {
      vars <- apply(count_matrix, 1, var, na.rm = TRUE)
      top_features <- names(sort(vars, decreasing = TRUE)[1:1000])
      count_matrix <- count_matrix[top_features, ]
    }
    
    vars <- apply(count_matrix, 1, var, na.rm = TRUE)
    count_matrix <- count_matrix[vars > 0 & !is.na(vars), ]
    
    if (nrow(count_matrix) == 0) stop("No valid features remaining for PCA.")
    
    pca_data <- t(count_matrix)
    pca_result <- prcomp(pca_data, scale. = TRUE, center = TRUE)
    
    pca_df <- data.frame(
      PC1 = pca_result$x[, 1],
      PC2 = pca_result$x[, 2],
      sample_name = colnames(count_matrix),
      batch = as.factor(sample_data$batch)
    )
    
    if (has_bio_var && !"no_bio_var" %in% sample_data$biological_var) {
      pca_df$biological_var <- as.factor(sample_data$biological_var)
    } else {
      pca_df$biological_var <- as.factor(rep("No biological variable", nrow(pca_df)))
    }
    
    variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
    list(pca_df = pca_df, variance_explained = variance_explained, has_bio_var = has_bio_var)
  }
  
  # ===========================
  # PCA Plots
  # ===========================
  create_pca_plots <- function(pca_result, title) {
    req(pca_result)
    pca_df <- pca_result$pca_df
    var_exp <- pca_result$variance_explained
    has_bio_var <- pca_result$has_bio_var
    
    # Distinct color palettes for each plot
    batch_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(pca_df$batch)))
    bio_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                    "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
    
    p_batch <- plot_ly(pca_df, x = ~PC1, y = ~PC2, 
                       color = ~batch,
                       colors = batch_colors,
                       text = ~sample_name, 
                       hovertemplate = "%{text}<extra></extra>",
                       type = "scatter", mode = "markers",
                       legendgroup = "batch",
                       showlegend = TRUE) %>%
      layout(
        title = "Colored by Batch",
        xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
        yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")),
        legend = list(title = list(text = "<b>Batch</b>"))
      )
    
    if (has_bio_var && !all(pca_df$biological_var == "No biological variable")) {
      p_bio <- plot_ly(pca_df, x = ~PC1, y = ~PC2, 
                       color = ~biological_var,
                       colors = bio_colors[1:length(unique(pca_df$biological_var))],
                       text = ~sample_name, 
                       hovertemplate = "%{text}<extra></extra>",
                       type = "scatter", mode = "markers",
                       legendgroup = "bio",
                       showlegend = TRUE) %>%
        layout(
          title = "Colored by Biological Variable",
          xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
          yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")),
          legend = list(title = list(text = "<b>Biological Variable</b>"))
        )
      
      subplot(p_bio, p_batch, nrows = 1, margin = 0.08,
              shareX = FALSE, shareY = FALSE,
              titleX = TRUE, titleY = TRUE) %>%
        layout(
          title = list(text = title, x = 0.5),
          annotations = list(
            list(x = 0.2, y = 1.05, text = "<b>Biological Variable</b>",
                 xref = "paper", yref = "paper", showarrow = FALSE,
                 font = list(size = 13)),
            list(x = 0.8, y = 1.05, text = "<b>Batch</b>",
                 xref = "paper", yref = "paper", showarrow = FALSE,
                 font = list(size = 13))
          )
        )
    } else {
      p_batch %>% layout(title = list(text = title, x = 0.5))
    }
  }
  
  output$pca_before <- renderPlotly({
    req(values$pca_before)
    create_pca_plots(values$pca_before, "Data before batch correction")
  })
  
  # ===========================
  # Configuration Summaries
  # ===========================
  output$count_config_summary <- renderText({
    req(values$count_data, input$feature_col)
    dropped <- if (is.null(input$drop_cols)) "None" else paste(input$drop_cols, collapse = ", ")
    paste0("Feature column: ", input$feature_col, "\n",
           "Dropped columns: ", dropped, "\n",
           "Note: Only columns matching sample metadata will be retained as sample data")
  })
  
  output$sample_config_summary <- renderText({
    req(values$sample_data, input$batch_col, input$sample_name_col, input$has_bio_var)
    bio_text <- if (input$has_bio_var && !is.null(input$bio_col)) {
      paste("Biological variable column:", input$bio_col)
    } else {
      "No biological variable selected"
    }
    paste0("Sample name column: ", input$sample_name_col, "\n",
           "Batch column: ", input$batch_col, "\n",
           bio_text)
  })
  
  # ===========================
  # Batch Correction
  # ===========================
  observeEvent(input$run_combat, {
    req(values$preprocessed_data, values$matched_data, input$batch_method)
    
    method_name <- ifelse(input$batch_method == "combat_seq", "ComBat-Seq", "ComBat")
    showModal(modalDialog(paste("Running", method_name, "batch correction..."), footer = NULL))
    
    tryCatch({
      count_matrix <- values$preprocessed_data
      sample_data <- values$matched_data$sample_data
      has_bio_var <- values$matched_data$has_bio_var
      batch_factor <- as.factor(sample_data$batch)
      
      if (input$batch_method == "combat_seq") {
        batch_numeric <- as.numeric(batch_factor)
        if (has_bio_var && !"no_bio_var" %in% sample_data$biological_var) {
          bio_numeric <- as.numeric(as.factor(sample_data$biological_var))
          corrected_matrix <- ComBat_seq(counts = count_matrix,
                                         batch = batch_numeric,
                                         group = bio_numeric)
        } else {
          corrected_matrix <- ComBat_seq(counts = count_matrix, batch = batch_numeric)
        }
        
      } else if (input$batch_method == "combat_original") {
        if (any(count_matrix <= 0, na.rm = TRUE)) {
          count_matrix[count_matrix <= 0] <- 0.01
        }
        log_matrix <- log2(count_matrix)
        
        mod <- if (has_bio_var && !"no_bio_var" %in% sample_data$biological_var) {
          model.matrix(~as.factor(sample_data$biological_var))
        } else {
          NULL
        }
        
        corrected_log_matrix <- ComBat(dat = log_matrix, batch = batch_factor,
                                       mod = mod, par.prior = TRUE)
        corrected_matrix <- 2^corrected_log_matrix
      }
      
      values$corrected_data <- corrected_matrix
      values$pca_after <- perform_pca(corrected_matrix, sample_data, has_bio_var)
      
      removeModal()
      showNotification(paste(method_name, "completed successfully!"), type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error in", method_name, "batch correction:", e$message), type = "error")
    })
  })
  
  output$pca_after <- renderPlotly({
    req(values$pca_after)
    create_pca_plots(values$pca_after, "Data after batch correction")
  })
  
  # ===========================
  # Download Handlers
  # ===========================
  output$download_corrected <- downloadHandler(
    filename = function() paste0("batch_corrected_data_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$corrected_data, values$matched_data)
      corrected_df <- tibble::rownames_to_column(as.data.frame(values$corrected_data), var = "feature")
      bio_var_note <- if (values$matched_data$has_bio_var) {
        "biological variation preserved"
      } else {
        "no biological variation specified"
      }
      header <- c(
        "# Batch-corrected data",
        paste("# Method:", ifelse(input$batch_method == "combat_seq", "ComBat-Seq", "ComBat")),
        paste("# Biological variable:", bio_var_note),
        paste("# Date:", Sys.time()),
        "# Use for visualization and clustering only — not for differential analysis",
        ""
      )
      writeLines(header, file)
      write.csv(corrected_df, file, row.names = FALSE, append = TRUE)
    }
  )
  
  output$download_metadata <- downloadHandler(
    filename = function() paste0("sample_metadata_", Sys.Date(), ".csv"),
    content = function(file) {
      req(values$matched_data)
      metadata_df <- values$matched_data$sample_data
      metadata_df$batch_correction_factor <- metadata_df$batch
      header <- c(
        "# Sample metadata with batch correction factor",
        paste("# Date:", Sys.time()),
        ""
      )
      writeLines(header, file)
      write.csv(metadata_df, file, row.names = FALSE, append = TRUE)
    }
  )
}

shinyApp(ui = ui, server = server)
