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
library(openxlsx)
library(rmarkdown)
library(ggplot2)
library(RColorBrewer)
library(limma)

options(shiny.maxRequestSize = 30*1024^3)

# Modules. Sourced explicitly rather than relying on auto-loading so the load
# order is visible and the app behaves the same however it is launched.
source("modules/report_utils.R")
source("modules/preprocessing.R")
source("modules/report.R")
source("modules/de_analysis.R")
source("modules/mod_pipeline.R")
source("modules/mod_differential.R")

# ===========================
# Demo datasets
# ===========================
# Each entry describes one or both ionization modes. Two-mode datasets carry a
# `neg` block; single-mode ones do not.
DEMO_DATASETS <- list(
  simulated = list(
    pos = list(feature_file = "test-data/peak-data.csv",
               sample_file = "test-data/sample-data.csv",
               feature_col = "feature_id", batch_col = "batch",
               bio_col = "biological_var", sample_name_col = "sample_name"),
    labels = c(pos = "Simulated")
  ),
  meningioma = list(
    pos = list(feature_file = "test-data/Meningiomas85-Met-POS-BatchEffect_feature_data.csv",
               sample_file = "test-data/Meningiomas85-Met-POS-BatchEffect_sample_data.csv",
               feature_col = NULL, batch_col = "Batch",
               bio_col = "Group", sample_name_col = "Sample"),
    neg = list(feature_file = "test-data/Meningiomas85-Met-NEG-BatchEffect_feature_data.csv",
               sample_file = "test-data/Meningiomas85-Met-NEG-BatchEffect_sample_data.csv",
               feature_col = NULL, batch_col = "Batch",
               bio_col = "Group", sample_name_col = "Sample"),
    labels = c(pos = "Positive", neg = "Negative")
  ),
  procmet = list(
    pos = list(feature_file = "test-data/PROCMET_Batch1and2_RP_POS_feature_data.csv",
               sample_file = "test-data/PROCMET_Batch1and2_RP_POS_sample_data.csv",
               feature_col = "feature_id", batch_col = "Batch",
               bio_col = "Group", sample_name_col = "Sample"),
    neg = list(feature_file = "test-data/PROCMET_Batch1and2_RP_NEG_feature_data.csv",
               sample_file = "test-data/PROCMET_Batch1and2_RP_NEG_sample_data.csv",
               feature_col = "feature_id", batch_col = "Batch",
               bio_col = "Group", sample_name_col = "Sample"),
    labels = c(pos = "Positive", neg = "Negative")
  )
)

# ===========================
# UI
# ===========================
ui <- dashboardPage(
  dashboardHeader(title = "Metabo Tools"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Upload Data", tabName = "upload", icon = icon("upload")),
      menuItem("Preprocess Data", tabName = "configure", icon = icon("cogs")),
      menuItem("Batch Correction & Scaling", tabName = "correction", icon = icon("magic")),
      menuItem("Differential Analysis", tabName = "differential", icon = icon("flask")),
      menuItem("Download Results", tabName = "download", icon = icon("download")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  dashboardBody(
    useShinyjs(),
    tabItems(

      # Upload tab
      tabItem(tabName = "upload",
              fluidRow(
                box(title = "Try the Demo", status = "success", solidHeader = TRUE, width = 12,
                    p("Load one of the example datasets below to explore the app before uploading your own data."),
                    fluidRow(
                      column(4,
                             div(style = "display: flex; align-items: center; gap: 8px;",
                                 actionButton("load_demo", "Load Simulated Data",
                                              class = "btn-success btn-lg", icon = icon("play")),
                                 actionButton("info_simulated", icon("info-circle"),
                                              class = "btn-default btn-sm", title = "About this dataset")
                             )
                      ),
                      column(4,
                             div(style = "display: flex; align-items: center; gap: 8px;",
                                 actionButton("load_meningioma", "Load Meningioma Data (POS + NEG)",
                                              class = "btn-info btn-lg", icon = icon("play")),
                                 actionButton("info_meningioma", icon("info-circle"),
                                              class = "btn-default btn-sm", title = "About this dataset")
                             )
                      ),
                      column(4,
                             div(style = "display: flex; align-items: center; gap: 8px;",
                                 actionButton("load_procmet", "Load PROCMET Data (POS + NEG)",
                                              class = "btn-warning btn-lg", icon = icon("play")),
                                 actionButton("info_procmet", icon("info-circle"),
                                              class = "btn-default btn-sm", title = "About this dataset")
                             )
                      )
                    ),
                    br(),
                    div(id = "demo_info", style = "display: none;",
                        div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; border: 1px solid #c3e6cb;",
                            p(strong("Demo data loaded!"), "Go to the Preprocess Data tab to proceed."))
                    )
                )
              ),
              fluidRow(
                box(title = "Ionization Modes", status = "primary", solidHeader = TRUE, width = 12,
                    p("Positive and negative mode acquisitions are separate experiments: intensities are not comparable between them, so each is normalized, filtered and batch corrected on its own. They are brought together only for differential analysis, where a single FDR correction is applied across both."),
                    fluidRow(
                      column(4, textInput("label_pos", "Name for the first mode:", value = "Positive")),
                      column(4,
                             br(),
                             checkboxInput("use_second_mode", "I have a second ionization mode", value = FALSE)),
                      column(4,
                             conditionalPanel(
                               condition = "input.use_second_mode == true",
                               textInput("label_neg", "Name for the second mode:", value = "Negative")
                             ))
                    )
                )
              ),
              h3("Mode 1"),
              pipelineUploadUI("pos", "Mode 1"),
              conditionalPanel(
                condition = "input.use_second_mode == true",
                h3("Mode 2"),
                pipelineUploadUI("neg", "Mode 2")
              )
      ),

      # Configure tab
      tabItem(tabName = "configure",
              h3(textOutput("configure_heading_pos", inline = TRUE)),
              pipelineConfigureUI("pos", "Mode 1"),
              conditionalPanel(
                condition = "input.use_second_mode == true",
                hr(),
                h3(textOutput("configure_heading_neg", inline = TRUE)),
                pipelineConfigureUI("neg", "Mode 2")
              )
      ),

      # Batch Correction & Scaling tab
      tabItem(tabName = "correction",
              fluidRow(
                box(title = "About this step", status = "info", solidHeader = TRUE, width = 12,
                    p("Batch correction is optional. If the PCA on the previous tab shows no separation by batch, there is nothing to remove — skip ComBat and apply scaling only. Differential analysis can still account for batch by including it as a covariate, which is the safer choice when a batch effect is small or uncertain."),
                    p(style = "color: #666; font-size: 12px;",
                      "Scaling is applied here either way, because the scaled matrices are what PCA, clustering and heatmaps use.")
                )
              ),
              pipelineCorrectionUI("pos", "Mode 1"),
              conditionalPanel(
                condition = "input.use_second_mode == true",
                pipelineCorrectionUI("neg", "Mode 2")
              )
      ),

      # Differential Analysis tab
      tabItem(tabName = "differential",
              differentialUI("de")
      ),

      # About tab
      tabItem(tabName = "about",
              fluidRow(
                box(title = "About Metabo Tools", status = "primary", solidHeader = TRUE, width = 12,
                    p(strong("Metabo Tools"), "is a web application for metabolomics data preprocessing, batch correction, and differential analysis."),
                    p(tags$em("This tool is currently under active development and testing. Please use results with appropriate caution and report any issues.")),
                    br(),
                    h4("Source Code"),
                    p("The source code is available on GitHub:"),
                    p(tags$a(href = "https://github.com/UFHCC-BCBSR/app-metab-tools",
                             "https://github.com/UFHCC-BCBSR/app-metab-tools", target = "_blank")),
                    br(),
                    h4("Contact"),
                    p("For questions, feedback, or to report issues, contact:",
                      tags$a(href = "mailto:hkates@ufl.edu", "hkates@ufl.edu")),
                    br(),
                    h4("Developed in partnership with SECIM"),
                    p("This tool was developed in partnership with the ",
                      tags$a(href = "https://secim.ufl.edu/",
                             "University of Florida Southeast Center for Integrated Metabolomics (SECIM)",
                             target = "_blank"), "."),
                    br(),
                    div(style = "font-size: 12px; color: #888;",
                        p("University of Florida Health Cancer Center — Biostatistics, Computational Biology, and Bioinformatics Shared Resource (BCBSR)"))
                )
              )
      ),

      # Download tab
      tabItem(tabName = "download",
              fluidRow(
                box(title = "Download Results", status = "info", solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.downloads_ready",
                      p("Download your complete results package. The Excel file contains every processed data version for each ionization mode, a README sheet, and one sheet per differential analysis comparison. The HTML report contains a methods paragraph, processing log, file guide and PCA plots for each mode, plus a full section for every comparison you ran."),
                      br(),
                      downloadButton("download_excel", "Download Excel Data Package", class = "btn-primary btn-lg"),
                      br(), br(),
                      downloadButton("download_html", "Download HTML Report", class = "btn-info btn-lg")
                    ),
                    conditionalPanel(
                      condition = "!output.downloads_ready",
                      p("Complete preprocessing and the scaling step for at least one mode to enable downloads.")
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

  # ---- demo datasets ---------------------------------------------------
  preset_pos <- reactiveVal(NULL)
  preset_neg <- reactiveVal(NULL)

  load_demo_dataset <- function(key) {
    ds <- DEMO_DATASETS[[key]]
    has_neg <- !is.null(ds$neg)
    updateCheckboxInput(session, "use_second_mode", value = has_neg)
    updateTextInput(session, "label_pos", value = unname(ds$labels["pos"]))
    if (has_neg) updateTextInput(session, "label_neg", value = unname(ds$labels["neg"]))
    preset_pos(ds$pos)
    preset_neg(if (has_neg) ds$neg else NULL)
    shinyjs::show("demo_info")
    showNotification(
      paste0("Data loaded", if (has_neg) " for both modes" else "", ". Go to the Preprocess Data tab."),
      type = "message")
  }

  observeEvent(input$load_demo, load_demo_dataset("simulated"))
  observeEvent(input$load_meningioma, load_demo_dataset("meningioma"))
  observeEvent(input$load_procmet, load_demo_dataset("procmet"))

  observeEvent(input$info_simulated, {
    showModal(modalDialog(
      title = "Simulated Dataset",
      p("60 simulated samples across 3 batches with 2 biological groups (Control and Treatment) and 1000 features, single mode."),
      p("Batch effects were spiked in as large additive shifts on the log2 scale to ensure clear separation by batch before correction. A biological effect was spiked into a subset of 300 features with variable effect sizes."),
      p("Effect sizes are larger than typically observed in real metabolomics data — this dataset is for demonstration only."),
      easyClose = TRUE, footer = modalButton("Close")))
  })

  observeEvent(input$info_meningioma, {
    showModal(modalDialog(
      title = "Meningioma Dataset (positive and negative mode)",
      p("Real untargeted metabolomics peak intensity tables acquired in both positive and negative ionization mode."),
      p("85 samples across two batches and three levels of a biological variable, making it suitable for testing batch correction with a biological covariate."),
      p("Feature columns contain raw peak intensities. Missing values are present and will be handled during preprocessing."),
      easyClose = TRUE, footer = modalButton("Close")))
  })

  observeEvent(input$info_procmet, {
    showModal(modalDialog(
      title = "PROCMET Dataset (positive and negative mode)",
      p("Real reverse-phase untargeted metabolomics data in both ionization modes: 189 injections per mode, 3160 positive and 1249 negative features."),
      p("Every specimen was injected in both Batch 1 and Batch 2, so the batches are technical repeats of the same samples. The PCA is the check for whether a batch effect exists at all — if it does not, ComBat can be skipped."),
      p("Groups include NORM (normal), HCY, IVA, C3G, EB and LTRS, plus process blanks that carry an empty group label and are excluded from differential analysis automatically."),
      easyClose = TRUE, footer = modalButton("Close")))
  })

  # ---- per-mode pipelines ----------------------------------------------
  label_pos <- reactive({ if (isTruthy(input$label_pos)) input$label_pos else "Mode 1" })
  label_neg <- reactive({ if (isTruthy(input$label_neg)) input$label_neg else "Mode 2" })

  mode_pos <- pipelineServer("pos", label_pos, preset_pos)
  mode_neg <- pipelineServer("neg", label_neg, preset_neg)

  output$configure_heading_pos <- renderText({ paste0("Mode 1 — ", label_pos()) })
  output$configure_heading_neg <- renderText({ paste0("Mode 2 — ", label_neg()) })

  # Every mode that has been preprocessed, in display order. The second mode
  # only counts when the user has turned it on.
  modes <- reactive({
    out <- list()
    if (!is.null(mode_pos())) out[[length(out) + 1]] <- mode_pos()
    if (isTRUE(input$use_second_mode) && !is.null(mode_neg())) {
      out[[length(out) + 1]] <- mode_neg()
    }
    out
  })

  # ---- differential analysis -------------------------------------------
  de_runs <- differentialServer("de", modes)

  # ---- downloads --------------------------------------------------------
  ready_modes <- reactive({
    Filter(function(m) isTRUE(m$correction_complete), modes())
  })

  output$downloads_ready <- reactive({ length(ready_modes()) > 0 })
  outputOptions(output, "downloads_ready", suspendWhenHidden = FALSE)

  download_basename <- reactive({
    "metabo_tools"
  })

  output$download_excel <- downloadHandler(
    filename = function() {
      paste0(download_basename(), "_metabo_tools_results_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      ms <- ready_modes()
      req(length(ms) > 0)
      build_excel_package(ms, de_runs(), file)
    }
  )

  output$download_html <- downloadHandler(
    filename = function() {
      paste0(download_basename(), "_metabo_tools_report_", Sys.Date(), ".html")
    },
    content = function(file) {
      ms <- ready_modes()
      req(length(ms) > 0)
      generate_html_report(ms, de_runs(), file)
    }
  )
}

shinyApp(ui = ui, server = server)
