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
  threshold <- quantile(iqrs, top_percent / 100)
  count_matrix[iqrs >= threshold, ]
}

sample_normalize <- function(mat, method) {
  if (method == "SumNorm") {
    col_sums <- colSums(mat, na.rm = TRUE)
    mat <- sweep(mat, 2, col_sums, "/") * median(col_sums)
  } else if (method == "MedianNorm") {
    col_medians <- apply(mat, 2, median, na.rm = TRUE)
    mat <- sweep(mat, 2, col_medians, "/") * median(col_medians)
  }
  mat
}

scale_matrix <- function(mat, method) {
  if (method == "ParetoNorm") {
    row_means <- rowMeans(mat, na.rm = TRUE)
    row_sds <- apply(mat, 1, sd, na.rm = TRUE)
    mat <- sweep(mat, 1, row_means, "-")
    mat <- sweep(mat, 1, sqrt(row_sds), "/")
  } else if (method == "AutoNorm") {
    row_means <- rowMeans(mat, na.rm = TRUE)
    row_sds <- apply(mat, 1, sd, na.rm = TRUE)
    mat <- sweep(mat, 1, row_means, "-")
    mat <- sweep(mat, 1, row_sds, "/")
  }
  mat
}

perform_pca <- function(count_matrix, sample_data, has_bio_var = TRUE) {
  plot_matrix <- log2(abs(count_matrix) + 1)
  plot_matrix <- plot_matrix[apply(plot_matrix, 1, function(x) all(is.finite(x))), ]
  if (nrow(plot_matrix) > 1000) {
    vars <- apply(plot_matrix, 1, var, na.rm = TRUE)
    top_features <- names(sort(vars, decreasing = TRUE)[1:1000])
    plot_matrix <- plot_matrix[top_features, ]
  }
  vars <- apply(plot_matrix, 1, var, na.rm = TRUE)
  plot_matrix <- plot_matrix[vars > 0 & !is.na(vars), ]
  if (nrow(plot_matrix) == 0) stop("No valid features remaining for PCA.")
  pca_data <- t(plot_matrix)
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

create_pca_plots <- function(pca_result, title) {
  req(pca_result)
  pca_df <- pca_result$pca_df
  var_exp <- pca_result$variance_explained
  has_bio_var <- pca_result$has_bio_var
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

# ===========================
# Generate README content
# ===========================
generate_readme <- function(params) {
  data.frame(
    Section = c(
      "ABOUT THIS FILE",
      "ABOUT THIS FILE",
      "PIPELINE OVERVIEW",
      "PIPELINE OVERVIEW",
      "PIPELINE OVERVIEW",
      "PIPELINE OVERVIEW",
      "PIPELINE OVERVIEW",
      "PIPELINE OVERVIEW",
      "SHEET: data_original",
      "SHEET: data_preprocessed",
      "SHEET: data_normalized",
      "SHEET: data_normalized_scaled",
      "SHEET: data_batch_corrected_scaled",
      "SHEET: metadata",
      "PROCESSING LOG",
      "PROCESSING LOG",
      "PROCESSING LOG",
      "PROCESSING LOG",
      "PROCESSING LOG",
      "PROCESSING LOG",
      "PROCESSING LOG",
      "PROCESSING LOG"
    ),
    Detail = c(
      paste("Generated by Metabo Tools on", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "This file contains all processed versions of your data. See each sheet description below to understand which to use for your analysis.",
      "Step 1: Sample matching — count matrix columns matched to metadata sample names",
      if (params$do_missing_filter) "Step 2: Missing value filter — features with >50% missing values removed" else "Step 2: Missing value filter — skipped",
      if (params$do_imputation) "Step 3: KNN imputation — remaining missing/zero values imputed" else "Step 3: KNN imputation — skipped",
      if (params$do_iqr_filter) paste0("Step 4: IQR filter — bottom ", params$iqr_threshold, "% lowest variance features removed") else "Step 4: IQR filter — skipped",
      paste("Step 5: Sample normalization —", params$row_norm),
      "Step 6: Log2 transformation",
      paste("Step 7: Batch correction — ComBat (biological variable:", params$bio_var, ")"),
      paste("Step 8: Scaling —", params$scale_norm),
      "Use this if you want to start the entire processing pipeline over from scratch. This is the raw matched matrix with no processing applied.",
      "Use this if you want to apply your own normalization, transformation, and scaling. Missing values have been filtered and imputed but nothing else has been applied.",
      "Use this if you want to apply your own scaling and/or batch correction method. Sample normalization and log2 transformation have been applied.",
      paste0("Use this for statistical modeling where you want to account for batch effects by including the batch column from the metadata file as a covariate in your model (e.g. limma). ",
             "Batch effects are still present in this data — they are controlled for in the model, not removed from the matrix. ",
             "Scaling applied: ", params$scale_norm, "."),
      paste0("Use this for PCA, clustering, heatmaps, and statistical analysis where batch effects have already been removed from the data matrix. ",
             "Do not include batch as a covariate in your model when using this data — it has already been corrected. ",
             "Scaling applied: ", params$scale_norm, "."),
      "Use the batch column when building your model design if using data_normalized_scaled. Includes all matched samples and their group assignments.",
      paste("Features in raw matched matrix:", params$n_original),
      if (params$do_missing_filter) paste("Features after missing value filter:", params$n_after_missing) else "Missing value filter: not applied",
      if (params$do_imputation) paste("Missing/zero values imputed:", params$n_imputed) else "KNN imputation: not applied",
      if (params$do_iqr_filter) paste("Features after IQR filter:", params$n_after_iqr) else "IQR filter: not applied",
      paste("Features in final matrix:", params$n_final),
      paste("Samples matched:", params$n_samples)
    )
  )
}

# ===========================
# Generate HTML Report
# ===========================
embed_plotly <- function(fig) {
  if (is.null(fig)) return("")
  div_id <- paste0("plotly_", paste(sample(letters, 10, replace = TRUE), collapse = ""))
  fig_json <- plotly::plotly_json(fig, jsonedit = FALSE)
  paste0(
    '<div id="', div_id, '" style="width:100%;height:400px;"></div>',
    '<script>',
    'Plotly.newPlot("', div_id, '",',
    fig_json,
    ');</script>'
  )
}

generate_html_report <- function(params, pca_before_result, pca_after_result, output_path) {

  pca_before_batch_fig <- NULL
  pca_before_bio_fig <- NULL
  pca_after_batch_fig <- NULL
  pca_after_bio_fig <- NULL

  if (!is.null(pca_before_result)) {
    pca_df <- pca_before_result$pca_df
    var_exp <- pca_before_result$variance_explained
    batch_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(pca_df$batch)))
    pca_before_batch_fig <- plot_ly(pca_df, x = ~PC1, y = ~PC2,
                                    color = ~batch, colors = batch_colors,
                                    text = ~sample_name,
                                    hovertemplate = "%{text}<extra></extra>",
                                    type = "scatter", mode = "markers") %>%
      layout(title = "Before Batch Correction — Colored by Batch",
             xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
             yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")))

    if (pca_before_result$has_bio_var && !all(pca_df$biological_var == "No biological variable")) {
      bio_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                      "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
      pca_before_bio_fig <- plot_ly(pca_df, x = ~PC1, y = ~PC2,
                                    color = ~biological_var,
                                    colors = bio_colors[1:length(unique(pca_df$biological_var))],
                                    text = ~sample_name,
                                    hovertemplate = "%{text}<extra></extra>",
                                    type = "scatter", mode = "markers") %>%
        layout(title = "Before Batch Correction — Colored by Biological Variable",
               xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
               yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")))
    }
  }

  if (!is.null(pca_after_result)) {
    pca_df <- pca_after_result$pca_df
    var_exp <- pca_after_result$variance_explained
    batch_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(pca_df$batch)))
    pca_after_batch_fig <- plot_ly(pca_df, x = ~PC1, y = ~PC2,
                                   color = ~batch, colors = batch_colors,
                                   text = ~sample_name,
                                   hovertemplate = "%{text}<extra></extra>",
                                   type = "scatter", mode = "markers") %>%
      layout(title = "After Batch Correction — Colored by Batch",
             xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
             yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")))

    if (pca_after_result$has_bio_var && !all(pca_df$biological_var == "No biological variable")) {
      bio_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                      "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
      pca_after_bio_fig <- plot_ly(pca_df, x = ~PC1, y = ~PC2,
                                   color = ~biological_var,
                                   colors = bio_colors[1:length(unique(pca_df$biological_var))],
                                   text = ~sample_name,
                                   hovertemplate = "%{text}<extra></extra>",
                                   type = "scatter", mode = "markers") %>%
        layout(title = "After Batch Correction — Colored by Biological Variable",
               xaxis = list(title = paste0("PC1 (", round(var_exp[1], 1), "%)")),
               yaxis = list(title = paste0("PC2 (", round(var_exp[2], 1), "%)")))
    }
  }

  scale_label <- switch(params$scale_norm,
                        "ParetoNorm" = "Pareto scaling",
                        "AutoNorm" = "Auto scaling (mean-centering and division by standard deviation)",
                        "None" = "no scaling")

  norm_label <- switch(params$row_norm,
                       "SumNorm" = "sum normalization",
                       "MedianNorm" = "median normalization",
                       "None" = "no sample normalization")

sva_version <- as.character(packageVersion("sva"))
impute_version <- as.character(packageVersion("impute"))

scaling_sentence <- switch(params$scale_norm,
  "ParetoNorm" = "Data were subsequently Pareto scaled (mean-centered and divided by the square root of the standard deviation of each feature). ",
  "AutoNorm"   = "Data were subsequently auto scaled (mean-centered and divided by the standard deviation of each feature). ",
  "None"       = ""
)

methods_paragraph <- paste0(
  "Metabolomics data processing was performed using Metabo Tools, a web application developed by the ",
  "University of Florida Southeast Center for Integrated Metabolomics (SECIM; ",
  "<a href='https://secim.ufl.edu/'>https://secim.ufl.edu/</a>). ",
  "All analyses were conducted in R using the sva package (v", sva_version, "; Johnson et al., 2007) ",
  "for batch correction and the impute package (v", impute_version, ") for missing value imputation. ",
  "Raw peak intensity data from ", params$n_samples, " samples and ", params$n_original, " features were imported. ",
  if (params$do_missing_filter) paste0(
    "Features with greater than 50% missing or zero values were removed, retaining ",
    params$n_after_missing, " features. ") else "",
  if (params$do_imputation) paste0(
    "Remaining missing and zero values (n = ", params$n_imputed, ") were imputed using ",
    "k-nearest neighbours (KNN, k = 10) as implemented in the Bioconductor impute package. ") else "",
  if (params$do_iqr_filter) paste0(
    "Low-variance features were removed by filtering the bottom ", params$iqr_threshold,
    "% by interquartile range, retaining ", params$n_after_iqr, " features. ") else "",
  "Samples were normalized by ", norm_label, " followed by log2 transformation. ",
  "Batch effects were removed using the ComBat empirical Bayes method",
  if (params$bio_var != "None") paste0(
    " with biological group (", params$bio_var, ") included as a covariate to preserve ",
    "biological variation during correction (Johnson et al., 2007). ") else " (Johnson et al., 2007). ",
  scaling_sentence,
  "The final processed matrix contained ", params$n_final, " features across ", params$n_samples, " samples.",
  "<br><br>",
  "<strong>Reference:</strong> Johnson WE, Li C, Rabinovic A. ",
  "Adjusting batch effects in microarray expression data using empirical Bayes methods. ",
  "<em>Biostatistics</em>. 2007;8(1):118-127. doi:10.1093/biostatistics/kxj037"
)

  html_content <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Metabo Tools Processing Report</title>
  <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <style>
    body { font-family: Arial, sans-serif; max-width: 1100px; margin: 40px auto; padding: 0 20px; color: #333; }
    h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
    h2 { color: #2980b9; margin-top: 40px; }
    h3 { color: #16a085; }
    .meta { background: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 30px; font-size: 14px; }
    .methods { background: #eaf4fb; padding: 20px; border-radius: 5px; border-left: 5px solid #3498db; margin-bottom: 30px; line-height: 1.7; }
    .step-log { background: #f8f9fa; padding: 15px; border-radius: 5px; font-family: monospace; font-size: 13px; line-height: 1.8; }
    .step-log .removed { color: #e74c3c; }
    .step-log .kept { color: #27ae60; }
    .file-table { width: 100%; border-collapse: collapse; margin-top: 15px; }
    .file-table th { background: #2980b9; color: white; padding: 10px 15px; text-align: left; }
    .file-table td { padding: 10px 15px; border-bottom: 1px solid #ddd; vertical-align: top; }
    .file-table tr:nth-child(even) { background: #f8f9fa; }
    .file-table .filename { font-family: monospace; font-weight: bold; color: #2c3e50; white-space: nowrap; }
    .highlight-a { background: #d5f5e3 !important; }
    .highlight-b { background: #d6eaf8 !important; }
    .tag { display: inline-block; padding: 2px 8px; border-radius: 3px; font-size: 11px; font-weight: bold; margin-left: 5px; }
    .tag-primary { background: #27ae60; color: white; }
    .tag-secondary { background: #2980b9; color: white; }
    .pca-section { margin-top: 20px; }
    .pca-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 15px; }
    .warning { background: #fef9e7; border-left: 5px solid #f39c12; padding: 15px; border-radius: 5px; margin-top: 20px; }
  </style>
</head>
<body>

<h1>Metabo Tools — Processing Report</h1>

<div class="meta">
  <strong>Generated:</strong> ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '<br>
  <strong>Samples:</strong> ', params$n_samples, '<br>
  <strong>Features (input):</strong> ', params$n_original, '<br>
  <strong>Features (final):</strong> ', params$n_final, '<br>
  <strong>Batch correction:</strong> ComBat<br>
  <strong>Scaling:</strong> ', scale_label, '
</div>

<h2>Methods Paragraph</h2>
<div class="methods">
  <p>', methods_paragraph, '</p>
  <p style="font-size:12px; color:#888; margin-top:10px;">
    This paragraph is provided as a draft for your methods section. Please review and adjust as needed for your specific analysis context.
  </p>
</div>

<h2>Processing Steps</h2>
<div class="step-log">
  <strong>Input:</strong> ', params$n_original, ' features, ', params$n_samples, ' samples<br>',
  if (params$do_missing_filter) paste0('  &#8594; Missing value filter (&gt;50%): <span class="removed">', params$n_original - params$n_after_missing, ' features removed</span>, <span class="kept">', params$n_after_missing, ' retained</span><br>') else '  &#8594; Missing value filter: skipped<br>',
  if (params$do_imputation) paste0('  &#8594; KNN imputation: <span class="kept">', params$n_imputed, ' values imputed</span><br>') else '  &#8594; KNN imputation: skipped<br>',
  if (params$do_iqr_filter) paste0('  &#8594; IQR filter (bottom ', params$iqr_threshold, '%): <span class="removed">', params$n_after_missing - params$n_after_iqr, ' features removed</span>, <span class="kept">', params$n_after_iqr, ' retained</span><br>') else '  &#8594; IQR filter: skipped<br>',
  '  &#8594; Sample normalization: ', norm_label, '<br>
  &#8594; Log2 transformation applied<br>
  &#8594; ComBat batch correction applied',
  if (params$bio_var != "None") paste0(' (biological variable: ', params$bio_var, ')') else '',
  '<br>
  &#8594; Scaling: ', scale_label, '<br>
  <strong>Output:</strong> ', params$n_final, ' features, ', params$n_samples, ' samples
</div>

<h2>Output Files Guide</h2>
<p>Your download contains an Excel file with the following sheets. Use the table below to identify which sheet to use for your analysis.</p>

<table class="file-table">
  <tr>
    <th>Sheet</th>
    <th>Contents</th>
    <th>Use this for</th>
  </tr>
  <tr>
    <td class="filename">data_original</td>
    <td>Raw matched matrix, no processing</td>
    <td>Archive. Use this if you want to start the entire processing pipeline over from scratch.</td>
  </tr>
  <tr>
    <td class="filename">data_preprocessed</td>
    <td>After missing value filtering and KNN imputation</td>
    <td>Use this if you want to apply your own normalization, transformation, and scaling.</td>
  </tr>
  <tr>
    <td class="filename">data_normalized</td>
    <td>After sample normalization and log2 transformation</td>
    <td>Use this if you want to apply your own scaling and/or batch correction method.</td>
  </tr>
  <tr class="highlight-b">
    <td class="filename">data_normalized_scaled <span class="tag tag-secondary">limma + batch covariate</span></td>
    <td>After sample normalization, log2 transformation, and scaling (', scale_label, ')</td>
    <td>Use this for statistical modeling where you want to account for batch effects by including the batch column from the metadata file as a covariate in your model (e.g. limma). Batch effects are still present in this data — they are controlled for in the model, not removed from the matrix.</td>
  </tr>
  <tr class="highlight-a">
    <td class="filename">data_batch_corrected_scaled <span class="tag tag-primary">PCA / clustering / limma</span></td>
    <td>After sample normalization, log2 transformation, ComBat batch correction, and scaling (', scale_label, ')</td>
    <td>Use this for PCA, clustering, heatmaps, and statistical analysis where batch effects have already been removed from the data matrix. Do not include batch as a covariate in your model when using this data.</td>
  </tr>
  <tr>
    <td class="filename">metadata</td>
    <td>Sample metadata with batch column</td>
    <td>Use the batch column when building your model design if using data_normalized_scaled. Includes all matched samples and their group assignments.</td>
  </tr>
</table>

<div class="warning">
  <strong>Note on batch correction and statistical modeling:</strong> These two approaches are mutually exclusive.
  If you use <code>data_batch_corrected_scaled</code>, do not include batch as a covariate in your model.
  If you use <code>data_normalized_scaled</code>, include batch as a covariate. Do not do both.
</div>

<h2>PCA Plots</h2>
<p>PCA plots are shown before and after batch correction. Points are colored by batch and by biological variable to help assess whether batch effects have been successfully removed while biological signal is preserved.</p>
')

  html_content <- paste0(html_content, '<div class="pca-section"><h3>Before Batch Correction</h3><div class="pca-grid">')

  if (!is.null(pca_before_batch_fig)) {
    html_content <- paste0(html_content, '<div>', embed_plotly(pca_before_batch_fig), '</div>')
  }
  if (!is.null(pca_before_bio_fig)) {
    html_content <- paste0(html_content, '<div>', embed_plotly(pca_before_bio_fig), '</div>')
  }

  html_content <- paste0(html_content, '</div><h3>After Batch Correction</h3><div class="pca-grid">')

  if (!is.null(pca_after_batch_fig)) {
    html_content <- paste0(html_content, '<div>', embed_plotly(pca_after_batch_fig), '</div>')
  }
  if (!is.null(pca_after_bio_fig)) {
    html_content <- paste0(html_content, '<div>', embed_plotly(pca_after_bio_fig), '</div>')
  }

  html_content <- paste0(html_content, '
</div></div>

</body>
</html>
')

  writeLines(html_content, output_path)
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
      menuItem("Batch Correction & Scaling", tabName = "correction", icon = icon("magic")),
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
                    p("Want to try the app with example data first?"),
                    actionButton("load_demo", "Load Test Data",
                                 class = "btn-success btn-lg",
                                 icon = icon("play")),
                    br(), br(),
                    div(id = "demo_info", style = "display: none;",
                        div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; border: 1px solid #c3e6cb;",
                            p(strong("Demo data loaded!"),
                              "Go to the Preprocess Data tab to proceed."))
                    )
                )
              ),
              fluidRow(
                box(title = "Upload Count Data", status = "primary", solidHeader = TRUE, width = 6,
                    fileInput("count_file", "Choose Count Data File (.csv or .xlsx)",
                              accept = c(".csv", ".xlsx")),
                    conditionalPanel(
                      condition = "output.count_uploaded",
                      h4("Preview:"),
                      DT::dataTableOutput("count_preview")
                    )
                ),
                box(title = "Upload Sample Metadata", status = "primary", solidHeader = TRUE, width = 6,
                    fileInput("sample_file", "Choose Sample Metadata File (.csv or .xlsx)",
                              accept = c(".csv", ".xlsx")),
                    conditionalPanel(
                      condition = "output.sample_uploaded",
                      h4("Preview:"),
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
                      checkboxInput("has_bio_var", "My data has a biological variable of interest", value = TRUE),
                      conditionalPanel(
                        condition = "input.has_bio_var == true",
                        selectInput("bio_col", "Select Biological Variable Column:", choices = NULL)
                      ),
                      br(),
                      verbatimTextOutput("sample_config_summary")
                    )
                )
              ),
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
              fluidRow(
                box(title = "Step 2: Preprocessing", status = "warning",
                    solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.matching_previewed",
                      fluidRow(
                        column(6,
                               div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
                                   p(strong("Filtering and imputation:")),
                                   checkboxInput("do_missing_filter",
                                                 "Filter features with >50% missing values",
                                                 value = TRUE),
                                   checkboxInput("do_imputation",
                                                 "Impute remaining missing values (KNN)",
                                                 value = TRUE),
                                   checkboxInput("do_iqr_filter",
                                                 "Filter low-variance features (IQR)",
                                                 value = FALSE),
                                   conditionalPanel(
                                     condition = "input.do_iqr_filter == true",
                                     numericInput("iqr_threshold",
                                                  "Remove bottom X% lowest variance features:",
                                                  value = 10, min = 1, max = 50, step = 5)
                                   )
                               )
                        ),
                        column(6,
                               div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px;",
                                   p(strong("Sample normalization:")),
                                   selectInput("row_norm", "Method:",
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

      # Batch Correction & Scaling tab
      tabItem(tabName = "correction",
              fluidRow(
                box(title = "Batch Correction + Scaling", status = "success",
                    solidHeader = TRUE, width = 12,
                    conditionalPanel(
                      condition = "output.pca_before_complete",
                      div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
                          p("ComBat will be applied to your log2-normalized data. Select a scaling method to apply after batch correction.",
                            "The same scaling will also be applied to the non-batch-corrected data so you have both versions available in your download."),
                          br(),
                          selectInput("scale_norm", "Scaling method:",
                                      choices = list(
                                        "Pareto Scaling" = "ParetoNorm",
                                        "Auto Scaling" = "AutoNorm",
                                        "None" = "None"
                                      ), selected = "ParetoNorm")
                      ),
                      actionButton("run_combat", "Run Batch Correction + Scale",
                                   class = "btn-success btn-lg"),
                      br(), br()
                    ),
                    conditionalPanel(
                      condition = "output.correction_complete",
                      h4("PCA after batch correction and scaling:"),
                      plotlyOutput("pca_after", height = "400px")
                    )
                )
              )
      ),

      # About tab
      tabItem(tabName = "about",
              fluidRow(
                box(title = "About Metabo Tools", status = "primary", solidHeader = TRUE, width = 12,
                    p(strong("Metabo Tools"), "is a web application for metabolomics data preprocessing and batch correction."),
                    p(tags$em("This tool is currently under active development and testing. Please use results with appropriate caution and report any issues.")),
                    br(),
                    h4("Source Code"),
                    p("The source code is available on GitHub:"),
                    p(tags$a(href = "https://github.com/UFHCC-BCBSR/app-metabo-tools",
                             "https://github.com/UFHCC-BCBSR/app-metabo-tools",
                             target = "_blank")),
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
                      condition = "output.correction_complete",
                      p("Download your complete results package. The Excel file contains all processed data versions and a README sheet. The HTML report contains a methods paragraph, processing log, file guide, and PCA plots for your records."),
                      br(),
                      downloadButton("download_excel", "Download Excel Data Package",
                                     class = "btn-primary btn-lg"),
                      br(), br(),
                      downloadButton("download_html", "Download HTML Report",
                                     class = "btn-info btn-lg")
                    ),
                    conditionalPanel(
                      condition = "!output.correction_complete",
                      p("Complete batch correction to enable downloads.")
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
    normalized_data = NULL,
    normalized_scaled_data = NULL,
    batch_corrected_data = NULL,
    batch_corrected_scaled_data = NULL,
    pca_before = NULL,
    pca_after = NULL,
    matching_previewed = FALSE,
    preprocessing_complete = FALSE,
    pca_before_complete = FALSE,
    preprocessing_log = NULL,
    processing_params = NULL
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
        full_sample_data <- values$sample_data[
          values$sample_data[[input$sample_name_col]] %in% sample_mapping, ]
        full_sample_data <- full_sample_data[
          match(sample_mapping, full_sample_data[[input$sample_name_col]]), ]

        values$matched_data <- list(
          count_matrix = count_matrix,
          sample_data = sample_data_matched,
          full_sample_data = full_sample_data,
          matched_samples = matched_samples,
          unmatched_samples = unmatched_samples,
          has_bio_var = has_bio_var
        )
        values$preprocessed_data <- NULL
        values$normalized_data <- NULL
        values$normalized_scaled_data <- NULL
        values$batch_corrected_data <- NULL
        values$batch_corrected_scaled_data <- NULL
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
      n_after_iqr <- n_original

      mat <- apply(mat, 2, as.numeric)
      rownames(mat) <- rownames(values$matched_data$count_matrix)

      if (input$do_missing_filter) {
        before <- nrow(mat)
        mat <- filter_missing(mat, threshold = 0.5)
        n_after_missing <- nrow(mat)
        log_lines <- c(log_lines,
                       paste0("Missing value filter: ", before, " -> ", nrow(mat),
                              " features (removed ", before - nrow(mat), ")"))
      }

      if (input$do_imputation) {
        n_imputed <- sum(is.na(mat) | mat == 0)
        mat[mat == 0] <- NA
        mat <- impute_knn(mat)
        log_lines <- c(log_lines,
                       paste0("KNN imputation: ", n_imputed, " values imputed"))
      }

      if (input$do_iqr_filter) {
        before <- nrow(mat)
        mat <- filter_iqr(mat, top_percent = input$iqr_threshold)
        n_after_iqr <- nrow(mat)
        log_lines <- c(log_lines,
                       paste0("IQR filter: ", before, " -> ", nrow(mat),
                              " features (removed ", before - nrow(mat), ")"))
      } else {
        n_after_iqr <- nrow(mat)
      }

      # Store preprocessed (filter + impute only)
      values$preprocessed_data <- mat

      # Sample normalization
      mat_norm <- sample_normalize(mat, input$row_norm)

      # Log2 transformation
      if (any(mat_norm <= 0, na.rm = TRUE)) {
        mat_norm[mat_norm <= 0] <- 0.01
      }
      mat_norm <- log2(mat_norm)

      # Store normalized (no scaling)
      values$normalized_data <- mat_norm

      log_lines <- c(
        paste0("Input features: ", n_original),
        log_lines,
        paste0("Sample normalization: ", input$row_norm),
        "Log2 transformation applied",
        paste0("Final features: ", nrow(mat_norm))
      )

      values$processing_params <- list(
        do_missing_filter = input$do_missing_filter,
        do_imputation = input$do_imputation,
        do_iqr_filter = input$do_iqr_filter,
        iqr_threshold = if (input$do_iqr_filter) input$iqr_threshold else NA,
        row_norm = input$row_norm,
        n_original = n_original,
        n_after_missing = n_after_missing,
        n_after_iqr = n_after_iqr,
        n_imputed = n_imputed,
        n_final = nrow(mat_norm),
        n_samples = ncol(mat_norm),
        bio_var = if (values$matched_data$has_bio_var) input$bio_col else "None"
      )

      values$preprocessing_log <- paste(log_lines, collapse = "\n")
      values$preprocessing_complete <- TRUE
      values$pca_before_complete <- FALSE
      showNotification("Preprocessing complete!", type = "message")
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
        values$matched_data$has_bio_var
      )
      values$pca_before_complete <- TRUE
      showNotification("PCA complete!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error in PCA:", e$message), type = "error")
    })
  })

  observeEvent(input$run_combat, {
    req(values$normalized_data, values$matched_data, input$scale_norm)
    shinyjs::disable("run_combat")
    showModal(modalDialog("Running ComBat batch correction...", footer = NULL))
    tryCatch({
      mat_norm <- values$normalized_data
      sample_data <- values$matched_data$sample_data
      has_bio_var <- values$matched_data$has_bio_var
      batch_factor <- as.factor(sample_data$batch)

      mod <- if (has_bio_var && !"no_bio_var" %in% sample_data$biological_var) {
        model.matrix(~as.factor(sample_data$biological_var))
      } else {
        NULL
      }

      mat_corrected <- ComBat(dat = mat_norm, batch = batch_factor,
                              mod = mod, par.prior = TRUE)

      # Scale both corrected and non-corrected
      mat_norm_scaled <- scale_matrix(mat_norm, input$scale_norm)
      mat_corrected_scaled <- scale_matrix(mat_corrected, input$scale_norm)

      values$normalized_scaled_data <- mat_norm_scaled
      values$batch_corrected_data <- mat_corrected
      values$batch_corrected_scaled_data <- mat_corrected_scaled

      # Update params with scale
      values$processing_params$scale_norm <- input$scale_norm

      values$pca_after <- perform_pca(mat_corrected_scaled, sample_data, has_bio_var)

      removeModal()
      shinyjs::enable("run_combat")
      showNotification("Batch correction and scaling complete!", type = "message")
    }, error = function(e) {
      removeModal()
      shinyjs::enable("run_combat")
      showNotification(paste("Error in batch correction:", e$message), type = "error")
    })
  }, ignoreInit = TRUE)

  output$pca_before <- renderPlotly({
    req(values$pca_before)
    create_pca_plots(values$pca_before, "Before batch correction")
  })

  output$pca_after <- renderPlotly({
    req(values$pca_after)
    create_pca_plots(values$pca_after, "After batch correction and scaling")
  })

  output$matching_previewed <- reactive({ values$matching_previewed })
  output$preprocessing_complete <- reactive({ values$preprocessing_complete })
  output$pca_before_complete <- reactive({ values$pca_before_complete })
  output$correction_complete <- reactive({ !is.null(values$batch_corrected_scaled_data) })

  outputOptions(output, "count_uploaded", suspendWhenHidden = FALSE)
  outputOptions(output, "sample_uploaded", suspendWhenHidden = FALSE)
  outputOptions(output, "ready_to_preview", suspendWhenHidden = FALSE)
  outputOptions(output, "matching_previewed", suspendWhenHidden = FALSE)
  outputOptions(output, "preprocessing_complete", suspendWhenHidden = FALSE)
  outputOptions(output, "pca_before_complete", suspendWhenHidden = FALSE)
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

  output$preprocessing_summary <- renderText({
    req(values$preprocessing_log)
    values$preprocessing_log
  })

  output$count_config_summary <- renderText({
    req(values$count_data, input$feature_col)
    dropped <- if (is.null(input$drop_cols)) "None" else paste(input$drop_cols, collapse = ", ")
    paste0("Feature column: ", input$feature_col, "\n",
           "Dropped columns: ", dropped)
  })

  output$sample_config_summary <- renderText({
    req(values$sample_data, input$batch_col, input$sample_name_col)
    bio_text <- if (!is.null(input$has_bio_var) && input$has_bio_var && !is.null(input$bio_col)) {
      paste("Biological variable column:", input$bio_col)
    } else {
      "No biological variable selected"
    }
    paste0("Sample name column: ", input$sample_name_col, "\n",
           "Batch column: ", input$batch_col, "\n",
           bio_text)
  })

  # ===========================
  # Downloads
  # ===========================
  output$download_excel <- downloadHandler(
    filename = function() paste0("metabo_tools_results_", Sys.Date(), ".xlsx"),
    content = function(file) {
      req(values$matched_data, values$preprocessed_data, values$normalized_data,
          values$normalized_scaled_data, values$batch_corrected_scaled_data,
          values$processing_params)

      params <- values$processing_params
      readme_df <- generate_readme(params)

      mat_to_df <- function(mat) {
        tibble::rownames_to_column(as.data.frame(mat), var = "feature")
      }

      wb <- createWorkbook()

      addWorksheet(wb, "README")
      addWorksheet(wb, "data_original")
      addWorksheet(wb, "data_preprocessed")
      addWorksheet(wb, "data_normalized")
      addWorksheet(wb, "data_normalized_scaled")
      addWorksheet(wb, "data_batch_corrected_scaled")
      addWorksheet(wb, "metadata")

      writeData(wb, "README", readme_df)
      writeData(wb, "data_original", mat_to_df(values$matched_data$count_matrix))
      writeData(wb, "data_preprocessed", mat_to_df(values$preprocessed_data))
      writeData(wb, "data_normalized", mat_to_df(values$normalized_data))
      writeData(wb, "data_normalized_scaled", mat_to_df(values$normalized_scaled_data))
      writeData(wb, "data_batch_corrected_scaled", mat_to_df(values$batch_corrected_scaled_data))
      writeData(wb, "metadata", values$matched_data$full_sample_data)

      # Style README
      headerStyle <- createStyle(textDecoration = "bold", fgFill = "#2980b9",
                                 fontColour = "white", wrapText = TRUE)
      addStyle(wb, "README", headerStyle, rows = 1, cols = 1:2, gridExpand = TRUE)
      setColWidths(wb, "README", cols = 1:2, widths = c(35, 90))

      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

  output$download_html <- downloadHandler(
    filename = function() paste0("metabo_tools_report_", Sys.Date(), ".html"),
    content = function(file) {
      req(values$processing_params)
      generate_html_report(
        params = values$processing_params,
        pca_before_result = values$pca_before,
        pca_after_result = values$pca_after,
        output_path = file
      )
    }
  )
}

shinyApp(ui = ui, server = server)
