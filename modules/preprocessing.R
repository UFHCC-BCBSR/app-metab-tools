# ===========================
# Preprocessing and PCA
# ===========================
# Pure functions shared by every ionization mode: filtering, imputation,
# normalization, scaling, sample matching and PCA. No Shiny reactivity, so each
# can be run and checked from a plain R session.

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

perform_pca <- function(count_matrix, sample_data, has_bio_var = TRUE,
                       bio_var_name = "Biological Variable", has_batch = TRUE) {
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
  # Samples with no group or batch label (process blanks, for example) would
  # otherwise show up in the legend as "trace 0".
  label_blanks <- function(x) {
    x <- as.character(x)
    ifelse(is.na(x) | !nzchar(trimws(x)), "(unlabeled)", x)
  }
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample_name = colnames(count_matrix),
    batch = as.factor(label_blanks(sample_data$batch))
  )
  if (has_bio_var && !"no_bio_var" %in% sample_data$biological_var) {
    pca_df$biological_var <- as.factor(label_blanks(sample_data$biological_var))
  } else {
    pca_df$biological_var <- as.factor(rep("No biological variable", nrow(pca_df)))
  }
  variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
  list(pca_df = pca_df, variance_explained = variance_explained,
       has_bio_var = has_bio_var, bio_var_name = bio_var_name,
       has_batch = has_batch && length(unique(pca_df$batch)) > 1)
}

create_pca_plots <- function(pca_result, title) {
  if (is.null(pca_result)) return(NULL)
  pca_df       <- pca_result$pca_df
  var_exp      <- pca_result$variance_explained
  bio_var_name <- if (!is.null(pca_result$bio_var_name)) pca_result$bio_var_name else "Biological Variable"

  show_bio <- isTRUE(pca_result$has_bio_var) &&
    !all(pca_df$biological_var == "No biological variable")
  show_batch <- isTRUE(pca_result$has_batch)

  axis_titles <- list(
    x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2], 1), "%)")
  )

  bio_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                  "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")

  p_bio <- NULL
  if (show_bio) {
    n_bio <- length(unique(pca_df$biological_var))
    bio_pal <- if (n_bio <= length(bio_colors)) {
      bio_colors[seq_len(n_bio)]
    } else {
      colorRampPalette(bio_colors)(n_bio)
    }
    p_bio <- plot_ly(pca_df, x = ~PC1, y = ~PC2,
                     color = ~biological_var,
                     colors = bio_pal,
                     text = ~sample_name,
                     hovertemplate = "%{text}<extra></extra>",
                     type = "scatter", mode = "markers",
                     legendgrouptitle = list(text = paste0("<b>", bio_var_name, "</b>")),
                     legendgroup = "bio",
                     showlegend = TRUE) %>%
      layout(xaxis = list(title = axis_titles$x),
             yaxis = list(title = axis_titles$y))
  }

  p_batch <- NULL
  if (show_batch) {
    batch_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(pca_df$batch)))
    p_batch <- plot_ly(pca_df, x = ~PC1, y = ~PC2,
                       color = ~batch,
                       colors = batch_colors,
                       text = ~sample_name,
                       hovertemplate = "%{text}<extra></extra>",
                       type = "scatter", mode = "markers",
                       legendgrouptitle = list(text = "<b>Batch</b>"),
                       legendgroup = "batch",
                       showlegend = TRUE) %>%
      layout(xaxis = list(title = axis_titles$x),
             yaxis = list(title = axis_titles$y))
  }

  if (show_bio && show_batch) {
    subplot(p_bio, p_batch, nrows = 1, margin = 0.08,
            shareX = FALSE, shareY = FALSE,
            titleX = TRUE, titleY = TRUE) %>%
      layout(
        title = list(text = title, x = 0.5, y = 0.98, yanchor = "top"),
        legend = list(tracegroupgap = 40),
        margin = list(t = 90),
        annotations = list(
          list(x = 0.2, y = 1.02, text = paste0("<b>", bio_var_name, "</b>"),
               xref = "paper", yref = "paper", xanchor = "center",
               yanchor = "bottom", showarrow = FALSE, font = list(size = 13)),
          list(x = 0.8, y = 1.02, text = "<b>Batch</b>",
               xref = "paper", yref = "paper", xanchor = "center",
               yanchor = "bottom", showarrow = FALSE, font = list(size = 13))
        )
      )
  } else if (show_bio) {
    p_bio %>% layout(title = list(text = title, x = 0.5))
  } else if (show_batch) {
    p_batch %>% layout(title = list(text = title, x = 0.5))
  } else {
    plot_ly(pca_df, x = ~PC1, y = ~PC2, text = ~sample_name,
            hovertemplate = "%{text}<extra></extra>",
            type = "scatter", mode = "markers",
            marker = list(color = "#2980b9")) %>%
      layout(title = list(text = title, x = 0.5),
             xaxis = list(title = axis_titles$x),
             yaxis = list(title = axis_titles$y))
  }
}


# Match metadata sample names to feature-table columns. Names are compared with
# '-' and '.' treated alike because spreadsheet round-trips rewrite them, and a
# substring match is the last resort for run IDs carrying extra suffixes.
match_samples <- function(count_data, sample_data, feature_col, sample_name_col,
                          batch_col = NULL, bio_col = NULL, drop_cols = NULL) {
  count_processed <- count_data
  if (!is.null(drop_cols) && length(drop_cols) > 0) {
    count_processed <- count_processed[, !colnames(count_processed) %in% drop_cols, drop = FALSE]
  }
  potential_sample_cols <- setdiff(colnames(count_processed), feature_col)

  has_bio_var <- !is.null(bio_col) && nzchar(bio_col)
  has_batch <- !is.null(batch_col) && nzchar(batch_col)

  cols <- c(sample_name_col, if (has_batch) batch_col, if (has_bio_var) bio_col)
  sample_processed <- sample_data[, cols, drop = FALSE]
  colnames(sample_processed) <- c("sample_name",
                                  if (has_batch) "batch",
                                  if (has_bio_var) "biological_var")
  if (!has_batch) sample_processed$batch <- "no_batch"
  if (!has_bio_var) sample_processed$biological_var <- "no_bio_var"

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
  if (length(matched_samples) == 0) return(NULL)

  count_matrix <- as.matrix(count_processed[, matched_samples, drop = FALSE])
  # Feature names must be unique. limma's topTable() throws away every rowname
  # and substitutes row numbers as soon as it sees one duplicate, which would
  # silently rename all results, so duplicates are suffixed here instead.
  feature_names <- as.character(count_processed[[feature_col]])
  n_duplicate_features <- sum(duplicated(feature_names))
  rownames(count_matrix) <- make.unique(feature_names)

  sample_data_matched <- sample_processed[sample_processed$sample_name %in% sample_mapping, ]
  sample_data_matched <- sample_data_matched[match(sample_mapping, sample_data_matched$sample_name), ]

  full_sample_data <- sample_data[sample_data[[sample_name_col]] %in% sample_mapping, ]
  full_sample_data <- full_sample_data[match(sample_mapping, full_sample_data[[sample_name_col]]), ]

  list(
    count_matrix = count_matrix,
    sample_data = sample_data_matched,
    full_sample_data = full_sample_data,
    matched_samples = matched_samples,
    unmatched_samples = unmatched_samples,
    has_bio_var = has_bio_var,
    has_batch = has_batch,
    n_duplicate_features = n_duplicate_features
  )
}

# Does this mode have a batch variable with more than one level? Everything
# downstream — ComBat, the batch covariate option, the PCA coloring — depends
# on this one question, so it is asked in exactly one place.
has_usable_batch <- function(matched) {
  if (is.null(matched) || !isTRUE(matched$has_batch)) return(FALSE)
  length(unique(matched$sample_data$batch)) > 1
}
