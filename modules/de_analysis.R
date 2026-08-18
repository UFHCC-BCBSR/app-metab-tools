# ===========================
# Differential Analysis — core logic
# ===========================
# Pure functions: model fitting, plots, and report/Excel content for the
# differential analysis module. Nothing here touches Shiny reactivity, so every
# piece can be run and checked from a plain R session.
#
# Modelling follows limma: one linear model per feature on log2 intensities,
# with variance moderated across features by empirical Bayes. Two matrices can
# be modelled and the two are mutually exclusive, exactly as the processing
# report warns:
#   "normalized"      log2-normalized data, batch included as a model covariate
#   "batch_corrected" ComBat-corrected data, no batch term in the design
# Scaled matrices are never modelled — Pareto/auto scaling would leave the
# coefficients in scaled units and destroy the fold change interpretation.

DE_MAX_HEATMAP_FEATURES <- 50
DE_REPORT_BOXPLOT_FEATURES <- 6
DE_REPORT_TABLE_ROWS <- 25

DE_COLORS <- list(
  up = "#c0392b",
  down = "#2980b9",
  ns = "#bdc3c7"
)

de_group_colors <- function(levels) {
  base <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
            "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  if (length(levels) <= length(base)) {
    base[seq_along(levels)]
  } else {
    colorRampPalette(base)(length(levels))
  }
}

# ---------------------------
# Design and model fitting
# ---------------------------

# Metadata columns usable as a grouping factor: at least two levels, and few
# enough levels that the column is a grouping and not a per-sample identifier.
de_group_candidates <- function(sample_df) {
  if (is.null(sample_df) || ncol(sample_df) == 0) return(character(0))
  n <- nrow(sample_df)
  keep <- vapply(sample_df, function(col) {
    lv <- unique(col[!is.na(col)])
    length(lv) >= 2 && length(lv) <= max(2, floor(n / 2))
  }, logical(1))
  colnames(sample_df)[keep]
}

# Design matrix with one column per group (no intercept) plus optional batch
# dummies. Group columns are named from the factor levels so contrasts can be
# written in terms of the levels the user selected.
de_build_design <- function(groups, batch = NULL) {
  groups <- droplevels(as.factor(groups))
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- de_level_keys(groups)
  if (!is.null(batch)) {
    batch <- droplevels(as.factor(batch))
    if (nlevels(batch) > 1) {
      bmat <- model.matrix(~ batch)[, -1, drop = FALSE]
      colnames(bmat) <- make.names(paste0("batch_", levels(batch)[-1]), unique = TRUE)
      design <- cbind(design, bmat)
    }
  }
  design
}

# Map group levels to syntactically valid design column names.
de_level_keys <- function(groups) {
  groups <- droplevels(as.factor(groups))
  setNames(make.names(levels(groups), unique = TRUE), levels(groups))
}

# Fit the model and return a tidy results table. `mode` is "pairwise" for a
# moderated t-test between two groups, or "global" for a moderated F-test
# asking whether a feature differs across any of the groups.
de_run <- function(mat, groups, batch = NULL,
                   mode = c("pairwise", "global"),
                   group_a = NULL, group_b = NULL,
                   matrix_source = "normalized",
                   use_batch = FALSE,
                   group_var = "group",
                   fdr_cutoff = 0.05,
                   lfc_cutoff = 1) {
  mode <- match.arg(mode)
  if (is.null(mat) || nrow(mat) == 0) stop("No features available to test.")
  if (ncol(mat) != length(groups)) {
    stop("Sample metadata and data matrix are out of sync — re-run sample matching.")
  }

  # Drop samples with no group (and no batch, when batch is in the model)
  keep <- !is.na(groups) & groups != ""
  if (use_batch && !is.null(batch)) keep <- keep & !is.na(batch)
  mat <- mat[, keep, drop = FALSE]
  groups <- droplevels(as.factor(groups[keep]))
  if (!is.null(batch)) batch <- droplevels(as.factor(batch[keep]))

  if (nlevels(groups) < 2) {
    stop("The grouping variable has fewer than two groups among the matched samples.")
  }

  if (mode == "pairwise") {
    if (is.null(group_a) || is.null(group_b)) stop("Select two groups to compare.")
    if (identical(group_a, group_b)) stop("Select two different groups to compare.")
    missing <- setdiff(c(group_a, group_b), levels(groups))
    if (length(missing) > 0) {
      stop(paste0("Group not found in the matched samples: ", paste(missing, collapse = ", ")))
    }
    small <- c(group_a, group_b)[table(groups)[c(group_a, group_b)] < 2]
    if (length(small) > 0) {
      stop(paste0("Need at least 2 samples per group; too few in: ", paste(small, collapse = ", ")))
    }
  }

  design <- de_build_design(groups, if (use_batch) batch else NULL)
  if (qr(design)$rank < ncol(design)) {
    stop(paste0("The model is not identifiable: '", group_var,
                "' and batch are confounded, so their effects cannot be separated. ",
                "Use the ComBat-corrected matrix instead, or choose a grouping variable ",
                "that varies within batches."))
  }
  if (nrow(design) - ncol(design) < 1) {
    stop("Not enough samples to estimate the model — no residual degrees of freedom remain.")
  }

  keys <- de_level_keys(groups)
  contrast_strings <- if (mode == "pairwise") {
    setNames(paste0(keys[[group_a]], "-", keys[[group_b]]),
             paste0(group_a, " vs ", group_b))
  } else {
    ref <- levels(groups)[1]
    others <- levels(groups)[-1]
    setNames(paste0(keys[others], "-", keys[[ref]]),
             paste0(others, " vs ", ref))
  }

  fit <- limma::lmFit(mat, design)
  cm <- limma::makeContrasts(contrasts = unname(contrast_strings), levels = design)
  # topTable() runs its columns through check.names, so the contrast columns are
  # named with syntactic keys here and the readable labels kept alongside.
  contrast_keys <- make.names(names(contrast_strings), unique = TRUE)
  colnames(cm) <- contrast_keys
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = TRUE)

  if (mode == "pairwise") {
    tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "P")
    res <- data.frame(
      feature = rownames(tt),
      log2FC = tt$logFC,
      fold_change = 2^tt$logFC,
      mean_intensity = tt$AveExpr,
      t = tt$t,
      p_value = tt$P.Value,
      fdr = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )
    res$significant <- res$fdr <= fdr_cutoff & abs(res$log2FC) >= lfc_cutoff
    n_up <- sum(res$significant & res$log2FC > 0)
    n_down <- sum(res$significant & res$log2FC < 0)
  } else {
    # With exactly two groups the F-test reduces to the t-test on the single
    # contrast, and limma has no F column to sort on — F is then t squared.
    multi <- ncol(cm) > 1
    tt <- if (multi) {
      limma::topTable(fit2, number = Inf, sort.by = "F")
    } else {
      limma::topTable(fit2, coef = 1, number = Inf, sort.by = "P")
    }
    lfc_mat <- if (multi) {
      as.matrix(tt[, contrast_keys, drop = FALSE])
    } else {
      as.matrix(tt[, "logFC", drop = FALSE])
    }
    colnames(lfc_mat) <- paste0(
      "log2FC_", gsub("[^A-Za-z0-9]+", "_", names(contrast_strings)))
    res <- data.frame(
      feature = rownames(tt),
      mean_intensity = tt$AveExpr,
      F = if (multi) tt$F else tt$t^2,
      p_value = tt$P.Value,
      fdr = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )
    res <- cbind(res, as.data.frame(lfc_mat))
    res$max_abs_log2FC <- apply(abs(lfc_mat), 1, max)
    res$significant <- res$fdr <= fdr_cutoff
    n_up <- NA_integer_
    n_down <- NA_integer_
  }

  label <- if (mode == "pairwise") {
    paste0(group_a, " vs ", group_b)
  } else {
    paste0("Any difference across ", group_var)
  }

  list(
    label = label,
    mode = mode,
    group_var = group_var,
    group_a = group_a,
    group_b = group_b,
    group_levels = levels(groups),
    groups = groups,
    matrix = mat,
    sample_names = colnames(mat),
    batch = batch,
    matrix_source = matrix_source,
    use_batch = use_batch,
    contrasts = contrast_strings,
    fdr_cutoff = fdr_cutoff,
    lfc_cutoff = lfc_cutoff,
    results = res,
    n_features = nrow(res),
    n_samples = ncol(mat),
    n_sig = sum(res$significant),
    n_up = n_up,
    n_down = n_down,
    limma_version = as.character(packageVersion("limma")),
    run_time = Sys.time()
  )
}

# Identifies the matrix a comparison was fitted on, so saved comparisons can be
# dropped when — and only when — the data underneath them changes.
de_mat_fingerprint <- function(m) {
  if (is.null(m)) return(NA_character_)
  paste(nrow(m), ncol(m), format(sum(m, na.rm = TRUE), digits = 15), sep = ":")
}

de_matrix_label <- function(run) {
  if (identical(run$matrix_source, "batch_corrected")) {
    "ComBat batch-corrected data (no batch term in the model)"
  } else if (isTRUE(run$use_batch)) {
    "log2-normalized data with batch as a model covariate"
  } else {
    "log2-normalized data (no batch adjustment)"
  }
}

de_design_summary <- function(run) {
  test <- if (run$mode == "pairwise") "moderated t-test" else "moderated F-test"
  paste0(
    "Comparison: ", run$label, "\n",
    "Test: ", test, " (limma, empirical Bayes)\n",
    "Data: ", de_matrix_label(run), "\n",
    "Samples: ", run$n_samples, " | Features tested: ", run$n_features, "\n",
    "Thresholds: FDR <= ", run$fdr_cutoff,
    if (run$mode == "pairwise") paste0(", |log2FC| >= ", run$lfc_cutoff) else "", "\n",
    "Significant features: ", run$n_sig,
    if (run$mode == "pairwise") paste0(" (", run$n_up, " up, ", run$n_down, " down)") else ""
  )
}

# The raw p-value at which the FDR cutoff bites — drawn on volcano and
# significance plots so the threshold line means what the table means.
de_p_threshold <- function(run) {
  passing <- run$results$p_value[run$results$fdr <= run$fdr_cutoff]
  if (length(passing) == 0) return(NA_real_)
  max(passing, na.rm = TRUE)
}

# ---------------------------
# Plots
# ---------------------------

de_volcano_plot <- function(run, title = NULL) {
  if (run$mode != "pairwise") return(NULL)
  res <- run$results
  status <- ifelse(!res$significant, "Not significant",
                   ifelse(res$log2FC > 0,
                          paste0("Up in ", run$group_a),
                          paste0("Up in ", run$group_b)))
  status <- factor(status, levels = c(paste0("Up in ", run$group_a),
                                      paste0("Up in ", run$group_b),
                                      "Not significant"))
  df <- data.frame(
    log2FC = res$log2FC,
    neg_log10_p = -log10(res$p_value),
    status = status,
    feature = res$feature,
    fdr = res$fdr,
    stringsAsFactors = FALSE
  )
  p_thresh <- de_p_threshold(run)

  shapes <- list(
    list(type = "line", x0 = run$lfc_cutoff, x1 = run$lfc_cutoff,
         y0 = 0, y1 = max(df$neg_log10_p, na.rm = TRUE) * 1.05,
         line = list(dash = "dot", color = "#95a5a6", width = 1)),
    list(type = "line", x0 = -run$lfc_cutoff, x1 = -run$lfc_cutoff,
         y0 = 0, y1 = max(df$neg_log10_p, na.rm = TRUE) * 1.05,
         line = list(dash = "dot", color = "#95a5a6", width = 1))
  )
  annotations <- list()
  if (!is.na(p_thresh)) {
    rng <- range(df$log2FC, na.rm = TRUE)
    shapes[[length(shapes) + 1]] <- list(
      type = "line", x0 = rng[1], x1 = rng[2],
      y0 = -log10(p_thresh), y1 = -log10(p_thresh),
      line = list(dash = "dot", color = "#95a5a6", width = 1))
    annotations[[1]] <- list(
      x = rng[2], y = -log10(p_thresh), xanchor = "right", yanchor = "bottom",
      text = paste0("FDR = ", run$fdr_cutoff), showarrow = FALSE,
      font = list(size = 11, color = "#7f8c8d"))
  }

  plot_ly(df, x = ~log2FC, y = ~neg_log10_p,
          color = ~status,
          colors = c(DE_COLORS$up, DE_COLORS$down, DE_COLORS$ns),
          type = "scatter", mode = "markers",
          marker = list(size = 6, opacity = 0.75),
          text = ~paste0("<b>", feature, "</b><br>log2FC: ", round(log2FC, 3),
                         "<br>p: ", signif(10^(-neg_log10_p), 3),
                         "<br>FDR: ", signif(fdr, 3)),
          hovertemplate = "%{text}<extra></extra>") %>%
    layout(
      title = list(text = if (is.null(title)) run$label else title, x = 0.5),
      xaxis = list(title = "log2 fold change", zeroline = FALSE),
      yaxis = list(title = "-log10(p-value)"),
      shapes = shapes,
      annotations = annotations,
      legend = list(orientation = "h", y = -0.18, x = 0.5, xanchor = "center"),
      margin = list(b = 90)
    )
}

# The F-test has no single fold change to put on an x-axis, so significance is
# shown against mean intensity instead of inventing a volcano.
de_significance_plot <- function(run, title = NULL) {
  res <- run$results
  df <- data.frame(
    mean_intensity = res$mean_intensity,
    neg_log10_p = -log10(res$p_value),
    status = ifelse(res$significant, "Significant", "Not significant"),
    feature = res$feature,
    fdr = res$fdr,
    stringsAsFactors = FALSE
  )
  df$status <- factor(df$status, levels = c("Significant", "Not significant"))
  p_thresh <- de_p_threshold(run)
  shapes <- list()
  annotations <- list()
  if (!is.na(p_thresh)) {
    rng <- range(df$mean_intensity, na.rm = TRUE)
    shapes[[1]] <- list(type = "line", x0 = rng[1], x1 = rng[2],
                        y0 = -log10(p_thresh), y1 = -log10(p_thresh),
                        line = list(dash = "dot", color = "#95a5a6", width = 1))
    annotations[[1]] <- list(
      x = rng[2], y = -log10(p_thresh), xanchor = "right", yanchor = "bottom",
      text = paste0("FDR = ", run$fdr_cutoff), showarrow = FALSE,
      font = list(size = 11, color = "#7f8c8d"))
  }

  plot_ly(df, x = ~mean_intensity, y = ~neg_log10_p,
          color = ~status,
          colors = c(DE_COLORS$up, DE_COLORS$ns),
          type = "scatter", mode = "markers",
          marker = list(size = 6, opacity = 0.75),
          text = ~paste0("<b>", feature, "</b><br>mean log2 intensity: ",
                         round(mean_intensity, 3),
                         "<br>FDR: ", signif(fdr, 3)),
          hovertemplate = "%{text}<extra></extra>") %>%
    layout(
      title = list(text = if (is.null(title)) run$label else title, x = 0.5),
      xaxis = list(title = "Mean log2 intensity"),
      yaxis = list(title = "-log10(p-value)"),
      shapes = shapes,
      annotations = annotations,
      legend = list(orientation = "h", y = -0.18, x = 0.5, xanchor = "center"),
      margin = list(b = 90)
    )
}

# Features worth showing: significant ones, or the top of the table when
# nothing clears the threshold (with the caller told which it got).
de_display_features <- function(run, top_n = DE_MAX_HEATMAP_FEATURES) {
  res <- run$results
  sig <- res$feature[res$significant]
  used_fallback <- length(sig) == 0
  feats <- if (used_fallback) head(res$feature, top_n) else head(sig, top_n)
  list(features = feats, used_fallback = used_fallback,
       n_available = if (used_fallback) nrow(res) else length(sig))
}

de_heatmap_plot <- function(run, top_n = DE_MAX_HEATMAP_FEATURES) {
  sel <- de_display_features(run, top_n)
  feats <- sel$features
  if (length(feats) < 2) return(NULL)

  mat <- run$matrix[feats, , drop = FALSE]
  # Row z-scores: the heatmap shows each feature's pattern across samples, not
  # absolute intensity, which is what makes features comparable in one panel.
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  mat <- sweep(mat, 1, rowMeans(mat, na.rm = TRUE), "-")
  mat <- sweep(mat, 1, ifelse(row_sd > 0, row_sd, 1), "/")

  ord <- order(run$groups, run$sample_names)
  mat <- mat[, ord, drop = FALSE]
  groups_ord <- run$groups[ord]

  if (nrow(mat) > 2) {
    hc <- hclust(dist(mat))
    mat <- mat[hc$order, , drop = FALSE]
  }

  limit <- max(abs(mat), na.rm = TRUE)
  # Group block boundaries and label positions. Names are stripped because a
  # named length-1 numeric serialises to a JSON object, which plotly ignores.
  sizes <- as.numeric(table(groups_ord))
  group_names <- names(table(groups_ord))
  ends <- cumsum(sizes)
  shapes <- lapply(head(ends, -1), function(b) {
    list(type = "line", x0 = b - 0.5, x1 = b - 0.5,
         y0 = -0.5, y1 = nrow(mat) - 0.5,
         line = list(color = "#2c3e50", width = 1.5))
  })
  centers <- ends - sizes / 2
  annotations <- lapply(seq_along(centers), function(i) {
    list(x = centers[i] - 0.5, y = nrow(mat) - 0.5, yanchor = "bottom",
         text = paste0("<b>", group_names[i], "</b>"), showarrow = FALSE,
         font = list(size = 12))
  })

  plot_ly(
    x = colnames(mat), y = rownames(mat), z = mat,
    type = "heatmap",
    colors = colorRamp(c(DE_COLORS$down, "#ffffff", DE_COLORS$up)),
    zmin = -limit, zmax = limit,
    colorbar = list(title = "z-score", len = 0.6),
    hovertemplate = "%{y}<br>%{x}<br>z = %{z:.2f}<extra></extra>"
  ) %>%
    layout(
      title = list(text = paste0(
        if (sel$used_fallback) "Top features by p-value — " else "Significant features — ",
        run$label), x = 0.5),
      xaxis = list(title = "", tickangle = -45, showticklabels = ncol(mat) <= 40),
      yaxis = list(title = "", showticklabels = nrow(mat) <= 60),
      shapes = shapes,
      annotations = annotations,
      margin = list(t = 60)
    )
}

de_feature_boxplot <- function(run, feature, showlegend = TRUE,
                               title = feature, show_y_title = TRUE) {
  if (is.null(feature) || !feature %in% rownames(run$matrix)) return(NULL)
  df <- data.frame(
    value = as.numeric(run$matrix[feature, ]),
    group = run$groups,
    sample = run$sample_names,
    stringsAsFactors = FALSE
  )
  colors <- de_group_colors(levels(run$groups))
  y_title <- if (identical(run$matrix_source, "batch_corrected")) {
    "log2 intensity (batch-corrected)"
  } else {
    "log2 intensity"
  }

  plot_ly(df, x = ~group, y = ~value, color = ~group, colors = colors,
          type = "box", boxpoints = "all", jitter = 0.4, pointpos = 0,
          marker = list(size = 5, opacity = 0.7),
          text = ~sample, hoverinfo = "y+text",
          showlegend = showlegend) %>%
    layout(
      title = if (is.null(title)) NULL else list(text = title, x = 0.5),
      xaxis = list(title = ""),
      yaxis = list(title = if (show_y_title) y_title else "")
    )
}

# Long feature names blow out a small panel title.
de_trunc <- function(x, n = 28) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "\u2026"), x)
}

# A grid of the top features, used in the report so a reader sees the effects
# behind the table without clicking anything.
de_top_boxplots <- function(run, n = DE_REPORT_BOXPLOT_FEATURES) {
  sel <- de_display_features(run, n)
  feats <- head(sel$features, n)
  if (length(feats) == 0) return(NULL)
  # subplot() keeps only one layout title, so each panel is labelled with a
  # paper-referenced annotation, which subplot remaps into that panel's domain.
  plots <- lapply(feats, function(f) {
    de_feature_boxplot(run, f, showlegend = FALSE, title = NULL,
                       show_y_title = FALSE) %>%
      add_annotations(text = de_trunc(f), x = 0.5, y = 1.04,
                      xref = "paper", yref = "paper",
                      xanchor = "center", yanchor = "bottom",
                      showarrow = FALSE, font = list(size = 12))
  })
  n_rows <- ceiling(length(plots) / 3)
  subplot(plots, nrows = n_rows, margin = c(0.03, 0.03, 0.07, 0.07),
          titleY = FALSE) %>%
    layout(showlegend = FALSE, margin = list(t = 40))
}

# ---------------------------
# Report and download content
# ---------------------------

de_methods_paragraph <- function(run) {
  test_sentence <- if (run$mode == "pairwise") {
    paste0("Differential abundance between ", html_escape(run$group_a), " and ",
           html_escape(run$group_b), " was assessed with a moderated t-test.")
  } else {
    paste0("Features differing across any level of ", html_escape(run$group_var),
           " were identified with a moderated F-test.")
  }
  data_sentence <- if (identical(run$matrix_source, "batch_corrected")) {
    paste0("Models were fitted to ComBat batch-corrected log2 intensities; ",
           "batch was therefore not included as a covariate.")
  } else if (isTRUE(run$use_batch)) {
    paste0("Models were fitted to log2-normalized intensities with batch included ",
           "as a covariate, so batch effects were controlled for in the model ",
           "rather than removed from the matrix.")
  } else {
    "Models were fitted to log2-normalized intensities with no batch adjustment."
  }
  sig_sentence <- if (run$mode == "pairwise") {
    paste0("Features were called differentially abundant at a Benjamini-Hochberg FDR of ",
           run$fdr_cutoff, " and an absolute log2 fold change of at least ", run$lfc_cutoff,
           ", giving ", run$n_sig, " features (", run$n_up, " higher in ",
           html_escape(run$group_a), ", ", run$n_down, " higher in ",
           html_escape(run$group_b), ").")
  } else {
    paste0("Features were called significant at a Benjamini-Hochberg FDR of ",
           run$fdr_cutoff, ", giving ", run$n_sig, " features.")
  }

  paste0(
    "Differential analysis was performed in R using the limma package (v",
    run$limma_version, "; Ritchie et al., 2015). ", test_sentence, " ",
    "A linear model was fitted per feature across all ", run$n_samples,
    " samples with one coefficient per level of ", html_escape(run$group_var), ". ",
    data_sentence, " ",
    "Feature-level variances were moderated towards an intensity-dependent trend by ",
    "empirical Bayes shrinkage, and p-values were adjusted for multiple testing across ",
    run$n_features, " features using the Benjamini-Hochberg procedure. ",
    sig_sentence,
    "<br><br>",
    "<strong>Reference:</strong> Ritchie ME, Phipson B, Wu D, et al. ",
    "limma powers differential expression analyses for RNA-sequencing and microarray studies. ",
    "<em>Nucleic Acids Research</em>. 2015;43(7):e47. doi:10.1093/nar/gkv007"
  )
}

# Top rows of a results table, formatted for HTML display.
de_report_table <- function(run, n = DE_REPORT_TABLE_ROWS) {
  res <- head(run$results, n)
  if (run$mode == "pairwise") {
    df <- data.frame(
      Feature = res$feature,
      `log2FC` = fmt_num(res$log2FC),
      `Fold change` = fmt_num(res$fold_change),
      `Mean log2 intensity` = fmt_num(res$mean_intensity),
      `t` = fmt_num(res$t),
      `p-value` = fmt_pval(res$p_value),
      FDR = fmt_pval(res$fdr),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Feature = res$feature,
      `Mean log2 intensity` = fmt_num(res$mean_intensity),
      `F` = fmt_num(res$F),
      `Largest |log2FC|` = fmt_num(res$max_abs_log2FC),
      `p-value` = fmt_pval(res$p_value),
      FDR = fmt_pval(res$fdr),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  }
  df
}

de_report_section <- function(run, index = 1) {
  main_plot <- if (run$mode == "pairwise") de_volcano_plot(run) else de_significance_plot(run)
  heatmap <- de_heatmap_plot(run)
  boxplots <- de_top_boxplots(run)
  sel <- de_display_features(run)

  summary_row <- if (run$mode == "pairwise") {
    paste0("<strong>Significant:</strong> ", run$n_sig, " of ", run$n_features,
           " features (", run$n_up, " higher in ", html_escape(run$group_a),
           ", ", run$n_down, " higher in ", html_escape(run$group_b), ")")
  } else {
    paste0("<strong>Significant:</strong> ", run$n_sig, " of ", run$n_features, " features")
  }

  paste0(
    '<h2>Differential Analysis ', index, ': ', html_escape(run$label), '</h2>',
    '<div class="meta">',
    '<strong>Test:</strong> ',
    if (run$mode == "pairwise") "moderated t-test (limma)" else "moderated F-test (limma)", '<br>',
    '<strong>Grouping variable:</strong> ', html_escape(run$group_var), '<br>',
    '<strong>Data modelled:</strong> ', de_matrix_label(run), '<br>',
    '<strong>Thresholds:</strong> FDR &le; ', run$fdr_cutoff,
    if (run$mode == "pairwise") paste0(", |log2FC| &ge; ", run$lfc_cutoff) else "", '<br>',
    summary_row,
    '</div>',

    '<h3>Methods</h3>',
    '<div class="methods"><p>', de_methods_paragraph(run), '</p></div>',

    '<h3>Top features</h3>',
    '<p>The ', min(DE_REPORT_TABLE_ROWS, nrow(run$results)),
    ' features with the smallest p-values are shown below. The full table is in the ',
    '<code>', html_escape(de_sheet_name(run, index)), '</code> sheet of the Excel download.</p>',
    html_table(de_report_table(run)),

    '<h3>', if (run$mode == "pairwise") "Volcano plot" else "Significance plot", '</h3>',
    '<p>',
    if (run$mode == "pairwise") {
      paste0("Each point is a feature. Coloured points clear both thresholds. ",
             "Positive log2 fold changes are higher in ", html_escape(run$group_a), ".")
    } else {
      "Each point is a feature, plotted against its mean intensity. Coloured points clear the FDR threshold."
    },
    '</p>',
    embed_plotly(main_plot, height = "450px"),

    if (!is.null(heatmap)) paste0(
      '<h3>Heatmap</h3>',
      '<p>',
      if (sel$used_fallback) {
        "No feature cleared the thresholds, so the strongest features by p-value are shown. "
      } else {
        paste0("Showing ", min(length(sel$features), DE_MAX_HEATMAP_FEATURES), " of ",
               sel$n_available, " significant features. ")
      },
      'Values are per-feature z-scores; samples are grouped by ', html_escape(run$group_var), '.</p>',
      embed_plotly(heatmap, height = "600px")
    ) else "",

    if (!is.null(boxplots)) paste0(
      '<h3>Top feature abundances</h3>',
      '<p>Individual sample values for the ',
      min(DE_REPORT_BOXPLOT_FEATURES, length(sel$features)),
      ' strongest features.</p>',
      embed_plotly(boxplots, height = "500px")
    ) else ""
  )
}

# Excel sheet names are capped at 31 characters by the format itself.
de_sheet_name <- function(run, index) {
  safe <- gsub("[^A-Za-z0-9]+", "_", run$label)
  safe <- gsub("^_|_$", "", safe)
  substr(paste0("de", index, "_", safe), 1, 31)
}

de_readme_rows <- function(runs) {
  if (length(runs) == 0) return(NULL)
  rows <- lapply(seq_along(runs), function(i) {
    run <- runs[[i]]
    test <- if (run$mode == "pairwise") "moderated t-test" else "moderated F-test"
    data.frame(
      Section = paste0("SHEET: ", de_sheet_name(run, i)),
      Detail = paste0(
        "Differential analysis results for ", run$label, ". ",
        "limma ", test, " on ", de_matrix_label(run), ". ",
        "Thresholds: FDR <= ", run$fdr_cutoff,
        if (run$mode == "pairwise") paste0(", |log2FC| >= ", run$lfc_cutoff) else "",
        ". Significant features: ", run$n_sig, " of ", run$n_features, ". ",
        "Columns: feature, ",
        if (run$mode == "pairwise") {
          "log2FC (log2 fold change, positive = higher in the first group), fold_change, mean_intensity, t, p_value, fdr, significant."
        } else {
          "mean_intensity, F, p_value, fdr, one log2FC column per group versus the reference level, max_abs_log2FC, significant."
        }
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

de_excel_sheets <- function(runs) {
  if (length(runs) == 0) return(list())
  lapply(seq_along(runs), function(i) {
    list(name = de_sheet_name(runs[[i]], i), data = runs[[i]]$results)
  })
}
