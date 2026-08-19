# ===========================
# PLS-DA with validation
# ===========================
# PLS-DA is supervised: it is handed the class labels and finds the projection
# that best reproduces them. With thousands of features and tens of samples it
# will separate randomly assigned labels, so a scores plot on its own is not
# evidence of anything. Everything here is built so the validation travels with
# the picture and cannot be quietly dropped:
#
#   R2Y   how much of the class information the model captures (fit)
#   Q2    how well it predicts samples held out in cross-validation
#   pQ2   how often models on shuffled labels reach that Q2 (permutation test)
#
# ropls does the cross-validation and permutation as part of fitting, and
# declines to build a model at all when no component survives cross-validation.
# That refusal is a result, not an error, and is reported as one.

PLSDA_DEFAULT_PERMUTATIONS <- 200
PLSDA_VIP_SHOWN <- 20
PLSDA_Q2_GOOD <- 0.4      # conventional "acceptable predictive ability" in metabolomics
PLSDA_P_CUTOFF <- 0.05

# Fit one mode of one comparison. Returns a list describing what happened,
# including the case where no model could be built.
plsda_fit <- function(run, mode_label = NULL,
                      n_perm = PLSDA_DEFAULT_PERMUTATIONS, seed = 1) {
  if (is.null(mode_label)) mode_label <- run$mode_labels[1]
  fit <- de_fit_for_mode(run, mode_label)
  if (is.null(fit)) return(NULL)

  X <- t(fit$matrix)
  y <- droplevels(fit$groups)
  if (nrow(X) < 6 || nlevels(y) < 2) {
    return(list(ok = FALSE, status = "too_small", mode_label = mode_label,
                n_samples = nrow(X), n_features = ncol(X), groups = levels(y),
                n_perm = n_perm))
  }

  set.seed(seed)
  model <- tryCatch(
    suppressWarnings(suppressMessages(
      ropls::opls(X, y, predI = NA, permI = n_perm,
                  crossvalI = min(7, nrow(X)),
                  fig.pdfC = "none", info.txtC = "none")
    )),
    error = function(e) e
  )

  base <- list(mode_label = mode_label, n_samples = nrow(X), n_features = ncol(X),
               groups = levels(y), group_vector = y, sample_names = rownames(X),
               n_perm = n_perm, group_var = run$group_var, label = run$label)

  if (inherits(model, "error")) {
    return(c(base, list(ok = FALSE, status = "error", message = conditionMessage(model))))
  }
  summ <- ropls::getSummaryDF(model)
  # ropls returns an empty summary when cross-validation kept no component
  if (is.null(summ) || nrow(summ) == 0) {
    return(c(base, list(ok = FALSE, status = "no_component")))
  }

  scores <- ropls::getScoreMN(model)
  vip <- tryCatch(ropls::getVipVn(model), error = function(e) numeric(0))

  q2 <- summ[["Q2(cum)"]]
  pq2 <- summ[["pQ2"]]
  status <- if (is.na(q2) || q2 <= 0 || (!is.na(pq2) && pq2 > PLSDA_P_CUTOFF)) {
    "not_predictive"
  } else if (q2 < PLSDA_Q2_GOOD) {
    "weak"
  } else {
    "predictive"
  }

  c(base, list(
    ok = TRUE, status = status,
    n_components = summ[["pre"]],
    r2x = summ[["R2X(cum)"]], r2y = summ[["R2Y(cum)"]], q2 = q2,
    p_r2y = summ[["pR2Y"]], p_q2 = pq2,
    scores = scores, vip = vip
  ))
}

# Permutation p-values cannot resolve below 1/(n+1).
plsda_fmt_p <- function(p, n_perm) {
  if (is.null(p) || is.na(p)) return("n/a")
  floor_p <- 1 / (n_perm + 1)
  if (p < floor_p) paste0("< ", signif(floor_p, 2)) else formatC(p, format = "g", digits = 2)
}

# Plain language, because "Q2 = -0.02" is not a conclusion most readers can act on.
plsda_verdict <- function(fit) {
  if (is.null(fit)) return(list(level = "bad", headline = "PLS-DA not available.", detail = ""))
  if (isTRUE(fit$ok)) {
    p_txt <- plsda_fmt_p(fit$p_q2, fit$n_perm)
    stats <- paste0("R2Y = ", signif(fit$r2y, 3), ", Q2 = ", signif(fit$q2, 3),
                    ", permutation p = ", p_txt, " (", fit$n_components,
                    if (fit$n_components == 1) " component, " else " components, ",
                    fit$n_perm, " permutations).")
    switch(fit$status,
      "predictive" = list(level = "good",
        headline = "The model predicts held-out samples and beats shuffled labels.",
        detail = paste0(stats, " The separation below reflects real structure, not overfitting.")),
      "weak" = list(level = "warn",
        headline = "The model predicts better than chance, but weakly.",
        detail = paste0(stats, " Q2 below ", PLSDA_Q2_GOOD,
                        " is usually treated as suggestive rather than conclusive; the separation is real but modest.")),
      "not_predictive" = list(level = "bad",
        headline = "This model does not predict better than chance.",
        detail = paste0(stats, " The separation in the scores plot is the algorithm fitting noise. ",
                        "Do not read a group difference from it."))
    )
  } else {
    switch(fit$status,
      "no_component" = list(level = "bad",
        headline = "Cross-validation kept no component: these groups are not separable beyond chance.",
        detail = paste0("PLS-DA was fitted on ", fit$n_features, " features across ", fit$n_samples,
                        " samples and no component survived cross-validation. ",
                        "There is no scores plot to show, because any plot would be showing noise.")),
      "too_small" = list(level = "bad",
        headline = "Too few samples for a cross-validated PLS-DA.",
        detail = paste0("This comparison has ", fit$n_samples,
                        " samples; cross-validation and permutation testing need more.")),
      list(level = "bad", headline = "PLS-DA could not be fitted.",
           detail = if (!is.null(fit$message)) fit$message else ""))
  }
}

PLSDA_LEVEL_COLORS <- list(
  good = list(bg = "#d5f5e3", border = "#27ae60"),
  warn = list(bg = "#fef9e7", border = "#f39c12"),
  bad  = list(bg = "#fdedec", border = "#e74c3c")
)

plsda_verdict_html <- function(fit) {
  v <- plsda_verdict(fit)
  col <- PLSDA_LEVEL_COLORS[[v$level]]
  paste0('<div style="background:', col$bg, '; border-left: 5px solid ', col$border,
         '; padding: 12px 15px; border-radius: 4px; margin: 10px 0;">',
         '<strong>', html_escape(v$headline), '</strong><br>',
         '<span style="font-size: 0.92em;">', html_escape(v$detail), '</span></div>')
}

plsda_scores_plot <- function(fit) {
  if (is.null(fit) || !isTRUE(fit$ok)) return(NULL)
  scores <- fit$scores
  groups <- fit$group_vector
  colors <- de_group_colors(levels(groups))

  if (ncol(scores) >= 2) {
    df <- data.frame(t1 = scores[, 1], t2 = scores[, 2],
                     group = groups, sample = fit$sample_names)
    plot_ly(df, x = ~t1, y = ~t2, color = ~group, colors = colors,
            type = "scatter", mode = "markers",
            marker = list(size = 9, opacity = 0.85),
            text = ~paste0("<b>", sample, "</b><br>", fit$group_var, ": ", group),
            hovertemplate = "%{text}<extra></extra>") %>%
      layout(title = list(text = paste0("PLS-DA scores — ", fit$mode_label, " mode"), x = 0.5),
             xaxis = list(title = "Component 1"), yaxis = list(title = "Component 2"),
             legend = list(orientation = "h", y = -0.18, x = 0.5, xanchor = "center"),
             margin = list(b = 90))
  } else {
    # Cross-validation kept a single component; there is no second axis to plot
    # against, so the one component is shown by group.
    df <- data.frame(t1 = scores[, 1], group = groups, sample = fit$sample_names)
    plot_ly(df, x = ~group, y = ~t1, color = ~group, colors = colors,
            type = "box", boxpoints = "all", jitter = 0.4, pointpos = 0,
            marker = list(size = 6, opacity = 0.75),
            text = ~sample, hoverinfo = "y+text") %>%
      layout(title = list(text = paste0("PLS-DA component 1 by group — ",
                                        fit$mode_label, " mode"), x = 0.5),
             xaxis = list(title = ""), yaxis = list(title = "Component 1"),
             showlegend = FALSE)
  }
}

plsda_vip_plot <- function(fit, n = PLSDA_VIP_SHOWN) {
  if (is.null(fit) || !isTRUE(fit$ok) || length(fit$vip) == 0) return(NULL)
  v <- sort(fit$vip, decreasing = TRUE)
  v <- head(v, n)
  df <- data.frame(feature = factor(names(v), levels = rev(names(v))),
                   vip = as.numeric(v))
  plot_ly(df, x = ~vip, y = ~feature, type = "bar", orientation = "h",
          marker = list(color = "#2980b9"),
          hovertemplate = "%{y}<br>VIP = %{x:.2f}<extra></extra>") %>%
    layout(title = list(text = paste0("Top ", nrow(df), " features by VIP — ",
                                      fit$mode_label, " mode"), x = 0.5),
           xaxis = list(title = "VIP score"),
           yaxis = list(title = "", automargin = TRUE),
           shapes = list(list(type = "line", x0 = 1, x1 = 1, y0 = -0.5, y1 = nrow(df) - 0.5,
                              line = list(dash = "dot", color = "#7f8c8d", width = 1))),
           margin = list(l = 10))
}

plsda_vip_table <- function(fit, n = PLSDA_VIP_SHOWN) {
  if (is.null(fit) || !isTRUE(fit$ok) || length(fit$vip) == 0) return(NULL)
  v <- sort(fit$vip, decreasing = TRUE)
  data.frame(Feature = names(head(v, n)),
             VIP = as.numeric(head(v, n)),
             check.names = FALSE, stringsAsFactors = FALSE)
}

plsda_help_text <- function() {
  paste0(
    "PLS-DA is supervised: unlike PCA, it is told which group each sample belongs to and ",
    "finds the projection that best separates them. That makes it good at showing class ",
    "structure and prone to inventing it — with thousands of features it can separate ",
    "randomly assigned labels. The numbers above are what decide whether the picture means ",
    "anything: Q2 is how well the model predicts samples held out of fitting, and the ",
    "permutation p-value is how often shuffled labels do as well. VIP scores rank features ",
    "by their contribution to the separation; above 1 is the usual threshold."
  )
}

# Report section for one comparison, covering whichever modes were fitted.
plsda_report_section <- function(run) {
  fits <- run$plsda
  if (is.null(fits) || length(fits) == 0) return("")
  multi <- length(run$mode_labels) > 1
  paste0(
    "<h3>PLS-DA</h3>",
    "<p>", plsda_help_text(), "</p>",
    paste0(vapply(names(fits), function(ml) {
      fit <- fits[[ml]]
      body <- paste0(
        if (multi) paste0("<h4>", html_escape(ml), " mode</h4>") else "",
        plsda_verdict_html(fit))
      sp <- plsda_scores_plot(fit)
      vp <- plsda_vip_plot(fit)
      if (!is.null(sp)) body <- paste0(body, embed_plotly(sp, height = "480px"))
      if (!is.null(vp)) body <- paste0(body, embed_plotly(vp, height = "520px"))
      body
    }, character(1)), collapse = "")
  )
}
