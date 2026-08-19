#!/usr/bin/env Rscript
# ===========================
# Standalone slide deck generator
# ===========================
# Produces ONE self-contained .html file: reveal.js, plotly and every figure are
# inlined, so the deck needs no internet at presentation time and nothing beyond
# a static file server to host. Viewers need only a browser.
#
#   Rscript scripts/make_slides.R            # uses CONFIG below
#   Rscript scripts/make_slides.R out.html   # or name the output
#
# The science comes from the app's own modules, so the numbers on the slides are
# the numbers the app produces. Edit CONFIG to point at your data and to choose
# which comparisons become slides.
#
# Build requirements (on the machine rendering the deck, not the server):
#   the app's own packages, plus `revealjs` for its bundled reveal.js assets.

suppressPackageStartupMessages({
  library(shiny); library(plotly); library(limma); library(sva)
  library(impute); library(dplyr); library(tibble); library(RColorBrewer)
})

# Run from the app directory, or from scripts/ — either way, find modules/
find_app_root <- function() {
  for (p in c(getwd(), dirname(getwd()))) {
    if (dir.exists(file.path(p, "modules"))) return(p)
  }
  stop("Could not find modules/. Run this from the app directory: Rscript scripts/make_slides.R")
}
setwd(find_app_root())
for (f in c("report_utils.R", "preprocessing.R", "de_analysis.R")) {
  source(file.path("modules", f))
}

# ===========================
# CONFIG — edit this block
# ===========================
CONFIG <- list(
  title    = "PROCMET metabolomics",
  subtitle = "Processing and differential analysis",
  author   = "UFHCC BCBSR &middot; SECIM",
  output   = "procmet_slides.html",

  # One entry per ionization mode. Drop the second for a single-mode study.
  modes = list(
    list(label = "Positive",
         feature_file = "test-data/PROCMET_Batch1and2_RP_POS_feature_data.csv",
         sample_file  = "test-data/PROCMET_Batch1and2_RP_POS_sample_data.csv",
         feature_col = "feature_id", sample_name_col = "Sample",
         batch_col = "Batch", bio_col = "Group"),
    list(label = "Negative",
         feature_file = "test-data/PROCMET_Batch1and2_RP_NEG_feature_data.csv",
         sample_file  = "test-data/PROCMET_Batch1and2_RP_NEG_sample_data.csv",
         feature_col = "feature_id", sample_name_col = "Sample",
         batch_col = "Batch", bio_col = "Group")
  ),

  # Same options as the app's Preprocess and Batch Correction tabs
  preprocess = list(
    missing_filter = TRUE, imputation = TRUE,
    iqr_filter = FALSE, iqr_threshold = 10,
    row_norm = "SumNorm", scale_norm = "ParetoNorm",
    run_combat = FALSE
  ),

  group_var = "Group",
  # "covariate", "none", or "batch_corrected" (needs run_combat = TRUE)
  matrix_source = "covariate",
  fdr_cutoff = 0.05,
  lfc_cutoff = 1,

  # One slide group per comparison
  comparisons = list(
    list(type = "pairwise", a = "IVA", b = "NORM"),
    list(type = "pairwise", a = "HCY", b = "NORM"),
    list(type = "global", groups = c("NORM", "IVA", "HCY"))
  ),

  # Figures to include per comparison
  include_heatmaps = TRUE,
  include_boxplots = TRUE
)

# ===========================
# Run the pipeline
# ===========================
run_mode <- function(spec, pp) {
  message("  processing ", spec$label, " mode")
  feat <- read.csv(spec$feature_file, stringsAsFactors = FALSE)
  samp <- read.csv(spec$sample_file, stringsAsFactors = FALSE)
  matched <- match_samples(feat, samp, spec$feature_col, spec$sample_name_col,
                           spec$batch_col, spec$bio_col, NULL)
  if (is.null(matched)) stop(spec$label, ": no samples matched.")

  mat <- apply(matched$count_matrix, 2, as.numeric)
  rownames(mat) <- rownames(matched$count_matrix)
  n_original <- nrow(mat)
  if (isTRUE(pp$missing_filter)) mat <- filter_missing(mat, threshold = 0.5)
  if (isTRUE(pp$imputation)) { mat[mat == 0] <- NA; mat <- impute_knn(mat) }
  if (isTRUE(pp$iqr_filter)) mat <- filter_iqr(mat, top_percent = pp$iqr_threshold)
  mat <- sample_normalize(mat, pp$row_norm)
  if (any(mat <= 0, na.rm = TRUE)) mat[mat <= 0] <- 0.01
  normalized <- log2(mat)

  corrected <- NULL
  if (isTRUE(pp$run_combat) && has_usable_batch(matched)) {
    mod <- if (isTRUE(matched$has_bio_var)) {
      model.matrix(~as.factor(matched$sample_data$biological_var))
    } else NULL
    corrected <- ComBat(dat = normalized,
                        batch = as.factor(matched$sample_data$batch),
                        mod = mod, par.prior = TRUE)
  }

  pca_before <- perform_pca(normalized, matched$sample_data,
                            has_bio_var = isTRUE(matched$has_bio_var),
                            bio_var_name = spec$bio_col,
                            has_batch = has_usable_batch(matched))
  list(label = spec$label, matched_data = matched, normalized = normalized,
       batch_corrected = corrected, has_batch = has_usable_batch(matched),
       pca_before = pca_before, n_original = n_original, n_final = nrow(normalized),
       n_samples = ncol(normalized),
       params = list(bio_var = spec$bio_col, batch_var = spec$batch_col))
}

# ===========================
# Deck assembly
# ===========================
inline_js <- function(path) {
  js <- paste(readLines(path, warn = FALSE), collapse = "\n")
  # a literal </script> inside inlined JS would close the tag early
  gsub("</script", "<\\\\/script", js, fixed = TRUE)
}
inline_css <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")

# The embedding div carries the dimensions and plotly sizes itself to its
# container. Margins are left alone: the plot builders already reserve room for
# their legends, and overriding that pushes the legend onto the axis title. The
# figure's own title goes, though — the slide heading already says it.
slide_fig <- function(fig) {
  if (is.null(fig)) return(NULL)
  # A slide is wider than it is tall, so the legend goes to the right: the
  # bottom legend the app uses needs vertical room a slide does not have, and
  # collides with the axis title.
  layout(fig, title = list(text = ""),
         legend = list(orientation = "v", x = 1.02, xanchor = "left",
                       y = 1, yanchor = "top"),
         margin = list(r = 150, b = 70))
}

slide <- function(...) paste0("<section>", paste0(..., collapse = ""), "</section>")

fig_slide <- function(title, fig, note = NULL, height = 470) {
  if (is.null(fig)) return("")
  slide(
    "<h2>", title, "</h2>",
    embed_plotly(slide_fig(fig), height = paste0(height, "px")),
    if (!is.null(note)) paste0('<p class="note">', note, "</p>") else ""
  )
}

build_deck <- function(modes, runs, cfg) {
  rj <- system.file(package = "revealjs")
  rj_dir <- list.files(rj, pattern = "^reveal\\.js-", full.names = TRUE)[1]
  if (is.na(rj_dir)) stop("reveal.js assets not found; install.packages('revealjs')")
  plotly_js <- list.files(system.file("htmlwidgets/lib", package = "plotly"),
                          pattern = "plotly.*min\\.js$", recursive = TRUE, full.names = TRUE)[1]

  multi <- length(modes) > 1
  n_feat <- sum(vapply(modes, function(m) m$n_final, numeric(1)))

  overview <- slide(
    "<h2>What was done</h2>",
    "<table class='tbl'><tr><th>Mode</th><th>Samples</th><th>Features in</th><th>Features analysed</th></tr>",
    paste0(vapply(modes, function(m) {
      paste0("<tr><td>", m$label, "</td><td>", m$n_samples, "</td><td>",
             m$n_original, "</td><td>", m$n_final, "</td></tr>")
    }, character(1)), collapse = ""),
    "</table>",
    "<ul>",
    "<li>Each ionization mode filtered, imputed and normalized separately",
    if (multi) "<li>Modes combined only for statistics, with one FDR correction across both" else "",
    "<li>", if (isTRUE(cfg$preprocess$run_combat)) "ComBat batch correction applied" else
      "No batch correction applied", "</li>",
    "</ul>")

  pca_slides <- paste0(vapply(modes, function(m) {
    fig_slide(paste0("Quality control &mdash; ", m$label, " mode"),
              create_pca_plots(m$pca_before, ""),
              note = if (isTRUE(m$pca_before$has_batch))
                "Left: biological group. Right: batch. Overlapping batches mean no batch structure to remove."
              else NULL,
              height = 450)
  }, character(1)), collapse = "")

  comparison_slides <- paste0(vapply(seq_along(runs), function(i) {
    run <- runs[[i]]
    main <- if (run$mode == "pairwise") de_volcano_plot(run) else de_significance_plot(run)
    headline <- if (run$mode == "pairwise") {
      paste0("<strong>", run$n_sig, "</strong> features differ (", run$n_up, " higher in ",
             run$group_a, ", ", run$n_down, " higher in ", run$group_b, ")")
    } else {
      paste0("<strong>", run$n_sig, "</strong> of ", run$n_features,
             " features differ across ", paste(run$groups_included, collapse = ", "))
    }
    out <- fig_slide(run$label, main, note = headline)

    if (isTRUE(cfg$include_heatmaps)) {
      out <- paste0(out, paste0(vapply(run$mode_labels, function(ml) {
        fig_slide(paste0(run$label, if (multi) paste0(" &mdash; ", ml, " mode") else ""),
                  de_heatmap_plot(run, mode_label = ml), height = 500)
      }, character(1)), collapse = ""))
    }
    if (isTRUE(cfg$include_boxplots)) {
      out <- paste0(out, fig_slide(paste0(run$label, " &mdash; strongest features"),
                                   de_top_boxplots(run), height = 480))
    }
    out
  }, character(1)), collapse = "")

  closing <- slide(
    "<h2>Notes</h2><ul>",
    "<li>Counts are features, not compounds: a metabolite seen in both modes is tested twice</li>",
    "<li>Multi-group tests carry no fold change threshold, so their counts run higher</li>",
    "<li>Thresholds: FDR &le; ", cfg$fdr_cutoff, ", |log2FC| &ge; ", cfg$lfc_cutoff, "</li>",
    "<li>Full detail, methods paragraph and complete tables are in the HTML report and Excel package</li>",
    "</ul>")

  paste0(
    "<!DOCTYPE html>\n<html lang='en'><head><meta charset='UTF-8'>\n",
    "<title>", cfg$title, "</title>\n",
    "<style>", inline_css(file.path(rj_dir, "dist", "reset.css")), "</style>\n",
    "<style>", inline_css(file.path(rj_dir, "dist", "reveal.css")), "</style>\n",
    "<style>", inline_css(file.path(rj_dir, "dist", "theme", "white.css")), "</style>\n",
    "<style>
      /* the reveal theme uppercases headings; feature names and group labels
         are case-sensitive, so that is turned off */
      .reveal h1, .reveal h2, .reveal h3 { text-transform: none; }
      .reveal h1 { font-size: 1.9em; } .reveal h2 { font-size: 1.25em; color: #2980b9; }
      .reveal section { text-align: left; }
      .reveal .tbl { border-collapse: collapse; font-size: 0.6em; margin-bottom: 0.6em; }
      .reveal .tbl th { background: #2980b9; color: #fff; padding: 6px 12px; text-align: left; }
      .reveal .tbl td { padding: 5px 12px; border-bottom: 1px solid #ddd; }
      .reveal ul { font-size: 0.62em; }
      .reveal .note { font-size: 0.55em; color: #555; margin-top: 0.3em; }
      .reveal .title-slide { text-align: center; }
      .reveal .subtitle { color: #555; font-size: 0.8em; }
    </style>\n",
    "<script>", inline_js(plotly_js), "</script>\n",
    "</head><body><div class='reveal'><div class='slides'>\n",
    slide("<div class='title-slide'><h1>", cfg$title, "</h1>",
          "<p class='subtitle'>", cfg$subtitle, "</p>",
          "<p class='subtitle'>", cfg$author, "</p>",
          "<p class='subtitle'>", format(Sys.Date(), "%d %B %Y"), "</p></div>"),
    overview, pca_slides, comparison_slides, closing,
    "\n</div></div>\n",
    "<script>", inline_js(file.path(rj_dir, "dist", "reveal.js")), "</script>\n",
    "<script>
      Reveal.initialize({ hash: true, controls: true, progress: true,
                          center: false, transition: 'none', width: 1050, height: 700,
                          margin: 0.05, minScale: 0.2, maxScale: 1.0 });
      // reveal hides off-screen slides, so plots drawn while hidden need a nudge
      Reveal.on('slidechanged', function (e) {
        e.currentSlide.querySelectorAll('.js-plotly-plot').forEach(function (d) {
          Plotly.Plots.resize(d);
        });
      });
    </script>\n",
    "</body></html>\n")
}

# ===========================
# Main
# ===========================
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1) args[1] else CONFIG$output

message("Processing ", length(CONFIG$modes), " mode(s)...")
modes <- lapply(CONFIG$modes, run_mode, pp = CONFIG$preprocess)

message("Running ", length(CONFIG$comparisons), " comparison(s)...")
runs <- lapply(CONFIG$comparisons, function(cmp) {
  if (identical(cmp$type, "pairwise")) {
    message("  ", cmp$a, " vs ", cmp$b)
    de_run(modes, CONFIG$group_var, "pairwise", group_a = cmp$a, group_b = cmp$b,
           matrix_source = CONFIG$matrix_source,
           fdr_cutoff = CONFIG$fdr_cutoff, lfc_cutoff = CONFIG$lfc_cutoff)
  } else {
    message("  F-test: ", paste(cmp$groups, collapse = " / "))
    de_run(modes, CONFIG$group_var, "global", groups_included = cmp$groups,
           matrix_source = CONFIG$matrix_source,
           fdr_cutoff = CONFIG$fdr_cutoff, lfc_cutoff = 0)
  }
})

message("Building deck...")
writeLines(build_deck(modes, runs, CONFIG), out_file)
message("Wrote ", out_file, " (", round(file.size(out_file) / 1024^2, 1), " MB)")
