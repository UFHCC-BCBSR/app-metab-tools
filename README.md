# Metabo Tools — Batch Corrector

A Shiny web application for metabolomics data preprocessing and batch correction, developed by the University of Florida Health Cancer Center Biostatistics, Computational Biology, and Bioinformatics Shared Resource (BCBSR) in partnership with the [Southeast Center for Integrated Metabolomics (SECIM)](https://secim.ufl.edu/).

> **This tool is currently under active development and testing.**

## What it does

- Matches count/intensity data to sample metadata
- Filters and imputes missing values
- Applies sample normalization and log2 transformation
- Removes batch effects using ComBat
- Applies post-correction scaling (Pareto or Auto)
- Tests features for differential abundance between groups
- Downloads a complete Excel data package and HTML processing report

## Differential analysis

Comparisons are fitted with [limma](https://bioconductor.org/packages/limma):
one linear model per feature on log2 intensities, with feature variances
moderated across the dataset by empirical Bayes. Two kinds of comparison are
available:

| Comparison | Test | Reports |
|---|---|---|
| Two groups | Moderated t-test | log2 fold change, fold change, t, p-value, BH FDR |
| All groups | Moderated F-test | F, p-value, BH FDR, log2 fold change per group versus the reference |

Batch effects can be handled either way, and the two are mutually exclusive:

- **Log2-normalized data with batch as a covariate** (default) — batch is a term
  in the model, so its effect is controlled for rather than removed.
- **ComBat-corrected data with no batch term** — batch has already been removed
  from the matrix.

Either way the model is fitted to unscaled log2 intensities, so fold changes
stay in log2 intensity units. Pareto/auto scaled matrices are for PCA,
clustering and heatmap display, not for modelling.

Results are shown as a sortable table, a volcano plot (or a significance plot
for the F-test), a heatmap of significant features, and per-feature boxplots.
Every comparison you run is saved: each one becomes its own sheet in the Excel
download and its own section — methods paragraph, top features, and all plots —
in the HTML report. Comparisons are dropped automatically if you re-run
preprocessing or batch correction, so results never outlive the data they were
fitted on.

## Output

The Excel download contains the following data snapshots:

| Sheet | Contents |
|---|---|
| `data_original` | Raw matched matrix |
| `data_preprocessed` | After filtering and imputation |
| `data_normalized` | After sample normalization and log2 transformation |
| `data_normalized_scaled` | Normalized + scaled — use for limma with batch as covariate |
| `data_batch_corrected_scaled` | ComBat corrected + scaled — use for PCA, clustering, limma without batch covariate |
| `metadata` | Full sample metadata with batch column |
| `de1_...`, `de2_...` | One sheet per differential analysis comparison |

## Contact

Hannah Kates — [hkates@ufl.edu](mailto:hkates@ufl.edu)
