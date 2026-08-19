# Metabo Tools — Batch Corrector

A Shiny web application for metabolomics data preprocessing and batch correction, developed by the University of Florida Health Cancer Center Biostatistics, Computational Biology, and Bioinformatics Shared Resource (BCBSR) in partnership with the [Southeast Center for Integrated Metabolomics (SECIM)](https://secim.ufl.edu/).

> **This tool is currently under active development and testing.**

## What it does

- Matches count/intensity data to sample metadata
- Filters and imputes missing values
- Applies sample normalization and log2 transformation
- Removes batch effects using ComBat, when there are batch effects to remove
- Applies scaling (Pareto or Auto)
- Tests features for differential abundance between groups
- Handles positive and negative ionization mode together in one session
- Downloads a complete Excel data package and HTML processing report

## Ionization modes

Positive and negative mode acquisitions are separate experiments. Intensities
are not comparable between them, missing-value patterns differ, and each mode
has its own injection sequence — so each mode is matched, filtered, imputed,
normalized and (optionally) batch corrected entirely on its own. PCA is shown
per mode for the same reason: a combined PCA is dominated by whichever mode
contributes more features.

The modes come together only for differential analysis, where each is modeled
separately and the p-values from both are then adjusted **once**, together. That
keeps variance moderation inside a single noise regime while making the false
discovery rate a property of the experiment rather than of one acquisition.

Each mode has its own feature table and its own sample metadata, so the two
files need not share sample naming or batch spelling.

> Note: a metabolite that ionizes in both modes appears as a feature in each,
> so it is tested twice. Deduplicating those requires compound annotation the
> app does not have. Results carry a `mode` column; treat cross-mode duplicates
> as related rather than independent findings.

## Batch correction is optional

Batch correction is a step you should be able to decline. If the PCA shows no
separation by batch, there is nothing to remove, and running ComBat anyway
models structure that is not there.

- **No batch column at all** — uncheck "My data has a batch variable". The
  batch correction step then only applies scaling.
- **A batch column but no batch effect** — leave the batch variable in place and
  uncheck "Remove batch effects with ComBat". Scaling is still applied, and
  differential analysis can still include batch as a covariate.
- **A real batch effect** — run ComBat as before.

## Differential analysis

Comparisons are fitted with [limma](https://bioconductor.org/packages/limma):
one linear model per feature on log2 intensities, with feature variances
moderated across the dataset by empirical Bayes. Two kinds of comparison are
available:

| Comparison | Test | Reports |
|---|---|---|
| Two groups | Moderated t-test | log2 fold change, fold change, t, p-value, BH FDR |
| All groups | Moderated F-test | F, p-value, BH FDR, log2 fold change per group versus the reference |

Pick which groups go into a comparison: a two-group test selects the two, and
the F-test takes any subset of levels. Groups you leave out are excluded from
the model entirely, so they do not contribute to the variance estimate. Samples
whose group label is blank — process blanks, for instance — are always excluded.

Batch is handled in one of three ways, and the first and third are mutually
exclusive:

- **Log2-normalized data with batch as a covariate** (default when a batch
  variable exists) — batch is a term in the model, so its effect is controlled
  for rather than removed.
- **Log2-normalized data with no batch adjustment** — for data with no batch
  variable, or no batch effect worth modeling.
- **ComBat-corrected data with no batch term** — offered once ComBat has been
  run for every loaded mode; batch has already been removed from the matrix.

Either way the model is fitted to unscaled log2 intensities, so fold changes
stay in log2 intensity units. Pareto/auto scaled matrices are for PCA,
clustering and heatmap display, not for modeling.

Results are shown as a sortable table, a volcano plot (or a significance plot
for the F-test), a PCA of the samples in the comparison, a heatmap of
significant features, and per-feature boxplots.

The comparison PCA is separate from the quality-control PCA on the Preprocess
tab. That one covers every sample, which is what you want for spotting batch
structure; this one is restricted to the groups being compared, so it is not
dominated by blanks and untested groups. It is computed on all features, not
the significant ones — selecting features for differing between these groups
and then plotting those features would separate the groups no matter what the
data said.
With two modes loaded, the volcano marks each mode with its own symbol, and each
mode gets its own heatmap because positive and negative mode were run as two
separate experiments on different sample injections.

For a comparison across three or more groups there is no volcano plot to draw:
a volcano needs one fold change per feature, which only exists when two groups
are compared. The significance plot shown instead places each feature by its
average abundance and by the strength of the evidence that it differs somewhere
among the groups. The app explains this in full alongside the plot.
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
| `data_batch_corrected_scaled` | ComBat corrected + scaled — only present when ComBat was run |
| `metadata` | Full sample metadata with batch column |
| `de1_...`, `de2_...` | One sheet per differential analysis comparison |

With two modes loaded, each mode's sheets carry a prefix taken from its name
(`POS_data_normalized`, `NEG_metadata`, and so on). With a single mode the names
are unprefixed, as above.

## Contact

Hannah Kates — [hkates@ufl.edu](mailto:hkates@ufl.edu)
