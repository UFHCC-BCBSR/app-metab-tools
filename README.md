# Metabo Tools — Batch Corrector

A Shiny web application for metabolomics data preprocessing and batch correction, developed by the University of Florida Health Cancer Center Biostatistics, Computational Biology, and Bioinformatics Shared Resource (BCBSR) in partnership with the [Southeast Center for Integrated Metabolomics (SECIM)](https://secim.ufl.edu/).

> **This tool is currently under active development and testing.**

## What it does

- Matches count/intensity data to sample metadata
- Filters and imputes missing values
- Applies sample normalization and log2 transformation
- Removes batch effects using ComBat
- Applies post-correction scaling (Pareto or Auto)
- Downloads a complete Excel data package and HTML processing report

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

## Contact

Hannah Kates — [hkates@ufl.edu](mailto:hkates@ufl.edu)
