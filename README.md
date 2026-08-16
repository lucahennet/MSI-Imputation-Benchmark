# MSI-Imputation-Benchmark

A benchmarking framework for missing-value imputation in Mass Spectrometry Imaging (MSI) data. It simulates realistic missingness (MCAR, MNAR, spatial dropout, and intensity-coupled hybrids), runs a panel of imputation methods (classical, matrix-factorisation, and spatial), and scores each method on value-recovery accuracy, spatial-structure preservation, and downstream biological validity (PCA, differential analysis, clustering, co-abundance networks).

## Contents

- [Overview](#overview)
- [Pipeline](#pipeline)
- [Repository structure](#repository-structure)
- [Missingness mechanisms](#missingness-mechanisms)
- [Imputation methods](#imputation-methods)
- [Evaluation metrics](#evaluation-metrics)
- [Downstream validation](#downstream-validation)
- [Requirements](#requirements)
- [Usage](#usage)
- [Data layout](#data-layout)
- [Output](#output)
- [Reproducibility](#reproducibility)

## Overview

MSI produces high-dimensional pixel × feature matrices (ion intensities across a tissue section) that are commonly affected by missing values — either technical (below limit of detection) or structural (values dropped for zero-inflation, thresholding, etc.). Choosing an imputation strategy is non-trivial: methods that work well for bulk omics data don't necessarily preserve the *spatial* structure that makes MSI useful.

This project:

1. Loads real, fully-observed MSI datasets (spatial ROI exports or non-spatial Excel exports) to serve as ground truth.
2. Simulates missingness under several mechanisms and proportions.
3. Imputes the missing values with 10+ methods (R-native and external/Python).
4. Scores imputation quality via value-accuracy, spatial-fidelity, and spectral-preservation metrics.
5. Propagates each imputed matrix through downstream analyses (PCA, differential expression, clustering, co-abundance) and compares the results back to the ground-truth analysis.

## Pipeline

```
01_data(_nospatial).R   Load & assemble raw MSI data → wide matrix + coordinates
02_simulation.R         Simulate missingness (MCAR / MNAR / hybrid / spatial)
03_imputation.R         Imputation methods (classical, PCA-based, spatial)
04_preprocessing.R      Transform + scale before imputation, and reverse after
05_metrics.R            Accuracy & spatial-fidelity metrics
06_experiment.R         Orchestrates run → export → (external) import → score
07_visualisation.R      Benchmark plots (line/bar, runtime, spatial fidelity)
07.2_plots_manuscript.R Manuscript-styled variants of the above (multi-dataset)
08_downstream.R         PCA / DE / clustering / co-abundance vs. ground truth
09_downstream_visualisation.R  Plots for the downstream comparisons
10.1-10.3                Mixed-mechanism (MCAR+MNAR) benchmark on tissue zones
```

Each stage is a standalone script sourced by one of the three driver scripts (`main.R`, `main_nospatial.R`, `main_mix-missingness.R`) described in [Usage](#usage).

## Repository structure

```
MSI-Imputation-Benchmark/
├── main.R                       # Spatial pipeline (ROI-based MSI data)
├── main_nospatial.R             # Non-spatial pipeline (Excel MSI export)
├── main_mix-missingness.R       # Mixed MCAR/MNAR benchmark on tissue zones
├── scripts/R/
│   ├── 01_data.R                 # Assemble spatial MSI data (TXT coords + CSV features)
│   ├── 01_data_nospatial.R       # Load non-spatial MSI data from Excel
│   ├── 02_simulation.R           # Missingness simulators (MCAR/MNAR/hybrid/spatial)
│   ├── 03_imputation.R           # Imputation methods
│   ├── 04_preprocessing.R        # Transform/scale + inverse, applied around imputation
│   ├── 05_metrics.R              # Metric functions (RMSE, CCC, SAM, SSIM, Moran's I, ...)
│   ├── 06_experiment.R           # Experiment runner + Python interop (export/import)
│   ├── 07_visualisation.R        # Exploratory benchmark plots
│   ├── 07.2_plots_manuscript.R   # Publication-styled plots (multi-dataset facets)
│   ├── 08_downstream.R           # Downstream analyses + ground-truth comparison
│   ├── 09_downstream_visualisation.R  # Downstream comparison plots
│   ├── 10.1_mix_data_analysis.R  # Missingness structure analysis (MASH tissue zones)
│   ├── 10.2_mix_data_benchmark.R # RF vs QRILC vs Oracle benchmark on zones
│   └── 10.3_mix_data_plots.R     # Plots for the mixed-mechanism benchmark
├── data/                         # Input datasets (gitignored)
├── output/                       # Exported CSV/JSON for external methods (gitignored)
└── results/                      # Saved .rds benchmark outputs (gitignored)
```

## Missingness mechanisms

Implemented in [02_simulation.R](scripts/R/02_simulation.R):

| Function | Mechanism |
|---|---|
| `simulate_MCAR` | Missing completely at random, per feature |
| `simulate_MCAR_global` | MCAR sampled across the whole matrix at once |
| `simulate_MNAR` | Lowest-intensity values per feature are dropped (below-LOD proxy) |
| `simulate_hybrid_missing` | Weighted mix of MNAR + MCAR per feature |
| `simulate_msi_missing` | MNAR + spatial dropout (contiguous region) + MCAR remainder |

[10.2_mix_data_benchmark.R](scripts/R/10.2_mix_data_benchmark.R) adds a fourth, intensity-coupled variant: each *feature* is probabilistically assigned to MCAR or MNAR via a logistic function of its mean intensity, calibrated so the realised Spearman coupling matches an empirically observed target (`target_rho_from_summary()` in [10.1](scripts/R/10.1_mix_data_analysis.R)).

## Imputation methods

Implemented in [03_imputation.R](scripts/R/03_imputation.R):

| Method | Function | Notes |
|---|---|---|
| Zero / Mean / Median | `impute_zero`, `impute_mean`, `impute_median` | Naive baselines |
| Half-minimum | `impute_half_min` | Common LOD-substitution baseline |
| Random Forest | `impute_missForest` | via `missForest` |
| k-NN | `impute_knn` | via `VIM::kNN` |
| QRILC | `impute_qrilc` | Quantile regression, left-censored — MNAR-aware |
| PPCA / BPCA | `impute_ppca`, `impute_bpca` | Probabilistic / Bayesian PCA via `pcaMethods` |
| SVD | `impute_svd` | Iterative SVD via `pcaMethods` |
| Spatial k-NN | `impute_spatial_knn` | k nearest pixels by (X, Y) distance |
| Gaussian Process | `impute_gp` | Per-feature GP regression on coordinates |
| IDW | `impute_idw` | Inverse-distance weighting via `gstat` |
| GAIN | *external* | Python GAN-based imputer, run out-of-process (see below) |

Every method plugs into a common `list(fun, transform, scaling)` spec consumed by `run_experiment()` / `impute_with_prep()`, so adding a new method only requires a wrapper of the form `function(mat_na, coords) -> mat_imputed`.

Python-side methods (e.g. GAIN) are not part of this repo. `export_ext_dataset()` in [06_experiment.R](scripts/R/06_experiment.R) writes the pre-processed matrix, mask, coordinates and a `prep.rds` to `output/<run_id>/<method>/r<rep>/p<prop>/`; the external script is expected to write `imputed_r<rep>_p<prop>.csv` back to `output/results/<run_id>/<method>/`, which `import_ext_results()` then reads, reverses the pre-processing on, and scores identically to the R methods.

## Evaluation metrics

Implemented in [05_metrics.R](scripts/R/05_metrics.R), computed only on the masked (originally-missing) entries unless noted:

- **RMSE / NRMSE / MAE** — value-recovery error
- **CCC** — Concordance Correlation Coefficient
- **SAM** — Spectral Angle Mapper
- **VarRatio** — ratio of imputed to true per-feature variance (over/under-smoothing)
- **CorStruct** — preservation of the inter-feature correlation matrix
- **SSIM** — structural similarity of the reconstructed spatial image (requires coordinates)
- **MoranDiff** — absolute difference in Moran's I spatial autocorrelation (requires coordinates)

## Downstream validation

[08_downstream.R](scripts/R/08_downstream.R) runs a fixed module registry (`DOWNSTREAM_MODULES`) on every imputed matrix and compares it against the same analysis run on the ground truth:

- **PCA** — Procrustes SS + variance-explained difference on PC1/PC2
- **Differential expression** — t-test and Wilcoxon, scored by hit-list Jaccard overlap and p-value rank correlation
- **Clustering** — Adjusted Rand Index (vs. ground-truth labels and vs. known groups) + silhouette width
- **Co-abundance network** — Jaccard overlap of thresholded Spearman correlation edges
- **Feature variance / distribution match** — Spearman correlation of per-feature variances; mean KS statistic per feature
- **PLS-DA** — Procrustes SS between imputed and ground-truth discriminant subspaces

[09_downstream_visualisation.R](scripts/R/09_downstream_visualisation.R) provides the corresponding plots (PCA/PLS-DA overlays, Procrustes bars, DE/clustering/network summaries, a pooled intensity-distribution ridge plot).

## Requirements

R (≥ 4.x), managed via [`renv`](https://rstudio.github.io/renv/) (`renv.lock` in this repo — run `renv::restore()` to install exact versions).

Core packages used across the pipeline:

```r
# CRAN
tidyverse, patchwork, jsonlite, readxl, zoo, scales,
missForest, gstat, kernlab, VIM, spdep, vegan, cluster,
ggrepel, ggridges, RColorBrewer, pls

# Bioconductor (install via BiocManager::install())
pcaMethods, MsCoreUtils, SpatialPack, imputeLCMD
```

GAIN (or other external methods) requires a separate Python environment; this repo only handles the R-side export/import contract, not the Python implementation itself.

## Usage

Three independent entry points, depending on the dataset:

### `main.R` — spatial pipeline

For ROI-exported MSI data (a root folder with one `.txt` coordinate subfolder and N `.csv` feature subfolders — see [01_data.R](scripts/R/01_data.R)). Set `DATA`, `MISSING_PROPS`, `SEEDS`, and `missing_model` (`"mcar"`, `"mnar"`, `"hybrid"`, `"msi"`, `"mcar_g"`) at the top of the script, then run top to bottom. Produces `results/benchmark_<dataset>_<model>_<date>_<n>reps.rds` and, after `08_downstream.R`, `results/downstream/downstream_<...>.rds`.

```r
source("main.R")
```

### `main_nospatial.R` — non-spatial pipeline

For a single Excel export with metadata rows (acquisition mode, lipid name, m/z) followed by sample rows (see [01_data_nospatial.R](scripts/R/01_data_nospatial.R)). Same simulate → impute → score → downstream flow, but without spatial metrics (SSIM/Moran's I) or spatial imputation methods.

### `main_mix-missingness.R` — mixed-mechanism zone benchmark

A more targeted analysis on tissue-zone data (e.g. fibrotic vs. normal regions across multiple patients): characterises the *real* missingness structure (per-feature intensity/missingness coupling, fibrosis-marker subsetting), calibrates an intensity-coupled MCAR/MNAR simulator to match it, then benchmarks a global RF imputer against a global QRILC imputer and an "Oracle" that routes each feature to whichever method matches its true mechanism.

All three scripts follow the same pattern: **config block → load data → simulate/benchmark → save `.rds` → reload → plot**. The "reload" and "plot" sections contain multiple `readRDS(...)` calls for different dataset/mechanism combinations — edit these paths to point at your own saved results.

## Data layout

Spatial input (`main.R`), one folder per dataset:

```
data/<dataset_name>/
├── <txt_subfolder>/      # one .txt per ROI, pixel coordinates + scan count
└── <csv_subfolder(s)>/   # one .csv per sample, feature intensities keyed by ROI
```

Non-spatial input (`main_nospatial.R`): a single `.xlsx` with row 1 = acquisition mode (NEG/POS), row 2 = lipid name, row 3 = m/z, rows 4+ = samples (`SampleName` encodes individual + ROI).

`data/`, `output/`, and `results/` are gitignored — datasets and generated artefacts are not tracked in this repository.

## Output

- **`results/*.rds`** — full experiment bundle: config, tidy metrics, imputed-matrix store, coordinates, feature names.
- **`results/downstream/*.rds`** — downstream comparison metrics + ground-truth analysis objects.
- **`results/benchmark/*.rds`** — mixed-mechanism zone-benchmark bundles (results + zones, for reload/plotting).
- **`output/<run_id>/...`** — pre-processed data handed off to external (Python) methods, plus their results once imported back.

## Reproducibility

Each `(replicate, missingness proportion)` combination re-seeds immediately before simulation (`set.seed(seed)`), so results are deterministic given the same `SEEDS` vector. The mixed-mechanism benchmark additionally seeds per `(zone, proportion, replicate)` unit (`cfg$seed_base + rep * 1000 + prop * 1000`) and checkpoints each unit to `results/benchmark/checkpoints/`, so long runs can be resumed (`resume = TRUE`) without recomputation.
