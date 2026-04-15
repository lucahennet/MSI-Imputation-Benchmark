# =============================================================================
# main.R
# Purpose:  Runs the full MSI imputation benchmark pipeline
#
# Pipeline:
#   1. Load and assemble data          (01_data.R)
#   2. Simulate missing values         (02_simulation.R)
#   3. Export datasets for Python      (06_experiment.R)
#   4. Run R imputation methods        (06_experiment.R)
#   5. Import external results         (06_experiment.R)
#   6. Combine and plot metrics        (07_visualisation.R)
# =============================================================================


# Dependencies ------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(version = "3.22", update = FALSE)

BiocManager::install(c(
  "imputeLCMD",
  "pcaMethods",
  "msImpute"
), update = FALSE)

packages <- c(
  "tidyverse",
  "missForest",
  "pcaMethods",
  "msImpute",
  "gstat",
  "sp",
  "kernlab",
  "SpatialPack",
  "VIM",
  "spdep"
)

new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(stringr)
library(patchwork)
library(missForest) # RF
library(imputeLCMD) # QRILC
library(pcaMethods) # PPCA, BPCA, SVM
library(msImpute) # NNGP
library(gstat) # IDW
library(sp) # ?
library(kernlab) # GP
library(SpatialPack)
library(VIM) # kNN
library(spdep) # Moran’s I


# Source modules ----------------------------------------------------------

source("scripts/R/01_data.R")
source("scripts/R/02_simulation.R")
source("scripts/R/03_imputation.R")
source("scripts/R/04_preprocessing.R")
source("scripts/R/05_metrics.R")
source("scripts/R/06_experiment.R")
source("scripts/R/07_visualisation.R")


# Config ------------------------------------------------------------------

MISSING_PROPS   <- c(0.1, 0.4)        # missingness levels to benchmark
RANDOM_SEED     <- 42
missing_model   <- "msi"              # one of: "mcar", "mnar", "hybrid", "msi"

DATA = "data/MvaExport_10-04-2026_10.23.29.798"
# DATA = "data/ROI_Luca1303"

# Load the data
raw_wide <- assemble_msi_data(DATA) |>
  distinct(X, Y, .keep_all = TRUE)   # remove duplicate pixels

FEATURE_RANGE   <- 8:ncol(raw_wide)   # which features to load (min = 8, max = 1007)
# FEATURE_RANGE   <- 8:10

# TODO eventually add stats about the structure, missing values in the data


feature_cols <- colnames(raw_wide)[FEATURE_RANGE]

# Long-format, normalised
# df_visualisation <- raw_wide |>
#   pivot_msi_long(feature_indices = FEATURE_RANGE) |>
#   normalize_msi_data()
# generate_msi_plots(df_visualisation, ncol = 3)
# -> look where if I still use it or if I can remove it

# Wide matrix ready for imputation
df_impute <- raw_wide |>
  select(X, Y, RAW.TIC.OR.ROI.sum.peak, all_of(feature_cols)) |>
  select(
    X, Y, RAW.TIC.OR.ROI.sum.peak,
    where(~ !any(is.na(.)) && !any(. == 0))   # drop all-NA or all-zero features
  ) |>
  mutate(across(
    -c(X, Y, RAW.TIC.OR.ROI.sum.peak),
    ~ . / RAW.TIC.OR.ROI.sum.peak             # TIC normalisation
  )) |>
  select(-RAW.TIC.OR.ROI.sum.peak)

# Coordinates and numeric matrix
coords     <- df_impute |> select(X, Y)
mat_impute <- df_impute |> select(-X, -Y) |> as.matrix()

# Map string → simulation function
simulation_fn <- list(
  mcar   = simulate_MCAR,
  mnar   = simulate_MNAR,
  hybrid = simulate_hybrid_missing,
  msi    = simulate_msi_missing
)[[MISSING_MODEL]]

# Methods -----------------------------------------------------------------

imputation_methods <- list(
  
  # zero = list(
  #   fun = impute_zero,
  #   transform = "none",
  #   scaling = "none"
  # ),
  
  mean = list(
    fun = impute_mean,
    transform = "none",
    scaling = "none"
  ),
  
  median = list(
    fun = impute_median,
    transform = "none",
    scaling = "none"
  ),
  
  half_min = list(
    fun = impute_half_min,
    transform = "none",
    scaling = "none"
  )# ,
  
  # missForest = list(
  #   fun = impute_missForest,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  
  # knn = list(
  #   fun = impute_knn,
  #   transform = "log1p",
  #   scaling = "pareto"
  # ),
  
  # qrilc = list(
  #   fun = impute_qrilc,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  
  # ppca = list(
  #   fun = impute_ppca,
  #   transform = "log1p",
  #   scaling = "pareto"
  # ),
  
  # bpca = list(
  #   fun = impute_bpca,
  #   transform = "log1p",
  #   scaling = "pareto"
  # ),
  
  # svd = list(
  #   fun = impute_svd,
  #   transform = "log1p",
  #   scaling = "pareto"
  # ),
  
  # nngp = list(
  #   fun = impute_nngp,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  
  # sp_knn = list(
  #   fun = impute_spatial_knn,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  
  # gp = list(
  #   fun = impute_gp,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  
  # idw = list(
  #   fun = impute_idw,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  
)

external_methods <- list(
  
  GAIN = list(
    transform = "log1p",
    scaling   = "range"
  )# ,
  # 
  # UNET = list(
  #   transform = "log1p",
  #   scaling   = "zscore"
  # ),
  # 
  # AUTOENCODER = list(
  #   transform = "log1p",
  #   scaling   = "zscore"
  # ),
  # 
  # DIFFUSION = list(
  #   transform = "log1p",
  #   scaling   = "zscore"
  # ),
  # 
  # KRIGING = list(
  #   transform = "log1p",
  #   scaling   = "none"
  # )
)

# Run pipeline ------------------------------------------------------------

set.seed(RANDOM_SEED)

results_storage <- list()   # stores na_matrix + imputed matrices per prop

results <- map_dfr(MISSING_PROPS, function(p) {
  message("\n── Missing proportion: ", p, " ──────────────────────────")
  
  # 1. Simulate missing values
  mat_na <- simulation_fn(mat_impute, coords, prop = p)
  
  # 2. Export for Python methods
  walk(names(external_methods), function(method_name) {
    export_ext_dataset(
      mat_true     = mat_impute,
      mat_na       = mat_na,
      coords       = coords,
      method_name  = method_name,
      preprocessing = external_methods[[method_name]],
      prop         = p
    )
  })
  
  # 3. Run R imputation methods
  exp_output <- run_experiment(
    mat_true = mat_impute,
    mat_na   = mat_na,
    coords   = coords,
    methods  = imputation_methods,
    metrics  = metric_functions
  )
  
  # 4. Store matrices for visualisation
  results_storage[[as.character(p)]] <<- list(
    na_matrix = mat_na,
    imputed   = exp_output$matrices
  )
  
  exp_output$metrics |> mutate(prop = p)
})

results

# 5. Import and combine external (Python) results
external_results <- import_ext_results(
  mat_true         = mat_impute,
  coords           = coords,
  metric_functions = metric_functions
)

all_results <- bind_rows(results, external_results)


# Plots -------------------------------------------------------------------

source("scripts/R/07_visualisation.R")   # re-source if iterating on plots only

plot_metric(all_results, "NRMSE")
plot_metric(all_results, "MAE")
plot_metric(all_results, "CCC")
plot_metric(all_results, "SSIM")
plot_metric(all_results, "MoranDiff",  y_label = "Mean Abs. Diff. in Moran's I")
plot_metric(all_results, "VarRatio",   y_label = "Mean Variance Ratio")
plot_metric(all_results, "CorStruct",  y_label = "Correlation of Correlation Matrices")

plot_runtime_vs_nrmse(all_results)
plot_spatial_fidelity(all_results)
plot_spectral_preservation(all_results)

visualise_heatmaps(feature_idx = 10,  mode = "all_methods", target_prop = 0.1, ncol = 4)
visualise_heatmaps(feature_idx = 500, mode = "all_props",   target_method = "mean", ncol = 3)
