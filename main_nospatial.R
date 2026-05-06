# =============================================================================
# main_nospatial.R
# Purpose:  Runs the full MSI imputation benchmark pipeline without spatial methods
#
# Pipeline:
#   1. Load and assemble data          (01_data_nospatial.R)
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
  "readxl",
  "stringr",
  "patchwork",
  "missForest",
  "pcaMethods",
  "msImpute",
  "gstat",
  "sp",
  "kernlab",
  "SpatialPack",
  "VIM",
  "spdep",
  "jsonlite"
)

new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(readxl)
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
library(jsonlite)


# Source modules ----------------------------------------------------------

source("scripts/R/01_data_nospatial.R")
source("scripts/R/02_simulation.R")
source("scripts/R/03_imputation.R")
source("scripts/R/04_preprocessing.R")
source("scripts/R/05_metrics.R")
source("scripts/R/06_experiment.R")
source("scripts/R/07_visualisation.R")


# Config ------------------------------------------------------------------

MISSING_PROPS <- c(0.1, 0.4)  # missingness proportions to benchmark
# SEEDS <- c(42, 123, 456, 789, 1011)  # one replicate per seed
SEEDS <- c(42)
missing_model <- "mcar"  # one of: "mcar", "mnar", "hybrid", "mcar_g"

# Unique identifier for this pipeline run: <missing_model>_<date>, e.g. "msi_2026-04-22"
RUN_ID <- paste0(missing_model, "_", Sys.Date())

DATA = "data/20220118_POSNEG_combined_all_samples.xlsx"
ACQUISITION <- "BOTH"  # "NEG", "POS" or "BOTH"

# Load the data
df_all  <- load_nospatial_msi(DATA)

# Select metadata (no row filtering needed)
meta <- df_all |> select(Group, Regions, Individual, ROI_num)

# Choose columns based on ACQUISITION
if (ACQUISITION == "BOTH") {
  # Keep all NEG and POS columns
  mat_impute <- df_all |>
    select(starts_with("NEG_") | starts_with("POS_")) |>
    as.matrix()
} else if (ACQUISITION == "NEG") {
  # Keep only NEG columns
  mat_impute <- df_all |>
    select(starts_with("NEG_")) |>
    as.matrix()
} else if (ACQUISITION == "POS") {
  # Keep only POS columns
  mat_impute <- df_all |>
    select(starts_with("POS_")) |>
    as.matrix()
}

message("Matrix dimensions: ", nrow(mat_impute), " rows x ", ncol(mat_impute), " features")

dataset_summary_nospatial(mat_impute, meta)

# Remove features that already have Nas or all-zeros
valid_cols <- colSums(is.na(mat_impute)) == 0 & colSums(mat_impute == 0) == 0
mat_impute <- mat_impute[, valid_cols]

message(
  "Matrix dimensions after removing missing values: ", nrow(mat_impute),
  " rows x ", ncol(mat_impute), " features"
)

# Map string → simulation function
simulation_fn <- list(
  mcar   = function(mat, coords, prop) simulate_MCAR(mat, prop),
  mnar   = function(mat, coords, prop) simulate_MNAR(mat, prop),
  mcar_g = function(mat, coords, prop) simulate_MCAR_global(mat, prop),
  hybrid = function(mat, coords, prop) simulate_hybrid_missing(mat, total_prop = prop),
  msi    = function(mat, coords, prop) simulate_msi_missing(mat, coords, prop)
)[[missing_model]]


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
  ),
  missForest = list(
    fun = impute_missForest,
    transform = "log1p",
    scaling = "none"
  ),
  knn = list(
    fun = impute_knn,
    transform = "log1p",
    scaling = "pareto"
  ),
  qrilc = list(
    fun = impute_qrilc,
    transform = "none",
    scaling = "none"
  ),
  ppca = list(
    fun = impute_ppca,
    transform = "log1p",
    scaling = "pareto"
  ),
  bpca = list(
    fun = impute_bpca,
    transform = "log1p",
    scaling = "pareto"
  ),
  svd = list(
    fun = impute_svd,
    transform = "log1p",
    scaling = "pareto"
  ),
  nngp = list(
    fun = impute_nngp,
    transform = "none",
    scaling = "none"
  )
)

external_methods <- list(
  GAIN = list(
    transform = "log1p",
    scaling   = "range"
  )# ,
  # ConvGAIN = list(
  #   transform = "log1p",
  #   scaling   = "range"
  # ),
  # SpatialGAIN = list(
  #   transform = "log1p",
  #   scaling   = "range"
  # ),
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

# results_storage: stores na_matrix + imputed matrices per replicate × prop
# Indexed as results_storage[[ rep_key ]][[ p_str ]]
# where rep_key = "r1", "r2", ... (one per seed)
results_storage <- list()

results <- map_dfr(seq_along(SEEDS), function(rep_idx) {
  seed    <- SEEDS[[rep_idx]]
  rep_key <- paste0("r", rep_idx)  # e.g. "r1", "r2", ...

  message(
    "\n══ Replicate ", rep_idx, " / ", length(SEEDS),
    "  (seed = ", seed, ") ══════════════════════════"
  )

  results_storage[[rep_key]] <<- list()

  map_dfr(MISSING_PROPS, function(p) {
    message("  ── Missing proportion: ", p, " ──────────────────────────")

    # 1. Seed is fixed immediately before each simulation so each (rep, prop)
    #    combination is fully reproducible and independent.
    set.seed(seed)
    mat_na <- simulation_fn(mat_impute, coords, prop = p)  # simulate missing values

    p_str <- as.character(p)

    # 2. Export for Python methods (06_experiment)
    walk(names(external_methods), function(method_name) {
      export_ext_dataset(
        mat_true    = mat_impute,
        mat_na      = mat_na,
        coords      = coords,
        method_name = method_name,
        preprocessing = external_methods[[method_name]],
        prop        = p,
        ms_model    = missing_model,
        rep_idx     = rep_idx,
        run_id      = RUN_ID
      )
    })

    # 3. Run R imputation methods (06_experiment)
    exp_output <- run_experiment(
      mat_true = mat_impute,
      mat_na   = mat_na,
      coords   = coords,
      methods  = imputation_methods,
      metrics  = metric_functions
    )

    # 4. Store matrices for visualisation
    results_storage[[rep_key]][[p_str]] <<- list(
      na_matrix = mat_na,
      imputed   = exp_output$matrices
    )

    exp_output$metrics |>
      mutate(prop = p, seed = seed, rep = rep_idx)
  })
})

results

# 5. Import and combine external (Python) results
external_results <- import_ext_results(
  mat_true         = mat_impute,
  coords           = coords,
  metric_functions = metric_functions,
  run_id           = RUN_ID
)

# 6. Merge metrics from R and external methods into a single data frame
all_results <- bind_rows(results, external_results$metrics)

all_results

# 7. Merges matrices into results_storage
walk(names(external_results$matrices), function(rep_key) {
  walk(names(external_results$matrices[[rep_key]]), function(p_str) {
    results_storage[[rep_key]][[p_str]]$imputed <<- c(
      results_storage[[rep_key]][[p_str]]$imputed,
      external_results$matrices[[rep_key]][[p_str]]
    )
  })
})


# Save --------------------------------------------------------------------

experiment_output <- list(
  config = list(
    dataset       = DATA,
    missing_model = missing_model,
    missing_props = MISSING_PROPS,
    seeds         = SEEDS,
    n_replicates  = length(SEEDS),
    run_id        = RUN_ID,
    date          = Sys.time()
  ),
  results         = results,
  all_results     = all_results,
  results_storage = results_storage,
  feature_names   = colnames(mat_impute)
)

dir.create("results", showWarnings = FALSE)

# For future reproducibility and reference
saveRDS(
  experiment_output,
  file = file.path(
    "results",
    paste0(
      "benchmark_", basename(DATA), "_", RUN_ID,
      "_", length(SEEDS), "reps.rds"
    )
  )
)


# Load a previous experiment ----------------------------------------------

list.files(path = "results", pattern = "\\.rds")

exp <- readRDS("results/benchmark_20220118_POSNEG_combined_all_samples.xlsx_mcar_2026-05-06_1reps.rds")

all_results     <- exp$all_results
results_storage <- exp$results_storage
feature_names   <- exp$feature_names

plot_metric(exp$all_results, "NRMSE", type = "line")  # example


# Plots -------------------------------------------------------------------

source("scripts/R/07_visualisation.R")

# Aggregated metric plots (mean ± ribbon across replicates)
plot_metric(all_results, "NRMSE")
plot_metric(all_results, "MAE")
plot_metric(exp$all_results, "CCC", type='bar')
plot_metric(all_results, "VarRatio",   y_label = "Mean Variance Ratio")
plot_metric(exp$all_results, "CorStruct",  y_label = "Correlation of Correlation Matrices")

plot_runtime_vs_nrmse(exp$all_results)
plot_spectral_preservation(exp$all_results)