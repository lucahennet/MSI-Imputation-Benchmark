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
#   7. Run downstream analyses         (08_downstream.R). 
# =============================================================================


# Dependencies ------------------------------------------------------------

# Core Environment Setup
library(tidyverse)
library(patchwork)
library(jsonlite)

# Pipeline-specific
library(missForest) # RF
library(gstat) # IDW
library(kernlab) # GP
library(SpatialPack) # SSIM
library(VIM) # KNN
library(spdep) # Moran’s I
library(vegan) # Procrustes
library(cluster) # silhouette

# Pipeline-specific (BiocManager)
library(pcaMethods) # PPCA, BPCA, SVM
library(MsCoreUtils) # QRILC


# Source modules ----------------------------------------------------------

source("scripts/R/01_data.R")
source("scripts/R/02_simulation.R")
source("scripts/R/03_imputation.R")
source("scripts/R/04_preprocessing.R")
source("scripts/R/05_metrics.R")
source("scripts/R/06_experiment.R")


# Config ------------------------------------------------------------------

MISSING_PROPS <- c(0.1, 0.4)  # missingness proportions to benchmark
SEEDS <- 42
# SEEDS <- c(42, 123, 456, 789, 1011)  # one replicate per seed
missing_model <- "msi"  # one of: "mcar", "mnar", "hybrid", "msi", "mcar_g"

# DATA = "data/MvaExport_10-04-2026_10.23.29.798"
DATA = "data/MvaExport_10-04-2026_10.23.29.798_TIC"
DATA = "data/ROI_Luca1303"

# DATA = "data/MvaExport_30-06-2026_11.36.12.551_fib"
# DATA = "data/MvaExport_30-06-2026_12.04.46.174_norm"

# Unique identifier for this pipeline run: <missing_model>_<date>, e.g. "msi_2026-04-22"
RUN_ID <- paste0(missing_model, "_", Sys.Date())

# Load the data
raw_wide <- assemble_msi_data(DATA) |>
  distinct(X, Y, .keep_all = TRUE)  # remove duplicate pixels

FEATURE_RANGE <- 8:ncol(raw_wide)  # which features to load (min = 8, max = 1007)
# FEATURE_RANGE <- 8:10

dataset_summary(raw_wide, FEATURE_RANGE)

feature_cols <- colnames(raw_wide)[FEATURE_RANGE]

# Long-format, normalised
# df_visualisation <- raw_wide |>
#   pivot_msi_long(feature_indices = 1995:1997) |>
#   normalize_msi_data()
# generate_msi_plots(df_visualisation, ncol = 3)
# -> look where if I still use it or if I can remove it

# Wide matrix ready for imputation
df_impute <- raw_wide |>
  select(X, Y, RAW.TIC.OR.ROI.sum.peak, all_of(feature_cols)) |>
  select(
    X, Y, RAW.TIC.OR.ROI.sum.peak,
    where(~ !any(is.na(.)) && !any(. == 0)) # drop all-NA or all-zero features
  ) |>
  # - facultative part -
  # mutate(across(
  #   -c(X, Y, RAW.TIC.OR.ROI.sum.peak),
  #   ~ . / RAW.TIC.OR.ROI.sum.peak # TIC normalisation !depending on the dataset!
  # )) |>
  # # - end -
  select(-RAW.TIC.OR.ROI.sum.peak)

# Coordinates and numeric matrix
coords     <- df_impute |> select(X, Y)
mat_impute <- df_impute |> select(-X, -Y) |> as.matrix()

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
  # Zero = list(
  #   fun = impute_zero,
  #   transform = "none",
  #   scaling = "none"
  # ),
  # Mean = list(
  #   fun = impute_mean,
  #   transform = "none",
  #   scaling = "none"
  # ),
  Median = list(
    fun = impute_median,
    transform = "none",
    scaling = "none"
  ),
  HM = list(
    fun = impute_half_min,
    transform = "none",
    scaling = "none"
  )# ,
  # RF = list(
  #   fun = impute_missForest,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  # KNN = list(
  #   fun = impute_knn,
  #   transform = "log1p",
  #   scaling = "pareto"
  # ),
  # QRILC = list(
  #   fun = impute_qrilc,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  # PPCA = list(
  #   fun = impute_ppca,
  #   transform = "log1p",
  #   scaling = "pareto"
  # ),
  # BPCA = list(
  #   fun = impute_bpca,
  #   transform = "log1p",
  #   scaling = "pareto"
  # ),
  # SVD = list(
  #   fun = impute_svd,
  #   transform = "log1p",
  #   scaling = "pareto"
  # ),
  # NNGP = list(
  #   fun = impute_nngp,
  #   transform = "none",
  #   scaling = "none"
  # ),
  # spKNN = list(
  #   fun = impute_spatial_knn,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  # GP = list(
  #   fun = impute_gp,
  #   transform = "log1p",
  #   scaling = "none"
  # ),
  # IDW = list(
  #   fun = impute_idw,
  #   transform = "log1p",
  #   scaling = "none"
  # )
)

external_methods <- list(
  GAIN = list(
    transform = "log1p",
    scaling   = "range"
  )
)


# Run pipeline ------------------------------------------------------------

# results_storage: stores na_matrix + imputed matrices per replicate × prop
# Indexed as results_storage[[ rep_key ]][[ p_str ]]
# where rep_key = "r1", "r2", ... (one per seed)
results_storage <- list()

results <- map_dfr(seq_along(SEEDS), function(rep_idx) {
  seed <- SEEDS[[rep_idx]]
  rep_key <- paste0("r", rep_idx) # e.g. "r1", "r2", ...

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
    mat_na <- simulation_fn(mat_impute, coords, prop = p) # simulate missing values

    p_str <- as.character(p)

    # 2. Export for Python methods (06_experiment)
    walk(names(external_methods), function(method_name) {
      export_ext_dataset(
        mat_true = mat_impute,
        mat_na = mat_na,
        coords = coords,
        method_name = method_name,
        preprocessing = external_methods[[method_name]],
        prop = p,
        ms_model = missing_model,
        rep_idx = rep_idx,
        run_id = RUN_ID
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
  coords          = coords,
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

exp <- readRDS("results/benchmark_MvaExport_10-04-2026_10.23.29.798_TIC_mcar_2026-06-17_50reps.rds")
exp <- readRDS("results/benchmark_MvaExport_10-04-2026_10.23.29.798_TIC_mnar_2026-06-17_50reps.rds")
exp <- readRDS("results/benchmark_MvaExport_10-04-2026_10.23.29.798_TIC_hybrid_2026-06-17_50reps.rds")
exp <- readRDS("results/benchmark_MvaExport_10-04-2026_10.23.29.798_TIC_msi_2026-06-17_50reps.rds")

exp <- readRDS("results/benchmark_ROI_Luca1303_mcar_2026-06-16_50reps.rds")
exp <- readRDS("results/benchmark_ROI_Luca1303_mnar_2026-06-16_50reps.rds")
exp <- readRDS("results/benchmark_ROI_Luca1303_hybrid_2026-06-16_50reps.rds")
exp <- readRDS("results/benchmark_ROI_Luca1303_msi_2026-06-16_50reps.rds")

all_results     <- exp$all_results
results_storage <- exp$results_storage
coords          <- exp$coords
feature_names   <- exp$feature_names


# Plots -------------------------------------------------------------------

source("scripts/R/07_visualisation.R")

# Aggregated metric plots (mean ± ribbon across replicates)
plot_metric(all_results, "NRMSE")
plot_metric(all_results, "MAE")
plot_metric(all_results, "SAM")
plot_metric(exp$all_results, "CCC", type='bar')
plot_metric(exp$all_results, "SSIM", type='line')
plot_metric(exp$all_results, "MoranDiff",  y_label = "Mean Abs. Diff. in Moran's I")
plot_metric(all_results, "VarRatio",   y_label = "Mean Variance Ratio")
plot_metric(exp$all_results, "CorStruct",  y_label = "Correlation of Correlation Matrices")

plot_runtime_vs_nrmse(exp$all_results)
plot_spatial_fidelity(exp$all_results)
plot_spectral_preservation(exp$all_results)

visualise_heatmaps(feature_idx = 100, mode = "all_methods", target_prop = 0.1,     ncol = 4, rep_idx = 1)
visualise_heatmaps(feature_idx = 9,  mode = "all_props",   target_method = "PPCA", ncol = 3, rep_idx = 1)


# Downstream analyses -------------------------------------------------------

source("scripts/R/08_downstream.R")
source("scripts/R/09_downstream_visualisation.R")

# -- Config -----------------------------------------------------------------
# The spatial pipeline has no group/region labels, so DE is skipped.
# DS_CONTRAST <- NULL
# DS_LEVELS   <- NULL

# -- Ground truth -----------------------------------------------------------
message("\n══ Ground truth downstream analyses ══════════════════════════")
gt_downstream <- run_all_downstream(
  mat          = mat_impute,
  meta         = coords,         # coords (X, Y) passed as meta for PCA colouring
  contrast_col = "Group",        # will be silently skipped: "Group" not in coords
  de_levels    = NULL,
  cor_threshold = 0.7
)

# -- Comparison across all methods × props × replicates --------------------
message("\n══ Downstream comparisons ═════════════════════════════════════")
downstream_results <- compute_downstream_comparison(
  results_storage = results_storage,
  meta            = coords,
  ground_truth    = gt_downstream,
  contrast_col    = "Group",
  de_levels       = NULL
)

downstream_results$metrics  # inspect tidy table

# # -- Save -------------------------------------------------------------------
# dir.create("results/downstream", showWarnings = FALSE)
# 
# # downstream_results <- downstream_results
# saveRDS(
#   downstream_results,
#   file = file.path(
#     "results/downstream",
#     paste0(
#       "downstream_", basename(DATA), "_", RUN_ID,
#       "_", length(SEEDS), "reps.rds"
#     )
#   )
# )
# 
# # -- Load a previous experiment ---------------------------------------------
# list.files(path = "results/downstream", pattern = "\\.rds")
# 
# downstream_results <- readRDS("results/downstream/benchmark_ROI_Luca1303_mnar_2026-06-16_50reps_downstream.rds")
# downstream_results <- readRDS("results/downstream/benchmark_ROI_Luca1303_hybrid_2026-06-16_50reps_downstream.rds")

# -- Save Downstream Results --------------------------------------------------
dir.create("results/downstream", showWarnings = FALSE)

downstream_output <- list(
  metrics = downstream_results,
  gt      = gt_downstream
)

saveRDS(
  downstream_output,
  file = file.path(
    "results/downstream",
    paste0("downstream_", basename(DATA), "_", RUN_ID, "_", length(SEEDS), "reps.rds")
  )
)

# -- Load a previous experiment (Compatible with both local and cluster runs) -
list.files(path = "results/downstream", pattern = "\\.rds")

# Reading standardized uniform archive container
downstream_data <- readRDS("results/downstream/benchmark_MvaExport_10-04-2026_10.23.29.798_TIC_mcar_2026-06-17_50reps_downstream.rds")
downstream_data <- readRDS("results/downstream/benchmark_MvaExport_10-04-2026_10.23.29.798_TIC_mnar_2026-06-17_50reps_downstream.rds")
downstream_data <- readRDS("results/downstream/benchmark_MvaExport_10-04-2026_10.23.29.798_TIC_hybrid_2026-06-17_50reps_downstream.rds")
downstream_data <- readRDS("results/downstream/benchmark_MvaExport_10-04-2026_10.23.29.798_TIC_msi_2026-06-17_50reps_downstream.rds")

downstream_data <- readRDS("results/downstream/benchmark_ROI_Luca1303_mcar_2026-06-16_50reps_downstream.rds")
downstream_data <- readRDS("results/downstream/benchmark_ROI_Luca1303_mnar_2026-06-16_50reps_downstream.rds")
downstream_data <- readRDS("results/downstream/benchmark_ROI_Luca1303_hybrid_2026-06-16_50reps_downstream.rds")
downstream_data <- readRDS("results/downstream/benchmark_ROI_Luca1303_msi_2026-06-16_50reps_downstream.rds")

# Expose both environments directly to local visualization scripts seamlessly
downstream_results <- downstream_data$metrics
gt_downstream      <- downstream_data$gt

# -- Plots ------------------------------------------------------------------

# PCA
plot_pca_overlay(downstream_results, gt_downstream, target_prop = 0.1, rep_idx = 1)
plot_pca_overlay(downstream_results, gt_downstream, target_prop = 0.4, rep_idx = 1)
plot_procrustes(downstream_results)

# Differential
plot_de_comparisons(downstream_results)

# Clustering
plot_clustering_metrics(downstream_results)

# Network
plot_network_metrics(downstream_results)

# All-in-one summary
plot_downstream_summary(downstream_results, target_prop = 0.1)
plot_downstream_summary(downstream_results, target_prop = 0.4)

# Individual metrics via the existing plot_metric() from 07_visualisation.R
plot_metric(downstream_results$metrics, "ProcrustesSS",  y_label = "Procrustes SS (lower = better)")
plot_metric(downstream_results$metrics, "EdgeJaccard",   y_label = "Edge Jaccard (higher = better)")
plot_metric(downstream_results$metrics, "FeatureVarCor", y_label = "Per-feature Variance Correlation (Higher = Better)")
# plot_metric(downstream_results$metrics, "NMF_SpatialCor", y_label = "NMF Spatial Map Preservation (Higher = Better)")