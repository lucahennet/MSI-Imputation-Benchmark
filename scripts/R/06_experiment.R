# =============================================================================
# 06_experiment.R
# Purpose:  Run imputation methods, export datasets for external (Python)
#           methods, and import their results back into R.
# Inputs:   mat_true, mat_na, coords, method specs, metric functions
# Outputs:  metrics data frame, imputed matrices, exported CSV/RDS files
# Depends:  04_preprocessing.R, 05_metrics.R, tidyverse, stringr, jsonlite
# =============================================================================


# Functions ---------------------------------------------------------------

#' Run the full imputation and evaluation pipeline for a given dataset and set of methods.
#' 
#' For each method, applies the specified preprocessing, runs the imputation, 
#' reverses the preprocessing, and computes metrics. Also collects imputed matrices 
#' for later use (e.g. export).
#' @param mat_true The complete matrix without missing values (for metric computation)
#' @param mat_na   The matrix with missing values (NA) to be imputed
#' @param coords   Data frame of pixel coordinates (for spatial methods)
#' @param methods  A named list of method specifications, where each entry is a list with:
#'                  - fun: the imputation function to call
#'                  - transform: the name of the transformation to apply (e.g. "log1p")
#'                  - scaling: the name of the scaling to apply (e.g. "pareto")
#' @param metrics A named list of metric functions to compute (e.g. RMSE, MAE)
#' @return A list with:
#'  Returns a list with:
#'    $metrics  - tidy data frame with method, runtime, metric values, and
#'    $matrices - named list of imputed matrices (method name as key)
run_experiment <- function(mat_true, mat_na, coords, methods, metrics) {
  mask <- is.na(mat_na)
  imputed_matrices <- list()

  method_results <- map_dfr(names(methods), function(name) {
    cat("Running:", name, "\n")
    method_spec <- methods[[name]]
    prep <- apply_preprocessing(mat_na, method_spec)

    start_time <- Sys.time()
    mat_imp_raw <- tryCatch(
      {
        method_spec$fun(prep$data, coords)
      },
      error = function(e) {
        warning(paste("Error in method:", name, "|", e$message))
        return(NULL)
      }
    )
    runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    if (is.null(mat_imp_raw)) {
      return(tibble(method = name, success = FALSE))
    }

    mat_imp <- reverse_preprocessing(mat_imp_raw, prep)

    imputed_matrices[[name]] <<- mat_imp

    metric_vals <- compute_metrics(mat_true, mat_imp, mask, metrics)
    extra_vals <- extra_metrics(mat_true, mat_imp, coords, mask)

    tibble(
      method = name, success = TRUE, runtime_sec = runtime,
      !!!as.list(metric_vals), !!!extra_vals
    )
  })

  return(list(metrics = method_results, matrices = imputed_matrices))
}


#' Export a pre-processed dataset for one external method + one missingness
#' proportion + one replicate.
#'
#' Output folder structure:
#'   <root_dir>/<run_id>/<method_name>/r<rep_idx>/<p_str>/
#'     data.csv      – pre-processed matrix with NAs filled as 0
#'     mask.csv      – TRUE where observed, FALSE where missing
#'     coords.csv    – pixel coordinates
#'     prep.rds      – preprocessing object for reverse-transform in R
#'     metadata.json – run provenance
#'
#' @param run_id  String identifying this pipeline run, e.g. "msi_2026-04-22".
#'   Created once in main.R as: paste0(missing_model, "_", Sys.Date())
#' @param rep_idx Integer replicate index (1, 2, ...).
export_ext_dataset <- function(
  mat_true,
  mat_na,
  coords,
  method_name,
  preprocessing,
  prop,
  ms_model = "unknown",
  rep_idx = 1,
  run_id = paste0(ms_model, "_", Sys.Date()),
  root_dir = "output"
) {
  p_str <- paste0("p", round(prop * 100))
  rep_key <- paste0("r", rep_idx)

  # <root_dir>/<run_id>/<method_name>/r<N>/<p_str>/
  out_dir <- file.path(root_dir, run_id, method_name, rep_key, p_str)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  mask <- !is.na(mat_na)
  prep <- apply_preprocessing(mat_na, preprocessing)
  data_proc <- prep$data
  data_proc[is.na(data_proc)] <- 0

  write.csv(data_proc, file.path(out_dir, "data.csv"), row.names = FALSE)
  write.csv(mask, file.path(out_dir, "mask.csv"), row.names = FALSE)
  write.csv(coords, file.path(out_dir, "coords.csv"), row.names = FALSE)
  saveRDS(prep, file.path(out_dir, "prep.rds"))

  metadata <- list(
    method             = method_name,
    missing_proportion = prop,
    missing_model      = ms_model,
    replicate          = rep_idx,
    run_id             = run_id,
    transform          = preprocessing$transform,
    scaling            = preprocessing$scaling,
    n_pixels           = nrow(mat_na),
    n_features         = ncol(mat_na),
    export_time        = as.character(Sys.time())
  )

  writeLines(
    toJSON(metadata, pretty = TRUE, auto_unbox = TRUE),
    file.path(out_dir, "metadata.json")
  )

  message("Exported: ", run_id, " | ", method_name, " | ", rep_key, " | ", p_str)
}


#' Import imputed results produced by external (Python) methods.
#'
#' Expects Python to write output files as:
#'   <results_root>/<run_id>/<method_name>/imputed_r<rep_idx>_p<prop>.csv
#'
#' The matching IO folder (prep.rds, mask.csv) is read from:
#'   <io_root>/<run_id>/<method_name>/r<rep_idx>/<p_str>/
#'
#' @param run_id  The same run_id used during export (e.g. "msi_2026-04-22").
#'   If NULL, the function uses the most recently created sub-folder of
#'   <results_root>, so interactive use requires no extra argument.
#'
#' @return List with:
#'   $metrics  - tidy data frame: method, prop, rep + all metric columns
#'   $matrices - nested list [[ rep_key ]][[ p_str ]][[ method_name ]]
import_ext_results <- function(
  mat_true,
  coords,
  metric_functions,
  run_id = NULL,
  io_root = "output",
  results_root = file.path(io_root, "results")
) {
  empty <- list(metrics = tibble(), matrices = list())

  # Resolve run_id: default to the most recent sub-folder
  if (is.null(run_id)) {
    if (!dir.exists(results_root)) {
      message(
        "import_ext_results: results folder does not exist yet (", results_root,
        ") — skipping external results."
      )
      return(empty)
    }
    run_dirs <- list.dirs(results_root, recursive = FALSE, full.names = FALSE)
    if (length(run_dirs) == 0) {
      message(
        "import_ext_results: no run sub-folders found in ", results_root,
        " — skipping external results."
      )
      return(empty)
    }
    run_id <- tail(sort(run_dirs), 1) # alphabetical sort; ISO dates sort correctly
    message("import_ext_results: using run_id = '", run_id, "'")
  }

  run_results_dir <- file.path(results_root, run_id)
  run_io_dir <- file.path(io_root, run_id)

  if (!dir.exists(run_results_dir)) {
    message(
      "import_ext_results: results folder not found (", run_results_dir,
      ") — no Python results to import yet."
    )
    return(empty)
  }
  if (!dir.exists(run_io_dir)) {
    message(
      "import_ext_results: IO folder not found (", run_io_dir,
      ") — skipping external results."
    )
    return(empty)
  }

  method_dirs <- list.dirs(run_results_dir, recursive = FALSE, full.names = TRUE)

  if (length(method_dirs) == 0) {
    message("No method sub-folders found in: ", run_results_dir)
    return(list(metrics = tibble(), matrices = list()))
  }

  # ---- per-method ---------------------------------------------------------
  all_entries <- map(method_dirs, function(method_path) {
    method_name <- basename(method_path)

    imp_files <- list.files(
      method_path,
      pattern    = "^imputed_r\\d+_p\\d+\\.csv$",
      full.names = TRUE
    )

    if (length(imp_files) == 0) {
      message("  No result files for: ", method_name)
      return(NULL)
    }

    # ---- per-file ---------------------------------------------------------
    map(imp_files, function(file) {
      fname <- basename(file)
      rep_idx <- as.integer(str_extract(fname, "(?<=r)\\d+"))
      p_pct <- as.integer(str_extract(fname, "(?<=p)\\d+"))
      prop <- p_pct / 100
      rep_key <- paste0("r", rep_idx)
      p_str <- paste0("p", p_pct)

      p_dir <- file.path(run_io_dir, method_name, rep_key, p_str)

      if (!dir.exists(p_dir)) {
        message("  Missing IO folder: ", p_dir)
        return(NULL)
      }

      prep_path <- file.path(p_dir, "prep.rds")
      mask_path <- file.path(p_dir, "mask.csv")

      if (!file.exists(prep_path) || !file.exists(mask_path)) {
        message("  Missing prep.rds or mask.csv in: ", p_dir)
        return(NULL)
      }

      message("  Importing: ", method_name, " | ", rep_key, " | ", p_str)

      prep <- readRDS(prep_path)
      mask <- !as.matrix(read.csv(mask_path))
      mat_imp <- reverse_preprocessing(as.matrix(read.csv(file)), prep)

      metric_vals <- compute_metrics(mat_true, mat_imp, mask, metric_functions)
      extra_vals <- extra_metrics(mat_true, mat_imp, coords, mask)

      list(
        metric_row = tibble(
          method = method_name, prop = prop, rep = rep_idx,
          !!!as.list(metric_vals), !!!extra_vals
        ),
        rep_key = rep_key,
        p_str = p_str,
        mat = mat_imp
      )
    })
  }) |>
    unlist(recursive = FALSE) |> # flatten one level: list of per-file results
    Filter(Negate(is.null), x = _) # drop skipped files

  if (length(all_entries) == 0) {
    message("No results successfully imported.")
    return(list(metrics = tibble(), matrices = list()))
  }

  # ---- assemble outputs ---------------------------------------------------

  metrics_df <- map_dfr(all_entries, ~ .x$metric_row)

  # Build nested matrices list [[ rep_key ]][[ p_str ]][[ method_name ]]
  all_matrices <- list()
  walk(all_entries, function(e) {
    rk <- e$rep_key
    ps <- e$p_str
    mn <- e$metric_row$method

    if (is.null(all_matrices[[rk]])) all_matrices[[rk]] <<- list()
    if (is.null(all_matrices[[rk]][[ps]])) all_matrices[[rk]][[ps]] <<- list()
    all_matrices[[rk]][[ps]][[mn]] <<- e$mat
  })

  list(metrics = metrics_df, matrices = all_matrices)
}