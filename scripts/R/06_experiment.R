# =============================================================================
# 06_experiment.R
# Purpose:  Load, assemble and reshape MSI data from raw ROI files
# Inputs:   - A root directory containing one TXT subfolder and N CSV subfolders
# Outputs:  - raw_wide       : wide-format data frame (pixels × features)
#           - df_visualisation: long-format, TIC-normalised, for plotting
#           - df_impute      : clean wide matrix ready for imputation
#           - coords         : data frame of X/Y pixel coordinates
#           - mat_impute     : numeric matrix (pixels × features), no metadata
# Depends:  tidyverse, stringr
# =============================================================================


# Functions ---------------------------------------------------------------

run_experiment <- function(mat_true, mat_na, coords, methods, metrics) {
  mask <- is.na(mat_na)
  imputed_matrices <- list()
  
  method_results <- map_dfr(names(methods), function(name) {
    cat("Running:", name, "\n")
    method_spec <- methods[[name]]
    prep <- apply_preprocessing(mat_na, method_spec)
    
    start_time <- Sys.time()
    mat_imp_raw <- tryCatch({
      method_spec$fun(prep$data, coords)
    }, error = function(e) {
      warning(paste("Error in method:", name, "|", e$message))
      return(NULL)
    })
    runtime <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
    
    if (is.null(mat_imp_raw)) return(tibble(method = name, success = FALSE))
    
    mat_imp <- reverse_preprocessing(mat_imp_raw, prep)
    
    # STORE the matrix for visualization later
    imputed_matrices[[name]] <<- mat_imp
    
    metric_vals <- compute_metrics(mat_true, mat_imp, mask, metrics)
    extra_vals  <- extra_metrics(mat_true, mat_imp, coords, mask)
    
    tibble(
      method = name, success = TRUE, runtime_sec = runtime,
      !!!as.list(metric_vals), !!!extra_vals
    )
  })
  
  return(list(metrics = method_results, matrices = imputed_matrices))
}

export_ext_dataset <- function(
    mat_true,
    mat_na,
    coords,
    method_name,
    preprocessing,
    prop,
    root_dir = "output") {
  
  p_str <- paste0("p", round(prop * 100))
  
  out_dir <- file.path(root_dir, method_name, p_str)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  mask <- !is.na(mat_na)
  
  # Apply preprocessing
  prep <- apply_preprocessing(mat_na, preprocessing)
  
  data_proc <- prep$data
  data_proc[is.na(data_proc)] <- 0
  
  write.csv(data_proc,
            file.path(out_dir, "data.csv"),
            row.names = FALSE)
  
  write.csv(mask,
            file.path(out_dir, "mask.csv"),
            row.names = FALSE)
  
  write.csv(coords,
            file.path(out_dir, "coords.csv"),
            row.names = FALSE)
  
  saveRDS(prep,
          file.path(out_dir, "prep.rds"))
  
  # eventually add the missing model (MCARM MNAR, hybrid, ...)
  metadata <- list(
    method = method_name,
    missing_proportion = prop,
    transform = preprocessing$transform,
    scaling = preprocessing$scaling,
    n_pixels = nrow(mat_na),
    n_features = ncol(mat_na),
    export_time = as.character(Sys.time())
  )
  
  writeLines(
    jsonlite::toJSON(metadata, pretty = TRUE, auto_unbox = TRUE),
    file.path(out_dir, "metadata.json")
  )
  
  message("Dataset exported for ", method_name, " | ", p_str)
}

import_ext_results <- function(
    mat_true,
    coords,
    metric_functions,
    io_root = "output",
    results_root = file.path(io_root, "results")) {
  
  results <- list()
  
  method_dirs <- list.dirs(
    results_root,
    recursive = FALSE,
    full.names = TRUE
  )
  
  for (method_path in method_dirs) {
    
    method_name <- basename(method_path)
    
    imp_files <- list.files(
      method_path,
      pattern = "^imputed_.*\\.csv$",
      full.names = TRUE
    )
    
    if (length(imp_files) == 0) {
      cat("No files in", method_name, "\n")
      next
    }
    
    # load IO folder once per method (NOT per p_dir anymore)
    io_method_dir <- file.path(io_root, method_name)
    
    for (file in imp_files) {
      
      cat("Importing:", basename(file), "\n")
      
      # extract p from filename
      p_str <- stringr::str_extract(basename(file), "p\\d+")
      prop <- as.numeric(stringr::str_remove(p_str, "p")) / 100
      
      io_p_dir <- file.path(io_method_dir, p_str)
      
      if (!dir.exists(io_p_dir)) {
        cat("Missing IO folder:", io_p_dir, "\n")
        next
      }
      
      prep_path <- file.path(io_p_dir, "prep.rds")
      mask_path <- file.path(io_p_dir, "mask.csv")
      
      if (!file.exists(prep_path) || !file.exists(mask_path)) {
        cat("Missing prep or mask for", p_str, "\n")
        next
      }
      
      prep <- readRDS(prep_path)
      mask <- !as.matrix(read.csv(mask_path))
      
      mat_imp <- as.matrix(read.csv(file))
      mat_imp <- reverse_preprocessing(mat_imp, prep)
      
      metric_vals <- compute_metrics(
        mat_true,
        mat_imp,
        mask,
        metric_functions
      )
      
      extra_vals <- extra_metrics(
        mat_true,
        mat_imp,
        coords,
        mask
      )
      
      results[[length(results) + 1]] <- tibble::tibble(
        method = method_name,
        prop = prop,
        !!!as.list(metric_vals),
        !!!extra_vals
      )
    }
  }
  
  bind_rows(results)
}