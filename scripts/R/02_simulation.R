# =============================================================================
# 02_simulation.R
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

simulate_MCAR <- function(mat, prop = 0.1) {
  mat_na <- mat
  n <- nrow(mat)
  num_miss <- round(prop * n)

  for (j in seq_len(ncol(mat))) {
    # Ensure we don't try to sample more than available rows
    miss_idx <- sample(n, size = min(num_miss, n))
    mat_na[miss_idx, j] <- NA
  }
  mat_na
}

simulate_MNAR <- function(mat, prop = 0.1) {
  mat_na <- mat
  n <- nrow(mat)
  num_miss <- round(prop * n)
  
  for (j in seq_len(ncol(mat))) {
    # Find the indices of the smallest values in this column
    col_values <- mat[, j]
    threshold_idx <- order(col_values)[1:num_miss]
    mat_na[threshold_idx, j] <- NA
  }
  mat_na
}

simulate_hybrid_missing <- function(mat, total_prop = 0.1, mnar_weight = 0.7) {
  mat_na <- mat
  n <- nrow(mat)
  
  # Calculate how many NAs come from each logic
  prop_mnar <- total_prop * mnar_weight
  prop_mcar <- total_prop * (1 - mnar_weight)
  
  num_mnar <- round(prop_mnar * n)
  num_mcar <- round(prop_mcar * n)
  
  for (j in seq_len(ncol(mat))) {
    col_vals <- mat[, j]
    
    # 1. Apply MNAR (The lowest values)
    mnar_idx <- order(col_vals)[1:num_mnar]
    mat_na[mnar_idx, j] <- NA
    
    # 2. Apply MCAR (Randomly from the REMAINING indices)
    remaining_idx <- setdiff(seq_len(n), mnar_idx)
    mcar_idx <- sample(remaining_idx, size = min(num_mcar, length(remaining_idx)))
    mat_na[mcar_idx, j] <- NA
  }
  
  return(mat_na)
}

simulate_msi_missing <- function(
    mat,
    coords,
    prop = 0.2,
    spatial_fraction = 0.4,
    mnar_fraction = 0.5){
  
  mat_na <- mat
  n <- nrow(mat)
  
  total_miss <- round(prop * n)
  
  for(j in seq_len(ncol(mat))){
    
    values <- mat[,j]
    
    ## ---- 1 MNAR (low intensities)
    n_mnar <- round(total_miss * mnar_fraction)
    idx_mnar <- order(values)[1:n_mnar]
    
    ## ---- 2 spatial dropout
    n_spatial <- round(total_miss * spatial_fraction)
    
    center <- coords[sample(n,1),]
    
    d <- sqrt(
      (coords$X-center$X)^2 +
        (coords$Y-center$Y)^2
    )
    
    idx_spatial <- order(d)[1:n_spatial]
    
    ## ---- 3 MCAR remainder
    used <- unique(c(idx_mnar, idx_spatial))
    remaining <- setdiff(seq_len(n), used)
    
    n_mcar <- total_miss - length(used)
    # n_mcar <- max(0, total_miss - length(used))
    # if (n_mcar > 0) idx_mcar <- sample(remaining, min(n_mcar, length(remaining)))
    # else idx_mcar <- integer(0)
    
    idx_mcar <- sample(remaining, n_mcar)
    
    miss_idx <- unique(c(idx_mnar, idx_spatial, idx_mcar))
    
    mat_na[miss_idx, j] <- NA
  }
  
  mat_na
}

# need to clean it
simulate_MCAR_global <- function(data, prop = 0.1) {
  all_idx <- which(data != Inf, arr.ind = T)
  rdm_idx <- sample(1:nrow(all_idx), round(nrow(all_idx)*prop))
  slc_idx <- all_idx[rdm_idx, ]
  data_res <- data
  data_res[slc_idx] <- NA
  # return(list(data_res = data_res, mis_idx = slc_idx))
  return(data_res)
}

simulate_MNAR_global <- function (data_c, mis_var = 0.1, var_prop = seq(.3, .6, .1)) {
  data_mis <- data_c
  if (is.numeric(mis_var)) var_mis_list <- sample(1:ncol(data_c), round(ncol(data_c)*mis_var))
  else if (is.character(mis_var)) var_mis_list <- which(colnames(data_c) %in% mis_var)
  for (i in 1:length(var_mis_list)) {
    var_idx <- var_mis_list[i]
    cur_var <- data_mis[, var_idx]
    cutoff <- quantile(cur_var, sample(var_prop, 1))
    cur_var[cur_var < cutoff] <- NA
    data_mis[, var_idx] <- cur_var
  }
  mis_idx_df <- which(is.na(data_mis), arr.ind = T)
  # return (list(data_mis = data_mis, mis_idx_df = mis_idx_df))
  return(data_mis)
}