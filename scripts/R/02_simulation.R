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
    
    idx_mcar <- sample(remaining, n_mcar)
    
    miss_idx <- unique(c(idx_mnar, idx_spatial, idx_mcar))
    
    mat_na[miss_idx, j] <- NA
  }
  
  mat_na
}