# =============================================================================
# 02_simulation.R
# Purpose:  Simulate missing data in MSI datasets under different mechanisms
# Inputs:   - A complete numeric matrix (pixels × features)
#           - Proportion of missingness to simulate
# Outputs:  - A matrix with simulated missing values (NA)
# Depends:  
# =============================================================================


# Functions ---------------------------------------------------------------

#' Simulate MCAR by randomly setting a proportion of values to NA in each column
#' 
#' @param mat A complete numeric matrix (pixels × features)
#' @param prop Proportion of missing values to simulate (default 0.1 for 10%)
#' @return A matrix with simulated MCAR missing values (NA)
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

#' Simulate MNAR by setting the lowest values in each column to NA
#' 
#' @param mat A complete numeric matrix (pixels × features)
#' @param prop Proportion of missing values to simulate (default 0.1 for 10%)
#' @return A matrix with simulated MNAR missing values (NA)
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

#' Simulate a hybrid missingness mechanism combining MNAR and MCAR
#' 
#' @param mat A complete numeric matrix (pixels × features)
#' @param total_prop Total proportion of missing values to simulate (default 0.1
#' for 10%)
#' @param mnar_weight Proportion of missingness that is MNAR (default 0.7 means 
#' 70% MNAR and 30% MCAR)
#' @return A matrix with simulated hybrid missing values (NA)
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

#' Simulate a more complex missingness mechanism combining MNAR, spatial dropout, and MCAR
#' 
#' @param mat A complete numeric matrix (pixels × features)
#' @param coords A data frame with columns X and Y for pixel coordinates
#' @param prop Total proportion of missing values to simulate (default 0.2 for 20%)
#' @param spatial_fraction Proportion of missingness due to spatial dropout (default 0.4 
#' means 40% spatial, 60% MNAR/MCAR)
#' @param mnar_fraction Proportion of non-spatial missingness that is MNAR
#' (default 0.5 means 50% MNAR and 50% MCAR among the non-spatial missingness)
#' @return A matrix with simulated complex missing values (NA)
simulate_msi_missing <- function(
  mat,
  coords,
  prop = 0.2,
  spatial_fraction = 0.4,
  mnar_fraction = 0.5
) {
  mat_na <- mat
  n <- nrow(mat)

  total_miss <- round(prop * n)

  for (j in seq_len(ncol(mat))) {
    values <- mat[, j]

    # 1 MNAR (low intensities)
    n_mnar <- round(total_miss * mnar_fraction)
    idx_mnar <- order(values)[1:n_mnar]

    # 2 spatial dropout
    n_spatial <- round(total_miss * spatial_fraction)

    center <- coords[sample(n, 1), ]

    d <- sqrt(
      (coords$X - center$X)^2 +
        (coords$Y - center$Y)^2
    )

    idx_spatial <- order(d)[1:n_spatial]

    # 3 MCAR remainder
    used <- unique(c(idx_mnar, idx_spatial))
    remaining <- setdiff(seq_len(n), used)

    n_mcar <- total_miss - length(used)

    idx_mcar <- sample(remaining, n_mcar)

    miss_idx <- unique(c(idx_mnar, idx_spatial, idx_mcar))

    mat_na[miss_idx, j] <- NA
  }

  mat_na
}

#' Simulate MCAR missingness by randomly setting a proportion of values to NA across the entire matrix
#' 
#' @param data A complete numeric matrix (pixels × features)
#' @param prop Proportion of missing values to simulate (default 0.1 for 10%)
#' @return A matrix with simulated MCAR missing values (NA)
simulate_MCAR_global <- function(data, prop = 0.1) {
  all_idx <- which(data != Inf, arr.ind = T)
  rdm_idx <- sample(1:nrow(all_idx), round(nrow(all_idx) * prop))
  slc_idx <- all_idx[rdm_idx, ]
  data_res <- data
  data_res[slc_idx] <- NA
  return(data_res)
}