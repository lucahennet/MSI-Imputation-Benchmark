# =============================================================================
# 04_preprocessing.R
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

transformers <- list(
  none = list(
    forward  = function(x) x,
    backward = function(x) x
  ),
  
  log1p = list(
    forward  = function(x) log1p(x),
    backward = function(x) expm1(x)
  )
)

scalers <- list(
  none = list(
    forward = function(x) list(data = x, params = NULL),
    backward = function(x, params) x
  ),
  
  zscore = list(
    forward = function(x){
      mu <- colMeans(x, na.rm = TRUE)
      sd <- apply(x, 2, sd, na.rm = TRUE)
      sd[sd == 0 | is.na(sd)] <- 1
      list(data = sweep(sweep(x, 2, mu, "-"), 2, sd, "/"), 
           params = list(mu = mu, measure = sd))
    },
    backward = function(x, params) sweep(sweep(x, 2, params$measure, "*"), 2, params$mu, "+")
  ),
  
  pareto = list(
    # Pareto uses the square root of the SD as the scaling factor
    forward = function(x){
      mu <- colMeans(x, na.rm = TRUE)
      sd <- apply(x, 2, sd, na.rm = TRUE)
      sq_sd <- sqrt(sd)
      sq_sd[sq_sd == 0 | is.na(sq_sd)] <- 1
      list(data = sweep(sweep(x, 2, mu, "-"), 2, sq_sd, "/"), 
           params = list(mu = mu, measure = sq_sd))
    },
    backward = function(x, params) sweep(sweep(x, 2, params$measure, "*"), 2, params$mu, "+")
  ),
  
  range = list(
    # Scales everything between 0 and 1
    forward = function(x){
      min_val <- apply(x, 2, min, na.rm = TRUE)
      max_val <- apply(x, 2, max, na.rm = TRUE)
      diff_val <- max_val - min_val
      diff_val[diff_val == 0 | is.na(diff_val)] <- 1
      list(data = sweep(sweep(x, 2, min_val, "-"), 2, diff_val, "/"), 
           params = list(min = min_val, measure = diff_val))
    },
    backward = function(x, params) sweep(sweep(x, 2, params$measure, "*"), 2, params$min, "+")
  )
)

apply_preprocessing <- function(mat, method_spec){
  
  # Transform
  tr <- transformers[[method_spec$transform]]
  
  mat_t <- tr$forward(mat)
  
  # Scale
  sc <- scalers[[method_spec$scaling]]
  
  scaled <- sc$forward(mat_t)
  
  list(
    data = scaled$data,
    transform = tr,
    scaler = sc,
    scale_params = scaled$params
  )
}

reverse_preprocessing <- function(mat_imp, prep){
  
  # Undo scaling
  mat_imp <- prep$scaler$backward(
    mat_imp,
    prep$scale_params
  )
  
  # Undo transform
  mat_imp <- prep$transform$backward(mat_imp)
  
  mat_imp
}