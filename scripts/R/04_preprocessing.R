# =============================================================================
# 04_preprocessing.R
# Purpose:  Implement data transformation and scaling functions for MSI datasets, 
#           along with their inverses for post-imputation reversal.
# Inputs:   - A numeric matrix (pixels × features) and a method specification list
# Outputs:  - A list containing the pre-processed matrix, the forward and backward functions
# Depends:  
# =============================================================================


# Functions ---------------------------------------------------------------

#' Apply specified transformation and scaling to the input matrix, and return the
#' pre-processed data along with the functions needed to reverse these steps after imputation.
#' 
#' Supported transformations: "none", "log1p"
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

#' Apply specified scaling to the input matrix, and return the scaled data along 
#' with the parameters needed to reverse the scaling after imputation.
#' 
#' Supported scalings: "none", "zscore", "pareto", "range"
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

#' Main function to apply the specified transformation and scaling to the input matrix.
#' 
#' @param mat         Numeric matrix (pixels × features) to be pre-processed
#' @param method_spec List with elements 'transform' and 'scaling' specifying
#'                    the desired transformation and scaling methods 
#'                    (e.g. list(transform = "log1p", scaling = "zscore"))
#'  @return A list containing:
#'        - data: The pre-processed matrix ready for imputation
#'        - transform: The transformation functions (forward and backward)
#'        - scaler: The scaling functions (forward and backward)
#'        - scale_params: The parameters needed to reverse the scaling after imputation
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

#' Reverse the transformation and scaling applied to the imputed matrix, using the
#' preprocessing object returned by apply_preprocessing() to access the necessary parameters and functions.
#' 
#' @param mat_imp  Numeric matrix (pixels × features) with imputed values, still in pre-processed space
#' @param prep     The list returned by apply_preprocessing() containing the transform and scaler
#' @return The imputed matrix transformed back to the original data space, ready for metric computation
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