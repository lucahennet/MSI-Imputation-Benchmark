# =============================================================================
# 05_metrics.R
# Purpose:  Load, assemble and reshape MSI data from raw ROI files
# Inputs:   - A root directory containing one TXT subfolder and N CSV subfolders
# Outputs:  - raw_wide       : wide-format data frame (pixels × features)
#           - df_visualisation: long-format, TIC-normalised, for plotting
#           - df_impute      : clean wide matrix ready for imputation
#           - coords         : data frame of X/Y pixel coordinates
#           - mat_impute     : numeric matrix (pixels × features), no metadata
# Depends:  tidyverse, stringr
# =============================================================================

# Value accuracy → RMSE, NRMSE, MAE
# Spatial structure → SSIM

# Functions ---------------------------------------------------------------

# Convert long format (values + coordinates) back to image matrix
to_image <- function(values, coords) {
  
  # Unique sorted coordinates
  x_vals <- sort(unique(coords$X))
  y_vals <- sort(unique(coords$Y))
  
  # Create empty grid
  mat <- matrix(NA, nrow = length(y_vals), ncol = length(x_vals))
  
  # Map each (X,Y) to matrix position
  for (i in seq_along(values)) {
    x_idx <- match(coords$X[i], x_vals)
    y_idx <- match(coords$Y[i], y_vals)
    
    mat[y_idx, x_idx] <- values[i]
  }
  
  # Replace remaining NA (if any) with 0
  mat[is.na(mat)] <- 0
  
  mat
}

compute_ssim <- function(true, pred, coords) {
  
  scale_image <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(x * 0)
    (x - rng[1]) / diff(rng) * 255
  }
  
  vals <- numeric(ncol(true))
  
  for (j in seq_len(ncol(true))) {
    
    img_true <- scale_image(to_image(true[, j], coords))
    img_pred <- scale_image(to_image(pred[, j], coords))
    
    vals[j] <- SpatialPack::SSIM(img_true, img_pred)$SSIM
  }
  
  mean(vals, na.rm = TRUE)
}

# Concordance Correlation Coefficient (CCC)
compute_ccc <- function(true, pred, mask) {
  x <- true[mask]
  y <- pred[mask]
  
  mu_x <- mean(x)
  mu_y <- mean(y)
  
  var_x <- var(x)
  var_y <- var(y)
  
  cov_xy <- cov(x, y)
  
  (2 * cov_xy) / (var_x + var_y + (mu_x - mu_y)^2)
}

# Spectral Angle Mapper (SAM)
compute_sam <- function(true, pred, mask){
  
  t <- true[mask]
  p <- pred[mask]
  
  acos(sum(t*p) / (sqrt(sum(t^2))*sqrt(sum(p^2))))
}

# rebuilds spatial weights on every call -> slow?
compute_moran_preservation <- function(true, pred, coords, k = 6){
  
  # Build spatial weights
  listw <- nb2listw(knn2nb(knearneigh(coords, k = k)), style = "W")
  
  moran_vals <- numeric(ncol(true))
  
  for(j in seq_len(ncol(true))){
    
    mt <- moran.test(true[,j], listw)$estimate[1]
    mp <- moran.test(pred[,j], listw)$estimate[1]
    
    moran_vals[j] <- abs(mt - mp)
  }
  
  mean(moran_vals, na.rm = TRUE)
}

compute_variance_ratio <- function(true, pred, mask = NULL){
  
  ratios <- sapply(seq_len(ncol(true)), function(j){
    var(pred[,j]) / var(true[,j])
  })
  
  mean(ratios, na.rm = TRUE)
}

compute_correlation_preservation <- function(true, pred, mask = NULL){
  
  cor_true <- cor(true)
  cor_pred <- cor(pred)
  
  idx <- upper.tri(cor_true)
  
  cor(cor_true[idx], cor_pred[idx])
}

metric_functions <- list(
  RMSE = function(true, pred, mask) {
    sqrt(mean((true[mask] - pred[mask])^2))
  },
  
  NRMSE = function(true, pred, mask) {
    sqrt(mean((true[mask] - pred[mask])^2)) / sd(true[mask])
  },
  
  MAE = function(true, pred, mask) {
    mean(abs(true[mask] - pred[mask]))
  },
  
  CCC = compute_ccc,
  
  SAM = compute_sam,
  
  VarRatio = compute_variance_ratio,
  
  CorStruct = compute_correlation_preservation
)

extra_metrics <- function(true, pred, coords, mask) {
  list(
    SSIM = compute_ssim(true, pred, coords),
    MoranDiff = compute_moran_preservation(true, pred, coords)
  )
}

# Computes all metrics dynamically
compute_metrics <- function(true, pred, mask, metrics) {
  
  map_dbl(metrics, function(f) {
    f(true, pred, mask)
  })
}