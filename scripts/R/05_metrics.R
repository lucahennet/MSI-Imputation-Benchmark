# =============================================================================
# 05_metrics.R
# Purpose:  Implement evaluation metrics to assess imputation performance, 
#           including both value accuracy and spatial structure preservation.
# Inputs:   - True complete matrix, imputed matrix, spatial coordinates, and mask of missing values
# Outputs:  - A list of computed metric values for each imputation method
# Depends:  SpatialPack, spdep, tidyverse
# =============================================================================


# Functions ---------------------------------------------------------------

#' Convert a vector of values and corresponding coordinates into an image matrix.
#' 
#' This function takes a vector of values and a data frame of coordinates (X, Y) 
#' and reconstructs the original image matrix. It assumes that the coordinates are 
#' integer pixel positions and that the values correspond to those positions. 
#' The resulting matrix will have rows corresponding to Y and columns corresponding 
#' to X, with NA for any missing positions.
#' 
#' @param values A numeric vector of values corresponding to pixel intensities.
#' @param coords A data frame with columns X and Y indicating the pixel coordinates for each
#' value in the 'values' vector. The length of 'values' should match the number of rows in 'coords'.
#' @return A matrix where the entry at (Y, X) corresponds to the value
#' provided in 'values' for that coordinate. The matrix will be filled with NA for any coordinates not provided.
#' 
#' Example usage:
#' coords <- data.frame(X = c(1, 2, 1), Y = c(1, 1, 2))
#' values <- c(10, 20, 30)
#' to_image(values, coords)
#' #' This will produce a matrix:
#' #      [,1] [,2]
#' # [1,]   10   20
#' # [2,]   30   NA
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

#' Compute the Structural Similarity Index (SSIM) between two matrices, treating them as images.
#' 
#' This function first scales the input matrices to the range [0, 255] to mimic 
#' typical image pixel values, then converts them into image matrices using the 
#' provided coordinates, and finally computes the SSIM for each feature (column) 
#' and averages the results. The SSIM is a perceptual metric that considers luminance, 
#' contrast, and structure to assess the similarity between two images.
#' 
#' @param true A numeric matrix of true values (pixels × features).
#' @param pred A numeric matrix of predicted values (pixels × features).
#' @param coords A data frame with columns X and Y indicating the pixel coordinates 
#' for each row in the 'true' and 'pred' matrices. The number of rows in 'true' 
#' and 'pred' should match the number of rows in 'coords'.
#' @return A single numeric value representing the average SSIM across all features. 
#' Higher values indicate greater similarity between the true and predicted images, 
#' with a maximum of 1 for identical images.
compute_ssim <- function(true, pred, coords) {
  scale_image <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) {
      return(x * 0)
    }
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

#' Compute the Concordance Correlation Coefficient (CCC) between true and predicted 
#' values, considering only the masked entries.
#' 
#' The CCC is a measure of agreement between two sets of data, combining elements 
#' of both precision and accuracy. It ranges from -1 to 1, where 1 indicates perfect 
#' concordance, 0 indicates no concordance, and -1 indicates perfect discordance. 
#' This function computes the CCC for the entries specified by the mask, allowing 
#' for an assessment of how well the predicted values match the true values in the context of imputation.
#' 
#' @param true A numeric matrix of true values (pixels × features).
#' @param pred A numeric matrix of predicted values (pixels × features).
#' @param mask A logical matrix of the same dimensions as 'true' and 'pred indicating 
#' which entries to include in the CCC calculation (typically the positions of missing values).
#' @return A single numeric value representing the CCC for the masked entries. 
#' Values closer to 1 indicate better agreement between the true and predicted 
#' values, while values closer to -1 indicate worse agreement.
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

#' Compute the Spectral Angle Mapper (SAM) between true and predicted values,
#' considering only the masked entries.
#' 
#' SAM is a measure of similarity between two vectors, often used in remote sensing
#' to compare spectral signatures. It calculates the angle between the true and
#' predicted vectors in a multi-dimensional space. A smaller angle indicates greater
#' similarity, with 0 degrees indicating identical vectors. 
#' 
#' @param true A numeric matrix of true values (pixels × features).
#' @param pred A numeric matrix of predicted values (pixels × features).
#' @param mask A logical matrix of the same dimensions as 'true' and 'pred
#' indicating which entries to include in the SAM calculation (typically the positions of missing values).
#' @return A single numeric value representing the SAM in radians for the masked entries.
compute_sam <- function(true, pred, mask) {
  t <- true[mask]
  p <- pred[mask]

  acos(sum(t * p) / (sqrt(sum(t^2)) * sqrt(sum(p^2))))
}

#' Compute the preservation of spatial autocorrelation using Moran's I between true and predicted values.
#' 
#' @param true A numeric matrix of true values (pixels × features).
#' @param pred A numeric matrix of predicted values (pixels × features).
#' @param coords A data frame with columns X and Y indicating the pixel coordinates 
#' for each row in the 'true' and 'pred' matrices. The number of rows in 'true' 
#' and 'pred' should match the number of rows in 'coords'.
#' @param k The number of nearest neighbours to use for constructing spatial weights 
#' (default is 6). This determines how many neighbouring pixels are considered 
#' when calculating Moran's I, with higher values capturing broader spatial relationships.
#' @return A single numeric value representing the average absolute difference in 
#' Moran's I between the true and predicted data across all features. Values closer 
#' to 0 indicate better preservation of spatial autocorrelation in the imputed data.
compute_moran_preservation <- function(true, pred, coords, k = 6) {
  # Build spatial weights
  listw <- nb2listw(knn2nb(knearneigh(coords, k = k)), style = "W")

  moran_vals <- numeric(ncol(true))

  for (j in seq_len(ncol(true))) {
    mt <- moran.test(true[, j], listw)$estimate[1]
    mp <- moran.test(pred[, j], listw)$estimate[1]

    moran_vals[j] <- abs(mt - mp)
  }

  mean(moran_vals, na.rm = TRUE)
}

#' Compute the variance ratio between predicted and true values for each feature,
#' and average across features.
#' 
#' The variance ratio is calculated as the variance of the predicted values divided by
#' the variance of the true values for each feature. A ratio close to 1 indicates
#' that the imputed data has a similar level of variability as the true data, while
#' ratios significantly greater than 1 suggest overestimation of variance, and ratios
#' significantly less than 1 suggest underestimation.
#' 
#' @param true A numeric matrix of true values (pixels × features).
#' @param pred A numeric matrix of predicted values (pixels × features).
#' @param mask A logical matrix of the same dimensions as 'true' and 'pred
#' indicating which entries to include in the variance ratio calculation (typically 
#' the positions of missing values for which we want to assess the imputation). 
#' If NULL, all entries are used.
#' @return A single numeric value representing the average variance ratio across all features.
compute_variance_ratio <- function(true, pred, mask = NULL) {
  ratios <- sapply(seq_len(ncol(true)), function(j) {
    var(pred[, j]) / var(true[, j])
  })

  mean(ratios, na.rm = TRUE)
}

#' Compute the preservation of correlation structure between true and predicted values.
#' 
#' This function calculates the correlation matrix for both the true and predicted data,
#' then extracts the upper triangular elements (excluding the diagonal) to compare the
#' pairwise correlations between features.
#' 
#' @param true A numeric matrix of true values (pixels × features).
#' @param pred A numeric matrix of predicted values (pixels × features).
#' @param mask A logical matrix of the same dimensions as 'true' and 'pred
#' indicating which entries to include in the correlation preservation calculation 
#' (typically the positions of missing values for which we want to assess the imputation). 
#' If NULL, all entries are used.
#' @return A single numeric value representing the correlation between the upper 
#' triangular elements of the true and predicted correlation matrices. Values 
#' closer to 1 indicate better preservation of the correlation structure in the 
#' imputed data, while values closer to 0 indicate poor preservation.
compute_correlation_preservation <- function(true, pred, mask = NULL) {
  cor_true <- cor(true)
  cor_pred <- cor(pred)

  idx <- upper.tri(cor_true)

  cor(cor_true[idx], cor_pred[idx])
}

#' A list of metric functions to compute various evaluation metrics for imputation performance.
#' 
#' Each function takes the true values, predicted values, and a mask indicating 
#' which entries to include in the calculation.
#' 
#' - RMSE: Root Mean Squared Error, measuring the average magnitude of the errors 
#' between true and predicted values.
#' - NRMSE: Normalized RMSE, which scales the RMSE by the standard deviation of 
#' the true values to provide a relative measure of error.
#' - MAE: Mean Absolute Error, measuring the average absolute difference between true and predicted
#' values, providing a more interpretable measure of error in the same units as the data.
#' - CCC: Concordance Correlation Coefficient
#' - VarRatio: Variance Ratio
#' - CorStruct: Correlation Structure Preservation
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

#' Compute additional metrics that require spatial coordinates
#' 
#' This function computes metrics that assess not only the accuracy of the imputed 
#' values but also the preservation of spatial structure. It includes:
#' - SSIM: Structural Similarity Index
#' - MoranDiff: The absolute difference in Moran's I between the true and predicted data
#' 
#' @param true A numeric matrix of true values (pixels × features).
#' @param pred A numeric matrix of predicted values (pixels × features).
#' @param coords A data frame with columns X and Y indicating the pixel coordinates 
#' for each row in the 'true' and 'pred' matrices. The number of rows in 'true' 
#' and 'pred' should match the number of rows in 'coords'.
#' @param mask A logical matrix of the same dimensions as 'true' and 'pred
#' indicating which entries to include in the calculation of metrics that require 
#' masking (if applicable). This parameter is not used for SSIM and MoranDiff but 
#' is included for consistency with other metric functions.
#' @return A list containing the computed SSIM and MoranDiff values.
extra_metrics <- function(true, pred, coords, mask) {
  has_spatial_coords <- !is.null(coords) &&
    is.data.frame(coords) &&
    all(c("X", "Y") %in% names(coords)) &&
    is.numeric(coords$X) && is.numeric(coords$Y)

  result <- list()

  if (has_spatial_coords) {
    result$SSIM <- compute_ssim(true, pred, coords)
    result$MoranDiff <- compute_moran_preservation(true, pred, coords)
  }

  result
}

#' Compute a set of specified metrics for the true and predicted matrices, using 
#' the provided mask to determine which entries to include in the calculations.
#' 
#' This function iterates over the list of metric functions defined in 'metric_functions' 
#' and applies each one to the true and predicted matrices, using the mask to focus 
#' on the relevant entries (typically the positions of missing values). The results 
#' are returned as a named vector, where each name corresponds to a metric and each 
#' value is the computed metric for the given true and predicted data.
#' 
#' @param true A numeric matrix of true values (pixels × features).
#' @param pred A numeric matrix of predicted values (pixels × features).
#' @param mask A logical matrix of the same dimensions as 'true' and 'pred
#' indicating which entries to include in the metric calculations (typically the 
#' positions of missing values). 
#' @param metrics A character vector of metric names to compute, corresponding to the keys in
#' the 'metric_functions' list. If NULL, all metrics in the list will be computed.
#' @return A named numeric vector where each name corresponds to a metric and each value is
#' the computed metric for the given true and predicted data. The metrics are 
#' computed only for the entries specified by the mask.
compute_metrics <- function(true, pred, mask, metrics) {
  map_dbl(metrics, function(f) {
    f(true, pred, mask)
  })
}