# =============================================================================
# 03_imputation.R
# Purpose:  Implement various imputation methods for missing values in MSI data
# Inputs:   - A numeric matrix with missing values (NA) and optional spatial coordinates
# Outputs:  - A complete numeric matrix with imputed values replacing NA
# Depends:  missForest, gstat, kernlab, VIM, imputeLCMD, pcaMethods, MsCoreUtils
# =============================================================================


# Functions ---------------------------------------------------------------

impute_zero <- function(mat_na, coords = NULL) {
  mat_na[is.na(mat_na)] <- 0
  mat_na
}

impute_mean <- function(mat_na, coords = NULL) {
  apply(mat_na, 2, function(col) {
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  }) %>% as.matrix()
}

impute_median <- function(mat_na, coords = NULL) {
  apply(mat_na, 2, function(col) {
    col[is.na(col)] <- median(col, na.rm = TRUE)
    col
  }) %>% as.matrix()
}

impute_half_min <- function(mat_na, coords = NULL) {
  mat_imp <- mat_na
  
  for (j in seq_len(ncol(mat_na))) {
    col <- mat_na[, j]
    min_val <- min(col, na.rm = TRUE)
    col[is.na(col)] <- min_val / 2
    mat_imp[, j] <- col
  }
  mat_imp
}

impute_missForest <- function(mat_na, coords = NULL) {
  res <- missForest(mat_na)
  return(res$ximp)
}

impute_knn <- function(mat_na, coords = NULL){
  df <- as.data.frame(mat_na)
  res <- kNN(df, k = 5, imp_var = FALSE)
  as.matrix(res)
}

# impute_qrilc <- function(mat_na, coords = NULL){
#   res <- impute.QRILC(t(mat_na))
#   return(t(res[[1]]))
# }

impute_qrilc <- function(mat_na, coords = NULL) {
  res <- MsCoreUtils::impute_QRILC(mat_na, MARGIN = 2L)
  return(res)
}

impute_ppca <- function(mat_na, coords = NULL){
  res <- pca(t(mat_na), method = "ppca", nPcs = 2)
  t(completeObs(res))
}

impute_bpca <- function(mat_na, coords = NULL){
  res <- pca(t(mat_na), method = "bpca", nPcs = 2)
  t(completeObs(res))
}

# impute_bpca <- function(mat_na, coords = NULL) {
#   res <- MsCoreUtils::impute_bpca(mat_na, MARGIN = 2L)
#   return(res)
# }

impute_svd <- function(mat_na, coords = NULL){
  res <- pca(t(mat_na), method = "svdImpute", nPcs = 5)
  return(t(completeObs(res)))
}

impute_nngp <- function(data, coords = NULL){
  res <- msImpute(t(data), method = "v2") # method!
  return(t(res))
}

#' Impute missing values using spatial KNN
#' 
#' For each missing value, find the k nearest observed pixels based on spatial 
#' coordinates and impute using their mean.
impute_spatial_knn <- function(mat_na, coords, k = 5) {
  mat_imp <- mat_na

  for (j in seq_len(ncol(mat_na))) {
    missing_idx <- which(is.na(mat_na[, j]))
    observed_idx <- which(!is.na(mat_na[, j]))

    for (i in missing_idx) {
      # Euclidean distances to observed points
      dists <- sqrt(
        (coords$X[observed_idx] - coords$X[i])^2 +
          (coords$Y[observed_idx] - coords$Y[i])^2
      )

      nn <- observed_idx[order(dists)][1:k]
      mat_imp[i, j] <- mean(mat_na[nn, j], na.rm = TRUE)
    }
  }

  mat_imp
}

# Gaussian Process regression for spatial imputation

#' Impute missing values using Gaussian Process regression based on spatial coordinates
#' 
#' For each feature, fit a Gaussian Process model to the observed values and 
#' predict missing ones based on their spatial coordinates. This method captures 
#' spatial correlations.
impute_gp <- function(mat_na, coords) {
  mat_imp <- mat_na

  for (j in seq_len(ncol(mat_na))) {
    y <- mat_na[, j]
    obs <- !is.na(y)

    if (sum(obs) < 5) next

    model <- suppressMessages(
      gausspr(
        x = as.matrix(coords[obs, ]),
        y = y[obs]
      )
    )

    pred <- predict(model, as.matrix(coords[!obs, ]))

    mat_imp[!obs, j] <- pred
  }

  mat_imp
}

#' Impute missing values using Inverse Distance Weighting (IDW) based on spatial coordinates
#' 
#' For each feature, use the observed values to predict missing ones by weighting them
#' inversely by their distance to the missing point.
impute_idw <- function(mat_na, coords, idp = 2) {
  mat_imp <- mat_na

  for (j in seq_len(ncol(mat_na))) {
    df <- data.frame(
      x = coords$X,
      y = coords$Y,
      z = mat_na[, j]
    )

    observed <- df[!is.na(df$z), ]
    missing <- df[is.na(df$z), ]

    if (nrow(missing) == 0) next

    coordinates(observed) <- ~ x + y
    coordinates(missing) <- ~ x + y

    idw_res <- idw(z ~ 1, observed, missing, idp = idp)

    mat_imp[is.na(mat_na[, j]), j] <- idw_res$var1.pred
  }

  mat_imp
}