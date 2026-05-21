# =============================================================================
# 08_downstream.R
# Purpose:  Run downstream MSI analyses on imputed matrices and compare
#           results against the ground truth to assess biological validity.
# Inputs:   imputed matrices from results_storage, meta data frame, ground truth
# Outputs:  tidy downstream_results data frame + PCA comparison plots
# Depends:  tidyverse, vegan, ggrepel
# =============================================================================


# =========================================================================
# 1. PCA
# =========================================================================

#' Run PCA on a single matrix and return scores, variance explained, loadings.
#'
#' @param mat   Numeric matrix (ROIs x features), fully observed.
#' @param meta  Data frame with one row per ROI: Group, Regions, Individual, ROI_num.
#' @param n_pcs Number of PCs to retain (capped at matrix rank).
#' @return Named list: $scores (data frame, PC cols + meta), $variance_explained,
#'         $loadings (features x n_pcs matrix).
run_pca_downstream <- function(mat, meta = NULL, n_pcs = 10) {
  mat_log <- log1p(mat)
  mat_scaled <- scale(mat_log)

  pca <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
  n_pcs <- min(n_pcs, ncol(pca$x))

  scores <- as.data.frame(pca$x[, 1:n_pcs])
  if (!is.null(meta)) scores <- bind_cols(scores, meta)

  list(
    scores             = scores,
    variance_explained = summary(pca)$importance[2, 1:n_pcs],
    loadings           = pca$rotation[, 1:n_pcs]
  )
}

#' Compute symmetric Procrustes sum-of-squares between two PC score matrices.
#'
#' @param scores_true   Scores data frame from run_pca_downstream() on ground truth.
#' @param scores_method Scores data frame from run_pca_downstream() on imputed mat.
#' @param n_pcs         Number of PCs to include in the comparison.
#' @return Scalar: Procrustes SS (0 = identical, 1 = maximally different).
procrustes_dist <- function(scores_true, scores_method, n_pcs = 5) {
  pc_cols <- paste0("PC", seq_len(n_pcs))
  pc_cols <- intersect(pc_cols, intersect(colnames(scores_true), colnames(scores_method)))
  if (length(pc_cols) < 2) {
    return(NA_real_)
  }

  X <- as.matrix(scores_true[, pc_cols])
  Y <- as.matrix(scores_method[, pc_cols])

  vegan::procrustes(X, Y, symmetric = TRUE)$ss
}

#' Run PCA and return comparison metrics against ground truth.
#'
#' @param mat Imputed matrix.
#' @param meta Metadata data frame.
#' @param gt_pca Output of run_pca_downstream() on ground truth.
#' @return Named list: $result (full PCA output), scalar metrics ProcrustesSS,
#'         VarDiffPC1, VarDiffPC2.
compare_pca <- function(mat, meta, gt_pca) {
  res <- run_pca_downstream(mat, meta)

  list(
    result = res,
    ProcrustesSS = procrustes_dist(gt_pca$scores, res$scores),
    VarDiffPC1 = abs(gt_pca$variance_explained[["PC1"]] -
      res$variance_explained[["PC1"]]),
    VarDiffPC2 = abs(gt_pca$variance_explained[["PC2"]] -
      res$variance_explained[["PC2"]])
  )
}


#' # =========================================================================
#' # 2. Differential lipid analysis
#' # =========================================================================
#' 
#' #' Run differential analysis (limma-trend) between two groups.
#' #'
#' #' @param mat          Numeric matrix (ROIs x features), fully observed.
#' #' @param meta         Metadata data frame; must contain `contrast_col`.
#' #' @param contrast_col Column name in meta defining the grouping (e.g. "Group").
#' #' @param levels       Character vector of length 2: c(reference, test).
#' #'                     If NULL, uses the first two factor levels.
#' #' @return Data frame: feature, logFC, AveExpr, t, P.Value, adj.P.Val.




#' Run feature-by-feature t-tests or Wilcoxon tests between two groups.
run_de_downstream <- function(mat, meta, contrast_col = "Group", levels = NULL,
                              method = c("t.test", "wilcoxon")) {
  method <- match.arg(method)

  # Spatial pipeline: no group labels available → skip DE
  if (is.null(meta) || !contrast_col %in% colnames(meta)) {
    return(tibble())
  }

  group_vec <- meta[[contrast_col]]

  if (!is.null(levels)) {
    keep <- group_vec %in% levels
    mat <- mat[keep, , drop = FALSE]
    group_vec <- group_vec[keep]
  }

  group_fac <- factor(group_vec)
  if (nlevels(group_fac) < 2) {
    return(tibble())
  }

  mat_log <- log1p(mat)

  res <- purrr::map_dfr(colnames(mat_log), function(feat) {
    vals <- split(mat_log[, feat], group_fac)
    g_ref <- vals[[1]]
    g_test <- vals[[2]]

    logFC <- mean(g_test, na.rm = TRUE) - mean(g_ref, na.rm = TRUE)

    p_val <- tryCatch(
      {
        if (method == "t.test") {
          t.test(g_test, g_ref)$p.value
        } else {
          wilcox.test(g_test, g_ref, exact = FALSE)$p.value
        }
      },
      error = function(e) NA_real_
    )

    tibble(feature = feat, logFC = logFC, P.Value = p_val)
  }) |>
    mutate(adj.P.Val = p.adjust(P.Value, method = "BH"))

  return(res)
}

#' Assess how well the imputation preserves p-value ranking and significance
compare_de <- function(mat, meta, gt_de, contrast_col = "Group", levels = NULL,
                       method = c("t.test", "wilcoxon"), fdr_thresh = 0.05) {
  method <- match.arg(method)
  if (nrow(gt_de) == 0) {
    return(list(result = tibble(), Jaccard = NA_real_, RankCor = NA_real_))
  }

  res <- run_de_downstream(mat, meta, contrast_col = contrast_col, levels = levels, method = method)
  if (nrow(res) == 0) {
    return(list(result = res, Jaccard = NA_real_, RankCor = NA_real_))
  }

  # 1. Jaccard similarity of declared significant biomarkers
  sig_true <- gt_de$feature[gt_de$adj.P.Val < fdr_thresh & !is.na(gt_de$adj.P.Val)]
  sig_imp <- res$feature[res$adj.P.Val < fdr_thresh & !is.na(res$adj.P.Val)]

  n_intersect <- length(intersect(sig_true, sig_imp))
  n_union <- length(union(sig_true, sig_imp))
  jaccard <- if (n_union == 0) NA_real_ else n_intersect / n_union

  # 2. Spearman Correlation of raw p-values across all shared features
  common <- intersect(gt_de$feature, res$feature)
  rank_cor <- if (length(common) > 2) {
    cor(gt_de$P.Value[match(common, gt_de$feature)],
      res$P.Value[match(common, res$feature)],
      method = "spearman", use = "complete.obs"
    )
  } else {
    NA_real_
  }

  list(result = res, Jaccard = jaccard, RankCor = rank_cor)
}


# =========================================================================
# 3. Clustering
# =========================================================================

#' Cluster ROIs via k-means (Ward hclust as fallback option).
#'
#' @param mat    Numeric matrix (ROIs x features), fully observed.
#' @param k      Number of clusters. Defaults to number of unique Groups.
#' @param meta   Metadata (used only if k is NULL).
#' @param method "kmeans" or "hclust".
#' @return Named list: $labels (integer vector of cluster assignments),
#'         $silhouette_avg (mean silhouette width).

run_clustering_downstream <- function(mat, meta = NULL, k = NULL) {
  if (is.null(k)) {
    if (!is.null(meta) && "Group" %in% colnames(meta)) {
      k <- length(unique(meta$Group))
    } else {
      k <- 3L # sensible spatial default; override via k argument
    }
  }

  mat_scaled <- scale(log1p(mat))

  set.seed(42)


  labels <- kmeans(mat_scaled, centers = k, nstart = 25, iter.max = 100)$cluster
  sil <- tryCatch(
    mean(cluster::silhouette(labels, dist(mat_scaled))[, 3]),
    error = function(e) NA_real_
  )
  list(labels = labels, silhouette_avg = sil)
}


#' Compare clustering from an imputed matrix against ground truth clustering.
#'
#' Returns:
#'   ARI           — Adjusted Rand Index vs ground truth cluster labels
#'   ARI_Group     — ARI of imputed labels vs known meta$Group (biological)
#'   Silhouette    — mean silhouette width of imputed clustering
#'
#' @param mat         Imputed matrix.
#' @param meta        Metadata data frame.
#' @param gt_clust    Output of run_clustering_downstream() on ground truth.
compare_clustering <- function(mat, meta, gt_clust) {
  res <- run_clustering_downstream(mat, meta)

  # ari_gt    <- mclust::adjustedRandIndex(gt_clust$labels, res$labels)
  # ari_group <- mclust::adjustedRandIndex(as.integer(factor(meta$Group)), res$labels)

  list(
    result = res,
    # ARI        = ari_gt,
    # ARI_Group  = ari_group,
    Silhouette = res$silhouette_avg
  )
}


# =========================================================================
# 4. Co-abundance / correlation network
# =========================================================================

#' Build a Spearman co-abundance network.
#'
#' @param mat           Numeric matrix (ROIs x features), fully observed.
#' @param cor_threshold Minimum |correlation| to include an edge.
#' @return Named list: $cor_matrix, $edges (tibble: feat_a, feat_b, cor).
run_coabundance_downstream <- function(mat, cor_threshold = 0.7) {
  mat_log <- log1p(mat)
  cor_mat <- cor(mat_log, method = "spearman")

  idx <- which(upper.tri(cor_mat) & abs(cor_mat) >= cor_threshold, arr.ind = TRUE)
  edges <- tibble(
    feat_a = colnames(mat)[idx[, 1]],
    feat_b = colnames(mat)[idx[, 2]],
    cor    = cor_mat[idx]
  )

  list(cor_matrix = cor_mat, edges = edges)
}

#' Compare co-abundance network from an imputed matrix against ground truth.
#'
#' Returns:
#'   EdgeJaccard   — Jaccard similarity of edge sets (both thresholded)
#'   CorrFrobenius — Frobenius norm of |cor_imputed − cor_true| / n_features²
#'                   (normalised so it is comparable across feature sets)
#'
#' @param mat       Imputed matrix.
#' @param gt_coab   Output of run_coabundance_downstream() on ground truth.
compare_coabundance <- function(mat, gt_coab) {
  res <- run_coabundance_downstream(mat, cor_threshold = attr(gt_coab, "threshold") %||% 0.7)

  # Edge Jaccard
  edges_true <- paste(gt_coab$edges$feat_a, gt_coab$edges$feat_b, sep = "||")
  edges_imp <- paste(res$edges$feat_a, res$edges$feat_b, sep = "||")
  n_intersect <- length(intersect(edges_true, edges_imp))
  n_union <- length(union(edges_true, edges_imp))
  jaccard <- if (n_union == 0) NA_real_ else n_intersect / n_union

  list(
    result        = res,
    EdgeJaccard   = jaccard
  )
}

# Helper: %||% (NULL coalescing) if not already available
`%||%` <- function(a, b) if (!is.null(a)) a else b


# =========================================================================
# Orchestration
# =========================================================================

#' Run all downstream analyses on a single fully-observed matrix.
#'
#' @param mat            Numeric matrix (ROIs x features), no NAs.
#' @param meta           Metadata data frame.
#' @param contrast_col   Column for DE analysis (default "Group").
#' @param de_levels      Two-element vector of group names to contrast (or NULL
#'                       to use the first two factor levels automatically).
#' @param cor_threshold  Edge threshold for co-abundance network.
#' @return Named list: $pca, $de, $clustering, $coabundance.


#' Apply all downstream comparisons to every method x prop x replicate in
#' results_storage, returning a tidy metrics data frame and stored results
#' for each analysis type (for use in 09_downstream_visualisation.R).
#'
#' @param results_storage Nested list [[rep_key]][[p_str]]$imputed[[method]].
#' @param meta            Metadata data frame (same row order as mat_impute).
#' @param ground_truth    Output of run_all_downstream(mat_impute, meta).
#' @param contrast_col    Column for DE comparison (default "Group").
#' @param de_levels       Two-level vector for DE contrast (or NULL).
#' @return List:
#'   $metrics     — tidy data frame (method, prop, rep, all scalar metrics)
#'   $pca_store   — [[rep_key]][[p_str]][[method]] → full PCA result
#'   $de_store    — [[rep_key]][[p_str]][[method]] → full DE result
#'   $clust_store — [[rep_key]][[p_str]][[method]] → full clustering result
#'   $coab_store  — [[rep_key]][[p_str]][[method]] → full co-abundance result


run_all_downstream <- function(mat, meta = NULL, contrast_col = "Group",
                               de_levels = NULL, cor_threshold = 0.7) {
  has_de <- !is.null(meta) && contrast_col %in% colnames(meta)

  message("  [downstream] PCA ...")
  pca <- run_pca_downstream(mat, meta)

  if (has_de) {
    message("  [downstream] Differential analysis (", contrast_col, ") ...")
  } else {
    message("  [downstream] Differential analysis — skipped (no group labels in meta)")
  }

  de_t <- run_de_downstream(mat, meta,
    contrast_col = contrast_col,
    levels = de_levels, method = "t.test"
  )

  de_w <- run_de_downstream(mat, meta,
    contrast_col = contrast_col,
    levels = de_levels, method = "wilcoxon"
  )

  message("  [downstream] Clustering ...")
  clust <- run_clustering_downstream(mat, meta)

  message("  [downstream] Co-abundance network ...")
  coab <- run_coabundance_downstream(mat, cor_threshold = cor_threshold)
  attr(coab, "threshold") <- cor_threshold

  list(pca = pca, de_t = de_t, de_w = de_w, clustering = clust, coabundance = coab)
}

compute_downstream_comparison <- function(results_storage, meta = NULL, ground_truth,
                                          contrast_col = "Group", de_levels = NULL) {
  has_de <- !is.null(meta) && contrast_col %in% colnames(meta)
  pca_store <- list()
  de_t_store <- list()
  de_w_store <- list()
  clust_store <- list()
  coab_store <- list()

  metrics_df <- map_dfr(names(results_storage), function(rep_key) {
    rep_idx <- as.integer(str_extract(rep_key, "\\d+"))

    map_dfr(names(results_storage[[rep_key]]), function(p_str) {
      prop <- as.numeric(p_str)
      imputed_list <- results_storage[[rep_key]][[p_str]]$imputed
      if (is.null(imputed_list) || length(imputed_list) == 0) {
        return(tibble())
      }

      if (is.null(pca_store[[rep_key]])) pca_store[[rep_key]] <<- list()
      if (is.null(de_t_store[[rep_key]])) de_t_store[[rep_key]] <<- list()
      if (is.null(de_w_store[[rep_key]])) de_w_store[[rep_key]] <<- list()
      if (is.null(clust_store[[rep_key]])) clust_store[[rep_key]] <<- list()
      if (is.null(coab_store[[rep_key]])) coab_store[[rep_key]] <<- list()

      pca_store[[rep_key]][[p_str]] <<- list()
      de_t_store[[rep_key]][[p_str]] <<- list()
      de_w_store[[rep_key]][[p_str]] <<- list()
      clust_store[[rep_key]][[p_str]] <<- list()
      coab_store[[rep_key]][[p_str]] <<- list()

      map_dfr(names(imputed_list), function(method_name) {
        mat_imp <- imputed_list[[method_name]]

        message("  Downstream [", rep_key, "][p=", p_str, "] – ", method_name)

        tryCatch(
          {
            pca_res <- compare_pca(
              mat    = mat_imp,
              meta   = meta,
              gt_pca = ground_truth$pca
            )

            de_t_res <- compare_de(
              mat          = mat_imp,
              meta         = meta,
              gt_de        = ground_truth$de_t,
              contrast_col = contrast_col,
              levels       = de_levels,
              method       = "t.test"
            )

            de_w_res <- compare_de(
              mat          = mat_imp,
              meta         = meta,
              gt_de        = ground_truth$de_w,
              contrast_col = contrast_col,
              levels       = de_levels,
              method       = "wilcoxon"
            )

            cl_res <- compare_clustering(
              mat      = mat_imp,
              meta     = meta,
              gt_clust = ground_truth$clustering
            )

            cb_res <- compare_coabundance(
              mat     = mat_imp,
              gt_coab = ground_truth$coabundance
            )

            # Store full objects for visualization
            pca_store[[rep_key]][[p_str]][[method_name]] <<- pca_res$result
            de_t_store[[rep_key]][[p_str]][[method_name]] <<- de_t_res$result
            de_w_store[[rep_key]][[p_str]][[method_name]] <<- de_w_res$result
            clust_store[[rep_key]][[p_str]][[method_name]] <<- cl_res$result
            coab_store[[rep_key]][[p_str]][[method_name]] <<- cb_res$result

            tibble(
              method         = method_name,
              prop           = prop,
              rep            = rep_idx,
              ProcrustesSS   = pca_res$ProcrustesSS,
              TTest_Jaccard  = de_t_res$Jaccard,
              TTest_RankCor  = de_t_res$RankCor,
              Wilcox_Jaccard = de_w_res$Jaccard,
              Wilcox_RankCor = de_w_res$RankCor,
              Silhouette     = cl_res$Silhouette,
              EdgeJaccard    = cb_res$EdgeJaccard
            )
          },
          error = function(e) {
            warning("Downstream failed for ", method_name, " [", rep_key, "][", p_str, "]: ", e$message)
            tibble(
              method         = method_name,
              prop           = prop,
              rep            = rep_idx,
              ProcrustesSS   = NA_real_,
              TTest_Jaccard  = NA_real_,
              TTest_RankCor  = NA_real_,
              Wilcox_Jaccard = NA_real_,
              Wilcox_RankCor = NA_real_,
              Silhouette     = NA_real_,
              EdgeJaccard    = NA_real_
            )
          }
        )
      })
    })
  })

  list(
    metrics     = metrics_df,
    pca_store   = pca_store,
    de_t_store  = de_t_store,
    de_w_store  = de_w_store,
    clust_store = clust_store,
    coab_store  = coab_store
  )
}