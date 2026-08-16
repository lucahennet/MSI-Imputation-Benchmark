# =============================================================================
# 08_downstream.R
# Purpose:  Run downstream MSI analyses on imputed matrices and compare
#           results against the ground truth to assess biological validity.
# Inputs:   imputed matrices from results_storage, meta data frame, ground truth
# Outputs:  downstream_results list — $metrics (tidy) + per-analysis stores
# Depends:  tidyverse, vegan, cluster
# =============================================================================


# ═════════════════════════════════════════════════════════════════════════════
# Helpers
# ═════════════════════════════════════════════════════════════════════════════

# NULL-coalescing (must precede first use)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Adjusted Rand Index — standalone, no extra dependency.
.ari <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n <= 1) {
    return(0)
  }
  sum_nij2 <- sum(choose(tab, 2))
  sum_ai2 <- sum(choose(rowSums(tab), 2))
  sum_bj2 <- sum(choose(colSums(tab), 2))
  expected <- (sum_ai2 * sum_bj2) / choose(n, 2)
  max_index <- 0.5 * (sum_ai2 + sum_bj2)
  if (max_index == expected) {
    return(0)
  }
  (sum_nij2 - expected) / (max_index - expected)
}

#' Standard log1p + z-score preprocessing shared by all analyses.
.preprocess <- function(mat) scale(log1p(mat))

#' Jaccard similarity between two character vectors.
.jaccard <- function(a, b) {
  n_union <- length(union(a, b))
  if (n_union == 0) NA_real_ else length(intersect(a, b)) / n_union
}


# ═════════════════════════════════════════════════════════════════════════════
# Analyses  (run_*  →  compute ground-truth object)
# ═════════════════════════════════════════════════════════════════════════════

#' Run PCA and return scores, variance explained, and loadings.
#'
#' @param mat   Numeric matrix (ROIs × features), fully observed.
#' @param meta  Optional data frame appended column-wise to scores.
#' @param n_pcs Number of PCs to retain (capped at matrix rank).
#' @return List: $scores, $variance_explained, $loadings.
run_pca_downstream <- function(mat, meta = NULL, n_pcs = 10) {
  pca <- prcomp(.preprocess(mat), center = FALSE, scale. = FALSE)
  n_pcs <- min(n_pcs, ncol(pca$x))
  scores <- as.data.frame(pca$x[, 1:n_pcs, drop = FALSE])
  if (!is.null(meta)) scores <- bind_cols(scores, meta)
  list(
    scores             = scores,
    variance_explained = summary(pca)$importance[2, 1:n_pcs],
    loadings           = pca$rotation[, 1:n_pcs, drop = FALSE]
  )
}

#' Run feature-by-feature t-tests or Wilcoxon tests between two groups.
#'
#' Returns an empty tibble when meta lacks the contrast column, so spatial
#' (label-free) pipelines work without any special-casing by callers.
#'
#' @param mat          Numeric matrix (ROIs × features), fully observed.
#' @param meta         Metadata data frame, or NULL.
#' @param contrast_col Column name in meta that encodes the grouping.
#' @param levels       c(reference, test) — NULL uses first two factor levels.
#' @param method       "t.test" or "wilcoxon".
#' @return Tibble: feature, logFC, P.Value, adj.P.Val.
run_de_downstream <- function(mat, meta, contrast_col = "Group", levels = NULL,
                              method = c("t.test", "wilcoxon")) {
  method <- match.arg(method)
  if (is.null(meta) || !contrast_col %in% colnames(meta)) {
    return(tibble())
  }
  
  # Subset to specified levels if provided, and ensure the grouping factor is properly formed
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
  
  # Log-transform once for all tests; t-tests on log-transformed data approximate log-fold changes better
  mat_log <- log1p(mat)
  purrr::map_dfr(colnames(mat_log), function(feat) {
    vals <- split(mat_log[, feat], group_fac)
    logFC <- mean(vals[[2]], na.rm = TRUE) - mean(vals[[1]], na.rm = TRUE)
    p_val <- tryCatch(
      if (method == "t.test") {
        t.test(vals[[2]], vals[[1]])$p.value
      } else {
        wilcox.test(vals[[2]], vals[[1]], exact = FALSE)$p.value
      },
      error = function(e) NA_real_
    )
    tibble(feature = feat, logFC = logFC, P.Value = p_val)
  }) |>
    mutate(adj.P.Val = p.adjust(P.Value, method = "BH"))
}

#' Cluster ROIs via k-means; k defaults to the number of unique groups in meta.
#'
#' @param mat  Numeric matrix (ROIs × features), fully observed.
#' @param meta Optional metadata — "Group" column sets k when present.
#' @param k    Explicit cluster count (overrides meta-derived k).
#' @return List: $labels (integer vector), $silhouette_avg.
run_clustering_downstream <- function(mat, meta = NULL, k = NULL) {
  k <- k %||%
    if (!is.null(meta) && "Group" %in% colnames(meta)) {
      length(unique(meta$Group))
    } else {
      3L
    }
  
  mat_s <- .preprocess(mat)
  set.seed(42)
  labels <- kmeans(mat_s, centers = k, nstart = 25, iter.max = 100)$cluster
  sil <- tryCatch(
    mean(cluster::silhouette(labels, dist(mat_s))[, 3]),
    error = function(e) NA_real_
  )
  list(labels = labels, silhouette_avg = sil)
}

#' Build a Spearman co-abundance network (upper-triangle edges above threshold).
#'
#' @param mat           Numeric matrix (ROIs × features), fully observed.
#' @param cor_threshold Minimum |correlation| to retain an edge.
#' @return List: $cor_matrix, $edges (tibble: feat_a, feat_b, cor), $threshold.
run_coabundance_downstream <- function(mat, cor_threshold = 0.7) {
  cor_mat <- cor(log1p(mat), method = "spearman")
  idx <- which(upper.tri(cor_mat) & abs(cor_mat) >= cor_threshold, arr.ind = TRUE)
  edges <- tibble(
    feat_a = colnames(mat)[idx[, 1]],
    feat_b = colnames(mat)[idx[, 2]],
    cor    = cor_mat[idx]
  )
  list(cor_matrix = cor_mat, edges = edges, threshold = cor_threshold)
}


# ═════════════════════════════════════════════════════════════════════════════
# Comparisons  (compare_*  →  scalar metrics + stored full result)
# ═════════════════════════════════════════════════════════════════════════════

#' Procrustes SS + variance-explained differences between imputed and GT PCA.
#'
#' @param mat    Imputed matrix.
#' @param meta   Metadata data frame.
#' @param gt     Ground-truth object ($pca from run_all_downstream()).
#' @return List: $result, $ProcrustesSS, $VarDiffPC1, $VarDiffPC2.
compare_pca <- function(mat, meta, gt) {
  res <- run_pca_downstream(mat, meta)
  
  # Procrustes analysis on the common PCs (up to 5)
  pc_cols <- intersect(
    paste0("PC", 1:5),
    intersect(colnames(gt$scores), colnames(res$scores))
  )
  # Only compute Procrustes SS if we have at least 2 common PCs; not meaningful otherwise
  ss <- if (length(pc_cols) >= 2) {
    vegan::procrustes(as.matrix(gt$scores[, pc_cols]),
                      as.matrix(res$scores[, pc_cols]),
                      symmetric = TRUE
    )$ss
  } else {
    NA_real_
  }
  
  list(
    # Store only the scores data frame — the plots never use loadings or
    # variance_explained from the store (scalars are already in $metrics).
    result       = list(scores = res$scores),
    ProcrustesSS = ss,
    VarDiffPC1   = abs(gt$variance_explained[["PC1"]] - res$variance_explained[["PC1"]]),
    VarDiffPC2   = abs(gt$variance_explained[["PC2"]] - res$variance_explained[["PC2"]])
  )
}

#' Jaccard hit-list overlap + Spearman rank correlation of p-values vs. GT DE
#'
#' @param mat          Imputed matrix.
#' @param meta         Metadata data frame.
#' @param gt           Ground-truth DE tibble (from run_de_downstream()).
#' @param contrast_col,levels,method,fdr_thresh  Forwarded to run_de_downstream().
#' @return List: $result, $Jaccard, $RankCor.
compare_de <- function(mat, meta, gt, contrast_col = "Group", levels = NULL,
                       method = c("t.test", "wilcoxon"), fdr_thresh = 0.05) {
  method <- match.arg(method)
  if (nrow(gt) == 0) {
    return(list(result = tibble(), Jaccard = NA_real_, RankCor = NA_real_))
  }
  
  # Run DE on the imputed matrix with the same parameters as the GT to ensure a fair comparison
  res <- run_de_downstream(mat, meta,
                           contrast_col = contrast_col,
                           levels = levels, method = method
  )
  if (nrow(res) == 0) {
    return(list(result = res, Jaccard = NA_real_, RankCor = NA_real_))
  }
  
  # Define significant features based on adjusted p-value threshold
  sig_gt <- gt$feature[!is.na(gt$adj.P.Val) & gt$adj.P.Val < fdr_thresh]
  sig_imp <- res$feature[!is.na(res$adj.P.Val) & res$adj.P.Val < fdr_thresh]
  
  # Jaccard similarity of significant feature sets
  common <- intersect(gt$feature, res$feature)
  rank_cor <- if (length(common) > 2) {
    cor(gt$P.Value[match(common, gt$feature)],
        res$P.Value[match(common, res$feature)],
        method = "spearman", use = "complete.obs"
    )
  } else {
    NA_real_
  }
  
  list(result = res, Jaccard = .jaccard(sig_gt, sig_imp), RankCor = rank_cor)
}

#' ARI (vs GT labels + vs known Group) and silhouette width for imputed clustering.
#'
#' @param mat  Imputed matrix.
#' @param meta Metadata data frame.
#' @param gt   Ground-truth clustering object ($clustering from run_all_downstream()).
#' @return List: $result, $ARI, $ARI_Group, $Silhouette.
compare_clustering <- function(mat, meta, gt) {
  res <- run_clustering_downstream(mat, meta)
  ari_group <- if (!is.null(meta) && "Group" %in% colnames(meta)) {
    .ari(as.integer(factor(meta$Group)), res$labels)
  } else {
    NA_real_
  }
  list(
    result     = res,
    ARI        = .ari(gt$labels, res$labels),
    ARI_Group  = ari_group,
    Silhouette = res$silhouette_avg
  )
}

#' Edge-set Jaccard between imputed and GT co-abundance networks.
#'
#' Re-uses the exact threshold stored in the GT object so the comparison is fair.
#'
#' @param mat Imputed matrix.
#' @param gt  Ground-truth co-abundance object ($coabundance from run_all_downstream()).
#' @return List: $result, $EdgeJaccard.
compare_coabundance <- function(mat, gt) {
  res <- run_coabundance_downstream(mat, cor_threshold = gt$threshold %||% 0.7)
  edges_gt  <- paste(gt$edges$feat_a,  gt$edges$feat_b,  sep = "||")
  edges_imp <- paste(res$edges$feat_a, res$edges$feat_b, sep = "||")
  list(
    # Drop cor_matrix (features × features) — far too large to store per
    # method × prop × rep, and no plot reads it back.  Keep only the edge
    # tibble in case future plots need edge-level detail.
    result      = list(edges = res$edges, threshold = res$threshold),
    EdgeJaccard = .jaccard(edges_gt, edges_imp)
  )
}

#' Spearman correlation of per-feature variances (detects over-smoothing).
#'
#' @param mat  Imputed matrix.
#' @param gt   Ground-truth matrix (mat_true from run_all_downstream()).
#' @return Scalar, or NA when either variance vector is constant.
compare_feature_variance <- function(mat, gt) {
  v_gt <- apply(gt, 2, var, na.rm = TRUE)
  v_imp <- apply(mat, 2, var, na.rm = TRUE)
  if (sd(v_gt) == 0 || sd(v_imp) == 0) {
    return(NA_real_)
  }
  cor(v_gt, v_imp, method = "spearman", use = "complete.obs")
}

#' Compute mean Kolmogorov-Smirnov statistic across all features to check distribution shift.
#'
#' Both matrices are log1p-transformed before comparison so the scale matches
#' the rest of the pipeline.  Ties are broken by jittering (suppresses the
#' ks.test "cannot compute exact p-value with ties" warning that is expected
#' for integer/zero-inflated MS data and does not affect the D statistic).
#'
#' @param mat Imputed matrix (ROIs × features).
#' @param gt  Ground-truth matrix (ROIs × features).
#' @return Scalar — mean KS D-statistic across features (lower = better).
compare_feature_distributions <- function(mat, gt) {
  common_feats <- intersect(colnames(gt), colnames(mat))
  if (length(common_feats) == 0) return(NA_real_)
  
  gt_log  <- log1p(gt[,  common_feats, drop = FALSE])
  mat_log <- log1p(mat[, common_feats, drop = FALSE])
  
  ks_stats <- vapply(common_feats, function(feat) {
    x <- gt_log[, feat]
    y <- mat_log[, feat]
    # suppressWarnings: ties warning is cosmetic; D statistic is unaffected
    res <- tryCatch(
      suppressWarnings(ks.test(x, y)),
      error = function(e) NULL
    )
    if (is.null(res)) NA_real_ else res$statistic[[1L]]
  }, numeric(1))
  
  mean(ks_stats, na.rm = TRUE)
}

#' Compare PLS-DA subspace tracking vs Ground Truth
#'
#' Fits a PLS-DA model on the imputed matrix and measures how well its
#' discriminant subspace matches the ground-truth model via Procrustes SS.
#'
#' @param mat          Imputed matrix (ROIs × features).
#' @param meta         Metadata data frame.
#' @param gt           Ground-truth PLS-DA object ($plsda_gt from run_all_downstream()).
#' @param contrast_col Column in meta encoding the grouping factor.
#' @return List: $result (plsr model or NULL), $PLS_ProcrustesSS (scalar).
compare_plsda <- function(mat, meta, gt, contrast_col = "Group") {
  if (is.null(meta) || !contrast_col %in% colnames(meta) ||
      !requireNamespace("pls", quietly = TRUE)) {
    return(list(result = NULL, PLS_ProcrustesSS = NA_real_))
  }
  
  group_factor <- factor(meta[[contrast_col]])
  if (nlevels(group_factor) < 2 || is.null(gt$result)) {
    return(list(result = NULL, PLS_ProcrustesSS = NA_real_))
  }
  
  Y_dummy  <- model.matrix(~ group_factor - 1)
  X_scaled <- .preprocess(mat)
  
  # ncomp must be < nrow(X) and <= ncol(X); cap at 2 discriminant components
  ncomp_safe <- min(2L, nrow(X_scaled) - 1L, ncol(X_scaled))
  if (ncomp_safe < 1L) return(list(result = NULL, PLS_ProcrustesSS = NA_real_))
  
  mod_res <- tryCatch(
    pls::plsr(Y_dummy ~ X_scaled, ncomp = ncomp_safe, method = "oscorespls"),
    error = function(e) NULL
  )
  if (is.null(mod_res)) return(list(result = NULL, PLS_ProcrustesSS = NA_real_))
  
  # pls::scores() returns a 3-D array [obs × comp × responses] for plsr;
  # drop the third dimension to get the [obs × comp] score matrix.
  extract_scores <- function(model, n) {
    s <- pls::scores(model)
    if (length(dim(s)) == 3L) s <- s[, seq_len(n), 1L, drop = FALSE]
    else                       s <- s[, seq_len(n), drop = FALSE]
    s
  }
  
  n_comp_common <- min(ncomp_safe, dim(pls::scores(gt$result))[2L])
  gt_scores  <- extract_scores(gt$result, n_comp_common)
  imp_scores <- extract_scores(mod_res,   n_comp_common)
  
  ss <- tryCatch(
    vegan::procrustes(gt_scores, imp_scores, symmetric = TRUE)$ss,
    error = function(e) NA_real_
  )
  
  # Store only the score matrix — the full plsr model carries fitted values,
  # residuals, and the entire X matrix, making it very large.  plot_plsda_overlay
  # needs only the [obs × comp] scores to draw the scatter.
  imp_scores_df        <- as.data.frame(imp_scores)
  colnames(imp_scores_df) <- paste0("LV", seq_len(ncol(imp_scores_df)))
  
  list(result = imp_scores_df, PLS_ProcrustesSS = ss)
}


# ═════════════════════════════════════════════════════════════════════════════
# Module registry
# ═════════════════════════════════════════════════════════════════════════════
#
# Each entry is a named list describing one downstream analysis module.
#
# Fields
# ------
# store_key  chr   Key under which the full result list is stored in the output.
# gt_key     chr   Key used inside the ground_truth object (run_all_downstream()).
# compare_fn fn    compare_*(mat, [meta,] gt, ...) → list($result, $metric1, ...)
# metrics    chr   Names of the scalar metrics this module contributes to the
#                  tidy output table.  Must match names in compare_fn's return.
# needs_meta lgl   TRUE  → called as compare_fn(mat, meta, gt, ...)
#                  FALSE → called as compare_fn(mat, gt, ...)
# extra_args list  Static extra arguments forwarded to compare_fn (optional).
#
# ═════════════════════════════════════════════════════════════════════════════

DOWNSTREAM_MODULES <- list(
  pca = list(
    store_key  = "pca_store",
    gt_key     = "pca",
    compare_fn = compare_pca,
    metrics    = c("ProcrustesSS", "VarDiffPC1", "VarDiffPC2"),
    needs_meta = TRUE,
    extra_args = list()
  ),
  de_ttest = list(
    store_key  = NULL,   # full DE table not needed after scalars are extracted
    gt_key     = "de_t",
    compare_fn = compare_de,
    metrics    = c("TTest_Jaccard", "TTest_RankCor"),
    needs_meta = TRUE,
    extra_args = list(method = "t.test")
  ),
  de_wilcox = list(
    store_key  = NULL,   # full DE table not needed after scalars are extracted
    gt_key     = "de_w",
    compare_fn = compare_de,
    metrics    = c("Wilcox_Jaccard", "Wilcox_RankCor"),
    needs_meta = TRUE,
    extra_args = list(method = "wilcoxon")
  ),
  clustering = list(
    store_key  = "clust_store",
    gt_key     = "clustering",
    compare_fn = compare_clustering,
    metrics    = c("ARI", "ARI_Group", "Silhouette"),
    needs_meta = TRUE,
    extra_args = list()
  ),
  coabundance = list(
    store_key  = "coab_store",
    gt_key     = "coabundance",
    compare_fn = compare_coabundance,
    metrics    = "EdgeJaccard",
    needs_meta = FALSE,
    extra_args = list()
  ),
  feature_variance = list(
    store_key  = NULL,
    gt_key     = "mat_true",
    compare_fn = compare_feature_variance,
    metrics    = "FeatureVarCor",
    needs_meta = FALSE,
    extra_args = list()
  ),
  distribution_match = list(
    store_key  = NULL,
    gt_key     = "mat_true",
    compare_fn = compare_feature_distributions,
    metrics    = "MeanKSDistance",
    needs_meta = FALSE,
    extra_args = list()
  ),
  plsda = list(
    store_key  = "pls_store",
    gt_key     = "plsda_gt", # Will map to the run_all_downstream generation step
    compare_fn = compare_plsda,
    metrics    = "PLS_ProcrustesSS",
    needs_meta = TRUE,
    extra_args = list()
  )
)


# ═════════════════════════════════════════════════════════════════════════════
# Orchestration
# ═════════════════════════════════════════════════════════════════════════════

#' Run all downstream analyses on a single fully-observed matrix.
#'
#' Iterates DOWNSTREAM_MODULES to decide what to run
#'
#' @param mat           Numeric matrix (ROIs × features), no NAs.
#' @param meta          Metadata data frame, or NULL for spatial-only data.
#' @param contrast_col  Column for DE analysis.
#' @param de_levels     c(reference, test) for DE contrast, or NULL.
#' @param cor_threshold Edge threshold for the co-abundance network.
#' @return Named list: $mat_true, $pca, $de_t, $de_w,
#'                     $clustering, $coabundance, $nmf.
run_all_downstream <- function(mat, meta = NULL, contrast_col = "Group",
                               de_levels = NULL, cor_threshold = 0.7) {
  has_de <- !is.null(meta) && contrast_col %in% colnames(meta)
  
  message("  [downstream] PCA ...")
  pca <- run_pca_downstream(mat, meta)
  
  message(
    if (has_de) {
      paste0("  [downstream] DE (", contrast_col, ") ...")
    } else {
      "  [downstream] DE — skipped (no group labels in meta)"
    }
  )
  de_t <- run_de_downstream(mat, meta, contrast_col, de_levels, "t.test")
  de_w <- run_de_downstream(mat, meta, contrast_col, de_levels, "wilcoxon")
  
  message("  [downstream] Clustering ...")
  clust <- run_clustering_downstream(mat, meta)
  
  message("  [downstream] Co-abundance network ...")
  coab_full <- run_coabundance_downstream(mat, cor_threshold)
  # Drop the full correlation matrix from the stored GT object — it's
  # features × features and is only needed transiently inside compare_coabundance.
  coab <- list(
    edges     = coab_full$edges,
    threshold = coab_full$threshold
  )
  
  message("  [downstream] PLS-DA baseline structure ...")
  plsda_gt <- list(result = NULL)
  if (has_de && requireNamespace("pls", quietly = TRUE)) {
    group_factor <- factor(meta[[contrast_col]])
    if (nlevels(group_factor) >= 2) {
      Y_dummy  <- model.matrix(~ group_factor - 1)
      X_scaled <- .preprocess(mat)
      # ncomp must be < nrow(X) and <= ncol(X)
      ncomp_safe <- min(2L, nrow(X_scaled) - 1L, ncol(X_scaled))
      if (ncomp_safe >= 1L) {
        plsda_gt$result <- tryCatch(
          pls::plsr(Y_dummy ~ X_scaled, ncomp = ncomp_safe, method = "oscorespls"),
          error = function(e) NULL
        )
      }
    }
  }
  
  list(
    mat_true    = mat,
    pca         = pca,
    de_t        = de_t,
    de_w        = de_w,
    clustering  = clust,
    coabundance = coab,
    plsda_gt    = plsda_gt
  )
}


#' Compare every imputed matrix (all methods × props × replicates) against
#' the ground truth using all modules defined in DOWNSTREAM_MODULES.
compute_downstream_comparison <- function(results_storage, meta = NULL, ground_truth,
                                          contrast_col = "Group", de_levels = NULL) {
  stores <- Filter(
    Negate(is.null),
    setNames(lapply(DOWNSTREAM_MODULES, `[[`, "store_key"), names(DOWNSTREAM_MODULES))
  )
  store_env <- list2env(
    setNames(replicate(length(stores), list(), simplify = FALSE), unname(stores))
  )
  
  all_metrics <- unlist(lapply(DOWNSTREAM_MODULES, `[[`, "metrics"), use.names = FALSE)
  
  .na_row <- function(method, prop, rep_idx) {
    extra <- setNames(as.list(rep(NA_real_, length(all_metrics))), all_metrics)
    do.call(tibble, c(list(method = method, prop = prop, rep = rep_idx), extra))
  }
  
  metrics_df <- map_dfr(names(results_storage), function(rep_key) {
    rep_idx <- as.integer(str_extract(rep_key, "\\d+"))
    
    for (sk in ls(store_env)) store_env[[sk]][[rep_key]] <- list()
    
    map_dfr(names(results_storage[[rep_key]]), function(p_str) {
      prop <- as.numeric(p_str)
      imputed_list <- results_storage[[rep_key]][[p_str]]$imputed
      if (is.null(imputed_list) || length(imputed_list) == 0) {
        return(tibble())
      }
      
      for (sk in ls(store_env)) store_env[[sk]][[rep_key]][[p_str]] <- list()
      
      map_dfr(names(imputed_list), function(method_name) {
        mat_imp <- imputed_list[[method_name]]
        message("  Downstream [", rep_key, "][p=", p_str, "] – ", method_name)
        
        tryCatch(
          {
            module_results <- lapply(DOWNSTREAM_MODULES, function(mod) {
              gt_obj <- ground_truth[[mod$gt_key]]
              base <- if (isTRUE(mod$needs_meta)) {
                list(mat = mat_imp, meta = meta, gt = gt_obj)
              } else {
                list(mat = mat_imp, gt = gt_obj)
              }
              do.call(mod$compare_fn, c(base, mod$extra_args))
            })
            
            for (mod_name in names(DOWNSTREAM_MODULES)) {
              sk <- DOWNSTREAM_MODULES[[mod_name]]$store_key
              if (!is.null(sk)) {
                res_obj <- module_results[[mod_name]]
                store_env[[sk]][[rep_key]][[p_str]][[method_name]] <<- res_obj$result %||% res_obj
              }
            }
            
            # PATCHED: Maps returns to metric names if fields are missing
            scalars <- unlist(lapply(names(DOWNSTREAM_MODULES), function(mod_name) {
              mod <- DOWNSTREAM_MODULES[[mod_name]]
              res_obj <- module_results[[mod_name]]
              
              setNames(lapply(mod$metrics, function(m) {
                if (!is.list(res_obj)) {
                  return(res_obj)
                }
                if (m %in% names(res_obj)) {
                  return(res_obj[[m]])
                }
                
                # Check for un-prefixed fallback alternatives inside compare_de
                fallback <- str_remove(m, "^(TTest_|Wilcox_)")
                if (fallback %in% names(res_obj)) {
                  return(res_obj[[fallback]])
                }
                
                return(NA_real_)
              }), mod$metrics)
            }), recursive = FALSE)
            
            do.call(tibble, c(
              list(method = method_name, prop = prop, rep = rep_idx),
              scalars
            ))
          },
          error = function(e) {
            warning("Downstream failed for ", method_name, " [", rep_key, "][", p_str, "]: ", e$message)
            .na_row(method_name, prop, rep_idx)
          }
        )
      })
    })
  })
  
  c(list(metrics = metrics_df), as.list(store_env))
}