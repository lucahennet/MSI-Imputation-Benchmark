# =============================================================================
# 09_downstream_visualisation.R
# Purpose:  Visualise downstream analysis comparisons across imputation methods.
# Inputs:   downstream_results (from compute_downstream_comparison()),
#           gt_downstream      (from run_all_downstream() on ground truth),
#           meta               (ROI metadata)
# Depends:  ggplot2, patchwork, RColorBrewer, ggridges, pls
# =============================================================================


# ── Shared helpers ────────────────────────────────────────────────────────────

#' Map a vector of method names to a stable colour palette.
#' Handles any number of methods by interpolating Set1 beyond 9 colours.
.method_colours <- function(methods) {
  unique_methods <- unique(methods)
  n <- length(unique_methods)
  pal <- if (n <= 9) {
    RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n)
  }
  setNames(pal, unique_methods)
}

.prop_label <- function(prop) paste0(as.numeric(prop) * 100, "% missing")

.base_theme <- function() {
  theme_minimal(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.4),
      legend.position = "right"
    )
}

#' Summarise a single metric column to mean ± SD per method × prop.
.summarise_metric <- function(df, metric) {
  df |>
    group_by(method, prop) |>
    summarise(
      mean_val = mean(.data[[metric]], na.rm = TRUE),
      sd_val   = sd(.data[[metric]], na.rm = TRUE),
      .groups  = "drop"
    )
}


# ── 1. PCA overlay ────────────────────────────────────────────────────────────

#' Faceted PCA score plots — imputed (colour) over ground truth (grey).
#'
#' Falls back to cluster labels or a plain "ROI" group when the requested
#' colour_by column is absent from the scores (e.g. spatial-only data).
#'
#' @param downstream_results Output of compute_downstream_comparison().
#' @param gt_downstream      Output of run_all_downstream(mat_impute, meta).
#' @param target_prop        Missingness proportion to display (e.g. 0.4).
#' @param rep_idx            Replicate index (default 1).
#' @param pc_x,pc_y          PCs to plot on each axis.
#' @param colour_by          Meta column to use for colour (default "Group").
plot_pca_overlay <- function(downstream_results,
                             gt_downstream,
                             target_prop = 0.4,
                             rep_idx = 1,
                             pc_x = "PC1",
                             pc_y = "PC2",
                             colour_by = "Group",
                             meta = NULL) {
  rep_key <- paste0("r", rep_idx)
  p_str <- as.character(target_prop)
  pca_list <- downstream_results$pca_store[[rep_key]][[p_str]]

  if (is.null(pca_list) || length(pca_list) == 0) {
    stop("No PCA results for rep=", rep_idx, ", prop=", target_prop)
  }

  # Resolve colour column.
  # Priority: (1) already in PCA scores, (2) supplied via meta, (3) fallbacks.
  gt_scores <- gt_downstream$pca$scores
  if (!colour_by %in% colnames(gt_scores) && !is.null(meta) && colour_by %in% colnames(meta)) {
    gt_scores[[colour_by]] <- meta[[colour_by]]
    pca_list <- lapply(pca_list, function(p) {
      p$scores[[colour_by]] <- meta[[colour_by]]
      p
    })
  }
  if (!colour_by %in% colnames(gt_scores)) {
    if (!is.null(gt_downstream$clustering$labels)) {
      cluster_col <- factor(paste0("Cluster_", gt_downstream$clustering$labels))
      gt_scores$SpatialCluster <- cluster_col
      pca_list <- lapply(pca_list, function(p) {
        p$scores$SpatialCluster <- cluster_col
        p
      })
      colour_by <- "SpatialCluster"
    } else {
      gt_scores$ROI <- "ROI"
      pca_list <- lapply(pca_list, function(p) {
        p$scores$ROI <- "ROI"
        p
      })
      colour_by <- "ROI"
    }
  }

  gt_xy <- gt_scores |> rename(GT_x = !!pc_x, GT_y = !!pc_y)

  method_scores <- map_dfr(names(pca_list), function(m) {
    pca_list[[m]]$scores |>
      select(all_of(c(pc_x, pc_y, colour_by))) |>
      mutate(method = m)
  })

  ve <- gt_downstream$pca$variance_explained
  xlab <- paste0(pc_x, " (", round(ve[[pc_x]] * 100, 1), "% var, GT)")
  ylab <- paste0(pc_y, " (", round(ve[[pc_y]] * 100, 1), "% var, GT)")

  ggplot(method_scores, aes(.data[[pc_x]], .data[[pc_y]], colour = .data[[colour_by]])) +
    geom_point(
      data = gt_xy,
      aes(x = GT_x, y = GT_y),
      colour = "grey80", size = 1.8, inherit.aes = FALSE
    ) +
    geom_point(size = 2, alpha = 0.8) +
    facet_wrap(~method) +
    scale_colour_viridis_d(option = "D") +
    .base_theme() +
    labs(
      title = paste0(
        "PCA score overlay  |  ", .prop_label(target_prop),
        "  |  rep ", rep_idx
      ),
      subtitle = "Grey = ground truth; coloured = imputed",
      x = xlab, y = ylab, colour = colour_by
    )
}

#' Bar chart of Procrustes SS across methods × proportions.
#' Lower = imputed PC space closer to ground truth.
plot_procrustes <- function(downstream_results) {
  df <- .summarise_metric(downstream_results$metrics, "ProcrustesSS")

  ggplot(df, aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(0, mean_val - sd_val), ymax = mean_val + sd_val),
      width = 0.25
    ) +
    facet_wrap(~prop, labeller = as_labeller(.prop_label)) +
    coord_flip() +
    scale_fill_manual(values = .method_colours(df$method)) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title = "Procrustes SS: PC space vs ground truth",
      subtitle = "Lower = imputed PC structure closer to true",
      x = NULL, y = "Procrustes SS"
    )
}


# ── 2. Differential analysis ──────────────────────────────────────────────────

#' Bar charts for t-test and Wilcoxon DE metrics (Jaccard hit-list overlap,
#' rank correlation of p-values) per method × prop.


plot_de_comparisons <- function(downstream_results) {
  req_cols <- c("TTest_Jaccard", "TTest_RankCor", "Wilcox_Jaccard", "Wilcox_RankCor")
  available_cols <- intersect(req_cols, colnames(downstream_results$metrics))

  if (length(available_cols) == 0 || all(is.na(downstream_results$metrics[, available_cols, drop = FALSE]))) {
    message("Skipping DE plot: all metrics are missing or NA (no group labels in meta).")
    return(invisible(NULL))
  }

  metrics_long <- downstream_results$metrics |>
    select(method, prop, rep, all_of(available_cols)) |>
    pivot_longer(all_of(available_cols), names_to = "metric", values_to = "value") |>
    group_by(method, prop, metric) |>
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val   = sd(value, na.rm = TRUE),
      .groups  = "drop"
    )

  ggplot(metrics_long, aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(0, mean_val - sd_val), ymax = pmin(1, mean_val + sd_val)), width = 0.25) +
    facet_grid(metric ~ prop, labeller = labeller(prop = as_labeller(.prop_label)), scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = .method_colours(metrics_long$method)) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title = "Classical statistical test tracking vs ground truth",
      subtitle = "Jaccard: hit-list overlap  |  RankCor: p-value rank correlation (higher = better)",
      x = NULL, y = "Value"
    )
}


# ── 3. Clustering ─────────────────────────────────────────────────────────────

#' Grouped bar chart of clustering metrics (ARI, ARI_Group, Silhouette).
#'
#' ARI columns are shown only when they are non-NA (i.e. group labels existed).
plot_clustering_metrics <- function(downstream_results) {
  available <- c("ARI", "ARI_Group", "Silhouette")
  available <- available[vapply(
    available, function(m) {
      m %in% colnames(downstream_results$metrics) &&
        any(!is.na(downstream_results$metrics[[m]]))
    },
    logical(1)
  )]

  if (length(available) == 0) {
    message("No clustering metrics to plot.")
    return(invisible(NULL))
  }

  metrics_long <- downstream_results$metrics |>
    select(method, prop, rep, all_of(available)) |>
    pivot_longer(all_of(available), names_to = "metric", values_to = "value") |>
    group_by(method, prop, metric) |>
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val   = sd(value, na.rm = TRUE),
      .groups  = "drop"
    )

  ggplot(metrics_long, aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(
      ymin = pmax(-1, mean_val - sd_val),
      ymax = pmin(1, mean_val + sd_val)
    ), width = 0.25) +
    facet_grid(metric ~ prop,
      labeller = labeller(prop = as_labeller(.prop_label)),
      scales = "free_y"
    ) +
    coord_flip() +
    scale_fill_manual(values = .method_colours(metrics_long$method)) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title = "Clustering metrics vs ground truth",
      subtitle = "ARI: vs GT labels  |  Silhouette: cluster compactness (higher = better)",
      x = NULL, y = "Value"
    )
}


# ── 4. Co-abundance network ───────────────────────────────────────────────────

#' Bar chart of EdgeJaccard (edge-set overlap) per method × prop.
plot_network_metrics <- function(downstream_results) {
  df <- .summarise_metric(downstream_results$metrics, "EdgeJaccard")

  ggplot(df, aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(
      ymin = pmax(0, mean_val - sd_val),
      ymax = pmin(1, mean_val + sd_val)
    ), width = 0.25) +
    facet_wrap(~prop, labeller = as_labeller(.prop_label)) +
    coord_flip() +
    scale_fill_manual(values = .method_colours(df$method)) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title = "Co-abundance network topology match",
      subtitle = "EdgeJaccard: edge overlap (higher = better)",
      x = NULL, y = "Jaccard index"
    )
}

#' Whole-dataset intensity distribution overlay: all features, all replicates.
#'
#' Pools every (ROI × feature) log1p-intensity value across *all* replicates
#' for a given missingness proportion into a single distribution per method,
#' then plots them as overlapping density curves (one curve = one method,
#' faceted by missingness proportion).  The ground-truth distribution is drawn
#' as a black reference line so deviations are immediately visible.
#'
#' Averaging across replicates is done implicitly by pooling: values from
#' rep 1, rep 2, … are concatenated into one vector per method × prop before
#' density estimation.  This is appropriate here because we are characterising
#' the marginal intensity distribution — a shape property — not a per-replicate
#' metric.  Individual replicates are not shown separately because the
#' replicate-to-replicate variation of the global distribution is negligible
#' compared to method-to-method differences.
#'
#' @param results_storage  The *original* results_storage object produced by
#'   the imputation loop: structure [[rep_key]][[p_str]]$imputed[[method]].
#'   Do NOT pass downstream_results here.
#' @param mat_true   Fully-observed ground-truth matrix (ROIs × features).
#' @param props      Missingness proportions to display.  NULL = all available.
plot_distribution_overlay <- function(results_storage, mat_true, props = NULL) {
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("Package 'ggridges' is required: install.packages('ggridges')")
  }

  rep_keys <- names(results_storage)
  if (length(rep_keys) == 0) stop("results_storage is empty.")

  # Determine proportions to show
  all_props <- names(results_storage[[rep_keys[1]]])
  props <- if (is.null(props)) all_props else intersect(as.character(props), all_props)
  if (length(props) == 0) stop("None of the requested props found in results_storage.")

  # Ground-truth vector (log1p of every value) — same reference for all props
  gt_vals <- log1p(as.vector(mat_true))
  gt_vals <- gt_vals[is.finite(gt_vals)]

  rows <- list()

  for (p_str in props) {
    # Collect method names from the first replicate that has data
    method_names <- NULL
    for (rk in rep_keys) {
      imp <- results_storage[[rk]][[p_str]]$imputed
      if (!is.null(imp) && length(imp) > 0) {
        method_names <- names(imp)
        break
      }
    }
    if (is.null(method_names)) next

    # Pool values across all replicates for each method
    for (m in method_names) {
      pooled <- unlist(lapply(rep_keys, function(rk) {
        mat <- results_storage[[rk]][[p_str]]$imputed[[m]]
        if (is.null(mat)) {
          return(NULL)
        }
        v <- log1p(as.vector(mat))
        v[is.finite(v)]
      }), use.names = FALSE)

      if (length(pooled) == 0) next
      rows[[length(rows) + 1L]] <- tibble(
        prop   = as.numeric(p_str),
        method = m,
        value  = pooled
      )
    }

    # Ground truth row (repeated per prop so it appears in every facet)
    rows[[length(rows) + 1L]] <- tibble(
      prop   = as.numeric(p_str),
      method = "Ground Truth",
      value  = gt_vals
    )
  }

  if (length(rows) == 0) stop("No data collected — check results_storage structure.")
  combined <- bind_rows(rows)

  # Ensure Ground Truth is always the topmost ridge
  method_order <- c("Ground Truth", sort(setdiff(unique(combined$method), "Ground Truth")))
  combined$method <- factor(combined$method, levels = rev(method_order))

  n_methods <- length(method_order) - 1L # exclude GT
  custom_palette <- c("Ground Truth" = "#111111", .method_colours(method_order[-1]))

  n_reps <- length(rep_keys)
  n_feat <- ncol(mat_true)
  subtitle_txt <- paste0(
    "All ", n_feat, " features × ", n_reps,
    if (n_reps == 1) " replicate" else " replicates (pooled)",
    " — log\u2081p(intensity)"
  )

  ggplot(combined, aes(x = value, y = method, fill = method, colour = method)) +
    ggridges::geom_density_ridges(
      alpha          = 0.45,
      scale          = 1.15,
      linewidth      = 0.6,
      rel_min_height = 0.005
    ) +
    facet_wrap(~prop, labeller = as_labeller(.prop_label)) +
    scale_fill_manual(values = custom_palette) +
    scale_colour_manual(values = custom_palette) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title    = "Global intensity distribution: imputed vs ground truth",
      subtitle = subtitle_txt,
      x        = "log\u2081p(Intensity)",
      y        = NULL
    )
}


# ── 5. PLS-DA overlay ─────────────────────────────────────────────────────────

#' Faceted PLS-DA score plots — imputed (colour) over ground truth (grey).
#'
#' Mirrors plot_pca_overlay in structure.  Scores are extracted from the
#' plsr models stored in downstream_results$pls_store.  The ground-truth
#' model comes from gt_downstream$plsda_gt$result.
#'
#' @param downstream_results Output of compute_downstream_comparison().
#' @param gt_downstream      Output of run_all_downstream().
#' @param target_prop        Missingness proportion to display (e.g. 0.4).
#' @param rep_idx            Replicate index (default 1).
#' @param colour_by          Column in the scores data frame used for colour.
#'                           Must have been joined into the PCA scores (i.e. a
#'                           meta column); falls back to "Group" or "ROI".
plot_plsda_overlay <- function(downstream_results,
                               gt_downstream,
                               target_prop = 0.4,
                               rep_idx = 1,
                               colour_by = "Group",
                               meta = NULL) {
  rep_key <- paste0("r", rep_idx)
  p_str <- as.character(target_prop)
  pls_list <- downstream_results$pls_store[[rep_key]][[p_str]]

  if (is.null(pls_list) || length(pls_list) == 0) {
    stop(
      "No PLS-DA results for rep=", rep_idx, ", prop=", target_prop,
      ". Check that the plsda module ran and that meta has the contrast column."
    )
  }

  # Ground-truth scores: extract from the GT plsr model (still a full model in gt_downstream)
  gt_mod <- gt_downstream$plsda_gt$result
  if (is.null(gt_mod)) {
    stop("Ground-truth PLS-DA model is NULL. Re-run run_all_downstream() with a valid contrast column.")
  }
  .pls_scores_2d <- function(mod) {
    s <- pls::scores(mod)
    if (length(dim(s)) == 3L) {
      s <- s[, seq_len(min(2L, dim(s)[2L])), 1L, drop = FALSE]
    } else {
      s <- s[, seq_len(min(2L, ncol(s))), drop = FALSE]
    }
    df <- as.data.frame(s)
    colnames(df) <- paste0("LV", seq_len(ncol(df)))
    df
  }
  gt_scores <- .pls_scores_2d(gt_mod)

  # Resolve colour column: check gt PCA scores first (meta was joined there),
  # then fall back to the explicit meta argument, then to a plain ROI index.
  gt_pca_scores <- gt_downstream$pca$scores
  if (colour_by %in% colnames(gt_pca_scores)) {
    gt_scores[[colour_by]] <- gt_pca_scores[[colour_by]]
  } else if (!is.null(meta) && colour_by %in% colnames(meta)) {
    gt_scores[[colour_by]] <- meta[[colour_by]]
  } else {
    colour_by <- "ROI"
    gt_scores$ROI <- seq_len(nrow(gt_scores))
  }

  # pls_store now holds pre-extracted score data frames (LV1, LV2), not plsr models
  method_scores <- map_dfr(names(pls_list), function(m) {
    sc <- pls_list[[m]]
    if (is.null(sc) || !is.data.frame(sc)) {
      return(tibble())
    }
    sc[[colour_by]] <- gt_scores[[colour_by]]
    sc$method <- m
    sc
  })

  if (nrow(method_scores) == 0) {
    message("All PLS-DA score matrices are NULL for this prop/rep — nothing to plot.")
    return(invisible(NULL))
  }

  gt_xy <- gt_scores |> rename(GT_x = LV1, GT_y = LV2)

  ggplot(method_scores, aes(LV1, LV2, colour = .data[[colour_by]])) +
    geom_point(
      data = gt_xy,
      aes(x = GT_x, y = GT_y),
      colour = "grey80",
      size = 1.8,
      inherit.aes = FALSE
    ) +
    geom_point(size = 2, alpha = 0.8) +
    facet_wrap(~method) +
    scale_colour_viridis_d(option = "D") +
    .base_theme() +
    labs(
      title = paste0(
        "PLS-DA score overlay  |  ", .prop_label(target_prop),
        "  |  rep ", rep_idx
      ),
      subtitle = "Grey = ground truth; coloured = imputed",
      x = "LV1",
      y = "LV2",
      colour = colour_by
    )
}


# ── 6. Summary: all downstream metrics in one view ────────────────────────────

#' Faceted bar chart showing all downstream scalar metrics for a single
#' missingness proportion. Useful as a single-page overview figure.
#'
#' Metrics that are entirely NA for this prop are silently dropped,
#' so the spatial pipeline (no DE) and the group pipeline both work unchanged.
#'
#' @param downstream_results Output of compute_downstream_comparison().
#' @param target_prop        Missingness proportion to display.
#' @param higher_is_better   Character vector of metric names where higher is
#'                           better; all others are assumed lower = better.
plot_downstream_summary <- function(downstream_results,
                                    target_prop = 0.4,
                                    higher_is_better = c(
                                      "FeatureVarCor", "NMF_SpatialCor",
                                      "ARI", "ARI_Group", "Silhouette",
                                      "TTest_Jaccard", "TTest_RankCor",
                                      "Wilcox_Jaccard", "Wilcox_RankCor",
                                      "EdgeJaccard"
                                    )) {
  candidate_metrics <- c(
    "ProcrustesSS", "PLS_ProcrustesSS", "MeanKSDistance",
    "FeatureVarCor", "NMF_SpatialCor",
    "ARI", "ARI_Group", "Silhouette",
    "TTest_Jaccard", "TTest_RankCor",
    "Wilcox_Jaccard", "Wilcox_RankCor",
    "EdgeJaccard"
  )

  # Keep only metrics present in the data and not entirely NA for this prop
  sub <- downstream_results$metrics |> filter(prop == target_prop)
  available <- intersect(candidate_metrics, colnames(sub))
  not_all_na <- available[vapply(available, function(m) {
    any(!is.na(sub[[m]]))
  }, logical(1))]

  if (length(not_all_na) == 0) {
    message("plot_downstream_summary: no non-NA metrics for prop=", target_prop)
    return(invisible(NULL))
  }

  df <- sub |>
    select(method, rep, all_of(not_all_na)) |>
    pivot_longer(all_of(not_all_na), names_to = "metric", values_to = "value") |>
    mutate(direction = ifelse(metric %in% higher_is_better,
      "higher = better", "lower = better"
    )) |>
    group_by(method, metric, direction) |>
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      sd_val   = sd(value, na.rm = TRUE),
      .groups  = "drop"
    )

  ggplot(df, aes(x = method, y = mean_val, fill = method)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(
      ymin = pmax(0, mean_val - sd_val),
      ymax = mean_val + sd_val
    ), width = 0.25) +
    facet_wrap(~metric, scales = "free_y", ncol = 4) +
    coord_flip() +
    scale_fill_manual(values = .method_colours(df$method)) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title = paste0("Downstream analysis summary  |  ", .prop_label(target_prop)),
      subtitle = "All available metrics vs ground truth (mean ± SD across replicates)",
      x = NULL, y = "Mean value"
    )
}