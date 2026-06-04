# =============================================================================
# 09_downstream_visualisation.R
# Purpose:  Visualise downstream analysis comparisons across imputation methods.
# Inputs:   downstream_results (from compute_downstream_comparison()),
#           gt_downstream      (from run_all_downstream() on ground truth),
#           meta               (ROI metadata)
# Depends:  ggplot2, patchwork, RColorBrewer
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
                             colour_by = "Group") {
  rep_key <- paste0("r", rep_idx)
  p_str <- as.character(target_prop)
  pca_list <- downstream_results$pca_store[[rep_key]][[p_str]]

  if (is.null(pca_list) || length(pca_list) == 0) {
    stop("No PCA results for rep=", rep_idx, ", prop=", target_prop)
  }

  # Resolve colour column, falling back gracefully for spatial-only data
  gt_scores <- gt_downstream$pca$scores
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


# ── 5. Summary: all downstream metrics in one view ────────────────────────────

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
    "ProcrustesSS",
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