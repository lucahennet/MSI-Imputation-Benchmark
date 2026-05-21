# =============================================================================
# 09_downstream_visualisation.R
# Purpose:  Visualise downstream analysis comparisons across imputation methods.
# Inputs:   downstream_results (from compute_downstream_comparison()),
#           gt_downstream (from run_all_downstream() on ground truth),
#           meta (ROI metadata)
# Depends:  ggplot2, patchwork, ggrepel
# =============================================================================


# Shared helpers ----------------------------------------------------------

# .method_colours <- function(methods) {
#   pal <- RColorBrewer::brewer.pal(max(3, length(methods)), "Set1")
#   setNames(pal[seq_along(methods)], methods)
# }

.method_colours <- function(methods) {
  unique_methods <- unique(methods)
  n <- length(unique_methods)
  if (n <= 9) {
    pal <- RColorBrewer::brewer.pal(max(3, n), "Set1")[1:n]
  } else {
    # Smoothly expands Set1 to handle 10+ methods without warnings
    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n)
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

# Summarise metrics_df → mean ± SD per method × prop
.summarise_metric <- function(df, metric) {
  df |>
    group_by(method, prop) |>
    summarise(
      mean_val = mean(.data[[metric]], na.rm = TRUE),
      sd_val   = sd(.data[[metric]], na.rm = TRUE),
      .groups  = "drop"
    )
}


# =============================================================================
# 1. PCA overlay
# =============================================================================

#' Faceted PCA score plots — imputed (colour) over ground truth (grey).
#'
#' @param downstream_results  Output of compute_downstream_comparison().
#' @param gt_downstream       Output of run_all_downstream(mat_impute, meta).
#' @param target_prop         Missingness proportion to display (e.g. 0.4).
#' @param rep_idx             Replicate index (default 1).
#' @param pc_x / pc_y         PCs to plot on each axis.
#' @param colour_by           Meta column to use for colour (default "Group").
plot_pca_overlay <- function(downstream_results,
                             gt_downstream,
                             target_prop = 0.4,
                             rep_idx     = 1,
                             pc_x        = "PC1",
                             pc_y        = "PC2",
                             colour_by   = "Group") {
  rep_key  <- paste0("r", rep_idx)
  p_str    <- as.character(target_prop)
  pca_list <- downstream_results$pca_store[[rep_key]][[p_str]]
  
  if (is.null(pca_list) || length(pca_list) == 0)
    stop("No PCA results for rep=", rep_idx, ", prop=", target_prop)
  
  # Dynamic Metadata Verification
  if (!colour_by %in% colnames(gt_downstream$pca$scores)) {
    if (!is.null(gt_downstream$clustering$labels)) {
      gt_downstream$pca$scores$SpatialCluster <- factor(paste0("Cluster_", gt_downstream$clustering$labels))
      for (m in names(pca_list)) {
        pca_list[[m]]$scores$SpatialCluster <- factor(paste0("Cluster_", gt_downstream$clustering$labels))
      }
      colour_by <- "SpatialCluster"
    } else {
      gt_downstream$pca$scores$ROI <- "ROI"
      for (m in names(pca_list)) {
        pca_list[[m]]$scores$ROI <- "ROI"
      }
      colour_by <- "ROI"
    }
  }
  
  gt_scores <- gt_downstream$pca$scores |>
    rename(GT_x = !!pc_x, GT_y = !!pc_y)
  
  method_scores <- map_dfr(names(pca_list), function(m) {
    pca_list[[m]]$scores |>
      select(all_of(c(pc_x, pc_y, colour_by))) |>
      mutate(method = m)
  })
  
  # Variance explained from ground truth
  ve   <- gt_downstream$pca$variance_explained
  xlab <- paste0(pc_x, " (", round(ve[[pc_x]] * 100, 1), "% var, GT)")
  ylab <- paste0(pc_y, " (", round(ve[[pc_y]] * 100, 1), "% var, GT)")
  
  ggplot(method_scores, aes(.data[[pc_x]], .data[[pc_y]], colour = .data[[colour_by]])) +
    geom_point(
      data         = gt_scores,
      aes(x = GT_x, y = GT_y),
      colour       = "grey80",
      size         = 1.8,
      inherit.aes  = FALSE
    ) +
    geom_point(size = 2, alpha = 0.8) +
    facet_wrap(~ method) +
    scale_colour_viridis_d(option = "D") +
    .base_theme() +
    labs(
      title    = paste0("PCA score overlay  |  ", .prop_label(target_prop),
                        "  |  rep ", rep_idx),
      subtitle = "Grey = ground truth; coloured = imputed",
      x        = xlab, y = ylab, colour = colour_by
    )
}

#' Bar chart of Procrustes SS across methods × proportions.
#' Lower = imputed PC space closer to ground truth.
plot_procrustes <- function(downstream_results) {
  df <- .summarise_metric(downstream_results$metrics, "ProcrustesSS")

  ggplot(df, aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(
      ymin = pmax(0, mean_val - sd_val),
      ymax = mean_val + sd_val
    ), width = 0.25) +
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


#' # =============================================================================
#' # 2. Differential analysis
#' # =============================================================================
#' 
#' #' Overlaid volcano plots — ground truth in grey, imputed in colour.
#' #'
#' #' @param downstream_results Output of compute_downstream_comparison().
#' #' @param gt_downstream      Output of run_all_downstream().
#' #' @param target_prop        Missingness proportion to display.
#' #' @param rep_idx            Replicate index.
#' #' @param fdr_thresh         FDR threshold for colouring significant features.
#' #' @param n_label            Number of top features to label per panel.
#' plot_volcano_overlap <- function(downstream_results,
#'                                  gt_downstream,
#'                                  target_prop = 0.4,
#'                                  rep_idx     = 1,
#'                                  fdr_thresh  = 0.05,
#'                                  n_label     = 10) {
#'   rep_key <- paste0("r", rep_idx)
#'   p_str   <- as.character(target_prop)
#'   de_list <- downstream_results$de_store[[rep_key]][[p_str]]
#'   
#'   if (is.null(de_list) || length(de_list) == 0)
#'     stop("No DE results for rep=", rep_idx, ", prop=", target_prop)
#'   
#'   gt_de <- gt_downstream$de |>
#'     mutate(neg_log_p = -log10(pmax(P.Value, 1e-300)),
#'            sig       = adj.P.Val < fdr_thresh)
#'   
#'   method_de <- map_dfr(names(de_list), function(m) {
#'     r <- de_list[[m]]
#'     if (nrow(r) == 0) return(tibble())
#'     r |>
#'       mutate(neg_log_p = -log10(pmax(P.Value, 1e-300)),
#'              sig       = adj.P.Val < fdr_thresh,
#'              method    = m)
#'   })
#'   
#'   # Top features to label (by ground truth significance)
#'   top_feats <- gt_de |>
#'     filter(sig) |>
#'     arrange(adj.P.Val) |>
#'     slice_head(n = n_label) |>
#'     pull(feature)
#'   
#'   ggplot(method_de, aes(logFC, neg_log_p)) +
#'     # GT background
#'     geom_point(
#'       data        = gt_de,
#'       aes(logFC, neg_log_p),
#'       colour      = "grey80", size = 0.8, inherit.aes = FALSE
#'     ) +
#'     geom_point(aes(colour = sig), size = 0.9, alpha = 0.7) +
#'     geom_hline(yintercept = -log10(fdr_thresh), linetype = "dashed",
#'                colour = "grey50", linewidth = 0.4) +
#'     ggrepel::geom_text_repel(
#'       data        = method_de |> filter(feature %in% top_feats & sig),
#'       aes(label   = feature),
#'       size        = 2.5, max.overlaps = 15, colour = "black"
#'     ) +
#'     scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "#e41a1c"),
#'                         labels = c("ns", paste0("FDR<", fdr_thresh))) +
#'     facet_wrap(~ method) +
#'     .base_theme() +
#'     labs(
#'       title    = paste0("Volcano overlap  |  ", .prop_label(target_prop),
#'                         "  |  rep ", rep_idx),
#'       subtitle = "Grey = ground truth; coloured = imputed",
#'       x        = "log2 FC", y = "-log10(p-value)",
#'       colour   = NULL
#'     )
#' }
#' 
#' #' Bar charts for DE scalar metrics (Recall, Jaccard, FC_Cor) per method × prop.
#' plot_del_metrics <- function(downstream_results) {
#'   metrics_long <- downstream_results$metrics |>
#'     select(method, prop, rep, DEL_Recall, DEL_Jaccard, FC_Cor) |>
#'     pivot_longer(c(DEL_Recall, DEL_Jaccard, FC_Cor),
#'                  names_to = "metric", values_to = "value") |>
#'     group_by(method, prop, metric) |>
#'     summarise(
#'       mean_val = mean(value, na.rm = TRUE),
#'       sd_val   = sd(value,   na.rm = TRUE),
#'       .groups  = "drop"
#'     )
#'   
#'   ggplot(metrics_long,
#'          aes(x = reorder(method, -mean_val), y = mean_val, fill = method)) +
#'     geom_col(width = 0.65, alpha = 0.85) +
#'     geom_errorbar(aes(ymin = pmax(0, mean_val - sd_val),
#'                       ymax = pmin(1, mean_val + sd_val)), width = 0.25) +
#'     facet_grid(metric ~ prop,
#'                labeller = labeller(prop = as_labeller(.prop_label)),
#'                scales   = "free_y") +
#'     coord_flip() +
#'     scale_fill_brewer(palette = "Set1") +
#'     .base_theme() +
#'     theme(legend.position = "none") +
#'     labs(
#'       title = "Differential analysis metrics vs ground truth",
#'       subtitle = "DEL_Recall: fraction of true DELs recovered  |  DEL_Jaccard: set overlap  |  FC_Cor: log-FC agreement",
#'       x = NULL, y = "Value"
#'     )
#' }

plot_de_comparisons <- function(downstream_results) {
  metrics_long <- downstream_results$metrics |>
    select(method, prop, rep, TTest_Jaccard, TTest_RankCor, Wilcox_Jaccard, Wilcox_RankCor) |>
    pivot_longer(c(TTest_Jaccard, TTest_RankCor, Wilcox_Jaccard, Wilcox_RankCor), names_to = "metric", values_to = "value") |>
    group_by(method, prop, metric) |>
    summarise(mean_val = mean(value, na.rm = TRUE), sd_val = sd(value, na.rm = TRUE), .groups = "drop")

  ggplot(metrics_long, aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(0, mean_val - sd_val), ymax = pmin(1, mean_val + sd_val)), width = 0.25) +
    facet_grid(metric ~ prop, labeller = labeller(prop = as_labeller(.prop_label)), scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = .method_colours(metrics_long$method)) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title = "Classical Statistical Test Tracking vs Ground Truth",
      subtitle = "Jaccard: Hit-list overlap  |  RankCor: P-value correlation (Higher = Better)",
      x = NULL, y = "Value"
    )
}


# =============================================================================
# 3. Clustering
# =============================================================================

#' #' Heatmap of ARI vs ground truth per method × prop.
#' #' Rows = methods, columns = missing proportions.
#' plot_ari_heatmap <- function(downstream_results) {
#'   df <- downstream_results$metrics |>
#'     group_by(method, prop) |>
#'     summarise(ARI = mean(ARI, na.rm = TRUE), .groups = "drop") |>
#'     mutate(prop_label = .prop_label(prop))
#'   
#'   ggplot(df, aes(x = prop_label, y = method, fill = ARI)) +
#'     geom_tile(colour = "white", linewidth = 0.5) +
#'     geom_text(aes(label = round(ARI, 2)), size = 3.5) +
#'     scale_fill_gradient2(low = "#d73027", mid = "#ffffbf", high = "#1a9850",
#'                          midpoint = 0.5, limits = c(0, 1),
#'                          name = "ARI") +
#'     .base_theme() +
#'     theme(panel.border = element_blank()) +
#'     labs(
#'       title    = "Clustering ARI vs ground truth",
#'       subtitle = "1 = identical assignment, 0 = random",
#'       x = "Missing proportion", y = NULL
#'     )
#' }

#' Grouped bar chart of clustering metrics (ARI, ARI_Group, Silhouette).
# plot_clustering_metrics <- function(downstream_results) {
#   metrics_long <- downstream_results$metrics |>
#     select(method, prop, rep, ARI, ARI_Group, Silhouette) |>
#     pivot_longer(c(ARI, ARI_Group, Silhouette),
#                  names_to = "metric", values_to = "value") |>
#     group_by(method, prop, metric) |>
#     summarise(
#       mean_val = mean(value, na.rm = TRUE),
#       sd_val   = sd(value,   na.rm = TRUE),
#       .groups  = "drop"
#     )
#   
#   ggplot(metrics_long,
#          aes(x = reorder(method, -mean_val), y = mean_val, fill = method)) +
#     geom_col(width = 0.65, alpha = 0.85) +
#     geom_errorbar(aes(ymin = pmax(-1, mean_val - sd_val),
#                       ymax = pmin(1,  mean_val + sd_val)), width = 0.25) +
#     facet_grid(metric ~ prop,
#                labeller = labeller(prop = as_labeller(.prop_label)),
#                scales   = "free_y") +
#     coord_flip() +
#     scale_fill_manual(values = .method_colours(metrics_long$method)) +
#     .base_theme() +
#     theme(legend.position = "none") +
#     labs(
#       title    = "Clustering metrics vs ground truth",
#       subtitle = "ARI: vs GT labels  |  ARI_Group: vs known groups  |  Silhouette: cluster compactness",
#       x = NULL, y = "Value"
#     )
# }
plot_clustering_metrics <- function(downstream_results) {
  df <- .summarise_metric(downstream_results$metrics, "Silhouette")
  
  ggplot(df, aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(-1, mean_val - sd_val), ymax = pmin(1, mean_val + sd_val)), width = 0.25) +
    facet_wrap(~prop, labeller = as_labeller(.prop_label)) +
    coord_flip() +
    scale_fill_manual(values = .method_colours(df$method)) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title = "Clustering: Mean Silhouette Width",
      subtitle = "Higher = Better, tighter and more distinct clusters",
      x = NULL, y = "Silhouette Width"
    )
}


# =============================================================================
# 4. Co-abundance / network
# =============================================================================

#' Bar charts of EdgeJaccard and CorrFrobenius per method × prop.
# plot_network_metrics <- function(downstream_results) {
#   metrics_long <- downstream_results$metrics |>
#     select(method, prop, rep, EdgeJaccard, CorrFrobenius) |>
#     pivot_longer(c(EdgeJaccard, CorrFrobenius),
#                  names_to = "metric", values_to = "value") |>
#     group_by(method, prop, metric) |>
#     summarise(
#       mean_val = mean(value, na.rm = TRUE),
#       sd_val   = sd(value,   na.rm = TRUE),
#       .groups  = "drop"
#     )
#   
#   ggplot(metrics_long,
#          aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
#     geom_col(width = 0.65, alpha = 0.85) +
#     geom_errorbar(aes(ymin = pmax(0, mean_val - sd_val),
#                       ymax = mean_val + sd_val), width = 0.25) +
#     facet_grid(metric ~ prop,
#                labeller = labeller(prop = as_labeller(.prop_label)),
#                scales   = "free_y") +
#     coord_flip() +
#     scale_fill_manual(values = .method_colours(metrics_long$method)) +
#     .base_theme() +
#     theme(legend.position = "none") +
#     labs(
#       title    = "Co-abundance network metrics vs ground truth",
#       subtitle = "EdgeJaccard: edge set overlap (higher = better)  |  CorrFrobenius: correlation matrix distance (lower = better)",
#       x = NULL, y = "Value"
#     )
# }

plot_network_metrics <- function(downstream_results) {
  df <- .summarise_metric(downstream_results$metrics, "EdgeJaccard")

  ggplot(df, aes(x = reorder(method, mean_val), y = mean_val, fill = method)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(aes(ymin = pmax(0, mean_val - sd_val), ymax = pmin(1, mean_val + sd_val)), width = 0.25) +
    facet_wrap(~prop, labeller = as_labeller(.prop_label)) +
    coord_flip() +
    scale_fill_manual(values = .method_colours(df$method)) +
    .base_theme() +
    theme(legend.position = "none") +
    labs(
      title = "Co-abundance Network Topology Match",
      subtitle = "EdgeJaccard: Edge Overlap (Higher = Better)",
      x = NULL, y = "Jaccard Index"
    )
}


# =============================================================================
# 5. Summary: all downstream metrics in one view
# =============================================================================

#' Faceted bar chart showing all downstream scalar metrics for a single prop.

#' Useful as a single-page overview figure.

#'

#' @param downstream_results Output of compute_downstream_comparison().

#' @param target_prop        Missingness proportion to display.

#' @param higher_is_better   Character vector of metric names where higher = better.

#'                           All others are assumed lower = better.

plot_downstream_summary <- function(downstream_results,
                                    target_prop = 0.4,
                                    higher_is_better = c(
                                      "TTest_Jaccard", "TTest_RankCor",
                                      "Wilcox_Jaccard", "Wilcox_RankCor",
                                      "Silhouette", "EdgeJaccard"
                                    )) {
  candidate_metrics <- c(
    "ProcrustesSS", "TTest_Jaccard", "TTest_RankCor",
    "Wilcox_Jaccard", "Wilcox_RankCor", "Silhouette", "EdgeJaccard"
  )


  # Keep only metrics present in the data and not entirely NA

  # (DE metrics are all-NA in the spatial pipeline → silently dropped)

  available <- intersect(candidate_metrics, colnames(downstream_results$metrics))

  not_all_na <- available[vapply(available, function(m) {
    any(!is.na(downstream_results$metrics[[m]]))
  }, logical(1))]


  if (length(not_all_na) == 0) {
    message("plot_downstream_summary: no non-NA metrics for prop=", target_prop)

    return(invisible(NULL))
  }


  df <- downstream_results$metrics |>
    filter(prop == target_prop) |>
    select(method, rep, all_of(not_all_na)) |>
    pivot_longer(all_of(not_all_na), names_to = "metric", values_to = "value") |>
    mutate(
      direction = ifelse(metric %in% higher_is_better, "higher = better", "lower = better")
    ) |>
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
      subtitle = "All verified metrics vs ground truth",
      x = NULL, y = "Mean value (across replicates)"
    )
}