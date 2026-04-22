# =============================================================================
# 07_visualisation.R
# Purpose:  Visualise benchmark results: aggregated metric plots, runtime
#           scatter plots, and per-feature heatmaps.
# Inputs:   all_results (tidy data frame with cols method, prop, rep, metrics)
#           results_storage (nested list [[ rep_key ]][[ p_str ]])
#           mat_impute, coords (from main.R environment)
# Depends:  ggplot2, patchwork
# =============================================================================


# Functions ---------------------------------------------------------------

generate_msi_plots <- function(data_long, ncol = 5) {
  plot_list <- data_long %>%
    group_split(Feature) %>%
    map(~ {
      ggplot(.x, aes(X, Y, fill = NormIntensity)) +
        geom_tile() +
        scale_fill_viridis_c(option = "turbo", name = "Rel. Int.", na.value = "#FF00FF") +
        coord_fixed() +
        scale_y_reverse() +
        theme_minimal() +
        labs(title = unique(.x$Feature)) +
        theme(
          plot.title = element_text(size = 8, face = "bold"),
          axis.title = element_blank(),
          axis.text = element_text(size = 6),
          legend.key.height = unit(0.4, "cm")
        )
    })

  wrap_plots(plot_list, ncol = ncol)
}


# Generic function to plot any metric over missingness proportions

#' Plot a single metric aggregated across replicates.
#'
#' When multiple replicates are present (col `rep` in df), draws mean ± 1 SD
#' ribbon plus a mean line. Falls back to a plain line + points when only one
#' replicate exists (or for runtime, which only makes sense as a scatter).
#'
#' @param df        Tidy results data frame (all_results).
#' @param metric_name  Column name of the metric to plot.
#' @param y_label   Y-axis label; defaults to metric_name.
plot_metric <- function(df, metric_name, y_label = NULL) {
  if (is.null(y_label)) y_label <- metric_name
  
  has_replicates <- "rep" %in% colnames(df) && length(unique(df$rep)) > 1
  
  if (has_replicates) {
    # Aggregate across replicates
    agg <- df |>
      group_by(method, prop) |>
      summarise(
        mean_val = mean(.data[[metric_name]], na.rm = TRUE),
        sd_val   = sd(.data[[metric_name]],   na.rm = TRUE),
        .groups  = "drop"
      )
    
    ggplot(agg, aes(x = prop, y = mean_val, color = method, fill = method)) +
      geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
                  alpha = 0.15, color = NA) +
      geom_line(alpha = 0.8) +
      geom_point(size = 2) +
      theme_minimal() +
      labs(
        title    = paste("Benchmark:", metric_name),
        subtitle = paste0("Mean ± SD across ", length(unique(df$rep)), " replicates"),
        x        = "Proportion of Missing Values",
        y        = y_label,
        color    = "Method",
        fill     = "Method"
      ) +
      theme(legend.position = "right")
    
  } else {
    # Single replicate — plain lines
    ggplot(df, aes(x = prop, y = .data[[metric_name]], color = method)) +
      geom_line(alpha = 0.7) +
      geom_point(size = 2) +
      theme_minimal() +
      labs(
        title = paste("Benchmark:", metric_name),
        x     = "Proportion of Missing Values",
        y     = y_label,
        color = "Method"
      ) +
      theme(legend.position = "right")
  }
}


#' Accuracy vs. computational complexity scatter.
#' Uses mean NRMSE and mean runtime across replicates (if present).
plot_runtime_vs_nrmse <- function(df) {
  
  has_replicates <- "rep" %in% colnames(df) && length(unique(df$rep)) > 1
  
  plot_df <- if (has_replicates) {
    df |>
      group_by(method, prop) |>
      summarise(
        NRMSE       = mean(NRMSE,       na.rm = TRUE),
        runtime_sec = mean(runtime_sec, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    df
  }
  
  ggplot(plot_df, aes(x = runtime_sec, y = NRMSE, color = method, shape = factor(prop))) +
    geom_point(size = 3, alpha = 0.8) +
    scale_x_log10() +
    theme_minimal() +
    labs(
      title    = "Accuracy vs. Computational Complexity",
      subtitle = if (has_replicates) "Points show mean across replicates" else NULL,
      x        = "Runtime (seconds, log scale)",
      y        = "NRMSE (lower is better)",
      shape    = "Prop. missing"
    )
}


#' Spatial fidelity trade-off: SSIM vs. Moran's I difference.
plot_spatial_fidelity <- function(df) {
  
  has_replicates <- "rep" %in% colnames(df) && length(unique(df$rep)) > 1
  
  plot_df <- if (has_replicates) {
    df |>
      group_by(method, prop) |>
      summarise(
        MoranDiff = mean(MoranDiff, na.rm = TRUE),
        SSIM      = mean(SSIM,      na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    df
  }
  
  ggplot(plot_df, aes(x = MoranDiff, y = SSIM, color = method, shape = factor(prop))) +
    geom_point(size = 4, alpha = 0.8) +
    theme_minimal() +
    labs(
      title    = "Spatial Fidelity Trade-off",
      subtitle = if (has_replicates) "Points show mean across replicates" else NULL,
      x        = "Moran's I Difference (lower is better)",
      y        = "SSIM (higher is better)",
      shape    = "Prop. missing"
    )
}


#' Biochemical / spectral preservation: variance ratio vs. correlation structure.
plot_spectral_preservation <- function(df) {
  
  has_replicates <- "rep" %in% colnames(df) && length(unique(df$rep)) > 1
  
  plot_df <- if (has_replicates) {
    df |>
      group_by(method, prop) |>
      summarise(
        VarRatio  = mean(VarRatio,  na.rm = TRUE),
        CorStruct = mean(CorStruct, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    df
  }
  
  ggplot(plot_df, aes(x = VarRatio, y = CorStruct, color = method, shape = factor(prop))) +
    geom_point(size = 5, alpha = 0.7) +
    theme_minimal() +
    labs(
      title    = "Biochemical & Spectral Preservation",
      subtitle = if (has_replicates) "Points show mean across replicates" else NULL,
      x        = "Variance Ratio (target: 1.0)",
      y        = "Correlation Structure Preservation",
      shape    = "Prop. missing"
    )
}


#' Visualise per-feature heatmaps for one replicate.
#'
#' @param feature_idx    Column index of the feature in mat_impute.
#' @param mode           "all_methods" – one prop, all methods.
#'                       "all_props"   – one method, all props.
#' @param target_prop    Required when mode = "all_methods".
#' @param target_method  Required when mode = "all_props".
#' @param rep_idx        Which replicate to visualise (default: 1).
#' @param ncol           Number of columns in the facet grid.
#' @param results_storage, mat_impute, coords  Passed explicitly to avoid
#'   relying on global environment; fall back to globals if NULL.
visualise_heatmaps <- function(
    feature_idx,
    mode            = c("all_methods", "all_props"),
    target_prop     = NULL,
    target_method   = NULL,
    rep_idx         = 1,
    ncol            = 5,
    .results_storage = NULL,
    .mat_impute      = NULL,
    .coords          = NULL) {
  
  mode <- match.arg(mode)
  
  # Fall back to globals (kept for interactive convenience, but explicit is better)
  rs  <- .results_storage %||% results_storage
  mat <- .mat_impute      %||% mat_impute
  crd <- .coords          %||% coords
  
  rep_key      <- paste0("r", rep_idx)
  feature_name <- colnames(mat)[feature_idx]
  
  if (!(rep_key %in% names(rs))) {
    stop("Replicate '", rep_key, "' not found in results_storage. ",
         "Available: ", paste(names(rs), collapse = ", "))
  }
  
  plots_data <- list()
  plots_data[["Original"]] <- mat[, feature_idx]
  
  if (mode == "all_methods") {
    if (is.null(target_prop)) stop("target_prop is required for mode = 'all_methods'")
    p_str <- as.character(target_prop)
    if (!(p_str %in% names(rs[[rep_key]]))) {
      stop("Proportion '", p_str, "' not found in results_storage[['", rep_key, "']].")
    }
    plots_data[["Missing"]] <- rs[[rep_key]][[p_str]]$na_matrix[, feature_idx]
    for (m_name in names(rs[[rep_key]][[p_str]]$imputed)) {
      plots_data[[m_name]] <- rs[[rep_key]][[p_str]]$imputed[[m_name]][, feature_idx]
    }
    
  } else {
    if (is.null(target_method)) stop("target_method is required for mode = 'all_props'")
    for (p_str in names(rs[[rep_key]])) {
      plots_data[[paste0("Missing_", p_str)]] <- rs[[rep_key]][[p_str]]$na_matrix[, feature_idx]
      imp <- rs[[rep_key]][[p_str]]$imputed[[target_method]]
      if (!is.null(imp)) {
        plots_data[[paste0(target_method, "_", p_str)]] <- imp[, feature_idx]
      }
    }
  }
  
  # Global colour scale across all panels
  global_limits <- range(unlist(plots_data), na.rm = TRUE)
  
  make_base_plot <- function(val, title, na_magenta = FALSE) {
    df <- crd |> mutate(Intensity = val)
    ggplot(df, aes(X, Y, fill = Intensity)) +
      geom_tile() +
      coord_fixed() +
      scale_y_reverse() +
      theme_minimal() +
      theme(
        axis.text    = element_blank(),
        axis.title   = element_blank(),
        panel.grid   = element_blank()
      ) +
      labs(title = title) +
      scale_fill_viridis_c(
        option = "turbo",
        limits = global_limits,
        na.value = if (na_magenta) "#FF00FF" else "grey90"
      )
  }
  
  final_plots <- map2(plots_data, names(plots_data), function(val, name) {
    make_base_plot(val, name, na_magenta = grepl("Missing", name))
  })
  
  wrap_plots(final_plots, ncol = ncol) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = paste("Feature:", feature_name),
      subtitle = paste0("Replicate ", rep_idx,
                        " | Unified intensity scale across all panels"),
      caption  = "Magenta = missing values"
    )
}