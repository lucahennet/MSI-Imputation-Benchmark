# =============================================================================
# 07_visualisation.R
# Purpose:  Visualise benchmark results: aggregated metric plots, runtime
#           scatter plots, and per-feature heatmaps.
# Inputs:   all_results (tidy data frame with cols method, prop, rep, metrics)
#           results_storage (nested list [[ rep_key ]][[ p_str ]])
#           mat_impute, coords (from main.R environment)
# Depends:  ggplot2, patchwork, tidyverse
# =============================================================================


# Functions ---------------------------------------------------------------

#' Generate per-feature heatmaps for all methods or all proportions.
#' 
#' Assumes data_long has columns Feature, X, Y, NormIntensity. Plots one heatmap per feature,
#' with a unified colour scale across all panels. Missing values (NA) are shown in
#' magenta to highlight them. The resulting plots are arranged in a grid with ncol columns.
#' @param data_long A long-format data frame with columns Feature, X, Y, NormIntensity.
#' @param ncol Number of columns in the facet grid (default: 5).
#' @return A patchwork object containing the arranged heatmaps.
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


#' Plot a single metric aggregated across replicates, with optional line or bar style.
#' 
#' When multiple replicates are present (col `rep` in df), draws mean ± 1 SD ribbon 
#' plus a mean line. Falls back to a plain line + points when only one replicate 
#' exists (or for runtime which only makes sense as a scatter). When type = "bar", 
#' shows bars with error bars instead of lines.
#' @param df Tidy results data frame (all_results).
#' @param metric_name Column name of the metric to plot.
#' @param y_label Y-axis label; defaults to metric_name.
#' @param type "line" (default) or "bar" to choose the plot style.
plot_metric <- function(df, metric_name, y_label = NULL, type = c("line", "bar")) {
  type <- match.arg(type)
  if (is.null(y_label)) y_label <- metric_name

  has_replicates <- "rep" %in% colnames(df) && length(unique(df$rep)) > 1
  n_reps <- length(unique(df$rep))
  multi_prop <- length(unique(df$prop)) > 1

  # Always aggregate — SD is 0 / NA for a single replicate
  agg <- df |>
    group_by(method, prop) |>
    summarise(
      mean_val = mean(.data[[metric_name]], na.rm = TRUE),
      sd_val   = if (has_replicates) sd(.data[[metric_name]], na.rm = TRUE) else 0,
      .groups  = "drop"
    )

  subtitle <- if (has_replicates) {
    paste0("Mean \u00b1 SD across ", n_reps, " replicates")
  } else {
    "Single replicate"
  }

  # ---- Line plot ----------------------------------------------------------
  if (type == "line") {
    p <- ggplot(agg, aes(x = prop, y = mean_val, color = method, fill = method)) +
      theme_minimal() +
      labs(
        title = paste("Benchmark:", metric_name),
        subtitle = subtitle,
        x = "Proportion of Missing Values",
        y = y_label,
        color = "Method", fill = "Method"
      ) +
      theme(legend.position = "right")

    if (has_replicates) {
      p <- p + geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
        alpha = 0.15, color = NA
      )
    }

    p + geom_line(alpha = 0.8) + geom_point(size = 2)

    # ---- Bar plot -----------------------------------------------------------
  } else {
    # Represent proportion as a factor label so bars are evenly spaced
    agg <- agg |>
      mutate(prop_label = paste0(prop * 100, "% missing"))

    # When only one proportion exists, skip the facet
    p <- ggplot(agg, aes(x = method, y = mean_val, fill = method)) +
      geom_col(width = 0.7, alpha = 0.85) +
      geom_errorbar(
        aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val),
        width = 0.25, linewidth = 0.6
      ) +
      theme_minimal() +
      labs(
        title    = paste("Benchmark:", metric_name),
        subtitle = subtitle,
        x        = NULL,
        y        = y_label,
        fill     = "Method"
      ) +
      theme(
        legend.position = "none", # colour is already on the x-axis
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.major.x = element_blank()
      )

    if (multi_prop) {
      p <- p + facet_wrap(~prop_label, scales = "fixed")
    }

    p
  }
}


#' Plot runtime vs. NRMSE, averaging across replicates and proportions to get one 
#' point per method.
plot_runtime_vs_nrmse <- function(df) {
  library(scales)
  
  # 1. Detect data structure for the subtitle
  n_reps  <- if("rep" %in% colnames(df)) length(unique(df$rep)) else 1
  n_props <- length(unique(df$prop))
  
  # Create a dynamic subtitle
  sub_text <- paste0(
    "Averaged over ", n_reps, " replicate(s) ",
    "across ", n_props, " missingness level(s) (", 
    paste(percent(unique(df$prop)), collapse = ", "), ")"
  )
  
  # 2. Process data: Average replicates, then average proportions
  plot_df <- df |>
    group_by(method, prop) |>
    summarise(
      m_nrmse   = mean(NRMSE, na.rm = TRUE),
      m_time    = mean(runtime_sec, na.rm = TRUE),
      .groups   = "drop"
    ) |>
    group_by(method) |>
    summarise(
      NRMSE        = mean(m_nrmse),
      runtime_sec  = mean(m_time),
      # Calculate range for "stability" visual
      nrmse_sd     = sd(m_nrmse), 
      .groups      = "drop"
    )
  
  # 3. Plotting
  ggplot(plot_df, aes(x = runtime_sec, y = NRMSE, color = method)) +
    # Add a horizontal "stability line" to show NRMSE variation across props
    geom_errorbar(aes(xmin = runtime_sec * 0.8, xmax = runtime_sec * 1.2), 
                   alpha = 0.3, width = 0) +
    geom_point(size = 5, alpha = 0.9) +
    geom_text(aes(label = method), vjust = -1.5, size = 3.5, fontface = "bold") +
    
    # Improved Log Scale: consistent spacing, readable time units
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 60, 300, 1200),
      labels = c("1ms", "10ms", "0.1s", "1s", "10s", "1m", "5m", "20m")
    ) +
    
    # Aesthetics
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(), # Remove messy minor lines
      # axis.title = element_text(face = "bold")
    ) +
    labs(
      title    = "Performance vs. Computational Cost",
      subtitle = sub_text,
      x        = "Average Runtime (Log Scale)",
      y        = "Mean NRMSE (Lower is Better)"
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