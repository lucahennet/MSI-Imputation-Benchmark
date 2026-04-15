# =============================================================================
# 07_visualisation.R
# Purpose:  Load, assemble and reshape MSI data from raw ROI files
# Inputs:   - A root directory containing one TXT subfolder and N CSV subfolders
# Outputs:  - raw_wide       : wide-format data frame (pixels × features)
#           - df_visualisation: long-format, TIC-normalised, for plotting
#           - df_impute      : clean wide matrix ready for imputation
#           - coords         : data frame of X/Y pixel coordinates
#           - mat_impute     : numeric matrix (pixels × features), no metadata
# Depends:  tidyverse, stringr
# =============================================================================


# Functions ---------------------------------------------------------------

generate_msi_plots <- function(data_long, ncol = 5) {
  plot_list <- data_long %>%
    group_split(Feature) %>%
    map(~{
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
plot_metric <- function(df, metric_name, y_label = NULL) {
  # Default label to metric name if not provided
  if (is.null(y_label)) y_label <- metric_name
  
  ggplot(df, aes(x = prop, y = !!sym(metric_name), color = method)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(
      title = paste("Benchmark:", metric_name),
      x = "Proportion of Missing Values",
      y = y_label,
      color = "Method"
    ) +
    theme(legend.position = "right")
}

# Complexity: Accuracy vs Time
plot_runtime_vs_nrmse <- function(df) {
  ggplot(df, aes(x = runtime_sec, y = NRMSE, color = method)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_x_log10() +
    theme_minimal() +
    labs(
      title = "Accuracy vs. Computational Complexity",
      x = "Runtime (seconds, log scale)",
      y = "NRMSE (Lower is better)"
    )
}

# Spatial Fidelity: Structure vs Spatial Correlation
plot_spatial_fidelity <- function(df) {
  ggplot(df, aes(x = MoranDiff, y = SSIM, color = method)) +
    geom_point(size = 4, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = "Spatial Fidelity Trade-off",
      x = "Moran's I Difference (Lower is better)",
      y = "SSIM (Higher is better)"
    )
}

# Spectral/Biochemical Validity
plot_spectral_preservation <- function(df) {
  ggplot(df, aes(x = VarRatio, y = CorStruct, color = method)) +
    geom_point(size = 5, alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Biochemical & Spectral Preservation",
      x = "Variance Ratio (Target: 1.0)",
      y = "Correlation Structure Preservation"
    )
}

visualise_heatmaps <- function(feature_idx, 
                               mode = c("all_methods", "all_props"), 
                               target_prop = NULL, 
                               target_method = NULL,
                               ncol = 5) {
  
  feature_name <- colnames(mat_impute)[feature_idx]
  plots_data <- list()
  
  # 1. Collect all data first to find the global min/max
  # Base: Original
  plots_data[["Original"]] <- mat_impute[, feature_idx]
  
  if (mode == "all_methods") {
    p_str <- as.character(target_prop)
    plots_data[["Missing"]] <- results_storage[[p_str]]$na_matrix[, feature_idx]
    for (m_name in names(results_storage[[p_str]]$imputed)) {
      plots_data[[m_name]] <- results_storage[[p_str]]$imputed[[m_name]][, feature_idx]
    }
  } else {
    for (p_str in names(results_storage)) {
      plots_data[[paste0("Missing_", p_str)]] <- results_storage[[p_str]]$na_matrix[, feature_idx]
      if (!is.null(results_storage[[p_str]]$imputed[[target_method]])) {
        plots_data[[paste0(target_method, "_", p_str)]] <- results_storage[[p_str]]$imputed[[target_method]][, feature_idx]
      }
    }
  }
  
  # 2. Calculate Global Limits (ignoring NAs)
  all_values <- unlist(plots_data)
  global_limits <- range(all_values, na.rm = TRUE)
  
  # 3. Create plots using the global limits
  make_base_plot <- function(val, title, na_magenta = FALSE) {
    df <- coords %>% mutate(Intensity = val)
    p <- ggplot(df, aes(X, Y, fill = Intensity)) +
      geom_tile() + 
      coord_fixed() + 
      scale_y_reverse() + 
      theme_minimal() +
      theme(axis.text = element_blank(), 
            axis.title = element_blank(), 
            panel.grid = element_blank()) +
      labs(title = title) +
      # Enforce global limits and guide
      scale_fill_viridis_c(option = "turbo", 
                           limits = global_limits, 
                           na.value = if(na_magenta) "#FF00FF" else "grey90")
    return(p)
  }
  
  # Build the plot list
  final_plots <- map2(plots_data, names(plots_data), function(val, name) {
    is_missing_plot <- grepl("Missing", name)
    make_base_plot(val, name, na_magenta = is_missing_plot)
  })
  
  # 4. Wrap with a collected guide
  wrap_plots(final_plots, ncol = ncol) + 
    plot_layout(guides = "collect") +
    plot_annotation(title = paste("Feature:", feature_name),
                    subtitle = "Unified intensity scale across all panels",
                    caption = "Magenta indicates missing values")
}