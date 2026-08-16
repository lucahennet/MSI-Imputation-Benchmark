library(tidyverse)


# ---- Aesthetics -----------------------------------------------------------

METHOD_COLORS <- c(
  "Median" = "black",
  "HM" = "#7570b3",
  "QRILC" = "#66a61e",
  "RF" = "#e6ab02",
  "KNN" = "#d95f02",
  "PPCA" = "#e7298a",
  "BPCA" = "#1b9e77",
  "spKNN" = "red",
  "GP" = "blue",
  "IDW" = "turquoise",
  "Oracle" = "purple"
)
METHOD_COLORS <- c(
  "Median" = "black",
  "HM" = "grey50",
  "QRILC" = "#1D9E75",
  "RF" = "#E24B4A",
  "KNN" = "#E07B39",
  "PPCA" = "#FF95CB",
  "BPCA" = "#EF9F27",
  "spKNN" = "#378ADD",
  "GP" = "#495BAF",
  "IDW" = "#B26CF7",
  "Oracle" = "#F1FD61"
)
METHOD_COLORS <- c(
  "Median" = "black",
  "HM" = "grey50",
  "QRILC" = "#1D9E75",
  "RF" = "#E24B4A",
  "KNN" = "#E07B39",
  "PPCA" = "#FF95CB",
  "BPCA" = "#EF9F27",
  "spKNN" = "#B26CF7",
  "GP" = "#495BAF",
  "IDW" = "#378ADD",
  "Oracle" = "#615CD6"
)
#"#BDF7F3"
#"#D69EC4"
#"#243159"

METHOD_ORDER <- names(METHOD_COLORS)

PROP_SHAPES <- c(
  "0.05" = 18,
  "0.1" = 16,
  "0.2" = 3,
  "0.3" = 17,
  "0.4" = 15,
  "0.5" = 8
)

theme_manuscript <- function() {
  theme_minimal() +
    theme(
      text = element_text(color = "black"),
      axis.title.x = element_text(size = 11, margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, margin = margin(r = 10)),
      axis.text = element_text(size = 9.5, color = "grey20"),
      axis.text.x = element_text(vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.position = "right"
    )
}


# ---- Shared helpers -------------------------------------------------------

# Combine a named list of result data frames, lock method order + factor levels.
.prep_df <- function(df_list, numeric_prop = TRUE) {
  df <- bind_rows(df_list, .id = "dataset") |>
    mutate(dataset = factor(dataset, levels = names(df_list)))
  df <- if (numeric_prop) {
    df |> mutate(prop = as.numeric(as.character(prop)))
  } else {
    df |> mutate(prop_factor = factor(as.character(prop)))
  }
  avail <- unique(df$method)
  ord   <- c(intersect(METHOD_ORDER, avail), setdiff(avail, METHOD_ORDER))
  df |> mutate(method = factor(method, levels = ord))
}

.has_reps <- function(df) "rep" %in% names(df) && length(unique(df$rep)) > 1

.export <- function(p, export_prefix, width, height) {
  if (!is.null(export_prefix)) {
    dir.create(dirname(export_prefix), recursive = TRUE, showWarnings = FALSE)
    ggsave(paste0(export_prefix, ".pdf"), p, width = width, height = height,
           device = cairo_pdf)
    ggsave(paste0(export_prefix, ".png"), p, width = width, height = height, dpi = 300)
  }
  invisible(p)
}


# ===========================================================================
# LINE METRICS  (NRMSE, SAM, RMSE, MAE, CCC, ...)
# Replaces the old plot_nrmse / plot_sam (which differed only in ylim).
# Works for main (2 datasets) and mix (1 dataset "MASH") without change.
# ===========================================================================

#' @param metric_name  column to plot (e.g. "NRMSE", "SAM", "RMSE")
#' @param ylim         y-axis limits, e.g. c(0, 0.8). NULL = auto.
#' @param subset       optional: filter to a `subset` value ("overall","MCAR",
#'                     "MNAR") when the results carry a subset column (mix data).
#' @param error        "sd" (default) or "se" for the ribbon.
plot_metric_line <- function(df_list,
                             metric_name,
                             y_label       = NULL,
                             ylim          = NULL,
                             subset        = NULL,
                             error         = c("sd", "se"),
                             method_colors = METHOD_COLORS,
                             export_prefix = NULL,
                             width = 8, height = 3.5) {
  
  error   <- match.arg(error)
  y_label <- y_label %||% metric_name
  
  df <- .prep_df(df_list, numeric_prop = TRUE)
  
  # Optional subset filter (mix-missingness results have overall/MCAR/MNAR)
  if (!is.null(subset) && "subset" %in% names(df)) {
    df <- df |> filter(subset == !!subset)
  }
  
  reps <- .has_reps(df)
  agg <- df |>
    group_by(dataset, method, prop) |>
    summarise(
      mean_val = mean(.data[[metric_name]], na.rm = TRUE),
      spread   = if (reps) {
        s <- sd(.data[[metric_name]], na.rm = TRUE)
        if (error == "se") s / sqrt(dplyr::n()) else s
      } else 0,
      .groups = "drop"
    ) |>
    mutate(spread = coalesce(spread, 0))
  
  p <- ggplot(agg, aes(prop, mean_val, color = method, fill = method)) +
    facet_wrap(~dataset, ncol = 2, scales = "fixed") +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    scale_x_continuous(breaks = unique(agg$prop),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Missing Proportion", y = y_label, color = "Method", fill = "Method") +
    theme_manuscript()
  
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  
  if (reps) {
    p <- p + geom_ribbon(aes(ymin = mean_val - spread, ymax = mean_val + spread),
                         alpha = 0.15, color = NA)
  }
  p <- p + geom_line(alpha = 0.8) + geom_point(size = 2)
  
  .export(p, export_prefix, width, height)
  p
}


# ===========================================================================
# SPATIAL FIDELITY  (SSIM vs Moran's I difference) — unchanged behaviour
# ===========================================================================
# 3.5 or 5 for the height
plot_spatial_fidelity <- function(df_list,
                                  x_label = expression(paste("Moran's ", italic("I"), " Difference")),
                                  y_label = "SSIM",
                                  xlim    = c(0, NA),
                                  ylim    = c(NA, 1.0),
                                  subset  = NULL,
                                  method_colors = METHOD_COLORS,
                                  export_prefix = NULL,
                                  width = 8, height = 5) {

  df <- .prep_df(df_list, numeric_prop = FALSE)
  if (!is.null(subset) && "subset" %in% names(df)) df <- df |> filter(subset == !!subset)

  agg <- df |>
    group_by(dataset, method, prop_factor) |>
    summarise(MoranDiff = mean(MoranDiff, na.rm = TRUE),
              SSIM      = mean(SSIM, na.rm = TRUE), .groups = "drop")

  p <- ggplot(agg, aes(MoranDiff, SSIM, color = method, shape = prop_factor)) +
    facet_wrap(~dataset, ncol = 2, scales = "fixed") +
    geom_point(size = 2.8, alpha = 0.85, stroke = 1) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_shape_manual(values = PROP_SHAPES, drop = FALSE) +
    scale_x_continuous(
      breaks = function(limits) unique(c(0, pretty(limits)[pretty(limits) >= 0]))
    ) +
    scale_y_continuous(
      breaks = function(limits) unique(c(pretty(limits)[pretty(limits) <= 1.0], 1.0))
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x = x_label, y = y_label, color = "Method", shape = "Missing Proportion") +
    theme_manuscript() +
    theme(panel.spacing = unit(1, "lines"))

  .export(p, export_prefix, width, height)
  p
}

# ===========================================================================
# SPECTRAL PRESERVATION  (Variance ratio vs Correlation structure) — unchanged
# ===========================================================================

plot_spectral_preservation <- function(df_list,
                                       x_label = "Variance Ratio (log scale)",
                                       y_label = expression(paste("", italic("R")[cor], "")),
                                       subset  = NULL,
                                       method_colors = METHOD_COLORS,
                                       export_prefix = NULL,
                                       width = 8, height = 5) {
  
  df <- .prep_df(df_list, numeric_prop = FALSE)
  if (!is.null(subset) && "subset" %in% names(df)) df <- df |> filter(subset == !!subset)
  
  agg <- df |>
    group_by(dataset, method, prop_factor) |>
    summarise(VarRatio  = mean(VarRatio, na.rm = TRUE),
              CorStruct = mean(CorStruct, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(agg, aes(VarRatio, CorStruct, color = method, shape = prop_factor)) +
    facet_wrap(~dataset, ncol = 2, scales = "fixed") +
    geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(size = 2.8, alpha = 0.85, stroke = 1) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_shape_manual(values = PROP_SHAPES, drop = FALSE) +
    scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1.0, 2.0, 3.0),
                  labels = c("0.1", "0.25", "0.5", "1.0", "2.0", "3.0")) +
    labs(x = x_label, y = y_label, color = "Method", shape = "Missing Proportion") +
    coord_cartesian(xlim = c(0.25, 3), ylim = c(NA, 1)) +
    theme_manuscript() +
    theme(panel.spacing = unit(1.5, "lines"))
  
  .export(p, export_prefix, width, height)
  p
}


# ===========================================================================
# ORACLE GAIN OVER BASELINE  (Bootstrap ribbon using nested resampling)
# ===========================================================================

plot_gain <- function(df_list,
                      base_method   = "RF",
                      target_method = "Oracle",
                      metric_name   = "RMSE",
                      subset_sel    = "overall",
                      n_boot        = 500,
                      seed          = 1,
                      y_label       = NULL,
                      export_prefix = NULL,
                      width = 8, height = 3.5) {
  
  y_label <- y_label %||% paste0(target_method, " gain over ", base_method, " (%)")
  
  df <- .prep_df(df_list, numeric_prop = TRUE)
  if ("subset" %in% names(df) && !is.null(subset_sel)) {
    df <- df |> filter(subset == !!subset_sel)
  }
  
  set.seed(seed)
  gain_color <- METHOD_COLORS[[target_method]] %||% "#1D9E75"
  
  res_df <- df |>
    group_by(dataset) |>
    group_modify(~ {
      d <- .x
      # Point estimate: ratio of means
      point <- d |>
        filter(method %in% c(base_method, target_method)) |>
        group_by(prop, method) |>
        summarise(m_val = mean(.data[[metric_name]], na.rm = TRUE), .groups = "drop") |>
        pivot_wider(names_from = method, values_from = m_val) |>
        mutate(gain = 100 * (.data[[base_method]] - .data[[target_method]]) / .data[[base_method]]) |>
        select(prop, gain)
      
      # Resample nested unit blocks (zones x reps) without left_join warnings
      has_units <- all(c("zone", "rep") %in% names(d))
      if (has_units) {
        nested_units <- d |>
          filter(method %in% c(base_method, target_method)) |>
          group_by(prop, zone, rep) |>
          nest() |>
          ungroup()
        
        boot <- map_dfr(seq_len(n_boot), function(b) {
          nested_units |>
            group_by(prop) |>
            slice_sample(prop = 1, replace = TRUE) |>
            unnest(cols = c(data)) |>
            group_by(prop, method) |>
            summarise(m_val = mean(.data[[metric_name]], na.rm = TRUE), .groups = "drop") |>
            pivot_wider(names_from = method, values_from = m_val) |>
            mutate(gain = 100 * (.data[[base_method]] - .data[[target_method]]) / .data[[base_method]], boot = b) |>
            select(prop, gain, boot)
        })
        band <- boot |>
          group_by(prop) |>
          summarise(lo = quantile(gain, 0.025, na.rm = TRUE),
                    hi = quantile(gain, 0.975, na.rm = TRUE), .groups = "drop")
        point <- left_join(point, band, by = "prop")
      } else {
        point <- point |> mutate(lo = gain, hi = gain)
      }
      point
    }) |>
    ungroup()
  
  p <- ggplot(res_df, aes(prop, gain)) +
    facet_wrap(~dataset, ncol = length(df_list), scales = "fixed") +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = gain_color) +
    geom_line(color = gain_color, linewidth = 0.8) +
    geom_point(color = gain_color, size = 2.5) +
    scale_x_continuous(breaks = unique(res_df$prop),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Missing Proportion", y = y_label) +
    theme_manuscript()
  
  .export(p, export_prefix, width, height)
  p
}


# ===========================================================================
# SUBSET BAR COMPARISON  (MCAR vs MNAR side-by-side, ordered by proportion)
# ===========================================================================

plot_subset_bar <- function(df_list,
                            metric_name   = "RMSE",
                            subsets       = c("MCAR", "MNAR"),
                            y_label       = NULL,
                            error         = c("se", "sd"),
                            method_colors = METHOD_COLORS,
                            export_prefix = NULL,
                            width = 8, height = 4) {
  
  error   <- match.arg(error)
  y_label <- y_label %||% metric_name
  
  df <- .prep_df(df_list, numeric_prop = FALSE)
  if (!("subset" %in% names(df))) {
    stop("Column 'subset' not found in dataset. Required for MCAR/MNAR bar comparison.")
  }
  
  df <- df |>
    filter(subset %in% subsets) |>
    mutate(subset = factor(subset, levels = subsets))
  
  reps <- .has_reps(df)
  
  # Compute metrics and enforce strict numerical ordering of proportions
  agg <- df |>
    mutate(prop_num = as.numeric(as.character(prop_factor))) |>
    group_by(dataset, prop_num, subset, method) |>
    summarise(
      mean_val = mean(.data[[metric_name]], na.rm = TRUE),
      spread   = if (reps) {
        s <- sd(.data[[metric_name]], na.rm = TRUE)
        if (error == "se") s / sqrt(dplyr::n()) else s
      } else 0,
      .groups = "drop"
    ) |>
    arrange(prop_num) |>
    mutate(
      spread = coalesce(spread, 0),
      prop_label = factor(
        paste0(prop_num * 100, "% missing"),
        levels = paste0(sort(unique(prop_num)) * 100, "% missing")
      )
    )
  
  p <- ggplot(agg, aes(x = method, y = mean_val, fill = method)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = mean_val - spread, ymax = mean_val + spread),
                  width = 0.25, linewidth = 0.5) +
    facet_grid(subset ~ prop_label, scales = "fixed") +
    scale_fill_manual(values = method_colors, drop = FALSE) +
    labs(x = NULL, y = y_label, fill = "Method") +
    theme_manuscript() +
    theme(
      axis.text.x = element_blank(),   # Removes x-axis labels under the bars
      axis.ticks.x = element_blank(),  # Removes x-axis tick marks
      panel.grid.major.x = element_blank()
    )
  
  .export(p, export_prefix, width, height)
  p
}


# ===========================================================================
# Back-compat shims so existing driver scripts keep working
# ===========================================================================
plot_nrmse <- function(df_list, metric_name = "NRMSE", y_label = "NRMSE",
                       ylim = c(0, 0.8), ...) {
  plot_metric_line(df_list, metric_name, y_label = y_label, ylim = ylim, ...)
}
plot_sam <- function(df_list, metric_name = "SAM", y_label = "SAM",
                     ylim = c(0, 0.3), ...) {
  plot_metric_line(df_list, metric_name, y_label = y_label, ylim = ylim, ...)
}


# ===========================================================================
# Histograms
# ===========================================================================

plot_zone_summary <- function(zone_df, patient = NULL,
                              export_prefix = NULL, width = 8, height = 11.5) {
  df <- if (!is.null(patient)) filter(zone_df, patient == !!patient) else zone_df
  y_max <- max(df$pct_missing, na.rm = TRUE) * 1.20
  
  p <- df |>
    mutate(zone_label = paste0(zone_id, "\n(", tissue_type, ")")) |>
    ggplot(aes(x = zone_label, y = pct_missing, fill = tissue_type)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.1f%%", pct_missing)),
              vjust = -0.4, size = 2.8, colour = "grey20") +
    scale_fill_manual(
      values = c(fib = "#E24B4A", norm = "#378ADD"),
      labels = c(fib = "Fibrosis", norm = "Normal"), name = "Tissue type"
    ) +
    scale_y_continuous(
      labels = \(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.02)), limits = c(0, y_max)
    ) +
    facet_grid(patient ~ mode) +
    labs(x = NULL, y = "Missing Proportion") +
    theme_manuscript() +
    theme(
      legend.position    = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size = 8)
    )
  
  if (!is.null(export_prefix)) {
    dir.create(dirname(export_prefix), recursive = TRUE, showWarnings = FALSE)
    ggsave(paste0(export_prefix, ".pdf"), p, width = width, height = height,
           device = cairo_pdf)
    ggsave(paste0(export_prefix, ".png"), p, width = width, height = height, dpi = 300)
  }
  p
}

plot_feature_classes <- function(class_df, patient = NULL,
                                 export_prefix = NULL, width = 8, height = 11.5) {
  df <- if (!is.null(patient)) filter(class_df, patient == !!patient) else class_df
  
  p <- df |>
    ggplot(aes(x = class, y = n, fill = class)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = paste0(n, "\n(", sprintf("%.1f%%", pct), ")")),
              vjust = -0.25, size = 2.6, colour = "grey20") +
    scale_fill_manual(values = CLASS_COLOURS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22)), labels = scales::comma) +
    facet_grid(patient ~ mode) +
    labs(x = NULL, y = "Number of features") +
    theme_manuscript() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(angle = 20, hjust = 1, vjust = 1, margin = margin(t = 2), size = 8)
    )
  
  if (!is.null(export_prefix)) {
    dir.create(dirname(export_prefix), recursive = TRUE, showWarnings = FALSE)
    ggsave(paste0(export_prefix, ".pdf"), p, width = width, height = height,
           device = cairo_pdf)
    ggsave(paste0(export_prefix, ".png"), p, width = width, height = height, dpi = 300)
  }
  p
}


# ===========================================================================
# Approach A plots
# ===========================================================================

# Approach A combined figure: scatter (A) + intensity-bin bars (B) side by side.
# Styled with theme_manuscript() to match the other manuscript figures.
# Requires: fim_all, threshold_table (from 10.1), patchwork, theme_manuscript().

# Manuscript palette for the two modes (from METHOD_COLORS family)
MODE_COLORS <- c(NEG = "#FF95CB", POS = "#68E29D")   # amber / orange

# ---- Panel A: missingness vs mean intensity (scatter + LOESS) -------------
pA1 <- fim_all |>
  ggplot(aes(x = mean_obs_int, y = 100 * prop_missing)) +
  geom_jitter(alpha = 0.22, size = 0.9, height = 0.15, width = 0,
              colour = "#378ADD") +
  geom_smooth(method = "loess", se = TRUE, colour = "#E24B4A",
              fill = "#E24B4A", alpha = 0.15, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 1, y = Inf, label = "intensity = 1",
           hjust = -0.08, vjust = 1.5, size = 2.8, colour = "grey40") +
  scale_x_log10() +
  coord_cartesian(ylim = c(0, NA)) +
  facet_wrap(~mode) +
  labs(x = "Mean observed intensity (log scale)", y = "Missing Proportion (%)") +
  theme_manuscript() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 2)),   # pull y-title closer
    axis.title.x = element_text(margin = margin(t = 2))    # (optional) same for x
  )

# ---- Panel B: mean missingness by intensity bin ---------------------------
pA2 <- threshold_table |>
  ggplot(aes(x = int_bin, y = mean_missingness, fill = mode)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.6) +
  geom_text(aes(label = n_features),
            position = position_dodge(width = 0.85),
            angle = 90, hjust = -0.2, vjust = 0.5, size = 2.4, colour = "grey20") +
  scale_fill_manual(values = MODE_COLORS, name = "Mode") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Mean observed intensity bin", y = "Mean missingness (%)") +
  theme_manuscript() +
  theme(
    legend.position = c(0.98, 0.98),        # inside plot, top-right
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
    legend.key.size = unit(0.35, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.margin = margin(t = 1, r = 2, b = 1, l = 2),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(margin = margin(r = 2)),
    axis.title.x = element_text(margin = margin(t = 2))
  )

# ---- Combine ---------------------------------------------------------------
# tag_levels adds (A) / (B) labels in the manuscript style.
pA_combined <- (pA1 | pA2) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))

# ---- Export ----------------------------------------------------------------
save_approachA <- function(prefix = "figures/approachA_combined",
                           width = 8, height = 5) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(prefix, ".pdf"), pA_combined, width = width, height = height,
         device = cairo_pdf)
  ggsave(paste0(prefix, ".png"), pA_combined, width = width, height = height, dpi = 300)
}


# Runtime vs. NRMSE -------------------------------------------------------

theme_manuscript <- function() {
  theme_minimal() +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(size = 14, face = "bold", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10.5, color = "grey30", margin = margin(b = 8)),
      axis.title.x = element_text(size = 11, margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, margin = margin(r = 10)),
      axis.text = element_text(size = 9.5, color = "grey20"),
      axis.text.x = element_text(vjust = 0.5),
      axis.text.y = element_text(hjust = 1),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.position = "right"
    )
}

#' Plot runtime vs. NRMSE, averaged across replicates and proportions to get
#' one point per method. Styled to match the other manuscript figures
#' (theme_manuscript(), METHOD_COLORS/METHOD_ORDER).
#'
#' @param df            combined data frame with columns: method, prop,
#'                      NRMSE, runtime_sec, dataset (e.g. "Zebrafish"/"Mouse"),
#'                      and optionally rep, mechanism.
#' @param mechanism     optional: filter to one mechanism (e.g. "MCAR") if
#'                      df contains multiple. NULL = pool across all present.
#' @param dataset_labels named vector mapping raw dataset values to facet
#'                      strip labels, e.g. c(Zebrafish = "Zebrafish Dataset").
#' @param error         "sd" (default) or "se" for the error bars, both axes
#' @param method_colors named colour vector, defaults to METHOD_COLORS
#' @param export_prefix if given, saves .pdf + .png via .export()
plot_runtime_vs_nrmse <- function(df,
                                  mechanism      = NULL,
                                  dataset_labels = c(Zebrafish = "Zebrafish Dataset",
                                                     Mouse     = "Mouse Dataset"),
                                  error          = c("sd", "se"),
                                  method_colors  = METHOD_COLORS,
                                  export_prefix  = NULL,
                                  width = 8, height = 3.5) {
  error <- match.arg(error)
  
  if (!is.null(mechanism) && "mechanism" %in% names(df)) {
    df <- df |> filter(mechanism == !!mechanism)
  }
  
  # Lock dataset facet order (Zebrafish left, Mouse right) from first
  # appearance in the data, same idea as .prep_df()'s dataset factor
  df <- df |> mutate(dataset = factor(dataset, levels = unique(dataset)))
  
  # Lock method order/colour, same convention as .prep_df() elsewhere
  avail <- unique(df$method)
  ord   <- c(intersect(METHOD_ORDER, avail), setdiff(avail, METHOD_ORDER))
  df    <- df |> mutate(method = factor(method, levels = ord))
  
  n_reps  <- if ("rep" %in% colnames(df)) length(unique(df$rep)) else 1
  n_props <- length(unique(df$prop))
  
  group_vars <- c("dataset", "method")
  
  # Step 1: average over reps within each (dataset, method, prop)
  per_prop <- df |>
    group_by(across(all_of(c(group_vars, "prop")))) |>
    summarise(
      m_nrmse = mean(NRMSE, na.rm = TRUE),
      m_time  = mean(runtime_sec, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Step 2: average across props; sd/se of those per-prop means = error bars
  plot_df <- per_prop |>
    group_by(across(all_of(group_vars))) |>
    summarise(
      NRMSE       = mean(m_nrmse),
      runtime_sec = mean(m_time),
      nrmse_sd    = sd(m_nrmse),
      time_sd     = sd(m_time),
      n           = n(),
      .groups     = "drop"
    ) |>
    mutate(
      nrmse_err = if (error == "se") nrmse_sd / sqrt(n) else nrmse_sd,
      time_err  = if (error == "se") time_sd  / sqrt(n) else time_sd
    )
  
  p <- ggplot(plot_df, aes(x = runtime_sec, y = NRMSE, color = method)) +
    facet_wrap(~dataset, ncol = 2, scales = "fixed",
               labeller = as_labeller(dataset_labels)) +
    geom_errorbar(aes(ymin = pmax(NRMSE - nrmse_err, 0), ymax = NRMSE + nrmse_err),
                  width = 0, alpha = 0.35, linewidth = 0.5) +
    geom_errorbar(aes(xmin = pmax(runtime_sec - time_err, 1e-4),
                      xmax = runtime_sec + time_err),
                  width = 0, alpha = 0.35, linewidth = 0.5,
                  orientation = "y") +
    geom_point(size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = method), size = 3, fontface = "bold",
                             show.legend = FALSE, seed = 42,
                             min.segment.length = 0.3) +
    scale_x_log10(
      breaks = c(0.001, 0.01, 0.1, 1, 10, 60, 300, 1200),
      labels = c("1ms", "10ms", "0.1s", "1s", "10s", "1m", "5m", "20m")
    ) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    theme_manuscript() +
    theme(
      panel.spacing    = unit(1, "lines"),
      legend.position  = "none",
      panel.grid.minor = element_blank()
    ) +
    labs(
      x     = "Average Runtime (log scale)",
      y     = paste0("Mean NRMSE \u00b1 ", toupper(error))
    )
  
  .export(p, export_prefix, width, height)
  p
}
