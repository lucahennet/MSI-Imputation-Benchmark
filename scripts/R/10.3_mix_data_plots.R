# =============================================================================
# 10.3_mix_missingness_plots.R
# Purpose:  Plotting functions for the mixed missingness benchmark results
# Inputs:   - 
# Outputs:  - 
# Depends:  tidyverse, patchwork, scales
# =============================================================================


# ---- Paired significance test (methods share the same masked matrix) -------
# Wilcoxon signed-rank on per-unit differences; reports win-rate + median gap.
paired_test <- function(res, m1, m2, metric = "RMSE", subset_sel = "overall") {
  res |>
    filter(subset == subset_sel, method %in% c(m1, m2)) |>
    select(zone, prop, rep, method, value = all_of(metric)) |>
    pivot_wider(names_from = method, values_from = value) |>
    group_by(prop) |>
    summarise(
      n = n(),
      median_diff = median(.data[[m1]] - .data[[m2]], na.rm = TRUE),
      p_value = wilcox.test(.data[[m1]], .data[[m2]], paired = TRUE)$p.value,
      m1_wins_pct = 100 * mean(.data[[m1]] < .data[[m2]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      comparison = paste(m1, "vs", m2), metric = metric,
      subset = subset_sel, p_adj = p.adjust(p_value, method = "BH")
    )
}

#' Run the full set of comparisons that matter.
all_tests <- function(res, metric = "RMSE") {
  bind_rows(
    paired_test(res, "RF", "QRILC", metric, "overall"), # which global method wins
    paired_test(res, "Oracle", "RF", metric, "overall"), # value of matching vs RF
    paired_test(res, "Oracle", "QRILC", metric, "overall"), # value of matching vs QRILC
    paired_test(res, "QRILC", "RF", metric, "MNAR"), # crossover, MNAR side
    paired_test(res, "RF", "QRILC", metric, "MCAR") # crossover, MCAR side
  )
}


# ---- Metric vs proportion, mean +/- SE across reps x zones ------------------
plot_metric <- function(res, metric_name, subset = "overall",
                        y_label = NULL, error = c("se", "sd"),
                        type = c("line", "bar")) {
  error <- match.arg(error)
  type <- match.arg(type)
  if (is.null(y_label)) y_label <- metric_name
  df <- res |> filter(subset == !!subset)

  agg <- df |>
    group_by(method, prop) |>
    summarise(
      mean_val = mean(.data[[metric_name]], na.rm = TRUE),
      spread = if (error == "se") {
        sd(.data[[metric_name]], na.rm = TRUE) / sqrt(n())
      } else {
        sd(.data[[metric_name]], na.rm = TRUE)
      },
      .groups = "drop"
    )

  subtitle <- sprintf(
    "%s subset | mean +/- %s across %d rep x %d zones",
    subset, toupper(error),
    length(unique(df$rep)), length(unique(df$zone))
  )

  # ---- Line: metric vs proportion ----------------------------------------
  if (type == "line") {
    ggplot(agg, aes(prop, mean_val, colour = method, fill = method)) +
      geom_ribbon(aes(ymin = mean_val - spread, ymax = mean_val + spread),
        alpha = 0.15, colour = NA
      ) +
      geom_line() +
      geom_point(size = 2) +
      scale_x_continuous(labels = percent) +
      theme_minimal() +
      labs(
        title = paste("Benchmark:", metric_name), subtitle = subtitle,
        x = "Proportion of missing values", y = y_label,
        colour = "Method", fill = "Method"
      )

    # ---- Bar: method bars, faceted by proportion ---------------------------
  } else {
    agg |>
      mutate(prop_label = paste0(prop * 100, "% missing")) |>
      ggplot(aes(method, mean_val, fill = method)) +
      geom_col(width = 0.7, alpha = 0.85) +
      geom_errorbar(aes(ymin = mean_val - spread, ymax = mean_val + spread),
        width = 0.25, linewidth = 0.6
      ) +
      facet_wrap(~prop_label) +
      theme_minimal() +
      labs(
        title = paste("Benchmark:", metric_name), subtitle = subtitle,
        x = NULL, y = y_label, fill = "Method"
      ) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid.major.x = element_blank()
      )
  }
}

# ---- Oracle gain: value of perfect mechanism knowledge --------------------
# plot_gain <- function(res) {
#   res |>
#     filter(subset == "overall") |>
#     group_by(prop, rep, zone) |>
#     summarise(
#       best_global = min(RMSE[method %in% c("RF", "QRILC")]),
#       oracle = RMSE[method == "Oracle"], .groups = "drop"
#     ) |>
#     mutate(gain = 100 * (best_global - oracle) / best_global) |>
#     group_by(prop) |>
#     summarise(mean_gain = mean(gain), se = sd(gain) / sqrt(n()), .groups = "drop") |>
#     ggplot(aes(prop, mean_gain)) +
#     geom_ribbon(aes(ymin = mean_gain - se, ymax = mean_gain + se),
#       alpha = 0.15, fill = "#1D9E75"
#     ) +
#     geom_line(colour = "#1D9E75") +
#     geom_point(colour = "#1D9E75", size = 2.5) +
#     scale_x_continuous(labels = percent) +
#     theme_minimal() +
#     labs(
#       title = "Value of perfect mechanism knowledge (ceiling)",
#       subtitle = "Oracle RMSE gain over best single global method (mean +/- SE)",
#       x = "Proportion of missing values", y = "Oracle gain (%)"
#     )
# }

# Oracle gain vs RF — computed to MATCH the text (ratio of means, baseline RF).
# Aggregates RMSE across units FIRST, then takes the gain, so the plotted
# value equals 100*(mean_RF - mean_Oracle)/mean_RF at each prop.
#
# For an error band, we bootstrap over units (resample zones x reps) so the
# ribbon reflects uncertainty in the *ratio of means*, consistent with the
# point estimate. (You can't put an SE on a single ratio-of-means directly.)

plot_gain <- function(res, n_boot = 500, seed = 1) {
  ov <- res |> filter(subset == "overall")

  # point estimate: ratio of means, per prop
  point <- ov |>
    group_by(prop, method) |>
    summarise(RMSE = mean(RMSE), .groups = "drop") |>
    pivot_wider(names_from = method, values_from = RMSE) |>
    mutate(gain = 100 * (RF - Oracle) / RF) |>
    select(prop, gain)

  # bootstrap band: resample units within each prop
  set.seed(seed)
  units <- ov |> distinct(prop, zone, rep)
  boot <- map_dfr(seq_len(n_boot), function(b) {
    units |>
      group_by(prop) |>
      slice_sample(prop = 1, replace = TRUE) |>
      ungroup() |>
      left_join(ov, by = c("prop", "zone", "rep")) |>
      group_by(prop, method) |>
      summarise(RMSE = mean(RMSE), .groups = "drop") |>
      pivot_wider(names_from = method, values_from = RMSE) |>
      mutate(gain = 100 * (RF - Oracle) / RF, boot = b) |>
      select(prop, gain, boot)
  })

  band <- boot |>
    group_by(prop) |>
    summarise(
      lo = quantile(gain, 0.025),
      hi = quantile(gain, 0.975), .groups = "drop"
    )

  point |>
    left_join(band, by = "prop") |>
    ggplot(aes(prop, gain)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#1D9E75") +
    geom_line(colour = "#1D9E75", linewidth = 0.8) +
    geom_point(colour = "#1D9E75", size = 2.5) +
    scale_x_continuous(labels = scales::percent) +
    theme_minimal() +
    labs(
      title = "Oracle RMSE gain over RF",
      subtitle = "Ratio of mean RMSE (matches reported values); 95% bootstrap band",
      x = "Proportion of missing values", y = "Oracle gain over RF (%)"
    )
}

# ---- Spatial fidelity: SSIM vs Moran's I difference ------------------------
# Uses the SSIM / MoranDiff columns (present on the `overall` subset rows,
# where extra_metrics() attached them).
plot_spatial_fidelity <- function(res) {
  if (!all(c("SSIM", "MoranDiff") %in% names(res))) {
    stop(
      "SSIM / MoranDiff not in res. Re-run the benchmark with coords so ",
      "extra_metrics() is computed."
    )
  }
  res |>
    filter(subset == "overall") |>
    group_by(method, prop) |>
    summarise(
      SSIM = mean(SSIM, na.rm = TRUE),
      MoranDiff = mean(MoranDiff, na.rm = TRUE), .groups = "drop"
    ) |>
    ggplot(aes(MoranDiff, SSIM, colour = method, shape = factor(prop))) +
    geom_point(size = 4, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = "Spatial fidelity trade-off",
      subtitle = "Mean across replicates x zones (overall subset)",
      x = "Moran's I difference (lower is better)",
      y = "SSIM (higher is better)", shape = "Prop. missing"
    )
}


# ---- Spectral preservation: variance ratio vs correlation structure --------
plot_spectral_preservation <- function(res, subset = "overall") {
  res |>
    filter(subset == !!subset) |>
    group_by(method, prop) |>
    summarise(
      VarRatio = mean(VarRatio, na.rm = TRUE),
      CorStruct = mean(CorStruct, na.rm = TRUE), .groups = "drop"
    ) |>
    ggplot(aes(VarRatio, CorStruct, colour = method, shape = factor(prop))) +
    geom_point(size = 5, alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Biochemical & spectral preservation",
      subtitle = sprintf("%s subset | mean across replicates x zones", subset),
      x = "Variance ratio (target: 1.0)",
      y = "Correlation structure preservation", shape = "Prop. missing"
    )
}