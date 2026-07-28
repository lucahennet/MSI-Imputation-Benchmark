# Core functions ----------------------------------------------------------

load_mat <- function(data_path) {
  raw_wide <- assemble_msi_data(data_path) |> distinct(X, Y, .keep_all = TRUE)
  feature_cols <- colnames(raw_wide)[8:ncol(raw_wide)]
  df <- raw_wide |> select(X, Y, ROI, all_of(feature_cols))
  list(
    coords = df |> select(X, Y, ROI),
    mat = df |> select(-X, -Y, -ROI) |> as.matrix()
  )
}

annotate_zones <- function(coords, zone_meta) {
  coords |>
    left_join(zone_meta, by = join_by(ROI >= ROI_min, ROI <= ROI_max)) |>
    select(-ROI_min, -ROI_max)
}

build_long <- function(mat, coords, patient, mode) {
  as_tibble(mat) |>
    bind_cols(coords) |>
    pivot_longer(
      cols = -c(X, Y, ROI, zone_id, tissue_type),
      names_to = "feature", values_to = "intensity"
    ) |>
    mutate(missing = intensity == 0, patient = patient, mode = mode)
}

classify_feature <- function(prop) {
  case_when(
    prop == 0 ~ "No missing",
    prop < THRESH_LOW ~ "Low [<5%]",
    prop < THRESH_MODERATE ~ "Moderate [5-20%]",
    prop <= THRESH_HIGH ~ "High [20-50%]",
    TRUE ~ ">50%"
  )
}

# m/z from a feature column name ("X452.2813_.1" -> 452.2813)
feature_mz <- function(feature_names) {
  as.numeric(str_extract(feature_names, "(?<=^X)\\d+\\.?\\d*"))
}


# Summary functions -------------------------------------------------------

summarise_dataset <- function(long_df, patient, mode) {
  tibble(
    patient = patient, mode = mode,
    n_pixels = n_distinct(paste(long_df$X, long_df$Y)),
    n_features = n_distinct(long_df$feature),
    n_values = nrow(long_df),
    n_missing = sum(long_df$missing),
    pct_missing = 100 * mean(long_df$missing),
    n_features_w_missing = long_df |>
      group_by(feature) |> summarise(a = any(missing)) |> pull(a) |> sum(),
    pct_features_missing = 100 * (long_df |>
                                    group_by(feature) |> summarise(a = any(missing)) |> pull(a) |> mean()),
    n_pixels_w_missing = long_df |>
      group_by(X, Y) |> summarise(a = any(missing), .groups = "drop") |> pull(a) |> sum(),
    pct_pixels_missing = 100 * (long_df |>
                                  group_by(X, Y) |> summarise(a = any(missing), .groups = "drop") |> pull(a) |> mean())
  )
}

summarise_zones <- function(long_df, patient, mode) {
  long_df |>
    filter(!is.na(zone_id)) |>
    group_by(zone_id, tissue_type) |>
    summarise(
      n_pixels = n_distinct(paste(X, Y)), n_features = n_distinct(feature),
      n_missing = sum(missing), pct_missing = 100 * mean(missing),
      n_features_any_zero = n_distinct(feature[missing]),
      pct_features_zero = 100 * n_distinct(feature[missing]) / n_distinct(feature),
      .groups = "drop"
    ) |>
    mutate(patient = patient, mode = mode)
}

summarise_features <- function(long_df, patient, mode) {
  long_df |>
    group_by(feature) |>
    summarise(
      n_pixels = n(), n_missing = sum(missing), n_observed = n() - sum(missing),
      prop_missing = mean(missing), pct_missing = 100 * mean(missing),
      .groups = "drop"
    ) |>
    mutate(
      patient = patient, mode = mode,
      class = factor(classify_feature(prop_missing), levels = CLASS_LEVELS)
    ) |>
    arrange(desc(prop_missing))
}


# Terminal output ---------------------------------------------------------

print_freq_table <- function(feat_df, patient, mode) {
  bin_labels <- paste0(seq(0, 95, BIN_WIDTH), "-", seq(5, 100, BIN_WIDTH), "%")
  tbl <- feat_df |>
    mutate(bin = cut(pct_missing,
                     breaks = seq(0, 100, BIN_WIDTH),
                     include.lowest = TRUE, right = FALSE, labels = bin_labels
    )) |>
    count(bin, .drop = FALSE) |>
    mutate(pct = 100 * n / nrow(feat_df)) |>
    filter(n > 0)
  message(sprintf("\n  Frequency table: %s %s  (%d features)", patient, mode, nrow(feat_df)))
  message(sprintf("  %-12s  %5s  %6s", "Range", "n", "%"))
  message(sprintf("  %-12s  %5s  %6s", "-----", "-", "------"))
  for (i in seq_len(nrow(tbl))) {
    message(sprintf(
      "  %-12s  %5d  %5.1f%%",
      as.character(tbl$bin[i]), tbl$n[i], tbl$pct[i]
    ))
  }
}

print_class_summary <- function(feat_df, patient, mode) {
  tbl <- feat_df |>
    count(class, .drop = FALSE) |>
    mutate(pct = 100 * n / sum(n))
  message(sprintf("\n  Class summary: %s %s", patient, mode))
  for (i in seq_len(nrow(tbl))) {
    message(sprintf(
      "  %-22s  %5d  (%4.1f%%)",
      as.character(tbl$class[i]), tbl$n[i], tbl$pct[i]
    ))
  }
}


# Plot functions ----------------------------------------------------------

plot_histogram <- function(feat_df, mode_label, patient, y_limit = NULL, save = FALSE) {
  y_lim <- y_limit %||% {
    feat_df |>
      mutate(bin = floor(pct_missing / BIN_WIDTH) * BIN_WIDTH) |>
      count(bin) |>
      pull(n) |>
      max() |>
      (\(x) ceiling(x * 1.12))()
  }
  p <- feat_df |>
    mutate(
      bin_left = floor(pct_missing / BIN_WIDTH) * BIN_WIDTH,
      bin_class = factor(classify_feature(bin_left / 100), levels = CLASS_LEVELS)
    ) |>
    ggplot(aes(x = pct_missing, fill = bin_class)) +
    geom_histogram(
      breaks = seq(0, 100, BIN_WIDTH), closed = "left",
      colour = "white", linewidth = 0.25
    ) +
    geom_vline(
      xintercept = c(THRESH_LOW, THRESH_MODERATE, THRESH_HIGH) * 100,
      linetype = "dashed", colour = "grey35", linewidth = 0.5
    ) +
    annotate("text",
             x = c(THRESH_LOW, THRESH_MODERATE, THRESH_HIGH) * 100 + 0.8,
             y = y_lim * 0.96, label = c("5%", "20%", "50%"),
             hjust = 0, size = 2.8, colour = "grey35"
    ) +
    scale_fill_manual(values = CLASS_COLOURS, drop = FALSE, name = NULL) +
    scale_x_continuous(
      limits = c(0, 100), breaks = seq(0, 100, 10),
      expand = expansion(mult = c(0, 0.01)), labels = \(x) paste0(x, "%")
    ) +
    scale_y_continuous(
      limits = c(0, y_lim), expand = expansion(mult = c(0, 0)),
      labels = scales::comma
    ) +
    labs(
      title = paste0(patient, " — ", mode_label),
      subtitle = paste0(
        nrow(feat_df), " features  |  ",
        sprintf("%.1f%%", 100 * mean(feat_df$prop_missing)),
        " overall missing rate"
      ),
      x = "Missing values per feature (%)", y = "Number of features"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
      panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(colour = "grey40", size = 9)
    ) +
    guides(fill = guide_legend(nrow = 1))
  if (save) {
    ggsave(file.path(OUT_DIR, paste0(patient, "_histogram_", mode_label, ".pdf")),
           p,
           width = 8, height = 4.5
    )
  }
  p
}

plot_zone_summary <- function(zone_df, patient = NULL) {
  df <- if (!is.null(patient)) filter(zone_df, patient == !!patient) else zone_df
  y_max <- max(df$pct_missing, na.rm = TRUE) * 1.20
  df |>
    mutate(zone_label = paste0(zone_id, "\n(", tissue_type, ")")) |>
    ggplot(aes(x = zone_label, y = pct_missing, fill = tissue_type)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.1f%%", pct_missing)), vjust = -0.4, size = 3) +
    scale_fill_manual(
      values = c(fib = "#E24B4A", norm = "#378ADD"),
      labels = c(fib = "Fibrosis", norm = "Normal"), name = NULL
    ) +
    scale_y_continuous(
      limits = c(0, y_max), labels = \(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0))
    ) +
    facet_grid(patient ~ mode) +
    labs(
      title = if (!is.null(patient)) {
        paste0("Zone missingness — ", patient)
      } else {
        "Zone missingness — all patients"
      },
      subtitle = "Each bar = one 5x5 zone.", x = NULL, y = "Missing values (%)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom", panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(colour = "grey40", size = 9)
    )
}

plot_feature_classes <- function(class_df, patient = NULL) {
  df <- if (!is.null(patient)) filter(class_df, patient == !!patient) else class_df
  df |>
    ggplot(aes(x = class, y = n, fill = class)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    geom_text(aes(label = paste0(n, "\n(", sprintf("%.1f%%", pct), ")")),
              vjust = -0.25, size = 3
    ) +
    scale_fill_manual(values = CLASS_COLOURS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22)), labels = scales::comma) +
    facet_grid(patient ~ mode) +
    labs(
      title = if (!is.null(patient)) {
        paste0("Feature classes — ", patient)
      } else {
        "Feature classes — all patients"
      },
      x = NULL, y = "Number of features"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
      plot.title = element_text(face = "bold", size = 12)
    )
}


# Build-all helper: turn long_data (list of patients) into the 3 summaries --
build_summaries <- function(long_data) {
  ds <- list()
  zn <- list()
  ft <- list()
  for (pid in names(long_data)) {
    ln <- long_data[[pid]]$NEG
    lp <- long_data[[pid]]$POS
    ds[[pid]] <- bind_rows(
      summarise_dataset(ln, pid, "NEG"),
      summarise_dataset(lp, pid, "POS")
    )
    zn[[pid]] <- bind_rows(
      summarise_zones(ln, pid, "NEG"),
      summarise_zones(lp, pid, "POS")
    )
    ft[[pid]] <- bind_rows(
      summarise_features(ln, pid, "NEG"),
      summarise_features(lp, pid, "POS")
    )
  }
  feat <- bind_rows(ft)
  list(
    ds_summary = bind_rows(ds),
    zone_summary = bind_rows(zn),
    feat_summary = feat,
    class_counts = feat |>
      count(patient, mode, class, .drop = FALSE) |>
      group_by(patient, mode) |> mutate(pct = 100 * n / sum(n)) |> ungroup()
  )
}

