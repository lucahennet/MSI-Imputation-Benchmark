# =============================================================================
# main_mix-missingness.R
# Purpose:  
#
# Structure:
#   PART 1 — 
#   PART 2 — 
#
#
# Usage:
#
# =============================================================================

library(tidyverse)
library(patchwork)
library(readxl)
library(scales)
library(missForest) # RF
library(MsCoreUtils) # QRILC
library(SpatialPack) # SSIM
library(spdep) # Moran’s I

source("scripts/R/01_data.R")
source("scripts/R/03_imputation.R")
source("scripts/R/04_preprocessing.R")
source("scripts/R/05_metrics.R")
source("scripts/R/10.1_mix_data_analysis.R")
source("scripts/R/10.2_mix_data_benchmark.R")
source("scripts/R/10.3_mix_data_plots.R")


# Config ------------------------------------------------------------------

OUT_DIR     <- "results/missing_structure"
MARKER_XLSX <- "data/20260722_Fibrosis_markers.xlsx"
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE, recursive = TRUE)

THRESH_LOW      <- 0.05
THRESH_MODERATE <- 0.20
THRESH_HIGH     <- 0.50
BIN_WIDTH       <- 5

CLASS_LEVELS <- c("No missing", "Low [<5%]", "Moderate [5-20%]",
                  "High [20-50%]", ">50%")
CLASS_COLOURS <- c(
  "No missing"       = "#1D9E75",
  "Low [<5%]"        = "#378ADD",
  "Moderate [5-20%]" = "#EF9F27",
  "High [20-50%]"    = "#E07B39",
  ">50%"             = "#E24B4A"
)


# Patient list ------------------------------------------------------------

PATIENTS <- list(
  P1 = list(
    id = "P1",
    neg_path = "data/NEG_8976",
    pos_path = "data/POS_8976",
    zone_meta_NEG = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      154, 178, "N1", "norm",
      78, 102, "N2", "norm",
      129, 153, "N3", "norm",
      3, 27, "F1", "fib",
      53, 77, "F2", "fib",
      103, 127, "F3", "fib"
    ),
    zone_meta_POS = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      26, 50, "N1", "norm",
      79, 103, "N2", "norm",
      130, 154, "N3", "norm",
      1, 25, "F1", "fib",
      52, 76, "F2", "fib",
      155, 156, "F2", "fib", # small extension
      104, 128, "F3", "fib"
    )
  ),
  P2 = list(
    id = "P2",
    neg_path = "data/NEG_32244",
    pos_path = "data/POS_32244",
    zone_meta_NEG = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      27, 51, "N1", "norm",
      77, 100, "N2", "norm",
      151, 152, "N2", "norm", # small extension
      126, 150, "N3", "norm",
      1, 26, "F1", "fib",
      52, 76, "F2", "fib",
      101, 125, "F3", "fib"
    ),
    zone_meta_POS = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      76, 100, "N1", "norm",
      101, 125, "N2", "norm",
      126, 150, "N3", "norm",
      1, 25, "F1", "fib",
      26, 50, "F2", "fib",
      151, 151, "F2", "fib", # single-ROI extension
      51, 75, "F3", "fib"
    )
  ),
  P3 = list(
    id = "P3",
    neg_path = "data/NEG_32216",
    pos_path = "data/POS_32216",
    zone_meta_NEG = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      76, 100, "N1", "norm",
      101, 125, "N2", "norm",
      126, 150, "N3", "norm",
      1, 25, "F1", "fib",
      26, 50, "F2", "fib",
      51, 75, "F3", "fib"
    ),
    zone_meta_POS = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      52, 76, "N1", "norm",
      77, 101, "N2", "norm",
      102, 126, "N3", "norm",
      1, 26, "F1", "fib",
      27, 51, "F2", "fib",
      127, 151, "F3", "fib"
    )
  ),
  P4 = list(
    id = "P4",
    neg_path = "data/NEG_17342",
    pos_path = "data/POS_17342",
    zone_meta_NEG = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      76, 100, "N1", "norm",
      101, 125, "N2", "norm",
      126, 150, "N3", "norm",
      1, 25, "F1", "fib",
      26, 50, "F2", "fib",
      51, 75, "F3", "fib"
    ),
    zone_meta_POS = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      76, 100, "N1", "norm",
      101, 125, "N2", "norm",
      126, 150, "N3", "norm",
      1, 25, "F1", "fib",
      26, 50, "F2", "fib",
      51, 75, "F3", "fib"
    )
  ),
  P5 = list(
    id = "P5",
    neg_path = "data/NEG_17139",
    pos_path = "data/POS_17139",
    zone_meta_NEG = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      76, 100, "N1", "norm",
      101, 125, "N2", "norm",
      126, 150, "N3", "norm",
      1, 25, "F1", "fib",
      26, 50, "F2", "fib",
      51, 75, "F3", "fib"
    ),
    zone_meta_POS = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      76, 100, "N1", "norm",
      101, 125, "N2", "norm",
      126, 150, "N3", "norm",
      1, 25, "F1", "fib",
      26, 50, "F2", "fib",
      51, 75, "F3", "fib"
    )
  ),
  P6 = list(
    id = "P6",
    neg_path = "data/NEG_20458",
    pos_path = "data/POS_20458",
    zone_meta_NEG = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      76, 100, "N1", "norm",
      101, 125, "N2", "norm",
      126, 150, "N3", "norm",
      1, 25, "F1", "fib",
      26, 50, "F2", "fib",
      51, 75, "F3", "fib"
    ),
    zone_meta_POS = tribble(
      ~ROI_min, ~ROI_max, ~zone_id, ~tissue_type,
      76, 100, "N1", "norm",
      101, 125, "N2", "norm",
      126, 150, "N3", "norm",
      1, 25, "F1", "fib",
      26, 50, "F2", "fib",
      51, 75, "F3", "fib"
    )
  )
)


# =============================================================================
# DATASET MISSINGNESS ANALYSIS
# =============================================================================

# ================= PART 1: load data, full-feature objects ===================

# Pre-load all raw matrices
raw_mats <- list()
for (pt in PATIENTS) {
  raw_mats[[pt$id]] <- list(NEG = load_mat(pt$neg_path), POS = load_mat(pt$pos_path))
}

# Common features across all patients (per mode)
common_neg <- Reduce(intersect, map(raw_mats, \(pt) colnames(pt$NEG$mat)))
common_pos <- Reduce(intersect, map(raw_mats, \(pt) colnames(pt$POS$mat)))
message(sprintf("NEG aligned: %d common features", length(common_neg)))
message(sprintf("POS aligned: %d common features", length(common_pos)))

# Build long_data on the FULL common feature set
long_data <- list()
for (pt in PATIENTS) {
  message("Processing patient ", pt$id)
  neg <- raw_mats[[pt$id]]$NEG
  pos <- raw_mats[[pt$id]]$POS
  neg$mat <- neg$mat[, common_neg, drop = FALSE]
  pos$mat <- pos$mat[, common_pos, drop = FALSE]
  coords_neg <- annotate_zones(neg$coords, pt$zone_meta_NEG)
  coords_pos <- annotate_zones(pos$coords, pt$zone_meta_POS)
  long_data[[pt$id]] <- list(
    NEG = build_long(neg$mat, coords_neg, pt$id, "NEG"),
    POS = build_long(pos$mat, coords_pos, pt$id, "POS")
  )
}
message("Done loading.")

# Full-feature summaries
full <- build_summaries(long_data)
ds_summary   <- full$ds_summary
zone_summary <- full$zone_summary
feat_summary <- full$feat_summary
class_counts <- full$class_counts


# ================= PART 2: select fibrosis-marker features ===================

markers_raw <- read_excel(MARKER_XLSX)

markers <- markers_raw |>
  mutate(
    mode       = str_extract(MARKER, "(NEG|POS)"),
    mz         = as.numeric(str_extract(MARKER, "(?<=(NEG|POS)_)\\d+\\.?\\d*")),
    annotation = str_trim(str_remove(MARKER, "\\s*\\((NEG|POS)_[\\d.]+\\)$"))
  ) |>
  select(annotation,
         marker_raw = MARKER, mode, mz,
         p_value = P_VALUE, FDR, FWER,
         mean_fib = Fibrosis, mean_nonfib = `Non-Fibrosis`, effect_D = D
  )

message(sprintf(
  "\nMarkers in Excel: %d (%d NEG, %d POS)",
  nrow(markers), sum(markers$mode == "NEG"), sum(markers$mode == "POS")
))

# Exact m/z match against our common features (mode-aware)
feat_lookup <- bind_rows(
  tibble(mode = "NEG", feature = common_neg, feat_mz = feature_mz(common_neg)),
  tibble(mode = "POS", feature = common_pos, feat_mz = feature_mz(common_pos))
) |> filter(!is.na(feat_mz))

marker_matches <- markers |>
  left_join(feat_lookup, by = c("mode" = "mode", "mz" = "feat_mz")) |>
  mutate(matched = !is.na(feature))

n_matched <- sum(marker_matches$matched)
n_unmatched <- nrow(marker_matches) - n_matched

message("\n=== Marker matching summary ===")
message(sprintf("  Markers total:     %d", nrow(marker_matches)))
message("\n  By mode:")
print(marker_matches |> group_by(mode) |>
        summarise(
          n_markers = n(), n_matched = sum(matched),
          n_unmatched = sum(!matched), .groups = "drop"
        ))

if (n_unmatched > 0) {
  message("\n=== Unmatched markers (m/z absent from our common features) ===")
  print(marker_matches |> filter(!matched) |> select(annotation, mode, mz))
}

selected_features <- marker_matches |>
  filter(matched) |>
  select(annotation, mode, mz, feature, p_value, effect_D, mean_fib, mean_nonfib)
selected_neg <- selected_features |>
  filter(mode == "NEG") |>
  pull(feature)
selected_pos <- selected_features |>
  filter(mode == "POS") |>
  pull(feature)

message(sprintf(
  "\nSelected marker features present in data: %d NEG, %d POS",
  length(selected_neg), length(selected_pos)
))
# write_csv(selected_features, file.path(OUT_DIR, "tables", "selected_fibrosis_markers.csv"))


# ================= PART 3: marker-only objects ===============================

long_data_markers <- map(long_data, function(pt) {
  list(
    NEG = pt$NEG |> filter(feature %in% selected_neg),
    POS = pt$POS |> filter(feature %in% selected_pos)
  )
})

# Same summaries, computed on the marker subset
markers_sum <- build_summaries(long_data_markers)
ds_summary_markers <- markers_sum$ds_summary
zone_summary_markers <- markers_sum$zone_summary
feat_summary_markers <- markers_sum$feat_summary
class_counts_markers <- markers_sum$class_counts


# ================= PART 4: summaries & plots =================================

# --- GLOBAL summary (full feature set) ---
message("\n=== Global dataset summary (ALL features) ===")
ds_summary |>
  select(
    patient, mode, n_pixels, n_features, n_missing, pct_missing,
    n_features_w_missing, pct_features_missing,
    n_pixels_w_missing, pct_pixels_missing
  ) |>
  print(n = Inf)

# --- GLOBAL summary (marker subset) ---
message("\n=== Dataset summary (FIBROSIS MARKERS only) ===")
ds_summary_markers |>
  select(
    patient, mode, n_features, n_missing, pct_missing,
    n_features_w_missing, pct_features_missing
  ) |>
  print(n = Inf)


# --- Example specific summaries (swap feat_summary <-> feat_summary_markers) ---

print_freq_table(feat_summary_markers |> filter(patient=="P2", mode=="NEG"), "P2", "NEG")
print_class_summary(feat_summary_markers |> filter(patient=="P2", mode=="NEG"), "P2", "NEG")

plot_histogram(feat_summary_markers |> filter(patient=="P1", mode=="NEG"), "NEG", "P1")
plot_zone_summary(zone_summary_markers, patient = "P1")
plot_zone_summary(zone_summary_markers)              # all patients, markers only
plot_feature_classes(class_counts_markers)           # all patients, markers 


# =============================================================================
# FEATURE MISSINGNESS ANALYSIS
# =============================================================================

# Purpose:  Investigate whether LOW-INTENSITY features are more prone to
#           missingness (MNAR / below-LOD) than high-intensity features
#
# Two approaches:
#   A (primary, threshold): per-FEATURE relationship between mean observed
#     intensity and missingness rate. If low mean intensity predicts high
#     missingness -> MNAR-consistent. Fit a threshold on mean intensity.
#   B (confirmation): WITHIN-feature, are observed intensities lower in the
#     zones where the feature also goes missing? Direct below-LOD signature.
#
# Levels:
#   - Per patient x mode (check the relationship is consistent)
#   - Pooled across all 6 patients (final threshold per mode)

# Approach A gives a global trend
# Approach B gives local evidence

 
if (!exists("long_data_markers")) stop("Run missing_structure.R first (need long_data).")


# ================= APPROACH A ================================================
# Are features that have low signal intensity generally more affected by missing values?

# per-feature: mean observed intensity vs. missingness rate

# For each feature (within a patient x mode), summarise:
#   - mean/median OBSERVED intensity (only non-missing pixels)
#   - missingness rate
# Then relate the two.

feature_intensity_missing <- function(long_df, patient, mode) {
  long_df |>
    group_by(feature) |>
    summarise(
      n_pixels = n(),
      n_missing = sum(missing),
      prop_missing = mean(missing),
      # observed-only intensity summaries
      mean_obs_int = mean(intensity[!missing]),
      median_obs_int = median(intensity[!missing]),
      min_obs_int = ifelse(any(!missing), min(intensity[!missing]), NA_real_),
      .groups = "drop"
    ) |>
    filter(n_missing < n_pixels) |> # need at least some observed values
    mutate(patient = patient, mode = mode)
}

fim_all <- map_dfr(names(long_data), function(pid) {
  bind_rows(
    feature_intensity_missing(long_data[[pid]]$NEG, pid, "NEG"),
    feature_intensity_missing(long_data[[pid]]$POS, pid, "POS")
  )
})


# --- Correlation: does low intensity go with high missingness? -----------
# Spearman (rank) correlation is robust and doesn't assume linearity.
# NEGATIVE correlation = lower intensity -> more missing = MNAR-consistent.

corr_by_group <- fim_all |>
  group_by(patient, mode) |>
  summarise(
    n_features = n(),
    # correlate log intensity with missingness rate
    spearman_rho = cor(log1p(mean_obs_int), prop_missing,
      method = "spearman", use = "complete.obs"
    ),
    .groups = "drop"
  )

# Pooled (all patients) per mode
corr_pooled <- fim_all |>
  group_by(mode) |>
  summarise(
    n_features = n(),
    spearman_rho = cor(log1p(mean_obs_int), prop_missing,
      method = "spearman", use = "complete.obs"
    ),
    .groups = "drop"
  )

message("\n=== Approach A: intensity vs missingness (Spearman rho) ===")
message("  (negative rho = low intensity -> more missing = MNAR-consistent)")
message("\n  Per patient x mode:")
print(corr_by_group)
message("\n  Pooled per mode:")
print(corr_pooled)


# --- Threshold derivation -------------------------------------------------
# Group features into intensity bins and show missingness rate per bin.
# The intensity level where missingness climbs sharply is the MNAR/MCAR split.

INTENSITY_BREAKS <- c(0, 0.5, 1, 2, 5, 10, 50, 100, Inf)

threshold_table <- fim_all |>
  mutate(int_bin = cut(mean_obs_int,
    breaks = INTENSITY_BREAKS,
    include.lowest = TRUE, right = FALSE
  )) |>
  group_by(mode, int_bin) |>
  summarise(
    n_features = n(),
    mean_missingness = 100 * mean(prop_missing),
    median_missingness = 100 * median(prop_missing),
    pct_high_missing = 100 * mean(prop_missing > 0.2), # frac of features >20% missing
    .groups = "drop"
  )

message("\n=== Approach A: missingness by intensity bin (pooled) ===")
message("  (if missingness is high at low intensity and drops off, that")
message("   transition point is your MNAR/MCAR threshold)")
print(threshold_table, n = Inf)


# ================= APPROACH B ================================================
# Does the feature disappear preferentially in regions where it is already weak?

# within-feature: observed intensity where feature goes missing

# For each feature, compare the observed intensities in zones that HAVE
# missing pixels vs zones with NO missing pixels for that feature.
# If observed values are lower in missing-affected zones, that is the
# below-LOD (MNAR) signature at the within-feature level.

within_feature_shift <- function(long_df, patient, mode) {
  # Per feature x zone: missing count and observed-intensity median
  fz <- long_df |>
    filter(!is.na(zone_id)) |>
    group_by(feature, zone_id) |>
    summarise(
      n_missing = sum(missing),
      median_obs_int = ifelse(any(!missing), median(intensity[!missing]), NA_real_),
      .groups = "drop"
    ) |>
    mutate(zone_has_missing = n_missing > 0)

  # Per feature: compare median observed intensity in missing vs non-missing zones
  fz |>
    filter(!is.na(median_obs_int)) |>
    group_by(feature) |>
    summarise(
      n_zones = n(),
      n_zones_missing = sum(zone_has_missing),
      n_zones_clean = sum(!zone_has_missing),
      med_int_missing_zones = median(median_obs_int[zone_has_missing]),
      med_int_clean_zones = median(median_obs_int[!zone_has_missing]),
      .groups = "drop"
    ) |>
    filter(n_zones_missing > 0, n_zones_clean > 0) |> # need both to compare
    mutate(
      # negative = observed intensity LOWER in missing-affected zones = MNAR-like
      log2fc = log2(med_int_missing_zones / med_int_clean_zones),
      patient = patient, mode = mode
    )
}

wfs_all <- map_dfr(names(long_data), function(pid) {
  bind_rows(
    within_feature_shift(long_data[[pid]]$NEG, pid, "NEG"),
    within_feature_shift(long_data[[pid]]$POS, pid, "POS")
  )
})

# Summary: what fraction of features show LOWER observed intensity in
# missing-affected zones (MNAR-consistent)?
wfs_summary <- wfs_all |>
  group_by(mode) |>
  summarise(
    n_features = n(),
    frac_mnar_like = mean(log2fc < 0), # lower intensity where missing
    median_log2fc = median(log2fc),
    # Wilcoxon signed-rank: is log2fc systematically < 0?
    wilcox_p = wilcox.test(log2fc, mu = 0)$p.value,
    .groups = "drop"
  )

message("\n=== Approach B: within-feature intensity shift (confirmation) ===")
message("  frac_mnar_like = fraction of features with LOWER observed intensity")
message("  in zones where they also go missing (MNAR-consistent).")
print(wfs_summary)


# ================= PLOTS =====================================================

# Plot A1: scatter of mean intensity vs missingness (pooled, per mode)
pA1 <- fim_all |>
  ggplot(aes(x = mean_obs_int, y = 100 * prop_missing)) +
  geom_point(alpha = 0.25, size = 1, colour = "#378ADD") +
  geom_smooth(method = "loess", se = TRUE, colour = "#E24B4A", linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  annotate("text",
    x = 1, y = 95, label = "intensity = 1", hjust = -0.1,
    size = 2.8, colour = "grey40"
  ) +
  scale_x_log10() +
  facet_wrap(~mode) +
  labs(
    title = "Approach A: feature missingness vs mean observed intensity",
    subtitle = "Each point = one feature (all patients pooled). x on log scale.",
    x = "Mean observed intensity (log scale)",
    y = "Missingness (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 9)
  )
pA1

# Plot A2: missingness by intensity bin (the threshold view)
pA2 <- threshold_table |>
  ggplot(aes(x = int_bin, y = mean_missingness, fill = mode)) +
  geom_col(position = position_dodge(), width = 0.7) +
  geom_text(aes(label = n_features),
    position = position_dodge(width = 0.7),
    vjust = -0.3, size = 2.5
  ) +
  scale_fill_manual(values = c(NEG = "#378ADD", POS = "#E07B39"), name = NULL) +
  labs(
    title = "Approach A: mean missingness by intensity bin",
    subtitle = "Bar labels = number of features in each bin. High missingness at low intensity = MNAR-consistent.",
    x = "Mean observed intensity bin", y = "Mean missingness (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 8)
  )
pA2

# Plot B: distribution of within-feature log2fc
pB <- wfs_all |>
  ggplot(aes(x = log2fc, fill = mode)) +
  geom_histogram(
    bins = 40, colour = "white", linewidth = 0.15,
    position = "identity", alpha = 0.6
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
  annotate("text",
    x = 0, y = Inf, label = "  MNAR-like <-- | --> not MNAR",
    hjust = 0.5, vjust = 1.5, size = 2.8, colour = "grey40"
  ) +
  scale_fill_manual(values = c(NEG = "#378ADD", POS = "#E07B39"), name = NULL) +
  facet_wrap(~mode, ncol = 1) +
  labs(
    title = "Approach B: within-feature intensity shift",
    subtitle = "log2(median intensity in missing-affected zones / clean zones). Negative = MNAR-consistent.",
    x = "log2 fold-change", y = "Number of features"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none", strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(colour = "grey40", size = 8)
  )
pB

ggsave(file.path(OUT_DIR, "intensity_A1_scatter.pdf"),   pA1, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "intensity_A2_threshold.pdf"), pA2, width = 10, height = 5)
ggsave(file.path(OUT_DIR, "intensity_B_shift.pdf"),      pB,  width = 8,  height = 6)

message("\nPlots saved to: ", OUT_DIR)


# =============================================================================
# MIX DATA BENCHMARK
# =============================================================================

CFG <- list(
  # props        = c(0.05, 0.10, 0.20, 0.30, 0.40),  # 0.01 dropped: 1 px at n=25
  props        = c(0.05),
  # reps         = 3,
  reps         = 1,
  mnar_frac    = 0.5,      # fraction of features assigned MNAR (at random)
  min_features = 20,       # skip a zone with fewer complete features
  min_observed = 5,        # keep at least this many observed pixels per feature
  seed_base    = 42,
  checkpoint_dir = "results/benchmark/checkpoints",
  out_dir        = "results/benchmark"
)

SPECS <- list(
  RF    = list(fun = impute_missForest, transform = "log1p", scaling = "none"),
  QRILC = list(fun = impute_qrilc,      transform = "log1p", scaling = "none")
)

validate_setup()
zones <- zones_from_long(long_data_markers)
dry_run(zones)
res  <- run_benchmark(zones, resume = FALSE)
summ <- summarise_benchmark(res)
save_benchmark(res, zones)        # add zones to enable heatmaps on reload


# =============================================================================
# MIX DATA PLOTS
# =============================================================================

res <- readRDS("results/benchmark/zone_benchmark_3reps_2026-07-27.rds")$results
res <- readRDS("results/benchmark/zone_benchmark_1reps_2026-07-27.rds")$results
all_tests(res)
plot_metric(res, "RMSE", "MNAR")
plot_gain(res)
ins <- inspect_unit(zones, "P1_NEG_N1", prop = 0.05); visualise_heatmaps(ins, 1)

plot_metric(res, "RMSE",  subset = "overall")              # line
plot_metric(res, "RMSE",  subset = "MNAR",    type = "bar") # bar, faceted by prop
plot_metric(res, "NRMSE", subset = "MCAR",    type = "bar")
plot_metric(res, "NRMSE", subset = "overall")
plot_metric(res, "SAM",   subset = "overall")
plot_gain(res)          # alias -> plot_gain
plot_spatial_fidelity(res)     # uses SSIM + MoranDiff
plot_spectral_preservation(res) # uses VarRatio + CorStruct
