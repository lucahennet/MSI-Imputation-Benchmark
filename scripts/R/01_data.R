# =============================================================================
# 01_data.R
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

#' Extract ROI coordinates from a single TXT file
#'
#' @param file_path Path to a .txt ROI file
#' @return A tibble with columns ROI, X, Y — or NULL if file is skipped
read_roi_txt <- function(file_path) {
  lines <- readLines(file_path)

  scans <- str_extract(lines[grep("\\*NO Of SCANS:", lines)], "\\d+") |> as.numeric()
  if (is.na(scans) || scans != 1) {
    return(NULL)
  }

  roi_id <- str_extract(lines[grep("\\*NAME:", lines)], "ROI\\s*\\d+") |>
    str_extract("\\d+") |>
    as.numeric()

  coord_line <- lines[grep("\\*REGION SHAPE:", lines)]
  coords <- str_match_all(coord_line, "X=(\\d+), Y=(\\d+)")[[1]]

  if (nrow(coords) == 0) {
    return(NULL)
  }

  tibble(
    ROI = roi_id,
    X   = as.numeric(coords[, 2]),
    Y   = as.numeric(coords[, 3])
  )
}

#' Read one CSV folder and join its data with pixel coordinates
#'
#' @param folder_path Path to a folder containing exactly one .csv file
#' @param coords_df   Tibble of ROI coordinates from read_roi_txt()
#' @return A tibble with sample, X, Y, ROI and all feature columns
process_csv_folder <- function(folder_path, coords_df) {
  csv_file <- list.files(folder_path, "\\.csv$", full.names = TRUE)[1]
  if (is.na(csv_file)) {
    return(NULL)
  }

  read.csv(csv_file) |>
    mutate(ROI = as.numeric(str_extract(ROI, "\\d+$"))) |>
    inner_join(coords_df, by = "ROI") |>
    mutate(sample = basename(folder_path))
}

#' Assemble the full MSI dataset from a root directory
#'
#' Scans for a TXT subfolder (coordinates) and CSV subfolders (features),
#' then joins them and reorders columns.
#'
#' @param root_path Path to the root data directory
#' @return A wide tibble: sample, X, Y, ROI, TIC, then feature columns
assemble_msi_data <- function(root_path) {
  folders <- list.dirs(root_path, recursive = FALSE)

  txt_folder <- folders[map_lgl(folders, ~ any(str_detect(list.files(.x), "\\.txt$")))]
  txt_files <- list.files(txt_folder, "\\.txt$", full.names = TRUE)
  coords_ref <- map_dfr(txt_files, read_roi_txt)

  csv_folders <- folders[map_lgl(folders, ~ any(str_detect(list.files(.x), "\\.csv$")))]
  raw_joined <- map_dfr(csv_folders, process_csv_folder, coords_df = coords_ref)

  raw_joined |>
    select(sample, X, Y, ROI, RAW.TIC.OR.ROI.sum.peak, everything())
}

#' Pivot selected features to long format
#'
#' @param data_wide     Wide tibble from assemble_msi_data()
#' @param feature_indices Integer vector of column indices to pivot
#' @return Long tibble with columns X, Y, TIC, Feature, Intensity
pivot_msi_long <- function(data_wide, feature_indices) {
  target_cols <- colnames(data_wide)[feature_indices]

  data_wide |>
    select(X, Y, RAW.TIC.OR.ROI.sum.peak, all_of(target_cols)) |>
    pivot_longer(
      cols      = all_of(target_cols),
      names_to  = "Feature",
      values_to = "Intensity"
    ) |>
    mutate(Feature = factor(Feature, levels = target_cols))
}

#' Apply TIC normalisation to long-format MSI data
#'
#' @param data_long Long tibble from pivot_msi_long()
#' @return Same tibble with NormIntensity column; TIC column dropped
normalize_msi_data <- function(data_long) {
  data_long |>
    mutate(NormIntensity = Intensity / RAW.TIC.OR.ROI.sum.peak) |>
    select(-RAW.TIC.OR.ROI.sum.peak)
}

#' Summarise the dataset with key metrics
#' 
#' @param raw_wide     Wide tibble from assemble_msi_data()
#' @param feature_range Integer vector of column indices corresponding to features
#' @return A tibble summarising dataset dimensions, NA/zero counts and percentages
dataset_summary <- function(raw_wide, feature_range) {
  n_pixels <- nrow(raw_wide)
  n_features <- length(feature_range) # eventually remove the first one that do not count

  raw_subset <- raw_wide[, feature_range]

  total_cells <- n_pixels * n_features
  na_count <- sum(is.na(raw_subset))
  na_perc <- (na_count / total_cells) * 100

  cols_with_na <- colSums(is.na(raw_subset)) > 0
  num_cols_with_na <- sum(cols_with_na)
  perc_cols_with_na <- (num_cols_with_na / n_features) * 100

  zero_count <- sum(raw_subset == 0, na.rm = TRUE)
  zero_perc <- (zero_count / total_cells) * 100

  cols_with_zero <- colSums(raw_subset == 0, na.rm = TRUE) > 0
  num_cols_with_zero <- sum(cols_with_zero)
  perc_cols_with_zero <- (num_cols_with_zero / n_features) * 100

  rows_with_zero <- sum(apply(raw_subset, 1, function(row) any(row == 0, na.rm = TRUE)))
  perc_rows_with_zero <- (rows_with_zero / n_pixels) * 100

  tibble(
    Metric = c(
      "Total Pixels", "Total Features",
      "Total NAs", "NA Percentage",
      "Features with NAs", "Features (cols) with NAs (%)",
      "Total Zeros", "Zero Percentage", "Features with Zeros",
      "Features (cols) with Zeros (%)", "Rows with Zeros", "Rows with Zeros (%)"
    ),
    Value = c(
      as.character(n_pixels), as.character(n_features),
      as.character(na_count), sprintf("%.2f%%", na_perc),
      as.character(num_cols_with_na), sprintf("%.2f%%", perc_cols_with_na),
      as.character(zero_count), sprintf("%.2f%%", zero_perc), as.character(num_cols_with_zero),
      sprintf("%.2f%%", perc_cols_with_zero), as.character(rows_with_zero),
      sprintf("%.2f%%", perc_rows_with_zero)
    )
  )
}