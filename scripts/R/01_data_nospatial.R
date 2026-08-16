# =============================================================================
# 01_data_nospatial.R
# Purpose:  Load and process non-spatial MSI datasets from Excel files
# Inputs:   - A root directory containing one Excel file
# Outputs:  - 
# Depends:  tidyverse, stringr, readxl, zoo
# =============================================================================


# Functions ---------------------------------------------------------------

#' Load a non-spatial MSI dataset from an Excel file, and reshape it into a tidy format.
#' 
#' Assumes the first 3 rows contain metadata (acquisition mode, lipid names, m/z values) 
#' and the rest are data rows. The function extracts acquisition mode and m/z to 
#' create feature column names, and parses sample metadata from the SampleName column.
#' @param path Path to the Excel file containing the MSI data
#' @return A tidy data frame with columns for Group, Regions, Individual, Acquisition,
#' ROI_num, SampleName, and one column per feature (e.g. "NEG_253.216", "POS_279.332", etc.)
load_nospatial_msi <- function(path) {
  raw <- suppressMessages(read_excel(path, col_names = FALSE))

  # Row 1: acquisition mode (NEG/POS), propagate forward across merged cells
  acq_row <- as.character(raw[1, ])
  acq_filled <- zoo::na.locf(ifelse(acq_row == "NEG" | acq_row == "POS", acq_row, NA), na.rm = FALSE)
  acq_filled[1:3] <- NA # first 3 cols are metadata

  # Row 2: lipid names; Row 3: m/z values (also used as col IDs)
  lipid_names <- as.character(raw[2, ])
  mz_row <- as.character(raw[3, ]) # "Group", "Regions", "m/z", then numeric m/z

  # Build column names: "NEG_253.216", "POS_279.332", etc.
  mz_numeric <- suppressWarnings(as.numeric(mz_row))
  mz_formatted <- ifelse(is.na(mz_numeric), mz_row, sprintf("%.4f", mz_numeric))
  feature_cols <- paste0(acq_filled, "_", mz_formatted)
  feature_cols[1:3] <- c("Group", "Regions", "SampleName")

  # Data rows
  data_raw <- raw[4:nrow(raw), ]
  colnames(data_raw) <- feature_cols

  # Parse sample metadata from SampleName (e.g. "NEG_M17_ROI 10")
  data_raw |>
    mutate(across(starts_with("NEG_") | starts_with("POS_"), as.numeric)) |>
    mutate(
      Individual  = str_extract(SampleName, "M_?\\d+"),
      ROI_num     = as.integer(str_extract(SampleName, "(?<=ROI )\\d+")),
      Acquisition = str_extract(SampleName, "NEG|POS")
    ) |>
    select(Group, Regions, Individual, Acquisition, ROI_num, SampleName, everything())
}

dataset_summary_nospatial <- function(mat, meta) {
  n_rois <- nrow(mat)
  n_features <- ncol(mat)

  total_cells <- n_rois * n_features
  na_count <- sum(is.na(mat))
  na_perc <- na_count / total_cells * 100
  cols_with_na <- sum(colSums(is.na(mat)) > 0)
  zero_count <- sum(mat == 0, na.rm = TRUE)
  zero_perc <- zero_count / total_cells * 100
  cols_with_zero <- sum(colSums(mat == 0, na.rm = TRUE) > 0)
  rows_with_zero <- sum(apply(mat, 1, function(r) any(r == 0, na.rm = TRUE)))

  tibble(
    Metric = c(
      "Total ROIs", "Total Features",
      "Groups", "Regions", "Individuals",
      "Total NAs", "NA (%)",
      "Features with NAs",
      "Total Zeros", "Zero (%)",
      "Features with Zeros", "ROIs with Zeros"
    ),
    Value = c(
      n_rois, n_features,
      paste(sort(unique(meta$Group)), collapse = ", "),
      paste(sort(unique(meta$Regions)), collapse = ", "),
      paste(sort(unique(meta$Individual)), collapse = ", "),
      na_count, sprintf("%.2f%%", na_perc),
      cols_with_na,
      zero_count, sprintf("%.2f%%", zero_perc),
      cols_with_zero, rows_with_zero
    )
  )
}