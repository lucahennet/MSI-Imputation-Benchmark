# =============================================================================
# 03_benchmark.R
# Imputation benchmark on small zones with MIXED missing mechanisms.
#
# Methods compared:
#   RF     - missForest on every feature
#   QRILC  - QRILC on every feature
#   Oracle - route each feature to the method matching its TRUE simulated
#            mechanism (RF for MCAR, QRILC for MNAR). Upper bound on what
#            perfect mechanism knowledge is worth.
#
# NOTE — no intensity-threshold method here, deliberately. The simulation
# assigns MNAR at random with respect to intensity, so an intensity-based
# routing rule cannot be evaluated fairly against it (its premise and the
# generator disagree). The intensity finding lives in 02_intensity.R as a
# DESCRIPTIVE result about the real data, not as an imputation rule tested
# on synthetic missingness.
#
# INPUT CONTRACT (reusable across experiments):
#   run_benchmark(zones): named list, each element
#       list(mat = <pixels x features, no NA/zeros>, coords = <X,Y per pixel>)
#
# Prereq: 00_setup.R (functions), 03_imputation.R, 04_preprocessing.R,
#         05_metrics.R.
# =============================================================================


# ============================ SETUP CHECK ====================================

validate_setup <- function() {
  ok <- TRUE
  for (p in c("missForest", "imputeLCMD", "MsCoreUtils", "tidyverse")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("MISSING package: ", p)
      ok <- FALSE
    }
  }
  for (f in c(
    "impute_missForest", "impute_qrilc", "apply_preprocessing",
    "reverse_preprocessing", "compute_metrics"
  )) {
    if (!exists(f)) {
      message("MISSING function: ", f)
      ok <- FALSE
    }
  }
  if (!exists("metric_functions")) {
    message("MISSING: metric_functions")
    ok <- FALSE
  }
  if (ok) {
    m <- matrix(runif(25 * 30, 1, 10), 25, 30, dimnames = list(NULL, paste0("f", 1:30)))
    m[sample(750, 60)] <- NA
    ok <- tryCatch(
      {
        impute_qrilc(log1p(m))
        TRUE
      },
      error = function(e) {
        message(
          "QRILC test failed: ",
          conditionMessage(e)
        )
        FALSE
      }
    )
  }
  message(if (ok) "Setup OK." else "Setup INCOMPLETE - fix the above first.")
  invisible(ok)
}


# ============================== HELPERS ======================================

impute_with_prep <- function(mat_na, spec) {
  prep <- apply_preprocessing(mat_na, spec)
  reverse_preprocessing(spec$fun(prep$data, NULL), prep)
}

n_to_mask <- function(n, prop, min_observed) {
  k <- max(1L, as.integer(ceiling(prop * n)))
  min(k, n - min_observed)
}

#' Assign each feature MCAR or MNAR at random, then apply missingness.
simulate_missing <- function(mat, prop, mnar_frac, min_observed) {
  n <- nrow(mat)
  p <- ncol(mat)
  k <- n_to_mask(n, prop, min_observed)
  if (k < 1) {
    return(NULL)
  }

  mnar_idx <- sample.int(p, round(mnar_frac * p))
  is_mnar <- seq_len(p) %in% mnar_idx

  mat_na <- mat
  for (j in seq_len(p)) {
    idx <- if (is_mnar[j]) order(mat[, j])[seq_len(k)] else sample.int(n, k)
    mat_na[idx, j] <- NA
  }
  list(
    mat_na = mat_na,
    mechanism = if_else(is_mnar, "MNAR", "MCAR"),
    n_masked = k
  )
}

#' Score one prediction: overall + per true mechanism.
#' Spatial metrics (need coords + full matrix) attach to the `overall` row only.
score <- function(truth, pred, mask, mechanism, method, coords = NULL) {
  strata <- c(
    list(overall = mask),
    map(set_names(c("MCAR", "MNAR")), function(mech) {
      m <- mask
      m[, mechanism != mech] <- FALSE
      m
    })
  )
  spatial <- if (!is.null(coords)) extra_metrics(truth, pred, coords, mask) else list()

  imap_dfr(strata, function(msk, nm) {
    if (!any(msk)) {
      return(NULL)
    }
    row <- as_tibble_row(compute_metrics(truth, pred, msk, metric_functions)) |>
      mutate(subset = nm)
    if (nm == "overall" && length(spatial) > 0) {
      row <- bind_cols(row, as_tibble_row(unlist(spatial)))
    }
    row
  }) |>
    mutate(method = method)
}


# ============================= CORE UNIT =====================================

run_unit <- function(zone, zone_id, prop, rep, cfg = CFG, return_matrices = FALSE) {
  mat <- zone$mat
  set.seed(cfg$seed_base + rep * 1000L + as.integer(prop * 1000))

  sim <- simulate_missing(mat, prop, cfg$mnar_frac, cfg$min_observed)
  if (is.null(sim)) {
    return(NULL)
  }

  mask <- is.na(sim$mat_na)

  rf <- impute_with_prep(sim$mat_na, SPECS$RF)
  qr <- impute_with_prep(sim$mat_na, SPECS$QRILC)

  # Oracle: route by TRUE mechanism (upper bound)
  oracle <- rf
  oracle[, sim$mechanism == "MNAR"] <- qr[, sim$mechanism == "MNAR"]

  scores <- bind_rows(
    score(mat, rf, mask, sim$mechanism, "RF", zone$coords),
    score(mat, qr, mask, sim$mechanism, "QRILC", zone$coords),
    score(mat, oracle, mask, sim$mechanism, "Oracle", zone$coords)
  ) |>
    mutate(
      zone = zone_id, prop = prop, rep = rep,
      n_pixels = nrow(mat), n_features = ncol(mat),
      n_masked_per_feature = sim$n_masked,
      pct_masked_actual = 100 * sim$n_masked / nrow(mat),
      n_mnar = sum(sim$mechanism == "MNAR")
    )

  if (!return_matrices) {
    return(scores)
  }

  list(
    scores = scores, coords = zone$coords, truth = mat,
    masked = sim$mat_na,
    imputed = list(RF = rf, QRILC = qr, Oracle = oracle),
    mechanism = sim$mechanism
  )
}


# ============================== DRIVER =======================================

dry_run <- function(zones, cfg = CFG) {
  z <- zones[[1]]
  p_max <- max(cfg$props)
  message(sprintf(
    "Test unit: %s (%d px x %d features)",
    names(zones)[1], nrow(z$mat), ncol(z$mat)
  ))
  t <- system.time(run_unit(z, names(zones)[1], p_max, 1, cfg))[["elapsed"]]
  n_units <- length(zones) * length(cfg$props) * cfg$reps
  message(sprintf(
    "One unit at prop=%.2f: %.1f s  ->  %d units ~ %.1f min",
    p_max, t, n_units, t * n_units / 60
  ))
  invisible(t)
}

run_benchmark <- function(zones, cfg = CFG, resume = TRUE) {
  dir.create(cfg$checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  grid <- expand_grid(zone = names(zones), prop = cfg$props, rep = seq_len(cfg$reps))
  message(sprintf("%d units to run.", nrow(grid)))
  results <- vector("list", nrow(grid))
  t0 <- Sys.time()

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    key <- sprintf("%s_p%03d_r%d", g$zone, round(g$prop * 100), g$rep)
    ckpt <- file.path(cfg$checkpoint_dir, paste0(key, ".rds"))
    if (resume && file.exists(ckpt)) {
      results[[i]] <- readRDS(ckpt)
      next
    }

    z <- zones[[g$zone]]
    if (ncol(z$mat) < cfg$min_features) {
      message("skip ", g$zone, ": only ", ncol(z$mat), " features")
      next
    }
    res <- tryCatch(run_unit(z, g$zone, g$prop, g$rep, cfg),
      error = function(e) {
        message("FAILED ", key, ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(res)) {
      saveRDS(res, ckpt)
      results[[i]] <- res
    }

    if (i %% 10 == 0 || i == nrow(grid)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message(sprintf(
        "  %d/%d  (%.1f min, ~%.1f left)",
        i, nrow(grid), el, el / i * (nrow(grid) - i)
      ))
    }
  }
  out <- bind_rows(results)
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(out, file.path(cfg$out_dir, "benchmark_results.csv"))
  out
}


# ============================= SUMMARISE =====================================

summarise_benchmark <- function(res) {
  by_subset <- res |>
    group_by(prop, subset, method) |>
    summarise(
      RMSE_mean = mean(RMSE), RMSE_sd = sd(RMSE),
      RMSE_se = sd(RMSE) / sqrt(n()),
      MAE_mean = mean(MAE), NRMSE_mean = mean(NRMSE),
      n = n(), .groups = "drop"
    ) |>
    arrange(subset, prop, RMSE_mean)

  headline <- res |>
    filter(subset == "overall") |>
    group_by(prop, method) |>
    summarise(RMSE = mean(RMSE), .groups = "drop") |>
    pivot_wider(names_from = method, values_from = RMSE) |>
    mutate(
      best_global = pmin(RF, QRILC),
      oracle_gain_pct = 100 * (best_global - Oracle) / best_global
    )

  message("\n=== Method RMSE by subset ===")
  print(by_subset, n = Inf)
  message("\n=== Oracle gain over best global ===")
  print(headline)

  list(by_subset = by_subset, headline = headline)
}


# ====================== ZONE EXTRACTION (adapter) ============================

zones_from_long <- function(long_data) {
  out <- list()
  for (pid in names(long_data)) {
    for (md in names(long_data[[pid]])) {
      ld <- long_data[[pid]][[md]] |> filter(!is.na(zone_id))
      for (z in unique(ld$zone_id)) {
        w <- ld |>
          filter(zone_id == z) |>
          select(X, Y, feature, intensity) |>
          pivot_wider(names_from = feature, values_from = intensity) |>
          arrange(X, Y)
        m <- as.matrix(w |> select(-X, -Y))
        keep <- colSums(m == 0) == 0 & colSums(is.na(m)) == 0
        out[[paste(pid, md, z, sep = "_")]] <- list(
          mat = m[, keep, drop = FALSE], coords = w |> select(X, Y)
        )
      }
    }
  }
  out
}


# ============================= SAVE / LOAD ===================================

save_benchmark <- function(res, zones = NULL, cfg = CFG, file = NULL) {
  bundle <- list(
    config = cfg, results = res,
    summary = summarise_benchmark(res), zones = zones
  )
  if (is.null(file)) {
    file <- file.path(
      cfg$out_dir,
      sprintf("zone_benchmark_%dreps_%s.rds", cfg$reps, Sys.Date())
    )
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(bundle, file)
  message("Saved: ", file)
  invisible(file)
}