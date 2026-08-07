## ---------------------------------------------------------
## create_analysis_data.R
## ---------------------------------------------------------
## Optional data-preparation script for one or more analysis setups.
## Set the parameters in scripts/analysis_setup.R, then source this script to create the
## setup-specific data files used by scripts/create_tables_and_figures.R.

library(tidyverse)
library(foreach)
library(doParallel)
library(here)

source(here("scripts", "analysis_setup.R"))

## Missing counterfactual z-/p-value files are created separately for each
## estimator. Set this to TRUE (for example in scripts/main.R) to overwrite
## existing matching files. These steps can be very time consuming with
## n_iterations <- 1000.
recreate_counterfactuals <- if (exists("recreate_counterfactuals")) recreate_counterfactuals else FALSE

save_counterfactual <- function(path, calculate, progress_bar = NULL,
                                progress_value = NULL, expected_iterations = NULL) {
  cached_result_is_stale <- FALSE
  if (file.exists(path) && !is.null(expected_iterations)) {
    cached_result <- readRDS(path)
    cached_result_is_stale <- !is.list(cached_result) ||
      length(cached_result) == 0 ||
      !isTRUE(nrow(cached_result[[1]]) == expected_iterations)
  }
  if (recreate_counterfactuals || !file.exists(path) || cached_result_is_stale) {
    saveRDS(calculate(), path)
  }
  if (!is.null(progress_value)) update_progress_bar(progress_bar, progress_value)
}

ensure_output_dirs <- function() {
  invisible(lapply(
    list(
      here("results", "main"),
      here("data", "derived_data", "pet_peese_rstandard"),
      here("data", "derived_data", "pet_peese_all_data"),
      here("data", "derived_data", "multilevel_random"),
      here("data", "derived_data", "multilevel_random_all_data"),
      here("data", "derived_data")
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

load_estimator_data <- function(estimator, outlier_variant = "outliers_removed") {
  estimator_dir <- switch(
    paste(estimator, outlier_variant, sep = "_"),
    pet_peese_outliers_removed = here("data", "derived_data", "pet_peese_rstandard"),
    pet_peese_all_data = here("data", "derived_data", "pet_peese_all_data"),
    multilevel_random_outliers_removed = here("data", "derived_data", "multilevel_random"),
    multilevel_random_all_data = here("data", "derived_data", "multilevel_random_all_data"),
    stop("Unknown estimator/outlier variant: ", estimator, "/", outlier_variant)
  )

  estimator_files <- list.files(estimator_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(estimator_files) == 0) {
    stop("No ", estimator, " RDS files found in ", estimator_dir,
         ". Run scripts/compute_meta_estimates.R first.")
  }

  estimator_files %>%
    map_dfr(readRDS) %>%
    mutate(
      sei = sqrt(vi),
      sse_yn = ifelse(!is.na(small_study_effect_pval) &
                        small_study_effect_pval <= 0.05, "yes", "no")
    )
}

add_power_variables <- function(dat, meta_average_multiplier = 0.5) {
  alpha <- 0.05
  q1 <- qnorm(1 - alpha / 2)
  q2 <- qnorm(alpha / 2)

  dat %>%
    mutate(
      GE = meta_average_multiplier * GE,
      lambda = abs(GE) / sei,
      power = 1 - pnorm(q1 - lambda) + pnorm(q2 - lambda),
      yn80 = ifelse(power >= 0.8, "yes", "no"),
      yn20 = ifelse(power <= 0.2, "yes", "no")
    )
}

split_meta_analyses <- function(dat) {
  split(dat, dat$cID)
}

make_grids <- function() {
  list(
    p_grid_tab = c(
      -Inf,
      qnorm(0.001 / 2, lower.tail = TRUE), qnorm(0.01 / 2, lower.tail = TRUE),
      qnorm(0.05 / 2, lower.tail = TRUE), qnorm(0.1 / 2, lower.tail = TRUE),
      qnorm(0.2 / 2, lower.tail = TRUE), qnorm(0.3 / 2, lower.tail = TRUE),
      qnorm(0.4 / 2, lower.tail = TRUE), qnorm(0.5 / 2, lower.tail = TRUE),
      qnorm(0.6 / 2, lower.tail = TRUE), qnorm(0.7 / 2, lower.tail = TRUE),
      qnorm(0.8 / 2, lower.tail = TRUE), qnorm(0.9 / 2, lower.tail = TRUE),
      0,
      qnorm(0.9 / 2, lower.tail = FALSE), qnorm(0.8 / 2, lower.tail = FALSE),
      qnorm(0.7 / 2, lower.tail = FALSE), qnorm(0.6 / 2, lower.tail = FALSE),
      qnorm(0.5 / 2, lower.tail = FALSE), qnorm(0.4 / 2, lower.tail = FALSE),
      qnorm(0.3 / 2, lower.tail = FALSE), qnorm(0.2 / 2, lower.tail = FALSE),
      qnorm(0.1 / 2, lower.tail = FALSE), qnorm(0.05 / 2, lower.tail = FALSE),
      qnorm(0.01 / 2, lower.tail = FALSE), qnorm(0.001 / 2, lower.tail = FALSE),
      Inf
    ),
    z_grid_plot = seq(-10.25, 10.25, 0.1025)
  )
}

run_parallel_cf <- function(dat, grid, heterogeneity_multiplier) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  cf(dat = dat, z.grid = grid, heterogeneity_multiplier = heterogeneity_multiplier)
}

run_parallel_cf_ci <- function(dat, grid, cluster, heterogeneity_multiplier) {
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  cf.ci.cluster(
    dat = dat,
    z.grid = grid,
    iters = n_iterations,
    cluster = cluster,
    heterogeneity_multiplier = heterogeneity_multiplier
  )
}

write_analysis_setup <- function(estimator_raw, grids, meta_average_multiplier,
                                 heterogeneity_multiplier, setup_label,
                                 estimator, outlier_variant = "outliers_removed",
                                 progress_bar = NULL, progress_offset = 0) {
  variant_suffix <- if (outlier_variant == "all_data") "_all_data" else ""
  result_suffix <- paste(setup_label, estimator, outlier_variant, sep = "_")
  estimator_power <- add_power_variables(estimator_raw, meta_average_multiplier)
  counterfactual_data <- estimator_raw %>% mutate(GE = meta_average_multiplier * GE)
  myDat_counterfactual <- split_meta_analyses(counterfactual_data)

  saveRDS(
    estimator_raw,
    here(
      "data", "derived_data",
      paste0("pps_rstandard_raw_", setup_label, "_", estimator, variant_suffix, ".rds")
    )
  )
  if (outlier_variant == "outliers_removed") {
    saveRDS(
      estimator_power,
      here(
        "data", "derived_data",
        paste0("pps_rstandard_power_", setup_label, "_", estimator, ".rds")
      )
    )
  }
  saveRDS(
    list(
      n_cores = n_cores,
      n_iterations = n_iterations,
      meta_average_multiplier = meta_average_multiplier,
      heterogeneity_multiplier = heterogeneity_multiplier,
      setup_label = setup_label,
      estimator = estimator,
      outlier_variant = outlier_variant
    ),
    here("data", "derived_data", paste0("analysis_settings_", result_suffix, ".rds"))
  )

  z_plot_path <- here("data", "derived_data", paste0("z_plot_", setup_label, "_", estimator, ".rds"))
  z_plot_ci_path <- here("data", "derived_data", paste0("z_plot_ci_", setup_label, "_", estimator, ".rds"))
  p_tab_path <- here("data", "derived_data", paste0("p_tab_", result_suffix, ".rds"))
  p_tab_ci_path <- here("data", "derived_data", paste0("p_tab_ci_", result_suffix, ".rds"))

  if (outlier_variant == "outliers_removed") {
    save_counterfactual(z_plot_path, function() run_parallel_cf(myDat_counterfactual, grids$z_grid_plot, heterogeneity_multiplier), progress_bar, progress_offset + 1)
    save_counterfactual(z_plot_ci_path, function() run_parallel_cf_ci(myDat_counterfactual, grids$z_grid_plot, unique(estimator_raw$cID), heterogeneity_multiplier), progress_bar, progress_offset + 2, n_iterations)
  } else {
    update_progress_bar(progress_bar, progress_offset + 1)
    update_progress_bar(progress_bar, progress_offset + 2)
  }
  save_counterfactual(p_tab_path, function() run_parallel_cf(myDat_counterfactual, grids$p_grid_tab, heterogeneity_multiplier), progress_bar, progress_offset + 3)
  save_counterfactual(p_tab_ci_path, function() run_parallel_cf_ci(myDat_counterfactual, grids$p_grid_tab, unique(estimator_raw$cID), heterogeneity_multiplier), progress_bar, progress_offset + 4, n_iterations)
}

ensure_output_dirs()

grids <- make_grids()
analysis_data_parameters <- tidyr::crossing(
  analysis_setups,
  estimator = meta_analysis_estimators,
  outlier_variant = "outliers_removed"
)

counterfactual_progress <- new_progress_bar(
  4 * nrow(analysis_data_parameters),
  "Counterfactual data progress"
)
for (parameter_index in seq_len(nrow(analysis_data_parameters))) {
  parameters <- analysis_data_parameters[parameter_index, ]
  estimator_raw <- load_estimator_data(parameters$estimator, parameters$outlier_variant)
  write_analysis_setup(
    estimator_raw, grids,
    parameters$meta_average_multiplier, parameters$heterogeneity_multiplier,
    parameters$setup_label, parameters$estimator, parameters$outlier_variant,
    counterfactual_progress, 4 * (parameter_index - 1)
  )
}
close_progress_bar(counterfactual_progress)
