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
      sse_yn = case_when(
        is.na(small_study_effect_pval) ~ NA_character_,
        small_study_effect_pval <= 0.05 ~ "yes",
        TRUE ~ "no"
      )
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

write_analysis_setup <- function(estimator_raw, grids, meta_average_multiplier,
                                 heterogeneity_multiplier, setup_label,
                                 estimator, outlier_variant = "outliers_removed",
                                 progress_bar = NULL, progress_offset = 0) {
  variant_suffix <- if (outlier_variant == "all_data") "_all_data" else ""
  result_suffix <- paste(setup_label, estimator, outlier_variant, sep = "_")
  estimator_power <- add_power_variables(estimator_raw, meta_average_multiplier)
  counterfactual_data <- estimator_raw %>%
    mutate(GE = meta_average_multiplier * GE) %>%
    filter_counterfactual_data(
      heterogeneity_multiplier,
      context = paste(setup_label, estimator, outlier_variant)
    )
  myDat_counterfactual <- split_meta_analyses(counterfactual_data)
  component_cache <- new.env(parent = emptyenv())
  get_components <- function(grid_name, grid) {
    if (!exists(grid_name, envir = component_cache, inherits = FALSE)) {
      assign(
        grid_name,
        counterfactual_components(
          myDat_counterfactual, grid, heterogeneity_multiplier
        ),
        envir = component_cache
      )
    }
    get(grid_name, envir = component_cache, inherits = FALSE)
  }

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
    save_counterfactual(z_plot_path, function() cf(
      myDat_counterfactual, grids$z_grid_plot, heterogeneity_multiplier,
      components = get_components("z", grids$z_grid_plot)
    ), progress_bar, progress_offset + 1)
    save_counterfactual(z_plot_ci_path, function() cf.ci.cluster(
      myDat_counterfactual, grids$z_grid_plot, n_iterations,
      unique(counterfactual_data$cID), heterogeneity_multiplier,
      components = get_components("z", grids$z_grid_plot)
    ), progress_bar, progress_offset + 2, n_iterations)
  } else {
    update_progress_bar(progress_bar, progress_offset + 1)
    update_progress_bar(progress_bar, progress_offset + 2)
  }
  save_counterfactual(p_tab_path, function() cf(
    myDat_counterfactual, grids$p_grid_tab, heterogeneity_multiplier,
    components = get_components("p", grids$p_grid_tab)
  ), progress_bar, progress_offset + 3)
  save_counterfactual(p_tab_ci_path, function() cf.ci.cluster(
    myDat_counterfactual, grids$p_grid_tab, n_iterations,
    unique(counterfactual_data$cID), heterogeneity_multiplier,
    components = get_components("p", grids$p_grid_tab)
  ), progress_bar, progress_offset + 4, n_iterations)
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
counterfactual_cluster <- makeCluster(n_cores)
registerDoParallel(counterfactual_cluster)
estimator_data_cache <- new.env(parent = emptyenv())
tryCatch(
  {
    for (parameter_index in seq_len(nrow(analysis_data_parameters))) {
      parameters <- analysis_data_parameters[parameter_index, ]
      cache_key <- paste(
        parameters$estimator, parameters$outlier_variant, sep = "_"
      )
      if (!exists(cache_key, envir = estimator_data_cache, inherits = FALSE)) {
        assign(
          cache_key,
          load_estimator_data(
            parameters$estimator, parameters$outlier_variant
          ),
          envir = estimator_data_cache
        )
      }
      estimator_raw <- get(
        cache_key, envir = estimator_data_cache, inherits = FALSE
      )
      write_analysis_setup(
        estimator_raw, grids,
        parameters$meta_average_multiplier, parameters$heterogeneity_multiplier,
        parameters$setup_label, parameters$estimator,
        parameters$outlier_variant,
        counterfactual_progress, 4 * (parameter_index - 1)
      )
    }
  },
  finally = {
    stopCluster(counterfactual_cluster)
    close_progress_bar(counterfactual_progress)
  }
)

## Derive the multilevel random-effects regression inputs consumed by Table 4.
## Keep these outputs with the other generated analysis data rather than in the
## manuscript-output directories.
regression_data_dir <- here("data", "derived_data", "regression_data")
dir.create(regression_data_dir, recursive = TRUE, showWarnings = FALSE)

regression_setup_label <- "meta_0p5_heterogeneity_0"
regression_estimator <- "multilevel_random"
random_effects_path <- here(
  "data", "derived_data",
  paste0(
    "pps_rstandard_raw_", regression_setup_label, "_",
    regression_estimator, ".rds"
  )
)
random_effects_data <- readRDS(random_effects_path) %>%
  add_power_variables(0.5)

critical_value <- qnorm(0.975)
esr_data <- random_effects_data %>%
  mutate(
    observed_z = abs(yi / sei),
    expected_significant_probability =
      pnorm(-critical_value, mean = GE / sei) +
      pnorm(critical_value, mean = GE / sei, lower.tail = FALSE)
  ) %>%
  group_by(cID) %>%
  summarise(
    tot.all = n(),
    tot.sig = sum(observed_z >= critical_value),
    expected.sig = sum(expected_significant_probability),
    esr.sig.count = tot.sig - expected.sig,
    esr.all.count = esr.sig.count,
    esr.sig = if_else(tot.sig > 0, esr.sig.count / tot.sig, 0),
    esr.all = esr.all.count / tot.all,
    dif = esr.sig.count,
    .groups = "drop"
  ) %>%
  dplyr::select(cID, esr.all, esr.sig, esr.all.count, esr.sig.count,
    dif, tot.all, tot.sig)
saveRDS(
  esr_data,
  file.path(regression_data_dir, "esr05_multilevel_random.rds")
)

## Create the per-meta-analysis power workbook used by the exploratory
## regressions. Publication year comes from the primary input data; the
## separately supplied five-year journal impact factors are joined by cID.
regression_power_data <- random_effects_data %>%
  group_by(cID) %>%
  summarise(
    metaID = first(metaID),
    median = median(power, na.rm = TRUE),
    sape = mean(power >= 0.8, na.rm = TRUE),
    nips = n_distinct(sID),
    esty = first(etype),
    guid = first(guide),
    prer = first(prere),
    subf = first(subfd),
    sdes = first(sdesn),
    .groups = "drop"
  ) %>%
  mutate(cID = as.character(cID))

journal_impact_factors <- readxl::read_excel(
  here("data", "journal_impact_factors.xlsx"))
 
regression_power_data <- regression_power_data %>%
    left_join(journal_impact_factors, by = c("metaID", "cID"))  

if (anyNA(regression_power_data$pyear) ||
    anyNA(regression_power_data$jif_5yr_wos)) {
  stop(
    "Publication year or journal impact factor is missing for one or more ",
    "meta-analyses."
  )
}

power_summary_path <- file.path(
  regression_data_dir,
  "power_summary_multilevel_random_half_meta_average.xlsx"
)
openxlsx::write.xlsx(
  regression_power_data, power_summary_path, overwrite = TRUE
)

regression_covariates <- regression_power_data %>%
  dplyr::select(cID, pyear, jif_5yr_wos)
saveRDS(
  regression_covariates,
  file.path(regression_data_dir, "regression_covariates.rds")
)
