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

## Missing counterfactual z-/p-value files are always created because the tables
## and plots consume them. Set this to TRUE to overwrite existing matching files.
## These steps can be very time consuming with n_iterations <- 1000.
recreate_counterfactuals <- if (exists("recreate_counterfactuals")) recreate_counterfactuals else FALSE

save_counterfactual <- function(path, calculate) {
  if (recreate_counterfactuals || !file.exists(path)) {
    saveRDS(calculate(), path)
  }
}

ensure_output_dirs <- function() {
  invisible(lapply(
    list(
      here("results", "main"),
      here("results", "main", "pet_peese_rstandard"),
      here("results", "main", "derived_data")
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

load_pet_peese_data <- function() {
  pet_peese_files <- list.files(
    here("results", "main", "pet_peese_rstandard"),
    pattern = "\\.rds$",
    full.names = TRUE
  )

  if (length(pet_peese_files) == 0) {
    stop(
      "No PET-PEESE RDS files found in results/main/pet_peese_rstandard. ",
      "Run the PET-PEESE estimation block in scripts/main.R first."
    )
  }

  pet_peese_files %>%
    map_dfr(readRDS) %>%
    mutate(
      sei = sqrt(vi),
      sse_yn = ifelse(small_study_effect_pval <= 0.05, "yes", "no")
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

write_analysis_setup <- function(pps_rstandard_raw, grids, meta_average_multiplier, heterogeneity_multiplier, setup_label) {
  pps_rstandard_power <- add_power_variables(pps_rstandard_raw, meta_average_multiplier)
  pps_counterfactual <- pps_rstandard_raw %>% mutate(GE = meta_average_multiplier * GE)
  myDat_counterfactual <- split_meta_analyses(pps_counterfactual)

  saveRDS(
    pps_rstandard_raw,
    here("results", "main", "derived_data", paste0("pps_rstandard_raw_", setup_label, ".rds"))
  )
  saveRDS(
    pps_rstandard_power,
    here("results", "main", "derived_data", paste0("pps_rstandard_power_", setup_label, ".rds"))
  )
  saveRDS(
    list(
      n_cores = n_cores,
      n_iterations = n_iterations,
      meta_average_multiplier = meta_average_multiplier,
      heterogeneity_multiplier = heterogeneity_multiplier,
      setup_label = setup_label
    ),
    here("results", "main", "derived_data", paste0("analysis_settings_", setup_label, ".rds"))
  )

  z_plot_path <- here("results", "main", paste0("z_plot_pet_peese_rstandard_", setup_label, ".rds"))
  z_plot_ci_path <- here("results", "main", paste0("z_plot_ci_pet_peese_rstandard_", setup_label, ".rds"))
  p_tab_path <- here("results", "main", paste0("p_tab_pps_rstandard_", setup_label, ".rds"))
  p_tab_ci_path <- here("results", "main", paste0("p_tab_ci_pps_rstandard_", setup_label, ".rds"))

  save_counterfactual(z_plot_path, function() run_parallel_cf(myDat_counterfactual, grids$z_grid_plot, heterogeneity_multiplier))
  save_counterfactual(z_plot_ci_path, function() run_parallel_cf_ci(myDat_counterfactual, grids$z_grid_plot, unique(pps_rstandard_raw$cID), heterogeneity_multiplier))
  save_counterfactual(p_tab_path, function() run_parallel_cf(myDat_counterfactual, grids$p_grid_tab, heterogeneity_multiplier))
  save_counterfactual(p_tab_ci_path, function() run_parallel_cf_ci(myDat_counterfactual, grids$p_grid_tab, unique(pps_rstandard_raw$cID), heterogeneity_multiplier))
}

ensure_output_dirs()

pps_rstandard_raw <- load_pet_peese_data()
grids <- make_grids()

analysis_setups %>%
  pwalk(function(meta_average_multiplier, heterogeneity_multiplier, setup_label, ...) {
    write_analysis_setup(pps_rstandard_raw, grids, meta_average_multiplier, heterogeneity_multiplier, setup_label)
  })
