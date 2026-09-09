## Shared packages, runtime settings, helper functions, and output directories.
## Define any of these settings before sourcing this file to override its default.

library(metafor); library(clubSandwich)
library(tidyverse); library(xtable)
library(foreach); library(doParallel); library(readxl); library(openxlsx)
library(MASS); library(car); library(lmtest); library(sandwich)
library(stargazer)
library(gridExtra); library(ggeasy)
library(orchaRd)
library(scales); library(here)

source(here("scripts", "runtime_settings.R"))

## Console progress bars use base R so they also work in a restored renv and in
## non-interactive batch jobs.

new_progress_bar <- function(total, label) {
  message(label)
  progress_bar <- utils::txtProgressBar(
    min = 0, max = max(1, total), initial = 0, style = 3,
    file = stderr()
  )
  flush.console()
  progress_bar
}

update_progress_bar <- function(progress_bar, value) {
  if (!is.null(progress_bar)) {
    utils::setTxtProgressBar(progress_bar, value)
    flush.console()
  }
  invisible(value)
}

close_progress_bar <- function(progress_bar) {
  if (!is.null(progress_bar)) {
    close(progress_bar)
    flush.console()
  }
  invisible(NULL)
}

## When the setup is sourced again in the same workflow, recover the configured
## multiplier vectors from the setup grid. Rendering scripts use the scalar
## names for the setup currently being rendered, so those names no longer
## contain the original configuration by the time a later script sources this
## file (for example, before creating the supplementary figures).
if (exists("analysis_setups") &&
    all(c("meta_average_multiplier", "heterogeneity_multiplier") %in%
        names(analysis_setups))) {
  meta_average_multiplier <- unique(analysis_setups$meta_average_multiplier)
  heterogeneity_multiplier <- unique(analysis_setups$heterogeneity_multiplier)
}

## Set defaults without replacing values supplied before source("scripts/main.R").
## This makes short test runs and custom setups possible from the calling script.
resolve_runtime_setting("n_cores", 6L)
resolve_runtime_setting("n_iterations", 1000L)
if (length(n_iterations) != 1L || is.na(n_iterations) ||
    n_iterations < 1 || n_iterations != as.integer(n_iterations)) {
  stop("n_iterations must be one positive whole number.")
}
n_iterations <- as.integer(n_iterations)
if (!exists("meta_average_multiplier")) {
  meta_average_multiplier <- c(0.25, 0.5, 1)
}
if (!exists("heterogeneity_multiplier")) {
  heterogeneity_multiplier <- c(0, 0.25, 0.5, 0.75)
}

make_setup_label <- function(meta_average_multiplier, heterogeneity_multiplier) {
  paste0(
    "meta_", gsub("\\.", "p", as.character(meta_average_multiplier)),
    "_heterogeneity_", gsub("\\.", "p", as.character(heterogeneity_multiplier))
  )
}

analysis_setups <- tidyr::expand_grid(
  meta_average_multiplier = meta_average_multiplier,
  heterogeneity_multiplier = heterogeneity_multiplier
) %>%
  dplyr::mutate(setup_label = make_setup_label(meta_average_multiplier, heterogeneity_multiplier))

## Estimators used throughout the manuscript-output workflow. The identifiers
## are also used as filename suffixes, so keep them filesystem friendly.
meta_analysis_estimators <- c("pet_peese", "multilevel_random")

source(here("scripts", "functions.R"))

output_dirs <- list(
  here("results", "main"),
  here("results", "supplement"),
  here("data", "derived_data")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
