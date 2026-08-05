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

## Console progress bars use base R so they also work in a restored renv and in
## non-interactive batch jobs. Callers can set show_progress <- FALSE before
## sourcing a script when its output is being redirected to a log.
if (!exists("show_progress")) show_progress <- TRUE

new_progress_bar <- function(total, label) {
  if (!isTRUE(show_progress)) return(NULL)
  message(label)
  utils::txtProgressBar(min = 0, max = max(1, total), style = 3)
}

update_progress_bar <- function(progress_bar, value) {
  if (!is.null(progress_bar)) utils::setTxtProgressBar(progress_bar, value)
  invisible(value)
}

close_progress_bar <- function(progress_bar) {
  if (!is.null(progress_bar)) close(progress_bar)
  invisible(NULL)
}

## Set analysis parameters here. All downstream scripts consume these settings.
n_cores <- 4
n_iterations <- 10
meta_average_multiplier <- c(0.25, 0.5, 1)
heterogeneity_multiplier <- c(0, 0.25, 0.5, 0.75)

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
  here("results", "main", "pet_peese_rstandard"),
  here("results", "main", "multilevel_random"),
  here("results", "robustness")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
