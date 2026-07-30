## Shared packages, runtime settings, helper functions, and output directories.
rm(list = ls(all = TRUE))

library(metafor); library(clubSandwich)
library(tidyverse); library(xtable)
library(foreach); library(doParallel); library(readxl); library(openxlsx)
library(MASS); library(car); library(lmtest); library(sandwich)
library(stargazer)
library(gridExtra); library(ggeasy)
library(gt); library(gtExtras); library(orchaRd)
library(scales); library(here)

n_cores <- 7
n_iterations <- 1000
meta_average_multiplier <- 0.5
heterogeneity_multiplier <- 0.25

source(here("scripts", "functions.R"))

output_dirs <- list(
  here("results", "main"),
  here("results", "main", "pet_peese_rstandard"),
  here("results", "robustness")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
