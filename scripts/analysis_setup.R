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

n_cores <- if (exists("n_cores", inherits = FALSE)) n_cores else 7
n_iterations <- if (exists("n_iterations", inherits = FALSE)) n_iterations else 1000
meta_average_multiplier <- if (exists("meta_average_multiplier", inherits = FALSE)) meta_average_multiplier else c(0.5, 1)
heterogeneity_multiplier <- if (exists("heterogeneity_multiplier", inherits = FALSE)) heterogeneity_multiplier else c(0.25, 0.5, 0.75)

source(here("scripts", "functions.R"))

output_dirs <- list(
  here("results", "main"),
  here("results", "main", "pet_peese_rstandard"),
  here("results", "robustness")
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
