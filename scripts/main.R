## -------
## main.R
## -------
## This entry point only sets global parameters and sources the scripts needed to
## reproduce the analysis. Uncomment the optional scripts below if the underlying
## analysis data or journal classification files need to be rebuilt.

## User-adjustable analysis settings.
## n_cores controls parallel analyses: use 1 on low-memory laptops, 2-4 on
## typical 4-8 core laptops/desktops, and 6-8 or more on workstations/HPC nodes.
## Keep at least one core free for the operating system.
n_cores <- 7

## Number of Monte Carlo/bootstrap iterations for confidence intervals. The paper
## uses 1000; smaller values are useful only for quick code checks.
n_iterations <- 1000

## Multiplier applied to genuine-effect estimates before power analyses. The
## default 0.5 uses half of the meta-average; change to 1 for the full
## meta-average or 0.25 for one fourth of the meta-average.
meta_average_multiplier <- 0.5

## Multiplier applied to the between-study variance (tau^2) in counterfactual
## z-/p-value simulations. The manuscript sensitivity analyses use 0, 0.25,
## and 0.5; set one value here to choose the degree of heterogeneity.
heterogeneity_multiplier <- 0.25

## Label used when reading/writing setup-specific intermediate outputs.
setup_label <- "half"

## Set to TRUE only when the counterfactual z-/p-value files should be rebuilt.
## These steps can be very time consuming with n_iterations <- 1000.
recreate_counterfactuals <- FALSE

## Optional: rebuild the data needed to recreate tables and figures.
# source(here::here("scripts", "create_analysis_data.R"))

## Optional: rebuild the journal classification file.
# source(here::here("scripts", "classification.R"))

## Always recreate manuscript tables and figures from the available data.
source(here::here("scripts", "create_tables_and_figures.R"))
