## -------
## main.R
## -------
## Orchestrate the analysis from the repository root. Data-generation scripts are
## optional because the repository already contains the derived data consumed by
## the table, figure, and regression scripts.

source(here::here("scripts", "analysis_setup.R"))

## Set to TRUE to overwrite and rebuild counterfactual simulation files.
recreate_counterfactuals <- FALSE

## Optional: Fit PET-PEESE and multilevel random-effects models, save their effect-level
## analysis data, and write the meta-analysis estimates used below.
source(here::here("scripts", "compute_meta_estimates.R"))

## Build setup-specific datasets and any missing counterfactual simulations before
## rendering. Existing counterfactual files are reused unless explicitly rebuilt.
source(here::here("scripts", "create_analysis_data.R"))

## Fit the negative-binomial models needed by Table 4 and save their diagnostic
## plots. The fitted objects are consumed by create_tables_and_figures.R.
source(here::here("scripts", "run_exploratory_regressions.R"))

## Recreate the ordered manuscript tables and figures from available data.
source(here::here("scripts", "create_tables_and_figures.R"))

## Recreate supplementary sensitivity and excess-significance figures.
source(here::here("scripts", "create_supplement_tables_and_figures.R"))

## Optional: recreate the original counterfactual, subfield, and robustness data
## and render the complete legacy collection of tables and figures.
source(here::here("scripts", "create_full_tables_figures_and_robustness.R"))
