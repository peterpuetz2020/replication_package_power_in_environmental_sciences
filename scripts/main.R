## -------
## main.R
## -------
## Orchestrate the analysis from the repository root. Data-generation scripts are
## optional because the repository already contains the derived data consumed by
## the table, figure, and regression scripts.

source(here::here("scripts", "analysis_setup.R"))

## Optional: Fit PET-PEESE and multilevel random-effects models, save their effect-level
## analysis data, and write the meta-analysis estimates used below.
source(here::here("scripts", "compute_meta_estimates.R"))

## Build setup-specific datasets and any missing counterfactual simulations before
## rendering. Existing counterfactual files are reused unless explicitly rebuilt.
source(here::here("scripts", "create_analysis_data.R"))

## Recreate the ordered manuscript tables and figures from available data.
source(here::here("scripts", "create_tables_and_figures.R"))

## Recreate supplementary sensitivity and excess-significance figures.
source(here::here("scripts", "create_supplement_tables_and_figures.R"))

## Optional: recreate the original counterfactual, subfield, and robustness data
## and render the complete legacy collection of tables and figures.
source(here::here("scripts", "create_full_tables_figures_and_robustness.R"))

## Optional: after running the complete legacy workflow above, recreate the ESR
## workbook used by the regression analysis.
# source(here::here("scripts", "create_esr_data.R"))

## Recreate exploratory regression models and diagnostic figures. This script
## loads the supplied ESR and power workbooks when regenerated copies are absent.
source(here::here("scripts", "run_exploratory_regressions.R"))
