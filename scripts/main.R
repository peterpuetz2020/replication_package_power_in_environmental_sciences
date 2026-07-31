## -------
## main.R
## -------
## Orchestrate the analysis from the repository root. Data-generation scripts are
## optional because the repository already contains the derived data consumed by
## the table, figure, and regression scripts.

source(here::here("scripts", "analysis_setup.R"))

## Optional: refit PET-PEESE models and replace the per-meta-analysis RDS files.
# source(here::here("scripts", "fit_pet_peese_models.R"))

## Build setup-specific datasets and any missing counterfactual simulations before
## rendering. Existing counterfactual files are reused unless explicitly rebuilt.
source(here::here("scripts", "create_analysis_data.R"))

## Recreate the ordered manuscript tables and figures from available data.
source(here::here("scripts", "create_tables_and_figures.R"))

## Optional: recreate the original counterfactual, subfield, and robustness data
## and render the complete legacy collection of tables and figures.
source(here::here("scripts", "create_full_tables_figures_and_robustness.R"))

## Optional: after running the complete legacy workflow above, recreate the ESR
## workbook used by the regression analysis.
# source(here::here("scripts", "create_esr_data.R"))

## Recreate exploratory regression models and diagnostic figures. This script
## loads the supplied ESR and power workbooks when regenerated copies are absent.
source(here::here("scripts", "run_exploratory_regressions.R"))
