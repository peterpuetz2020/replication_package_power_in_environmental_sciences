## -------
## main.R
## -------
## Orchestrate the analysis from the repository root. Data-generation scripts are
## optional because the repository already contains the derived data consumed by
## the table, figure, and regression scripts.

source(here::here("scripts", "analysis_setup.R"))

## Set to TRUE to overwrite and rebuild counterfactual simulation files.
recreate_counterfactuals <- FALSE

analysis_steps <- c(
  "Fit meta-analysis models",
  "Create analysis data",
  "Create main tables and figures",
  "Create supplementary outputs",
  "Create legacy and robustness outputs",
  "Run exploratory regressions"
)
workflow_progress <- new_progress_bar(length(analysis_steps), "Overall analysis progress")
workflow_step <- 0
run_analysis_step <- function(path) {
  source(here::here("scripts", path), local = parent.frame())
  workflow_step <<- workflow_step + 1
  update_progress_bar(workflow_progress, workflow_step)
}

## Optional: Fit PET-PEESE and multilevel random-effects models, save their effect-level
## analysis data, and write the meta-analysis estimates used below.
run_analysis_step("compute_meta_estimates.R")

## Build setup-specific datasets and any missing counterfactual simulations before
## rendering. Existing counterfactual files are reused unless explicitly rebuilt.
run_analysis_step("create_analysis_data.R")

## Recreate the ordered manuscript tables and figures from available data.
run_analysis_step("create_tables_and_figures.R")

## Recreate supplementary sensitivity and excess-significance figures.
run_analysis_step("create_supplement_tables_and_figures.R")

## Optional: recreate the original counterfactual, subfield, and robustness data
## and render the complete legacy collection of tables and figures.
run_analysis_step("create_full_tables_figures_and_robustness.R")

## Optional: after running the complete legacy workflow above, recreate the ESR
## workbook used by the regression analysis.
# source(here::here("scripts", "create_esr_data.R"))

## Recreate exploratory regression models and diagnostic figures. This script
## loads the supplied ESR and power workbooks when regenerated copies are absent.
run_analysis_step("run_exploratory_regressions.R")
close_progress_bar(workflow_progress)
