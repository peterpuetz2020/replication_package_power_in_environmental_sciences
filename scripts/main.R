## -------
## main.R
## -------
## Orchestrate the complete analysis from the repository root. The two slow raw-
## data preparation stages are controlled by explicit switches and are disabled
## by default. Define either switch before sourcing main.R to override it.
run_classification <- if (exists("run_classification")) run_classification else FALSE
run_meta_analysis_fitting <- if (exists("run_meta_analysis_fitting")) {
  run_meta_analysis_fitting
} else {
  FALSE
}
recreate_counterfactuals <- if (exists("recreate_counterfactuals")) {
  recreate_counterfactuals
} else {
  FALSE
}
recreate_meta_analysis_estimates <- if (exists("recreate_meta_analysis_estimates")) {
  recreate_meta_analysis_estimates
} else {
  FALSE
}

## Recreate the article/subfield classification only when explicitly requested.
## The supplied classification CSV files are used on ordinary replication runs.
if (run_classification) {
  source(here::here("scripts", "classification.R"))
}

source(here::here("scripts", "analysis_setup.R"))

## Fitting all 708 meta-analyses is the longest stage. Run it when explicitly
## requested or automatically when any required estimator directory has no RDS
## inputs. The fitting caches make an interrupted automatic run resumable.
estimator_input_dirs <- c(
  here::here("data", "derived_data", "pet_peese_rstandard"),
  here::here("data", "derived_data", "pet_peese_all_data"),
  here::here("data", "derived_data", "multilevel_random"),
  here::here("data", "derived_data", "multilevel_random_all_data")
)
missing_estimator_inputs <- estimator_input_dirs[vapply(
  estimator_input_dirs,
  function(path) {
    !dir.exists(path) || length(list.files(path, pattern = "\\.rds$")) == 0
  },
  logical(1)
)]

if (length(missing_estimator_inputs) > 0 && !run_meta_analysis_fitting) {
  message(
    "Required estimator RDS files are missing from: ",
    paste(missing_estimator_inputs, collapse = ", "),
    ". Running scripts/compute_meta_estimates.R automatically."
  )
}

if (run_meta_analysis_fitting || length(missing_estimator_inputs) > 0) {
  source(here::here("scripts", "compute_meta_estimates.R"))
}

## Build setup-specific datasets and any missing counterfactual simulations before
## rendering. Existing counterfactual files are reused unless explicitly rebuilt.
source(here::here("scripts", "create_analysis_data.R"))

## Fit the negative-binomial models needed by Table 4. The fitted objects are
## consumed by the main and supplementary output scripts below.
source(here::here("scripts", "run_exploratory_regressions.R"))

## Recreate the ordered manuscript tables and figures from available data.
source(here::here("scripts", "create_tables_and_figures.R"))

## Recreate supplementary sensitivity and excess-significance figures.
source(here::here("scripts", "create_supplement_tables_and_figures.R"))
