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

## Fitting all 708 meta-analyses is the longest stage. Leave this disabled when
## the per-meta-analysis files already exist under data/derived_data, and enable
## it for a clean, from-source rebuild.
if (run_meta_analysis_fitting) {
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
