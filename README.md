# Replication package: power in environmental sciences

Replication materials for **“Low statistical power and overrepresentation of
statistically significant findings in the environmental sciences”** by Teshome
K. Deressa, Peter Pütz, David I. Stern, Jaco Vangronsveld, Jan Minx, Sebastien
Lizin, Robert Malina, Stephan B. Bruns.

The repository contains the input data, R environment, analysis scripts, and
manuscript outputs needed to reproduce the paper and Supporting Information.

## Requirements

- R 4.5.0 (the version recorded in `renv.lock`)
- An internet connection for the initial package restore
- RStudio is optional; open `power and bias.Rproj` when using it
- On Windows, Rtools45 may be needed if packages must be compiled from source

Start R in the repository root and restore the project library:

```r
install.packages("renv")
renv::restore()
```

The lockfile pins the CRAN dependencies and the GitHub package `orchaRd`.

## Reproduce the results

Run the complete workflow from the repository root:

```r
source("scripts/main.R")
```

`main.R` performs the steps in dependency order:

1. optionally recreates the subfield classifications;
2. fits fixed-effects, PET-PEESE and multilevel random-effects estimates when their derived
   RDS inputs are missing;
3. creates setup-specific analysis data and counterfactual simulations;
4. runs the regression analyses;
5. writes the main and supplementary tables and figures.

The supplied meta-estimate RDS files make refitting optional. To rerun that
long stage, set `run_meta_analysis_fitting <- TRUE`; completed meta-analyses are
cached individually so an interrupted fit can resume.

### Force a complete rebuild

The supplied classification CSV files are used by default. To recreate them and
refit every meta-analysis from the original inputs, run:

```r
run_classification <- TRUE
run_meta_analysis_fitting <- TRUE
recreate_meta_analysis_estimates <- TRUE
recreate_counterfactuals <- TRUE
source("scripts/main.R")
```

### Run individual stages

After restoring the environment, stages can also be sourced separately in this
order:

```r
source("scripts/classification.R")                    # optional
source("scripts/compute_meta_estimates.R")            # model fitting
source("scripts/create_analysis_data.R")              # derived analysis data
source("scripts/run_exploratory_regressions.R")       # regression models
source("scripts/create_tables_and_figures.R")         # main outputs
source("scripts/create_supplement_tables_and_figures.R") # supplement
```

The rendering scripts can reuse existing derived files. If the required files
are absent, use `scripts/main.R`, which runs the prerequisite stages.

## Analysis options

Define options before sourcing `scripts/main.R` or
`scripts/create_analysis_data.R`:

| Option | Default | Purpose |
|---|---:|---|
| `n_cores` | `6` | Parallel workers. Leave at least one CPU core free. |
| `n_iterations` | `1000` | Bootstrap iterations for confidence intervals. |
| `meta_average_multiplier` | `c(0.25, 0.5, 1)` | Multipliers applied to the meta-analytic average in power calculations. |
| `heterogeneity_multiplier` | `c(0, 0.25, 0.5, 0.75)` | Multipliers applied to between-effect heterogeneity in counterfactual calculations. |
| `run_classification` | `FALSE` | Recreate the classification CSV files. |
| `run_meta_analysis_fitting` | `FALSE` | Run model fitting even when estimator inputs exist. |
| `recreate_meta_analysis_estimates` | `FALSE` | Ignore fitting caches and refit all meta-analyses. |
| `recreate_counterfactuals` | `FALSE` | Replace matching counterfactual caches. |

For example, a small exploratory run can use:

```r
n_cores <- 2
n_iterations <- 100
meta_average_multiplier <- 0.5
heterogeneity_multiplier <- c(0, 0.5)
source("scripts/main.R")
```

Use `n_iterations <- 1000` for the manuscript results. Both multiplier options
accept vectors; the data-creation workflow evaluates every combination.
Custom setup names can be supplied through an `analysis_setups` data frame with
`meta_average_multiplier`, `heterogeneity_multiplier`, and unique `setup_label`
columns.

To reproduce all main-text outputs, retain the manuscript’s `0.5` meta-average 
setups with heterogeneity multipliers `0` and `0.5`.

## Repository structure

```text
data/
  MasterData.xlsx                         main effect-size data
  meta-articles.xlsx                      article metadata
  meta-classified into subfields.csv      supplied classifications
  ref_env_percentage.csv                  supplied classification input
  journal_impact_factors.xlsx             journal metadata
  scimago_sjr.xlsx                        classification source data
scripts/
  main.R                                  workflow orchestrator
  analysis_setup.R                        packages, settings, and shared setup
  runtime_settings.R                      handling of option defaults
  classification.R                        optional classification recreation
  compute_meta_estimates.R                fixed-effects, PET-PEESE and random-effects fitting
  create_analysis_data.R                  setup-specific derived data
  run_exploratory_regressions.R           regression analyses
  create_tables_and_figures.R             main manuscript outputs
  create_supplement_tables_and_figures.R  supplementary outputs
  functions.R                             shared analysis functions
results/
  main/                                   main tables and figures
  supplement/                             supplementary tables and figures
renv.lock                                 pinned R package environment
```

Generated files are written to `data/derived_data/`. Only the supplied
meta-estimate RDS files in its four estimator directories are versioned; all
other derived files are ignored and recreated as needed.


## Outputs

Main outputs are written to `results/main/` and Supporting Information outputs
to `results/supplement/`. Tables are produced in editable manuscript formats
(including DOCX, TeX, or XLSX as applicable), and figures are rendered as
scalable PDF and SVG files. Generated PDFs are excluded from version control but
are recreated by the scripts.

Setup-specific intermediate filenames contain their multiplier choices, for
example `meta_0p5_heterogeneity_0p25`. Cache keys include the data, grid,
settings, bootstrap clusters, and iteration count, so incompatible cached
counterfactual results are rebuilt automatically.
