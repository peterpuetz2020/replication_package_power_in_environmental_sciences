# Teshome K. Deressa, Peter Pütz, David I. Stern, Jaco Vangronsveld, Jan Minx, Sebastien Lizin, Robert Malina, Stephan B. Bruns

## Overview
This is the replication package of “Low statistical power and overrepresentation of statistically
significant findings in the environmental sciences” (including the Supporting Information).

The repository contains the data and R scripts needed to reproduce the subfield classification, tables, figures, and robustness checks.

## Computational environment

The project uses [`renv`](https://rstudio.github.io/renv/) to restore the package versions recorded in `renv.lock`. The lockfile records **R 4.5.0** and the CRAN repository snapshot configured through Posit Package Manager. Most packages are installed from CRAN; the non-CRAN dependency `orchaRd` is pinned in `renv.lock` to the GitHub repository `daniel1noble/orchaRd` at commit `5e9ac55cd28d717681bcbcf17527bace42500903`, so `renv::restore()` can install it reproducibly.

Chrome or Chromium is **not required**. The analysis writes tables directly as CSV/XLSX files and renders its PNG/PDF outputs through R graphics devices (including the heterogeneity summary table); it does not take browser screenshots. Consequently, the browser automation packages `webshot2`, `chromote`, and `websocket` are not included in the lockfile. Some installed packages list browser tooling only as an optional suggested dependency, but none of the replication scripts use that feature.

### RStudio

RStudio is optional but recommended for interactive reproduction. Install a current release of RStudio Desktop after installing R 4.5.0. Then open the repository's RStudio project file, `power and bias.Rproj`, so that paths resolve from the repository root and `renv` activates automatically.

Without RStudio, open R in the repository root before running the commands below. The `.Rprofile` file activates `renv` when R starts in this project.

### Rtools on Windows

Rtools is not used directly by these analysis scripts. On Windows, it may still be needed while restoring packages if R has to compile packages from source, including GitHub or CRAN packages without matching binaries. For R 4.5.0, install **Rtools45** if compilation is required.

## Restoring packages with renv

1. Install R 4.5.0.
2. On Windows, install Rtools45 if `renv::restore()` needs to compile packages from source.
3. Open the project:
   - With RStudio: open `power and bias.Rproj`.
   - Without RStudio: start R from the repository root.
4. Restore the project library:

```r
install.packages("renv")
renv::restore()
```

## Reproducing the results

Run the scripts from the repository root in this order for the complete legacy workflow:

```r
source("scripts/classification.R")
source("scripts/main.R")
```

Alternatively, open `scripts/classification.R` and run it first, then open `scripts/main.R` and run it second. In RStudio, this can be done by opening each file and choosing **Source**. Without RStudio, paste or source the same commands in an R session started from the repository root.

### Ordered manuscript tables and figures

To recreate the manuscript tables and figures in publication order, use the dedicated script:

```r
source("scripts/create_tables_and_figures.R")
```

This script creates outputs in numeric order: Table 1, Figure 1, Table 2, Table 3, and Figure 2. Each table/figure section reloads the data, settings, and grid definitions it needs, so a single section can be run independently in a fresh R session after the helper setup at the beginning of the file has been sourced. Setup-specific filenames are generated from the selected multipliers (for example, `meta_0p5_heterogeneity_0p25`). Figure 1 and Table 2 are written for every combination of `meta_average_multiplier` and `heterogeneity_multiplier`. Because heterogeneity is not used by Table 1, Table 3, or Figure 2, those outputs are written once per meta-average multiplier, using the first heterogeneity value only. A supplied `analysis_setups` data frame can instead assign a custom, unique `setup_label` to each combination.

The workflow renders each setup for PET-PEESE and multilevel random effects.
Before the first run, create the random-effects estimator datasets:

```r
source("scripts/compare_meta_analysis_estimators.R")
source("scripts/create_tables_and_figures.R")
```

The comparison script applies the PET residual screen consistently, then saves
effect-level inputs under `results/main/multilevel_random/`. Output filenames
end in `_pet_peese` or `_multilevel_random`, so results from each estimator
remain separate.

Publication outputs use the consistent names `Table_<number>_<setup_label>.<ext>` and `Figure_<number>_<setup_label>.<ext>`. Figures are saved in both PDF and EPS format. The workbook underlying Figure 2 is named `Figure_2_data_<setup_label>.xlsx`. Supplementary analyses use separately numbered, descriptive `Robustness_Table_<number>_...` and `Robustness_Figure_<number>_...` filenames.

### Optional setup-specific data creation

The data-preparation step is separated from table/figure rendering. To generate one or more alternative setups, define `meta_average_multipliers`, `heterogeneity_multipliers`, or `n_iterations` before running:

```r
source("scripts/create_analysis_data.R")
```

This writes setup-specific derived datasets under `results/main/derived_data/`. `meta_average_multipliers` and `heterogeneity_multipliers` can each contain one or more values; `scripts/create_analysis_data.R` computes and stores outputs for every combination. To use custom labels, define an `analysis_setups` tibble/data frame with `meta_average_multiplier`, `heterogeneity_multiplier`, and `setup_label` columns before sourcing the script. The script also creates any missing setup-specific counterfactual z-value and p-value RDS files needed by the tables and plots, while reusing files that already exist. Define `recreate_counterfactuals <- TRUE` before sourcing the script to overwrite and rebuild all matching counterfactual files. `scripts/create_tables_and_figures.R` uses these derived datasets when available; otherwise, it falls back to the existing PET-PEESE RDS files under `results/main/pet_peese_rstandard/`.

### Optional data-recreation step

The optional model-fitting script, `scripts/fit_pet_peese_models.R`, reads `data/MasterData.xlsx`. Running `scripts/classification.R` is only necessary if you want to recreate the classification files from the raw Scopus and Scimago inputs. `scripts/main.R` is now an orchestrator: it loads shared setup, renders outputs from the supplied derived data, and runs the exploratory regressions. Expensive scripts that replace supplied derived data are listed as commented, optional `source()` calls in `scripts/main.R`.

```r
source("scripts/main.R")
```

## Runtime settings in `scripts/analysis_setup.R`

In `scripts/analysis_setup.R`, adjust:

- `n_cores`: number of parallel worker cores.
- `n_iterations`: number of Monte Carlo/bootstrap iterations for confidence intervals.
- `meta_average_multiplier`: vector of multipliers applied to the meta-analytic average when calculating power (default `c(0.5, 1)`).
- `heterogeneity_multiplier`: vector of multipliers applied to the between-effect heterogeneity in the counterfactual calculations (default `c(0.25, 0.5, 0.75)`).

The two multiplier vectors are independent. Setup-specific data generation, Figure 1, and Table 2 evaluate every combination of their values. Table 1, Table 3, and Figure 2 are saved once for each meta-average multiplier, using the first supplied heterogeneity multiplier in their filenames. Define custom vectors before sourcing `scripts/main.R`:

```r
meta_average_multiplier <- c(0.5, 0.75, 1)
heterogeneity_multiplier <- c(0.25, 0.5)
source("scripts/main.R")
```

The paper uses `n_iterations <- 1000`. Smaller values are useful for quick checks only and should not be used for final replication.

Useful `n_cores` choices depend on the computer:

- `n_cores <- 1`: low-memory laptops or troubleshooting.
- `n_cores <- 2` to `4`: typical 4- to 8-core laptops/desktops.
- `n_cores <- 6` to `8`: stronger desktops or workstations.
- Higher values: high-performance computing nodes, if enough memory is available.

Leave at least one core free for the operating system and other applications.

## Output folders

All generated results are saved under `results/`:

- `results/main/` contains the main paper tables, figures, and intermediate result files.
- `results/robustness/` contains robustness-check tables, figures, diagnostics, and intermediate result files.

The analysis scripts create their output folders automatically if they do not already exist.

## Repository structure

- `data/`: `MasterData.xlsx` for the main analysis, plus raw and derived files used by the optional classification workflow.
- `scripts/classification.R`: optional script that recreates subfield classifications from raw inputs.
- `scripts/functions.R`: helper functions used by the analysis.
- `scripts/analysis_setup.R`: shared packages, runtime settings, helpers, and output-directory setup.
- `scripts/fit_pet_peese_models.R`: optional PET-PEESE fitting and per-meta-analysis RDS generation.
- `scripts/create_analysis_data.R`: optional setup-specific derived-data creation for table/figure rendering.
- `scripts/create_tables_and_figures.R`: creates manuscript tables and figures in numeric order, with independently rerunnable sections.
- `scripts/create_full_tables_figures_and_robustness.R`: optional complete legacy counterfactual, supplementary, and robustness workflow.
- `scripts/create_esr_data.R`: optional recreation of the ESR workbook after the complete legacy workflow.
- `scripts/run_exploratory_regressions.R`: exploratory regression models and diagnostic figures, loading supplied data when regenerated files are absent.
- `scripts/main.R`: short orchestrator that sources the analysis steps.
- `results/main/`: generated main results.
- `results/robustness/`: generated robustness-check results.
- `power and bias.Rproj`: RStudio project file for opening the repository in RStudio.
- `renv.lock`, `renv/activate.R`, `.Rprofile`: `renv` project files for dependency restoration.
