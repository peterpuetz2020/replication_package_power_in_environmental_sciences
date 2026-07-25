# Teshome K. Deressa, Peter Pütz, David I. Stern, Jaco Vangronsveld, Jan Minx, Sebastien Lizin, Robert Malina, Stephan B. Bruns

## Overview
This is the replication package of “Low statistical power and overrepresentation of statistically
significant findings in the environmental sciences” (including the Supporting Information).

The repository contains the data and R scripts needed to reproduce the subfield classification, tables, figures, and robustness checks.

## Computational environment

The project uses [`renv`](https://rstudio.github.io/renv/) to restore the package versions recorded in `renv.lock`. The lockfile records **R 4.5.0** and the CRAN repository snapshot configured through Posit Package Manager. Most packages are installed from CRAN; the non-CRAN dependency `orchaRd` is pinned in `renv.lock` to the GitHub repository `daniel1noble/orchaRd` at commit `5e9ac55cd28d717681bcbcf17527bace42500903`, so `renv::restore()` can install it reproducibly.

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

This script creates outputs in numeric order: Table 1, Figure 1, Table 2, Table 3, and Figure 2. Each table/figure section reloads the data, settings, and grid definitions it needs, so a single section can be run independently in a fresh R session after the helper setup at the beginning of the file has been sourced. Table 2 is written once for every combination of `meta_average_multipliers` and `heterogeneity_multipliers`; with multiple combinations, filenames use labels such as `meta_0p5_heterogeneity_0p25`. A supplied `analysis_setups` data frame can instead assign a custom, unique `setup_label` to each combination.

### Optional setup-specific data creation

The data-preparation step is separated from table/figure rendering. If you change the setup parameters at the beginning of `scripts/create_analysis_data.R` (for example `meta_average_multipliers`, `heterogeneity_multipliers`, `n_iterations`, or `setup_label`), run:

```r
source("scripts/create_analysis_data.R")
```

This writes setup-specific derived datasets under `results/main/derived_data/`. `meta_average_multipliers` and `heterogeneity_multipliers` can each contain one or more values; `scripts/create_analysis_data.R` computes and stores outputs for every combination. To use custom labels, define an `analysis_setups` tibble/data frame with `meta_average_multiplier`, `heterogeneity_multiplier`, and `setup_label` columns before sourcing the script. By default, it creates the derived analysis datasets only and leaves the expensive counterfactual simulations untouched. Set `recreate_counterfactuals <- TRUE` in `scripts/create_analysis_data.R` only when you also want to rebuild the setup-specific counterfactual z-value and p-value RDS files for each selected combination. `scripts/create_tables_and_figures.R` uses these derived datasets when available; otherwise, it falls back to the existing PET-PEESE RDS files under `results/main/pet_peese_rstandard/`.

### Optional data-recreation step

The main analysis in `scripts/main.R` reads `data/MasterData.xlsx`. Other files in `data/` support the optional classification workflow or are derived/intermediate files, but they are not read by `scripts/main.R`. Running `scripts/classification.R` is only necessary if you want to recreate the classification files from the raw Scopus and Scimago inputs. If you only want to reproduce the tables, figures, and robustness checks from the available main-analysis data, you can skip `scripts/classification.R` and run only:

```r
source("scripts/main.R")
```

## Runtime settings in `scripts/main.R`

At the beginning of `scripts/main.R`, adjust:

- `n_cores`: number of parallel worker cores.
- `n_iterations`: number of Monte Carlo/bootstrap iterations for confidence intervals.

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
- `scripts/create_analysis_data.R`: optional setup-specific derived-data creation for table/figure rendering.
- `scripts/create_tables_and_figures.R`: creates manuscript tables and figures in numeric order, with independently rerunnable sections.
- `scripts/main.R`: reproduces main and supplementary results.
- `results/main/`: generated main results.
- `results/robustness/`: generated robustness-check results.
- `power and bias.Rproj`: RStudio project file for opening the repository in RStudio.
- `renv.lock`, `renv/activate.R`, `.Rprofile`: `renv` project files for dependency restoration.
