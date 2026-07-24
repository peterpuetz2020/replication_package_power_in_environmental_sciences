# Replication package: Power in environmental sciences

This repository contains the data and R scripts needed to reproduce the subfield classification, tables, figures, and robustness checks for the paper.

## Computational environment

The project uses [`renv`](https://rstudio.github.io/renv/) to restore the package versions recorded in `renv.lock`. The lockfile records **R 4.5.0** and the CRAN repository snapshot configured through Posit Package Manager. Most packages are installed from CRAN; the non-CRAN dependency `orchaRd` is pinned in `renv.lock` to the GitHub repository `daniel1noble/orchaRd` at commit `5e9ac55cd28d717681bcbcf17527bace42500903`, so `renv::restore()` can install it reproducibly.

> Note: the historical comments in the scripts mention R 3.5.2, but this repository is configured for restoration and replication with the R version recorded in `renv.lock`.

## Restoring packages with renv

1. Install R 4.5.0.
2. Open R in the repository root.
3. Restore the project library:

```r
install.packages("renv")
renv::restore()
```

If GitHub rate limits prevent installation of the pinned `orchaRd` dependency, set a GitHub personal access token before restoring:

```r
Sys.setenv(GITHUB_PAT = "<your-token>")
renv::restore()
```

## Reproducing the results

Run the scripts from the repository root in this order:

```r
source("scripts/classification.R")
source("scripts/main.R")
```

Alternatively, from a shell with R available:

```sh
Rscript scripts/classification.R
Rscript scripts/main.R
```

The classification script classifies meta-analyses into environmental-science subfields and writes the derived classification files to `data/`. The main script uses `scripts/functions.R` to reproduce the main-paper and supplementary tables and figures.

## Output folders

All generated results are saved under `results/`:

- `results/main/` contains the main paper tables, figures, and intermediate result files.
- `results/robustness/` contains robustness-check tables, figures, diagnostics, and intermediate result files.

The main script creates these output folders automatically if they do not already exist.

## Repository structure

- `data/`: source data and derived classification inputs.
- `scripts/classification.R`: classifies meta-analyses into subfields.
- `scripts/functions.R`: helper functions used by the analysis.
- `scripts/main.R`: reproduces main and supplementary results.
- `results/main/`: generated main results.
- `results/robustness/`: generated robustness-check results.
- `renv.lock`, `renv/activate.R`, `.Rprofile`: `renv` project files for dependency restoration.
