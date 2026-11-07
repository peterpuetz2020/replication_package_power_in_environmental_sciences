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
2. fits fixed-effects, PET-PEESE, WAAP-WLS and multilevel random-effects estimates when their derived
   RDS inputs are missing;
3. creates setup-specific analysis data and counterfactual simulations;
4. writes the main and supplementary tables and figures.

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
| `meta_average_multiplier` | `c(0.25, 0.5, 0.75, 1)` | Multipliers applied to the meta-analytic average in power calculations. |
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
  compute_meta_estimates.R                fixed-effects, PET-PEESE, WAAP-WLS and random-effects fitting
  create_analysis_data.R                  setup-specific derived data
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
scalable PDF and SVG files. PDF, DOCX and XLSX manuscript outputs are versioned;
SVG and LaTeX files are ignored by Git and recreated by the scripts.

Setup-specific intermediate filenames contain their multiplier choices, for
example `meta_0p5_heterogeneity_0p25`. Cache keys include the data, grid,
settings, bootstrap clusters, and iteration count, so incompatible cached
counterfactual results are rebuilt automatically.

## Current output numbering

The exploratory negative-binomial/OLS regressions, their input workbooks,
main Table 4, and their supplementary tables and diagnostic figures have been
removed. PET-PEESE and Egger small-study-effect tests remain part of the
meta-analysis methods.

The default analysis evaluates **16 parameter combinations**: meta-average
multipliers 0.25, 0.5, 0.75 and 1, crossed with four heterogeneity multipliers.
Figure S5 includes all four meta-average values in a two-column, two-row layout. Main Tables 1–3 and Figures
1–2 retain the manuscript's half-meta-average specification.

Figure S5 additionally includes WAAP-WLS, with primary-study clustering and the
same CR2/Satterthwaite correction as fixed effects. Both estimators are shown
only at zero genuine heterogeneity; neither borrows random-effects variance estimates. See
[WAAP-WLS estimation details](docs/waap_wls.md). All WAAP estimates and fitting
summaries are saved under `data/derived_data/`; Figure S5 displays the results.

| Supplement output | Contents |
|---|---|
| Table S1 | Excess significance by subfield |
| Tables S2–S4 | Power summaries at meta-average multipliers 0.25, 0.75, 1 |
| Table S5 | Power after excluding non-significant meta-analyses |
| Table S6 | Counterfactual p-value distributions at multipliers 0.5 versus 1 |
| Table S7 | Effect-size measures by subfield |
| Table S8 | Small-study effects by subfield |
| Figure S1 | Heterogeneity distributions |
| Figures S2–S4 | Counterfactual sensitivities at multipliers 0.25, 0.75, 1 |
| Figure S5 | ESR comparisons for all four multipliers, arranged in two columns |

Table S6 has exactly Table 2's layout: 13 p-value-bin rows, the two
`ESR_0.05` measures (among all results and among significant results), and
meta-analysis/test sample sizes, with genuine heterogeneity as columns.
It retains the 0.5 counterfactual and replaces Table 2's empirical reference
with the counterfactual at multiplier 1. Bin differences are therefore
**reference (1) minus counterfactual (0.5)**. If `q_m` is the counterfactual
probability of `p < 0.05` at multiplier `m`, the ESR rows are
`ESR_all = q_1 - q_0.5` and `ESR_sig = (q_1 - q_0.5) / q_1`.

Both distributions use identical valid effects and paired article-cluster
bootstrap draws. Each significant-result ESR draw uses its own reference
significant probability as denominator. The table reports estimates and
95% percentile intervals, matching Table 2's formatting. Full-precision rows,
including both counterfactual values and sample sizes, are saved to
`data/derived_data/counterfactual_comparison_0p5_vs_1.csv`.

Tables are written as DOCX and LaTeX, and figures as PDF and SVG. Use 1,000
bootstrap iterations for final manuscript outputs. A quick numerical check is:

```sh
Rscript tests/test_counterfactual_comparison.R
```

`tests/test_bootstrap_parallel.R` checks that sequential, parallel and batched
execution preserve the same bootstrap draws. After a default analysis run,
`tests/test_analysis_outputs.R` checks the 16-setup output inventory and verifies
Table S6's intervals against the paired counterfactual caches.
This output check defaults to 1,000 iterations; pass an exploratory run's
iteration count as its first argument when checking a smaller run.

## Separate empirical discontinuity tests

Run `Rscript scripts/caliper_tests.R` to test the observed distribution
near two-sided p = 0.05 using calipers and the Kranz–Pütz (2022) rounding-omission
method. The script reads reported `yi`/`vi`/`sei` precision from the workbook entries,
uses the smaller `vi`/`sei` decimal count for each variance rounding interval,
propagates variance rounding through `sqrt(vi)`, and applies an explicitly
labelled variance-based extension of the omission rule. It reports
article-clustered inference and writes separate diagnostics under
`results/discontinuity_tests/`. It is independent of `scripts/main.R`.
Caliper half-widths are 0.05, 0.1, 0.15 and 0.2. Summary tables are XLSX
workbooks; ggplot diagnostics saved as SVG and PDF compare full and omission-filtered samples
using only pointwise 95% article-cluster bootstrap percentile intervals from 1,000 resamples.
Caliper output is pooled only and contains confidence intervals without p-values.
See [the method, precision conventions and input options](docs/discontinuity_tests.md).
