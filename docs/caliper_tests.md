# Empirical calipers and Figure S6

The standard workflow computes caliper diagnostics in
`scripts/create_analysis_data.R`, using shared functions in `scripts/functions.R`.
`scripts/create_supplement_tables_and_figures.R` renders **Figure S6**. Running
`scripts/main.R` runs both stages. There is no separate caliper analysis script.

The analysis uses every valid observed row of `data/MasterData.xlsx`, with
`abs(yi/sqrt(vi))` as its normal-approximation statistic. No reporting-precision
selection or omission adjustment is applied. Invalid/nonpositive variances and
nonfinite estimates are excluded and recorded in the row audit.

The cutoff is `qnorm(0.975)`, consistent with the main workflow's two-sided
5% normal threshold. Half-widths are **0.05, 0.10, 0.15 and 0.20** in z units.
Both window endpoints are included; exact cutoff ties count as above.
The point estimate is the observed above-cutoff count divided by the observed
window count. Each estimate has equal weight, so this is a test-level share,
not an equal-weighted mean of article shares or a mean of bootstrap estimates.

Only pooled calipers are calculated, without p-values. Confidence intervals are
**pointwise 95% article-cluster bootstrap percentile intervals**, using **1,000
resamples** and seed **20261107**. Whole `metaID` articles are resampled, including
those with no observations in a window. Each draw recomputes its own numerator
and denominator, with the same article draws used across window widths. Empty
bootstrap windows are counted; bounds require at least 95% valid draws. The
caller's random-number state is preserved.

`metaID` identifies meta-analysis articles. This matches the counterfactual
bootstrap's clustering unit, but does not identify primary studies shared across
different meta-analysis articles. Finite-window probability balance is stronger
than density continuity; a caliper imbalance is not proof of publication bias.

Derived outputs are stored under `data/derived_data/`:

- `caliper_tests.xlsx`: four pooled rows of counts, shares and bootstrap bounds,
  with a `Notes` sheet describing the calculation;
- `caliper_coverage.xlsx`: input, valid/excluded and article counts;
- `caliper_analysis.rds`: the results, row audit and settings used by the renderer.

Figure S6 is saved as PDF and SVG under `results/supplement/`. Its left panel
shows the full-sample histogram within the largest caliper; its right panel shows
caliper shares and their bootstrap intervals. Each x-axis tick gives `(n = ...)`
for that window. The y range is 0.40–0.60; complete interval bounds are retained
in the workbook when this zoom clips a tail. The earlier standalone diagnostic
outputs and precision-filtering workbooks have been removed.

Shared helpers also accept directly supplied z values, or t values with actual
degrees of freedom (`statistic_col`, `df_col`); t values are transformed to
equivalent normal z values. No degrees of freedom are inferred from study counts.

```sh
Rscript tests/test_calipers.R
```

Tests cover window counts, article dependence under within-article duplication,
bootstrap percentiles, unchanged numeric precision, t calibration, invalid rows,
derived-data handoff and the two-panel figure inputs.
