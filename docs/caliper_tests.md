# Empirical calipers and Figure S6

The standard workflow computes caliper diagnostics in
`scripts/create_analysis_data.R`, using shared functions in `scripts/functions.R`.
`scripts/create_supplement_tables_and_figures.R` renders **Figure S6**. Running
`scripts/main.R` runs both stages. There is no separate caliper analysis script.

The full empirical sample uses every valid observed row of `data/MasterData.xlsx`,
with `abs(yi/sqrt(vi))` as its normal-approximation statistic. Invalid/nonpositive
variances and nonfinite estimates are excluded and recorded in the row audit.
A second sample applies a Kranz–Pütz omission adjustment with the authors’
[`RoundingMatters`](https://github.com/skranz/RoundingMatters) package, version
0.1.0, pinned in `renv.lock` to commit
`ca9b1e2a086d7107167092bc112e26eb968de4cc`.

## Reporting precision and the omission sample

Read the decimal counts of `yi`, `vi` and `sei` separately from the workbook’s
original entries and explicit number formats. This preserves trailing zeros in
text cells and fixed decimal formats; digits are not inferred from computed
`sqrt(vi)` values or numeric-to-text conversions in R. For General numeric cells,
the stored XML spelling supplies the available precision; trailing zeros already
discarded by Excel cannot be recovered. Formula cells with General formatting
have unknown precision.

For each observation, pass separate precisions to `RoundingMatters::min.max.z()`:

- `mu.dec = decimals(yi)`;
- `sigma.dec = min(decimals(vi), decimals(sei))`.

The center of the coefficient interval is the stored `yi`, and the center of the
SE interval is `sqrt(vi)`, consistent with the rest of the analysis. The requested
minimum decimal count is used **on the SE scale**; this replaces the earlier
custom variance-to-SE precision propagation. Neither center is rounded to its
selected precision. Consequently, the adjusted sample can differ from earlier
runs of that custom implementation.

The omission criterion is the **unchanged published rule**:

```r
s <- RoundingMatters::significand(sqrt(vi), num.deci = sigma.dec)
keep <- s >= 37
```

The package computes the significand as `round(sqrt(vi) * 10^sigma.dec)`.
The threshold remains exactly 37, independently of caliper cutoff or half-width.
There is no hypothetical-z criterion and no extension of the rule to unequal
coefficient/SE precision.

The rule's common-precision assumption is checked separately: it is applied only
where `mu.dec == sigma.dec` under the selected precision convention above.
Unequal-precision observations remain in the full empirical sample and receive
rounding bounds, but have `kp_status = "unequal_precision_not_applicable"` and do
not enter the omission panels. They are counted separately from observations
omitted by the published rule. Missing precision likewise makes the omission
analysis unavailable for that row. The plots therefore compare the full sample
with the **equal-precision eligible subset after the original omission rule**;
differences can reflect both eligibility and omission.

`RoundingMatters::min.max.z()` supplies diagnostic bounds for every valid pair
with known precision, including unequal precision. Its separate `mu.dec` and
`sigma.dec` arguments do not justify altering the paper's omission criterion.
Those bounds do not determine `kp_keep`. The audit now records the significand,
common-precision eligibility and status; the custom `kp_reference_z_lower`
column has been removed.

The rule follows Kranz and Pütz
([doi:10.1257/aer.20210121](https://doi.org/10.1257/aer.20210121)) and the authors'
package implementation. It reduces rounding problems but does not guarantee
that every retained observation is correctly classified. Its interpretation
requires reporting precision not to be systematically related to the underlying
z distribution. Observed `yi/sqrt(vi)` values are unchanged; no random derounding
or simulated replacements are used.

## Caliper estimates and intervals

The cutoff is `qnorm(0.975)`, consistent with the main workflow's two-sided
5% normal threshold. Half-widths are **0.05, 0.10, 0.15 and 0.20** in z units.
Both window endpoints are included; exact cutoff ties count as above.
The point estimate is the observed above-cutoff count divided by the observed
window count. Each estimate has equal weight, so this is a test-level share,
not an equal-weighted mean of article shares or a mean of bootstrap estimates.

Only pooled calipers are calculated, separately for the full and omission
samples, without p-values. Confidence intervals are
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

- `caliper_tests.xlsx`: eight pooled rows (two samples × four half-widths) of
  counts, observed shares, bootstrap bounds and retention counts,
  with a `Notes` sheet describing the calculation;
- `caliper_coverage.xlsx`: input, valid/excluded and article counts, plus known/
  unknown reporting precision, equal/unequal precision eligibility and omission
  counts;
- `caliper_rounding_risks.xlsx`: possible significance switches and uncertain
  window membership under the package bounds, counted over each entire sample
  rather than just its observed window members; these are not confirmed errors;
- `caliper_analysis.rds`: the results, settings and row audit, including original
  entries, decimal counts, precision sources, ratio bounds, SE significands,
  common-precision eligibility and omission status.

Figure S6 is saved as PDF and SVG under `results/supplement/`. It has a two-column,
two-row layout: full sample above, Kranz–Pütz equal-precision omission sample below. The left
panels show histograms within the largest caliper, on shared x/y scales; the
right panels show observed caliper shares and their bootstrap intervals. Each
x-axis tick gives `(n = ...)` for that window. The caliper y range is 0.45–0.65;
complete interval bounds are retained in the workbook when this zoom clips a
tail. The calculation remains integrated into the main workflow, with numeric
workbooks in derived data rather than a separate results directory.

Shared helpers also accept directly supplied z values, or t values with actual
degrees of freedom (`statistic_col`, `df_col`); t values are transformed to
equivalent normal z values. No degrees of freedom are inferred from study counts.

```sh
Rscript tests/test_calipers.R
```

Tests cover window counts, article dependence under within-article duplication,
bootstrap percentiles, unchanged numeric values, t calibration, invalid rows,
Excel precision extraction, unequal-precision ratio bounds, unchanged original
36/37 omission boundary, equal-precision eligibility, unknown precision, option
preservation, derived-data handoff and the four-panel figure.
