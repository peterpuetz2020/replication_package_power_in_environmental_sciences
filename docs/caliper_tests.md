# Empirical calipers and Figure S5

The standard workflow computes caliper diagnostics in
`scripts/create_analysis_data.R`, using shared functions in `scripts/functions.R`.
`scripts/create_supplement_tables_and_figures.R` renders **Figure S5**. Running
`scripts/main.R` runs both stages.

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

For each observation, let `d_yi`, `d_sei` and `d_vi` be the separately extracted
reported decimal counts. Eligibility is

```r
(d_yi == d_sei) | (d_yi == d_vi)
```

Each equality uses known decimal counts. Compare `yi` against **either column
separately**, rather than comparing it against the minimum of `vi` and `sei`.
For an eligible observation, the matching count is necessarily `d_yi`; set
`sigma.dec = d_yi` for the KP calculation. When both inputs match, they
supply the same count. A match in one column remains eligible if the other
column's precision is missing.

For example, `d_yi = 4`, `d_sei = 4`, `d_vi = 2` is eligible through `sei` and
uses four decimals. Likewise, `d_yi = 4`, `d_sei = 3`, `d_vi = 4` is eligible
through `vi` and uses four decimals. Both matches are eligible even if the
other input is coarser.

The SE value remains `sqrt(vi)`, consistent with the rest of the analysis.
Matching the reported decimal count of `vi` is the requested precision convention
**on the SE scale**, not a variance-to-SE precision transformation. Observed
point statistics are not rounded. `kp_sigma_decimals` and `kp_se_significand`
are calculated only when the KP rule applies; they are `NA` for nonmatching or
unknown-precision observations, which are retained unchanged.

The omission criterion is the **unchanged published rule**:

```r
s <- RoundingMatters::significand(sqrt(vi), num.deci = sigma.dec)
keep <- s >= 37
```

The package computes the significand as `round(sqrt(vi) * 10^sigma.dec)`.
The threshold remains exactly 37, independently of caliper cutoff or half-width.

The rule is applied only where `d_yi == d_sei` **or** `d_yi == d_vi`, using
the matching count. A matching observation is omitted only when it fails the
unchanged threshold. **All valid observations without a precision match remain
in the adjusted sample unchanged.** Unknown reporting precision likewise means
that the rule is not applied and the observation is retained. Invalid observed
estimates or variances remain excluded from both samples.

The audit distinguishes rule application from final sample inclusion:

- `kp_rule_applied`: a valid observation has a known match through sei or vi;
- `kp_keep`: the observation is included in the adjusted sample, either because
  it passes KP or because KP was not applied;
- `kp_status = "retained"`: matched and passed KP;
- `kp_status = "omitted_by_precision_rule"`: matched and failed KP;
- `kp_status = "retained_no_precision_match"`: no match; retained unchanged;
- `kp_status = "retained_unknown_precision"`: precision unavailable; retained
  unchanged;
- `kp_status = "invalid_observation"`: invalid point data; excluded.

Thus the adjusted sample equals the full valid sample **minus matched KP
failures**. It is not restricted to precision-matched observations.

The audit records the significand, match flags (`kp_matches_sei_precision`,
`kp_matches_vi_precision`), match source (`kp_precision_match`: `sei`, `vi`,
`both`, `none` or `unknown`), rule application and inclusion status. No interval
bounds or rounding-risk summaries are calculated: neither is needed for the
current KP rule or Figure S5.

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

Only pooled calipers are calculated, separately for the full and KP-adjusted samples, without p-values. Confidence intervals are
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

Numerical results are saved only in `data/derived_data/caliper_analysis.rds`:

- `tests`: eight pooled rows (two samples × four half-widths), with counts,
  observed shares, bootstrap bounds, total `n_retained` and `n_omitted`;
- `coverage`: input and valid/excluded counts, reporting-precision matches,
  counts retained after KP or without applying KP, and omission counts;
- `audit`: observed coefficients, variances and z values, original reported
  entries, decimal counts, precision sources, applicable SE significands,
  match source, `kp_rule_applied`, and KP inclusion/status;
- `settings`: cutoff, half-widths, bootstrap settings, data source and KP rule.

No caliper XLSX workbooks are written.

Figure S5 is saved as PDF and SVG under `results/supplement/`. Its two columns
are **Full empirical sample** and **Kranz–Pütz adjustment**, with one title per
column. Histograms are above the corresponding caliper estimates. Both
histograms use the same bins and x/y scales within the largest caliper.
The caliper panels have no titles or subtitles, and there is no bottom caption
or legend. The interval method and resample count are recorded above and in
the RDS settings. Each x-axis tick gives `(n = ...)` for that window. The caliper y
range is 0.45–0.65; complete interval bounds are retained in the RDS when
this zoom clips a tail. Only the full sample and KP-adjusted sample are calculated.

Figures S2–S4 show counterfactual sensitivities at multipliers 0.25, 0.75 and 1;
Figure S6 compares the meta-analytic estimators.

```sh
Rscript tests/test_calipers.R
```

Tests cover window counts, article dependence under within-article duplication,
bootstrap percentiles, unchanged numeric values, invalid rows, Excel precision
extraction, the original 36/37 omission boundary, matches through sei or vi independently, unknown
precision in the nonmatching input, retention of nonmatches and unknown
precision, option preservation, derived-data handoff and the two sample columns
in Figure S5.
