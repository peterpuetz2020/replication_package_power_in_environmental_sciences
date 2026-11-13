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

For each observation, let `d_yi`, `d_sei` and `d_vi` be the separately extracted
reported decimal counts. Eligibility is

```r
(d_yi == d_sei) | (d_yi == d_vi)
```

Each equality uses known decimal counts. Compare `yi` against **either column
separately**, rather than comparing it against the minimum of `vi` and `sei`.
For an eligible observation, the matching count is necessarily `d_yi`; set
`mu.dec = sigma.dec = d_yi` for the KP calculation. When both inputs match, they
supply the same count. A match in one column remains eligible if the other
column's precision is missing.

For example, `d_yi = 4`, `d_sei = 4`, `d_vi = 2` is eligible through `sei` and
uses four decimals. Likewise, `d_yi = 4`, `d_sei = 3`, `d_vi = 4` is eligible
through `vi` and uses four decimals. The earlier minimum-based check incorrectly
rejected both of these examples.

`RoundingMatters::min.max.z()` receives `mu.dec` and `sigma.dec` separately.
For unmatched rows, its diagnostic bounds use `mu.dec = d_yi` and
`sigma.dec = min(d_vi, d_sei)` over the available counts; that fallback does not
make the row eligible for the omission analysis.

The coefficient center remains stored `yi`, and the SE center remains `sqrt(vi)`.
Matching the reported decimal count of `vi` is the requested precision convention
**on the SE scale**, not a variance-to-SE precision transformation. Neither
center is rounded for the observed point statistics.

The omission criterion is the **unchanged published rule**:

```r
s <- RoundingMatters::significand(sqrt(vi), num.deci = sigma.dec)
keep <- s >= 37
```

The package computes the significand as `round(sqrt(vi) * 10^sigma.dec)`.
The threshold remains exactly 37, independently of caliper cutoff or half-width.
There is no hypothetical-z criterion and no extension of the rule to unequal
coefficient/SE precision.

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

`RoundingMatters::min.max.z()` supplies diagnostic bounds for every valid pair
with known precision, including unequal precision. Its separate `mu.dec` and
`sigma.dec` arguments do not justify altering the paper's omission criterion.
Those bounds do not determine `kp_keep`. The audit now records the significand,
match flags (`kp_matches_sei_precision`, `kp_matches_vi_precision`), match source
(`kp_precision_match`: `sei`, `vi`, `both`, `none` or `unknown`), eligibility and
status; the custom `kp_reference_z_lower` column has been removed.

The rule follows Kranz and Pütz
([doi:10.1257/aer.20210121](https://doi.org/10.1257/aer.20210121)) and the authors'
package implementation. It reduces rounding problems but does not guarantee
that every retained observation is correctly classified. Its interpretation
requires reporting precision not to be systematically related to the underlying
z distribution. Observed `yi/sqrt(vi)` values are unchanged; no random derounding
or simulated replacements are used.

## Separate rounding-ambiguity sensitivity

The third Figure S6 row starts independently from the **full valid sample**.
It does not first apply the KP omissions and does not require equal coefficient/
SE precision. For every observation, including matches and nonmatches, use

```r
mu.dec <- reported_yi_decimals
sigma.dec <- min(reported_vi_decimals, reported_sei_decimals)
```

As above, `sigma.dec` is applied on the SE scale centered at `sqrt(vi)`. The
coefficient is centered at `abs(yi)`. `RoundingMatters::min.max.z()` gives the
possible ratio bounds, which are converted to absolute z bounds. Let `L` and
`U` denote those lower and upper absolute bounds. Remove the observation when

```r
L < qnorm(0.975) & U >= qnorm(0.975)
```

Both endpoints of the rounding intervals are considered. A value exactly at the
cutoff belongs to the above-cutoff side, consistently with the calipers. An
interval entirely below the cutoff, or entirely at/above it, is retained. No
additional restriction on possible caliper-window membership is imposed.

The `rounding_filter_*` audit columns record this calculation separately from
KP: `sigma_decimals`, `precision_known`, `abs_z_lower`, `abs_z_upper`,
`straddles_cutoff`, `keep` and `status`. In contrast, `rounding_abs_z_*` uses the
KP matching/fallback precision convention described above. For the new filter,
`status` is `omitted_straddles_cutoff`, `retained_one_side`,
`retained_unknown_precision` or `invalid_observation`. Missing precision in one
of vi/sei uses the other available count. Without usable precision the
observation remains unchanged and is flagged; it is not claimed to have a
certain classification.

This is a separate sensitivity analysis of rounding ambiguity, rather than the
published KP omission rule. Removing observations based on whether their bounds
cross the significance cutoff can itself change the local distribution and
caliper shares; the filtered result is not a corrected estimate of publication
bias.

## Caliper estimates and intervals

The cutoff is `qnorm(0.975)`, consistent with the main workflow's two-sided
5% normal threshold. Half-widths are **0.05, 0.10, 0.15 and 0.20** in z units.
Both window endpoints are included; exact cutoff ties count as above.
The point estimate is the observed above-cutoff count divided by the observed
window count. Each estimate has equal weight, so this is a test-level share,
not an equal-weighted mean of article shares or a mean of bootstrap estimates.

Only pooled calipers are calculated, separately for the full, KP-adjusted and rounding-ambiguity
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

- `caliper_tests.xlsx`: twelve pooled rows (three samples × four half-widths) of
  counts, observed shares, bootstrap bounds, `n_retained` (total sample size
  before window restriction) and `n_omitted` (removed from the valid sample),
  with a `Notes` sheet describing the calculation;
- `caliper_coverage.xlsx`: input, valid/excluded and article counts, plus known/
  unknown reporting precision, matching/no-match eligibility, counts of
  sei-only, vi-only and both matches, counts retained after KP or without
  applying KP, KP omission counts and rounding-ambiguity omission/retention
  counts;
- `caliper_rounding_risks.xlsx`: possible significance switches and uncertain
  window membership under the package bounds, counted over each entire sample
  rather than just its observed window members; `precision_basis` identifies
  the bounds convention used by each method. These are not confirmed errors;
- `caliper_analysis.rds`: the results, settings and row audit, including original
  entries, decimal counts, precision sources, ratio bounds, SE significands,
  match source, `kp_rule_applied`, KP inclusion/status, and the separate
  `rounding_filter_*` bounds and inclusion/status.

Figure S6 is saved as PDF and SVG under `results/supplement/`. It has a two-column,
three-row layout: full sample, KP-adjusted sample, then the independent
rounding-ambiguity sample. The KP sample retains nonmatching observations
unchanged. The left panels show histograms within the largest caliper, on shared
x/y scales; the
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
36/37 omission boundary, matches through sei or vi independently, unknown
precision in the nonmatching input, retention of nonmatches and unknown
precision, independent minimum-precision ambiguity filtering, cutoff ties,
option preservation, derived-data handoff and the six-panel figure.
