# Empirical discontinuity diagnostics

Run separately from the manuscript workflow:

```sh
Rscript scripts/discontinuity_tests.R
```

This uses **all observed rows** of `data/MasterData.xlsx`, with `yi` as the
coefficient and **`vi` as the reported variance**. By default, the variance
rounding resolution uses the smaller of the reported decimal counts in `vi`
and `sei`, separately for each observation. The SE point value is computed
as `sqrt(vi)`.
It does not use counterfactual
distributions, fitted meta-means, or the model-based outlier screen.

## Tests and estimands

The null value is zero. The default statistic is `abs(yi / sqrt(vi))`, with cutoff
`qnorm(0.975)`, the exact normal threshold for a two-sided p-value of 0.05,
consistent with the main workflow. Caliper
half-widths are 0.05, 0.1, 0.15 and 0.2 in **z units**. The runner rejects
half-widths larger than 0.2.
Window endpoints are included, and exact cutoff ties count as above, matching
the author replication's binomial/randomization procedure. Supply
`cutoff = 1.96` only when the paper's rounded threshold is specifically needed.

For each caliper, the script reports the fraction above the cutoff, one-sided
and two-sided exact binomial tests against 0.5, and binomial confidence
intervals. Those nominal tests assume independent observations. The additional
cluster-robust results use meta-analysis article `metaID` by default, a
`G/(G-1)` sandwich correction for the estimated share, and a `t(G-1)` reference.
Missing clusters, fewer than two clusters, and degenerate cluster variance are
explicitly flagged instead of producing spurious precision. Primary studies
may overlap across meta-analysis articles; this default clustering does not
identify those overlaps.

Pooled and subfield-specific tests are returned. Each p-value column also has
a Holm adjustment across all available windows, groups and methods. Local
probability balance in a finite window is stronger than density continuity;
these are **discontinuity diagnostics**, not proof of p-hacking or publication
bias. The relevant coefficient/SE boundary is
`abs(coefficient) = cutoff * SE`; there is no single corresponding cutoff in
the marginal coefficient or SE distribution by itself.

## Kranz–Pütz rule and its variance-based extension

The primary omission rule from
[Kranz and Pütz (2022), *American Economic Review* 112(9), 3124–3136](https://www.aeaweb.org/articles?id=10.1257/aer.20210121)
retains coefficient/SE observations with **SE significand at least 37**. For
example, an SE of `0.04` reported to two decimal places has significand 4,
whereas `0.0400` reported to four decimal places has significand 400.

The **default extends that rule to reported variances**, as requested by the
data owner. A decimal place in a variance is not a decimal place in an SE, so
the code does not copy the selected decimal count to `sqrt(vi)` or round the
SE to that many decimals. The user's requested conservative precision rule is:

```text
d = min(reported decimals of vi, reported decimals of sei)
Delta_v = 10^(-d)
v_L = max(0, v - Delta_v/2)       v_U = v + Delta_v/2
SE_L = sqrt(v_L)                SE_U = sqrt(v_U)
SE = sqrt(v)
Delta_SE_equiv = 2 * max(SE - SE_L, SE_U - SE)
s_equiv = SE / Delta_SE_equiv
```

The stored variance and coefficient values are not rounded to this selected
resolution. When `sei` has fewer decimals, the selected resolution describes
a wider uncertainty interval around the stored variance, not its original
rounding grid. Both original decimal counts and the selected minimum are
recorded in the audit. Both must be known for reported-precision mode; missing
metadata can only enter the separately labelled inferred-digit sensitivity.

The variance-based omission retains `s_equiv >= 37`. `Delta_SE_equiv` is a
conservative symmetric SE rounding step enclosing the exact transformed
variance interval. This preserves the relative-SE-precision interpretation of
the original criterion, but **is an extension, not the paper's literal
reported-SE-significand rule**. The score is not an integer reported SE
significand. If `k = v / Delta_v` is at least 0.5, then
`s_equiv = k * (1 + sqrt(1 - 1/(2*k)))`, approximately `2*k` at high
precision. With the minimum-decimal rule, `k` need not be an integer; the
audit records it as `variance_resolution_ratio`. The reported variance's own
integer significand is recorded separately.

For example, variance `0.04` reported to two decimals implies an SE interval
`[sqrt(0.035), sqrt(0.045)]`, giving `s_equiv` about 7.74 and omission.
Variance `0.0400` reported to four decimals gives a score about 799.75 and
retention, provided `sei` has at least the corresponding number of decimals.
The computed SE is 0.2 in both cases.

The default reads the workbook's `yi`, `vi` and `sei` entries in their original XML
spellings, preserves string-cell trailing zeros, and recognizes simple fixed
decimal and scientific number formats. Explicit display rounding is reflected
in the values being tested. Formula cells with General format and unsupported
formats are flagged as having unknown reporting precision. General numeric
cells cannot recover trailing zeros already lost when Excel stored a number;
their precision source is therefore recorded as `excel_stored_numeric_entry`.
The workbook itself is never modified.

The audit records variance decimal places and significands, equivalent SE
precision scores, keep/drop decisions, exact variance/SE bounds, original
reported entries and precision provenance. `se_decimals` and `se_significand`
are unavailable in variance mode, because an SE is computed rather than
directly reported. Coefficient bounds use the separate `yi` resolution:

```text
z_L = max(0, abs(yi) - 10^(-d_y)/2) / SE_U
z_U = (abs(yi) + 10^(-d_y)/2) / SE_L
```

An SE lower bound of zero gives an infinite upper z bound. These intervals
flag possible changes in significance classification or caliper membership.
The paper's Lemma 1 assumes coefficients and SEs reported at the same decimal
resolution; that assumption does not describe the variance-based extension.
The precision-score rule
does not guarantee every possible rounding problem is eliminated.
Its interpretation also assumes that reporting coarseness is not systematically
related to the unrounded z distribution; selection on precision can otherwise
change the share of significant results.

For inputs lacking recoverable reporting precision, an optional
`kp_inferred_variance_digits_sensitivity` explicitly assumes the numeric representation
at 15 significant digits, with trailing fractional zeros removed, gives the
precision. It is labelled **unverified sensitivity**,
not pooled with the reported-precision correction, and is not generated when
all valid rows have known precision. Disable it with `infer_precision = FALSE`.
The default corrected results are labelled `kp_reported_min_variance_se_precision`.
Use `precision_basis = "variance"` for the earlier variance-only precision rule
and its label `kp_reported_variance_precision`.
For directly reported SE inputs, `precision_basis = "se"` restores the original
rule and the labels `kp_reported_precision` /
`kp_inferred_digits_sensitivity`.
Directly reported z/t statistics are retained by the omission rule, as in the
paper: its coefficient/SE significand rule does not apply to direct reports.

The implementation follows the binomial/randomization diagnostic and either
the original omission criterion or the explicitly labelled variance extension,
not the paper's design-comparison probit specifications or its
separate latent-density analyses. Author sources:

- [Published-comment replication, pinned commit](https://github.com/peterpuetz2020/replication_repo_methods_matter_comment/tree/7220ff41dcf6b9c7648171c6d864c7871669841e): `Scripts/make_new_data.R` and `Scripts/deround_bch.R`.
- [RoundingMatters reference implementation](https://github.com/skranz/RoundingMatters/tree/ca9b1e2a086d7107167092bc112e26eb968de4cc): `R/window_funs.R`, `R/risks.R`, and `R/digit_tools.R`.

## Other inputs and original strings

Source the file to access functions without automatically running an analysis:

```r
source("scripts/discontinuity_tests.R")

# A data frame with original coefficient/variance/SE strings, preserving zeros.
d$coefficient_decimals <- reported_decimal_places(d$coefficient_as_reported)
d$variance_decimals <- reported_decimal_places(d$variance_as_reported)
d$se_decimals <- reported_decimal_places(d$se_as_reported)
run_discontinuity_tests(
  d, estimate_col = "coefficient_as_reported", variance_col = "variance_as_reported",
  se_col = "se_as_reported", se_digits_col = "se_decimals",
  estimate_digits_col = "coefficient_decimals", variance_digits_col = "variance_decimals",
  cluster_col = "article_id", subgroup_col = NULL, infer_precision = FALSE
)

# Use the literal published criterion for directly reported SEs.
run_discontinuity_tests(
  "data/MasterData.xlsx", precision_basis = "se",
  output_dir = "results/discontinuity_tests_reported_se",
  audit_dir = "data/derived_data/discontinuity_tests_reported_se"
)

# Direct t values require their actual degrees of freedom.
run_discontinuity_tests(
  t_data, statistic_col = "t_value", df_col = "df",
  cluster_col = "article_id", subgroup_col = NULL
)
```

t inputs are converted to equivalent absolute z values using their two-sided
Student-t p-values, so the common cutoff still means p = 0.05. The workbook
contains no reported t values or degrees of freedom; its default ratio must
be described as a **normal approximation**, not as an exact t test. No df is
invented from the number of studies.

Explicit precision columns override automatic workbook extraction. In default
minimum-decimal mode (`precision_basis = "min_variance_se"`), `variance_col` /
`variance_digits_col` and `se_col` / `se_digits_col` supply the two reported
decimal counts. `se_col` is not used to compute the point statistic. In SE mode, `se_col` and
`se_digits_col` control precision. `use_workbook_entries = FALSE` disables
automatic workbook-entry extraction.

## Outputs and validation

`results/discontinuity_tests/` contains:

- `caliper_tests.xlsx`: counts, nominal and cluster-robust inference, Holm
  adjustments, retention counts and status flags;
- `precision_coverage.xlsx`: valid/excluded rows and precision coverage;
- `rounding_risks.xlsx`: classification/window-membership risk counts;
- `empirical_discontinuity_diagnostics.pdf`: a two-column comparison of full
  and omission-filtered empirical samples, as described below;
- `empirical_discontinuity_diagnostics.svg`: the same ggplot figure in SVG format;
- `METHODS.txt`: conventions and limitations for interpreting the run.

Each XLSX file has a `Results` sheet containing numeric results and a `Notes`
sheet documenting the input, precision rule, confidence intervals and Holm
correction. Obsolete CSV versions of these three generated outputs are removed
when their XLSX replacements are written. The default run returns **64 test
rows**: 8 groups (pooled plus 7 subfields), 2 methods and 4 widths. Holm
adjustments are recalculated for these 64 rows separately per p-value column.

The ggplot diagnostic figure has empirical histograms on the left and caliper-share
plots on the right. The upper row uses all valid observations; the lower row
uses the Kranz–Puetz omission sample, meaning the filtered observations from
`MasterData.xlsx`, not the economics dataset in the authors' original paper.
The default lower row uses the variance-based extension. If a missing-precision
sensitivity is needed, that sample is explicitly labelled instead.

Both caliper panels use **y limits 0.35–0.70**. Their error bars are **pointwise
95% article-clustered Student-t intervals**, based on `metaID`, the sandwich
variance correction `G/(G-1)` and `t(G-1)` quantiles. They are not the exact
binomial intervals also saved in the workbook, and are not simultaneous or
Holm-adjusted intervals. Each panel identifies its input sample, and
each caliper's own `n_window` is printed as `(n = ...)` just below its half-width
on the x-axis. Workbook-name and total-N subtitle lines are omitted. The
histograms filter the data to **qnorm(0.975) +/- 0.2**, approximately **[1.76, 2.16]**, by default and
use exactly those x limits without padding. They have shared bins and vertical
scales and show their own local sample sizes. The default tests use the same
exact normal threshold as the main workflow, so the largest caliper and histogram include the same observations.

Clustering is appropriate when multiple estimates share a source of dependence.
In publication-bias applications that source is often the original research
article. Here `metaID` identifies a **meta-analysis article**, consistent with
the main workflow's article-cluster bootstrap. It is a practical cluster level
for this repository, but cannot account for primary studies shared across
different meta-analysis articles. Globally matched primary-article identifiers
would support more direct checks of that dependence. The Kranz–Puetz replication
uses article clustering in its design-comparison probits; its exact-binomial
window figures assume independent observations. Our clustered share intervals
are an additional dependence adjustment, not a literal reproduction of those
binomial figures.

The row-level audit is saved to the ignored
`data/derived_data/discontinuity_tests/observation_audit.rds`. Alternative output
and audit directories can be supplied to the runner.

```sh
Rscript tests/test_discontinuity_tests.R
```

Tests cover hand-calculated binomial and clustered inference, duplicated
observations within clusters, unavailable/degenerate cases, SE significands,
propagated variance bounds and their agreement with an independent rounding
grid, stable precision scores at high precision, preserved Excel trailing zeros
and display precision, selecting the minimum of `vi`/`sei` decimals per row,
keeping point values unchanged, ignoring `sei` point values in variance mode, Student-t calibration,
and a planted discontinuity.

The default workbook run includes 64,810 valid observations. Retention counts
under the minimum-decimal variance extension are saved in the coverage workbook.
The underlying one-sided binomial procedure
was also checked against
the pinned author implementation for every caliper width on both full and
original-SE-omission samples at the paper's literal cutoff of 1.96. The original
SE-mode significands and rounding bounds were checked against the author
functions; the variance extension is validated separately by its mathematical
interval calculations and regression tests.
