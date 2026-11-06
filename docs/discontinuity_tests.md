# Empirical discontinuity diagnostics

Run separately from the manuscript workflow:

```sh
Rscript scripts/discontinuity_tests.R
```

This uses **all observed rows** of `data/MasterData.xlsx`, with `yi` as the
coefficient and `sei` as its standard error. It does not use counterfactual
distributions, fitted meta-means, or the model-based outlier screen.

## Tests and estimands

The null value is zero. The default statistic is `abs(yi / sei)`, with cutoff
`qnorm(0.975)` for a two-sided normal-approximation p-value of 0.05. Caliper
half-widths are 0.05, 0.075, 0.1, 0.2, 0.3, 0.4 and 0.5 in **z units**.
Window endpoints are included, and exact cutoff ties count as above, matching
the author replication's binomial/randomization procedure. `cutoff = 1.96`
can be supplied to reproduce the paper's rounded normal cutoff literally.

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
bias. The joint coefficient/SE plot shows the relevant boundary
`abs(coefficient) = cutoff * SE`. There is no single corresponding cutoff in
the marginal coefficient or SE distribution by itself.

## Kranz–Pütz rounding adjustment

The primary omission rule from
[Kranz and Pütz (2022), *American Economic Review* 112(9), 3124–3136](https://www.aeaweb.org/articles?id=10.1257/aer.20210121)
retains coefficient/SE observations with **SE significand at least 37**. For
example, an SE of `0.04` reported to two decimal places has significand 4,
whereas `0.0400` reported to four decimal places has significand 400.

Following the data owner's instruction, the script treats the workbook's
`yi` and `sei` entries as the reported values. It reads their original XML
spellings, preserves string-cell trailing zeros, and recognizes simple fixed
decimal and scientific number formats. Explicit display rounding is reflected
in the values being tested. Formula cells with General format and unsupported
formats are flagged as having unknown reporting precision. General numeric
cells cannot recover trailing zeros already lost when Excel stored a number;
their precision source is therefore recorded as `excel_stored_numeric_entry`.
The workbook itself is never modified.

The audit records decimal places, SE significands, keep/drop decisions and
precision provenance. It also records whether the coefficient and SE share a
decimal resolution, as assumed by the paper's Lemma 1. Rounding intervals use
the individual coefficient/SE resolutions and flag possible changes in
significance classification or caliper membership. The significand-37 rule
does not guarantee every possible rounding problem is eliminated.
Its interpretation also assumes that reporting coarseness is not systematically
related to the unrounded z distribution; selection on precision can otherwise
change the share of significant results.

For inputs lacking recoverable reporting precision, an optional
`kp_inferred_digits_sensitivity` explicitly assumes the numeric representation
at 15 significant digits, with trailing fractional zeros removed, gives the
precision. It is labelled **unverified sensitivity**,
not pooled with the reported-precision correction, and is not generated when
all valid rows have known precision. Disable it with `infer_precision = FALSE`.
Directly reported z/t statistics are retained by the omission rule, as in the
paper: its coefficient/SE significand rule does not apply to direct reports.

The implementation follows the omission criterion and binomial/randomization
diagnostic, not the paper's design-comparison probit specifications or its
separate latent-density analyses. Author sources:

- [Published-comment replication, pinned commit](https://github.com/peterpuetz2020/replication_repo_methods_matter_comment/tree/7220ff41dcf6b9c7648171c6d864c7871669841e): `Scripts/make_new_data.R` and `Scripts/deround_bch.R`.
- [RoundingMatters reference implementation](https://github.com/skranz/RoundingMatters/tree/ca9b1e2a086d7107167092bc112e26eb968de4cc): `R/window_funs.R`, `R/risks.R`, and `R/digit_tools.R`.

## Other inputs and original strings

Source the file to access functions without automatically running an analysis:

```r
source("scripts/discontinuity_tests.R")

# A data frame with original coefficient/SE strings, preserving zeros.
d$coefficient_decimals <- reported_decimal_places(d$coefficient_as_reported)
d$se_decimals <- reported_decimal_places(d$se_as_reported)
run_discontinuity_tests(
  d, estimate_col = "coefficient_as_reported", se_col = "se_as_reported",
  estimate_digits_col = "coefficient_decimals", se_digits_col = "se_decimals",
  cluster_col = "article_id", subgroup_col = NULL, infer_precision = FALSE
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

Explicit precision columns override automatic workbook extraction. Use
`se_col = NULL, variance_col = "vi"` to derive SEs from variances if necessary;
this cannot establish the original SE reporting precision. `use_workbook_entries
= FALSE` disables automatic workbook-entry extraction.

## Outputs and validation

`results/discontinuity_tests/` contains:

- `caliper_tests.csv`: counts, nominal and cluster-robust inference, Holm
  adjustments, retention counts and status flags;
- `precision_coverage.csv`: valid/excluded rows and precision coverage;
- `rounding_risks.csv`: classification/window-membership risk counts;
- `empirical_discontinuity_diagnostics.pdf`: empirical local distribution,
  caliper shares, the joint coefficient/SE distribution, and the adjusted
  distribution;
- `METHODS.txt`: conventions and limitations for interpreting the run.

The row-level audit is saved to the ignored
`data/derived_data/discontinuity_tests/observation_audit.rds`. Alternative output
and audit directories can be supplied to the runner.

```sh
Rscript tests/test_discontinuity_tests.R
```

Tests cover hand-calculated binomial and clustered inference, duplicated
observations within clusters, unavailable/degenerate cases, SE significands
and rounding bounds, preserved Excel trailing zeros and display precision,
Student-t calibration, and a planted discontinuity.

The default workbook run includes 64,810 valid observations and retains 61,839
under the omission rule. Only 7,286 rows have the same reported coefficient/SE
decimal resolution. The one-sided binomial results were also checked against
the pinned author implementation for every caliper width on both full and
omission-adjusted samples at the paper's literal cutoff of 1.96; SE significands
and rounding bounds were checked against the author functions on the workbook.
