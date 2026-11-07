# WAAP-WLS and Figure S5

WAAP-WLS is an inverse-variance weighted mean, not a multilevel model. It is
included alongside PET–PEESE, random effects and the existing zero-heterogeneity
fixed-effect comparison in Figure S5. Main-text estimator choices are unchanged.

For every `cID` and both shared samples (all data and the existing RE outlier
screen), the estimator first calculates the full-sample inverse-variance mean.
Using its magnitude as the assumed true effect, it calculates two-sided normal
power at alpha 0.05 with the same exact formula as the main power workflow.
Estimates with power >= 0.80 qualify. With at least two qualifying estimates,
the mean is fitted to that subset; otherwise full-sample WLS is used. Signed
means are retained. This selection precedes every meta-average multiplier.

The selected fit uses exactly the fixed-effect implementation:
`metafor::rma.mv(..., method = "FE")`, study-level `sID` clusters,
`clubSandwich::vcovCR(type = "CR2")` and `coef_test()` Satterthwaite inference.
There are no random effects in this model. The selected subset can contain
fewer primary studies than the full sample; with fewer than two study clusters,
the WAAP point estimate is retained and its inference is flagged unavailable.
It is not silently replaced with independent-observation standard errors.

The adequate subset determines the fitted mean only. All effects in the shared
variant sample remain available for the empirical ESR reference and denominators.
WAAP has no estimator of `tau2`; Figure S5's nonzero-heterogeneity scenarios
therefore combine the WAAP mean with the corresponding existing random-effects
heterogeneity estimate on that same variant sample. This combination is labelled
in the figure and audit workbook; it is not a multilevel WAAP model.

Run after the existing meta-analysis inputs are available:

```sh
Rscript scripts/compute_waap_wls.R
```

`scripts/main.R` runs this stage automatically. The existing fitting stage also
creates WAAP inputs after producing the common samples. Individual results are
cached under the ignored `data/derived_data/waap_wls/` and
`waap_wls_all_data/` directories; each cache includes its source-file hash and
method settings. `results/supplement/WAAP_WLS_estimates.xlsx` records means,
CR2 standard errors, degrees of freedom, confidence intervals, adequate/selected
counts, independent study counts, inference status and WLS fallback use.
Rebuilding the standard table/figure workflow adds WAAP to Figure S5 for all
four meta-average and all four genuine-heterogeneity shares.
The updated Figure S5 uses 1,000 article-bootstrap replications, matching the
existing figure's caches. Set `n_iterations <- 1000` before the full manuscript
workflow for that resolution.

```sh
Rscript tests/test_waap_wls.R
Rscript tests/test_waap_outputs.R
```

Tests cover adequate-subset selection, WLS fallback, sign symmetry, unavailable
single-cluster inference, and numerical equivalence of standard errors, p-values
and degrees of freedom to the existing fixed-effect CR2/Satterthwaite calculation.
