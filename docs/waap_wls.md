# WAAP-WLS and Figure S5

WAAP-WLS is an inverse-variance weighted mean, not a multilevel model. It is
included alongside PET–PEESE, random effects and fixed effects in Figure S5.
Both WAAP-WLS and fixed effects appear **only at zero genuine heterogeneity**.
Main-text estimator choices are unchanged.

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
WAAP has no estimator of `tau2`. Its counterfactuals therefore assume zero genuine
heterogeneity, just like the fixed-effect comparison. No random-effects variance
components are borrowed. Stored `tau2` and `isq` are zero by assumption, with
`heterogeneity_source = "not_estimated_zero_assumed"`; they are not WAAP estimates
of these quantities. Figure generation rejects nonzero-heterogeneity scenarios
for either estimator, including stale figure inputs.

WAAP-WLS is implemented in the common fitting stage alongside the other estimators:

```sh
Rscript scripts/compute_meta_estimates.R
```

`scripts/main.R` runs the common fitting stage automatically. When the other
estimator inputs already exist, it uses the same script in its cached-input
mode to update WAAP without refitting the other models. There is no separate
WAAP calculation script. Individual results are
cached under the ignored `data/derived_data/waap_wls/` and
`waap_wls_all_data/` directories; each cache includes its source-file hash and
method settings. The ignored `data/derived_data/WAAP_WLS_estimates.xlsx` records means,
CR2 standard errors, degrees of freedom, confidence intervals, adequate/selected
counts, independent study counts, inference status and WLS fallback use.
WAAP estimates are saved only under derived data; Figure S5 is the manuscript
output displaying them. An obsolete copy of the workbook in `results/supplement/`
is removed after the derived-data workbook has been written successfully.
Rebuilding the standard table/figure workflow adds WAAP to Figure S5 for all
four meta-average multipliers at zero genuine heterogeneity only.
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
