## ---------------------------------------------------------------------------
## compare_meta_analysis_estimators.R
## ---------------------------------------------------------------------------
## Re-estimate every meta-analysis using the extended PET-PEESE procedure used
## in the original analysis and fixed- and random-effects meta-analysis
## estimators. The random-effects estimator uses the same effect- and
## study-level random intercepts as PET-PEESE. The resulting estimates are then
## compared on the same, PET-screened set of observations.

library(metafor)
library(clubSandwich)
library(tidyverse)
library(readxl)
library(foreach)
library(doParallel)
library(here)
n_cores <- 4
## Use a value supplied by the calling session, if present.
n_cores <- if (exists("n_cores")) n_cores else 7

input_file <- here("data", "MasterData.xlsx")
output_dir <- here("results", "estimator_comparison")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
estimator_data_dirs <- c(
  multilevel_random = here("results", "main", "multilevel_random"),
  fixed = here("results", "main", "fixed")
)
invisible(lapply(estimator_data_dirs, dir.create, recursive = TRUE,
                 showWarnings = FALSE))

meta <- read_excel(input_file) %>%
  ## These models also fail in the original PET-PEESE analysis because their
  ## PET/PEESE model matrices are rank deficient.
  filter(!cID %in% c("259-1", "259-2", "2204-2", "850-1"))

meta_analyses <- split(meta, meta$cID)

## Fit PET first, remove observations with an absolute standardized residual of
## at least three, and refit PET. PEESE is selected when the PET intercept is
## significant at 10%, as in the original analysis.
fit_pet_peese <- function(dat) {
  pet_initial <- rma.mv(
    yi, vi,
    mods = ~ 1 + sei,
    random = list(~ 1 | eID, ~ 1 | sID),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  )
  
  standardized_residual <- as.data.frame(rstandard.rma.mv(pet_initial))$resid
  keep <- abs(standardized_residual) < 3
  screened_data <- dat[keep, , drop = FALSE]
  
  ## A very small meta-analysis can lose nearly all observations to the screen.
  ## Refitting PET then fails before the conventional-model fallback is ever
  ## reached (metafor: "Processing terminated since k <= 1"). Retain the
  ## unscreened data when the screened PET model cannot support its intercept,
  ## slope, and cluster-robust inference.
  residual_screen_fallback <- nrow(screened_data) < 3 ||
    dplyr::n_distinct(screened_data$sei) < 2 ||
    dplyr::n_distinct(screened_data$sID) < 2
  analysis_data <- if (residual_screen_fallback) dat else screened_data
  
  pet <- rma.mv(
    yi, vi,
    mods = ~ 1 + sei,
    random = list(~ 1 | eID, ~ 1 | sID),
    method = "REML", test = "t", data = analysis_data,
    control = list(rel.tol = 1e-8)
  )
  pet_test <- coef_test(pet, vcov = vcovCR(pet, type = "CR2"))
  
  if (pet_test$p_Satt[1] > 0.10) {
    selected_model <- pet
    selected_test <- pet_test
    selected_method <- "PET"
  } else {
    selected_model <- rma.mv(
      yi, vi,
      mods = ~ 1 + vi,
      random = list(~ 1 | eID, ~ 1 | sID),
      method = "REML", test = "t", data = analysis_data,
      control = list(rel.tol = 1e-8)
    )
    selected_test <- coef_test(
      selected_model,
      vcov = vcovCR(selected_model, type = "CR2")
    )
    selected_method <- "PEESE"
  }
  
  list(
    data = analysis_data,
    method = selected_method,
    estimate = as.numeric(selected_test$beta[1]),
    standard_error = as.numeric(selected_test$SE[1]),
    p_value = as.numeric(selected_test$p_Satt[1]),
    ## Once PEESE is selected, its variance slope (rather than the PET standard
    ## error slope) is the small-study-effect test reported by the original
    ## implementation.
    small_study_effect_p_value = as.numeric(selected_test$p_Satt[2]),
    n_outliers_identified = sum(!keep),
    n_outliers_removed = if (residual_screen_fallback) 0L else sum(!keep),
    residual_screen_fallback = residual_screen_fallback
  )
}

## Fit the comparison estimators to exactly the observations retained by the
## PET residual screen. The random-effects model reproduces PET-PEESE's two
## random intercepts. A fixed-effects model cannot contain random effects by
## definition, but rma.mv still permits the same multivariate data interface;
## CR2 inference is used for both models, as it is for PET-PEESE.
fit_one_meta_analysis <- function(dat) {
  pet_peese <- fit_pet_peese(dat)
  analysis_data <- pet_peese$data
  
  ## Do this check before calling rma.mv(): metafor terminates immediately for
  ## some one-effect inputs, so inspecting fixed$k after fitting is too late.
  conventional_data <- analysis_data %>%
    filter(is.finite(yi), is.finite(vi), vi >= 0)
  conventional_k <- nrow(conventional_data)
  
  if (conventional_k == 0) {
    stop("No finite yi/vi pairs remain for cID ", analysis_data$cID[1])
  } else if (conventional_k == 1) {
    ## With one effect, both conventional estimators have the closed-form
    ## inverse-variance result. Between-effect heterogeneity is unidentifiable,
    ## so report the boundary value tau^2 = 0 without invoking rma.mv().
    fixed_effect_estimate <- conventional_data$yi[1]
    fixed_effect_se <- sqrt(conventional_data$vi[1])
    fixed_effect_p_value <- if (fixed_effect_se > 0) {
      2 * pnorm(-abs(fixed_effect_estimate / fixed_effect_se))
    } else if (fixed_effect_estimate == 0) {
      1
    } else {
      0
    }
    random_effect_estimate <- fixed_effect_estimate
    random_effect_se <- fixed_effect_se
    random_effect_p_value <- fixed_effect_p_value
    random_effect_tau2 <- 0
    random_effect_isq <- 0
    random_effect_fallback <- TRUE
  } else {
    fixed <- rma.mv(
      yi, vi,
      mods = ~ 1,
      data = conventional_data,
      method = "FE", test = "t",
      control = list(rel.tol = 1e-8)
    )
    fixed_test <- coef_test(
      fixed,
      vcov = vcovCR(fixed, cluster = conventional_data$sID, type = "CR2")
    )

    random <- rma.mv(
      yi, vi,
      ## This is a separate intercept-only meta-analysis: unlike PET and
      ## PEESE, it does not regress yi on sei or vi.
      mods = ~ 1,
      random = list(~ 1 | eID, ~ 1 | sID),
      data = conventional_data,
      method = "REML", test = "t",
      control = list(rel.tol = 1e-8)
    )
    random_test <- coef_test(random, vcov = vcovCR(random, type = "CR2"))

    fixed_effect_estimate <- as.numeric(fixed_test$beta[1])
    fixed_effect_se <- as.numeric(fixed_test$SE[1])
    fixed_effect_p_value <- as.numeric(fixed_test$p_Satt[1])
    random_effect_estimate <- as.numeric(random_test$beta[1])
    random_effect_se <- as.numeric(random_test$SE[1])
    random_effect_p_value <- as.numeric(random_test$p_Satt[1])
    ## sigma2[2] is the between-study variance, matching the PET-PEESE code.
    random_effect_tau2 <- as.numeric(random$sigma2[2])
    random_effect_isq <- as.numeric(orchaRd::i2_ml(random, method = "matrix")[1])
    random_effect_fallback <- FALSE
  }

  ## Preserve the effect-level structure consumed by the downstream power and
  ## counterfactual analyses. Each meta-analytic estimate is repeated for all
  ## observations retained by the common PET residual screen, just as in the
  ## original PET-PEESE RDS files.
  make_estimator_data <- function(estimate, tau2, isq, p_value) {
    conventional_data %>%
      transmute(
        metaID, cID, sID, eID, yi, vi,
        GE = estimate,
        tau2 = tau2,
        isq = isq,
        sig_overall = p_value,
        small_study_effect_pval = NA_real_,
        etype = estype,
        guide, prere = prereg, subfd = subfield, sdesn = sdesign
      )
  }

  saveRDS(
    make_estimator_data(
      random_effect_estimate, random_effect_tau2,
      random_effect_isq, random_effect_p_value
    ),
    file.path(estimator_data_dirs[["multilevel_random"]],
              paste0("meta_", analysis_data$cID[1], ".rds"))
  )
  saveRDS(
    make_estimator_data(fixed_effect_estimate, 0, 0, fixed_effect_p_value),
    file.path(estimator_data_dirs[["fixed"]],
              paste0("meta_", analysis_data$cID[1], ".rds"))
  )
  
  tibble(
    cID = as.character(analysis_data$cID[1]),
    k_original = nrow(dat),
    k_analyzed = nrow(analysis_data),
    k_conventional = conventional_k,
    n_outliers_identified = pet_peese$n_outliers_identified,
    n_outliers_removed = pet_peese$n_outliers_removed,
    residual_screen_fallback = pet_peese$residual_screen_fallback,
    pet_peese_method = pet_peese$method,
    pet_peese_estimate = pet_peese$estimate,
    pet_peese_se = pet_peese$standard_error,
    pet_peese_p_value = pet_peese$p_value,
    small_study_effect_p_value = pet_peese$small_study_effect_p_value,
    fixed_effect_estimate = fixed_effect_estimate,
    fixed_effect_se = fixed_effect_se,
    fixed_effect_p_value = fixed_effect_p_value,
    random_effect_estimate = random_effect_estimate,
    random_effect_se = random_effect_se,
    random_effect_p_value = random_effect_p_value,
    random_effect_tau2 = random_effect_tau2,
    random_effect_fallback = random_effect_fallback
  )
}

run_meta_analyses <- function(meta_analyses, n_cores) {
  cluster <- makeCluster(n_cores)
  on.exit(stopCluster(cluster), add = TRUE)
  registerDoParallel(cluster)
  
  foreach(
    dat = meta_analyses,
    .packages = c("metafor", "clubSandwich", "dplyr", "tibble", "orchaRd"),
    ## The fitting functions are referenced indirectly from this wrapper, so
    ## foreach's automatic global detection does not reliably export them to
    ## PSOCK workers (notably on Windows). Export both functions explicitly.
    .export = c("fit_one_meta_analysis", "fit_pet_peese", "estimator_data_dirs"),
    .combine = bind_rows
  ) %dopar% fit_one_meta_analysis(dat)
}

estimates <- run_meta_analyses(meta_analyses, n_cores)

## Compare effect magnitudes rather than signed estimates. A positive deviation
## means that the comparison method produces a stronger effect (farther from
## zero), irrespective of whether the estimates are positive or negative. The
## percentage uses the absolute original estimate as its denominator; it is
## undefined when the original estimate is exactly zero.
estimates <- estimates %>%
  mutate(
    pet_peese_vs_fixed_deviation =
      abs(pet_peese_estimate) - abs(fixed_effect_estimate),
    pet_peese_vs_fixed_percentage_deviation = if_else(
      fixed_effect_estimate == 0,
      NA_real_,
      100 * pet_peese_vs_fixed_deviation / abs(fixed_effect_estimate)
    ),
    pet_peese_vs_random_deviation =
      abs(pet_peese_estimate) - abs(random_effect_estimate),
    pet_peese_vs_random_percentage_deviation = if_else(
      random_effect_estimate == 0,
      NA_real_,
      100 * pet_peese_vs_random_deviation / abs(random_effect_estimate)
    ),
    fixed_vs_random_deviation =
      abs(fixed_effect_estimate) - abs(random_effect_estimate),
    fixed_vs_random_percentage_deviation = if_else(
      random_effect_estimate == 0,
      NA_real_,
      100 * fixed_vs_random_deviation / abs(random_effect_estimate)
    )
  ) %>%
  arrange(cID)

summarize_comparison <- function(comparison, comparison_estimate,
                                 original_estimate) {
  magnitude_deviation <- abs(comparison_estimate) - abs(original_estimate)
  percentage_deviation <- if_else(
    original_estimate == 0,
    NA_real_,
    100 * magnitude_deviation / abs(original_estimate)
  )

  tibble(
    comparison = comparison,
    median_magnitude_deviation = median(magnitude_deviation, na.rm = TRUE),
    median_percentage_deviation = median(percentage_deviation, na.rm = TRUE),
    mean_absolute_difference = mean(
      abs(comparison_estimate - original_estimate), na.rm = TRUE
    ),
    rmse = sqrt(mean((comparison_estimate - original_estimate)^2, na.rm = TRUE)),
    correlation = cor(
      comparison_estimate, original_estimate,
      method = "kendall", use = "complete.obs"
    )
  )
}

comparison_summary <- bind_rows(
  summarize_comparison(
    "PET-PEESE vs fixed effect (original)",
    estimates$pet_peese_estimate, estimates$fixed_effect_estimate
  ),
  summarize_comparison(
    "PET-PEESE vs random effects (original)",
    estimates$pet_peese_estimate, estimates$random_effect_estimate
  ),
  summarize_comparison(
    "Fixed effect vs random effects (original)",
    estimates$fixed_effect_estimate, estimates$random_effect_estimate
  )
)

## Summarize changes in absolute magnitude in a form that can be reported as
## percentage attenuation. Ratios use the conventional estimator as the
## reference (comparison / original), so 100 * (1 - ratio) is positive when
## the comparison estimate is closer to zero. Exact zero comparison estimates
## are valid ratios and imply a 100% reduction. Sign reversals require two
## non-zero estimates; an estimate of zero has no direction to reverse.
summarize_attenuation <- function(comparison, comparison_estimate,
                                  original_estimate) {
  eligible <- is.finite(comparison_estimate) &
    is.finite(original_estimate) & original_estimate != 0
  comparison_estimate <- comparison_estimate[eligible]
  original_estimate <- original_estimate[eligible]
  magnitude_ratio <- abs(comparison_estimate) / abs(original_estimate)
  percentage_reduction <- 100 * (1 - magnitude_ratio)
  n_comparisons <- length(magnitude_ratio)

  ## exp(mean(log(0))) correctly returns zero if an estimate is attenuated all
  ## the way to zero. This is the continuous extension of the geometric mean.
  geometric_mean_ratio <- if (n_comparisons == 0) {
    NA_real_
  } else {
    exp(mean(log(magnitude_ratio)))
  }

  tibble(
    comparison = comparison,
    n_meta_analyses = n_comparisons,
    geometric_mean_percentage_reduction =
      100 * (1 - geometric_mean_ratio),
    median_percentage_reduction = median(percentage_reduction),
    percentage_attenuated = 100 * mean(magnitude_ratio < 1),
    percentage_amplified = 100 * mean(magnitude_ratio > 1),
    percentage_unchanged = 100 * mean(magnitude_ratio == 1),
    percentage_sign_reversal = 100 * mean(
      comparison_estimate != 0 &
        sign(comparison_estimate) != sign(original_estimate)
    )
  )
}

attenuation_summary <- bind_rows(
  summarize_attenuation(
    "PET-PEESE vs fixed effect (reference)",
    estimates$pet_peese_estimate, estimates$fixed_effect_estimate
  ),
  summarize_attenuation(
    "PET-PEESE vs random effects (reference)",
    estimates$pet_peese_estimate, estimates$random_effect_estimate
  ),
  summarize_attenuation(
    "Fixed effect vs random effects (reference)",
    estimates$fixed_effect_estimate, estimates$random_effect_estimate
  )
)

## Repeat the attenuation comparison after excluding pair-specific sign
## reversals. Zero estimates are retained because they have no direction and
## therefore do not constitute a sign reversal.
summarize_attenuation_without_sign_reversals <- function(
    comparison, comparison_estimate, original_estimate) {
  sign_reversal <- is.finite(comparison_estimate) &
    is.finite(original_estimate) &
    comparison_estimate != 0 &
    original_estimate != 0 &
    sign(comparison_estimate) != sign(original_estimate)
  
  summarize_attenuation(
    comparison,
    comparison_estimate[!sign_reversal],
    original_estimate[!sign_reversal]
  )
}

attenuation_summary_without_sign_reversals <- bind_rows(
  summarize_attenuation_without_sign_reversals(
    "PET-PEESE vs fixed effect (reference)",
    estimates$pet_peese_estimate, estimates$fixed_effect_estimate
  ),
  summarize_attenuation_without_sign_reversals(
    "PET-PEESE vs random effects (reference)",
    estimates$pet_peese_estimate, estimates$random_effect_estimate
  ),
  summarize_attenuation_without_sign_reversals(
    "Fixed effect vs random effects (reference)",
    estimates$fixed_effect_estimate, estimates$random_effect_estimate
  )
)

write_csv(estimates, file.path(output_dir, "meta_analysis_estimates.csv"))
write_csv(comparison_summary, file.path(output_dir, "estimator_comparison_summary.csv"))
write_csv(
  attenuation_summary,
  file.path(output_dir, "estimator_attenuation_summary.csv")
)
write_csv(
  attenuation_summary_without_sign_reversals,
  file.path(
    output_dir,
    "estimator_attenuation_summary_without_sign_reversals.csv"
  )
)

comparison_plot_data <- estimates %>%
  select(cID, pet_peese_estimate,
         `Fixed effect` = fixed_effect_estimate,
         `Random effects` = random_effect_estimate) %>%
  pivot_longer(
    c(`Fixed effect`, `Random effects`),
    names_to = "estimator", values_to = "estimate"
  )

x_limits <- c(-6, 6)
y_limits <- c(-6, 6)

r2_labels <- comparison_plot_data %>%
  filter(
    is.finite(pet_peese_estimate),
    is.finite(estimate),
    between(pet_peese_estimate, -6, 6),
    between(estimate, -6, 6)
  ) %>%
  group_by(estimator) %>%
  summarise(
    r2 = summary(lm(estimate ~ pet_peese_estimate))$r.squared,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("R² = %.3f", r2),
    x = x_limits[1] + 0.04 * diff(x_limits),
    y = y_limits[2] - 0.04 * diff(y_limits)
  )

comparison_plot <- ggplot(
  comparison_plot_data,
  aes(x = pet_peese_estimate, y = estimate)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    color = "grey55",
    linetype = 2
  ) +
  geom_point(alpha = 0.35, size = 1) +
  geom_text(
    data = r2_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1
  ) +
  facet_wrap(~ estimator) +
  coord_equal(
    xlim = x_limits,
    ylim = y_limits,
    expand = FALSE
  ) +
  labs(
    x = "PET-PEESE estimate",
    y = "Comparison estimate"
  ) +
  theme_bw()
comparison_plot

ggsave(
  file.path(output_dir, "estimator_comparison.pdf"),
  comparison_plot, width = 8, height = 4.5
)

print(comparison_summary)
print(attenuation_summary)
# how many sign changes from positive to negative by applying pet-peese
estimates %>% filter(random_effect_estimate > 0 & pet_peese_estimate < 0) %>% nrow()
# how many positive random effects at all 
estimates %>% filter(random_effect_estimate > 0) %>% nrow()
# vice verse
estimates %>% filter(random_effect_estimate < 0 & pet_peese_estimate > 0) %>% nrow()
estimates %>% filter(random_effect_estimate < 0) %>% nrow()

# seems not to be systematic, but remove estimates with sign changes beforehand for robustness check
print(attenuation_summary_without_sign_reversals)

temp <- estimates |> 
  dplyr::select(pet_peese_estimate, random_effect_estimate, pet_peese_vs_random_deviation)

# or the current 704 meta-analyses, PET–PEESE estimates were geometrically 24.05% 
# smaller than random-effects estimates, with a 7.82% median reduction; 57.10%
# were attenuated, 42.90% amplified, and 22.73% reversed direction.

# drop small effects which might affect results
estimates_drop_small <- estimates |> 
  filter(!between(random_effect_estimate, -0.1, 0.1))

attenuation_summary_drop_small <- bind_rows(
  summarize_attenuation(
    "PET-PEESE vs fixed effect (reference)",
    estimates_drop_small$pet_peese_estimate, estimates_drop_small$fixed_effect_estimate
  ),
  summarize_attenuation(
    "PET-PEESE vs random effects (reference)",
    estimates_drop_small$pet_peese_estimate, estimates_drop_small$random_effect_estimate
  ),
  summarize_attenuation(
    "Fixed effect vs random effects (reference)",
    estimates_drop_small$fixed_effect_estimate, estimates_drop_small$random_effect_estimate
  )
)
print(attenuation_summary)
print(attenuation_summary_drop_small)
# results remain quite stable (for median)
