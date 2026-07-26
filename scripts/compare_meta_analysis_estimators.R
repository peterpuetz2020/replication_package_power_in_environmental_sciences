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

## Use a value supplied by the calling session, if present.
n_cores <- if (exists("n_cores")) n_cores else 7

input_file <- here("data", "MasterData.xlsx")
output_dir <- here("results", "estimator_comparison")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
    random_effect_fallback <- TRUE
  } else {
    fixed <- rma.mv(
      yi, vi,
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
    random_effect_fallback <- FALSE
  }
  
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
    .packages = c("metafor", "clubSandwich", "dplyr", "tibble"),
    ## The fitting functions are referenced indirectly from this wrapper, so
    ## foreach's automatic global detection does not reliably export them to
    ## PSOCK workers (notably on Windows). Export both functions explicitly.
    .export = c("fit_one_meta_analysis", "fit_pet_peese"),
    .combine = bind_rows
  ) %dopar% fit_one_meta_analysis(dat)
}

estimates <- run_meta_analyses(meta_analyses, n_cores)

## Paired differences retain the direction of the discrepancy for every
## meta-analysis. Aggregate comparisons summarize its size and association.
estimates <- estimates %>%
  mutate(
    fixed_minus_pet_peese = fixed_effect_estimate - pet_peese_estimate,
    random_minus_pet_peese = random_effect_estimate - pet_peese_estimate,
    fixed_minus_random = fixed_effect_estimate - random_effect_estimate
  ) %>%
  arrange(cID)

comparison_summary <- tribble(
  ~comparison, ~mean_difference, ~median_difference, ~mean_absolute_difference, ~rmse, ~correlation,
  "Fixed effect minus PET-PEESE",
  mean(estimates$fixed_minus_pet_peese),
  median(estimates$fixed_minus_pet_peese),
  mean(abs(estimates$fixed_minus_pet_peese)),
  sqrt(mean(estimates$fixed_minus_pet_peese^2)),
  cor(estimates$fixed_effect_estimate, estimates$pet_peese_estimate, method = "kendall"),
  "Random effects minus PET-PEESE",
  mean(estimates$random_minus_pet_peese),
  median(estimates$random_minus_pet_peese),
  mean(abs(estimates$random_minus_pet_peese)),
  sqrt(mean(estimates$random_minus_pet_peese^2)),
  cor(estimates$random_effect_estimate, estimates$pet_peese_estimate, method = "kendall"),
  "Fixed effect minus random effects",
  mean(estimates$fixed_minus_random),
  median(estimates$fixed_minus_random),
  mean(abs(estimates$fixed_minus_random)),
  sqrt(mean(estimates$fixed_minus_random^2)),
  cor(estimates$fixed_effect_estimate, estimates$random_effect_estimate, method = "kendall")
)

write_csv(estimates, file.path(output_dir, "meta_analysis_estimates.csv"))
write_csv(comparison_summary, file.path(output_dir, "estimator_comparison_summary.csv"))

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
