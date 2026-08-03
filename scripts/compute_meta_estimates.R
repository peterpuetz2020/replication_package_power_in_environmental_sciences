## Fit PET-PEESE and multilevel random-effects models with and without outliers.

library(metafor)
library(clubSandwich)
library(tidyverse)
library(readxl)
library(foreach)
library(doParallel)
library(orchaRd)
library(here)

if (!exists("n_cores")) n_cores <- 4

meta <- read_excel(here("data", "MasterData.xlsx"))
meta_analyses <- split(meta, meta$cID)

output_dirs <- c(
  pet_peese_outlier_removed = here("results", "main", "pet_peese_rstandard"),
  pet_peese_all_data = here("results", "main", "pet_peese_all_data"),
  random_effect_outlier_removed = here("results", "main", "multilevel_random"),
  random_effect_all_data = here("results", "main", "multilevel_random_all_data")
)
derived_data_dir <- here("results", "main", "derived_data")
invisible(lapply(
  c(output_dirs, derived_data_dir), dir.create,
  recursive = TRUE, showWarnings = FALSE
))

extract_coefficient_statistic <- function(test, term, statistic) {
  as.numeric(test[term, statistic])
}

extract_between_study_variance <- function(model) {
  variance_components <- setNames(as.numeric(model$sigma2), model$s.names)
  as.numeric(variance_components[["sID"]])
}

extract_total_isq <- function(model) {
  isq_statistics <- i2_ml(model, method = "matrix")
  ## i2_ml() returns total I-squared as its first value. It is not named
  ## "I2_total", so indexing it by that name causes "subscript out of bounds".
  as.numeric(isq_statistics[1])
}

make_effect_data <- function(dat, result, small_study_effect_p_value = NA_real_) {
  dat %>%
    transmute(
      metaID, cID, sID, eID, yi, vi,
      GE = result$estimate, tau2 = result$tau2, isq = result$isq,
      sig_overall = result$p_value,
      small_study_effect_pval = small_study_effect_p_value,
      etype = estype, guide, prere = prereg, subfd = subfield, sdesn = sdesign
    )
}

fit_sparse_effects <- function(dat) {
  if (nrow(dat) == 0) {
    return(list(
      estimate = NA_real_, standard_error = NA_real_, p_value = NA_real_,
      tau2 = NA_real_, isq = NA_real_, fallback = TRUE
    ))
  }

  estimate <- dat$yi[[1]]
  standard_error <- sqrt(dat$vi[[1]])
  p_value <- if (standard_error > 0) {
    2 * pnorm(-abs(estimate / standard_error))
  } else if (estimate == 0) 1 else 0

  list(
    estimate = estimate, standard_error = standard_error,
    p_value = p_value, tau2 = 0, isq = 0, fallback = TRUE
  )
}

fit_pet_peese <- function(dat) {
  if (nrow(dat) <= 1) {
    result <- fit_sparse_effects(dat)
    result$method <- "PET-PEESE not estimable"
    result$small_study_effect_p_value <- NA_real_
    return(result)
  }

  pet <- rma.mv(
    yi, vi, mods = ~ 1 + sei,
    random = list(~ 1 | eID, ~ 1 | sID),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  )
  pet_test <- coef_test(pet, vcov = vcovCR(pet, type = "CR2"))
  pet_intercept_p_value <- extract_coefficient_statistic(
    pet_test, "intrcpt", "p_Satt"
  )

  if (pet_intercept_p_value > 0.10) {
    selected_model <- pet
    selected_test <- pet_test
    method <- "PET"
    slope_term <- "sei"
  } else {
    ## Multiplying vi by 100 for this PEESE fit avoids convergence issues. This
    ## constant rescaling does not alter the extracted PEESE results.
    peese_data <- dat %>% mutate(vi = 100 * vi)
    selected_model <- rma.mv(
      yi, vi, mods = ~ 1 + vi,
      random = list(~ 1 | eID, ~ 1 | sID),
      method = "REML", test = "t", data = peese_data,
      control = list(rel.tol = 1e-8)
    )
    selected_test <- coef_test(
      selected_model, vcov = vcovCR(selected_model, type = "CR2")
    )
    method <- "PEESE"
    slope_term <- "vi"
  }

  list(
    method = method,
    estimate = extract_coefficient_statistic(selected_test, "intrcpt", "beta"),
    standard_error = extract_coefficient_statistic(selected_test, "intrcpt", "SE"),
    p_value = extract_coefficient_statistic(selected_test, "intrcpt", "p_Satt"),
    small_study_effect_p_value = extract_coefficient_statistic(
      selected_test, slope_term, "p_Satt"
    ),
    tau2 = extract_between_study_variance(selected_model),
    isq = extract_total_isq(selected_model)
  )
}

fit_random_effects <- function(dat) {
  if (nrow(dat) <= 1) {
    return(fit_sparse_effects(dat))
  }

  model <- rma.mv(
    yi, vi, mods = ~ 1,
    random = list(~ 1 | eID, ~ 1 | sID),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  )
  model_test <- coef_test(model, vcov = vcovCR(model, type = "CR2"))
  list(
    estimate = extract_coefficient_statistic(model_test, "intrcpt", "beta"),
    standard_error = extract_coefficient_statistic(model_test, "intrcpt", "SE"),
    p_value = extract_coefficient_statistic(model_test, "intrcpt", "p_Satt"),
    tau2 = extract_between_study_variance(model),
    isq = extract_total_isq(model),
    fallback = FALSE
  )
}

fit_one_meta_analysis <- function(dat) {
  ## Identify outliers from the initial PET fit, as in the original workflow.
  outlier_model <- rma.mv(
    yi, vi, mods = ~ 1 + sei,
    random = list(~ 1 | eID, ~ 1 | sID),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  )
  standardized_residuals <- as.data.frame(rstandard.rma.mv(outlier_model))$resid
  outlier_removed_data <- dat[abs(standardized_residuals) < 3, , drop = FALSE]

  analyses <- list(all_data = dat, outlier_removed = outlier_removed_data)
  results <- lapply(analyses, function(analysis_data) {
    list(
      data = analysis_data,
      pet_peese = fit_pet_peese(analysis_data),
      random_effect = fit_random_effects(analysis_data)
    )
  })

  for (variant in names(results)) {
    result <- results[[variant]]
    saveRDS(
      make_effect_data(
        result$data, result$pet_peese,
        result$pet_peese$small_study_effect_p_value
      ),
      file.path(
        output_dirs[[paste0("pet_peese_", variant)]],
        paste0("meta_", dat$cID[[1]], ".rds")
      )
    )
    saveRDS(
      make_effect_data(result$data, result$random_effect),
      file.path(
        output_dirs[[paste0("random_effect_", variant)]],
        paste0("meta_", dat$cID[[1]], ".rds")
      )
    )
  }

  tibble(
    cID = as.character(dat$cID[[1]]),
    k_all_data = nrow(dat),
    k_outlier_removed = nrow(outlier_removed_data),
    n_outliers_removed = nrow(dat) - nrow(outlier_removed_data),
    pet_peese_method_all_data = results$all_data$pet_peese$method,
    pet_peese_method_outlier_removed = results$outlier_removed$pet_peese$method,
    pet_peese_estimate_all_data = results$all_data$pet_peese$estimate,
    pet_peese_estimate_outlier_removed = results$outlier_removed$pet_peese$estimate,
    random_effect_estimate_all_data = results$all_data$random_effect$estimate,
    random_effect_estimate_outlier_removed = results$outlier_removed$random_effect$estimate
  )
}

cluster <- makeCluster(n_cores)
registerDoParallel(cluster)
meta_analysis_estimates <- foreach(
  dat = meta_analyses,
  .packages = c("metafor", "clubSandwich", "dplyr", "tibble", "orchaRd"),
  .combine = bind_rows
) %dopar% fit_one_meta_analysis(dat)
stopCluster(cluster)

meta_analysis_estimates <- meta_analysis_estimates %>% arrange(cID)
write_csv(
  meta_analysis_estimates,
  file.path(derived_data_dir, "meta_analysis_estimates.csv")
)

summarize_outlier_difference <- function(estimator, all_data, outlier_removed) {
  eligible <- is.finite(all_data) & is.finite(outlier_removed)
  all_data <- all_data[eligible]
  outlier_removed <- outlier_removed[eligible]
  sign_change <- all_data != 0 & outlier_removed != 0 &
    sign(all_data) != sign(outlier_removed)

  tibble(
    estimator = estimator,
    n_meta_analyses = length(all_data),
    median_absolute_difference = median(abs(outlier_removed - all_data)),
    mean_absolute_difference = mean(abs(outlier_removed - all_data)),
    percentage_sign_changes = 100 * mean(sign_change)
  )
}

pet_peese_outlier_comparison <- summarize_outlier_difference(
  "PET-PEESE",
  meta_analysis_estimates$pet_peese_estimate_all_data,
  meta_analysis_estimates$pet_peese_estimate_outlier_removed
)
random_effect_outlier_comparison <- summarize_outlier_difference(
  "Multilevel random effects",
  meta_analysis_estimates$random_effect_estimate_all_data,
  meta_analysis_estimates$random_effect_estimate_outlier_removed
)

write_csv(
  pet_peese_outlier_comparison,
  file.path(derived_data_dir, "pet_peese_outlier_comparison.csv")
)
write_csv(
  random_effect_outlier_comparison,
  file.path(derived_data_dir, "random_effect_outlier_comparison.csv")
)
