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
if (!exists("show_progress")) show_progress <- TRUE
if (!exists("new_progress_bar")) {
  new_progress_bar <- function(total, label) {
    if (!isTRUE(show_progress)) return(NULL)
    message(label)
    utils::txtProgressBar(min = 0, max = max(1, total), style = 3)
  }
  update_progress_bar <- function(progress_bar, value) {
    if (!is.null(progress_bar)) utils::setTxtProgressBar(progress_bar, value)
  }
  close_progress_bar <- function(progress_bar) {
    if (!is.null(progress_bar)) close(progress_bar)
  }
}

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

random_effect_structure <- function(dat) {
  effect_sizes_by_study <- split(dat$eID, dat$sID, drop = TRUE)
  has_within_study_effect_sizes <- any(vapply(
    effect_sizes_by_study,
    function(effect_size_ids) length(unique(effect_size_ids)) > 1,
    logical(1)
  ))

  if (has_within_study_effect_sizes) {
    ~ 1 | sID/eID
  } else {
    ~ 1 | sID
  }
}

extract_variance_components <- function(model) {
  variance_components <- setNames(as.numeric(model$sigma2), model$s.names)
  list(
    within_study_effect_size = if ("sID/eID" %in% names(variance_components)) {
      as.numeric(variance_components[["sID/eID"]])
    } else {
      0
    },
    between_study = as.numeric(variance_components[["sID"]])
  )
}

extract_total_isq <- function(model) {
  isq_statistics <- i2_ml(model, method = "ratio")
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
      tau2 = NA_real_, within_study_tau2 = NA_real_, isq = NA_real_,
      fallback = TRUE
    ))
  }

  estimate <- dat$yi[[1]]
  standard_error <- sqrt(dat$vi[[1]])
  p_value <- if (is.finite(standard_error) && standard_error > 0) {
    2 * pnorm(-abs(estimate / standard_error))
  } else if (is.finite(estimate) && estimate == 0) 1 else NA_real_

  list(
    estimate = estimate, standard_error = standard_error,
    p_value = p_value, tau2 = 0, within_study_tau2 = 0, isq = 0,
    fallback = TRUE
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
    random = random_effect_structure(dat),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  )
  pet_test <- coef_test(pet, vcov = vcovCR(pet, type = "CR2"))
  pet_intercept_p_value <- extract_coefficient_statistic(
    pet_test, "intrcpt", "p_Satt"
  )

  ## CR2/Satterthwaite inference can return NA when a meta-analysis has too few
  ## independent clusters. In that case the PET-to-PEESE switch cannot be
  ## justified, so retain the already fitted PET model instead of evaluating an
  ## NA in `if` (which aborts the entire parallel foreach job).
  use_pet <- !is.finite(pet_intercept_p_value) ||
    pet_intercept_p_value > 0.10

  if (use_pet) {
    selected_model <- pet
    selected_test <- pet_test
    method <- "PET"
    slope_term <- "sei"
    outcome_scale <- 1
  } else {
    ## Put the complete PEESE model on a numerically more stable scale, rather
    ## than scaling vi alone (which would change the inverse-variance weights
    ## and the fitted variance components). If y* = c*y, then its sampling
    ## variance is c^2*vi. Fitting y* ~ c^2*vi is the same PEESE model in new
    ## units; the intercept/SE and variance components are transformed back
    ## below. The slope test and I-squared are invariant to this change of units.
    outcome_scale <- 10
    peese_data <- dat %>% mutate(
      yi = outcome_scale * yi,
      vi = outcome_scale^2 * vi
    )
    selected_model <- rma.mv(
      yi, vi, mods = ~ 1 + vi,
      random = random_effect_structure(peese_data),
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
    estimate = extract_coefficient_statistic(
      selected_test, "intrcpt", "beta"
    ) / outcome_scale,
    standard_error = extract_coefficient_statistic(
      selected_test, "intrcpt", "SE"
    ) / outcome_scale,
    p_value = extract_coefficient_statistic(selected_test, "intrcpt", "p_Satt"),
    small_study_effect_p_value = extract_coefficient_statistic(
      selected_test, slope_term, "p_Satt"
    ),
    tau2 = extract_variance_components(selected_model)$between_study /
      outcome_scale^2,
    within_study_tau2 =
      extract_variance_components(selected_model)$within_study_effect_size /
      outcome_scale^2,
    isq = extract_total_isq(selected_model)
  )
}

fit_random_effects <- function(dat, capture_model_warnings = FALSE) {
  if (nrow(dat) <= 1) {
    result <- fit_sparse_effects(dat)
    result$model_warnings <- character()
    return(result)
  }

  model_warnings <- character()
  fit_model <- function() {
    rma.mv(
      yi, vi, mods = ~ 1,
      random = random_effect_structure(dat),
      method = "REML", test = "t", data = dat,
      control = list(rel.tol = 1e-8)
    )
  }
  model <- if (capture_model_warnings) {
    withCallingHandlers(
      fit_model(),
      warning = function(warning_condition) {
        model_warnings <<- c(model_warnings, conditionMessage(warning_condition))
        invokeRestart("muffleWarning")
      }
    )
  } else {
    fit_model()
  }
  model_test <- coef_test(model, vcov = vcovCR(model, type = "CR2"))
  variance_components <- extract_variance_components(model)
  list(
    estimate = extract_coefficient_statistic(model_test, "intrcpt", "beta"),
    standard_error = extract_coefficient_statistic(model_test, "intrcpt", "SE"),
    p_value = extract_coefficient_statistic(model_test, "intrcpt", "p_Satt"),
    tau2 = variance_components$between_study,
    within_study_tau2 = variance_components$within_study_effect_size,
    isq = extract_total_isq(model),
    fallback = FALSE,
    model = model,
    model_warnings = model_warnings
  )
}

fit_random_effects_with_outlier_removal <- function(dat, cutoff = 7, n_rounds = 2) {
  analysis_data <- dat

  ## Screen and refit twice. Each screen is based on the random-effects model
  ## fitted to the observations retained by the preceding screen.
  for (round in seq_len(n_rounds)) {
    if (nrow(analysis_data) <= 1) break
    screening_fit <- fit_random_effects(analysis_data)
    standardized_residuals <- as.data.frame(
      rstandard.rma.mv(screening_fit$model)
    )$resid
    keep <- is.na(standardized_residuals) |
      abs(standardized_residuals) <= cutoff
    analysis_data <- analysis_data[keep, , drop = FALSE]
  }

  ## Capture warnings from the final refit so the affected cIDs can be reported
  ## after the parallel workers have returned.
  final_fit <- fit_random_effects(
    analysis_data,
    capture_model_warnings = TRUE
  )
  final_fit$model <- NULL

  list(data = analysis_data, result = final_fit)
}

fit_one_meta_analysis <- function(dat) {
  ## Identify outliers from the initial PET fit, as in the original workflow.
  outlier_model <- rma.mv(
    yi, vi, mods = ~ 1 + sei,
    random = random_effect_structure(dat),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  )
  standardized_residuals <- as.data.frame(rstandard.rma.mv(outlier_model))$resid
  outlier_removed_data <- dat[abs(standardized_residuals) < 3, , drop = FALSE]

  random_effect_outlier_removed <- fit_random_effects_with_outlier_removal(dat)

  analyses <- list(all_data = dat, outlier_removed = outlier_removed_data)
  results <- lapply(analyses, function(analysis_data) {
    list(
      pet_peese = fit_pet_peese(analysis_data),
      pet_peese_data = analysis_data
    )
  })
  results$all_data$random_effect <- fit_random_effects(dat)
  results$all_data$random_effect$model <- NULL
  results$all_data$random_effect_data <- dat
  results$outlier_removed$random_effect <- random_effect_outlier_removed$result
  results$outlier_removed$random_effect_data <- random_effect_outlier_removed$data

  for (variant in names(results)) {
    result <- results[[variant]]
    saveRDS(
      make_effect_data(
        result$pet_peese_data, result$pet_peese,
        result$pet_peese$small_study_effect_p_value
      ),
      file.path(
        output_dirs[[paste0("pet_peese_", variant)]],
        paste0("meta_", dat$cID[[1]], ".rds")
      )
    )
    saveRDS(
      make_effect_data(result$random_effect_data, result$random_effect),
      file.path(
        output_dirs[[paste0("random_effect_", variant)]],
        paste0("meta_", dat$cID[[1]], ".rds")
      )
    )
  }

  tibble(
    cID = as.character(dat$cID[[1]]),
    k_all_data = nrow(dat),
    ## Retain these two legacy columns for the PET-based outlier sample.
    k_outlier_removed = nrow(outlier_removed_data),
    n_outliers_removed = nrow(dat) - nrow(outlier_removed_data),
    k_random_effect_outlier_removed = nrow(random_effect_outlier_removed$data),
    n_random_effect_outliers_removed =
      nrow(dat) - nrow(random_effect_outlier_removed$data),
    pet_peese_method_all_data = results$all_data$pet_peese$method,
    pet_peese_method_outlier_removed = results$outlier_removed$pet_peese$method,
    pet_peese_estimate_all_data = results$all_data$pet_peese$estimate,
    pet_peese_estimate_outlier_removed = results$outlier_removed$pet_peese$estimate,
    random_effect_estimate_all_data = results$all_data$random_effect$estimate,
    random_effect_estimate_outlier_removed = results$outlier_removed$random_effect$estimate,
    random_effect_between_study_variance_all_data =
      results$all_data$random_effect$tau2,
    random_effect_within_study_effect_size_variance_all_data =
      results$all_data$random_effect$within_study_tau2,
    random_effect_between_study_variance_outlier_removed =
      results$outlier_removed$random_effect$tau2,
    random_effect_within_study_effect_size_variance_outlier_removed =
      results$outlier_removed$random_effect$within_study_tau2,
    random_effect_final_model_warnings = paste(
      results$outlier_removed$random_effect$model_warnings,
      collapse = " | "
    )
  )
}

cluster <- makeCluster(n_cores)
registerDoParallel(cluster)
model_progress <- new_progress_bar(length(meta_analyses), "Meta-analysis model progress")
meta_analysis_estimates <- vector("list", ceiling(length(meta_analyses) / n_cores))
analysis_batches <- split(meta_analyses, ceiling(seq_along(meta_analyses) / n_cores))
completed_analyses <- 0
for (batch_index in seq_along(analysis_batches)) {
  meta_analysis_estimates[[batch_index]] <- foreach(
    dat = analysis_batches[[batch_index]],
    .packages = c("metafor", "clubSandwich", "dplyr", "tibble", "orchaRd"),
    .combine = bind_rows
  ) %dopar% fit_one_meta_analysis(dat)
  completed_analyses <- completed_analyses + length(analysis_batches[[batch_index]])
  update_progress_bar(model_progress, completed_analyses)
}
meta_analysis_estimates <- bind_rows(meta_analysis_estimates)
stopCluster(cluster)
close_progress_bar(model_progress)

variance_ratio_warning <- grepl(
  "Ratio of largest to smallest sampling variance extremely large",
  meta_analysis_estimates$random_effect_final_model_warnings,
  fixed = TRUE
)
if (any(variance_ratio_warning)) {
  message(
    "Final random-effects model sampling-variance warning for cID(s): ",
    paste(meta_analysis_estimates$cID[variance_ratio_warning], collapse = ", ")
  )
}

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

summarize_heterogeneity_components <- function(variant, between_study, within_study) {
  eligible <- is.finite(between_study) & is.finite(within_study)
  between_study <- between_study[eligible]
  within_study <- within_study[eligible]
  total <- between_study + within_study
  positive_within <- within_study > 0
  positive_total <- total > 0

  tibble(
    sample = variant,
    n_meta_analyses = length(between_study),
    median_between_study_variance = median(between_study),
    median_within_study_effect_size_variance = median(within_study),
    ratio_of_median_variances = median(between_study) / median(within_study),
    median_within_meta_analysis_variance_ratio =
      median(between_study[positive_within] / within_study[positive_within]),
    median_between_study_percentage_of_total_heterogeneity =
      100 * median(between_study[positive_total] / total[positive_total]),
    percentage_with_zero_within_study_effect_size_variance =
      100 * mean(within_study == 0)
  )
}

random_effect_heterogeneity_comparison <- bind_rows(
  summarize_heterogeneity_components(
    "All data",
    meta_analysis_estimates$random_effect_between_study_variance_all_data,
    meta_analysis_estimates$random_effect_within_study_effect_size_variance_all_data
  ),
  summarize_heterogeneity_components(
    "Outliers removed",
    meta_analysis_estimates$random_effect_between_study_variance_outlier_removed,
    meta_analysis_estimates$random_effect_within_study_effect_size_variance_outlier_removed
  )
)

write_csv(
  random_effect_heterogeneity_comparison,
  file.path(derived_data_dir, "random_effect_heterogeneity_comparison.csv")
)
