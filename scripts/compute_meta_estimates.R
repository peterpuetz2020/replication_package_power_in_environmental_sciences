## Fit PET-PEESE and multilevel random-effects models with and without outliers.

library(metafor)
library(clubSandwich)
library(tidyverse)
library(readxl)
library(foreach)
library(doParallel)
library(orchaRd)
library(here)

if (!exists("n_cores")) n_cores <- 6
if (!exists("show_progress")) show_progress <- TRUE
minimum_primary_studies <- 5L
if (!exists("new_progress_bar")) {
  new_progress_bar <- function(total, label) {
    if (!isTRUE(show_progress)) return(NULL)
    message(label)
    progress_bar <- utils::txtProgressBar(
      min = 0, max = max(1, total), initial = 0, style = 3,
      file = stderr()
    )
    flush.console()
    progress_bar
  }
  update_progress_bar <- function(progress_bar, value) {
    if (!is.null(progress_bar)) {
      utils::setTxtProgressBar(progress_bar, value)
      flush.console()
    }
  }
  close_progress_bar <- function(progress_bar) {
    if (!is.null(progress_bar)) {
      close(progress_bar)
      flush.console()
    }
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
derived_data_dir <- here("data", "derived_data")
invisible(lapply(
  c(output_dirs, derived_data_dir), dir.create,
  recursive = TRUE, showWarnings = FALSE
))

extract_coefficient_statistic <- function(test, term, statistic) {
  as.numeric(test[term, statistic])
}

coefficient_test <- function(model, dat) {
  study_cluster <- droplevels(factor(dat$sID))

  ## The workflow checks that at least five primary studies remain before any
  ## final model is fitted, so CR2 always has multiple independent clusters.
  coef_test(
    model,
    vcov = vcovCR(model, cluster = study_cluster, type = "CR2")
  )
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

primary_study_count <- function(dat) {
  dplyr::n_distinct(dat$sID, na.rm = TRUE)
}

drop_meta_analysis <- function(dat, reason, pet_studies = NA_integer_,
                               random_effect_studies = NA_integer_) {
  c_id <- as.character(dat$cID[[1]])

  ## A rerun must not leave output files from an earlier, less restrictive run.
  invisible(lapply(output_dirs, function(output_dir) {
    unlink(file.path(output_dir, paste0("meta_", c_id, ".rds")))
  }))

  tibble(
    cID = c_id,
    exclusion_reason = reason,
    n_primary_studies_pet_outlier_removed = as.integer(pet_studies),
    n_primary_studies_random_effect_outlier_removed =
      as.integer(random_effect_studies)
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

empty_meta_estimate <- function(method = NA_character_) {
  list(
    method = method, estimate = NA_real_, standard_error = NA_real_,
    p_value = NA_real_, small_study_effect_p_value = NA_real_,
    tau2 = NA_real_, within_study_tau2 = NA_real_, isq = NA_real_,
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

  pet <- suppressWarnings(rma.mv(
    yi, vi, mods = ~ 1 + sei,
    random = random_effect_structure(dat),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  ))
  pet_test <- coefficient_test(pet, dat)
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
    selected_model <- suppressWarnings(rma.mv(
      yi, vi, mods = ~ 1 + vi,
      random = random_effect_structure(peese_data),
      method = "REML", test = "t", data = peese_data,
      control = list(rel.tol = 1e-8)
    ))
    selected_test <- coefficient_test(selected_model, peese_data)
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

fit_random_effects <- function(dat) {
  if (nrow(dat) <= 1) {
    result <- fit_sparse_effects(dat)
    return(result)
  }

  fit_model <- function() {
    rma.mv(
      yi, vi, mods = ~ 1,
      random = random_effect_structure(dat),
      method = "REML", test = "t", data = dat,
      control = list(rel.tol = 1e-8)
    )
  }
  model <- suppressWarnings(fit_model())
  model_test <- coefficient_test(model, dat)
  variance_components <- extract_variance_components(model)
  list(
    estimate = extract_coefficient_statistic(model_test, "intrcpt", "beta"),
    standard_error = extract_coefficient_statistic(model_test, "intrcpt", "SE"),
    p_value = extract_coefficient_statistic(model_test, "intrcpt", "p_Satt"),
    tau2 = variance_components$between_study,
    within_study_tau2 = variance_components$within_study_effect_size,
    isq = extract_total_isq(model),
    fallback = FALSE,
    model = model
  )
}

fit_random_effects_with_outlier_removal <- function(dat, cutoff = 3) {
  analysis_data <- dat

  if (primary_study_count(analysis_data) < minimum_primary_studies) {
    return(list(data = analysis_data, result = NULL))
  }

  ## Fit once to identify outliers, remove them once, and then perform the
  ## single final refit below. Outlier detection is deliberately not iterated.
  screening_fit <- fit_random_effects(analysis_data)
  standardized_residuals <- as.data.frame(
    rstandard.rma.mv(screening_fit$model)
  )$z
  keep <- is.na(standardized_residuals) |
    abs(standardized_residuals) <= cutoff
  analysis_data <- analysis_data[keep, , drop = FALSE]

  ## Let the caller exclude the complete meta-analysis before attempting the
  ## final fit (and, in particular, before requesting CR2 inference).
  if (primary_study_count(analysis_data) < minimum_primary_studies) {
    return(list(data = analysis_data, result = NULL))
  }

  final_fit <- fit_random_effects(analysis_data)
  final_fit$model <- NULL

  list(data = analysis_data, result = final_fit)
}

fit_one_meta_analysis <- function(dat) {
  if (primary_study_count(dat) < minimum_primary_studies) {
    return(drop_meta_analysis(
      dat,
      paste0("Fewer than ", minimum_primary_studies, " primary studies")
    ))
  }

  ## Identify PET-PEESE outliers from the initial PET fit, as in the original
  ## workflow, but decide estimator eligibility separately after removal.
  outlier_model <- suppressWarnings(rma.mv(
    yi, vi, mods = ~ 1 + sei,
    random = random_effect_structure(dat),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  ))
  standardized_residuals <- as.data.frame(rstandard.rma.mv(outlier_model))$z
  pet_outlier_removed_data <- dat[abs(standardized_residuals) < 3, , drop = FALSE]
  pet_studies <- primary_study_count(pet_outlier_removed_data)

  random_effect_outlier_removed <- fit_random_effects_with_outlier_removal(dat)
  random_effect_studies <- primary_study_count(random_effect_outlier_removed$data)

  if (pet_studies < minimum_primary_studies &&
      random_effect_studies < minimum_primary_studies) {
    return(drop_meta_analysis(
      dat,
      paste0(
        "Both PET-PEESE and random-effects outlier removal left fewer than ",
        minimum_primary_studies, " primary studies"
      ),
      pet_studies = pet_studies,
      random_effect_studies = random_effect_studies
    ))
  }

  results <- list(
    all_data = list(
      pet_peese = fit_pet_peese(dat),
      pet_peese_data = dat,
      random_effect = fit_random_effects(dat),
      random_effect_data = dat
    ),
    outlier_removed = list(
      pet_peese = if (pet_studies >= minimum_primary_studies) {
        fit_pet_peese(pet_outlier_removed_data)
      } else {
        empty_meta_estimate("PET-PEESE not estimable after outlier removal")
      },
      pet_peese_data = pet_outlier_removed_data,
      random_effect = if (random_effect_studies >= minimum_primary_studies) {
        random_effect_outlier_removed$result
      } else {
        empty_meta_estimate("Random effects not estimable after outlier removal")
      },
      random_effect_data = random_effect_outlier_removed$data
    )
  )
  results$all_data$random_effect$model <- NULL
  if (!is.null(results$outlier_removed$random_effect$model)) {
    results$outlier_removed$random_effect$model <- NULL
  }

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
    exclusion_reason = NA_character_,
    n_primary_studies_pet_outlier_removed = pet_studies,
    n_primary_studies_random_effect_outlier_removed = random_effect_studies,
    k_all_data = nrow(dat),
    k_outlier_removed = nrow(pet_outlier_removed_data),
    n_outliers_removed = nrow(dat) - nrow(pet_outlier_removed_data),
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
      results$outlier_removed$random_effect$within_study_tau2
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

meta_analysis_estimates <- meta_analysis_estimates %>%
  filter(is.na(exclusion_reason)) %>%
  dplyr::select(-exclusion_reason)

meta_analysis_estimates <- meta_analysis_estimates %>% arrange(cID)
write_csv(
  meta_analysis_estimates,
  file.path(derived_data_dir, "meta_analysis_estimates.csv")
)
