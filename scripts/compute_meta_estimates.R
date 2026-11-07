## Fit PET-PEESE, fixed-effect, WAAP-WLS and multilevel random-effects models with and
## without outliers. Outliers are identified once, from the random-effects
## model, so every estimator uses the same screened sample.

## WAAP-WLS is part of the common fitting stage. These helpers are also used
## by the cached-input path in main.R; no separate WAAP script is required.
## WAAP-WLS point estimation with the same sID clustering, CR2 correction and
## Satterthwaite inference as the project's fixed-effect estimator.
fit_waap_wls <- function(dat, alpha = .05, adequate_power = .8) {
  if (!all(c("yi", "vi", "sID") %in% names(dat)) || nrow(dat) < 2L ||
      any(!is.finite(dat$yi) | !is.finite(dat$vi) | dat$vi <= 0) || anyNA(dat$sID)) {
    stop("WAAP-WLS requires at least two valid estimates, positive variances and study IDs.")
  }
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(adequate_power) != 1 || !is.finite(adequate_power) || adequate_power <= 0 || adequate_power >= 1) {
    stop("Invalid alpha or adequate-power threshold.")
  }
  wls_mean <- weighted.mean(dat$yi, 1 / dat$vi)
  lambda <- abs(wls_mean) / sqrt(dat$vi)
  power <- pnorm(-qnorm(1 - alpha / 2) - lambda) +
    pnorm(lambda - qnorm(1 - alpha / 2))
  adequate <- power >= adequate_power
  use_waap <- sum(adequate) >= 2L
  selected <- if (use_waap) dat[adequate, , drop = FALSE] else dat
  model <- suppressWarnings(metafor::rma.mv(yi, vi, mods = ~1,
    method = "FE", test = "t", data = selected))
  result <- list(method = if (use_waap) "WAAP" else "WLS fallback",
    estimate = unname(coef(model)[1]), standard_error = NA_real_, p_value = NA_real_,
    df = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
    n_adequate = sum(adequate), n_fit = nrow(selected),
    n_studies_fit = length(unique(selected$sID)), wls_mean = wls_mean,
    adequate = adequate, power = power, fallback = !use_waap,
    inference_status = "fewer_than_two_primary_study_clusters")
  if (result$n_studies_fit < 2L) return(result)
  test <- tryCatch(clubSandwich::coef_test(model, vcov = clubSandwich::vcovCR(
    model, cluster = droplevels(factor(selected$sID)), type = "CR2")), error = function(e) e)
  if (inherits(test, "error")) {
    result$inference_status <- paste0("CR2 unavailable: ", conditionMessage(test))
    return(result)
  }
  result$standard_error <- as.numeric(test$SE[1])
  result$p_value <- as.numeric(test$p_Satt[1])
  result$df <- as.numeric(test$df_Satt[1])
  if (is.finite(result$standard_error) && is.finite(result$df) && result$df > 0) {
    margin <- qt(.975, result$df) * result$standard_error
    result$ci_lower <- result$estimate - margin
    result$ci_upper <- result$estimate + margin
    result$inference_status <- "ok_CR2_Satterthwaite"
  } else result$inference_status <- "unavailable_CR2_Satterthwaite"
  result
}

compute_waap_wls <- function(root = ".", recreate = FALSE) {
  summaries <- list()
  for (variant in c("outliers_removed", "all_data")) {
    input_dir <- file.path(root, "data", "derived_data",
      if (variant == "all_data") "multilevel_random_all_data" else "multilevel_random")
    files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
    if (!length(files)) stop("Fit the shared random-effects samples before WAAP-WLS: ", input_dir)
    output_dir <- file.path(root, "data", "derived_data",
      if (variant == "all_data") "waap_wls_all_data" else "waap_wls")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    stale <- setdiff(list.files(output_dir, pattern = "\\.rds$"), basename(files))
    if (length(stale)) unlink(file.path(output_dir, stale))
    for (i in seq_along(files)) {
      dat <- readRDS(files[i])
      key <- list(version = 2L, source_md5 = unname(tools::md5sum(files[i])),
                  alpha = .05, adequate_power = .8)
      output <- file.path(output_dir, basename(files[i]))
      cached <- if (!recreate && file.exists(output)) readRDS(output) else NULL
      # Version 2 changes only heterogeneity storage, not WAAP fits or inference.
      previous_key <- key; previous_key$version <- 1L
      if (!recreate && !is.null(cached) && (identical(attr(cached, "waap_cache_key"), key) ||
          identical(attr(cached, "waap_cache_key"), previous_key))) {
        summary <- attr(cached, "waap_summary")
        if (!identical(attr(cached, "waap_cache_key"), key)) {
          cached$tau2 <- 0; cached$isq <- 0
          cached$heterogeneity_source <- "not_estimated_zero_assumed"
          summary$tau2_for_counterfactual <- 0
          summary$heterogeneity_source <- "not_estimated_zero_assumed"
          attr(cached, "waap_cache_key") <- key
          attr(cached, "waap_summary") <- summary
          saveRDS(cached, output)
        }
      } else {
        if (nrow(dat) < 2 || any(!is.finite(dat$yi) | !is.finite(dat$vi) | dat$vi <= 0) ||
            length(unique(dat$sID)) < 5L) {
          if (file.exists(output)) unlink(output)
          next
        }
        fit <- fit_waap_wls(dat)
        dat$GE <- fit$estimate
        dat$sig_overall <- fit$p_value
        # WAAP does not estimate variance components. Its counterfactuals use
        # zero genuine heterogeneity; no RE components are borrowed.
        dat$tau2 <- 0; dat$isq <- 0
        dat$heterogeneity_source <- "not_estimated_zero_assumed"
        dat$waap_selected_method <- fit$method
        dat$waap_n_adequate <- fit$n_adequate
        dat$waap_n_fit <- fit$n_fit
        dat$waap_n_studies_fit <- fit$n_studies_fit
        dat$waap_inference_status <- fit$inference_status
        summary <- data.frame(cID = as.character(dat$cID[1]), variant = variant,
          selected_method = fit$method, estimate = fit$estimate,
          standard_error = fit$standard_error, df = fit$df, p_value = fit$p_value,
          ci_lower = fit$ci_lower, ci_upper = fit$ci_upper,
          wls_reference_mean = fit$wls_mean, n_input = nrow(dat),
          n_adequate = fit$n_adequate, n_fit = fit$n_fit,
          n_studies_fit = fit$n_studies_fit, inference_status = fit$inference_status,
          tau2_for_counterfactual = unique(dat$tau2)[1],
          heterogeneity_source = dat$heterogeneity_source[1])
        attr(dat, "waap_cache_key") <- key
        attr(dat, "waap_summary") <- summary
        saveRDS(dat, output)
      }
      summaries[[length(summaries) + 1L]] <- summary
      if (i %% 50L == 0L) message("WAAP-WLS ", variant, ": ", i, "/", length(files))
    }
  }
  summary <- do.call(rbind, summaries)
  summary_path <- file.path(root, "data", "derived_data", "WAAP_WLS_estimates.xlsx")
  dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(list(Estimates = summary, Methods = data.frame(note = c(
    "Initial inverse-variance WLS mean; select estimates with exact normal two-sided power >= 0.80 at alpha 0.05.",
    "Preserve signed means; assess power using the mean's magnitude, before applying meta-average multipliers.",
    "At least two adequate estimates: WAAP; otherwise full-sample WLS fallback.",
    "Inference: primary-study sID clusters, CR2 and Satterthwaite, matching fixed effects. Unavailable inference is flagged, not replaced by independent-observation SEs.",
    "WAAP fits on the adequate subset; all shared-sample empirical effects remain in ESR denominators.",
    "WAAP does not estimate tau2. Like fixed effects, its counterfactuals are restricted to zero genuine heterogeneity; zero is assumed, not estimated."
  ))), summary_path, overwrite = TRUE)
  # Remove the obsolete generated supplement workbook after its replacement exists.
  legacy_path <- file.path(root, "results", "supplement", "WAAP_WLS_estimates.xlsx")
  if (file.exists(legacy_path)) unlink(legacy_path)
  invisible(summary)
}

## main.R uses this fast path when the other estimator inputs already exist.
## Directly sourcing this script runs the complete fitting stage by default.
if (!isTRUE(get0("waap_only", ifnotfound = FALSE))) {
library(metafor)
library(clubSandwich)
library(tidyverse)
library(readxl)
library(foreach)
library(doParallel)
library(orchaRd)
library(here)

if (!exists("n_cores")) n_cores <- 6
if (!exists("recreate_meta_analysis_estimates")) {
  recreate_meta_analysis_estimates <- FALSE
}
minimum_primary_studies <- 5L
meta_estimate_cache_version <- 3L
if (!exists("new_progress_bar")) {
  new_progress_bar <- function(total, label) {
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

master_data_path <- here("data", "MasterData.xlsx")
meta_analysis_data <- read_excel(master_data_path)
meta_analyses <- split(meta_analysis_data, meta_analysis_data$cID)

derived_data_dir <- here("data", "derived_data")
output_dirs <- c(
  pet_peese_outlier_removed = here(derived_data_dir, "pet_peese_rstandard"),
  pet_peese_all_data = here(derived_data_dir, "pet_peese_all_data"),
  fixed_effect_outlier_removed = here(derived_data_dir, "fixed_effect"),
  fixed_effect_all_data = here(derived_data_dir, "fixed_effect_all_data"),
  random_effect_outlier_removed = here(derived_data_dir, "multilevel_random"),
  random_effect_all_data = here(derived_data_dir, "multilevel_random_all_data")
)
invisible(lapply(
  c(output_dirs, derived_data_dir), dir.create,
  recursive = TRUE, showWarnings = FALSE
))

## Cache one compact summary per meta-analysis as soon as that analysis
## finishes. This makes interrupted runs resumable and avoids refitting all 708
## analyses on every invocation. The manifest invalidates the cache whenever
## the master workbook or cache schema changes; users can also force a complete
## rebuild with recreate_meta_analysis_estimates <- TRUE.
summary_cache_dir <- file.path(derived_data_dir, "meta_analysis_summary_cache")
dir.create(summary_cache_dir, recursive = TRUE, showWarnings = FALSE)
cache_manifest_path <- file.path(summary_cache_dir, "manifest.rds")
cache_manifest <- list(
  version = meta_estimate_cache_version,
  master_data_md5 = unname(tools::md5sum(master_data_path)),
  minimum_primary_studies = minimum_primary_studies
)
stored_manifest <- if (file.exists(cache_manifest_path)) {
  readRDS(cache_manifest_path)
} else {
  NULL
}
cache_is_current <- !recreate_meta_analysis_estimates &&
  identical(stored_manifest, cache_manifest)
if (!cache_is_current) {
  unlink(summary_cache_dir, recursive = TRUE)
  dir.create(summary_cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache_manifest, cache_manifest_path)
}

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

drop_meta_analysis <- function(dat, reason,
                               screened_studies = NA_integer_) {
  c_id <- as.character(dat$cID[[1]])

  ## A rerun must not leave output files from an earlier, less restrictive run.
  invisible(lapply(output_dirs, function(output_dir) {
    unlink(file.path(output_dir, paste0("meta_", c_id, ".rds")))
  }))

  tibble(
    cID = c_id,
    exclusion_reason = reason,
    n_primary_studies_outlier_removed = as.integer(screened_studies)
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
    p_value = NA_real_,
    tau2 = NA_real_, within_study_tau2 = NA_real_, isq = NA_real_,
    fallback = TRUE
  )
}

fit_pet_peese <- function(dat, pet = NULL, pet_test = NULL) {
  if (nrow(dat) <= 1) {
    result <- fit_sparse_effects(dat)
    result$method <- "PET-PEESE not estimable"
    return(result)
  }

  if (is.null(pet)) {
    pet <- suppressWarnings(rma.mv(
      yi, vi, mods = ~ 1 + sei,
      random = random_effect_structure(dat),
      method = "REML", test = "t", data = dat,
      control = list(rel.tol = 1e-8)
    ))
  }
  if (is.null(pet_test)) pet_test <- coefficient_test(pet, dat)
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
    tau2 = extract_variance_components(selected_model)$between_study /
      outcome_scale^2,
    within_study_tau2 =
      extract_variance_components(selected_model)$within_study_effect_size /
      outcome_scale^2,
    isq = extract_total_isq(selected_model)
  )
}

fit_egger <- function(dat) {
  egger <- suppressWarnings(rma.mv(
    yi, vi, mods = ~ sei,
    random = random_effect_structure(dat),
    method = "REML", test = "t", data = dat,
    control = list(rel.tol = 1e-8)
  ))
  egger_test <- coefficient_test(egger, dat)

  list(
    slope = extract_coefficient_statistic(egger_test, "sei", "beta"),
    standard_error = extract_coefficient_statistic(egger_test, "sei", "SE"),
    p_value = extract_coefficient_statistic(egger_test, "sei", "p_Satt")
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

fit_fixed_effects <- function(dat) {
  model <- suppressWarnings(rma.mv(
    yi, vi, mods = ~ 1, method = "FE", test = "t", data = dat
  ))
  model_test <- coefficient_test(model, dat)
  list(
    estimate = extract_coefficient_statistic(model_test, "intrcpt", "beta"),
    standard_error = extract_coefficient_statistic(model_test, "intrcpt", "SE"),
    p_value = extract_coefficient_statistic(model_test, "intrcpt", "p_Satt"),
    tau2 = 0,
    within_study_tau2 = 0,
    isq = 0,
    fallback = FALSE,
    model = model
  )
}

fit_random_effects_with_outlier_removal <- function(
    dat, cutoff = 3, screening_fit = NULL) {
  analysis_data <- dat

  if (primary_study_count(analysis_data) < minimum_primary_studies) {
    return(list(data = analysis_data, result = NULL))
  }

  ## Fit once to identify outliers, remove them once, and then perform the
  ## single final refit below. Outlier detection is deliberately not iterated.
  if (is.null(screening_fit)) screening_fit <- fit_random_effects(analysis_data)
  standardized_residuals <- as.data.frame(
    rstandard.rma.mv(screening_fit$model)
  )$z
  keep <- is.na(standardized_residuals) |
    abs(standardized_residuals) <= cutoff
  analysis_data <- analysis_data[keep, , drop = FALSE]

  ## When screening removes nothing, the screening model is already the exact
  ## final all-data model. Reusing it avoids an unnecessary REML fit and CR2
  ## calculation, which is a common case across hundreds of meta-analyses.
  if (all(keep)) {
    reused_fit <- screening_fit
    reused_fit$model <- NULL
    return(list(data = analysis_data, result = reused_fit))
  }

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

  random_effect_all_data <- fit_random_effects(dat)
  random_effect_outlier_removed <- fit_random_effects_with_outlier_removal(
    dat, screening_fit = random_effect_all_data
  )
  screened_data <- random_effect_outlier_removed$data
  screened_studies <- primary_study_count(screened_data)
  data_changed <- nrow(screened_data) != nrow(dat)

  pet_peese_all_data <- fit_pet_peese(dat)
  fixed_effect_all_data <- fit_fixed_effects(dat)

  if (screened_studies >= minimum_primary_studies) {
    pet_peese_outlier_removed <- if (data_changed) {
      fit_pet_peese(screened_data)
    } else {
      pet_peese_all_data
    }
    fixed_effect_outlier_removed <- if (data_changed) {
      fit_fixed_effects(screened_data)
    } else {
      fixed_effect_all_data
    }
    egger <- fit_egger(screened_data)
  } else {
    pet_peese_outlier_removed <- empty_meta_estimate(
      "PET-PEESE not estimable after outlier removal"
    )
    fixed_effect_outlier_removed <- empty_meta_estimate(
      "Fixed effect not estimable after outlier removal"
    )
    egger <- list(slope = NA_real_, standard_error = NA_real_, p_value = NA_real_)
  }

  results <- list(
    all_data = list(
      pet_peese = pet_peese_all_data,
      pet_peese_data = dat,
      fixed_effect = fixed_effect_all_data,
      fixed_effect_data = dat,
      random_effect = random_effect_all_data,
      random_effect_data = dat
    ),
    outlier_removed = list(
      pet_peese = pet_peese_outlier_removed,
      pet_peese_data = screened_data,
      fixed_effect = fixed_effect_outlier_removed,
      fixed_effect_data = screened_data,
      random_effect = if (screened_studies >= minimum_primary_studies) {
        random_effect_outlier_removed$result
      } else {
        empty_meta_estimate("Random effects not estimable after outlier removal")
      },
      random_effect_data = screened_data
    )
  )
  results$all_data$random_effect$model <- NULL
  results$all_data$fixed_effect$model <- NULL
  if (!is.null(results$outlier_removed$random_effect$model)) {
    results$outlier_removed$random_effect$model <- NULL
  }
  if (!is.null(results$outlier_removed$fixed_effect$model)) {
    results$outlier_removed$fixed_effect$model <- NULL
  }

  for (variant in names(results)) {
    result <- results[[variant]]
    saveRDS(
      make_effect_data(
        result$fixed_effect_data, result$fixed_effect
      ),
      file.path(
        output_dirs[[paste0("fixed_effect_", variant)]],
        paste0("meta_", dat$cID[[1]], ".rds")
      )
    )
    saveRDS(
      make_effect_data(
        result$pet_peese_data, result$pet_peese
      ),
      file.path(
        output_dirs[[paste0("pet_peese_", variant)]],
        paste0("meta_", dat$cID[[1]], ".rds")
      )
    )
    saveRDS(
      make_effect_data(
        result$random_effect_data, result$random_effect,
        if (variant == "outlier_removed") egger$p_value else NA_real_
      ),
      file.path(
        output_dirs[[paste0("random_effect_", variant)]],
        paste0("meta_", dat$cID[[1]], ".rds")
      )
    )
  }

  tibble(
    cID = as.character(dat$cID[[1]]),
    exclusion_reason = NA_character_,
    n_primary_studies_outlier_removed = screened_studies,
    k_all_data = nrow(dat),
    k_outlier_removed = nrow(screened_data),
    n_outliers_removed = nrow(dat) - nrow(screened_data),
    egger_slope_outlier_removed = egger$slope,
    egger_slope_se_outlier_removed = egger$standard_error,
    egger_p_value_outlier_removed = egger$p_value,
    pet_peese_method_all_data = results$all_data$pet_peese$method,
    pet_peese_method_outlier_removed = results$outlier_removed$pet_peese$method,
    pet_peese_estimate_all_data = results$all_data$pet_peese$estimate,
    pet_peese_estimate_outlier_removed = results$outlier_removed$pet_peese$estimate,
    random_effect_estimate_all_data = results$all_data$random_effect$estimate,
    random_effect_estimate_outlier_removed = results$outlier_removed$random_effect$estimate,
    fixed_effect_estimate_all_data = results$all_data$fixed_effect$estimate,
    fixed_effect_estimate_outlier_removed =
      results$outlier_removed$fixed_effect$estimate,
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

summary_cache_path <- function(c_id) {
  file.path(summary_cache_dir, paste0("meta_", c_id, ".rds"))
}

write_summary_cache <- function(result, path) {
  temporary_path <- tempfile(
    pattern = paste0(basename(path), "_"),
    tmpdir = dirname(path), fileext = ".tmp"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(result, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop("Could not move completed meta-analysis cache to ", path)
  }
  invisible(path)
}

expected_effect_paths <- function(c_id) {
  file.path(output_dirs, paste0("meta_", c_id, ".rds"))
}

cached_meta_analysis_is_complete <- function(c_id) {
  path <- summary_cache_path(c_id)
  if (!file.exists(path)) return(FALSE)

  summary <- tryCatch(readRDS(path), error = function(...) NULL)
  if (is.null(summary) || !"exclusion_reason" %in% names(summary)) {
    return(FALSE)
  }
  is_excluded <- !is.na(summary$exclusion_reason[[1]])
  is_excluded || all(file.exists(expected_effect_paths(c_id)))
}

fit_and_cache_meta_analysis <- function(dat) {
  result <- fit_one_meta_analysis(dat)
  write_summary_cache(
    result, summary_cache_path(as.character(dat$cID[[1]]))
  )
  result
}

meta_analysis_ids <- names(meta_analyses)
completed_ids <- meta_analysis_ids[vapply(
  meta_analysis_ids, cached_meta_analysis_is_complete, logical(1)
)]
pending_ids <- setdiff(meta_analysis_ids, completed_ids)
message(
  "Reusing ", length(completed_ids), " cached meta-analyses; fitting ",
  length(pending_ids), " of ", length(meta_analyses), "."
)

if (length(pending_ids) > 0L) {
  cluster <- makeCluster(min(n_cores, length(pending_ids)))
  registerDoParallel(cluster)
  invisible(tryCatch(
    foreach(
      dat = meta_analyses[pending_ids],
      .packages = c("metafor", "clubSandwich", "dplyr", "tibble", "orchaRd"),
      ## Dynamic scheduling prevents a few large meta-analyses from leaving
      ## otherwise idle workers near the end of the run. Each result is cached
      ## immediately; the small summaries are assembled from disk below.
      .options.snow = list(preschedule = FALSE)
    ) %dopar% fit_and_cache_meta_analysis(dat),
    finally = stopCluster(cluster)
  ))
}

## Read summaries from disk so cached and newly fitted analyses follow exactly
## the same aggregation path and an interrupted run needs no special handling.
all_meta_analysis_summaries <- purrr::map_dfr(
  meta_analysis_ids,
  ~ readRDS(summary_cache_path(.x))
)
excluded_meta_analyses <- all_meta_analysis_summaries %>%
  dplyr::filter(!is.na(exclusion_reason))
readr::write_csv(
  excluded_meta_analyses,
  file.path(derived_data_dir, "excluded_meta_analyses.csv")
)

meta_analysis_estimates <- all_meta_analysis_summaries %>%
  dplyr::filter(is.na(exclusion_reason)) %>%
  dplyr::select(-exclusion_reason)

meta_analysis_estimates <- meta_analysis_estimates %>% dplyr::arrange(cID)
readr::write_csv(
  meta_analysis_estimates,
  file.path(derived_data_dir, "meta_analysis_estimates.csv")
)

report_outlier_removal <- function(estimator, removed_column) {
  removed <- sum(meta_analysis_estimates[[removed_column]], na.rm = TRUE)
  all_estimates <- sum(meta_analysis_estimates$k_all_data, na.rm = TRUE)
  percentage <- if (all_estimates == 0) 0 else 100 * removed / all_estimates

  message(sprintf(
    "%s outlier removal: %d of %d primary estimates removed (%.2f%%).",
    estimator, removed, all_estimates, percentage
  ))
}

## The single random-effects screen supplies the data for every estimator.
report_outlier_removal("Random-effects model", "n_outliers_removed")


}

## Reuse the common all-data and outlier-screened samples for WAAP-WLS.
waap_wls_estimates <- compute_waap_wls(here::here(),
  recreate = isTRUE(get0("recreate_meta_analysis_estimates", ifnotfound = FALSE)))
