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
