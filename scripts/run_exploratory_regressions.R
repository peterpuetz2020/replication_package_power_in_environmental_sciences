## ---------------------------------------------------------
## run_exploratory_regressions.R
## ---------------------------------------------------------
## Fit the negative-binomial models reported in Table 4 and Tables S7-S8.
## Table rendering and supplementary diagnostic output are intentionally kept
## in their respective output scripts.

source(here::here("scripts", "analysis_setup.R"))

esr_path <- here::here(
  "data", "derived_data", "regression_data", "esr05_multilevel_random.rds"
)
power_path <- here::here(
  "data", "derived_data",
  "pps_rstandard_power_meta_0p5_heterogeneity_0_multilevel_random.rds"
)
covariate_path <- here::here(
  "data", "derived_data", "regression_data", "regression_covariates.rds"
)
if (!file.exists(esr_path) || !file.exists(power_path) || !file.exists(covariate_path)) {
  stop(
    "Missing random-effects regression inputs. Run scripts/create_analysis_data.R ",
    "before this script."
  )
}

esr <- readRDS(esr_path)
regression_covariates <- readRDS(covariate_path) %>%
  dplyr::mutate(cID = as.character(cID))
power <- readRDS(power_path) %>%
  dplyr::group_by(cID) %>%
  dplyr::summarise(
    metaID = dplyr::first(metaID),
    median = median(power, na.rm = TRUE),
    nips = dplyr::n_distinct(sID),
    esty = dplyr::first(etype),
    guid = dplyr::first(guide),
    prer = dplyr::first(prere),
    subf = dplyr::first(subfd),
    sdes = dplyr::first(sdesn),
    .groups = "drop"
  ) %>%
  dplyr::mutate(cID = as.character(cID)) %>%
  dplyr::left_join(regression_covariates, by = "cID")

final_nb <- dplyr::inner_join(
  dplyr::mutate(esr, cID = as.character(cID)), power, by = "cID"
) %>%
  dplyr::filter(esr.sig.count >= 0) %>%
  dplyr::mutate(
    med_perc = 100 * median,
    lognps = log(nips),
    logtotsig = log(tot.sig + 0.5),
    logtotall = log(tot.all),
    logjif = log(jif_5yr_wos),
    design_merged = ifelse(sdes == "experimental", "yes", "no"),
    metric = factor(dplyr::case_when(
      esty %in% c("lnRR", "log-mean ratio", "ratio") ~ "ratio",
      esty == "cohen's d" ~ "cohens_d",
      esty == "hedge's g" ~ "hedges_g",
      esty == "correlation" ~ "correlation",
      esty == "fisher's z" ~ "fishers_z",
      esty == "logOR" ~ "logOR",
      esty == "logRR" ~ "logRR",
      esty == "logHR" ~ "logHR",
      esty == "excess risk" ~ "excess_risk",
      esty == "regression coefficient" ~ "regression_coefficient",
      esty == "percentage change" ~ "percentage_change",
      esty %in% c("mean", "mean difference") ~ "mean",
      TRUE ~ NA_character_
    )),
    design_merged = factor(design_merged, levels = c("no", "yes")),
    subf = stats::relevel(factor(subf), ref = "Ecology")
  ) %>%
  tidyr::drop_na(esr.sig.count, med_perc, design_merged, guid, prer, lognps,
    logjif, pyear, metric, logtotsig, subf, metaID)

common_formula <- esr.sig.count ~ med_perc + design_merged + guid + prer +
  lognps + logjif + pyear + metric + offset(logtotsig)
nbMod1 <- MASS::glm.nb(common_formula, data = final_nb)
nbMod2 <- MASS::glm.nb(update(common_formula, . ~ . + subf), data = final_nb)

## Cluster-robust inference is reported in Table 4.
nbMod1.robu <- lmtest::coeftest(
  nbMod1, vcov. = sandwich::vcovCL(nbMod1, cluster = final_nb$metaID)
)
nbMod2.robu <- lmtest::coeftest(
  nbMod2, vcov. = sandwich::vcovCL(nbMod2, cluster = final_nb$metaID)
)

## Sensitivity specifications for Tables S7-S8. Excess-significance counts can
## be negative because they are observed minus expected counts. Negative values
## are not valid responses for a negative-binomial model, so use the same
## non-negative response restriction as the main specification above.
load_regression_multilevel_data <- function(setup_label) {
  derived_path <- here::here(
    "data", "derived_data",
    paste0("pps_rstandard_raw_", setup_label, "_multilevel_random.rds")
  )
  if (file.exists(derived_path)) {
    return(readRDS(derived_path))
  }

  input_files <- list.files(
    here::here("results", "main", "multilevel_random"),
    pattern = "\\.rds$",
    full.names = TRUE
  )
  if (length(input_files) == 0) {
    stop(
      "No multilevel random-effects inputs found. Run ",
      "compare_meta_analysis_estimators.R first."
    )
  }
  purrr::map_dfr(input_files, readRDS)
}

build_regression_data <- function(
    meta_multiplier,
    heterogeneity_multiplier,
    nonnegative_response = FALSE) {
  setup_label <- paste0(
    "meta_", gsub("\\.", "p", meta_multiplier),
    "_heterogeneity_", gsub("\\.", "p", heterogeneity_multiplier)
  )
  raw <- load_regression_multilevel_data(setup_label)
  crit <- stats::qnorm(.975)
  esr_alt <- raw %>%
    dplyr::mutate(
      mu = meta_multiplier * GE / sqrt(vi),
      sigma = sqrt(1 + heterogeneity_multiplier * tau2 / vi),
      z = abs(yi / sqrt(vi)),
      expected = stats::pnorm(-crit, mu, sigma) +
        stats::pnorm(crit, mu, sigma, lower.tail = FALSE)
    ) %>%
    dplyr::group_by(cID) %>%
    dplyr::summarise(
      tot.all = dplyr::n(),
      tot.sig = sum(z >= crit),
      esr.sig.count = tot.sig - sum(expected),
      esr.sig = dplyr::if_else(tot.sig > 0, esr.sig.count / tot.sig, 0),
      .groups = "drop"
    )
  power_alt <- raw %>%
    dplyr::mutate(
      power = 1 - stats::pnorm(crit - abs(meta_multiplier * GE) / sqrt(vi)) +
        stats::pnorm(-crit - abs(meta_multiplier * GE) / sqrt(vi))
    ) %>%
    dplyr::group_by(cID) %>%
    dplyr::summarise(
      metaID = dplyr::first(metaID),
      median = median(power, na.rm = TRUE),
      nips = dplyr::n_distinct(sID),
      esty = dplyr::first(etype),
      guid = dplyr::first(guide),
      prer = dplyr::first(prere),
      subf = dplyr::first(subfd),
      sdes = dplyr::first(sdesn),
      .groups = "drop"
    )
  covariates <- readRDS(covariate_path)
  result <- dplyr::inner_join(esr_alt, power_alt, by = "cID") %>%
    dplyr::left_join(covariates, by = "cID") %>%
    dplyr::mutate(
      med_perc = 100 * median,
      lognps = log(nips),
      logtotsig = log(tot.sig + .5),
      logtotall = log(tot.all),
      logjif = log(jif_5yr_wos),
      design_merged = factor(dplyr::if_else(sdes == "experimental", "yes", "no")),
      metric = factor(esty),
      subf = stats::relevel(factor(subf), ref = "Ecology")
    ) %>%
    tidyr::drop_na(
      esr.sig.count, med_perc, design_merged, guid, prer, lognps, logjif,
      pyear, metric, logtotsig, subf, metaID
    )
  if (nonnegative_response) {
    result <- dplyr::filter(result, esr.sig.count >= 0)
  }
  result
}

fit_nb_pair <- function(dat) {
  models <- list(
    MASS::glm.nb(common_formula, data = dat),
    MASS::glm.nb(update(common_formula, . ~ . + subf), data = dat)
  )
  robust <- purrr::map(
    models,
    ~ lmtest::coeftest(.x, vcov. = sandwich::vcovCL(.x, cluster = dat$metaID))
  )
  list(models = models, robust = robust)
}

nb_sensitivity_full <- fit_nb_pair(
  build_regression_data(1, 0, nonnegative_response = TRUE)
)
nb_sensitivity_heterogeneity <- fit_nb_pair(
  build_regression_data(.5, .25, nonnegative_response = TRUE)
)
