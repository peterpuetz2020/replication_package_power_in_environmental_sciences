## ---------------------------------------------------------
## run_exploratory_regressions.R
## ---------------------------------------------------------
## Fit the two negative-binomial models reported in Table 4 and save their
## diagnostics. Table rendering is intentionally kept in
## create_tables_and_figures.R, alongside the other manuscript tables.

source(here::here("scripts", "analysis_setup.R"))

first_existing <- function(...) {
  candidates <- c(...)
  existing <- candidates[file.exists(candidates)]
  if (!length(existing)) {
    stop(
      "None of the required regression input files exists: ",
      paste(candidates, collapse = ", "),
      ". Run scripts/create_esr_data.R and the main analysis-data workflow first."
    )
  }
  existing[[1]]
}

esr <- readxl::read_excel(first_existing(
  here::here("results", "main", "esr05_pet_peese_rstandard_half meta-average_704.xlsx"),
  here::here("data", "esr05_pet_peese_rstandard_half meta-average_704.xlsx")
))
power <- openxlsx::read.xlsx(first_existing(
  here::here("results", "main", "median_power_pps_rstandard_half meta-average_704.xlsx"),
  here::here("data", "median_power_pps_rstandard_half meta-average_704.xlsx")
))

final_nb <- merge(esr, power, by = "cID") %>%
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
  )

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

save_nb_diagnostics <- function(model, model_number) {
  diagnostic_data <- final_nb %>%
    dplyr::mutate(residual = stats::residuals(model), fitted = stats::fitted(model))
  stem <- here::here(
    "results", "supplement",
    paste0("Figure_S4_NB_Model_", model_number, "_diagnostics")
  )
  grDevices::pdf(paste0(stem, ".pdf"), width = 12, height = 12)
  old_par <- graphics::par(mfrow = c(3, 3))
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::plot(diagnostic_data$fitted, diagnostic_data$residual,
    xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs. fitted")
  graphics::abline(h = 0, col = "red", lty = 2)
  stats::qqnorm(diagnostic_data$residual, main = "Normal Q-Q plot")
  stats::qqline(diagnostic_data$residual, col = "red")
  for (variable in c("median", "lognps", "logtotall", "logjif", "pyear")) {
    graphics::plot(diagnostic_data[[variable]], diagnostic_data$residual,
      xlab = variable, ylab = "Residuals", main = paste("Residuals vs.", variable))
    graphics::abline(h = 0, col = "red", lty = 2)
  }
  for (variable in c("design_merged", "guid")) {
    graphics::boxplot(diagnostic_data$residual ~ diagnostic_data[[variable]],
      xlab = variable, ylab = "Residuals", main = paste("Residuals vs.", variable))
    graphics::abline(h = 0, col = "red", lty = 2)
  }
}

dir.create(here::here("results", "supplement"), recursive = TRUE, showWarnings = FALSE)
save_nb_diagnostics(nbMod1, 1)
save_nb_diagnostics(nbMod2, 2)
