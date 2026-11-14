## Saved-input and custom-grid regressions for all four Figure S6 panels.
library(dplyr)
library(tidyr)

assignment <- function(path, name) {
  matches <- Filter(function(expr) is.call(expr) &&
    identical(expr[[1]], as.name("<-")) && identical(expr[[2]], as.name(name)),
    as.list(parse(path)))
  stopifnot(length(matches) == 1L)
  matches[[1]]
}
eval(assignment("scripts/analysis_setup.R", "make_setup_label"))
parameter_expression <- assignment("scripts/create_tables_and_figures.R",
  "figure_s6_zero_heterogeneity_parameters")
table_parameter_expressions <- Filter(function(expr) is.call(expr) &&
  identical(expr[[1]], as.name("<-")) &&
  identical(expr[[2]], as.name("table_2_parameters")),
  as.list(parse("scripts/create_tables_and_figures.R")))
stopifnot(length(table_parameter_expressions) == 2L)
meta_analysis_estimators <- c("pet_peese", "multilevel_random")
for (heterogeneity in list(c(.25, .5, .75), c(0, .25, .5, .75))) {
  analysis_setups <- expand_grid(meta_average_multiplier = c(.5, .75, 1),
    heterogeneity_multiplier = heterogeneity) %>%
    mutate(setup_label = make_setup_label(meta_average_multiplier, heterogeneity_multiplier))
  eval(assignment("scripts/create_tables_and_figures.R", "figure_s6_setups"))
  stopifnot(setequal(figure_s6_setups$meta_average_multiplier, c(.25, .5, .75, 1)),
    setequal(figure_s6_setups$heterogeneity_multiplier, heterogeneity),
    nrow(figure_s6_setups) == 4L * length(heterogeneity))
  eval(parameter_expression)
  eval(table_parameter_expressions[[1]])
  eval(table_parameter_expressions[[2]])
  stopifnot(nrow(table_2_parameters) == 8L * length(heterogeneity) + 8L,
    setequal(table_2_parameters$meta_average_multiplier, c(.25, .5, .75, 1)),
    !anyDuplicated(table_2_parameters))
  p <- figure_s6_zero_heterogeneity_parameters
  stopifnot(nrow(p) == 8L, all(p$heterogeneity_multiplier == 0),
    !anyDuplicated(paste(p$meta_average_multiplier, p$estimator)),
    setequal(p$estimator, c("fixed_effect", "waap_wls")),
    identical(p$setup_label, make_setup_label(p$meta_average_multiplier, 0)))
}

## Execute the actual saved-input reader, without rendering or fitting models.
lines <- readLines("scripts/create_supplement_tables_and_figures.R")
start <- which(startsWith(lines, "## Figure S6:"))[1]
end <- which(startsWith(lines, "write.csv(") & seq_along(lines) > start)[1] - 1L
reader <- parse(text = lines[start:end])
root <- tempfile("figure-s6-inputs-")
dir.create(file.path(root, "data", "derived_data"), recursive = TRUE)
path <- file.path(root, "data", "derived_data", "Figure_S6_inputs.rds")
rows <- expand_grid(meta_average_multiplier = c(.25, .5, .75, 1),
  estimator = c("pet_peese", "multilevel_random", "fixed_effect", "waap_wls")) %>%
  mutate(heterogeneity_multiplier = if_else(estimator %in% c("fixed_effect", "waap_wls"), 0, .5),
    measure = "ESR_{0.05}^{sig}")
inputs <- list(all_combination_results = rows, all_esr_results = rows)
saveRDS(inputs, path)
session <- new.env(parent = .GlobalEnv)
session$here <- function(...) file.path(root, ...)
session$all_esr_results <- subset(rows, estimator != "waap_wls")
session$all_combination_results <- session$all_esr_results
eval(reader, envir = session)
stopifnot(identical(session$all_esr_results, rows),
  identical(session$all_combination_results, rows))

## A truly incomplete saved file must request rebuilding the figure inputs.
incomplete <- lapply(inputs, function(x) subset(x, estimator != "waap_wls"))
saveRDS(incomplete, path)
error_message <- tryCatch({eval(reader, envir = session); ""}, error = conditionMessage)
stopifnot(grepl("waap_wls", error_message, fixed = TRUE),
  grepl("create_tables_and_figures.R", error_message, fixed = TRUE),
  !grepl("compute_meta_estimates.R", error_message, fixed = TRUE))

## Reproduce the uploaded three-panel figure: missing 0.25 must be diagnosed.
incomplete_panels <- lapply(inputs, function(x) subset(x, meta_average_multiplier != .25))
saveRDS(incomplete_panels, path)
error_message <- tryCatch({eval(reader, envir = session); ""}, error = conditionMessage)
stopifnot(grepl("0.25", error_message, fixed = TRUE),
  grepl("create_tables_and_figures.R", error_message, fixed = TRUE))

## One missing estimator in a panel must not silently pass the global check.
incomplete_panels <- inputs
incomplete_panels$all_esr_results <- subset(rows,
  !(meta_average_multiplier == .25 & estimator == "waap_wls"))
saveRDS(incomplete_panels, path)
error_message <- tryCatch({eval(reader, envir = session); ""}, error = conditionMessage)
stopifnot(grepl("0.25 / waap_wls", error_message, fixed = TRUE))

## Preserve the zero-heterogeneity restriction even when session objects are valid.
invalid <- inputs
invalid$all_esr_results$heterogeneity_multiplier[which(rows$estimator == "waap_wls")[1]] <- .25
saveRDS(invalid, path)
error_message <- tryCatch({eval(reader, envir = session); ""}, error = conditionMessage)
stopifnot(grepl("zero heterogeneity", error_message, fixed = TRUE))
unlink(path)
error_message <- tryCatch({eval(reader, envir = session); ""}, error = conditionMessage)
stopifnot(grepl("inputs are unavailable", error_message, fixed = TRUE))
unlink(root, recursive = TRUE)
cat("PASS: all four Figure S6 multipliers survive partial grids; zero heterogeneity is preserved; stale objects and missing panel/estimator inputs are diagnosed.\n")
