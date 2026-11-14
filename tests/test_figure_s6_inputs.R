## Regressions for stale interactive objects and grids without zero heterogeneity.
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
for (heterogeneity in list(c(.25, .5, .75), c(0, .25, .5, .75))) {
  analysis_setups <- expand_grid(meta_average_multiplier = c(.5, .75),
    heterogeneity_multiplier = heterogeneity) %>%
    mutate(setup_label = make_setup_label(meta_average_multiplier, heterogeneity_multiplier))
  eval(parameter_expression)
  p <- figure_s6_zero_heterogeneity_parameters
  stopifnot(nrow(p) == 4L, all(p$heterogeneity_multiplier == 0),
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
rows <- data.frame(estimator = c("pet_peese", "multilevel_random", "fixed_effect", "waap_wls"),
  heterogeneity_multiplier = c(.5, .5, 0, 0))
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

## Preserve the zero-heterogeneity restriction even when session objects are valid.
invalid <- inputs
invalid$all_esr_results$heterogeneity_multiplier[4] <- .25
saveRDS(invalid, path)
error_message <- tryCatch({eval(reader, envir = session); ""}, error = conditionMessage)
stopifnot(grepl("zero heterogeneity", error_message, fixed = TRUE))
unlink(path)
error_message <- tryCatch({eval(reader, envir = session); ""}, error = conditionMessage)
stopifnot(grepl("inputs are unavailable", error_message, fixed = TRUE))
unlink(root, recursive = TRUE)
cat("PASS: zero-heterogeneity comparisons survive custom grids; saved Figure S6 inputs replace stale session objects; incomplete files are diagnosed.\n")
