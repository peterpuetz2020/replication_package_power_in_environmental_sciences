## Table S1 must follow the actual heterogeneity columns, including partial grids.
library(dplyr)
library(tidyr)
for (path in c("scripts/create_tables_and_figures.R", "scripts/create_supplement_tables_and_figures.R")) {
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
        as.character(expr[[2]]) %in% c("latex_escape", "write_subfield_esr_table")) {
      eval(expr, envir = .GlobalEnv)
    }
  }
}
lines <- readLines("scripts/create_supplement_tables_and_figures.R")
start <- which(startsWith(lines, "table_s1_columns <-"))[1]
end <- which(startsWith(lines, "## Table S6:"))[1] - 1L
assemble <- parse(text = lines[start:end])
subfield_levels <- c("Ecology", "Nature & Conservation")
root <- tempfile("subfield-esr-table-")
dir.create(root)
old_options <- options(warn = 2)
for (heterogeneity in list(c(.25, .5, .75), c(0, .25, .5, .75), .5, c(.125, .5))) {
  subfield_esr_results <- expand_grid(subfield = subfield_levels,
    heterogeneity_multiplier = rev(heterogeneity),
    measure = c("ESR_{0.05}^{sig}", "No. of meta-analysis", "No. of tests")) %>%
    mutate(meta_average_multiplier = .5,
      estimate = case_when(
        measure == "ESR_{0.05}^{sig}" ~ sprintf("%.3f", heterogeneity_multiplier),
        measure == "No. of meta-analysis" ~ as.character(10 + 1000 * heterogeneity_multiplier),
        TRUE ~ as.character(1000 + 10000 * heterogeneity_multiplier)),
      confidence_interval = if_else(measure == "ESR_{0.05}^{sig}", "[0.000, 1.000]", "0"))
  eval(assemble)
  headers <- paste0(100 * sort(heterogeneity), "%")
  stopifnot(identical(names(table_s1_columns), c("Subfield", "Measure", headers)))
  esr <- filter(table_s1_columns, Subfield == "Ecology", Measure == "ESR_{0.05}^{sig}")
  expected_values <- paste(sprintf("%.3f", sort(heterogeneity)), "[0.000, 1.000]")
  stopifnot(identical(unname(unlist(esr[headers])), expected_values))

  supplement_dir <- file.path(root, paste0("grid-", length(heterogeneity)))
  dir.create(supplement_dir)
  write_subfield_esr_table(table_s1_columns, number = 1)
  unzip(file.path(supplement_dir, "Table_S1.docx"), "word/document.xml", exdir = supplement_dir)
  doc <- xml2::read_xml(file.path(supplement_dir, "word", "document.xml"))
  word_rows <- xml2::xml_find_all(doc, ".//w:tbl[1]/w:tr")
  word <- do.call(rbind, lapply(word_rows, function(row) {
    xml2::xml_text(xml2::xml_find_all(row, "./w:tc"))
  }))
  expected_word <- table_s1_columns
  expected_word$Measure <- sub("^ESR", "ESO", expected_word$Measure)
  stopifnot(ncol(word) == 2 + length(heterogeneity),
    identical(as.character(word[1, ]), names(table_s1_columns)),
    all(unname(word[-1, ]) == unname(as.matrix(expected_word))))
  latex <- readLines(file.path(supplement_dir, "Table_S1.tex"))
  stopifnot(any(grepl("\\mathrm{ESO}", latex, fixed = TRUE)),
    !any(grepl("ESR", latex, fixed = TRUE)),
    any(grepl(paste0("*{", length(heterogeneity), "}{"), latex, fixed = TRUE)),
    any(grepl(paste0("\\multicolumn{", length(heterogeneity), "}{c}{Genuine heterogeneity}"), latex, fixed = TRUE)),
    any(grepl(paste0("\\cmidrule(lr){3-", length(heterogeneity) + 2L, "}"), latex, fixed = TRUE)),
    any(grepl(paste0("Subfield & Measure & ", paste(paste0(100 * sort(heterogeneity), "\\%"), collapse = " & ")),
      latex, fixed = TRUE)),
    any(grepl("Nature \\& Conservation", latex, fixed = TRUE)),
    any(endsWith(latex, paste0(" & ", paste(expected_values, collapse = " & "), " \\\\"))))
}
options(old_options)
unlink(root, recursive = TRUE)
cat("PASS: Table S1 Word/LaTeX headers and values match three-, four-, one- and two-column heterogeneity grids without warnings.\n")
