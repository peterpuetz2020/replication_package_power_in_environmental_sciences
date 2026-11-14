## Run after main.R, from a directory containing its generated outputs.
source("scripts/functions.R")
args <- commandArgs(trailingOnly = TRUE)
expected_iterations <- if (length(args)) as.integer(args[[1]]) else 1000L
stopifnot(length(expected_iterations) == 1L, !is.na(expected_iterations),
          expected_iterations > 0L)
for (heterogeneity in c(0, .25, .5, .75)) {
  label <- paste0("meta_0p5_heterogeneity_", gsub("\\.", "p", heterogeneity))
  cached <- readRDS(file.path("data", "derived_data", paste0("p_tab_ci_", label,
                        "_multilevel_random_outliers_removed.rds")))
  stopifnot(nrow(cached$result[[1]]) == expected_iterations)
}
esr <- readRDS("data/derived_data/Figure_S6_inputs.rds")$all_esr_results
stopifnot(setequal(unique(esr$meta_average_multiplier), c(.25, .5, .75, 1)))

power <- lapply(c(.5, .75, 1), function(m) {
  label <- paste0("meta_", gsub("\\.", "p", m), "_heterogeneity_0")
  readRDS(file.path("data", "derived_data", paste0("pps_rstandard_power_", label,
                   "_multilevel_random.rds")))
})
stopifnot(identical(power[[1]]$cID, power[[2]]$cID),
          identical(power[[2]]$cID, power[[3]]$cID),
          identical(is.finite(power[[1]]$power), is.finite(power[[2]]$power)),
          identical(is.finite(power[[2]]$power), is.finite(power[[3]]$power)),
          all(power[[2]]$power >= power[[1]]$power - 1e-12, na.rm = TRUE),
          all(power[[2]]$power <= power[[3]]$power + 1e-12, na.rm = TRUE),
          any(power[[2]]$power > power[[1]]$power + 1e-6, na.rm = TRUE))
required <- c(
  paste0("results/main/Table_", 1:3, ".docx"),
  paste0("results/supplement/Table_S", 1:7, ".docx"),
  paste0("results/supplement/Figure_S", 1:6, ".pdf"),
  paste0("results/supplement/Figure_S", 1:6, ".svg")
)
obsolete <- c("results/main/Table_4.docx", "results/main/Table_4.tex",
              paste0("results/supplement/Table_S", rep(8:9, each = 2), c(".docx", ".tex")),
              "data/derived_data/counterfactual_comparison_0p5_vs_1.csv",
              paste0("results/supplement/Figure_S", 7:8, ".pdf"),
              paste0("results/supplement/Figure_S", 7:8, ".svg"),
              "data/derived_data/regression_data")
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0),
          !any(file.exists(obsolete)))
read_word_table <- function(path) {
  doc_dir <- tempfile("analysis-docx-")
  dir.create(doc_dir)
  on.exit(unlink(doc_dir, recursive = TRUE))
  unzip(path, "word/document.xml", exdir = doc_dir)
  doc <- xml2::read_xml(file.path(doc_dir, "word/document.xml"))
  rows <- xml2::xml_find_all(doc, ".//w:tbl[1]/w:tr")
  do.call(rbind, lapply(rows, function(row) xml2::xml_text(xml2::xml_find_all(row,"./w:tc"))))
}
effect_sizes <- read_word_table("results/supplement/Table_S6.docx")
small_studies <- read_word_table("results/supplement/Table_S7.docx")
stopifnot(
  identical(dim(effect_sizes), c(13L, 9L)),
  identical(effect_sizes[1, ], c("No.", "Effect size", "ECO", "ENC", "ENE", "HTM", "MPL", "NLC", "WST")),
  identical(dim(small_studies), c(9L, 3L)),
  identical(small_studies[1, ], c("Subfield", "No. of meta-analyses", "Small-study effects (%)")),
  small_studies[2, 1] == "All meta-analyses",
  sum(as.numeric(effect_sizes[-1, 3:9])) == as.numeric(small_studies[2, 2]),
  all(colSums(matrix(as.numeric(effect_sizes[-1, 3:9]), ncol = 7)) ==
        as.numeric(small_studies[3:9, 2])))
effect_size_latex <- readLines("results/supplement/Table_S6.tex")
for (label in c("Cohen's $d$", "Correlation coefficient ($r$)",
                "Fisher's $z$", "Hedges's $g$",
                "Regression coefficient ($\\beta$)")) {
  stopifnot(any(grepl(label, effect_size_latex, fixed = TRUE)))
}
stopifnot(!any(grepl("β", effect_size_latex, fixed = TRUE)))
for (number in 2:4) {
  path <- paste0("results/supplement/Table_S",number,".docx")
  doc_dir <- tempfile("power-docx-"); dir.create(doc_dir)
  unzip(path,"word/document.xml",exdir=doc_dir)
  text <- xml2::xml_text(xml2::read_xml(file.path(doc_dir,"word/document.xml")))
  multiplier <- c(.25,.75,1)[number-1]
  stopifnot(grepl(paste0("Statistical power at ",multiplier," x the meta-average"),text,fixed=TRUE))
  unlink(doc_dir,recursive=TRUE)
}
cat("PASS: descriptive Tables S6/S7, LaTeX math, consecutive outputs, power sensitivities and counterfactual iterations verified\n")
