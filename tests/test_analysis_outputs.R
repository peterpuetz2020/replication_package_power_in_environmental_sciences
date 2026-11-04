## Run after main.R, from a directory containing its generated outputs.
source("scripts/functions.R")
args <- commandArgs(trailingOnly = TRUE)
expected_iterations <- if (length(args)) as.integer(args[[1]]) else 1000L
stopifnot(length(expected_iterations) == 1L, !is.na(expected_iterations),
          expected_iterations > 0L)
comparison <- read.csv("data/derived_data/counterfactual_comparison_0p5_vs_1.csv")
stopifnot(nrow(comparison) == 52L,
          setequal(unique(comparison$heterogeneity_multiplier), c(0, .25, .5, .75)),
          all(comparison$ci_lower <= comparison$ci_upper))
for (h in c(0, .25, .5, .75)) {
  rows <- comparison[comparison$heterogeneity_multiplier == h, ]
  stopifnot(identical(rows$p_value_interval, p_value_bin_labels()),
            abs(sum(rows$probability_0p5) - 1) < 1e-10,
            abs(sum(rows$probability_1) - 1) < 1e-10,
            max(abs(rows$difference - (rows$probability_0p5 - rows$probability_1))) < 1e-12)
  draws <- lapply(c(.5, 1), function(m) {
    label <- paste0("meta_", gsub("\\.", "p", m), "_heterogeneity_", gsub("\\.", "p", h))
    path <- file.path("data", "derived_data", paste0("p_tab_ci_", label,
                       "_multilevel_random_outliers_removed.rds"))
    readRDS(path)$result[[1]]
  })
  stopifnot(nrow(draws[[1]]) == expected_iterations,
            nrow(draws[[2]]) == expected_iterations)
  expected_ci <- apply(draws[[1]] - draws[[2]], 2, quantile, probs = c(.025, .975))
  stopifnot(max(abs(rows$ci_lower - expected_ci[1, ])) < 1e-12,
            max(abs(rows$ci_upper - expected_ci[2, ])) < 1e-12)
}
esr <- readRDS("data/derived_data/Figure_S4_inputs.rds")$all_esr_results
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
  paste0("results/supplement/Table_S", 1:8, ".docx"),
  paste0("results/supplement/Figure_S", 1:5, ".pdf"),
  paste0("results/supplement/Figure_S", 1:5, ".svg")
)
obsolete <- c("results/main/Table_4.docx", "results/main/Table_4.tex",
              "results/supplement/Table_S9.docx", "results/supplement/Table_S9.tex",
              paste0("results/supplement/Figure_S", 6:8, ".pdf"),
              paste0("results/supplement/Figure_S", 6:8, ".svg"),
              "data/derived_data/regression_data")
stopifnot(all(file.exists(required)), all(file.info(required)$size > 0),
          !any(file.exists(obsolete)))
doc_dir <- tempfile("counterfactual-docx-")
dir.create(doc_dir)
unzip("results/supplement/Table_S6.docx", "word/document.xml", exdir = doc_dir)
text <- xml2::xml_text(xml2::read_xml(file.path(doc_dir, "word/document.xml")))
stopifnot(grepl("0.5 x meta-average", text, fixed = TRUE),
          grepl("1 x meta-average", text, fixed = TRUE),
          grepl("0.001 < p < 0.01", text, fixed = TRUE))
unlink(doc_dir, recursive = TRUE)
dir.create(doc_dir)
unzip("results/supplement/Table_S5.docx", "word/document.xml", exdir = doc_dir)
text <- xml2::xml_text(xml2::read_xml(file.path(doc_dir, "word/document.xml")))
stopifnot(grepl("Statistical power at 0.75 x the meta-average", text, fixed = TRUE))
unlink(doc_dir, recursive = TRUE)
cat("PASS: 16 setups, all 52 comparison rows, paired cache intervals and current output inventory\n")
