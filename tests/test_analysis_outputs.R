## Run after main.R, from a directory containing its generated outputs.
source("scripts/functions.R")
args <- commandArgs(trailingOnly = TRUE)
expected_iterations <- if (length(args)) as.integer(args[[1]]) else 1000L
stopifnot(length(expected_iterations) == 1L, !is.na(expected_iterations),
          expected_iterations > 0L)
comparison <- read.csv("data/derived_data/counterfactual_comparison_0p5_vs_1.csv")
stopifnot(nrow(comparison) == 68L,
          setequal(unique(comparison$heterogeneity_multiplier), c(0, .25, .5, .75)))
for (h in c(0, .25, .5, .75)) {
  rows <- comparison[comparison$heterogeneity_multiplier == h, ]
  bins <- rows[1:13, ]
  stopifnot(identical(rows$measure, table_2_display_measures()),
            abs(sum(bins$reference_value) - 1) < 1e-10,
            abs(sum(bins$counterfactual_value) - 1) < 1e-10,
            max(abs(bins$estimate - (bins$reference_value - bins$counterfactual_value))) < 1e-12)
  draws <- lapply(c(.5, 1), function(m) {
    label <- paste0("meta_", gsub("\\.", "p", m), "_heterogeneity_", gsub("\\.", "p", h))
    path <- file.path("data", "derived_data", paste0("p_tab_ci_", label,
                       "_multilevel_random_outliers_removed.rds"))
    readRDS(path)$result[[1]]
  })
  stopifnot(nrow(draws[[1]]) == expected_iterations,
            nrow(draws[[2]]) == expected_iterations)
  difference_draws <- draws[[2]] - draws[[1]]
  expected_ci <- apply(difference_draws, 2, quantile, probs = c(.025, .975))
  ref_sig <- sum(bins$reference_value[11:13])
  target_sig <- sum(bins$counterfactual_value[11:13])
  bootstrap_all <- rowSums(difference_draws[, 11:13, drop = FALSE])
  bootstrap_sig <- bootstrap_all / rowSums(draws[[2]][, 11:13, drop = FALSE])
  stopifnot(max(abs(bins$ci_lower - expected_ci[1, ])) < 1e-12,
            max(abs(bins$ci_upper - expected_ci[2, ])) < 1e-12,
            abs(rows$estimate[14] - (ref_sig - target_sig)) < 1e-12,
            abs(rows$estimate[15] - (ref_sig - target_sig) / ref_sig) < 1e-12,
            max(abs(c(rows$ci_lower[14], rows$ci_upper[14]) - quantile(bootstrap_all, c(.025,.975)))) < 1e-12,
            max(abs(c(rows$ci_lower[15], rows$ci_upper[15]) - quantile(bootstrap_sig, c(.025,.975)))) < 1e-12,
            all(rows$estimate[16:17] == c(rows$n_meta_analyses[1], rows$n_effects[1])))
}
esr <- readRDS("data/derived_data/Figure_S5_inputs.rds")$all_esr_results
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
read_word_table <- function(path) {
  doc_dir <- tempfile("counterfactual-docx-")
  dir.create(doc_dir)
  on.exit(unlink(doc_dir, recursive = TRUE))
  unzip(path, "word/document.xml", exdir = doc_dir)
  doc <- xml2::read_xml(file.path(doc_dir, "word/document.xml"))
  rows <- xml2::xml_find_all(doc, ".//w:tbl[1]/w:tr")
  do.call(rbind, lapply(rows, function(row) xml2::xml_text(xml2::xml_find_all(row,"./w:tc"))))
}
main_table <- read_word_table("results/main/Table_2.docx")
supplement_table <- read_word_table("results/supplement/Table_S6.docx")
stopifnot(identical(dim(main_table), dim(supplement_table)),
          nrow(supplement_table) == 18L, ncol(supplement_table) == 5L,
          identical(main_table[,1], supplement_table[,1]),
          identical(main_table[1,], supplement_table[1,]))
for (number in 2:4) {
  path <- paste0("results/supplement/Table_S",number,".docx")
  doc_dir <- tempfile("power-docx-"); dir.create(doc_dir)
  unzip(path,"word/document.xml",exdir=doc_dir)
  text <- xml2::xml_text(xml2::read_xml(file.path(doc_dir,"word/document.xml")))
  multiplier <- c(.25,.75,1)[number-1]
  stopifnot(grepl(paste0("Statistical power at ",multiplier," x the meta-average"),text,fixed=TRUE))
  unlink(doc_dir,recursive=TRUE)
}
cat("PASS: Table S6 matches Table 2; paired bin and ESR intervals, output order and 16 setups verified\n")
