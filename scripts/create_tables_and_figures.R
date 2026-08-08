## ---------------------------------------------------------
## create_tables_and_figures.R
## ---------------------------------------------------------
## Recreate manuscript tables and figures in numeric order.
## Each section reinitializes its own inputs so it can be run independently
## (for example, by highlighting one section in RStudio) without relying on
## objects created by earlier table/figure sections.

library(tidyverse)
library(foreach)
library(doParallel)
library(openxlsx)
library(gridExtra)
library(ggeasy)
library(scales)
library(here)
library(officer)

source(here("scripts", "analysis_setup.R"))

required_setup_columns <- c("meta_average_multiplier", "heterogeneity_multiplier", "setup_label")
if (!all(required_setup_columns %in% names(analysis_setups))) {
  stop("analysis_setups must contain: ", paste(required_setup_columns, collapse = ", "))
}
if (anyDuplicated(analysis_setups$setup_label)) {
  stop("Each row of analysis_setups must have a unique setup_label.")
}

## Retain scalar defaults for helper calls made outside a setup-specific writer.
meta_average_multiplier <- analysis_setups$meta_average_multiplier[[1]]
heterogeneity_multiplier <- analysis_setups$heterogeneity_multiplier[[1]]
setup_label <- analysis_setups$setup_label[[1]]

## Heterogeneity affects only the counterfactual outputs (Figure 1 and Table 2).
## For all descriptive and power outputs, retain one row per meta-average value
## so identical files are not emitted once for every heterogeneity value.
heterogeneity_independent_setups <- analysis_setups %>%
  group_by(meta_average_multiplier) %>%
  slice(1) %>%
  ungroup()

safe_quantile <- function(x, probs, ...) {
  stats::quantile(x, probs = probs, na.rm = TRUE, ...)
}

ensure_output_dirs <- function() {
  invisible(lapply(
    list(
      here("results", "main"),
      here("results", "supplement"),
      here("results", "intermediate_results", "pet_peese_rstandard"),
      here("results", "intermediate_results", "multilevel_random"),
      here("data", "derived_data"),
      here("results", "robustness")
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

save_plot <- function(filename_stem, width, height, draw) {
  pdf(paste0(filename_stem, ".pdf"), width = width, height = height)
  draw()
  dev.off()

  cairo_ps(
    paste0(filename_stem, ".eps"),
    width = width,
    height = height,
    onefile = FALSE
  )
  draw()
  dev.off()

  svg(
    paste0(filename_stem, ".svg"),
    width = width,
    height = height
  )
  draw()
  dev.off()

  png(
    paste0(filename_stem, ".png"),
    width = width,
    height = height,
    units = "in",
    res = 300,
    type = "cairo"
  )
  draw()
  dev.off()
}

load_estimator_data <- function(estimator, setup_label_value = setup_label,
                                outlier_variant = "outliers_removed") {
  if (!outlier_variant %in% c("outliers_removed", "all_data")) {
    stop("Unknown outlier variant: ", outlier_variant)
  }

  variant_suffix <- if (outlier_variant == "all_data") "_all_data" else ""
  derived_path <- here(
    "data", "derived_data",
    paste0(
      "pps_rstandard_raw_", setup_label_value, "_", estimator,
      variant_suffix, ".rds"
    )
  )
  if (file.exists(derived_path)) {
    return(readRDS(derived_path))
  }

  estimator_dir <- switch(
    estimator,
    pet_peese = here(
      "results", "intermediate_results",
      if (outlier_variant == "all_data") "pet_peese_all_data" else "pet_peese_rstandard"
    ),
    multilevel_random = here(
      "results", "intermediate_results",
      if (outlier_variant == "all_data") "multilevel_random_all_data" else "multilevel_random"
    ),
    stop("Unknown estimator: ", estimator)
  )
  estimator_files <- list.files(
    estimator_dir,
    pattern = "\\.rds$",
    full.names = TRUE
  )

  if (length(estimator_files) == 0) {
    stop(
      "No ", estimator, " RDS files found in ", estimator_dir, ". ",
      "Run scripts/compute_meta_estimates.R first."
    )
  }

  estimator_files %>%
    map_dfr(readRDS) %>%
    mutate(
      sei = sqrt(vi),
      sse_yn = case_when(
        is.na(small_study_effect_pval) ~ NA_character_,
        small_study_effect_pval <= 0.05 ~ "yes",
        TRUE ~ "no"
      )
    )
}

add_power_variables <- function(dat, meta_average_multiplier = 0.5) {
  alpha <- 0.05
  q1 <- qnorm(1 - alpha / 2)
  q2 <- qnorm(alpha / 2)

  dat %>%
    mutate(
      GE = meta_average_multiplier * GE,
      lambda = abs(GE) / sei,
      power = 1 - pnorm(q1 - lambda) + pnorm(q2 - lambda),
      yn80 = ifelse(power >= 0.8, "yes", "no"),
      yn20 = ifelse(power <= 0.2, "yes", "no")
    )
}

load_power_data <- function(estimator, setup_label_value = setup_label, meta_average_multiplier_value = meta_average_multiplier) {
  derived_path <- here("data", "derived_data", paste0("pps_rstandard_power_", setup_label_value, "_", estimator, ".rds"))
  if (file.exists(derived_path)) {
    return(readRDS(derived_path))
  }
  load_estimator_data(estimator, setup_label_value) %>%
    add_power_variables(meta_average_multiplier_value)
}

split_meta_analyses <- function(dat, add_sape = FALSE) {
  split(dat, dat$cID) %>%
    lapply(function(x) {
      if (add_sape) {
        x$sape <- length(which(x$yn80 == "yes")) / length(x$vi)
      }
      x
    })
}

make_grids <- function() {
  list(
    p_grid_tab = c(
      -Inf,
      qnorm(0.001 / 2, lower.tail = TRUE), qnorm(0.01 / 2, lower.tail = TRUE),
      qnorm(0.05 / 2, lower.tail = TRUE), qnorm(0.1 / 2, lower.tail = TRUE),
      qnorm(0.2 / 2, lower.tail = TRUE), qnorm(0.3 / 2, lower.tail = TRUE),
      qnorm(0.4 / 2, lower.tail = TRUE), qnorm(0.5 / 2, lower.tail = TRUE),
      qnorm(0.6 / 2, lower.tail = TRUE), qnorm(0.7 / 2, lower.tail = TRUE),
      qnorm(0.8 / 2, lower.tail = TRUE), qnorm(0.9 / 2, lower.tail = TRUE),
      0,
      qnorm(0.9 / 2, lower.tail = FALSE), qnorm(0.8 / 2, lower.tail = FALSE),
      qnorm(0.7 / 2, lower.tail = FALSE), qnorm(0.6 / 2, lower.tail = FALSE),
      qnorm(0.5 / 2, lower.tail = FALSE), qnorm(0.4 / 2, lower.tail = FALSE),
      qnorm(0.3 / 2, lower.tail = FALSE), qnorm(0.2 / 2, lower.tail = FALSE),
      qnorm(0.1 / 2, lower.tail = FALSE), qnorm(0.05 / 2, lower.tail = FALSE),
      qnorm(0.01 / 2, lower.tail = FALSE), qnorm(0.001 / 2, lower.tail = FALSE),
      Inf
    ),
    p_grid_tab2 = c(
      0,
      qnorm(0.9 / 2, lower.tail = FALSE), qnorm(0.8 / 2, lower.tail = FALSE),
      qnorm(0.7 / 2, lower.tail = FALSE), qnorm(0.6 / 2, lower.tail = FALSE),
      qnorm(0.5 / 2, lower.tail = FALSE), qnorm(0.4 / 2, lower.tail = FALSE),
      qnorm(0.3 / 2, lower.tail = FALSE), qnorm(0.2 / 2, lower.tail = FALSE),
      qnorm(0.1 / 2, lower.tail = FALSE), qnorm(0.05 / 2, lower.tail = FALSE),
      qnorm(0.01 / 2, lower.tail = FALSE), qnorm(0.001 / 2, lower.tail = FALSE),
      Inf
    ),
    p_grid_plot = qnorm(seq(1, 0, -0.005), lower.tail = FALSE),
    z_grid_plot = seq(-10.25, 10.25, 0.1025),
    z_grid_plot2 = seq(0, 10.25, 0.1025)
  )
}

count_intervals <- function(values, grid) {
  vapply(seq_len(length(grid) - 1), function(a) {
    length(which(values >= grid[a] & values <= grid[a + 1]))
  }, numeric(1))
}

cf_ci_with_progress <- function(dat, grid, cluster, heterogeneity_multiplier, label) {
  components <- counterfactual_components(
    dat, grid, heterogeneity_multiplier
  )
  iteration_batches <- split(
    seq_len(n_iterations),
    ceiling(seq_len(n_iterations) / max(1, n_cores))
  )
  progress_bar <- new_progress_bar(
    n_iterations,
    paste0(label, " (bootstrap replications)")
  )
  on.exit(close_progress_bar(progress_bar), add = TRUE)

  batch_results <- lapply(iteration_batches, function(iteration_ids) {
    message(
      sprintf(
        "Running bootstrap replications %d-%d of %d",
        min(iteration_ids), max(iteration_ids), n_iterations
      )
    )
    result <- cf.ci.cluster(
      dat = dat,
      z.grid = grid,
      iters = length(iteration_ids),
      cluster = cluster,
      heterogeneity_multiplier = heterogeneity_multiplier,
      iteration_ids = iteration_ids,
      components = components
    )
    update_progress_bar(progress_bar, max(iteration_ids))
    result
  })

  lapply(seq_len(3), function(result_index) {
    do.call(rbind, lapply(batch_results, `[[`, result_index))
  })
}

get_counterfactual <- function(path, dat, grid, ci = FALSE, cluster = NULL, heterogeneity_multiplier_value = heterogeneity_multiplier) {
  cache_key <- counterfactual_cache_key(
    dat, grid, heterogeneity_multiplier_value, ci, cluster,
    if (ci) n_iterations else NULL
  )
  cached_result <- read_counterfactual_cache(path, cache_key)
  if (!is.null(cached_result)) {
    message("Using cached counterfactual: ", basename(path))
    return(cached_result)
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  result <- if (ci) {
    cf_ci_with_progress(
      dat = dat,
      grid = grid,
      cluster = cluster,
      heterogeneity_multiplier = heterogeneity_multiplier_value,
      label = paste0("Computing ", basename(path))
    )
  } else {
    message("Computing counterfactual: ", basename(path))
    cf(dat = dat, z.grid = grid, heterogeneity_multiplier = heterogeneity_multiplier_value)
  }
  write_counterfactual_cache(result, path, cache_key)
  result
}

ensure_output_dirs()

## -------------------------------
## Table 1
## -------------------------------
write_table_1 <- function() {
  ## Table 1 describes the sample rather than an outlier-screened estimate of
  ## a meta-average. Use the all-data random-effects estimates to calculate its
  ## power classification so that every meta-analysis and primary
  ## estimate is represented. Outlier-screened data are reserved for outputs
  ## whose estimands depend on the meta-average.
  pps_rstandard <- load_estimator_data(
    "multilevel_random", "meta_0p5_heterogeneity_0", "all_data"
  ) %>%
    add_power_variables(0.5)
myDat <- split_meta_analyses(pps_rstandard, add_sape = TRUE)
mss <- vapply(myDat, function(x) length(x$sei), numeric(1))

idx <- tibble(
  i = seq_along(myDat),
  sdesn = vapply(myDat, function(x) x$sdesn[1], character(1)),
  sape = vapply(myDat, function(x) x$sape[1], numeric(1)),
  guide = vapply(myDat, function(x) x$guide[1], character(1)),
  prere = vapply(myDat, function(x) x$prere[1], character(1))
)

fill_desc <- function(i) c(
  length(i), sum(mss[i]), round(mean(mss[i])), min(mss[i]),
  round(quantile(mss[i], 0.25)), round(median(mss[i])),
  round(quantile(mss[i], 0.75)), max(mss[i])
)

table_rows <- list(
  c("All meta-analyses", fill_desc(seq_along(myDat))),
  c("Research design", rep("", 8)),
  c("  Observational", fill_desc(idx$i[idx$sdesn %in% c("observational", "mixed")])),
  c("  Experimental", fill_desc(idx$i[idx$sdesn == "experimental"])),
  c("Statistical power", rep("", 8)),
  c("  SAPE > 0", fill_desc(idx$i[idx$sape > 0])),
  c("  SAPE = 0", fill_desc(idx$i[idx$sape == 0])),
  c("Followed guidelines", rep("", 8)),
  c("  Yes", fill_desc(idx$i[idx$guide == "yes"])),
  c("  No", fill_desc(idx$i[idx$guide == "no"])),
  c("Protocol registered", rep("", 8)),
  c("  Yes", fill_desc(idx$i[idx$prere == "yes"])),
  c("  No", fill_desc(idx$i[idx$prere == "no"]))
)
table_data <- as.data.frame(do.call(rbind, table_rows), stringsAsFactors = FALSE)
names(table_data) <- c("Group", "Meta-analyses", "Primary estimates",
                       "Mean", "Min", "Q25", "Q50", "Q75", "Max")
table_data[-1] <- lapply(table_data[-1], function(x) suppressWarnings(as.numeric(x)))
section_rows <- c(2, 5, 8, 11)

workbook <- createWorkbook()
addWorksheet(workbook, "Table 1", gridLines = FALSE)
writeData(workbook, "Table 1", "No. of", startCol = 2, startRow = 1)
writeData(workbook, "Table 1", "No. of", startCol = 3, startRow = 1)
writeData(workbook, "Table 1", "Primary estimates per meta-analysis", startCol = 4, startRow = 1)
writeData(workbook, "Table 1", names(table_data)[1:3], startRow = 2, colNames = FALSE)
writeData(workbook, "Table 1", names(table_data)[4:9], startCol = 4, startRow = 2, colNames = FALSE)
writeData(workbook, "Table 1", table_data, startRow = 3, colNames = FALSE)
mergeCells(workbook, "Table 1", cols = 4:9, rows = 1)
header_style <- createStyle(fontSize = 11, textDecoration = "bold", halign = "center",
                            valign = "center", border = "bottom", borderStyle = "medium")
section_style <- createStyle(textDecoration = "bold", fgFill = "#EAF0F6",
                             border = "bottom", borderColour = "#B7C3D0")
body_style <- createStyle(fontSize = 10, valign = "center")
addStyle(workbook, "Table 1", header_style, rows = 1:2, cols = 1:9, gridExpand = TRUE)
addStyle(workbook, "Table 1", body_style, rows = 3:(nrow(table_data) + 2), cols = 1:9, gridExpand = TRUE)
addStyle(workbook, "Table 1", section_style, rows = section_rows + 2, cols = 1:9, gridExpand = TRUE)
setColWidths(workbook, "Table 1", cols = 1, widths = 25)
setColWidths(workbook, "Table 1", cols = 2:3, widths = 17)
setColWidths(workbook, "Table 1", cols = 4:9, widths = 10)
freezePane(workbook, "Table 1", firstActiveRow = 3)
saveWorkbook(workbook, here("results", "main", "Table_1.xlsx"), overwrite = TRUE)

docx_table <- table_data
docx_table[section_rows, -1] <- ""
document <- officer::read_docx()
document <- officer::body_add_par(document, "Table 1. Characteristics of the meta-analyses",
                                  style = "heading 1")
## Do not request a template-specific Word style here. Minimal or customized
## reference documents do not necessarily define the built-in "Table Grid"
## style, and officer rejects a style that is absent from the document.
document <- officer::body_add_table(document, docx_table, style = NULL)
print(document, target = here("results", "main", "Table_1.docx"))
}

write_table_1()

## -------------------------------
## Figure 1
## -------------------------------
make_figure_1_panel <- function(meta_average_multiplier, heterogeneity_multiplier, setup_label, estimator) {
pps_rstandard <- load_estimator_data(estimator, setup_label)
grids <- make_grids()
my_dat <- pps_rstandard %>%
  mutate(GE = meta_average_multiplier * GE) %>%
  filter_counterfactual_data(
    heterogeneity_multiplier,
    context = paste("Figure 1", setup_label, estimator)
  )
myDat <- split_meta_analyses(my_dat)
facz <- abs(my_dat$yi / sqrt(my_dat$vi))
z.orig <- count_intervals(facz, grids$z_grid_plot2)
p.orig.plot <- count_intervals(facz, grids$p_grid_plot[which(grids$p_grid_plot >= 0)])
z.plot <- get_counterfactual(here("data", "derived_data", paste0("z_plot_", setup_label, "_", estimator, ".rds")), myDat, grids$z_grid_plot, heterogeneity_multiplier_value = heterogeneity_multiplier)
z.plot.ci <- get_counterfactual(here("data", "derived_data", paste0("z_plot_ci_", setup_label, "_", estimator, ".rds")), myDat, grids$z_grid_plot, ci = TRUE, cluster = unique(my_dat$cID), heterogeneity_multiplier_value = heterogeneity_multiplier)

xs <- as.vector(grids$z_grid_plot2[-length(grids$z_grid_plot2)] + (grids$z_grid_plot2[2] - grids$z_grid_plot2[1]) / 2)
## cf.ci.cluster() normalizes every bootstrap draw by the number of effects in
## that draw. Use the corresponding full-sample denominator for the point
## estimate; sum(p.orig.plot) can be smaller when an observed |z| falls beyond
## the finite plotting grid, shifting the orange curve above its bootstrap CI.
N <- nrow(my_dat)
datFull <- as.data.frame(cbind(
  xs = xs,
  q025 = as.vector(apply(z.plot.ci[[1]], 2, quantile, na.rm = TRUE, probs = c(0.025))),
  q975 = as.vector(apply(z.plot.ci[[1]], 2, quantile, na.rm = TRUE, probs = c(0.975))),
  n.f = as.vector(z.orig / N),
  n.cf = as.vector(z.plot / N)
))
ggplot(datFull) +
  geom_line(aes(xs, q025), color = "orange", lty = 3) +
  geom_line(aes(xs, n.cf), color = "orange", lty = 1) +
  geom_point(aes(xs, n.cf), shape = 20, fill = "orange", color = "orange", size = 1) +
  geom_line(aes(xs, q975), color = "orange", lty = 3) +
  geom_line(aes(xs, n.f), color = "blue", lty = 2) +
  geom_point(aes(xs, n.f), shape = 20, fill = "blue", color = "blue", size = 1) +
  coord_cartesian(xlim = c(0, 8)) +
  xlab("|z|-value") + ylab("Frequency") +
  ggtitle(paste0(heterogeneity_multiplier * 100, "% genuine heterogeneity")) +
  geom_vline(xintercept = c(1.64, 1.96, 2.58), lty = 2, color = c(3, 2, 6), linewidth = 0.5) +
  scale_x_continuous(
    breaks = c(0, 1.64, 1.96, 2.58, 4, 6, 8),
    guide = guide_axis(n.dodge = 2)
  ) +
  theme(panel.background = element_rect(fill = "gray100"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(linewidth = 0.5, color = "gray"))
}

figure_1_setups <- analysis_setups %>%
  filter(meta_average_multiplier == 0.5, heterogeneity_multiplier %in% c(0, 0.5)) %>%
  arrange(heterogeneity_multiplier)
figure_1_panels <- purrr::pmap(
  figure_1_setups,
  function(meta_average_multiplier, heterogeneity_multiplier, setup_label, ...) {
    make_figure_1_panel(meta_average_multiplier, heterogeneity_multiplier, setup_label, "multilevel_random")
  }
)
figure_1 <- arrangeGrob(grobs = figure_1_panels, ncol = 1)
save_plot(
  here("results", "main", "Figure_1_multilevel_random"),
  width = 10,
  height = 10,
  draw = function() grid::grid.draw(figure_1)
)

## -------------------------------
## Table 2
## -------------------------------
calculate_table_2 <- function(meta_average_multiplier, heterogeneity_multiplier,
                              setup_label, estimator, outlier_variant, ...) {
  pps_rstandard <- load_estimator_data(
    estimator, setup_label, outlier_variant
  )
  grids <- make_grids()
  my_dat <- pps_rstandard %>%
    mutate(GE = meta_average_multiplier * GE) %>%
    filter_counterfactual_data(
      heterogeneity_multiplier,
      context = paste("Table 2", setup_label, estimator)
    )
  myDat <- split_meta_analyses(my_dat)
  facz <- abs(my_dat$yi / sqrt(my_dat$vi))
  p.orig.tab <- count_intervals(facz, grids$p_grid_tab2)
  result_suffix <- paste(setup_label, estimator, outlier_variant, sep = "_")
  p.tab <- get_counterfactual(here("data", "derived_data", paste0("p_tab_", result_suffix, ".rds")), myDat, grids$p_grid_tab, heterogeneity_multiplier_value = heterogeneity_multiplier)
  p.tab.ci <- get_counterfactual(here("data", "derived_data", paste0("p_tab_ci_", result_suffix, ".rds")), myDat, grids$p_grid_tab, ci = TRUE, cluster = unique(my_dat$cID), heterogeneity_multiplier_value = heterogeneity_multiplier)

  include_p_value_intervals <- estimator == "multilevel_random" && meta_average_multiplier == 0.5
  p.table <- matrix(NA_character_, ncol = 2, nrow = length(grids$p_grid_tab2) - 1)
  colnames(p.table) <- c("estimate", "confidence_interval")
  N <- sum(p.orig.tab)
  if (include_p_value_intervals) {
    p.table[, 1] <- round((p.orig.tab - p.tab) / N, 3)
    q025 <- apply(matrix(p.orig.tab / N, nrow = nrow(p.tab.ci[[1]]), ncol = ncol(p.tab.ci[[1]]), byrow = TRUE) - p.tab.ci[[1]], 2, safe_quantile, probs = c(0.025))
    q975 <- apply(matrix(p.orig.tab / N, nrow = nrow(p.tab.ci[[1]]), ncol = ncol(p.tab.ci[[1]]), byrow = TRUE) - p.tab.ci[[1]], 2, safe_quantile, probs = c(0.975))
    p.table[, 2] <- paste("[", round(q025, 3), ", ", round(q975, 3), "]", sep = "")
  }
  for (level in list(c(10, 13, 1, "all"), c(11, 13, 1, "all"), c(10, 13, 2, "sig"), c(11, 13, 3, "sig"))) {
    lo <- as.integer(level[[1]]); hi <- as.integer(level[[2]]); ci_idx <- as.integer(level[[3]]); denom <- if (level[[4]] == "all") N else sum(p.orig.tab[lo:hi])
    point <- round(sum((p.orig.tab - p.tab)[lo:hi] / denom, na.rm = TRUE), 3)
    bs <- apply(p.tab.ci[[ci_idx]][, lo:hi], 1, sum)
    q <- round(safe_quantile(sum((p.orig.tab / denom)[lo:hi]) - bs, probs = c(0.025, 0.975)), 3)
    p.table <- rbind(p.table, c(point, paste("[", q[1], ", ", q[2], "]", sep = "")))
  }
  p.table <- rbind(p.table, c(length(myDat), 0), c(N, 0))
  rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7", "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3", "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01", "p < 0.001", "ESR_{0.1}^{all}", "ESR_{0.05}^{all}", "ESR_{0.1}^{sig}", "ESR_{0.05}^{sig}", "No. of meta-analysis", "No. of tests")
  summary_rows <- c("ESR_{0.05}^{all}", "ESR_{0.05}^{sig}", "No. of meta-analysis", "No. of tests")
  list(
    detailed = as.data.frame(p.table) %>% rownames_to_column("measure"),
    summary = as.data.frame(p.table) %>% rownames_to_column("measure") %>%
      filter(measure %in% summary_rows) %>%
      mutate(estimator = estimator, meta_average_multiplier = meta_average_multiplier,
             heterogeneity_multiplier = heterogeneity_multiplier, setup_label = setup_label,
             outlier_variant = outlier_variant,
             .before = 1)
  )
}

table_2_parameters <- tidyr::crossing(
  analysis_setups,
  estimator = meta_analysis_estimators,
  outlier_variant = "outliers_removed"
)
table_2_results <- table_2_parameters %>%
  pmap(calculate_table_2)

all_esr_results <- map_dfr(table_2_results, "summary")
all_combination_results <- map2_dfr(
  table_2_results,
  seq_len(nrow(table_2_parameters)),
  function(result, i) {
    parameters <- table_2_parameters[i, ]
    rows <- if (parameters$estimator == "multilevel_random" && parameters$meta_average_multiplier == 0.5) {
      result$detailed
    } else {
      result$summary
    }
    rows %>% mutate(
      estimator = parameters$estimator,
      meta_average_multiplier = parameters$meta_average_multiplier,
      heterogeneity_multiplier = parameters$heterogeneity_multiplier,
      setup_label = parameters$setup_label,
      outlier_variant = parameters$outlier_variant,
      .before = 1
    ) %>% dplyr::select(-any_of(c("estimator1", "meta_average_multiplier1", "heterogeneity_multiplier1", "setup_label1", "outlier_variant1")))
  }
)
figure_s3_inputs_path <- here(
  "data", "derived_data", "Figure_S3_inputs.rds"
)
saveRDS(
  list(
    all_combination_results = all_combination_results,
    all_esr_results = all_esr_results
  ),
  figure_s3_inputs_path
)
table_2_indices <- table_2_parameters %>%
  mutate(result_index = row_number()) %>%
  filter(estimator == "multilevel_random", meta_average_multiplier == 0.5) %>%
  arrange(heterogeneity_multiplier)
table_2_columns <- map2(
  table_2_indices$result_index,
  table_2_indices$heterogeneity_multiplier,
  function(i, h) table_2_results[[i]]$detailed %>%
    transmute(measure, !!paste0("heterogeneity_", h) := if_else(
      confidence_interval == "0", estimate, paste(estimate, confidence_interval)
    ))
)
table_2 <- reduce(table_2_columns, full_join, by = "measure") %>%
  filter(!measure %in% c("ESR_{0.1}^{all}", "ESR_{0.1}^{sig}"))
names(table_2) <- c(
  "p-value interval",
  "(1)\nHalf the meta-average\nDifference [95% CI]",
  "(2)\nHalf the meta-average and 25% genuine heterogeneity\nDifference [95% CI]",
  "(3)\nHalf the meta-average and 50% genuine heterogeneity\nDifference [95% CI]",
  "(4)\nHalf the meta-average and 75% genuine heterogeneity\nDifference [95% CI]"
)
table_2_document <- officer::read_docx()
table_2_document <- officer::body_add_table(
  table_2_document, table_2, style = NULL, header = TRUE,
  alignment = c("left", rep("center", 4)), align_table = "center"
)
print(table_2_document, target = here("results", "main", "Table_2.docx"))

## -------------------------------
## Table 3
## -------------------------------
write_table_3 <- function(meta_average_multiplier, heterogeneity_multiplier, setup_label, estimator, ...) {
pps_rstandard <- load_estimator_data(estimator, setup_label) %>% add_power_variables(meta_average_multiplier)
subf_desc <- pps_rstandard %>% group_by(subfd) %>% summarise(M = length(unique(cID)), N = length(power), median = round(median(power, na.rm = TRUE), 2), mean = round(mean(power, na.rm = TRUE), 2), Q25 = round(safe_quantile(power, 0.25), 2), Q75 = round(safe_quantile(power, 0.75), 2), sape = round(sum(power >= 0.8, na.rm = TRUE) / sum(!is.na(power)), 2), .groups = "drop")
med_med <- pps_rstandard %>% group_by(cID) %>% summarise(metaID = metaID[1], subfd = subfd[1], median = median(power, na.rm = TRUE), .groups = "drop")
med_med_subf <- med_med %>% group_by(subfd) %>% summarise(mmedian = round(median(median), 2), .groups = "drop")
power.tab3 <- matrix(nrow = 8, ncol = 8)
power.tab3[1, ] <- c(nrow(med_med), nrow(pps_rstandard), round(summary(med_med$median), 2)[3], round(summary(pps_rstandard$power), 2)[3], round(summary(pps_rstandard$power), 2)[4], round(summary(pps_rstandard$power), 2)[2], round(summary(pps_rstandard$power), 2)[5], 0.18)
for (i in seq_len(nrow(subf_desc))) power.tab3[i + 1, ] <- c(subf_desc$M[i], subf_desc$N[i], med_med_subf$mmedian[i], subf_desc$median[i], subf_desc$mean[i], subf_desc$Q25[i], subf_desc$Q75[i], subf_desc$sape[i])
colnames(power.tab3) <- c("No. of meta-\nanalyses", "No. of primary\nestimates", "Median of\nmedians", "Median", "Mean", "Q25", "Q75", "SAPE")
rownames(power.tab3) <- c("All meta-analyses", "Ecology", "Environmental Chemistry", "Environmental Engineering", "Health, Toxicology and Mutagenesis", "Management, Monitoring, Policy and Law", "Nature and Landscape Conservation", "Water Science and Technology")
power_table <- as.data.frame(power.tab3) %>% rownames_to_column("Subfield")
document <- officer::read_docx()
document <- officer::body_add_table(
  document, power_table, style = NULL, header = TRUE,
  alignment = c("left", rep("center", 8)), align_table = "center"
)
print(document, target = here("results", "main", "Table_3.docx"))
}

analysis_setups %>%
  filter(setup_label == "meta_0p5_heterogeneity_0") %>%
  mutate(estimator = "multilevel_random") %>%
  pwalk(write_table_3)

## -------------------------------
## Table 4
## -------------------------------
## run_exploratory_regressions.R creates these model and robust-inference
## objects. Keeping presentation here makes the numbering and output conventions
## consistent with the other manuscript tables.
required_nb_objects <- c(
  "nbMod1", "nbMod2", "nbMod3", "nbMod4",
  "nbMod1.robu", "nbMod2.robu", "nbMod3.robu", "nbMod4.robu"
)
if (!all(vapply(required_nb_objects, exists, logical(1), inherits = TRUE))) {
  stop("Run scripts/run_exploratory_regressions.R before creating Table 4.")
}

significance_stars <- function(p_value) {
  ifelse(p_value < 0.01, "***", ifelse(p_value < 0.05, "**",
    ifelse(p_value < 0.1, "*", "")))
}

format_nb_term <- function(robust_result, term) {
  if (is.na(term) || !term %in% rownames(robust_result)) return("")
  sprintf(
    "%.3f%s (%.3f)", robust_result[term, 1],
    significance_stars(robust_result[term, 4]), robust_result[term, 2]
  )
}

table_4_terms <- tibble::tribble(
  ~Variable, ~term_model_1, ~term_model_2,
  "Intercept", "(Intercept)", "(Intercept)",
  "Median power", "med_perc", "med_perc",
  "Experimental research design? (yes)", "design_mergedyes", "design_mergedyes",
  "Followed reporting guidelines? (yes)", "guidyes", "guidyes",
  "Protocol registered? (yes)", "preryes", "preryes",
  "Log number of independent studies", "lognps", "lognps",
  "Log journal impact factor", "logjif", "logjif",
  "Publication year", "pyear", "pyear",
  "Environmental Chemistry", NA_character_, "subfEnvironmental Chemistry",
  "Environmental Engineering", NA_character_, "subfEnvironmental Engineering",
  "Health, Toxicology and Mutagenesis", NA_character_, "subfHealth, Toxicology and Mutagenesis",
  "Management, Monitoring, Policy and Law", NA_character_, "subfManagement, Monitoring, Policy and Law",
  "Nature and Landscape Conservation", NA_character_, "subfNature and Landscape Conservation",
  "Water Science and Technology", NA_character_, "subfWater Science and Technology"
) %>%
  dplyr::mutate(
    `0% heterogeneity: Model 1 Estimate (SE)` = vapply(term_model_1, format_nb_term, character(1),
      robust_result = nbMod1.robu),
    `0% heterogeneity: Model 2 Estimate (SE)` = vapply(term_model_2, format_nb_term, character(1),
      robust_result = nbMod2.robu),
    `50% heterogeneity: Model 1 Estimate (SE)` = vapply(term_model_1, format_nb_term, character(1),
      robust_result = nbMod3.robu),
    `50% heterogeneity: Model 2 Estimate (SE)` = vapply(term_model_2, format_nb_term, character(1),
      robust_result = nbMod4.robu)
  ) %>%
  dplyr::select(-term_model_1, -term_model_2) %>%
  dplyr::bind_rows(tibble::tibble(
    Variable = c("Effect size type", "AIC", "No. of meta-analyses"),
    `0% heterogeneity: Model 1 Estimate (SE)` = c("Yes", sprintf("%.1f", AIC(nbMod1)), nobs(nbMod1)),
    `0% heterogeneity: Model 2 Estimate (SE)` = c("Yes", sprintf("%.1f", AIC(nbMod2)), nobs(nbMod2)),
    `50% heterogeneity: Model 1 Estimate (SE)` = c("Yes", sprintf("%.1f", AIC(nbMod3)), nobs(nbMod3)),
    `50% heterogeneity: Model 2 Estimate (SE)` = c("Yes", sprintf("%.1f", AIC(nbMod4)), nobs(nbMod4))
  ))

table_4_document <- officer::read_docx()
table_4_document <- officer::body_add_par(
  table_4_document,
  "Table 4. Negative binomial regression of excess significant results",
  style = "heading 1"
)
table_4_document <- officer::body_add_table(
  table_4_document, table_4_terms, style = NULL, header = TRUE,
  alignment = c("left", rep("center", 4)), align_table = "center"
)
table_4_document <- officer::body_add_par(
  table_4_document,
  "Note. The meta-average multiplier is 0.5. Cluster-robust standard errors in parentheses. * p < .10; ** p < .05; *** p < .01. Ecology is the reference subfield.",
  style = NULL
)
print(table_4_document, target = here("results", "main", "Table_4.docx"))

## -------------------------------
## Figure 2
## -------------------------------
write_figure_2 <- function(meta_average_multiplier, heterogeneity_multiplier, setup_label, estimator, ...) {
pps_rstandard <- load_power_data(estimator, setup_label, meta_average_multiplier)
pps_rstandard_median <- pps_rstandard %>% group_by(cID) %>% summarise(metaID = metaID[1], median = median(power, na.rm = TRUE), sape = sum(power >= 0.8, na.rm = TRUE) / sum(!is.na(power)), nips = length(unique(sID)), esty = unique(etype), guid = unique(guide), prer = unique(prere), subf = unique(subfd), sdes = unique(sdesn), .groups = "drop") %>% mutate(yn80 = ifelse(median >= 0.8, "yes", "no"), median100 = round(100 * median, 2), sape100 = round(100 * sape, 2))
write.xlsx(pps_rstandard_median, here("data", "derived_data", paste0("Figure_2_data_", setup_label, "_", estimator, ".xlsx")), overwrite = TRUE)

figure_2_summary <- tibble(
  result = c(
    "Median power of 20% or less",
    "Median power greater than 80%",
    "No adequately powered estimates",
    "Share of adequately powered estimates greater than 20%"
  ),
  count = c(
    sum(pps_rstandard_median$median <= 0.2, na.rm = TRUE),
    sum(pps_rstandard_median$median > 0.8, na.rm = TRUE),
    sum(pps_rstandard_median$sape == 0, na.rm = TRUE),
    sum(pps_rstandard_median$sape > 0.2, na.rm = TRUE)
  ),
  total_meta_analyses = nrow(pps_rstandard_median)
) %>%
  mutate(percentage = round(100 * count / total_meta_analyses, 1))

high_power_subfields <- pps_rstandard_median %>%
  filter(median > 0.8) %>%
  count(subf, name = "count") %>%
  mutate(
    total_high_power_meta_analyses = sum(count),
    percentage = round(100 * count / total_high_power_meta_analyses, 1)
  ) %>%
  arrange(desc(count)) %>%
  rename(subfield = subf)

write.xlsx(
  list(
    `Figure 2 summary` = figure_2_summary,
    `High-power subfields` = high_power_subfields
  ),
  here("results", "main", "Figure_2_summary.xlsx"),
  overwrite = TRUE
)
med_pwr <- pps_rstandard_median %>% ggplot(aes(x = median100, fill = as.factor(yn80))) + geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins = 30, alpha = I(0.6), linewidth = 0.1) + scale_fill_manual(values = c("brown2", "skyblue2")) + xlab("Median statistical power of primary estimates per meta-analysis") + ylab("Percentage") + ggtitle("(a)") + scale_x_continuous(breaks = breaks_width(20), labels = label_percent(scale = 1), expand = c(0, 0.5)) + scale_y_continuous(labels = label_percent(scale = 1), expand = c(0, 0.5)) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.5, color = "gray"))
sape <- pps_rstandard_median %>% ggplot(aes(x = sape100)) + geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins = 30, alpha = I(0.6), linewidth = 0.1, fill = "skyblue2") + xlab("Share of adequately powered primary estimates per meta-analysis") + ylab("Percentage") + ggtitle("(b)") + scale_x_continuous(breaks = breaks_width(20), labels = label_percent(scale = 1), expand = c(0, 0.5)) + scale_y_continuous(labels = label_percent(scale = 1), expand = c(0, 0.5)) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.5, color = "gray"))
figure_2 <- arrangeGrob(med_pwr, sape, ncol = 2)
save_plot(
  here("results", "main", "Figure_2"),
  width = 10,
  height = 4,
  draw = function() grid::grid.draw(figure_2)
)
}

analysis_setups %>%
  filter(
    setup_label == "meta_0p5_heterogeneity_0"
  ) %>%
  mutate(estimator = "multilevel_random") %>%
  pwalk(write_figure_2)
