## ---------------------------------------------------------
## create_supplement_tables_and_figures.R
## ---------------------------------------------------------
## Recreate all supplementary tables and figures from the analysis outputs. The
## fitted negative-binomial models used for Figure S4 are produced by
## run_exploratory_regressions.R.

library(tidyverse)
library(foreach)
library(doParallel)
library(gridExtra)
library(here)
library(openxlsx)

source(here("scripts", "analysis_setup.R"))

## Keep the supplementary renderer runnable in a fresh session after its cached
## analysis inputs have been created by create_tables_and_figures.R.
if (!exists("write_latex_table", mode = "function")) {
  latex_escape <- function(x) {
    replacements <- c(
      "\\" = "\\textbackslash{}", "&" = "\\&", "%" = "\\%", "#" = "\\#",
      "_" = "\\_", "$" = "\\$", "{" = "\\{", "}" = "\\}"
    )
    vapply(as.character(x), function(value) {
      characters <- strsplit(gsub("\n", " ", value, fixed = TRUE), "", fixed = TRUE)[[1]]
      paste0(ifelse(characters %in% names(replacements), replacements[characters], characters),
             collapse = "")
    }, character(1), USE.NAMES = FALSE)
  }
  write_latex_table <- function(x, path, alignment = NULL) {
    if (is.null(alignment)) alignment <- paste0("l", strrep("c", ncol(x) - 1))
    rows <- apply(x, 1, function(row) {
      paste0(paste(latex_escape(row), collapse = " & "), " \\\\")
    })
    writeLines(c(
      paste0("\\begin{tabular}{", alignment, "}"), "\\hline",
      paste0(paste(latex_escape(names(x)), collapse = " & "), " \\\\"),
      "\\hline", rows, "\\hline", "\\end{tabular}"
    ), path, useBytes = TRUE)
  }
}

supplement_dir <- here("results", "supplement")
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
derived_data_dir <- here("data", "derived_data")
dir.create(derived_data_dir, recursive = TRUE, showWarnings = FALSE)

save_supplement_plot <- function(filename_stem, width, height, draw) {
  save_with_device <- function(extension, open_device) {
    output_path <- paste0(filename_stem, extension)
    temporary_path <- tempfile(
      pattern = paste0(basename(filename_stem), "_"),
      tmpdir = dirname(filename_stem),
      fileext = extension
    )
    on.exit(unlink(temporary_path), add = TRUE)

    open_device(temporary_path)
    device <- dev.cur()
    on.exit({
      if (device %in% dev.list()) dev.off(device)
    }, add = TRUE)
    draw()
    dev.off(device)

    if (!file.rename(temporary_path, output_path)) {
      stop("Could not move completed plot to ", output_path)
    }
  }

  ## Write each format to a temporary file first. If drawing fails, this keeps
  ## a truncated device output from masquerading as a valid PDF (or image).
  save_with_device(".pdf", function(path) pdf(path, width = width, height = height))
  save_with_device(".eps", function(path) cairo_ps(
    path, width = width, height = height, onefile = FALSE
  ))
  save_with_device(".svg", function(path) svg(path, width = width, height = height))
  save_with_device(".png", function(path) png(
    path, width = width, height = height, units = "in", res = 300,
    type = "cairo"
  ))
}

load_multilevel_data <- function(setup_label) {
  derived_path <- here(
    "data", "derived_data",
    paste0("pps_rstandard_raw_", setup_label, "_multilevel_random.rds")
  )
  if (file.exists(derived_path)) {
    return(readRDS(derived_path))
  }

  input_files <- list.files(
    here("data", "derived_data", "multilevel_random"),
    pattern = "\\.rds$",
    full.names = TRUE
  )
  if (length(input_files) == 0) {
    stop("No multilevel random-effects inputs found. Run scripts/compute_meta_estimates.R first.")
  }
  purrr::map_dfr(input_files, readRDS)
}

load_multilevel_all_data <- function() {
  input_files <- list.files(
    here("data", "derived_data", "multilevel_random_all_data"),
    pattern = "\\.rds$", full.names = TRUE
  )
  if (length(input_files) == 0) {
    stop("No all-data random-effects inputs found. Run scripts/compute_meta_estimates.R first.")
  }
  purrr::map_dfr(input_files, readRDS)
}

## Small-study effects are tested by the standard-error slope in the multilevel
## random-effects Egger regression. Use the random-effects outlier-removed data,
## matching the primary estimator and its analysis sample.
load_small_study_effects <- function(setup_label) {
  derived_path <- here(
    "data", "derived_data",
    paste0("pps_rstandard_raw_", setup_label, "_multilevel_random.rds")
  )
  egger_data <- if (file.exists(derived_path)) {
    readRDS(derived_path)
  } else {
    input_files <- list.files(
      here("data", "derived_data", "multilevel_random"),
      pattern = "\\.rds$", full.names = TRUE
    )
    if (length(input_files) == 0) {
      stop("No random-effects Egger inputs found for the small-study-effect tests. ",
           "Run scripts/compute_meta_estimates.R first.")
    }
    purrr::map_dfr(input_files, readRDS)
  }

  egger_data %>%
    group_by(cID) %>%
    summarise(
      small_study_effect_pval = first(small_study_effect_pval),
      .groups = "drop"
    )
}

count_intervals <- function(values, grid) {
  vapply(seq_len(length(grid) - 1), function(i) {
    sum(values >= grid[i] & values <= grid[i + 1])
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

get_counterfactual <- function(path, dat, grid, ci = FALSE, cluster = NULL,
                               heterogeneity_multiplier) {
  cache_key <- counterfactual_cache_key(
    dat, grid, heterogeneity_multiplier, ci, cluster,
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
      dat = dat, grid = grid, cluster = cluster,
      heterogeneity_multiplier = heterogeneity_multiplier,
      label = paste0("Computing ", basename(path))
    )
  } else {
    message("Computing counterfactual: ", basename(path))
    cf(dat = dat, z.grid = grid, heterogeneity_multiplier = heterogeneity_multiplier)
  }
  write_counterfactual_cache(result, path, cache_key)
  result
}

make_figure_1_panel <- function(meta_average_multiplier, heterogeneity_multiplier,
                                setup_label) {
  dat <- load_multilevel_data(setup_label) %>%
    mutate(GE = meta_average_multiplier * GE) %>%
    filter_counterfactual_data(
      heterogeneity_multiplier,
      context = paste("supplement Figure 1", setup_label)
    )
  split_dat <- split(dat, dat$cID)
  z_grid <- seq(-10.25, 10.25, 0.1025)
  z_grid_positive <- seq(0, 10.25, 0.1025)
  observed_z <- abs(dat$yi / sqrt(dat$vi))
  observed_counts <- count_intervals(observed_z, z_grid_positive)

  filename_suffix <- paste0(setup_label, "_multilevel_random.rds")
  counterfactual <- get_counterfactual(
    here("data", "derived_data", paste0("z_plot_", filename_suffix)),
    split_dat, z_grid,
    heterogeneity_multiplier = heterogeneity_multiplier
  )
  counterfactual_ci <- get_counterfactual(
    here("data", "derived_data", paste0("z_plot_ci_", filename_suffix)),
    split_dat, z_grid, ci = TRUE, cluster = unique(dat$cID),
    heterogeneity_multiplier = heterogeneity_multiplier
  )

  bin_midpoints <- z_grid_positive[-length(z_grid_positive)] + diff(z_grid_positive)[1] / 2
  ## Match cf.ci.cluster(), which divides each bootstrap replication by its
  ## sampled number of effects. A count over the finite plotting grid omits any
  ## observed |z| above its upper boundary and gives the point curve a different
  ## denominator from its confidence interval.
  n_tests <- nrow(dat)
  plot_data <- tibble(
    z = bin_midpoints,
    ci_lower = apply(counterfactual_ci[[1]], 2, quantile, na.rm = TRUE, probs = 0.025),
    ci_upper = apply(counterfactual_ci[[1]], 2, quantile, na.rm = TRUE, probs = 0.975),
    observed = observed_counts / n_tests,
    counterfactual = as.vector(counterfactual) / n_tests
  )

  ggplot(plot_data) +
    geom_line(aes(z, ci_lower), color = "orange", linetype = 3) +
    geom_line(aes(z, counterfactual), color = "orange", linetype = 1) +
    geom_point(aes(z, counterfactual), shape = 20, fill = "orange", color = "orange", size = 1) +
    geom_line(aes(z, ci_upper), color = "orange", linetype = 3) +
    geom_line(aes(z, observed), color = "blue", linetype = 2) +
    geom_point(aes(z, observed), shape = 20, fill = "blue", color = "blue", size = 1) +
    coord_cartesian(xlim = c(0, 8)) +
    xlab("|z|-value") + ylab("Frequency") +
    ggtitle(paste0(heterogeneity_multiplier * 100, "% genuine heterogeneity")) +
    geom_vline(
      xintercept = c(1.64, 1.96, 2.58), linetype = 2,
      color = c(3, 2, 6), linewidth = 0.5
    ) +
    scale_x_continuous(
      breaks = c(0, 1.64, 1.96, 2.58, 4, 6, 8),
      guide = guide_axis(n.dodge = 2)
    ) +
    theme(
      panel.background = element_rect(fill = "gray100"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.5, color = "gray")
    )
}

write_counterfactual_figure <- function(meta_multiplier, figure_number) {
  figure_setups <- analysis_setups %>%
    filter(
      meta_average_multiplier == meta_multiplier,
      heterogeneity_multiplier %in% c(0, 0.5)
    ) %>%
    arrange(heterogeneity_multiplier)

  if (nrow(figure_setups) != 2) {
    stop("Expected heterogeneity multipliers 0 and 0.5 for meta-average multiplier ", meta_multiplier, ".")
  }
  panels <- pmap(figure_setups, function(meta_average_multiplier,
                                         heterogeneity_multiplier,
                                         setup_label, ...) {
    make_figure_1_panel(meta_average_multiplier, heterogeneity_multiplier, setup_label)
  })
  combined <- arrangeGrob(grobs = panels, ncol = 1)
  save_supplement_plot(
    file.path(supplement_dir, paste0("Figure_S", figure_number)),
    width = 10, height = 10,
    draw = function() grid::grid.draw(combined)
  )
}

## Figures S1 and S2: Figure 1 sensitivity analyses.
write_counterfactual_figure(meta_multiplier = 0.25, figure_number = 1)
write_counterfactual_figure(meta_multiplier = 1, figure_number = 2)

## Supplementary Figure 2 variants: power distributions for every estimator and
## meta-average multiplier not used for main-text Figure 2. Heterogeneity does
## not enter the power calculation, so these are intentionally generated only
## from the zero-heterogeneity setups.
load_power_data <- function(estimator, setup_label, meta_average_multiplier) {
  power_path <- here(
    "data", "derived_data",
    paste0("pps_rstandard_power_", setup_label, "_", estimator, ".rds")
  )
  if (!file.exists(power_path)) {
    stop("Missing ", power_path, ". Run create_tables_and_figures.R first.")
  }
  readRDS(power_path)
}

write_supplement_figure_2 <- function(meta_average_multiplier,
                                      heterogeneity_multiplier,
                                      setup_label, estimator, ...) {
  power_data <- load_power_data(estimator, setup_label, meta_average_multiplier)
  plot_data <- power_data %>%
    ## Meta-analyses without any finite power estimates produce NaN summaries.
    ## Exclude them before aggregation rather than letting each output device's
    ## stat_bin() remove the same non-finite histogram row with a warning.
    filter(is.finite(power)) %>%
    group_by(cID) %>%
    summarise(
      median = median(power),
      sape = mean(power >= 0.8),
      .groups = "drop"
    ) %>%
    mutate(
      adequately_powered = median >= 0.8,
      median = 100 * median,
      sape = 100 * sape
    )

  median_plot <- ggplot(plot_data, aes(median, fill = adequately_powered)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins = 30,
                   alpha = 0.6, linewidth = 0.1) +
    scale_fill_manual(values = c("brown2", "skyblue2")) +
    scale_x_continuous(breaks = scales::breaks_width(20),
                       labels = scales::label_percent(scale = 1)) +
    scale_y_continuous(labels = scales::label_percent(scale = 1)) +
    labs(x = "Median statistical power of primary estimates per meta-analysis",
         y = "Percentage", title = "(a)") +
    theme(legend.position = "none", panel.background = element_rect(fill = "white"))
  sape_plot <- ggplot(plot_data, aes(sape)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins = 30,
                   alpha = 0.6, linewidth = 0.1, fill = "skyblue2") +
    scale_x_continuous(breaks = scales::breaks_width(20),
                       labels = scales::label_percent(scale = 1)) +
    scale_y_continuous(labels = scales::label_percent(scale = 1)) +
    labs(x = "Share of adequately powered primary estimates per meta-analysis",
         y = "Percentage", title = "(b)") +
    theme(panel.background = element_rect(fill = "white"))

  stem <- paste0("Figure_2_", setup_label, "_", estimator)
  openxlsx::write.xlsx(plot_data, file.path(
    derived_data_dir, paste0("Figure_2_data_", setup_label, "_", estimator, ".xlsx")
  ),
                       overwrite = TRUE)
  combined <- arrangeGrob(median_plot, sape_plot, ncol = 2)
  save_supplement_plot(
    file.path(derived_data_dir, stem), width = 10, height = 4,
    draw = function() grid::grid.draw(combined)
  )
}

tidyr::crossing(
  analysis_setups %>% filter(heterogeneity_multiplier == 0),
  estimator = meta_analysis_estimators
) %>%
  filter(!(meta_average_multiplier == 0.5 & estimator == "multilevel_random")) %>%
  pwalk(write_supplement_figure_2)

## Figure S3: excess-significance estimates across all analysis setups. The
## estimates are calculated with Table 2 in create_tables_and_figures.R, while
## all supplementary output is deliberately written here.
if (!exists("all_combination_results") || !exists("all_esr_results")) {
  figure_s3_inputs_path <- here(
    "data", "derived_data", "Figure_S3_inputs.rds"
  )
  if (!file.exists(figure_s3_inputs_path)) {
    stop(
      "Figure S3 inputs are unavailable. Run create_tables_and_figures.R first ",
      "to create ", figure_s3_inputs_path, "."
    )
  }
  figure_s3_inputs <- readRDS(figure_s3_inputs_path)
  required_figure_s3_inputs <- c(
    "all_combination_results", "all_esr_results"
  )
  if (!all(required_figure_s3_inputs %in% names(figure_s3_inputs))) {
    stop(
      "Figure S3 input file is incomplete. Rerun create_tables_and_figures.R ",
      "to recreate ", figure_s3_inputs_path, "."
    )
  }
  all_combination_results <- figure_s3_inputs$all_combination_results
  all_esr_results <- figure_s3_inputs$all_esr_results
}
write.csv(
  all_combination_results,
  file.path(derived_data_dir, "Figure_S3_numbers.csv"),
  row.names = FALSE
)
esr_plot_data <- all_esr_results %>%
  filter(measure == "ESR_{0.05}^{sig}") %>%
  mutate(
    estimate = as.numeric(estimate),
    ci_lower = as.numeric(stringr::str_match(confidence_interval, "\\[([^,]+),")[, 2]),
    ci_upper = as.numeric(stringr::str_match(confidence_interval, ", ([^]]+)\\]")[, 2]),
    estimator = dplyr::recode(
      estimator,
      pet_peese = "PET-PEESE",
      multilevel_random = "Random effects"
    )
  ) %>%
  rename(`Meta average multiplier` = meta_average_multiplier)
esr_plot <- ggplot(
  esr_plot_data,
  aes(heterogeneity_multiplier, estimate, color = estimator)
) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), width = 0.03,
    position = position_dodge(width = 0.06)
  ) +
  geom_point(position = position_dodge(width = 0.06)) +
  facet_grid(. ~ `Meta average multiplier`, labeller = label_both) +
  labs(
    x = "Heterogeneity multiplier", y = expression(ESR[0.05]^sig),
    color = "Estimator"
  ) +
  theme_bw()
save_supplement_plot(
  file.path(supplement_dir, "Figure_S3"),
  width = 10, height = 4.5, draw = function() print(esr_plot)
)

## Figures S4-S10: the Figure S3 analysis repeated separately for all seven
## environmental-science subfields.
subfield_levels <- c(
  "Ecology", "Environmental Chemistry", "Environmental Engineering",
  "Health, Toxicology and Mutagenesis", "Management, Monitoring, Policy and Law",
  "Nature and Landscape Conservation", "Water Science and Technology"
)
subfield_esr_path <- here(
  "data", "derived_data", "Figure_S3_subfield_inputs.rds"
)
if (!file.exists(subfield_esr_path)) {
  stop("Subfield Figure S3 inputs are unavailable. Run create_tables_and_figures.R first.")
}
subfield_esr_results <- readRDS(subfield_esr_path)
walk2(subfield_levels, 4:10, function(subfield, figure_number) {
  plot_data <- subfield_esr_results %>%
    filter(.data$subfield == .env$subfield, measure == "ESR_{0.05}^{sig}") %>%
    mutate(
      estimate = as.numeric(estimate),
      ci_lower = as.numeric(str_match(confidence_interval, "\\[([^,]+),")[, 2]),
      ci_upper = as.numeric(str_match(confidence_interval, ", ([^]]+)\\]")[, 2]),
      estimator = dplyr::recode(estimator, pet_peese = "PET-PEESE",
                         multilevel_random = "Random effects")
    ) %>%
    rename(`Meta average multiplier` = meta_average_multiplier)
  plot <- ggplot(plot_data, aes(heterogeneity_multiplier, estimate, color = estimator)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = .03,
                  position = position_dodge(width = .06)) +
    geom_point(position = position_dodge(width = .06)) +
    facet_grid(. ~ `Meta average multiplier`, labeller = label_both) +
    labs(x = "Heterogeneity multiplier", y = expression(ESR[0.05]^sig),
         color = "Estimator", title = subfield) +
    theme_bw()
  save_supplement_plot(
    file.path(supplement_dir, paste0("Figure_S", figure_number)),
    width = 10, height = 4.5, draw = function() print(plot)
  )
})

## Tables and remaining figures follow the order of the supplementary material.
## Every table is emitted as both an editable Word document and copy-ready LaTeX.
library(officer)

write_word_table <- function(x, number, alignment = NULL) {
  if (is.null(alignment)) alignment <- c("left", rep("center", ncol(x) - 1))
  doc <- officer::read_docx()
  doc <- officer::body_add_table(doc, x, style = NULL, header = TRUE,
    alignment = alignment, align_table = "center")
  print(doc, target = file.path(supplement_dir, paste0("Table_S", number, ".docx")))
  write_latex_table(x, file.path(supplement_dir, paste0("Table_S", number, ".tex")))
}

make_power_table <- function(dat) {
  dat <- dat %>% mutate(
    sei = sqrt(vi), GE = 0.5 * GE,
    power = 1 - pnorm(qnorm(.975) - abs(GE) / sei) +
      pnorm(qnorm(.025) - abs(GE) / sei)
  )
  by_meta <- dat %>% group_by(cID) %>% summarise(
    Subfield = first(subfd), meta_median = median(power, na.rm = TRUE),
    .groups = "drop")
  detail <- dat %>% group_by(subfd) %>% summarise(
    `No. of meta-analyses` = n_distinct(cID),
    `No. of primary estimates` = n(), Median = median(power, na.rm = TRUE),
    Mean = mean(power, na.rm = TRUE), Q25 = quantile(power, .25, na.rm = TRUE),
    Q75 = quantile(power, .75, na.rm = TRUE), SAPE = mean(power >= .8, na.rm = TRUE),
    .groups = "drop") %>% rename(Subfield = subfd)
  bind_rows(
    tibble(Subfield = "All meta-analyses",
      `No. of meta-analyses` = n_distinct(dat$cID),
      `No. of primary estimates` = nrow(dat), Median = median(dat$power, na.rm = TRUE),
      Mean = mean(dat$power, na.rm = TRUE), Q25 = quantile(dat$power, .25, na.rm = TRUE),
      Q75 = quantile(dat$power, .75, na.rm = TRUE), SAPE = mean(dat$power >= .8, na.rm = TRUE)),
    detail %>% mutate(Subfield = factor(Subfield, subfield_levels)) %>% arrange(Subfield) %>%
      mutate(Subfield = as.character(Subfield))
  ) %>% left_join(
    bind_rows(
      tibble(Subfield = "All meta-analyses", `Median of medians` = median(by_meta$meta_median, na.rm = TRUE)),
      by_meta %>% group_by(Subfield) %>% summarise(`Median of medians` = median(meta_median, na.rm = TRUE), .groups = "drop")
    ), by = "Subfield"
  ) %>% dplyr::select(Subfield, `No. of meta-analyses`, `No. of primary estimates`,
    `Median of medians`, Median, Mean, Q25, Q75, SAPE) %>%
    mutate(across(`Median of medians`:SAPE, ~ sprintf("%.2f", .x)))
}

base_half <- load_multilevel_data("meta_0p5_heterogeneity_0")
## Table S1: remove primary estimates that are not significant at five percent.
write_word_table(make_power_table(base_half %>% filter(abs(yi / sqrt(vi)) > 1.96)), 1)
## Table S2: remove complete meta-analyses whose pooled effect is not significant.
significant_meta <- base_half %>% distinct(cID, sig_overall) %>%
  filter(!is.na(sig_overall), sig_overall < .05) %>% pull(cID)
write_word_table(make_power_table(base_half %>% filter(cID %in% significant_meta)), 2)

## Table S3: adequately powered meta-analyses and small-study effects by subfield.
small_study_effects <- load_small_study_effects("meta_0p5_heterogeneity_0")
table_s3_meta <- base_half %>%
  dplyr::select(-any_of(c("small_study_effect_pval", "sse_yn"))) %>%
  left_join(small_study_effects, by = "cID") %>% mutate(
  power = 1 - pnorm(qnorm(.975) - abs(.5 * GE) / sqrt(vi)) +
    pnorm(qnorm(.025) - abs(.5 * GE) / sqrt(vi))) %>%
  group_by(cID) %>% summarise(Subfield = first(subfd),
    adequately_powered = median(power, na.rm = TRUE) >= .8,
    small_study_effect = first(small_study_effect_pval) <= .05,
    .groups = "drop")
table_s3 <- table_s3_meta %>% group_by(Subfield) %>% summarise(
  `Median power >= 80% (%)` = 100 * mean(adequately_powered),
  `Small-study effects (%)` = 100 * mean(small_study_effect, na.rm = TRUE),
  .groups = "drop") %>% bind_rows(tibble(Subfield = "All meta-analyses",
    `Median power >= 80% (%)` = 100 * mean(table_s3_meta$adequately_powered),
    `Small-study effects (%)` = 100 * mean(table_s3_meta$small_study_effect, na.rm = TRUE))) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.1f", .x)))
write_word_table(table_s3, 3)

## Figure S11: heterogeneity distributions and the corresponding summaries.
heterogeneity_data <- base_half %>% distinct(cID, subfd, isq) %>%
  mutate(subfd = factor(subfd, subfield_levels))
heterogeneity_summary <- heterogeneity_data %>%
  group_by(subfd) %>%
  summarise(
    `No. of meta-analyses` = n(), Median = median(isq, na.rm = TRUE),
    Mean = mean(isq, na.rm = TRUE),
    Q25 = quantile(isq, .25, na.rm = TRUE),
    Q75 = quantile(isq, .75, na.rm = TRUE), .groups = "drop"
  ) %>%
  arrange(subfd)
openxlsx::write.xlsx(
  heterogeneity_summary %>% rename(Subfield = subfd),
  file.path(derived_data_dir, "Figure_S11_heterogeneity_by_subfield.xlsx"),
  overwrite = TRUE
)

figure_s4_rows <- heterogeneity_summary %>%
  mutate(row = rev(seq_len(n())))

## Draw the table and ridgelines in one coordinate system. Keeping every visual
## element in the same panel makes the row centres identical by construction;
## separate table and plot grobs can acquire different header and cell heights.
density_scale <- c(80, 99)
figure_s4_densities <- heterogeneity_data %>%
  filter(!is.na(isq)) %>%
  group_by(subfd) %>%
  group_modify(~ {
    curve <- density(.x$isq, from = 0, to = 100, adjust = .8, n = 256)
    tibble(
      x = scales::rescale(curve$x, to = density_scale, from = c(0, 100)),
      height = curve$y / max(curve$y)
    )
  }) %>%
  ungroup() %>%
  left_join(dplyr::select(figure_s4_rows, subfd, row), by = "subfd") %>%
  mutate(y = row + .34 * height)

column_positions <- c(Subfield = 1, `No. of meta-analyses` = 46,
                      Median = 55, Mean = 63,
                      Q25 = 70, Q75 = 77, Heterogeneity = 89.5)
heterogeneity_plot <- ggplot() +
  geom_hline(
    yintercept = c(.5, seq(1.5, nrow(figure_s4_rows) + .5),
                   nrow(figure_s4_rows) + 1.35),
    colour = "#d0d0d0", linewidth = .45
  ) +
  geom_ribbon(
    data = figure_s4_densities,
    aes(x = x, ymin = row, ymax = y, group = subfd),
    fill = "#66c2df", colour = "#b5b5b5", linewidth = .55
  ) +
  geom_text(
    data = figure_s4_rows,
    aes(x = column_positions[["Subfield"]], y = row, label = subfd),
    hjust = 0, size = 3.6
  ) +
  geom_text(
    data = figure_s4_rows,
    aes(x = column_positions[["No. of meta-analyses"]], y = row,
        label = `No. of meta-analyses`),
    hjust = 1, size = 3.6
  ) +
  geom_text(data = figure_s4_rows,
            aes(x = column_positions[["Median"]], y = row, label = sprintf("%.2f", Median)),
            hjust = 1, size = 3.6) +
  geom_text(data = figure_s4_rows,
            aes(x = column_positions[["Mean"]], y = row, label = sprintf("%.2f", Mean)),
            hjust = 1, size = 3.6) +
  geom_text(data = figure_s4_rows,
            aes(x = column_positions[["Q25"]], y = row, label = sprintf("%.2f", Q25)),
            hjust = 1, size = 3.6) +
  geom_text(data = figure_s4_rows,
            aes(x = column_positions[["Q75"]], y = row, label = sprintf("%.2f", Q75)),
            hjust = 1, size = 3.6) +
  annotate("text", x = column_positions, y = nrow(figure_s4_rows) + 1,
           label = names(column_positions),
           hjust = c(0, rep(.5, length(column_positions) - 1)), size = 3.6) +
  coord_cartesian(xlim = c(0, 100), ylim = c(.45, nrow(figure_s4_rows) + 1.35),
                  expand = FALSE, clip = "off") +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(6, 12, 6, 12, unit = "mm"),
    plot.caption = element_text(hjust = .5, margin = margin(t = 12))
  )

save_supplement_plot(file.path(supplement_dir, "Figure_S11"),
              9, 5.25, function() print(heterogeneity_plot))

## Figures S12-S13: subfield counterfactual distributions with zero and 50%
## genuine heterogeneity. The obsolete ESR Tables S4-S6 are no longer created.
subfield_grids <- list(
  p = c(-Inf, qnorm(c(.001, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9) / 2),
    0, qnorm(c(.9, .8, .7, .6, .5, .4, .3, .2, .1, .05, .01, .001) / 2,
      lower.tail = FALSE), Inf),
  p_absolute = c(0,
    qnorm(c(.9, .8, .7, .6, .5, .4, .3, .2, .1, .05, .01, .001) / 2,
      lower.tail = FALSE), Inf),
  z = seq(-10.25, 10.25, .1025),
  z_absolute = seq(0, 10.25, .1025)
)
legacy_codes <- c(eco = "Ecology", enc = "Environmental Chemistry",
  ene = "Environmental Engineering", htm = "Health, Toxicology and Mutagenesis",
  mpl = "Management, Monitoring, Policy and Law", nlc = "Nature and Landscape Conservation",
  wst = "Water Science and Technology")

prepare_subfield_data <- function(label, heterogeneity_multiplier) {
  setup_label <- if (heterogeneity_multiplier == 0) {
    "meta_0p5_heterogeneity_0"
  } else {
    "meta_0p5_heterogeneity_0p5"
  }
  load_multilevel_data(setup_label) %>%
    filter(subfd == label) %>%
    mutate(GE = .5 * GE) %>%
    filter_counterfactual_data(
      heterogeneity_multiplier, context = paste("supplement subfield", label)
    )
}

calculate_subfield_counterfactual <- function(label, code, grid, type,
                                               heterogeneity_multiplier,
                                               ci = FALSE) {
  dat <- prepare_subfield_data(label, heterogeneity_multiplier)
  heterogeneity_suffix <- paste0("_heterogeneity_", heterogeneity_multiplier)
  get_counterfactual(
    file.path(derived_data_dir, paste0("pet_peese_rstandard_", type,
      if (ci) "_ci" else "", ".", code, heterogeneity_suffix, ".rds")),
    split(dat, dat$cID), grid, ci = ci,
    cluster = if (ci) unique(dat$cID) else NULL,
    heterogeneity_multiplier = heterogeneity_multiplier
  )
}
make_subfield_counterfactual_plot <- function(heterogeneity_multiplier) {
 subfield_plot_data <- imap_dfr(legacy_codes, function(label, code) {
  dat <- prepare_subfield_data(label, heterogeneity_multiplier)
  point <- calculate_subfield_counterfactual(
    label, code, subfield_grids$z, "z_plot", heterogeneity_multiplier
  )
  ci <- calculate_subfield_counterfactual(
    label, code, subfield_grids$z, "z_plot", heterogeneity_multiplier, TRUE
  )
  observed <- abs(dat$yi / sqrt(dat$vi))
  tibble(Subfield = label,
    z = head(subfield_grids$z_absolute, -1) + diff(subfield_grids$z_absolute)[1] / 2,
    factual = count_intervals(observed, subfield_grids$z_absolute) / length(observed),
    counterfactual = as.vector(point) / length(observed),
    lower = apply(ci[[1]], 2, quantile, .025, na.rm = TRUE),
    upper = apply(ci[[1]], 2, quantile, .975, na.rm = TRUE))
 })
 ggplot(subfield_plot_data, aes(z)) +
  geom_line(aes(y = lower), colour = "orange", linetype = 3) +
  geom_line(aes(y = counterfactual), colour = "orange") +
  geom_point(aes(y = counterfactual), colour = "orange", size = 1) +
  geom_line(aes(y = upper), colour = "orange", linetype = 3) +
  geom_line(aes(y = factual), colour = "blue", linetype = 2) +
  geom_point(aes(y = factual), colour = "blue", size = 1) +
  geom_vline(
    data = tibble(
      xintercept = c(1.64, 1.96, 2.58),
      threshold_colour = palette()[c(3, 2, 6)]
    ),
    aes(xintercept = xintercept, colour = threshold_colour),
    linetype = 2, linewidth = .5
  ) +
  scale_colour_identity() +
  facet_wrap(~ Subfield, ncol = 2, scales = "free_y") + coord_cartesian(xlim = c(0, 8)) +
  scale_x_continuous(
    breaks = c(0, 1.64, 1.96, 2.58, 4, 6, 8),
    guide = guide_axis(n.dodge = 2)
  ) +
  labs(x = "|z|-value", y = "Frequency",
       title = paste0(heterogeneity_multiplier * 100, "% genuine heterogeneity")) +
  theme(
    panel.background = element_rect(fill = "gray100"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = .5, colour = "gray")
  )
}
walk2(c(0, .5), 12:13, function(heterogeneity_multiplier, figure_number) {
  plot <- make_subfield_counterfactual_plot(heterogeneity_multiplier)
  save_supplement_plot(file.path(supplement_dir, paste0("Figure_S", figure_number)),
    11, 10, function() print(plot))
})

## Regression table formatting for Tables S4-S6. The negative-binomial fits
## consumed by Tables S4-S5 are created in run_exploratory_regressions.R.
significance_stars <- function(p_value) {
  ifelse(p_value < .01, "***", ifelse(p_value < .05, "**",
    ifelse(p_value < .10, "*", "")))
}
format_model_table <- function(fit, include_adjusted_r2 = FALSE) {
  ## Force the fitted-model bundle before doing any validation. In particular,
  ## this makes the function safe to step through with debug()/debugonce()
  ## without repeatedly restarting evaluation of the lazy `fit` promise.
  models <- fit[["models"]]
  robust_tables <- fit[["robust"]]
  labels <- c("Intercept" = "(Intercept)", "Median power" = "med_perc",
    "Experimental research design? (yes)" = "design_mergedyes", "Followed reporting guidelines? (yes)" = "guidyes",
    "Protocol registered? (yes)" = "preryes", "Log number of independent studies" = "lognps",
    "Log journal impact factor" = "logjif", "Publication year" = "pyear",
    setNames(paste0("subf", subfield_levels[-1]), subfield_levels[-1]))
  if (length(models) != 4 || length(robust_tables) != 4) {
    stop("Expected four fitted models and four robust coefficient tables.")
  }
  format_coefficient <- function(coefficient_table, term) {
    coefficient_table <- as.matrix(coefficient_table)
    if (is.null(rownames(coefficient_table)) ||
        !term %in% rownames(coefficient_table)) {
      return("")
    }
    if (ncol(coefficient_table) < 4) {
      stop("Robust coefficient tables must contain estimate, SE, statistic, and p-value columns.")
    }
    estimate <- unname(coefficient_table[term, 1])
    standard_error <- unname(coefficient_table[term, 2])
    p_value <- unname(coefficient_table[term, 4])
    sprintf(
      "%.3f%s (%.3f)", estimate, significance_stars(p_value), standard_error
    )
  }
  result <- tibble(Variable = names(labels))
  column_names <- paste0(
    rep(c("0%", "50%"), each = 2), " heterogeneity: Model ",
    rep(1:2, 2), " Estimate (SE)"
  )
  for (i in seq_len(4)) {
    result[[column_names[[i]]]] <- vapply(
      unname(labels),
      function(term) format_coefficient(robust_tables[[i]], term),
      character(1),
      USE.NAMES = FALSE
    )
  }
  statistic <- if (include_adjusted_r2) "Adjusted R-squared" else "AIC"
  summary_rows <- tibble(Variable = c(
    "Effect size type", statistic, "No. of meta-analyses"
  ))
  for (i in seq_len(4)) {
    fit_statistic <- if (include_adjusted_r2) {
      summary(models[[i]])$adj.r.squared
    } else {
      stats::AIC(models[[i]])
    }
    summary_rows[[column_names[[i]]]] <- c(
      "Yes", sprintf("%.3f", fit_statistic), stats::nobs(models[[i]])
    )
  }
  bind_rows(result, summary_rows)
}
write_word_table(format_model_table(nb_sensitivity_full), 4)
write_word_table(format_model_table(nb_sensitivity_quarter), 5)

ols_formula <- esr_winsor ~ med_perc + design_merged + guid + prer + lognps + logjif + pyear + metric
ols_fits <- map(c(0, .5), function(heterogeneity_multiplier) {
  ols_data <- build_regression_data(.5, heterogeneity_multiplier) %>%
    mutate(esr_winsor = pmin(
      pmax(esr.sig, quantile(esr.sig, .05)), quantile(esr.sig, .95)
    ))
  models <- list(
    lm(ols_formula, ols_data),
    lm(update(ols_formula, . ~ . + subf), ols_data)
  )
  list(models = models, robust = map(models,
    ~ lmtest::coeftest(.x, vcov. = sandwich::vcovCL(.x, cluster = ols_data$metaID))))
})
ols_fit <- list(
  models = flatten(map(ols_fits, "models")),
  robust = flatten(map(ols_fits, "robust"))
)
write_word_table(format_model_table(ols_fit, TRUE), 6)

## Figures S14-S17: separate continuous and categorical diagnostics for both main
## negative-binomial specifications, matching the requested four-figure layout.
save_diagnostic_group <- function(model, model_number, kind, figure_number) {
  dat <- final_nb %>% mutate(residual = residuals(model), fitted_value = fitted(model))
  vars <- if (kind == "continuous") c("fitted_value", "med_perc", "lognps", "logtotall", "logjif", "pyear") else c("design_merged", "guid", "prer", "subf")
  plots <- map(vars, function(v) {
    if (is.numeric(dat[[v]])) ggplot(dat, aes(.data[[v]], residual)) + geom_point(colour = "skyblue3", shape = 1) + geom_hline(yintercept = 0, colour = "red") + theme_bw() + labs(x = v)
    else ggplot(dat, aes(.data[[v]], residual)) + geom_boxplot() + geom_hline(yintercept = 0, colour = "red") + theme_bw() + labs(x = v)
  })
  grob <- arrangeGrob(grobs = plots, ncol = 2)
  save_supplement_plot(file.path(supplement_dir, paste0("Figure_S", figure_number)),
    11, ifelse(kind == "continuous", 10, 7), function() grid::grid.draw(grob))
}
save_diagnostic_group(nbMod1, 1, "continuous", 14)
save_diagnostic_group(nbMod1, 1, "categorical", 15)
save_diagnostic_group(nbMod2, 2, "continuous", 16)
save_diagnostic_group(nbMod2, 2, "categorical", 17)

## Table S7 is descriptive and therefore uses all observations rather than an
## estimator-specific outlier-screened sample.
table_s7_all_data <- load_multilevel_all_data()
table_s7 <- table_s7_all_data %>% distinct(cID, etype, subfd) %>% count(etype, subfd) %>%
  complete(etype, subfd = subfield_levels, fill = list(n = 0)) %>%
  pivot_wider(names_from = subfd, values_from = n) %>% rename(`Effect size` = etype) %>%
  arrange(`Effect size`) %>% mutate(No. = row_number(), .before = 1)
write_word_table(table_s7, 7)
