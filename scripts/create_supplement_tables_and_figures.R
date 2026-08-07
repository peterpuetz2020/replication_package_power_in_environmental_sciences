## ---------------------------------------------------------
## create_supplement_tables_and_figures.R
## ---------------------------------------------------------
## Recreate the supplementary figures from the analysis outputs. Run
## create_tables_and_figures.R first so ESR_results_all_combinations.csv exists.

library(tidyverse)
library(foreach)
library(doParallel)
library(gridExtra)
library(here)
library(openxlsx)

source(here("scripts", "analysis_setup.R"))

supplement_dir <- here("results", "supplement")
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)

save_supplement_plot <- function(filename_stem, width, height, draw) {
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

load_multilevel_data <- function(setup_label) {
  derived_path <- here(
    "data", "derived_data",
    paste0("pps_rstandard_raw_", setup_label, "_multilevel_random.rds")
  )
  if (file.exists(derived_path)) {
    return(readRDS(derived_path))
  }

  input_files <- list.files(
    here("results", "main", "multilevel_random"),
    pattern = "\\.rds$",
    full.names = TRUE
  )
  if (length(input_files) == 0) {
    stop("No multilevel random-effects inputs found. Run compare_meta_analysis_estimators.R first.")
  }
  purrr::map_dfr(input_files, readRDS)
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
    file.path(supplement_dir, paste0("Figure_S", figure_number, "_multilevel_random")),
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
    group_by(cID) %>%
    summarise(
      median = median(power, na.rm = TRUE),
      sape = sum(power >= 0.8, na.rm = TRUE) / sum(!is.na(power)),
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
    supplement_dir, paste0("Figure_2_data_", setup_label, "_", estimator, ".xlsx")
  ),
                       overwrite = TRUE)
  combined <- arrangeGrob(median_plot, sape_plot, ncol = 2)
  save_supplement_plot(
    file.path(supplement_dir, stem), width = 10, height = 4,
    draw = function() grid::grid.draw(combined)
  )
}

tidyr::crossing(
  analysis_setups %>% filter(heterogeneity_multiplier == 0),
  estimator = meta_analysis_estimators
) %>%
  filter(!(meta_average_multiplier == 0.5 & estimator == "multilevel_random")) %>%
  pwalk(write_supplement_figure_2)

## Figure S3: excess significant results across every analysis combination.
esr_path <- here("results", "main", "ESR_results_all_combinations.csv")
if (!file.exists(esr_path)) {
  stop("Missing ", esr_path, ". Run create_tables_and_figures.R first.")
}

esr_plot_data_raw <- readr::read_csv(esr_path, show_col_types = FALSE) %>%
  filter(measure == "ESR_{0.05}^{sig}") %>%
  extract(confidence_interval, c("ci_lower", "ci_upper"),
          regex = "\\[([^,]+),\\s*([^]]+)\\]", convert = TRUE) %>%
  mutate(estimate = as.numeric(estimate))

if (nrow(esr_plot_data_raw) == 0) {
  stop("ESR results do not contain any ESR_{0.05}^{sig} rows.")
}

incomplete_esr_rows <- esr_plot_data_raw %>%
  filter(if_any(c(estimate, ci_lower, ci_upper), is.na))
if (nrow(incomplete_esr_rows) > 0) {
  warning(
    "Drawing ", nrow(incomplete_esr_rows),
    " ESR_{0.05}^{sig} point estimate(s) without confidence intervals in Figure S3. ",
    "Regenerate ESR_results_all_combinations.csv with create_tables_and_figures.R ",
    "to restore missing point estimates."
  )
}

inconsistent_esr_rows <- esr_plot_data_raw %>%
  filter(
    if_all(c(estimate, ci_lower, ci_upper), ~ !is.na(.x)),
    estimate < ci_lower | estimate > ci_upper
  )
if (nrow(inconsistent_esr_rows) > 0) {
  stop(
    "ESR_results_all_combinations.csv contains ", nrow(inconsistent_esr_rows),
    " point estimate(s) outside their confidence intervals. The counterfactual ",
    "cache that produced this file is stale; rerun create_tables_and_figures.R ",
    "to rebuild the cache and ESR results before creating supplement figures."
  )
}

esr_plot_data <- esr_plot_data_raw %>%
  mutate(
    estimator = dplyr::recode(
      estimator,
      multilevel_random = "Random effects",
      pet_peese = "PET-PEESE"
    ),
    meta_average_multiplier = factor(
      meta_average_multiplier,
      levels = sort(unique(meta_average_multiplier))
    )
  )

if (nrow(esr_plot_data) == 0) {
  stop("ESR results do not contain plottable ESR_{0.05}^{sig} estimates.")
}

esr_plot <- ggplot(
  esr_plot_data,
  aes(heterogeneity_multiplier, estimate, color = estimator, group = estimator)
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_line(position = position_dodge(width = 0.035), linewidth = 0.55) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), width = 0.025,
    data = ~ filter(.x, if_all(c(ci_lower, ci_upper), ~ !is.na(.x))),
    position = position_dodge(width = 0.035), linewidth = 0.5
  ) +
  geom_point(position = position_dodge(width = 0.035), size = 2) +
  facet_wrap(
    ~ meta_average_multiplier,
    labeller = labeller(meta_average_multiplier = function(x) paste("Meta-average multiplier =", x))
  ) +
  scale_color_manual(values = c("Random effects" = "#0072B2", "PET-PEESE" = "#D55E00")) +
  scale_x_continuous(breaks = sort(unique(esr_plot_data$heterogeneity_multiplier))) +
  labs(
    x = "Heterogeneity multiplier",
    y = expression(ESR[0.05]^sig),
    color = "Estimator"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  )

save_supplement_plot(
  file.path(supplement_dir, "Figure_S3_excess_p"),
  width = 10, height = 4.5,
  draw = function() print(esr_plot)
)
