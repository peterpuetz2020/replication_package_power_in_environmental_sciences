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

get_counterfactual <- function(path, dat, grid, ci = FALSE, cluster = NULL,
                               heterogeneity_multiplier) {
  if (file.exists(path)) {
    return(readRDS(path))
  }

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  result <- if (ci) {
    cf.ci.cluster(
      dat = dat, z.grid = grid, iters = n_iterations, cluster = cluster,
      heterogeneity_multiplier = heterogeneity_multiplier
    )
  } else {
    cf(dat = dat, z.grid = grid, heterogeneity_multiplier = heterogeneity_multiplier)
  }
  saveRDS(result, path)
  result
}

make_figure_1_panel <- function(meta_average_multiplier, heterogeneity_multiplier,
                                setup_label) {
  dat <- load_multilevel_data(setup_label) %>%
    mutate(GE = meta_average_multiplier * GE)
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
  n_tests <- sum(observed_counts)
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
    width = 5, height = 10,
    draw = function() grid::grid.draw(combined)
  )
}

## Figures S1 and S2: Figure 1 sensitivity analyses.
write_counterfactual_figure(meta_multiplier = 0.25, figure_number = 1)
write_counterfactual_figure(meta_multiplier = 1, figure_number = 2)

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
    "Omitting ", nrow(incomplete_esr_rows),
    " incomplete ESR_{0.05}^{sig} row(s) from Figure S3. ",
    "Regenerate ESR_results_all_combinations.csv with create_tables_and_figures.R ",
    "to restore missing point estimates."
  )
}

esr_plot_data <- esr_plot_data_raw %>%
  filter(if_all(c(estimate, ci_lower, ci_upper), ~ !is.na(.x))) %>%
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
  stop("ESR results do not contain plottable ESR_{0.05}^{sig} estimates and confidence intervals.")
}

esr_plot <- ggplot(
  esr_plot_data,
  aes(heterogeneity_multiplier, estimate, color = estimator, group = estimator)
) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_line(position = position_dodge(width = 0.035), linewidth = 0.55) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), width = 0.025,
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
