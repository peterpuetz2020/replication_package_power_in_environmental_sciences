## ---------------------------------------------------------
## create_supplement_tables_and_figures.R
## ---------------------------------------------------------
## Recreate all supplementary tables and figures from the analysis outputs. The
## fitted negative-binomial models used for the supplementary tables are produced by
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
    rows <- apply(x, 1, function(row) {
      paste0(paste(latex_escape(row), collapse = " & "), " \\\\")
    })
    writeLines(rows, path, useBytes = TRUE)
  }
}

supplement_dir <- here("results", "supplement")
dir.create(supplement_dir, recursive = TRUE, showWarnings = FALSE)
derived_data_dir <- here("data", "derived_data")
dir.create(derived_data_dir, recursive = TRUE, showWarnings = FALSE)

## Use one font size for every text element in the supplementary figures.
supplement_font_size <- 12
supplement_figure_theme <- function() {
  theme_minimal(base_size = supplement_font_size) +
    theme(text = element_text(size = supplement_font_size))
}

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

  ## Write each scalable format to a temporary file first. If drawing fails,
  ## this keeps a truncated device output from masquerading as a valid figure.
  save_with_device(".pdf", function(path) pdf(path, width = width, height = height))
  save_with_device(".svg", function(path) svg(path, width = width, height = height))
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
    split_dat, z_grid, ci = TRUE,
    cluster = counterfactual_meta_clusters(split_dat),
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
    supplement_figure_theme() +
    theme(
      panel.background = element_rect(fill = "gray100"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.5, color = "gray")
    )
}

write_counterfactual_figure <- function(meta_multiplier, figure_number) {
  figure_setups <- tibble(
    meta_average_multiplier = meta_multiplier,
    heterogeneity_multiplier = c(0, 0.5)
  ) %>%
    mutate(setup_label = make_setup_label(
      meta_average_multiplier, heterogeneity_multiplier
    ))
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

subfield_levels <- c(
  "Ecology", "Environmental Chemistry", "Environmental Engineering",
  "Health, Toxicology and Mutagenesis", "Management, Monitoring, Policy and Law",
  "Nature and Landscape Conservation", "Water Science and Technology"
)
base_half <- load_multilevel_data("meta_0p5_heterogeneity_0")

## Figure S1: heterogeneity distributions and the corresponding summaries.
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
  file.path(derived_data_dir, "Figure_S1_heterogeneity_by_subfield.xlsx"),
  overwrite = TRUE
)

figure_s1_rows <- heterogeneity_summary %>%
  mutate(row = rev(seq_len(n())))

## Draw the table and ridgelines in one coordinate system. Keeping every visual
## element in the same panel makes the row centres identical by construction;
## separate table and plot grobs can acquire different header and cell heights.
density_scale <- c(91, 111)
figure_s1_xlim <- c(0, max(density_scale) + 1)
density_height <- 0.34

figure_s1_densities <- heterogeneity_data %>%
  filter(!is.na(isq)) %>%
  group_by(subfd) %>%
  group_modify(~ {
    curve <- density(
      .x$isq,
      from = 0,
      to = 100,
      adjust = 0.8,
      n = 256
    )

    tibble(
      x = scales::rescale(
        curve$x,
        to = density_scale,
        from = c(0, 100)
      ),
      height = curve$y / max(curve$y)
    )
  }) %>%
  ungroup() %>%
  left_join(
    dplyr::select(figure_s1_rows, subfd, row),
    by = "subfd"
  ) %>%
  mutate(
    density_baseline = row - density_height / 2,
    density_top = density_baseline + density_height * height
  )


column_positions <- c(
  Subfield = 1,
  `No. of meta-analyses` = 48,
  Median = 60,
  Mean = 69,
  Q25 = 78,
  Q75 = 87,
  Heterogeneity = mean(density_scale)
)


header_df <- tibble(
  label = c(
    "Subfield",
    "No. of\nmeta-anal.",
    "Median",
    "Mean",
    "Q25",
    "Q75",
    "Heterogeneity"
  ),
  x = unname(column_positions),
  hjust = c(
    0,
    rep(0.5, 6)
  )
)


heterogeneity_plot <- ggplot() +

  # horizontal table lines
  geom_hline(
    yintercept = c(
      0.5,
      seq(
        1.5,
        nrow(figure_s1_rows) + 0.5,
        by = 1
      ),
      nrow(figure_s1_rows) + 1.35
    ),
    colour = "#d0d0d0",
    linewidth = 0.45
  ) +

  # heterogeneity distributions
  geom_ribbon(
    data = figure_s1_densities,
    aes(
      x = x,
      ymin = density_baseline,
      ymax = density_top,
      group = subfd
    ),
    fill = "#66c2df",
    colour = "#b5b5b5",
    linewidth = 0.55
  ) +

  # subfield names
  geom_text(
    data = figure_s1_rows,
    aes(
      x = column_positions[["Subfield"]],
      y = row,
      label = subfd
    ),
    hjust = 0,
    size = supplement_font_size / ggplot2::.pt
  ) +

  # number of meta-analyses
  geom_text(
    data = figure_s1_rows,
    aes(
      x = column_positions[["No. of meta-analyses"]],
      y = row,
      label = `No. of meta-analyses`
    ),
    hjust = 0.5,
    size = supplement_font_size / ggplot2::.pt
  ) +

  # median
  geom_text(
    data = figure_s1_rows,
    aes(
      x = column_positions[["Median"]],
      y = row,
      label = sprintf("%.2f", Median)
    ),
    hjust = 0.5,
    size = supplement_font_size / ggplot2::.pt
  ) +

  # mean
  geom_text(
    data = figure_s1_rows,
    aes(
      x = column_positions[["Mean"]],
      y = row,
      label = sprintf("%.2f", Mean)
    ),
    hjust = 0.5,
    size = supplement_font_size / ggplot2::.pt
  ) +

  # Q25
  geom_text(
    data = figure_s1_rows,
    aes(
      x = column_positions[["Q25"]],
      y = row,
      label = sprintf("%.2f", Q25)
    ),
    hjust = 0.5,
    size = supplement_font_size / ggplot2::.pt
  ) +

  # Q75
  geom_text(
    data = figure_s1_rows,
    aes(
      x = column_positions[["Q75"]],
      y = row,
      label = sprintf("%.2f", Q75)
    ),
    hjust = 0.5,
    size = supplement_font_size / ggplot2::.pt
  ) +

  # column headers
  geom_text(
    data = header_df,
    aes(
      x = x,
      y = nrow(figure_s1_rows) + 1.18,
      label = label,
      hjust = hjust
    ),
    vjust = 1,
    lineheight = 0.9,
    size = supplement_font_size / ggplot2::.pt
  ) +

  coord_cartesian(
    xlim = figure_s1_xlim,
    ylim = c(
      0.45,
      nrow(figure_s1_rows) + 1.35
    ),
    expand = FALSE,
    clip = "off"
  ) +

  theme_void(base_size = supplement_font_size) +

  theme(
    text = element_text(size = supplement_font_size),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(
      6, 12, 6, 12,
      unit = "mm"
    )
  )

save_supplement_plot(file.path(supplement_dir, "Figure_S1"),
              10.5, 5.25, function() print(heterogeneity_plot))

## Figures S2 and S3: Figure 1 sensitivity analyses.
write_counterfactual_figure(meta_multiplier = 0.25, figure_number = 2)
write_counterfactual_figure(meta_multiplier = 1, figure_number = 3)

## Figure S4: excess-significance estimates across all analysis setups. The
## estimates are calculated with Table 2 in create_tables_and_figures.R, while
## all supplementary output is deliberately written here.
if (!exists("all_combination_results") || !exists("all_esr_results")) {
  figure_s4_inputs_path <- here(
    "data", "derived_data", "Figure_S4_inputs.rds"
  )
  if (!file.exists(figure_s4_inputs_path)) {
    stop(
      "Figure S4 inputs are unavailable. Run create_tables_and_figures.R first ",
      "to create ", figure_s4_inputs_path, "."
    )
  }
  figure_s4_inputs <- readRDS(figure_s4_inputs_path)
  required_figure_s4_inputs <- c(
    "all_combination_results", "all_esr_results"
  )
  if (!all(required_figure_s4_inputs %in% names(figure_s4_inputs))) {
    stop(
      "Figure S4 input file is incomplete. Rerun create_tables_and_figures.R ",
      "to recreate ", figure_s4_inputs_path, "."
    )
  }
  all_combination_results <- figure_s4_inputs$all_combination_results
  all_esr_results <- figure_s4_inputs$all_esr_results
}
write.csv(
  all_combination_results,
  file.path(derived_data_dir, "Figure_S4_numbers.csv"),
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
  rename(`Meta-average multiplier` = meta_average_multiplier)
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
  facet_grid(. ~ `Meta-average multiplier`, labeller = label_both) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.5, 0.75)
  ) +
  labs(
    x = "Share of genuine heterogeneity", y = expression(ESR[0.05]^sig),
    color = "Estimator"
  ) +
  supplement_figure_theme() +
  theme(legend.position = "bottom")
save_supplement_plot(
  file.path(supplement_dir, "Figure_S4"),
  width = 10, height = 4.5, draw = function() print(esr_plot)
)

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

write_subfield_esr_table <- function(x, number, alignment = NULL) {

  # ============================================================
  # 1. WORD VERSION
  # ============================================================

  if (is.null(alignment)) {
    alignment <- c("left", rep("center", ncol(x) - 1))
  }

  doc <- officer::read_docx()

  doc <- officer::body_add_table(
    doc,
    x,
    style = NULL,
    header = TRUE,
    alignment = alignment,
    align_table = "center"
  )

  print(
    doc,
    target = file.path(
      supplement_dir,
      paste0("Table_S", number, ".docx")
    )
  )


  # ============================================================
  # 2. LATEX VERSION
  # ============================================================

  path <- file.path(
    supplement_dir,
    paste0("Table_S", number, ".tex")
  )

  subfields <- unique(x$Subfield)

  out <- c(
    "\\begin{tabularx}{\\textwidth}{",
    "  @{}",
    "  >{\\RaggedRight\\arraybackslash}p{0.13\\textwidth}",
    "  >{\\RaggedRight\\arraybackslash}p{0.145\\textwidth}",
    "  *{4}{>{\\centering\\arraybackslash}X}",
    "  @{}",
    "}",
    "\\toprule",
    "& & \\multicolumn{4}{c}{Genuine heterogeneity} \\\\",
    "\\cmidrule(lr){3-6}",
    "Subfield & Measure & 0\\% & 25\\% & 50\\% & 75\\% \\\\",
    "\\midrule"
  )


  # ============================================================
  # Table body
  # ============================================================

  for (s in subfields) {

    block <- x[
      x$Subfield == s,
      ,
      drop = FALSE
    ]

    esr <- block[
      block$Measure == "ESR_{0.05}^{sig}",
      ,
      drop = FALSE
    ]

    ma <- block[
      block$Measure == "No. of meta-analyses",
      ,
      drop = FALSE
    ]

    tests <- block[
      block$Measure == "No. of tests",
      ,
      drop = FALSE
    ]


    # Check that exactly one row of each type exists
    if (
      nrow(esr) != 1 ||
      nrow(ma) != 1 ||
      nrow(tests) != 1
    ) {
      stop(
        paste(
          "Unexpected table structure for subfield:",
          s
        )
      )
    }


    # Escape only the actual subfield text
    s_latex <- latex_escape(s)


    # ESR row
    row1 <- paste0(
      "\\multirow[t]{3}{0.13\\textwidth}{\\RaggedRight ",
      s_latex,
      "} & ",
      "$\\mathrm{ESR}_{0.05}^{\\mathrm{sig}}$",
      " & ", esr[[3]],
      " & ", esr[[4]],
      " & ", esr[[5]],
      " & ", esr[[6]],
      " \\\\"
    )


    # Number of meta-analyses
    row2 <- paste0(
      "& No. of meta-anal.",
      " & ", ma[[3]],
      " & ", ma[[4]],
      " & ", ma[[5]],
      " & ", ma[[6]],
      " \\\\"
    )


    # Number of tests
    row3 <- paste0(
      "& No. of tests",
      " & ", tests[[3]],
      " & ", tests[[4]],
      " & ", tests[[5]],
      " & ", tests[[6]],
      " \\\\"
    )


    out <- c(
      out,
      row1,
      row2,
      row3,
      "\\addlinespace[0.4em]"
    )
  }


  # Remove spacing after final subfield
  if (
    length(out) > 0 &&
    tail(out, 1) == "\\addlinespace[0.4em]"
  ) {
    out <- head(out, -1)
  }


  # Close tabularx
  out <- c(
    out,
    "\\bottomrule",
    "\\end{tabularx}"
  )


  # Write LaTeX file
  writeLines(
    out,
    path,
    useBytes = TRUE
  )


  # ============================================================
  # Useful confirmation
  # ============================================================

  message(
    "Written:\n",
    file.path(
      supplement_dir,
      paste0("Table_S", number, ".docx")
    ),
    "\n",
    path
  )

  invisible(path)
}



make_power_table <- function(dat, meta_average_multiplier = 0.5) {
  dat <- dat %>% mutate(
    sei = sqrt(vi), GE = meta_average_multiplier * GE,
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

## Tables S2-S3: sensitivity analyses at one-quarter and the full meta-average.
## Both use the full sample; only the assumed effect differs.
table_s2 <- make_power_table(base_half, 0.25)
table_s3 <- make_power_table(base_half, 1)

## Table S4: remove complete meta-analyses whose pooled effect is not significant.
significant_meta <- base_half %>% distinct(cID, sig_overall) %>%
  filter(!is.na(sig_overall), sig_overall < .05) %>% pull(cID)
table_s4 <- make_power_table(base_half %>% filter(cID %in% significant_meta))

## Table S9: share of meta-analyses with small-study effects by subfield.
small_study_effects <- load_small_study_effects("meta_0p5_heterogeneity_0")
table_s9_meta <- base_half %>%
  dplyr::select(-any_of(c("small_study_effect_pval", "sse_yn"))) %>%
  left_join(small_study_effects, by = "cID") %>%
  group_by(cID) %>% summarise(
    Subfield = first(subfd),
    small_study_effect = first(small_study_effect_pval) <= .05,
    .groups = "drop"
  )
table_s9_detail <- table_s9_meta %>%
  group_by(Subfield) %>%
  summarise(
    `No. of meta-analyses` = n_distinct(cID),
    `Small-study effects (%)` = 100 * mean(small_study_effect, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Subfield = factor(Subfield, subfield_levels)) %>%
  arrange(Subfield) %>%
  mutate(Subfield = as.character(Subfield))
table_s9 <- bind_rows(
  tibble(
    Subfield = "All meta-analyses",
    `No. of meta-analyses` = n_distinct(table_s9_meta$cID),
    `Small-study effects (%)` = 100 * mean(table_s9_meta$small_study_effect, na.rm = TRUE)
  ),
  table_s9_detail
) %>%
  mutate(`Small-study effects (%)` = sprintf("%.1f", `Small-study effects (%)`))

## Table S1: excess-significance results by subfield at half the meta-average
## and four degrees of genuine heterogeneity. Unlike Table 2, omit all p-value
## interval rows and retain only ESR_0.05^sig and the two sample-size rows.
subfield_esr_path <- here(
  "data", "derived_data", "Table_S6_subfield_inputs.rds"
)
if (!file.exists(subfield_esr_path)) {
  stop("Subfield ESR inputs are unavailable. Run create_tables_and_figures.R first.")
}
subfield_esr_results <- readRDS(subfield_esr_path)


table_s1_columns <- subfield_esr_results %>%
  filter(
    meta_average_multiplier == 0.5,
    measure %in% c(
      "ESR_{0.05}^{sig}",
      "No. of meta-analysis",
      "No. of tests"
    )
  ) %>%
  mutate(
    measure = dplyr::recode(
      measure,
      `No. of meta-analysis` = "No. of meta-analyses"
    ),
    value = if_else(
      confidence_interval == "0",
      estimate,
      paste(
        estimate,
        confidence_interval
      )
    ),
    heterogeneity_multiplier = factor(
      heterogeneity_multiplier,
      levels = c(
        0,
        0.25,
        0.5,
        0.75
      )
    )
  ) %>%
  dplyr::select(
    Subfield = subfield,
    Measure = measure,
    heterogeneity_multiplier,
    value
  ) %>%
  pivot_wider(
    names_from = heterogeneity_multiplier,
    values_from = value
  ) %>%
  arrange(
    factor(
      Subfield,
      levels = subfield_levels
    ),
    factor(
      Measure,
      levels = c(
        "ESR_{0.05}^{sig}",
        "No. of meta-analyses",
        "No. of tests"
      )
    )
  )


names(table_s1_columns) <- c(
  "Subfield",
  "Measure",
  "0%",
  "25%",
  "50%",
  "75%"
)

## Regression table formatting for Tables S6-S8. The negative-binomial fits
## consumed by Tables S7-S8 are created in run_exploratory_regressions.R.
significance_stars <- function(p_value) {
  ifelse(p_value < .01, "***", ifelse(p_value < .05, "**",
    ifelse(p_value < .10, "*", "")))
}
regression_labels <- c(
  "Intercept" = "(Intercept)",
  "Median power" = "med_perc",
  "Experimental research design? (yes)" = "design_mergedyes",
  "Followed reporting guidelines? (yes)" = "guidyes",
  "Protocol registered? (yes)" = "preryes",
  "Log number of independent studies" = "lognps",
  "Log journal impact factor" = "logjif",
  "Publication year" = "pyear"
)
format_model_table <- function(fit, include_adjusted_r2 = FALSE) {
  ## Force the fitted-model bundle before doing any validation. In particular,
  ## this makes the function safe to step through with debug()/debugonce()
  ## without repeatedly restarting evaluation of the lazy `fit` promise.
  models <- fit[["models"]]
  robust_tables <- fit[["robust"]]
  labels <- c(regression_labels,
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
    rep(c("0%", "50%"), each = 2), " genuine heterogeneity: Model ",
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
table_s6 <- format_model_table(nb_sensitivity_quarter)
table_s7 <- format_model_table(nb_sensitivity_full)

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
table_s5 <- format_model_table(ols_fit, TRUE)

## Figures S5-S8: separate continuous and categorical diagnostics for the two
## 50% genuine-heterogeneity negative-binomial specifications, matching the
## requested four-figure layout.
save_diagnostic_group <- function(model, model_number, kind, figure_number) {
  diagnostic_subfield_levels <- c(
    "Health, Toxicology and Mutagenesis", "Ecology",
    "Environmental Chemistry", "Environmental Engineering",
    "Management, Monitoring, Policy and Law",
    "Nature and Landscape Conservation", "Water Science and Technology"
  )
  dat <- final_nb_heterogeneity_0p5 %>%
    mutate(
      residual = residuals(model),
      fitted_value = fitted(model),
      subf = factor(subf, levels = diagnostic_subfield_levels, labels = 0:6)
    )
  vars <- if (kind == "continuous") c("fitted_value", "med_perc", "lognps", "logjif", "pyear") else c("design_merged", "guid", "prer", "subf")
  diagnostic_labels <- c(
    "fitted_value" = "Fitted value",
    setNames(
      sub(" \\(yes\\)$", "", names(regression_labels)),
      sub("yes$", "", unname(regression_labels))
    ),
    "subf" = "Subfield"
  )
  plots <- map(vars, function(v) {
    if (is.numeric(dat[[v]])) ggplot(dat, aes(.data[[v]], residual)) + geom_point(colour = "skyblue3", shape = 1) + geom_hline(yintercept = 0, colour = "red") + labs(x = diagnostic_labels[[v]]) + supplement_figure_theme()
    else ggplot(dat, aes(.data[[v]], residual)) + geom_boxplot() + geom_hline(yintercept = 0, colour = "red") + labs(x = diagnostic_labels[[v]]) + supplement_figure_theme()
  })
  grob <- arrangeGrob(grobs = plots, ncol = 2)
  save_supplement_plot(file.path(supplement_dir, paste0("Figure_S", figure_number)),
    11, ifelse(kind == "continuous", 10, 7), function() grid::grid.draw(grob))
}
save_diagnostic_group(nbMod3, 1, "continuous", 5)
save_diagnostic_group(nbMod3, 1, "categorical", 6)
save_diagnostic_group(nbMod4, 2, "continuous", 7)
save_diagnostic_group(nbMod4, 2, "categorical", 8)

## Table S8 is descriptive and therefore uses all observations rather than an
## estimator-specific outlier-screened sample.
table_s8_all_data <- load_multilevel_all_data()
table_s8 <- table_s8_all_data %>% distinct(cID, etype, subfd) %>%
  mutate(
    etype = dplyr::recode(
      etype,
      "cohen's d"              = "Cohen's $d$",
      "hedge's g"              = "Hedges's $g$",
      "correlation"            = "Correlation coefficient ($r$)",
      "fisher's z"             = "Fisher's $z$",
      "lnRR"                   = "Log-response ratio (lnRR)",
      "log-mean ratio"         = "Log-response ratio (lnRR)",
      "logOR"                  = "Log-odds ratio (logOR)",
      "logRR"                  = "Log-relative risk (logRR)",
      "logHR"                  = "Log-hazard ratio (logHR)",
      "mean"                   = "Raw mean or mean difference",
      "mean difference"        = "Raw mean or mean difference",
      "percentage change"      = "Percentage change",
      "excess risk"            = "Excess risk",
      "regression coefficient" = "Regression coefficient ($\\beta$)",
      "ratio"                  = "Log-response ratio (lnRR)"
    )
  ) %>%
 count(etype, subfd) %>%
  complete(etype, subfd = subfield_levels, fill = list(n = 0)) %>%
  pivot_wider(names_from = subfd, values_from = n) %>% rename(`Effect size` = etype) %>%
  arrange(`Effect size`) %>% mutate(No. = row_number(), .before = 1)

## Emit the supplementary tables in manuscript order. Keeping the writes in one
## place makes the filename-to-analysis mapping explicit and prevents a later
## analysis section from silently reusing a table number.
write_subfield_esr_table(table_s1_columns, number = 1)
write_word_table(table_s2, 2)
write_word_table(table_s3, 3)
write_word_table(table_s4, 4)
write_word_table(table_s5, 5)
write_word_table(table_s6, 6)
write_word_table(table_s7, 7)
write_word_table(table_s8, 8)
write_word_table(table_s9, 9)
