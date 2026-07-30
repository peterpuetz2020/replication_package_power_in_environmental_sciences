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

ensure_output_dirs <- function() {
  invisible(lapply(
    list(
      here("results", "main"),
      here("results", "main", "pet_peese_rstandard"),
      here("results", "main", "derived_data"),
      here("results", "robustness")
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
}

load_pet_peese_data <- function(setup_label_value = setup_label) {
  derived_path <- here("results", "main", "derived_data", paste0("pps_rstandard_raw_", setup_label_value, ".rds"))
  if (file.exists(derived_path)) {
    return(readRDS(derived_path))
  }

  pet_peese_files <- list.files(
    here("results", "main", "pet_peese_rstandard"),
    pattern = "\\.rds$",
    full.names = TRUE
  )

  if (length(pet_peese_files) == 0) {
    stop(
      "No PET-PEESE RDS files found in results/main/pet_peese_rstandard. ",
      "Run the PET-PEESE estimation block in scripts/main.R first."
    )
  }

  pet_peese_files %>%
    map_dfr(readRDS) %>%
    mutate(
      sei = sqrt(vi),
      sse_yn = ifelse(small_study_effect_pval <= 0.05, "yes", "no")
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

load_power_data <- function(setup_label_value = setup_label, meta_average_multiplier_value = meta_average_multiplier) {
  derived_path <- here("results", "main", "derived_data", paste0("pps_rstandard_power_", setup_label_value, ".rds"))
  if (file.exists(derived_path)) {
    return(readRDS(derived_path))
  }
  load_pet_peese_data(setup_label_value) %>%
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

get_counterfactual <- function(path, dat, grid, ci = FALSE, cluster = NULL, heterogeneity_multiplier_value = heterogeneity_multiplier) {
  if (file.exists(path)) {
    return(readRDS(path))
  }
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)
  result <- if (ci) {
    cf.ci.cluster(dat = dat, z.grid = grid, iters = n_iterations, cluster = cluster, heterogeneity_multiplier = heterogeneity_multiplier_value)
  } else {
    cf(dat = dat, z.grid = grid, heterogeneity_multiplier = heterogeneity_multiplier_value)
  }
  saveRDS(result, path)
  result
}

ensure_output_dirs()

## -------------------------------
## Table 1
## -------------------------------
write_table_1 <- function(meta_average_multiplier, heterogeneity_multiplier, setup_label, ...) {
pps_rstandard <- load_power_data(setup_label, meta_average_multiplier)
myDat <- split_meta_analyses(pps_rstandard, add_sape = TRUE)
mss <- vapply(myDat, function(x) length(x$sei), numeric(1))

idx <- tibble(
  i = seq_along(myDat),
  sdesn = vapply(myDat, function(x) x$sdesn[1], character(1)),
  sape = vapply(myDat, function(x) x$sape[1], numeric(1)),
  guide = vapply(myDat, function(x) x$guide[1], character(1)),
  prere = vapply(myDat, function(x) x$prere[1], character(1))
)

fill_desc <- function(i) c(length(i), sum(mss[i]), round(mean(mss[i])), round(median(mss[i])), min(mss[i]), round(quantile(mss[i], 0.25)), round(median(mss[i])), round(quantile(mss[i], 0.75)), max(mss[i]))
d.tab <- rbind(
  fill_desc(seq_along(myDat)),
  fill_desc(idx$i[idx$sdesn %in% c("observational", "mixed")]),
  fill_desc(idx$i[idx$sdesn == "experimental"]),
  fill_desc(idx$i[idx$sape > 0]),
  fill_desc(idx$i[idx$sape == 0]),
  fill_desc(idx$i[idx$guide == "yes"]),
  fill_desc(idx$i[idx$guide == "no"]),
  fill_desc(idx$i[idx$prere == "yes"]),
  fill_desc(idx$i[idx$prere == "no"])
)
rownames(d.tab) <- c("All meta-analyses", "Observational", "Experimental", "SAPE > 0", "SAPE = 0", "Yes", "No", "Yes", "No")
colnames(d.tab) <- c("M", "N", "Mean", "Median", "Min", "Q25", "Q50", "Q75", "Max")
write.csv(d.tab, here("results", "main", paste0("Descriptive_Table1_", setup_label, ".csv")))
}

analysis_setups %>% pwalk(write_table_1)

## -------------------------------
## Figure 1
## -------------------------------
write_figure_1 <- function(meta_average_multiplier, heterogeneity_multiplier, setup_label, ...) {
pps_rstandard <- load_pet_peese_data(setup_label)
grids <- make_grids()
my_dat <- pps_rstandard %>% mutate(GE = meta_average_multiplier * GE)
myDat <- split_meta_analyses(my_dat)
facz <- abs(my_dat$yi / sqrt(my_dat$vi))
z.orig <- count_intervals(facz, grids$z_grid_plot2)
p.orig.plot <- count_intervals(facz, grids$p_grid_plot[which(grids$p_grid_plot >= 0)])
z.plot <- get_counterfactual(here("results", "main", paste0("z_plot_pet_peese_rstandard_", setup_label, ".rds")), myDat, grids$z_grid_plot, heterogeneity_multiplier_value = heterogeneity_multiplier)
z.plot.ci <- get_counterfactual(here("results", "main", paste0("z_plot_ci_pet_peese_rstandard_", setup_label, ".rds")), myDat, grids$z_grid_plot, ci = TRUE, cluster = unique(pps_rstandard$cID), heterogeneity_multiplier_value = heterogeneity_multiplier)

xs <- as.vector(grids$z_grid_plot2[-length(grids$z_grid_plot2)] + (grids$z_grid_plot2[2] - grids$z_grid_plot2[1]) / 2)
N <- sum(p.orig.plot)
datFull <- as.data.frame(cbind(
  xs = xs,
  q025 = as.vector(apply(z.plot.ci[[1]], 2, quantile, na.rm = TRUE, probs = c(0.025))),
  q975 = as.vector(apply(z.plot.ci[[1]], 2, quantile, na.rm = TRUE, probs = c(0.975))),
  n.f = as.vector(z.orig / N),
  n.cf = as.vector(z.plot / N)
))

pdf(here("results", "main", paste0("zplot_pet_peese_rstandard_", setup_label, "_704_Fig1.pdf")), width = 10, height = 5)
ggplot(datFull) +
  geom_line(aes(xs, q025), color = "orange", lty = 3) +
  geom_line(aes(xs, n.cf), color = "orange", lty = 1) +
  geom_point(aes(xs, n.cf), shape = 20, fill = "orange", color = "orange", size = 1) +
  geom_line(aes(xs, q975), color = "orange", lty = 3) +
  geom_line(aes(xs, n.f), color = "blue", lty = 2) +
  geom_point(aes(xs, n.f), shape = 20, fill = "blue", color = "blue", size = 1) +
  coord_cartesian(xlim = c(0, 8)) +
  xlab("|z|-value") + ylab("Frequency") +
  geom_vline(xintercept = c(1.64, 1.96, 2.58), lty = 2, color = c(3, 2, 6), linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, 1.64, 1.96, 2.58, 4, 6, 8)) +
  theme(panel.background = element_rect(fill = "gray100"), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(linewidth = 0.5, color = "gray"))
dev.off()
}

analysis_setups %>% pwalk(write_figure_1)

## -------------------------------
## Table 2
## -------------------------------
write_table_2 <- function(meta_average_multiplier, heterogeneity_multiplier, setup_label, ...) {
  pps_rstandard <- load_pet_peese_data(setup_label)
  grids <- make_grids()
  my_dat <- pps_rstandard %>% mutate(GE = meta_average_multiplier * GE)
  myDat <- split_meta_analyses(my_dat)
  facz <- abs(my_dat$yi / sqrt(my_dat$vi))
  p.orig.tab <- count_intervals(facz, grids$p_grid_tab2)
  p.tab <- get_counterfactual(here("results", "main", paste0("p_tab_pps_rstandard_", setup_label, ".rds")), myDat, grids$p_grid_tab, heterogeneity_multiplier_value = heterogeneity_multiplier)
  p.tab.ci <- get_counterfactual(here("results", "main", paste0("p_tab_ci_pps_rstandard_", setup_label, ".rds")), myDat, grids$p_grid_tab, ci = TRUE, cluster = unique(pps_rstandard$cID), heterogeneity_multiplier_value = heterogeneity_multiplier)

  p.table <- matrix(ncol = 2, nrow = length(grids$p_grid_tab2) - 1)
  colnames(p.table) <- c("Difference", "0.95 CI")
  N <- sum(p.orig.tab)
  p.table[, 1] <- round((p.orig.tab - p.tab) / N, 3)
  q025 <- apply(matrix(p.orig.tab / N, nrow = nrow(p.tab.ci[[1]]), ncol = ncol(p.tab.ci[[1]]), byrow = TRUE) - p.tab.ci[[1]], 2, quantile, probs = c(0.025))
  q975 <- apply(matrix(p.orig.tab / N, nrow = nrow(p.tab.ci[[1]]), ncol = ncol(p.tab.ci[[1]]), byrow = TRUE) - p.tab.ci[[1]], 2, quantile, probs = c(0.975))
  p.table[, 2] <- paste("[", round(q025, 3), ", ", round(q975, 3), "]", sep = "")
  for (level in list(c(10, 13, 1, "all"), c(11, 13, 1, "all"), c(10, 13, 2, "sig"), c(11, 13, 3, "sig"))) {
    lo <- as.integer(level[[1]]); hi <- as.integer(level[[2]]); ci_idx <- as.integer(level[[3]]); denom <- if (level[[4]] == "all") N else sum(p.orig.tab[lo:hi])
    point <- round(sum((p.orig.tab - p.tab)[lo:hi] / denom), 3)
    bs <- apply(p.tab.ci[[ci_idx]][, lo:hi], 1, sum)
    q <- round(quantile(sum((p.orig.tab / denom)[lo:hi]) - bs, probs = c(0.025, 0.975)), 3)
    p.table <- rbind(p.table, c(point, paste("[", q[1], ", ", q[2], "]", sep = "")))
  }
  p.table <- rbind(p.table, c(length(myDat), 0), c(N, 0))
  rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7", "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3", "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01", "p < 0.001", "ESR_{0.1}^{all}", "ESR_{0.05}^{all}", "ESR_{0.1}^{sig}", "ESR_{0.05}^{sig}", "No. of meta-analysis", "No. of tests")
  write.csv(p.table, here("results", "main", paste0("p.table.ci_pet_peese_rstandard_", setup_label, "_Table 2.csv")))
}

analysis_setups %>% pwalk(write_table_2)

## -------------------------------
## Table 3
## -------------------------------
write_table_3 <- function(meta_average_multiplier, heterogeneity_multiplier, setup_label, ...) {
pps_rstandard <- load_pet_peese_data(setup_label) %>% add_power_variables(meta_average_multiplier)
subf_desc <- pps_rstandard %>% group_by(subfd) %>% summarise(M = length(unique(cID)), N = length(power), median = round(median(power), 2), mean = round(mean(power), 2), Q25 = round(quantile(power, 0.25), 2), Q75 = round(quantile(power, 0.75), 2), sape = round(length(which(power >= 0.8)) / length(power), 2), .groups = "drop")
med_med <- pps_rstandard %>% group_by(cID) %>% summarise(metaID = metaID[1], subfd = subfd[1], median = median(power), .groups = "drop")
med_med_subf <- med_med %>% group_by(subfd) %>% summarise(mmedian = round(median(median), 2), .groups = "drop")
power.tab3 <- matrix(nrow = 8, ncol = 8)
power.tab3[1, ] <- c(nrow(med_med), nrow(pps_rstandard), round(summary(med_med$median), 2)[3], round(summary(pps_rstandard$power), 2)[3], round(summary(pps_rstandard$power), 2)[4], round(summary(pps_rstandard$power), 2)[2], round(summary(pps_rstandard$power), 2)[5], 0.18)
for (i in seq_len(nrow(subf_desc))) power.tab3[i + 1, ] <- c(subf_desc$M[i], subf_desc$N[i], med_med_subf$mmedian[i], subf_desc$median[i], subf_desc$mean[i], subf_desc$Q25[i], subf_desc$Q75[i], subf_desc$sape[i])
colnames(power.tab3) <- c("M", "N", "mmedian", "median", "mean", "Q25", "Q75", "SAPE")
rownames(power.tab3) <- c("All meta-analyses", "Ecology", "Environmental Chemistry", "Environmental Engineering", "Health, Toxicology and Mutagenesis", "Management, Monitoring, Policy and Law", "Nature and Landscape Conservation", "Water Science and Technology")
write.csv(power.tab3, here("results", "main", paste0("Power_Table3_", setup_label, ".csv")))
}

analysis_setups %>% pwalk(write_table_3)

## -------------------------------
## Figure 2
## -------------------------------
write_figure_2 <- function(meta_average_multiplier, heterogeneity_multiplier, setup_label, ...) {
pps_rstandard <- load_power_data(setup_label, meta_average_multiplier)
pps_rstandard_median <- pps_rstandard %>% group_by(cID) %>% summarise(metaID = metaID[1], median = median(power), sape = length(which(power >= 0.8)) / length(power), nips = length(unique(sID)), esty = unique(etype), guid = unique(guide), prer = unique(prere), subf = unique(subfd), sdes = unique(sdesn), .groups = "drop") %>% mutate(yn80 = ifelse(median >= 0.8, "yes", "no"), median100 = round(100 * median, 2), sape100 = round(100 * sape, 2))
write.xlsx(pps_rstandard_median, here("results", "main", paste0("median_power_pps_rstandard_", setup_label, "_704.xlsx")), overwrite = TRUE)
med_pwr <- pps_rstandard_median %>% ggplot(aes(x = median100, fill = as.factor(yn80))) + geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins = 30, alpha = I(0.6), linewidth = 0.1) + scale_fill_manual(values = c("brown2", "skyblue2")) + xlab("Median statistical power of primary estimates per meta-analysis") + ylab("Percentage") + ggtitle("(a)") + scale_x_continuous(breaks = breaks_width(20), labels = label_percent(scale = 1), expand = c(0, 0.5)) + scale_y_continuous(labels = label_percent(scale = 1), expand = c(0, 0.5)) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.5, color = "gray"))
sape <- pps_rstandard_median %>% ggplot(aes(x = sape100)) + geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins = 30, alpha = I(0.6), linewidth = 0.1, fill = "skyblue2") + xlab("Share of adequately powered primary estimates per meta-analysis") + ylab("Percentage") + ggtitle("(b)") + scale_x_continuous(breaks = breaks_width(20), labels = label_percent(scale = 1), expand = c(0, 0.5)) + scale_y_continuous(labels = label_percent(scale = 1), expand = c(0, 0.5)) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), axis.line = element_line(linewidth = 0.5, color = "gray"))
pdf(here("results", "main", paste0("pps_rstandard_", setup_label, "_704_Fig2.pdf")), width = 10, height = 4)
grid.arrange(med_pwr, sape, ncol = 2)
dev.off()
}

analysis_setups %>% pwalk(write_figure_2)
