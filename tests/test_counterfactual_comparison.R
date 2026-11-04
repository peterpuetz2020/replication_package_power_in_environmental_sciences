library(foreach)
source("scripts/functions.R")
registerDoSEQ()

## Analytic normal probabilities give an independent reference for every bin.
grid <- counterfactual_p_grid()
dat <- data.frame(
  yi = c(.1, -.2, .3, -.4, .5, -.6), vi = 1, tau2 = 0, GE = 1,
  cID = c("a", "a", "b", "b", "c", "c"),
  metaID = c("one", "one", "two", "two", "two", "two")
)
half <- split(transform(dat, GE = .5), dat$cID)
full <- split(dat, dat$cID)
clusters <- counterfactual_meta_clusters(full)
half_draws <- cf.ci.cluster(half, grid, 40L, clusters, 0)[[1]]
full_draws <- cf.ci.cluster(full, grid, 40L, clusters, 0)[[1]]
result <- compare_counterfactual_bins(cf(half, grid, 0), cf(full, grid, 0),
                                      nrow(dat), half_draws, full_draws)
cutoffs <- c(1, .9, .8, .7, .6, .5, .4, .3, .2, .1, .05, .01, .001, 0)
analytic_bins <- function(mu) {
  z <- qnorm(cutoffs / 2, lower.tail = FALSE)
  tail_mass <- pnorm(-z, mean = mu) + pnorm(z, mean = mu, lower.tail = FALSE)
  -diff(tail_mass)
}
stopifnot(
  nrow(result) == 13L,
  identical(result$p_value_interval, p_value_bin_labels()),
  max(abs(result$reference_value - analytic_bins(.5))) < 1e-12,
  max(abs(result$counterfactual_value - analytic_bins(1))) < 1e-12,
  abs(sum(result$difference)) < 1e-12,
  result$difference[1] > 0, result$difference[13] < 0,
  max(abs(result$ci_lower - result$difference)) < 1e-12,
  max(abs(result$ci_upper - result$difference)) < 1e-12
)

## Perfectly shared bootstrap noise must cancel in a paired CI, even when
## each individual counterfactual has substantial uncertainty.
draws <- matrix(0, 4, 13)
draws[, 1] <- c(.2, .3, .4, .5)
draws[, 2] <- 1 - draws[, 1]
shift <- c(.02, -.02, rep(0, 11))
paired <- compare_counterfactual_bins(
  c(35, 65, rep(0, 11)), c(33, 67, rep(0, 11)), 100,
  draws, sweep(draws, 2, shift)
)
stopifnot(max(abs(paired$ci_lower - shift)) < 1e-12,
          max(abs(paired$ci_upper - shift)) < 1e-12)
invalid <- try(compare_counterfactual_bins(c(1, 2), c(1, 2), 3,
                                          half_draws, full_draws), silent = TRUE)
stopifnot(inherits(invalid, "try-error"))

## A simulated full-effect reference replaces Table 2's empirical reference.
tab <- compare_counterfactual_table(
  cf(full, grid, 0), cf(half, grid, 0), nrow(dat), full_draws, half_draws, 3L
)
ref_sig <- pnorm(-qnorm(.975), mean = 1) + pnorm(qnorm(.975), mean = 1, lower.tail = FALSE)
half_sig <- pnorm(-qnorm(.975), mean = .5) + pnorm(qnorm(.975), mean = .5, lower.tail = FALSE)
stopifnot(identical(tab$measure, table_2_display_measures()),
          abs(tab$estimate[14] - (ref_sig - half_sig)) < 1e-12,
          abs(tab$estimate[15] - (ref_sig - half_sig) / ref_sig) < 1e-12,
          identical(tail(tab$estimate, 2), c(3, nrow(dat))))
zero <- compare_counterfactual_table(
  cf(half, grid, 0), cf(half, grid, 0), nrow(dat), half_draws, half_draws, 3L
)
stopifnot(all(zero$estimate[1:15] == 0), all(zero$ci_lower[1:15] == 0),
          all(zero$ci_upper[1:15] == 0))

## Significant-result ESR divides by each reference draw's significant mass.
reference <- counterfactual <- matrix(0, 4, 13)
reference[, 13] <- c(.2, .4, .6, .8)
counterfactual[, 13] <- c(.05, .15, .35, .65)
reference[, 1] <- 1 - reference[, 13]
counterfactual[, 1] <- 1 - counterfactual[, 13]
ratio <- compare_counterfactual_table(
  c(50, rep(0, 11), 50), c(75, rep(0, 11), 25), 100,
  reference, counterfactual, 3L
)
stopifnot(abs(ratio$estimate[15] - .5) < 1e-12,
          abs(ratio$ci_lower[15] - .2046875) < 1e-12,
          abs(ratio$ci_upper[15] - .740625) < 1e-12)
cat("PASS: analytic bins, Table 2 row order, ESR denominators and paired intervals\n")
