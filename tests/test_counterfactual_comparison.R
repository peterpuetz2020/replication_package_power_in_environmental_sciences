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
  max(abs(result$probability_0p5 - analytic_bins(.5))) < 1e-12,
  max(abs(result$probability_1 - analytic_bins(1))) < 1e-12,
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
cat("PASS: analytic bin probabilities, tail direction, normalization and paired intervals\n")
