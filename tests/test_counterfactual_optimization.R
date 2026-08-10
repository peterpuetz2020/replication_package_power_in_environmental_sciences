library(foreach)
source(file.path("scripts", "functions.R"))
registerDoSEQ()

legacy_probabilities <- function(dat, grid, heterogeneity_multiplier) {
  lapply(dat, function(meta_data) {
    result <- matrix(NA_real_, nrow(meta_data), length(grid) - 1L)
    for (j in seq_len(nrow(meta_data))) {
      for (a in seq_len(length(grid) - 1L)) {
        mean_z <- meta_data$GE[j] / sqrt(meta_data$vi[j])
        sd_z <- sqrt(meta_data$vi[j] + heterogeneity_multiplier *
          meta_data$tau2[j]) / sqrt(meta_data$vi[j])
        result[j, a] <- pnorm(grid[a + 1L], mean_z, sd_z) -
          pnorm(grid[a], mean_z, sd_z)
      }
    }
    result
  })
}

legacy_fold <- function(values, grid) {
  half <- length(values) / 2L
  result <- numeric(half)
  for (index in 0:(half - 1L)) {
    result[half - index] <- values[1L + index] +
      values[length(values) - index]
  }
  names(result) <- grid[grid > 0]
  result
}

legacy_cf <- function(dat, grid, heterogeneity_multiplier) {
  probabilities <- legacy_probabilities(dat, grid, heterogeneity_multiplier)
  legacy_fold(colSums(do.call(rbind, probabilities)), grid)
}

legacy_cf_ci <- function(dat, grid, iters, cluster, heterogeneity_multiplier) {
  probabilities <- legacy_probabilities(dat, grid, heterogeneity_multiplier)
  output <- lapply(seq_len(iters), function(iteration) {
    set.seed(iteration + 12)
    sampled <- sample(unique(cluster), replace = TRUE)
    boots <- unlist(lapply(sampled, function(value) which(cluster == value)))
    numerator <- colSums(do.call(rbind, probabilities[boots]))
    observed_z <- lapply(dat[boots], function(meta_data) {
      abs(meta_data$yi / sqrt(meta_data$vi))
    })
    denominators <- c(
      sum(vapply(dat[boots], nrow, integer(1))),
      sum(vapply(observed_z, function(z) sum(z > abs(qnorm(0.1 / 2))), integer(1))),
      sum(vapply(observed_z, function(z) sum(z > abs(qnorm(0.05 / 2))), integer(1)))
    )
    lapply(denominators, function(denominator) {
      legacy_fold(numerator / denominator, grid)
    })
  })
  lapply(seq_len(3L), function(index) {
    do.call(rbind, lapply(output, `[[`, index))
  })
}

set.seed(2026)
dat <- lapply(c(3L, 4L, 2L), function(n) {
  vi <- runif(n, 0.02, 0.2)
  data.frame(
    GE = rnorm(n, 0.1, 0.05),
    yi = rnorm(n, 0.15, sqrt(vi)),
    vi = vi,
    tau2 = runif(n, 0, 0.08)
  )
})
grid <- seq(-3, 3, by = 0.5)
heterogeneity_multiplier <- 0.25
cluster <- c("shared", "shared", "other")
iters <- 12L

expected_point <- legacy_cf(dat, grid, heterogeneity_multiplier)
actual_point <- cf(dat, grid, heterogeneity_multiplier)
stopifnot(isTRUE(all.equal(actual_point, expected_point, tolerance = 1e-14)))

expected_ci <- legacy_cf_ci(
  dat, grid, iters, cluster, heterogeneity_multiplier
)
actual_ci <- cf.ci.cluster(
  dat, grid, iters, cluster, heterogeneity_multiplier
)
stopifnot(isTRUE(all.equal(actual_ci, expected_ci, tolerance = 1e-14)))

cache_path <- tempfile(fileext = ".rds")
cache_key <- counterfactual_cache_key(
  dat, grid, heterogeneity_multiplier, TRUE, cluster, iters
)
write_counterfactual_cache(actual_ci, cache_path, cache_key)
stopifnot(identical(read_counterfactual_cache(cache_path, cache_key), actual_ci))
different_iteration_key <- counterfactual_cache_key(
  dat, grid, heterogeneity_multiplier, TRUE, cluster, iters + 1L
)
stopifnot(is.null(read_counterfactual_cache(
  cache_path, different_iteration_key
)))
replacement <- lapply(actual_ci, function(result) result[1, , drop = FALSE])
write_counterfactual_cache(replacement, cache_path, different_iteration_key)
stopifnot(identical(
  read_counterfactual_cache(cache_path, different_iteration_key), replacement
))
unlink(cache_path)

components <- counterfactual_components(dat, grid, heterogeneity_multiplier)
stopifnot(isTRUE(all.equal(
  cf(dat, grid, components = components), expected_point, tolerance = 1e-14
)))

dat_with_missing_estimate <- dat
dat_with_missing_estimate[[1]]$GE[1] <- NA_real_
filtered_dat <- lapply(
  dat_with_missing_estimate,
  filter_counterfactual_data,
  heterogeneity_multiplier = heterogeneity_multiplier,
  context = "test"
)
stopifnot(nrow(filtered_dat[[1]]) == nrow(dat[[1]]) - 1L)
stopifnot(all(vapply(
  filtered_dat,
  function(meta_data) all(counterfactual_valid_rows(
    meta_data, heterogeneity_multiplier
  )),
  logical(1)
)))
stopifnot(all(is.finite(cf(filtered_dat, grid, heterogeneity_multiplier))))

invalid_error <- try(
  cf(dat_with_missing_estimate, grid, heterogeneity_multiplier),
  silent = TRUE
)
stopifnot(inherits(invalid_error, "try-error"))
cat("Counterfactual optimization equivalence checks passed.\n")
