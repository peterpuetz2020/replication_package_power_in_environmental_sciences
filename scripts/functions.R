
## ------------
## functions.R
## ------------

counterfactual_cache_version <- 3L

counterfactual_p_grid <- function() {
  cutoffs <- c(.001, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9)
  c(-Inf, qnorm(cutoffs / 2), 0, qnorm(rev(cutoffs) / 2, lower.tail = FALSE), Inf)
}

## Shared row order for Table 2 and counterfactual-distribution comparisons.
p_value_bin_labels <- function() {
  c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
    "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
    "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05",
    "0.001 < p < 0.01", "p < 0.001")
}

table_2_display_measures <- function() {
  c(p_value_bin_labels(), "ESR_{0.05}^{all}", "ESR_{0.05}^{sig}",
    "No. of meta-analysis", "No. of tests")
}

table_2_column_label <- function(meta_multiplier, heterogeneity) {
  paste0("(", meta_multiplier, " x meta-average; ", heterogeneity * 100,
         "% genuine heterogeneity)\nDifference [95% CI]")
}

## Both counterfactuals must use the same sample, article clusters and bootstrap
## iteration IDs. Subtract matched draws, rather than separate CI endpoints.
compare_counterfactual_bins <- function(reference_counts, counterfactual_counts, n_effects,
                                        reference_draws, counterfactual_draws) {
  n_bins <- length(p_value_bin_labels())
  if (length(n_effects) != 1L || !is.finite(n_effects) || n_effects <= 0 ||
      length(reference_counts) != n_bins || length(counterfactual_counts) != n_bins ||
      !is.matrix(reference_draws) || !is.matrix(counterfactual_draws) ||
      !identical(dim(reference_draws), dim(counterfactual_draws)) ||
      ncol(reference_draws) != n_bins || nrow(reference_draws) < 1L) {
    stop("Counterfactual comparison requires matching counts and paired bootstrap draws.")
  }
  if (any(!is.finite(c(reference_counts, counterfactual_counts, reference_draws, counterfactual_draws))) ||
      abs(sum(reference_counts) - n_effects) > 1e-8 * n_effects ||
      abs(sum(counterfactual_counts) - n_effects) > 1e-8 * n_effects ||
      any(abs(rowSums(reference_draws) - 1) > 1e-8) ||
      any(abs(rowSums(counterfactual_draws) - 1) > 1e-8)) {
    stop("Counterfactual probabilities must be finite and sum to one.")
  }
  difference_draws <- reference_draws - counterfactual_draws
  intervals <- apply(difference_draws, 2, quantile,
                     probs = c(0.025, 0.975), na.rm = TRUE)
  data.frame(
    p_value_interval = p_value_bin_labels(),
    reference_value = as.numeric(reference_counts) / n_effects,
    counterfactual_value = as.numeric(counterfactual_counts) / n_effects,
    difference = as.numeric(reference_counts - counterfactual_counts) / n_effects,
    ci_lower = intervals[1, ], ci_upper = intervals[2, ],
    row.names = NULL
  )
}

## Table 2's empirical-reference formula, with a simulated reference instead.
## Significant-result ESR uses the reference distribution's significant mass
## as denominator, including in each paired bootstrap replication.
compare_counterfactual_table <- function(reference_counts, counterfactual_counts,
                                         n_effects, reference_draws,
                                         counterfactual_draws, n_meta_analyses) {
  bins <- compare_counterfactual_bins(
    reference_counts, counterfactual_counts, n_effects,
    reference_draws, counterfactual_draws
  )
  result <- data.frame(
    measure = bins$p_value_interval,
    reference_value = bins$reference_value,
    counterfactual_value = bins$counterfactual_value,
    estimate = bins$difference,
    ci_lower = bins$ci_lower, ci_upper = bins$ci_upper
  )
  for (level in list(c(10L, "0.1"), c(11L, "0.05"))) {
    idx <- as.integer(level[[1]]):13L
    reference_mass <- sum(reference_counts[idx]) / n_effects
    counterfactual_mass <- sum(counterfactual_counts[idx]) / n_effects
    reference_boot <- rowSums(reference_draws[, idx, drop = FALSE])
    counterfactual_boot <- rowSums(counterfactual_draws[, idx, drop = FALSE])
    if (reference_mass <= 0 || any(reference_boot <= 0)) {
      stop("Significant-result ESR requires positive reference probability.")
    }
    for (denominator in c("all", "sig")) {
      scale <- if (denominator == "all") 1 else reference_mass
      boot_scale <- if (denominator == "all") 1 else reference_boot
      interval <- quantile((reference_boot - counterfactual_boot) / boot_scale,
                           probs = c(.025, .975), na.rm = TRUE)
      result <- rbind(result, data.frame(
        measure = paste0("ESR_{", level[[2]], "}^{", denominator, "}"),
        reference_value = reference_mass / scale,
        counterfactual_value = counterfactual_mass / scale,
        estimate = (reference_mass - counterfactual_mass) / scale,
        ci_lower = interval[[1]], ci_upper = interval[[2]]
      ))
    }
  }
  for (sample_size in c("No. of meta-analysis", "No. of tests")) {
    value <- if (sample_size == "No. of tests") n_effects else n_meta_analyses
    result <- rbind(result, data.frame(
      measure = sample_size, reference_value = value, counterfactual_value = value,
      estimate = value, ci_lower = NA_real_, ci_upper = NA_real_
    ))
  }
  result <- result[match(table_2_display_measures(), result$measure), ]
  rownames(result) <- NULL
  result
}

## Return one article-level cluster identifier for every meta-analysis in the
## split counterfactual data. Several meta-analyses (cID) can come from the same
## article (metaID), so bootstrapping cID would incorrectly treat them as
## independent.
counterfactual_meta_clusters <- function(dat) {
  vapply(dat, function(meta_data) {
    meta_ids <- unique(meta_data$metaID)
    if (length(meta_ids) != 1L || is.na(meta_ids)) {
      stop("Each meta-analysis must have exactly one non-missing metaID.")
    }
    as.character(meta_ids)
  }, character(1))
}

counterfactual_valid_rows <- function(dat, heterogeneity_multiplier = 0.25) {
  variance <- if (heterogeneity_multiplier == 0) {
    dat$vi
  } else {
    dat$vi + heterogeneity_multiplier * dat$tau2
  }
  is.finite(dat$GE) & is.finite(dat$yi) & is.finite(dat$vi) &
    dat$vi > 0 & is.finite(variance) & variance > 0
}

filter_counterfactual_data <- function(dat, heterogeneity_multiplier = 0.25,
                                       context = "counterfactual analysis") {
  valid <- counterfactual_valid_rows(dat, heterogeneity_multiplier)
  if (any(!valid)) {
    message(
      "Excluding ", sum(!valid), " of ", length(valid),
      " effects with non-finite counterfactual parameters from ", context, "."
    )
  }
  filtered <- dat[valid, , drop = FALSE]
  if (nrow(filtered) == 0L) {
    stop("No effects with finite counterfactual parameters remain for ", context, ".")
  }
  filtered
}

counterfactual_cache_key <- function(dat, z.grid, heterogeneity_multiplier,
                                     ci = FALSE, cluster = NULL, iters = NULL) {
  cache_input <- tempfile("counterfactual-cache-")
  on.exit(unlink(cache_input), add = TRUE)
  saveRDS(
    list(dat = dat, z.grid = z.grid,
         heterogeneity_multiplier = heterogeneity_multiplier,
         ci = ci, cluster = cluster, iters = iters),
    cache_input, version = 2
  )
  unname(tools::md5sum(cache_input))
}

read_counterfactual_cache <- function(path, key) {
  if (!file.exists(path)) return(NULL)
  cached <- readRDS(path)
  if (!is.list(cached) || !identical(cached$version, counterfactual_cache_version) ||
      !identical(cached$key, key) || is.null(cached$result)) {
    message("Ignoring stale counterfactual cache: ", basename(path))
    return(NULL)
  }
  cached$result
}

write_counterfactual_cache <- function(result, path, key) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(paste0(".", basename(path), "-"), dirname(path))
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(
    list(version = counterfactual_cache_version, key = key, result = result),
    temporary_path
  )
  if (!file.rename(temporary_path, path)) {
    if (!file.copy(temporary_path, path, overwrite = TRUE)) {
      stop("Could not replace counterfactual cache: ", path)
    }
    unlink(temporary_path)
  }
  invisible(path)
}

## Vectorized counterfactual calculations. Probability masses depend on an
## effect size and setup, but not on a bootstrap draw, so calculate and sum them
## once per meta-analysis. Bootstrap iterations then only resample rows of this
## much smaller contribution matrix.
counterfactual_components <- function(dat, z.grid,
                                      heterogeneity_multiplier = 0.25) {
  if (length(z.grid) < 3L || (length(z.grid) - 1L) %% 2L != 0L) {
    stop("z.grid must define an even number of intervals.")
  }

  contribution <- lapply(dat, function(meta_data) {
    if (!all(counterfactual_valid_rows(meta_data, heterogeneity_multiplier))) {
      stop(
        "Counterfactual data contain non-finite parameters. ",
        "Call filter_counterfactual_data() before splitting the data."
      )
    }
    mean_z <- meta_data$GE / sqrt(meta_data$vi)
    sd_z <- sqrt(
      meta_data$vi + heterogeneity_multiplier * meta_data$tau2
    ) / sqrt(meta_data$vi)
    cdf <- vapply(
      z.grid,
      function(boundary) pnorm(boundary, mean = mean_z, sd = sd_z),
      numeric(length(mean_z))
    )
    interval_probability <- cdf[, -1L, drop = FALSE] -
      cdf[, -ncol(cdf), drop = FALSE]
    colSums(interval_probability)
  })

  observed_z <- lapply(dat, function(meta_data) {
    abs(meta_data$yi / sqrt(meta_data$vi))
  })
  list(
    contribution = do.call(rbind, contribution),
    n_effects = vapply(dat, function(meta_data) nrow(meta_data), integer(1)),
    n_sig1 = vapply(
      observed_z, function(z) sum(z > abs(qnorm(0.1 / 2))), integer(1)
    ),
    n_sig5 = vapply(
      observed_z, function(z) sum(z > abs(qnorm(0.05 / 2))), integer(1)
    )
  )
}

fold_counterfactual_intervals <- function(freqs, z.grid) {
  freqs <- as.matrix(freqs)
  n_intervals <- ncol(freqs)
  half <- n_intervals / 2L
  folded <- freqs[, seq_len(half), drop = FALSE] +
    freqs[, n_intervals:as.integer(half + 1L), drop = FALSE]
  folded <- folded[, half:1L, drop = FALSE]
  colnames(folded) <- z.grid[z.grid > 0]
  folded
}

## Function for counterfactual z- or p-values.
cf <- function(dat, z.grid, heterogeneity_multiplier = 0.25,
               components = NULL) {
  if (is.null(components)) {
    components <- counterfactual_components(
      dat, z.grid, heterogeneity_multiplier
    )
  }
  interval_totals <- matrix(
    colSums(components$contribution), nrow = 1L
  )
  drop(fold_counterfactual_intervals(interval_totals, z.grid))
}

## Confidence intervals of counterfactual z-values (cluster bootstrap).
cf.ci.cluster <- function(dat, z.grid, iters, cluster,
                          heterogeneity_multiplier = 0.25,
                          iteration_ids = seq_len(iters), components = NULL) {
  if (length(iteration_ids) != iters) {
    stop("iteration_ids must contain exactly iters values.")
  }
  if (length(cluster) != length(dat)) {
    stop("cluster must have one value for each meta-analysis in dat.")
  }
  if (is.null(components)) {
    components <- counterfactual_components(
      dat, z.grid, heterogeneity_multiplier
    )
  }

  unique_clusters <- unique(cluster)
  meta_cluster_index <- match(cluster, unique_clusters)
  ## Each draw has its own seed, so prescheduling changes only communication
  ## overhead, not bootstrap samples or result ordering.
  bootstrap_results <- foreach(
    dd = iteration_ids, .options.snow = list(preschedule = TRUE)
  ) %dopar% {
    set.seed(dd + 12)
    sampled_clusters <- sample(unique_clusters, replace = TRUE)
    cluster_counts <- tabulate(
      match(sampled_clusters, unique_clusters),
      nbins = length(unique_clusters)
    )
    meta_counts <- cluster_counts[meta_cluster_index]
    numerator <- drop(meta_counts %*% components$contribution)
    list(
      numerator / sum(meta_counts * components$n_effects),
      numerator / sum(meta_counts * components$n_sig1),
      numerator / sum(meta_counts * components$n_sig5)
    )
  }

  bind_result <- function(index) {
    do.call(rbind, lapply(bootstrap_results, `[[`, index))
  }
  lapply(
    list(bind_result(1L), bind_result(2L), bind_result(3L)),
    fold_counterfactual_intervals,
    z.grid = z.grid
  )
}


## This is a separate function for the robustness check with increased standard errors. In the original function, we bootstrap shares,
#	this means we divide in each bs iterations by the number of all p-values or all significant p-values in the given iteration.
#	For the case with increased standard errors, it is important to divide by all significant p-values without inflating the standard error.

cf.ci.cluster.se <- function(dat, z.grid, iters, cluster) {

  # This is the loop over iterations indexed by dd
  freqs <- foreach(dd = 1:iters) %dopar% {

    # Within each foreach loop the random numbers are set in a reproducible way.
    set.seed(dd+12)

    #resampling of studies which may contain multiple meta-analyses
    bs.cluster <- sample(unique(cluster), replace=TRUE)

    boots <- NULL
    for (ii in 1:length(bs.cluster)) {
      boots <- c(boots,	which(cluster == bs.cluster[ii]))
    }

    # This is the loop over meta-analyses indexed by i
    freq.per.iter <- list()
    sig.1 <- NULL
    sig.5 <- NULL
    k <- 1

    for (i in boots) {

      freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$vi))

      # This is the loop over the observations per meta-analyses indexed by j
      for (j in 1:length(dat[[i]]$vi)) {

        # This is for each distribution/observation the loop over the grid indexed by a
        for (a in 1:(length(z.grid)-1)) {

          # This calculates the density between the intervals of the z.grid
          # freq is a matrix with the z.grid intervals as columns and
          # the observations/distributions per meta-analysis as rows
          freq[j,a] <- (pnorm(z.grid[a+1], mean = dat[[i]]$GE[j] / dat[[i]]$sei[j], sd = 1) -
                          pnorm(z.grid[a],   mean = dat[[i]]$GE[j] / dat[[i]]$sei[j], sd = 1))
          #freq[j,a] <- (pnorm(z.grid[a+1], mean = dat[[i]]$GE[j] / sqrt(dat[[i]]$vi[j]+dat[[i]]$tau2[j]), sd = 1) -
          #              pnorm(z.grid[a],   mean = dat[[i]]$GE[j] / sqrt(dat[[i]]$vi[j]+dat[[i]]$tau2[j]), sd = 1))
        }

        sig.1.temp <- length(which(abs(dat[[i]]$yi / sqrt(dat[[i]]$vi_orig)) > abs(qnorm(0.10/2))))
        sig.5.temp <- length(which(abs(dat[[i]]$yi / sqrt(dat[[i]]$vi_orig)) > abs(qnorm(0.05/2))))
      }

      sig.1 <- c(sig.1, sig.1.temp)
      sig.5 <- c(sig.5, sig.5.temp)
      freq.per.iter[[k]] <- freq #for each MA
      k <- k+1
    }
    # freq.per.iter is a list with each entry being the freq matrix (see above) for one
    # meta-analysis
    freq.per.iter <- do.call(rbind.data.frame, freq.per.iter)
    # sums the columns for each interval and divides by the number of observations
    freq.per.iter.all  <- apply(freq.per.iter, 2, sum) / dim(freq.per.iter)[1]
    freq.per.iter.sig1 <- apply(freq.per.iter, 2, sum) / sum(sig.1)
    freq.per.iter.sig5 <- apply(freq.per.iter, 2, sum) / sum(sig.5)

    names(freq.per.iter.all)  <- z.grid[-1] #some of the names change in an odd way.
    names(freq.per.iter.sig1) <- z.grid[-1]
    names(freq.per.iter.sig5) <- z.grid[-1]

    # freq.per.iter is one row with the counterfactual frequencies per grid
    return(list(freq.per.iter.all, freq.per.iter.sig1, freq.per.iter.sig5))
  }
  # freqs is a list with each element being the outcome of one iteration (freq.per.iter - see above)
  freqs.all  <- NULL
  freqs.sig1 <- NULL
  freqs.sig5 <- NULL

  for (jj in 1:iters) {
    freqs.all[jj]  <- freqs[[jj]][1]
    freqs.sig1[jj] <- freqs[[jj]][2]
    freqs.sig5[jj] <- freqs[[jj]][3]
  }

  freqs.all  <- do.call(rbind.data.frame, freqs.all)
  freqs.sig1 <- do.call(rbind.data.frame, freqs.sig1)
  freqs.sig5 <- do.call(rbind.data.frame, freqs.sig5)
  names(freqs.all)  <- z.grid[-1] #some of the names change in an odd way.
  names(freqs.sig1) <- z.grid[-1]
  names(freqs.sig5) <- z.grid[-1]

  # This takes the sum of the first and last column and the second and the second last
  # and so on. Therefore, it is essential that the z.grid is symmetric around zero.
  abs.freqs.all  <- matrix(ncol=(dim(freqs.all)[2] / 2),  nrow=iters)
  abs.freqs.sig1 <- matrix(ncol=(dim(freqs.sig1)[2] / 2), nrow=iters)
  abs.freqs.sig5 <- matrix(ncol=(dim(freqs.sig5)[2] / 2), nrow=iters)

  for (k in 0:((dim(freqs.all)[2] / 2) -1)) {
    abs.freqs.all[,(dim(freqs.all)[2] / 2)-k]   <- apply(matrix(c(freqs.all[,1+k],  freqs.all[,dim(freqs.all)[2]-k]),   ncol=2), 1, sum)
    abs.freqs.sig1[,(dim(freqs.sig1)[2] / 2)-k] <- apply(matrix(c(freqs.sig1[,1+k], freqs.sig1[,dim(freqs.sig1)[2]-k]), ncol=2), 1, sum)
    abs.freqs.sig5[,(dim(freqs.sig5)[2] / 2)-k] <- apply(matrix(c(freqs.sig5[,1+k], freqs.sig5[,dim(freqs.sig5)[2]-k]), ncol=2), 1, sum)
  }
  colnames(abs.freqs.all)  <- z.grid[which(z.grid > 0)] #some of the names change in an odd way.
  colnames(abs.freqs.sig1) <- z.grid[which(z.grid > 0)]
  colnames(abs.freqs.sig5) <- z.grid[which(z.grid > 0)]

  return(list(abs.freqs.all, abs.freqs.sig1, abs.freqs.sig5))
}


## This function creates the of expected significant for each meta-analysis

cf.disagg <- function(dat, z.grid) {

  # This is the loop over the meta-analyses indexed by i
  freqs <- foreach(i = 1:length(summary(dat)[,1])) %dopar% {

    # This is the loop over the observations per meta-analyses indexed by j
    freq <- matrix(ncol = (length(z.grid)-1), nrow = length(dat[[i]]$vi))
    for (j in 1:length(dat[[i]]$vi)) {

      # This is for each distribution/observation the loop over the grid indexed by a
      for (a in 1:(length(z.grid)-1)) {

        # This calculates the density between the intervals of the z.grid
        # freq is a matrix with the z.grid intervals as columns and
        # the observations/distributions per meta-analysis
        freq[j,a] <- (pnorm(z.grid[a+1], mean = dat[[i]]$GE[j] / dat[[i]]$sei[j], sd = 1) -
                        pnorm(z.grid[a], mean = dat[[i]]$GE[j] / dat[[i]]$sei[j], sd = 1))
      }
    }
    return(freq)
  }

  # freqs is a list with each element being the results for one meta-analysis. Each
  # of these lists contains the matrix "freq" (see above).
  #freqs <- do.call(rbind.data.frame, freqs)

  freqs.per.meta <- list()
  for (jjj in 1:length(summary(freqs)[,1])) {
    temp <- apply(freqs[[jjj]], 2, sum)
    temp2 <- NULL
    for (ii in 0:((length(temp) / 2) -1)) {
      temp2[(length(temp) / 2)-ii] <- sum(temp[1+ii], temp[length(temp)-ii])
    }
    freqs.per.meta[[jjj]] <- temp2
  }
  #freqs.per.meta[[ii]] <- apply(freqs[[ii]], 2, sum) # sums the columns for each interval.

  #names(freqs) <- z.grid[-1] #some of the names change in an odd way.

  # This takes the sum of the first and last column and the second and the second last
  # and so on. Therefore, it is essential that the z.grid is symmetric around zero.
  #abs.freqs <- NULL
  #for (ii in 0:((length(freqs) / 2) -1)) {

  #	abs.freqs[(length(freqs) / 2)-ii] <- sum(freqs[1+ii], freqs[length(freqs)-ii])
  #}
  #names(abs.freqs) <- z.grid[which(z.grid > 0)] #some of the names change in an odd way.

  #	[1] "0.0400000000000009" "0.0800000000000001" "0.120000000000001"  "0.16"               "0.200000000000001"
  # [6] "0.24"               "0.279999999999999"  "0.32"               "0.359999999999999"  "0.4"
  #[11] "0.44"               "0.48"               "0.52"               "0.56"               "0.6"
  #[16] "0.640000000000001"  "0.68"               "0.720000000000001"  "0.76"               "0.800000000000001"
  #[21] "0.84"               "0.880000000000001"  "0.92"               "0.960000000000001"  "1"
  #[26] "1.04"               "1.08"               "1.12"               "1.16"               "1.2"
  #[31] "1.24"               "1.28"               "1.32"               "1.36"               "1.4"
  #[36] "1.44"               "1.48"               "1.52"               "1.56"               "1.6"
  #[41] "1.64"               "1.68"               "1.72"               "1.76"               "1.8"
  #[46] "1.84"               "1.88"               "1.92"               "1.96"               "2"
  #[51] "2.04"               "2.08"               "2.12"               "2.16"               "2.2"
  #[56] "2.24"               "2.28"               "2.32"               "2.36"               "2.4"
  #[61] "2.44"               "2.48"               "2.52"               "2.56"               "2.6"

  # abs.freqs returns a vector with the cumulated probability mass for z.grid - mirrored at zero
  # The columns names are the upper limits of the intervals given by z.grid, but only for positive
  # values as the frequencies are mirrored at zero to obtain p-values for two sided tests
  # or absolute z-values for visualization.
  return(freqs.per.meta)
}


## Empirical caliper diagnostics: all observed estimates, article-cluster bootstrap.
cluster_bootstrap_calipers <- function(z, cluster, half_widths, cutoff = qnorm(.975),
                                       iterations = 1000L, seed = 20261107L) {
  if (length(iterations) != 1 || !is.finite(iterations) || iterations < 100 || iterations != trunc(iterations)) {
    stop("Use at least 100 integer bootstrap iterations.")
  }
  if (length(seed) != 1 || !is.finite(seed) || seed < 0 || seed > .Machine$integer.max || seed != trunc(seed)) stop("Invalid bootstrap seed.")
  if (!is.null(cluster) && length(cluster) != length(z)) stop("Cluster length mismatch.")
  out <- data.frame(half_width = half_widths, bootstrap_ci_lower = NA_real_,
    bootstrap_ci_upper = NA_real_, bootstrap_iterations = iterations,
    bootstrap_valid_iterations = 0L, bootstrap_clusters = NA_integer_,
    bootstrap_seed = seed, bootstrap_status = "no_cluster_column")
  if (is.null(cluster)) return(out)
  valid <- is.finite(z)
  z <- z[valid]; labels <- as.character(cluster[valid])
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) {
    out$bootstrap_status <- "missing_cluster_identifiers"; return(out)
  }
  g <- length(unique(labels)); out$bootstrap_clusters <- g
  if (g < 2L) {out$bootstrap_status <- "fewer_than_two_clusters"; return(out)}
  # Resample all articles in the retained sample, including articles with no
  # observations in a particular caliper. Recompute the ratio in every draw.
  totals <- lapply(half_widths, function(h) {
    inside <- z >= cutoff-h & z <= cutoff+h
    rowsum(cbind(n = as.numeric(inside), above = as.numeric(inside & z >= cutoff)),
      group = labels, reorder = FALSE)
  })
  n <- do.call(cbind,lapply(totals,function(x)x[,"n"]))
  above <- do.call(cbind,lapply(totals,function(x)x[,"above"]))
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else
    rm(".Random.seed",envir=.GlobalEnv))
  set.seed(seed)
  weights <- rmultinom(iterations,size=g,prob=rep(1/g,g))
  denominator <- t(weights) %*% n
  draws <- (t(weights) %*% above)/denominator
  for (j in seq_along(half_widths)) {
    finite <- is.finite(draws[,j])
    out$bootstrap_valid_iterations[j] <- sum(finite)
    if (sum(finite) < .95*iterations) {
      out$bootstrap_status[j] <- "too_many_empty_bootstrap_windows"
    } else {
      out[j,c("bootstrap_ci_lower","bootstrap_ci_upper")] <- as.list(quantile(draws[finite,j],c(.025,.975),names=FALSE))
      out$bootstrap_status[j] <- "ok_article_cluster_percentile"
    }
  }
  out
}

caliper_test <- function(z, cluster = NULL, cutoff = qnorm(.975), half_width = .1) {
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 ||
      length(half_width) != 1L || !is.finite(half_width) ||
      half_width <= 0 || half_width >= cutoff) stop("Invalid cutoff or caliper half-width.")
  if (!is.null(cluster) && length(cluster) != length(z)) stop("Cluster length mismatch.")
  inside <- is.finite(z) & z >= cutoff-half_width & z <= cutoff+half_width
  y <- z[inside] >= cutoff
  out <- data.frame(cutoff = cutoff, half_width = half_width,
    n_window = length(y), n_below = sum(!y), n_at_cutoff = sum(z[inside]==cutoff),
    n_above_including_ties = sum(y), share_above = if (length(y)) mean(y) else NA_real_,
    n_clusters = NA_integer_, cluster_status = "empty_window")
  if (!length(y)) return(out)
  if (is.null(cluster)) {out$cluster_status <- "no_cluster_column"; return(out)}
  labels <- as.character(cluster[inside])
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) {
    out$cluster_status <- "missing_cluster_identifiers"; return(out)
  }
  out$n_clusters <- length(unique(labels))
  out$cluster_status <- if (out$n_clusters < 2L) "fewer_than_two_clusters" else "ok_cluster_counts"
  out
}

prepare_caliper_data <- function(data, estimate_col = "yi", variance_col = "vi",
                                 statistic_col = NULL, df_col = NULL, cluster_col = "metaID") {
  n <- nrow(data)
  if (!n) stop("Caliper input data are empty.")
  numeric_col <- function(name) {
    if (!name %in% names(data)) stop("Missing column: ",name)
    values <- data[[name]]
    suppressWarnings(as.numeric(if (is.factor(values)) as.character(values) else values))
  }
  if (is.null(statistic_col)) {
    estimate <- numeric_col(estimate_col)
    variance <- numeric_col(variance_col)
    se <- suppressWarnings(sqrt(variance))
    statistic <- estimate/se
    valid <- is.finite(estimate) & is.finite(variance) & variance>0 & is.finite(statistic)
    reason <- ifelse(valid,"included","invalid_estimate_or_variance")
  } else {
    estimate <- variance <- se <- rep(NA_real_,n)
    statistic <- numeric_col(statistic_col)
    valid <- is.finite(statistic)
    reason <- ifelse(valid,"included","nonfinite_statistic")
  }
  z <- abs(statistic)
  df <- rep(NA_real_,n)
  if (!is.null(df_col)) {
    df <- numeric_col(df_col)
    bad <- !is.finite(df) | df<=0
    valid <- valid & !bad
    reason[bad] <- "invalid_degrees_of_freedom"
    z <- suppressWarnings(qnorm(pt(abs(statistic),df,lower.tail=FALSE,log.p=TRUE),
      lower.tail=FALSE,log.p=TRUE))
  }
  valid <- valid & is.finite(z)
  reason[reason=="included" & !valid] <- "nonfinite_equivalent_z"
  cluster <- rep(NA_character_,n)
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% names(data)) stop("Missing column: ",cluster_col)
    cluster <- as.character(data[[cluster_col]])
  }
  data.frame(row_id=seq_len(n),estimate=estimate,variance=variance,standard_error=se,
    raw_statistic=statistic,abs_z=z,df=df,cluster=cluster,valid=valid,
    exclusion_reason=reason,stringsAsFactors=FALSE)
}

compute_caliper_analysis <- function(data, half_widths = c(.05,.1,.15,.2),
                                     cutoff = qnorm(.975), bootstrap_iterations = 1000L,
                                     bootstrap_seed = 20261107L,
                                     estimate_col = "yi", variance_col = "vi",
                                     statistic_col = NULL, df_col = NULL, cluster_col = "metaID",
                                     source_data = "supplied data") {
  if (length(cutoff)!=1L || !is.finite(cutoff) || cutoff<=0 ||
      !length(half_widths) || any(!is.finite(half_widths) | half_widths<=0 |
                                 half_widths>.2 | half_widths>=cutoff)) {
    stop("Caliper half-widths must be positive, no larger than 0.2, and smaller than the cutoff.")
  }
  half_widths <- sort(unique(half_widths))
  audit <- prepare_caliper_data(data,estimate_col,variance_col,statistic_col,df_col,cluster_col)
  valid <- audit[audit$valid,,drop=FALSE]
  if (!nrow(valid)) stop("No valid observations remain for caliper analysis.")
  labels <- if (is.null(cluster_col)) NULL else valid$cluster
  bootstrap <- cluster_bootstrap_calipers(valid$abs_z,labels,half_widths,cutoff,
    bootstrap_iterations,bootstrap_seed)
  counts <- do.call(rbind,lapply(half_widths,function(h)caliper_test(valid$abs_z,labels,cutoff,h)))
  tests <- cbind(counts,bootstrap[,setdiff(names(bootstrap),"half_width"),drop=FALSE])
  tests$group <- "All"
  tests$method <- "empirical"
  coverage <- data.frame(n_input=nrow(audit),n_valid=nrow(valid),
    n_excluded_invalid=sum(!audit$valid),
    n_articles=length(unique(valid$cluster[!is.na(valid$cluster) & nzchar(trimws(valid$cluster))])),
    n_missing_cluster_ids=sum(is.na(valid$cluster) | !nzchar(trimws(valid$cluster))))
  list(tests=tests,coverage=coverage,audit=audit,
    settings=list(half_widths=half_widths,cutoff=cutoff,bootstrap_iterations=bootstrap_iterations,
      bootstrap_seed=bootstrap_seed,cluster_col=cluster_col,source_data=source_data,
      reference=if (is.null(df_col)) "normal_approximation" else "student_t_to_equivalent_z"))
}

save_caliper_analysis <- function(analysis, derived_dir) {
  dir.create(derived_dir,recursive=TRUE,showWarnings=FALSE)
  notes <- data.frame(note=c(
    paste0("Source: ",analysis$settings$source_data),
    "Pooled empirical estimates; every valid observation is retained. No reporting-precision selection is used.",
    paste0("Cutoff: ",format(analysis$settings$cutoff,digits=16),
      "; half-widths: ",paste(analysis$settings$half_widths,collapse=", ")),
    "Point estimate: observed fraction at/above cutoff among observations in the caliper, equally weighted per estimate.",
    paste0("Pointwise 95% article-cluster bootstrap percentile CIs; ",analysis$settings$bootstrap_iterations,
      " resamples; seed ",analysis$settings$bootstrap_seed,"; cluster identifier: ",
      if (is.null(analysis$settings$cluster_col)) "none" else analysis$settings$cluster_col),
    "Whole articles are resampled, including those with no observations in a given window; each draw recomputes its own ratio.",
    "No caliper p-values, analytical t intervals, binomial intervals or field-specific calipers are returned.",
    "Empty bootstrap windows are counted; intervals require at least 95% valid draws.",
    "metaID identifies meta-analysis articles; shared primary studies across articles are not captured by this clustering.",
    "Finite-window probability balance is stronger than density continuity; the diagnostics do not establish publication bias."
  ))
  openxlsx::write.xlsx(list(Results=analysis$tests,Notes=notes),
    file.path(derived_dir,"caliper_tests.xlsx"),asTable=TRUE,firstRow=TRUE,colWidths="auto",overwrite=TRUE)
  openxlsx::write.xlsx(list(Results=analysis$coverage,Notes=notes),
    file.path(derived_dir,"caliper_coverage.xlsx"),asTable=TRUE,firstRow=TRUE,colWidths="auto",overwrite=TRUE)
  saveRDS(analysis,file.path(derived_dir,"caliper_analysis.rds"))
  invisible(analysis)
}

make_caliper_plots <- function(analysis, font_size = 12) {
  settings <- analysis$settings
  valid <- analysis$audit[analysis$audit$valid,,drop=FALSE]
  limits <- settings$cutoff+c(-1,1)*max(settings$half_widths)
  z <- valid$abs_z[valid$abs_z>=limits[1] & valid$abs_z<=limits[2]]
  breaks <- seq(limits[1],limits[2],length.out=81)
  counts <- if (length(z)) hist(z,breaks=breaks,plot=FALSE)$counts else rep(0,80)
  bins <- data.frame(xmin=head(breaks,-1),xmax=tail(breaks,-1),count=counts)
  theme <- ggplot2::theme_minimal(base_size=font_size)+ggplot2::theme(
    text=ggplot2::element_text(size=font_size),
    axis.text=ggplot2::element_text(size=font_size),
    panel.grid.minor=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(face="bold"),
    plot.margin=ggplot2::margin(10,12,10,10))
  histogram <- ggplot2::ggplot(bins)+ggplot2::geom_rect(
    ggplot2::aes(xmin=xmin,xmax=xmax,ymin=0,ymax=count),fill="grey75",colour="grey20",linewidth=.25)+
    ggplot2::geom_vline(xintercept=settings$cutoff,colour="red3",linewidth=.6)+
    ggplot2::scale_x_continuous(limits=limits,expand=ggplot2::expansion(mult=0))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.05)))+
    ggplot2::labs(title="Observed |z| distribution",subtitle=paste0("n = ",format(length(z),big.mark=",",trim=TRUE)),
      x="|z|",y="Frequency")+theme
  tests <- analysis$tests
  ci <- is.finite(tests$bootstrap_ci_lower) & is.finite(tests$bootstrap_ci_upper)
  calipers <- ggplot2::ggplot(tests,ggplot2::aes(x=half_width,y=share_above))+
    ggplot2::geom_hline(yintercept=.5,colour="red3",linetype="dashed")+
    ggplot2::geom_errorbar(data=tests[ci,,drop=FALSE],
      ggplot2::aes(ymin=bootstrap_ci_lower,ymax=bootstrap_ci_upper),colour="#0072B2",width=.006)+
    ggplot2::geom_point(colour="#0072B2",size=2,na.rm=TRUE)+
    ggplot2::scale_x_continuous(breaks=tests$half_width,
      labels=paste0(sprintf("%.2f",tests$half_width),"\n(n = ",format(tests$n_window,big.mark=",",trim=TRUE),")"),
      limits=range(tests$half_width)+c(-.01,.01),expand=ggplot2::expansion(mult=0))+
    ggplot2::scale_y_continuous(breaks=seq(.45,.65,.05),labels=sprintf("%.2f",seq(.45,.65,.05)),
      expand=ggplot2::expansion(mult=0))+
    ggplot2::coord_cartesian(ylim=c(.45,.65))+
    ggplot2::labs(title="Caliper estimates",subtitle=paste0("95% cluster-bootstrap CIs; ",
      format(settings$bootstrap_iterations,big.mark=",",trim=TRUE)," resamples"),
      x="Caliper half-width (z units)",y="Share above cutoff")+theme
  if (nrow(tests)>1L) calipers <- calipers+ggplot2::geom_line(colour="#0072B2",na.rm=TRUE)
  list(histogram,calipers)
}
