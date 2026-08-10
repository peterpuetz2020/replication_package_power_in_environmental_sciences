
## ------------
## functions.R
## ------------

counterfactual_cache_version <- 3L

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
  bootstrap_results <- foreach(dd = iteration_ids) %dopar% {
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
