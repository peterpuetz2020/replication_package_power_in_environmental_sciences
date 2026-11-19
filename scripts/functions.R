
## ------------
## functions.R
## ------------

counterfactual_cache_version <- 3L

counterfactual_p_grid <- function() {
  cutoffs <- c(.001, .01, .05, .1, .2, .3, .4, .5, .6, .7, .8, .9)
  c(-Inf, qnorm(cutoffs / 2), 0, qnorm(rev(cutoffs) / 2, lower.tail = FALSE), Inf)
}

## P-value bin labels for Table 2.
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


## Use only with genuinely reported character strings.
##  Preserves zeros in "0.0400".
reported_decimal_places <- function(x) {
  if (!is.character(x)) stop("Use original reported strings to preserve trailing zeros.")
  vapply(trimws(x), function(text) {
    if (is.na(text) || !nzchar(text)) return(NA_real_)
    if (!grepl("^[+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?$", text)) {
      stop("Invalid reported numeric string: ", text)
    }
    parts <- strsplit(tolower(text), "e", fixed = TRUE)[[1]]
    pos <- regexpr(".", parts[[1]], fixed = TRUE)[[1]]
    decimals <- if (pos < 0) 0L else nchar(parts[[1]]) - pos
    decimals - if (length(parts) == 2L) as.integer(parts[[2]]) else 0L
  }, numeric(1), USE.NAMES = FALSE)
}

## Read the first sheet's reported entries before conversion to R doubles.
## General numeric cells retain their stored XML spelling; string cells retain
## trailing zeros. Simple fixed/scientific number formats retain displayed zeros.
excel_reporting_entries <- function(path, columns, n_rows) {
  if (!requireNamespace("xml2", quietly = TRUE)) stop("Install/restore xml2 first.")
  read_part <- function(member) {
    con <- unz(path, member, open = "rb")
    on.exit(close(con))
    xml2::read_xml(con)
  }
  ns <- c(s = "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
          r = "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
          p = "http://schemas.openxmlformats.org/package/2006/relationships")
  workbook <- read_part("xl/workbook.xml")
  sheet <- xml2::xml_find_first(workbook, ".//s:sheets/s:sheet[1]", ns)
  rid <- xml2::xml_attr(sheet, "r:id", ns)
  rels <- read_part("xl/_rels/workbook.xml.rels")
  target <- xml2::xml_attr(xml2::xml_find_first(rels,
    paste0("./p:Relationship[@Id='", rid, "']"), ns), "Target")
  member <- if (startsWith(target, "/")) sub("^/", "", target) else paste0("xl/", target)
  doc <- read_part(member)
  shared <- character()
  if ("xl/sharedStrings.xml" %in% unzip(path, list = TRUE)$Name) {
    shared <- xml2::xml_text(xml2::xml_find_all(read_part("xl/sharedStrings.xml"), ".//s:si", ns))
  }
  cell_text <- function(cells) {
    values <- xml2::xml_text(xml2::xml_find_first(cells, "./s:v", ns))
    type <- xml2::xml_attr(cells, "t")
    strings <- which(type == "s")
    values[strings] <- shared[as.integer(values[strings]) + 1L]
    inline <- which(type == "inlineStr")
    values[inline] <- xml2::xml_text(xml2::xml_find_first(cells[inline], "./s:is", ns))
    values
  }
  header <- xml2::xml_find_all(doc, ".//s:sheetData/s:row[@r='1']/s:c", ns)
  labels <- cell_text(header)
  letters <- sub("[0-9]+$", "", xml2::xml_attr(header, "r"))
  if (any(!columns %in% labels) || anyDuplicated(labels[labels %in% columns])) {
    stop("Excel reporting precision requires unique requested headers in worksheet row 1.")
  }
  styles <- read_part("xl/styles.xml")
  formats <- c(`0` = "General", `1` = "0", `2` = "0.00", `3` = "#,##0",
               `4` = "#,##0.00", `11` = "0.00E+00")
  custom <- xml2::xml_find_all(styles, ".//s:numFmts/s:numFmt", ns)
  formats[xml2::xml_attr(custom, "numFmtId")] <- xml2::xml_attr(custom, "formatCode")
  style_formats <- xml2::xml_attr(xml2::xml_find_all(styles, ".//s:cellXfs/s:xf", ns), "numFmtId")
  out <- lapply(columns, function(column) {
    letter <- letters[match(column, labels)]
    cells <- xml2::xml_find_all(doc, paste0(
      ".//s:sheetData/s:row/s:c[translate(@r,'0123456789','')='", letter, "']"), ns)
    row <- as.integer(sub("^[A-Z]+", "", xml2::xml_attr(cells, "r"))) - 1L
    keep <- row >= 1L & row <= n_rows
    cells <- cells[keep]; row <- row[keep]
    values <- cell_text(cells)
    type <- xml2::xml_attr(cells, "t")
    style <- suppressWarnings(as.integer(xml2::xml_attr(cells, "s")))
    style[is.na(style)] <- 0L
    codes <- unname(formats[style_formats[style + 1L]])
    formula <- !is.na(xml2::xml_text(xml2::xml_find_first(cells, "./s:f", ns)))
    decimals <- rep(NA_real_, length(values))
    source <- rep("unavailable_format_or_value", length(values))
    for (i in seq_along(values)) {
      value <- suppressWarnings(as.numeric(values[i]))
      if (!is.finite(value)) next
      string <- !is.na(type[i]) && type[i] %in% c("s", "inlineStr")
      code <- codes[i]
      if (string || (!formula[i] && !is.na(code) && code == "General")) {
        digits <- tryCatch(reported_decimal_places(values[i]), error = function(e) NA_real_)
        decimals[i] <- digits
        source[i] <- if (string) "excel_reported_string" else "excel_stored_numeric_entry"
      } else if (!is.na(code) && grepl("^[#0,]+(\\.0+)?$", code)) {
        d <- if (grepl(".", code, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", code)) else 0L
        values[i] <- formatC(value, format = "f", digits = d)
        decimals[i] <- d; source[i] <- "excel_fixed_decimal_format"
      } else if (!is.na(code) && grepl("^0\\.0+[eE][+-]0+$", code)) {
        d <- nchar(sub("^0\\.(0+)[eE].*$", "\\1", code))
        values[i] <- formatC(value, format = "e", digits = d)
        decimals[i] <- reported_decimal_places(values[i])
        source[i] <- "excel_scientific_format"
      }
    }
    result <- data.frame(entry = rep(NA_character_, n_rows), decimals = NA_real_,
                         source = "missing_cell")
    result[row, ] <- data.frame(entry = values, decimals = decimals, source = source)
    result
  })
  setNames(out, columns)
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

prepare_caliper_data <- function(data) {
  n <- nrow(data)
  if (!n) stop("Caliper input data are empty.")
  numeric_col <- function(name) {
    if (!name %in% names(data)) stop("Missing column: ",name)
    values <- data[[name]]
    suppressWarnings(as.numeric(if (is.factor(values)) as.character(values) else values))
  }
  estimate <- numeric_col("yi")
  variance <- numeric_col("vi")
  se <- suppressWarnings(sqrt(variance))
  statistic <- estimate/se
  valid <- is.finite(estimate) & is.finite(variance) & variance>0 & is.finite(statistic)
  reason <- ifelse(valid,"included","invalid_estimate_or_variance")
  if (!"metaID" %in% names(data)) stop("Missing column: metaID")
  data.frame(row_id=seq_len(n),estimate=estimate,variance=variance,standard_error=se,
    raw_statistic=statistic,abs_z=abs(statistic),cluster=as.character(data$metaID),
    valid=valid,exclusion_reason=reason,stringsAsFactors=FALSE)
}

## Keep the published Kranz-Puetz significand >= 37 rule unchanged. Its common-
## precision match is checked against sei OR vi separately. Use the matching
## decimal count for the rule. Keep observations without a match unchanged.
prepare_rounding_calipers <- function(audit, entries) {
  if (!requireNamespace("RoundingMatters", quietly=TRUE)) {
    stop("Restore the project library with renv::restore() to install RoundingMatters.")
  }
  for (column in c("yi","vi","sei")) {
    entry <- entries[[column]]
    if (is.null(entry) || nrow(entry)!=nrow(audit) ||
        !all(c("entry","decimals","source") %in% names(entry))) {
      stop("Reporting entries must contain aligned yi, vi and sei precision metadata.")
    }
    audit[[paste0("reported_",column,"_entry")]] <- entry$entry
    audit[[paste0("reported_",column,"_decimals")]] <- entry$decimals
    audit[[paste0(column,"_precision_source")]] <- entry$source
  }
  audit$kp_mu_decimals <- audit$reported_yi_decimals
  digits <- unlist(audit[c("reported_yi_decimals","reported_vi_decimals","reported_sei_decimals")],
    use.names=FALSE)
  if (any(is.finite(digits) & (digits!=trunc(digits) | abs(digits)>300))) {
    stop("Reporting decimal counts must be integers between -300 and 300.")
  }
  known <- is.finite(audit$kp_mu_decimals) &
    (is.finite(audit$reported_vi_decimals) | is.finite(audit$reported_sei_decimals))
  audit$kp_precision_known <- known
  audit$kp_matches_sei_precision <- is.finite(audit$kp_mu_decimals) &
    is.finite(audit$reported_sei_decimals) & audit$kp_mu_decimals==audit$reported_sei_decimals
  audit$kp_matches_vi_precision <- is.finite(audit$kp_mu_decimals) &
    is.finite(audit$reported_vi_decimals) & audit$kp_mu_decimals==audit$reported_vi_decimals
  audit$kp_common_precision <- audit$kp_matches_sei_precision | audit$kp_matches_vi_precision
  audit$kp_precision_match <- ifelse(!known,"unknown",
    ifelse(audit$kp_matches_sei_precision & audit$kp_matches_vi_precision,"both",
      ifelse(audit$kp_matches_sei_precision,"sei",
        ifelse(audit$kp_matches_vi_precision,"vi","none"))))
  audit$kp_rule_applied <- audit$valid & audit$kp_common_precision
  # Only applicable rows need a decimal count or significand for the KP rule.
  audit$kp_sigma_decimals <- ifelse(audit$kp_rule_applied,audit$kp_mu_decimals,NA_real_)
  audit$kp_se_significand <- NA_real_
  audit$kp_keep <- audit$valid
  rows <- which(audit$kp_rule_applied)
  if (length(rows)) {
    # The package's significand helper changes scipen. Restore the caller's
    # option so this diagnostic cannot affect later manuscript output formatting.
    old_options <- options("scipen")
    on.exit(options(old_options),add=TRUE)
    audit$kp_se_significand[rows] <- RoundingMatters::significand(
      audit$standard_error[rows],num.deci=audit$kp_sigma_decimals[rows])
    audit$kp_keep[rows] <- is.finite(audit$kp_se_significand[rows]) &
      audit$kp_se_significand[rows]>=37
  }
  audit$kp_status <- ifelse(!audit$valid,"invalid_observation",
    ifelse(!known,"retained_unknown_precision",
      ifelse(!audit$kp_common_precision,"retained_no_precision_match",
        ifelse(audit$kp_keep,"retained","omitted_by_precision_rule"))))
  audit
}

compute_caliper_analysis <- function(data, half_widths = c(.05,.1,.15,.2),
                                     cutoff = qnorm(.975), bootstrap_iterations = 1000L,
                                     bootstrap_seed = 20261107L,
                                     source_data = "supplied data", workbook_path = NULL,
                                     reporting_entries = NULL) {
  if (length(cutoff)!=1L || !is.finite(cutoff) || cutoff<=0 ||
      !length(half_widths) || any(!is.finite(half_widths) | half_widths<=0 |
                                 half_widths>.2 | half_widths>=cutoff)) {
    stop("Caliper half-widths must be positive, no larger than 0.2, and smaller than the cutoff.")
  }
  half_widths <- sort(unique(half_widths))
  audit <- prepare_caliper_data(data)
  if (!is.null(workbook_path)) {
    if (!is.null(reporting_entries)) stop("Supply either workbook_path or reporting_entries.")
    reporting_entries <- excel_reporting_entries(workbook_path,c("yi","vi","sei"),nrow(data))
  }
  use_rounding <- !is.null(reporting_entries)
  if (use_rounding) {
    audit <- prepare_rounding_calipers(audit,reporting_entries)
  }
  valid <- audit[audit$valid,,drop=FALSE]
  if (!nrow(valid)) stop("No valid observations remain for caliper analysis.")
  samples <- list(empirical=valid)
  if (use_rounding) samples$kranz_puetz <- valid[valid$kp_keep,,drop=FALSE]
  tests <- do.call(rbind,lapply(names(samples),function(method) {
    selected <- samples[[method]]
    labels <- selected$cluster
    bootstrap <- cluster_bootstrap_calipers(selected$abs_z,labels,half_widths,cutoff,
      bootstrap_iterations,bootstrap_seed)
    counts <- do.call(rbind,lapply(half_widths,function(h)caliper_test(selected$abs_z,labels,cutoff,h)))
    result <- cbind(counts,bootstrap[,setdiff(names(bootstrap),"half_width"),drop=FALSE])
    result$method <- method
    result$n_retained <- nrow(selected)
    result$n_omitted <- nrow(valid)-nrow(selected)
    result
  }))
  rownames(tests) <- NULL
  coverage <- data.frame(n_input=nrow(audit),n_valid=nrow(valid),
    n_excluded_invalid=sum(!audit$valid),
    n_articles=length(unique(valid$cluster[!is.na(valid$cluster) & nzchar(trimws(valid$cluster))])),
    n_missing_cluster_ids=sum(is.na(valid$cluster) | !nzchar(trimws(valid$cluster))))
  if (use_rounding) {
    coverage$n_known_reporting_precision <- sum(valid$kp_precision_known)
    coverage$n_unknown_reporting_precision <- sum(!valid$kp_precision_known)
    coverage$n_equal_reporting_precision <- sum(valid$kp_common_precision)
    coverage$n_unequal_reporting_precision <- sum(valid$kp_precision_known & !valid$kp_common_precision)
    coverage$n_matches_sei_only <- sum(valid$kp_precision_match=="sei")
    coverage$n_matches_vi_only <- sum(valid$kp_precision_match=="vi")
    coverage$n_matches_both <- sum(valid$kp_precision_match=="both")
    coverage$n_rule_applied <- sum(valid$kp_rule_applied)
    coverage$n_retained_after_kp_rule <- sum(valid$kp_rule_applied & valid$kp_keep)
    coverage$n_retained_without_kp_rule <- sum(!valid$kp_rule_applied & valid$kp_keep)
    coverage$n_retained_kranz_puetz <- sum(valid$kp_keep)
    coverage$n_omitted_by_precision_rule <- sum(valid$kp_rule_applied & !valid$kp_keep)
  }
  list(tests=tests,coverage=coverage,audit=audit,
    settings=list(half_widths=half_widths,cutoff=cutoff,bootstrap_iterations=bootstrap_iterations,
      bootstrap_seed=bootstrap_seed,cluster_col="metaID",source_data=source_data,
      reference="normal_approximation",
      rounding_adjustment=use_rounding,
      rounding_package=if (use_rounding) as.character(utils::packageVersion("RoundingMatters")) else NULL,
      omission_rule=if (use_rounding) "Apply SE significand >= 37 only when d_yi == d_sei OR d_yi == d_vi; retain all valid nonmatching/unknown-precision observations unchanged" else "none"))
}

save_caliper_analysis <- function(analysis, derived_dir) {
  dir.create(derived_dir,recursive=TRUE,showWarnings=FALSE)
  saveRDS(analysis,file.path(derived_dir,"caliper_analysis.rds"))
  # Retire the generated workbooks from earlier versions of this workflow.
  unlink(file.path(derived_dir,c("caliper_tests.xlsx","caliper_coverage.xlsx",
    "caliper_rounding_risks.xlsx")))
  invisible(analysis)
}

make_caliper_sample_plots <- function(analysis, font_size = 12, method = "empirical",
                                      histogram_ylim = NULL) {
  settings <- analysis$settings
  valid <- analysis$audit[analysis$audit$valid,,drop=FALSE]
  adjusted <- method=="kranz_puetz"
  if (adjusted) valid <- valid[valid$kp_keep,,drop=FALSE]
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
    ggplot2::scale_y_continuous(limits=histogram_ylim,
      expand=ggplot2::expansion(mult=if (is.null(histogram_ylim)) c(0,.05) else 0))+
    ggplot2::labs(title=if (adjusted) "Kranz-Pütz adjustment" else "Full empirical sample",
      subtitle=paste0("n = ",format(length(z),big.mark=",",trim=TRUE)),
      x="|z|",y="Frequency")+theme
  tests <- analysis$tests[analysis$tests$method==method,,drop=FALSE]
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
    ggplot2::labs(title=NULL,subtitle=NULL,
      x="Caliper half-width (z units)",y="Share above cutoff")+theme
  if (nrow(tests)>1L) calipers <- calipers+ggplot2::geom_line(colour="#0072B2",na.rm=TRUE)
  list(histogram,calipers)
}

make_caliper_plots <- function(analysis, font_size = 12) {
  full <- make_caliper_sample_plots(analysis,font_size)
  methods <- intersect(c("empirical","kranz_puetz"),analysis$tests$method)
  if (length(methods)==1L) return(full)
  histogram_ylim <- c(0,max(1,full[[1]]$data$count)*1.05)
  columns <- lapply(methods,function(method) make_caliper_sample_plots(
    analysis,font_size,method=method,histogram_ylim=histogram_ylim))
  # Histograms above calipers; each column represents one empirical sample.
  c(lapply(columns,`[[`,1L),lapply(columns,`[[`,2L))
}
