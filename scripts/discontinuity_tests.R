## Standalone empirical p = 0.05 diagnostics; not sourced by main.R.
## Kranz & Puetz (2022), AER 112(9):3124-3136, doi:10.1257/aer.20210121.
## Default: use min(reported vi decimals, reported sei decimals) as the
## variance rounding resolution, propagate through sqrt(vi), then apply a
## variance-based extension of their SE-significand >= 37 omission rule.
## precision_basis = "se" retains the original rule for directly reported SEs.
## Author replication: peterpuetz2020/replication_repo_methods_matter_comment,
## commit 7220ff41dcf6b9c7648171c6d864c7871669841e, make_new_data.R / deround_bch.R.

numeric_column <- function(data, name, default = NA_real_) {
  if (is.null(name)) return(rep(default, nrow(data)))
  if (!name %in% names(data)) stop("Missing column: ", name)
  suppressWarnings(as.numeric(as.character(data[[name]])))
}

## Use only with genuinely reported character strings, not numeric-to-text
## conversions of transformed estimates. Preserves zeros in "0.0400".
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
                         source = "missing_cell", formula = FALSE)
    result[row, ] <- data.frame(entry = values, decimals = decimals, source = source, formula = formula)
    result
  })
  setNames(out, columns)
}

## Shortest numeric representations cannot recover omitted trailing zeros.
## This is used ONLY for the explicitly labelled inferred-digit sensitivity.
inferred_decimal_places <- function(x) {
  vapply(x, function(value) {
    if (!is.finite(value)) return(NA_real_)
    text <- format(value, digits = 15, scientific = TRUE, trim = TRUE)
    parts <- strsplit(tolower(text), "e", fixed = TRUE)[[1]]
    mantissa <- sub("0+$", "", parts[[1]])
    fraction <- strsplit(mantissa, ".", fixed = TRUE)[[1]]
    digits <- if (length(fraction) == 2L) nchar(fraction[[2]]) else 0L
    max(0, digits - as.integer(parts[[2]]))
  }, numeric(1))
}

valid_precision <- function(value, digits) {
  supplied <- is.finite(digits)
  if (any(supplied & (digits != trunc(digits) | abs(digits) > 300))) {
    stop("Reporting decimal places must be integers between -300 and 300.")
  }
  rounded <- round(value, digits)
  supplied & is.finite(value) &
    abs(value - rounded) <= 32 * .Machine$double.eps *
      pmax(abs(value), abs(rounded), .Machine$double.xmin)
}

## Equivalent |z| makes a common cutoff correspond to two-sided p = 0.05
## even for reported t statistics with observation-specific degrees of freedom.
equivalent_z <- function(statistic, df = NULL) {
  if (is.null(df)) return(abs(statistic))
  qnorm(pt(abs(statistic), df, lower.tail = FALSE, log.p = TRUE),
        lower.tail = FALSE, log.p = TRUE)
}

## A variance reported to d decimals represents v +/- 0.5 * 10^(-d).
## Propagate that interval exactly, rather than treating the digits of sqrt(v)
## as reported SE precision. Rationalized differences avoid cancellation when
## the variance is very precisely reported. The equivalent symmetric SE step
## is twice the larger one-sided SE uncertainty; its precision score generalizes
## SE / reported_SE_step but is NOT a reported integer SE significand.
variance_precision_bounds <- function(variance, decimals) {
  step <- 10^(-decimals)
  se <- suppressWarnings(sqrt(variance))
  variance_lower <- pmax(0, variance - .5 * step)
  variance_upper <- variance + .5 * step
  se_lower <- suppressWarnings(sqrt(variance_lower))
  se_upper <- suppressWarnings(sqrt(variance_upper))
  lower_error <- ifelse(variance_lower > 0, .5 * step / (se + se_lower), se)
  upper_error <- .5 * step / (se_upper + se)
  equivalent_step <- 2 * pmax(lower_error, upper_error)
  data.frame(variance_lower = variance_lower, variance_upper = variance_upper,
    se_lower = se_lower, se_upper = se_upper,
    equivalent_se_step = equivalent_step,
    equivalent_se_precision_score = se / equivalent_step)
}

prepare_discontinuity_data <- function(data, estimate_col = "yi", se_col = "sei",
                                       variance_col = "vi", statistic_col = NULL,
                                       df_col = NULL, cluster_col = "metaID",
                                       subgroup_col = "subfield",
                                       estimate_digits_col = NULL, se_digits_col = NULL,
                                       variance_digits_col = NULL,
                                       precision_basis = c("min_variance_se", "variance", "se")) {
  precision_basis <- match.arg(precision_basis)
  n <- nrow(data)
  if (!n) stop("Input data are empty.")
  pair <- is.null(statistic_col)
  variance_based <- pair && precision_basis != "se"
  minimum_based <- pair && precision_basis == "min_variance_se"
  estimate <- if (pair) numeric_column(data, estimate_col) else rep(NA_real_, n)
  variance <- rep(NA_real_, n)
  if (variance_based || (pair && is.null(se_col))) {
    variance <- numeric_column(data, variance_col)
    se <- suppressWarnings(sqrt(variance))
    se_source <- paste0("derived_sqrt:", variance_col)
  } else if (pair) {
    se <- numeric_column(data, se_col)
    se_source <- paste0("column:", se_col)
  } else {
    se <- rep(NA_real_, n)
    se_source <- "not_applicable_direct_statistic"
  }
  statistic <- if (pair) estimate / se else numeric_column(data, statistic_col)
  df <- if (is.null(df_col)) NULL else numeric_column(data, df_col)
  valid <- is.finite(statistic)
  reason <- ifelse(valid, "included", "nonfinite_estimate_or_statistic")
  if (pair) {
    bad <- !is.finite(se) | se <= 0
    reason[bad] <- "nonpositive_or_nonfinite_se"
    valid <- valid & !bad & is.finite(estimate)
    if (variance_based) {
      bad_variance <- !is.finite(variance) | variance <= 0
      reason[bad_variance] <- "nonpositive_or_nonfinite_variance"
    }
  }
  if (!is.null(df)) {
    bad <- !is.finite(df) | df <= 0
    reason[bad] <- "invalid_degrees_of_freedom"
    valid <- valid & !bad
  }
  get_text <- function(name, fallback) {
    if (is.null(name)) return(rep(fallback, n))
    if (!name %in% names(data)) stop("Missing column: ", name)
    as.character(data[[name]])
  }
  cluster <- get_text(cluster_col, NA_character_)
  subgroup <- get_text(subgroup_col, NA_character_)
  md <- if (pair) numeric_column(data, estimate_digits_col) else rep(NA_real_, n)
  sd <- if (pair && (!variance_based || minimum_based)) numeric_column(data, se_digits_col) else rep(NA_real_, n)
  vd <- if (variance_based) numeric_column(data, variance_digits_col) else rep(NA_real_, n)
  reported_se <- if (minimum_based) numeric_column(data, se_col) else se
  if (any(is.finite(estimate) & is.finite(md) & !valid_precision(estimate, md)) ||
      any(is.finite(reported_se) & is.finite(sd) & !valid_precision(reported_se, sd)) ||
      any(is.finite(variance) & is.finite(vd) & !valid_precision(variance, vd))) {
    stop("Supplied reporting precision is inconsistent with the stored estimate/SE/variance.")
  }
  known <- pair & is.finite(if (variance_based) vd else sd)
  if (minimum_based) known <- known & is.finite(sd)
  inferred_sd <- inferred_md <- inferred_vd <- rep(NA_real_, n)
  if (pair) {
    if (variance_based) {
      unknown_vd <- !is.finite(vd)
      inferred_vd[unknown_vd] <- inferred_decimal_places(variance[unknown_vd])
      if (minimum_based) {
        unknown_sd <- !is.finite(sd)
        inferred_sd[unknown_sd] <- inferred_decimal_places(reported_se[unknown_sd])
      }
    } else inferred_sd[!known] <- inferred_decimal_places(se[!known])
    unknown_md <- !is.finite(md)
    inferred_md[unknown_md] <- inferred_decimal_places(estimate[unknown_md])
  }
  use_sd <- ifelse(is.finite(sd), sd, inferred_sd)
  original_vd <- ifelse(is.finite(vd), vd, inferred_vd)
  use_vd <- if (minimum_based) pmin(original_vd, use_sd) else original_vd
  use_md <- ifelse(is.finite(md), md, inferred_md)
  significand <- if (variance_based) rep(NA_real_, n) else round(se / 10^(-use_sd))
  variance_significand <- round(variance / 10^(-original_vd))
  variance_bounds <- variance_precision_bounds(variance, use_vd)
  se_lower <- if (variance_based) variance_bounds$se_lower else se - .5 * 10^(-use_sd)
  se_upper <- if (variance_based) variance_bounds$se_upper else se + .5 * 10^(-use_sd)
  omission_score <- if (variance_based) variance_bounds$equivalent_se_precision_score else significand
  ## Bounds for |estimate / SE| using the coefficient's own reported resolution.
  lower <- pmax(0, abs(estimate) - .5 * 10^(-use_md)) /
    se_upper
  upper <- ifelse(se_lower <= 0, Inf,
                  (abs(estimate) + .5 * 10^(-use_md)) / se_lower)
  z <- equivalent_z(statistic, df)
  lower <- equivalent_z(lower, df)
  upper <- equivalent_z(upper, df)
  valid <- valid & is.finite(z)
  reason[reason == "included" & !valid] <- "nonfinite_equivalent_z"
  data.frame(
    row_id = seq_len(n), estimate = estimate, standard_error = se,
    variance = variance, variance_decimals = use_vd,
    reported_variance_decimals = vd, reported_se_decimals = sd,
    variance_resolution_ratio = variance / 10^(-use_vd),
    raw_statistic = statistic, abs_z = z,
    p_two_sided = 2 * pnorm(z, lower.tail = FALSE),
    df = if (is.null(df)) rep(NA_real_, n) else df,
    reference = if (is.null(df)) "standard_normal_approximation" else "student_t_to_equivalent_z",
    cluster = cluster, subgroup = subgroup, valid = valid, exclusion_reason = reason,
    paired_input = pair, se_source = se_source,
    precision_basis = if (pair) precision_basis else "direct_statistic",
    estimate_decimals = use_md, se_decimals = if (variance_based) NA_real_ else use_sd,
    common_decimal_precision = if (variance_based) rep(NA, n) else pair & use_md == use_sd,
    precision_source = ifelse(!pair, "not_applicable_direct_statistic",
      ifelse(known, if (minimum_based) "reported_min_variance_se_precision" else if (variance_based) "reported_variance_precision" else "reported_se_precision",
             "inferred_digits_sensitivity")),
    se_significand = significand,
    variance_significand = variance_significand,
    equivalent_se_precision_score = variance_bounds$equivalent_se_precision_score,
    omission_precision_score = omission_score,
    rounding_se_lower = se_lower, rounding_se_upper = se_upper,
    rounding_variance_lower = variance_bounds$variance_lower,
    rounding_variance_upper = variance_bounds$variance_upper,
    kp_keep_reported = !pair | (known & omission_score >= 37),
    kp_keep_sensitivity = !pair | (!is.na(omission_score) & omission_score >= 37),
    rounding_z_lower = lower, rounding_z_upper = upper,
    rounding_interval_precision = ifelse(!pair, "unavailable",
      ifelse(is.finite(md) & known, "reported", "inferred_sensitivity")),
    stringsAsFactors = FALSE
  )
}

cluster_bootstrap_calipers <- function(z, cluster, half_widths, cutoff = qnorm(.975),
                                       iterations = 9999L, seed = 20261107L) {
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
  inside <- is.finite(z) & z >= cutoff - half_width & z <= cutoff + half_width
  y <- as.numeric(z[inside] >= cutoff)
  n <- length(y)
  n_above <- sum(y)
  out <- data.frame(
    cutoff = cutoff,
    half_width = half_width, n_window = n, n_below = n - n_above,
    n_at_cutoff = sum(z[inside] == cutoff), n_above_including_ties = n_above,
    share_above = if (n) mean(y) else NA_real_,
    binom_ci_lower = NA_real_, binom_ci_upper = NA_real_,
    n_clusters = NA_integer_, cluster_se = NA_real_,
    cluster_ci_lower = NA_real_, cluster_ci_upper = NA_real_,
    cluster_status = "empty_window"
  )
  if (!n) return(out)
  exact <- binom.test(n_above, n, p = .5)
  out$binom_ci_lower <- exact$conf.int[1]
  out$binom_ci_upper <- exact$conf.int[2]
  if (is.null(cluster)) {
    out$cluster_status <- "no_cluster_column"
    return(out)
  }
  labels <- as.character(cluster[inside])
  if (anyNA(labels) || any(!nzchar(trimws(labels)))) {
    out$cluster_status <- "missing_cluster_identifiers"
    return(out)
  }
  scores <- tapply(y - mean(y), labels, sum)
  g <- length(scores)
  out$n_clusters <- g
  if (g < 2L) {
    out$cluster_status <- "fewer_than_two_clusters"
    return(out)
  }
  se <- sqrt(g / (g - 1) * sum(scores^2) / n^2)
  out$cluster_se <- se
  if (!is.finite(se) || se <= sqrt(.Machine$double.eps)) {
    out$cluster_status <- "degenerate_cluster_variance"
    return(out)
  }
  margin <- qt(.975, g - 1) * se
  out$cluster_ci_lower <- max(0, mean(y) - margin)
  out$cluster_ci_upper <- min(1, mean(y) + margin)
  out$cluster_status <- "ok_cluster_t_approximation"
  out
}

run_discontinuity_tests <- function(input = "data/MasterData.xlsx",
                                     output_dir = "results/discontinuity_tests",
                                     audit_dir = "data/derived_data/discontinuity_tests",
                                     half_widths = c(.05, .1, .15, .2),
                                     cutoff = qnorm(.975), infer_precision = TRUE,
                                     subgroup_col = NULL,
                                     bootstrap_iterations = 9999L, bootstrap_seed = 20261107L,
                                     estimate_col = "yi", se_col = "sei", variance_col = "vi",
                                     statistic_col = NULL, df_col = NULL, cluster_col = "metaID",
                                     estimate_digits_col = NULL, se_digits_col = NULL,
                                     variance_digits_col = NULL,
                                     precision_basis = c("min_variance_se", "variance", "se"),
                                     use_workbook_entries = TRUE,
                                     ...) {
  precision_basis <- match.arg(precision_basis)
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Install/restore openxlsx first.")
  variance_based <- is.null(statistic_col) && precision_basis != "se"
  minimum_based <- is.null(statistic_col) && precision_basis == "min_variance_se"
  precision_col <- if (variance_based) variance_col else se_col
  data <- if (is.data.frame(input)) input else if (grepl("\\.csv$", input, ignore.case = TRUE)) {
    read.csv(input, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
  } else {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install/restore readxl first.")
    as.data.frame(readxl::read_excel(input))
  }
  if (!length(half_widths) || any(!is.finite(half_widths) | half_widths <= 0 |
                                half_widths > .2 | half_widths >= cutoff)) {
    stop("Caliper half-widths must be positive, no larger than 0.2, and smaller than the cutoff.")
  }
  entries <- NULL
  if (!is.data.frame(input) && grepl("\\.xlsx$", input, ignore.case = TRUE) &&
      use_workbook_entries && is.null(statistic_col) && !is.null(precision_col) &&
      is.null(estimate_digits_col) &&
      is.null(if (variance_based) variance_digits_col else se_digits_col) &&
      (!minimum_based || is.null(se_digits_col))) {
    entries <- excel_reporting_entries(input, unique(c(estimate_col, precision_col,
      if (minimum_based) se_col)), nrow(data))
    stored_estimate <- numeric_column(data, estimate_col)
    stored_precision_value <- numeric_column(data, precision_col)
    stored_se <- if (!is.null(se_col) && se_col %in% names(data)) numeric_column(data, se_col) else NULL
    data[[estimate_col]] <- entries[[estimate_col]]$entry
    data[[precision_col]] <- entries[[precision_col]]$entry
    data$.estimate_decimals <- entries[[estimate_col]]$decimals
    data$.precision_decimals <- entries[[precision_col]]$decimals
    estimate_digits_col <- ".estimate_decimals"
    if (variance_based) variance_digits_col <- ".precision_decimals" else se_digits_col <- ".precision_decimals"
    if (minimum_based) {
      data[[se_col]] <- entries[[se_col]]$entry
      data$.se_decimals <- entries[[se_col]]$decimals
      se_digits_col <- ".se_decimals"
    }
  }
  d <- prepare_discontinuity_data(data, subgroup_col = subgroup_col,
    estimate_col = estimate_col, se_col = se_col, variance_col = variance_col,
    statistic_col = statistic_col, df_col = df_col, cluster_col = cluster_col,
    estimate_digits_col = estimate_digits_col, se_digits_col = se_digits_col,
    variance_digits_col = variance_digits_col, precision_basis = precision_basis, ...)
  if (!is.null(entries)) {
    d$estimate_precision_source <- entries[[estimate_col]]$source
    d$reported_estimate_entry <- entries[[estimate_col]]$entry
    if (variance_based) {
      d$variance_precision_source <- entries[[variance_col]]$source
      d$reported_variance_entry <- entries[[variance_col]]$entry
      d$stored_workbook_variance <- stored_precision_value
      if (minimum_based) {
        d$se_precision_source <- entries[[se_col]]$source
        d$reported_se_entry <- entries[[se_col]]$entry
      }
    } else {
      d$se_precision_source <- entries[[se_col]]$source
      d$reported_se_entry <- entries[[se_col]]$entry
    }
    d$stored_workbook_estimate <- stored_estimate
    if (!is.null(stored_se)) d$stored_workbook_se <- stored_se
  }
  valid <- d[d$valid, , drop = FALSE]
  if (!nrow(valid)) stop("No valid observations remain; see input estimates, SEs and df.")
  # Caliper diagnostics are pooled only, irrespective of field metadata.
  groups <- list(All = rep(TRUE, nrow(d)))
  reported_source <- if (minimum_based) "reported_min_variance_se_precision" else if (variance_based) "reported_variance_precision" else "reported_se_precision"
  reported_mode <- if (minimum_based) "kp_reported_min_variance_se_precision" else if (variance_based) "kp_reported_variance_precision" else "kp_reported_precision"
  sensitivity_mode <- if (variance_based) "kp_inferred_variance_digits_sensitivity" else "kp_inferred_digits_sensitivity"
  omission_rule <- if (minimum_based) "variance_propagated_SE_precision_score>=37;d=min(d_vi,d_sei)" else if (variance_based) "variance_propagated_SE_precision_score>=37" else "reported_SE_significand>=37"
  modes <- c("unadjusted", reported_mode,
             if (infer_precision && any(valid$paired_input & valid$precision_source != reported_source)) sensitivity_mode)
  results <- list(); coverage <- list(); risks <- list()
  for (group in names(groups)) {
    rows <- groups[[group]]
    base <- d[rows & d$valid, , drop = FALSE]
    coverage[[group]] <- data.frame(
      group = group, n_input = sum(rows), n_valid = nrow(base),
      n_excluded_invalid = sum(rows & !d$valid),
      n_known_se_precision = sum(is.finite(base$reported_se_decimals)),
      n_known_variance_precision = sum(is.finite(base$reported_variance_decimals)),
      n_known_min_variance_se_precision = sum(base$precision_source == "reported_min_variance_se_precision"),
      n_unknown_pair_precision = sum(base$paired_input & base$precision_source != reported_source),
      n_common_reported_decimals = if (variance_based) NA_integer_ else sum(base$common_decimal_precision &
        base$rounding_interval_precision == "reported", na.rm = TRUE),
      n_direct_statistics = sum(!base$paired_input),
      n_retained_kp_reported = sum(base$kp_keep_reported),
      n_retained_kp_sensitivity = sum(base$kp_keep_sensitivity)
    )
    for (mode in modes) {
      keep <- if (mode == "unadjusted") rep(TRUE, nrow(base)) else if (mode == reported_mode) {
        base$kp_keep_reported
      } else base$kp_keep_sensitivity
      selected <- base[keep, , drop = FALSE]
      bootstrap <- cluster_bootstrap_calipers(selected$abs_z,
        if (is.null(cluster_col)) NULL else selected$cluster,
        half_widths, cutoff, bootstrap_iterations, bootstrap_seed)
      for (h in half_widths) {
        result <- caliper_test(selected$abs_z,
          if (is.null(cluster_col)) NULL else selected$cluster, cutoff, h)
        result <- cbind(result, bootstrap[match(h,bootstrap$half_width),
          setdiff(names(bootstrap),"half_width"),drop=FALSE])
        result$group <- group; result$method <- mode
        result$precision_basis <- if (is.null(statistic_col)) precision_basis else "direct_statistic"
        result$omission_rule <- if (mode == "unadjusted") "none" else omission_rule
        result$n_retained <- nrow(selected)
        result$n_omitted_or_unknown <- nrow(base) - nrow(selected)
        result$precision_status <- if (mode == "unadjusted") "not_required" else if (!nrow(selected)) {
          if (any(base$precision_source == reported_source)) "all_observations_omitted_by_precision_rule" else {
            "unavailable_no_reported_precision"
          }
        } else if (mode == sensitivity_mode) "assumed_not_verified" else if (variance_based) {
          "reported_variance_precision_extension"
        } else "known_precision_or_direct_statistic"
        results[[length(results) + 1L]] <- result
        known <- is.finite(selected$rounding_z_lower) & !is.na(selected$rounding_z_upper)
        lo <- selected$rounding_z_lower[known]; hi <- selected$rounding_z_upper[known]
        definitely_in <- lo >= cutoff - h & hi <= cutoff + h
        definitely_out <- hi < cutoff - h | lo > cutoff + h
        risks[[length(risks) + 1L]] <- data.frame(
          group = group, method = mode, half_width = h,
          n_with_rounding_bounds = sum(known),
          n_bounds_inferred = sum(known & selected$rounding_interval_precision != "reported"),
          n_potential_significance_switch = sum(lo < cutoff & hi >= cutoff),
          n_uncertain_window_membership = sum(!definitely_in & !definitely_out)
        )
      }
    }
  }
  results <- do.call(rbind, results)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  workbook_notes <- c(
    paste0("Input data: ", if (is.data.frame(input)) "supplied data frame" else input),
    paste0("Caliper half-widths in equivalent z units: ", paste(half_widths, collapse = ", ")),
    paste0("Precision basis: ", if (is.null(statistic_col)) precision_basis else "direct_statistic"),
    if (minimum_based) "Variance resolution uses d = min(reported vi decimals, reported sei decimals) per observation; point estimates are not rounded to the selected resolution.",
    paste0("Omission criterion: ", omission_rule,
           if (variance_based) " (variance-based extension of Kranz-Puetz, not the literal published SE rule)." else "."),
    paste0("Cluster identifier: ", if (is.null(cluster_col)) "none" else cluster_col),
    "Plotted CIs compare pointwise 95% article-clustered t(G-1) intervals (G/(G-1) sandwich correction) with cluster-bootstrap percentile intervals. Not simultaneous intervals.",
    "Binomial CIs in the test workbook are exact 95% intervals assuming independent observations; they are not the plotted intervals.",
    paste0("Pooled sample only; no caliper p-values. Cluster-bootstrap sensitivity: ", bootstrap_iterations,
      " whole-article resamples, seed ",bootstrap_seed,"; 2.5%/97.5% percentile intervals, with window counts recomputed in each draw."),
    "Blank cells indicate unavailable/not-applicable results; consult the status columns."
  )
  write_workbook <- function(table, stem, description) {
    openxlsx::write.xlsx(list(Results = table, Notes = data.frame(note = c(description, workbook_notes))),
      file = file.path(output_dir, paste0(stem, ".xlsx")), asTable = TRUE,
      firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  }
  write_workbook(results, "caliper_tests", "Pooled sample, one row per method and caliper width; counts, shares and confidence intervals; no p-values.")
  write_workbook(do.call(rbind, coverage), "precision_coverage", "One row per group; valid input and precision-filter coverage before caliper restriction.")
  write_workbook(do.call(rbind, risks), "rounding_risks", "Possible rounding risks across the entire retained group, not just observed window members; these are not confirmed errors.")
  unlink(file.path(output_dir, c("caliper_tests.csv", "precision_coverage.csv", "rounding_risks.csv")))
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(d, file.path(audit_dir, "observation_audit.rds"))
  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Install/restore ggplot2 and gridExtra first.")
  }
  sensitivity <- sensitivity_mode %in% modes
  adjusted_mode <- if (sensitivity) sensitivity_mode else reported_mode
  retained <- if (sensitivity) valid$kp_keep_sensitivity else valid$kp_keep_reported
  # The histogram and largest caliper use exactly the same window.
  histogram_limits <- cutoff + c(-1, 1) * max(half_widths)
  histogram_breaks <- seq(histogram_limits[1], histogram_limits[2], length.out = 81)
  in_display <- valid$abs_z >= histogram_limits[1] & valid$abs_z <= histogram_limits[2]
  hist_counts <- function(z) {
    counts <- if (length(z)) hist(z, breaks = histogram_breaks, plot = FALSE)$counts else rep(0, 80)
    data.frame(xmin = head(histogram_breaks, -1), xmax = tail(histogram_breaks, -1), count = counts)
  }
  full_hist <- hist_counts(valid$abs_z[in_display])
  adjusted_hist <- hist_counts(valid$abs_z[in_display & retained])
  histogram_ylim <- c(0, max(1, full_hist$count) * 1.05)
  plot_theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9),
      plot.margin = ggplot2::margin(10, 12, 10, 10))
  histogram_plot <- function(counts, title, n) {
    ggplot2::ggplot(counts) +
      ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count),
        fill = "grey75", colour = "grey20", linewidth = .25) +
      ggplot2::geom_vline(xintercept = cutoff, colour = "red3", linewidth = .6) +
      ggplot2::scale_x_continuous(limits = histogram_limits,
        breaks = pretty(histogram_limits, n = 4), expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_y_continuous(limits = histogram_ylim, expand = ggplot2::expansion(mult = 0)) +
      ggplot2::labs(title = title, subtitle = paste0("n = ", format(n, big.mark = ",", trim = TRUE)),
        x = "|z|", y = "Frequency") + plot_theme
  }
  caliper_plot <- function(method, title) {
    pooled <- results[results$group == "All" & results$method == method, ]
    pooled <- pooled[order(pooled$half_width), ]
    ci <- is.finite(pooled$cluster_ci_lower) & is.finite(pooled$cluster_ci_upper)
    boot_ci <- is.finite(pooled$bootstrap_ci_lower) & is.finite(pooled$bootstrap_ci_upper)
    intervals <- rbind(
      data.frame(x_plot=pooled$half_width[ci]-.003,
        lower=pooled$cluster_ci_lower[ci],upper=pooled$cluster_ci_upper[ci],ci_method=rep("Clustered t",sum(ci))),
      data.frame(x_plot=pooled$half_width[boot_ci]+.003,
        lower=pooled$bootstrap_ci_lower[boot_ci],upper=pooled$bootstrap_ci_upper[boot_ci],ci_method=rep("Cluster bootstrap",sum(boot_ci))))
    ci_label <- if (identical(cluster_col,"metaID")) {
      "Pointwise 95% CIs; meta-analysis article clusters"
    } else "Pointwise 95% confidence intervals"
    g <- ggplot2::ggplot(pooled, ggplot2::aes(x = half_width, y = share_above)) +
      ggplot2::geom_hline(yintercept = .5, colour = "red3", linetype = "dashed") +
      ggplot2::geom_errorbar(data = intervals,
        ggplot2::aes(x=x_plot,ymin=lower,ymax=upper,colour=ci_method),
        inherit.aes=FALSE,width=.004) +
      ggplot2::scale_colour_manual(values=c("Clustered t"="black","Cluster bootstrap"="#0072B2"),
        limits=c("Clustered t","Cluster bootstrap"),drop=FALSE,name=NULL) +
      ggplot2::geom_point(size = 2, na.rm = TRUE) +
      ggplot2::scale_x_continuous(breaks = pooled$half_width,
        labels = paste0(sprintf("%.2f", pooled$half_width), "\n(n = ",
                        format(pooled$n_window, big.mark = ",", trim = TRUE), ")"),
        limits = range(pooled$half_width) + c(-.01, .01), expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_y_continuous(breaks = seq(.35, .7, .05),
        labels = sprintf("%.2f", seq(.35, .7, .05)), expand = ggplot2::expansion(mult = 0)) +
      ggplot2::coord_cartesian(ylim = c(.35, .7)) +
      ggplot2::labs(title = title, subtitle = ci_label,
        x = "Caliper half-width (z units)", y = "Share above cutoff") + plot_theme +
      ggplot2::theme(legend.position="bottom",legend.text=ggplot2::element_text(size=9))
    if (nrow(pooled) > 1L) g <- g + ggplot2::geom_line(na.rm = TRUE)
    g
  }
  adjusted_title <- if (sensitivity) "Kranz-Puetz: precision sensitivity" else "Kranz-Puetz omission sample"
  plots <- list(
    histogram_plot(full_hist, "Full empirical sample", sum(in_display)),
    caliper_plot("unadjusted", "Calipers: full empirical sample"),
    histogram_plot(adjusted_hist, adjusted_title, sum(in_display & retained)),
    caliper_plot(adjusted_mode, paste0("Calipers: ", adjusted_title))
  )
  diagnostic <- local({
    # Grid needs a device for text metrics; avoid an implicit Rplots.pdf.
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    gridExtra::arrangeGrob(grobs = plots, ncol = 2)
  })
  ggplot2::ggsave(file.path(output_dir, "empirical_discontinuity_diagnostics.svg"),
    plot = diagnostic, device = grDevices::svg, width = 10, height = 8)
  ggplot2::ggsave(file.path(output_dir, "empirical_discontinuity_diagnostics.pdf"),
    plot = diagnostic, device = grDevices::pdf, width = 10, height = 8)
  writeLines(c(
    "Kranz & Puetz (2022), doi:10.1257/aer.20210121; published omission rule: reported SE significand >= 37.",
    if (variance_based) c(
      "THIS RUN EXTENDS THE RULE TO REPORTED VARIANCES; it does not apply the published SE rule literally.",
      paste0("Precision comes from ", variance_col, "; SE = sqrt(", variance_col, "); digits of a computed sqrt are not treated as reported precision."),
      if (minimum_based) paste0("Use d = min(reported decimals of ", variance_col, ", reported decimals of ", se_col, ") for each observation; this is a conservative resolution assumption around the stored variance, without rounding its value."),
      "For a variance v reported to d decimals, Delta_v = 10^(-d), v_L = max(0, v-Delta_v/2), v_U = v+Delta_v/2.",
      "SE bounds are sqrt(v_L) and sqrt(v_U). Equivalent SE step = 2*max(sqrt(v)-sqrt(v_L), sqrt(v_U)-sqrt(v)).",
      "The variance-based omission retains sqrt(v)/equivalent_SE_step >= 37; this score is not a reported integer SE significand."
    ) else "THIS RUN USES REPORTED SE PRECISION; no variance-to-SE extension is applied.",
    "Original, empirical observations only; no counterfactual simulations or model outlier screening.",
    paste0("Cutoff = ", format(cutoff, digits = 16), "; two-sided normal p = ", 2 * pnorm(cutoff, lower.tail = FALSE), "."),
    "Window includes both endpoints; exact cutoff ties count as above, following the author replication.",
    "Exact-binomial CIs assume independent observations. Caliper p-values are not computed or exported.",
    "Cluster inference is an intercept-share sandwich variance, small-sample G/(G-1), Student t(G-1).",
    "Default clusters are meta-analysis articles (metaID); unrelated articles may still share primary studies.",
    "Field-specific calipers are omitted. All caliper output concerns the pooled empirical sample.",
    paste0("Bootstrap sensitivity: ",bootstrap_iterations," whole-article resamples; seed ",bootstrap_seed,
      "; percentile 95% CIs. Empty windows are reported; bounds require >=95% valid draws."),
    paste0("This run tests half-widths ", paste(half_widths, collapse = ", "), " only; maximum permitted half-width is 0.2."),
    "Diagnostic layout: full empirical histogram/caliper panel above; Kranz-Puetz omission-sample histogram/caliper panel below.",
    "Figure uses ggplot2 and is saved as SVG and PDF. Histograms filter both data and x range to cutoff +/- max(half_widths), approximately 1.76 to 2.16 by default; the exact normal cutoff qnorm(0.975) matches the main workflow.",
    "Caliper x-axis labels give each window's own sample size as (n = ...) on the second line.",
    "Both caliper panels have y limits 0.35 to 0.70 and compare pointwise 95% article-clustered t(G-1) and cluster-bootstrap percentile CIs. These are not simultaneous intervals.",
    "Finite-width symmetry is stronger than density continuity: a rejection is not proof of publication bias.",
    paste0("Reporting precision basis: ", if (is.null(statistic_col)) precision_basis else "direct_statistic", ". Coefficient resolution uses its own reported entries."),
    "Excel precision uses original string/numeric XML entries or simple explicit number formats; formulas with General format have unknown precision.",
    "Excel General formatting cannot recover trailing zeros already lost when a number was stored.",
    "The paper's Lemma 1 assumes common coefficient/SE decimal places; the variance extension is not that setting. Individual coefficient and SE bounds are audited.",
    "Interpreting the omission adjustment requires that reporting coarseness is not systematically related to the unrounded z distribution.",
    "Inferred digits are ONLY a sensitivity assumption; trailing zeros and original rounding cannot be recovered.",
    "The omission rule reduces rounding problems but does not guarantee every retained observation is correctly classified.",
    "Rounding risk bounds are diagnostics conditional on the supplied or assumed decimal resolutions.",
    "Student t inputs require df_col; without df the estimate/SE ratio uses a normal approximation.",
    "The same calipers cannot diagnose separate marginal density jumps in estimates or SEs: the boundary is |estimate| = cutoff * SE."
  ), file.path(output_dir, "METHODS.txt"))
  message("Wrote ", nrow(results), " caliper results to ", output_dir)
  invisible(list(tests = results, coverage = do.call(rbind, coverage), audit = d))
}

if (sys.nframe() == 0L) run_discontinuity_tests()
