## Standalone empirical p = 0.05 diagnostics; not sourced by main.R.
## Kranz & Puetz (2022), AER 112(9):3124-3136, doi:10.1257/aer.20210121.
## Implements their primary omission rule: retain coefficient/SE reports with
## reported-SE significand >= 37. Reporting precision must not be confused with
## floating-point precision in a workbook of transformed effect sizes.
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

prepare_discontinuity_data <- function(data, estimate_col = "yi", se_col = "sei",
                                       variance_col = "vi", statistic_col = NULL,
                                       df_col = NULL, cluster_col = "metaID",
                                       subgroup_col = "subfield",
                                       estimate_digits_col = NULL, se_digits_col = NULL) {
  n <- nrow(data)
  if (!n) stop("Input data are empty.")
  pair <- is.null(statistic_col)
  estimate <- if (pair) numeric_column(data, estimate_col) else rep(NA_real_, n)
  if (pair && !is.null(se_col)) {
    se <- numeric_column(data, se_col)
    se_source <- paste0("column:", se_col)
  } else if (pair) {
    variance <- numeric_column(data, variance_col)
    se <- suppressWarnings(sqrt(variance))
    se_source <- paste0("derived_sqrt:", variance_col)
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
  sd <- if (pair) numeric_column(data, se_digits_col) else rep(NA_real_, n)
  if (any(is.finite(estimate) & is.finite(md) & !valid_precision(estimate, md)) ||
      any(is.finite(se) & is.finite(sd) & !valid_precision(se, sd))) {
    stop("Supplied reporting precision is inconsistent with the stored estimate/SE.")
  }
  known <- pair & is.finite(sd)
  inferred_sd <- inferred_md <- rep(NA_real_, n)
  if (pair) {
    inferred_sd[!known] <- inferred_decimal_places(se[!known])
    unknown_md <- !is.finite(md)
    inferred_md[unknown_md] <- inferred_decimal_places(estimate[unknown_md])
  }
  use_sd <- ifelse(known, sd, inferred_sd)
  use_md <- ifelse(is.finite(md), md, inferred_md)
  significand <- round(se / 10^(-use_sd))
  ## Bounds for |estimate / SE| under independent half-unit rounding errors.
  lower <- pmax(0, abs(estimate) - .5 * 10^(-use_md)) /
    (se + .5 * 10^(-use_sd))
  se_lower <- se - .5 * 10^(-use_sd)
  upper <- ifelse(se_lower <= 0, Inf,
                  (abs(estimate) + .5 * 10^(-use_md)) / se_lower)
  z <- equivalent_z(statistic, df)
  lower <- equivalent_z(lower, df)
  upper <- equivalent_z(upper, df)
  valid <- valid & is.finite(z)
  reason[reason == "included" & !valid] <- "nonfinite_equivalent_z"
  data.frame(
    row_id = seq_len(n), estimate = estimate, standard_error = se,
    raw_statistic = statistic, abs_z = z,
    p_two_sided = 2 * pnorm(z, lower.tail = FALSE),
    df = if (is.null(df)) rep(NA_real_, n) else df,
    reference = if (is.null(df)) "standard_normal_approximation" else "student_t_to_equivalent_z",
    cluster = cluster, subgroup = subgroup, valid = valid, exclusion_reason = reason,
    paired_input = pair, se_source = se_source,
    estimate_decimals = use_md, se_decimals = use_sd,
    common_decimal_precision = pair & use_md == use_sd,
    precision_source = ifelse(!pair, "not_applicable_direct_statistic",
                              ifelse(known, "reported_se_precision", "inferred_digits_sensitivity")),
    se_significand = significand,
    kp_keep_reported = !pair | (known & significand >= 37),
    kp_keep_sensitivity = !pair | (!is.na(significand) & significand >= 37),
    rounding_z_lower = lower, rounding_z_upper = upper,
    rounding_interval_precision = ifelse(!pair, "unavailable",
      ifelse(is.finite(md) & known, "reported", "inferred_sensitivity")),
    stringsAsFactors = FALSE
  )
}

caliper_test <- function(z, cluster = NULL, cutoff = qnorm(.975), half_width = .1) {
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 ||
      length(half_width) != 1L || !is.finite(half_width) ||
      half_width <= 0 || half_width >= cutoff) stop("Invalid cutoff or caliper half-width.")
  if (!is.null(cluster) && length(cluster) != length(z)) stop("Cluster length mismatch.")
  inside <- is.finite(z) & abs(z - cutoff) <= half_width
  y <- as.numeric(z[inside] >= cutoff)
  n <- length(y)
  n_above <- sum(y)
  out <- data.frame(
    cutoff = cutoff, p_cutoff = 2 * pnorm(cutoff, lower.tail = FALSE),
    half_width = half_width, n_window = n, n_below = n - n_above,
    n_at_cutoff = sum(z[inside] == cutoff), n_above_including_ties = n_above,
    share_above = if (n) mean(y) else NA_real_,
    binom_p_greater_iid = NA_real_, binom_p_two_sided_iid = NA_real_,
    binom_ci_lower = NA_real_, binom_ci_upper = NA_real_,
    n_clusters = NA_integer_, cluster_se = NA_real_,
    cluster_p_greater = NA_real_, cluster_p_two_sided = NA_real_,
    cluster_ci_lower = NA_real_, cluster_ci_upper = NA_real_,
    cluster_status = "empty_window"
  )
  if (!n) return(out)
  exact <- binom.test(n_above, n, p = .5)
  out$binom_p_greater_iid <- binom.test(n_above, n, p = .5, alternative = "greater")$p.value
  out$binom_p_two_sided_iid <- exact$p.value
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
  stat <- (mean(y) - .5) / se
  out$cluster_p_greater <- pt(stat, g - 1, lower.tail = FALSE)
  out$cluster_p_two_sided <- 2 * pt(abs(stat), g - 1, lower.tail = FALSE)
  margin <- qt(.975, g - 1) * se
  out$cluster_ci_lower <- max(0, mean(y) - margin)
  out$cluster_ci_upper <- min(1, mean(y) + margin)
  out$cluster_status <- "ok_cluster_t_approximation"
  out
}

run_discontinuity_tests <- function(input = "data/MasterData.xlsx",
                                     output_dir = "results/discontinuity_tests",
                                     audit_dir = "data/derived_data/discontinuity_tests",
                                     half_widths = c(.05, .075, .1, .2, .3, .4, .5),
                                     cutoff = qnorm(.975), infer_precision = TRUE,
                                     subgroup_col = "subfield",
                                     estimate_col = "yi", se_col = "sei", variance_col = "vi",
                                     statistic_col = NULL, df_col = NULL, cluster_col = "metaID",
                                     estimate_digits_col = NULL, se_digits_col = NULL,
                                     use_workbook_entries = TRUE,
                                     ...) {
  data <- if (is.data.frame(input)) input else if (grepl("\\.csv$", input, ignore.case = TRUE)) {
    read.csv(input, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
  } else {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install/restore readxl first.")
    as.data.frame(readxl::read_excel(input))
  }
  if (!length(half_widths) || any(!is.finite(half_widths) | half_widths <= 0 | half_widths >= cutoff)) {
    stop("Caliper half-widths must be positive and smaller than the cutoff.")
  }
  entries <- NULL
  if (!is.data.frame(input) && grepl("\\.xlsx$", input, ignore.case = TRUE) &&
      use_workbook_entries && is.null(statistic_col) && !is.null(se_col) &&
      is.null(estimate_digits_col) && is.null(se_digits_col)) {
    entries <- excel_reporting_entries(input, c(estimate_col, se_col), nrow(data))
    stored_estimate <- numeric_column(data, estimate_col)
    stored_se <- numeric_column(data, se_col)
    data[[estimate_col]] <- entries[[estimate_col]]$entry
    data[[se_col]] <- entries[[se_col]]$entry
    data$.estimate_decimals <- entries[[estimate_col]]$decimals
    data$.se_decimals <- entries[[se_col]]$decimals
    estimate_digits_col <- ".estimate_decimals"; se_digits_col <- ".se_decimals"
  }
  d <- prepare_discontinuity_data(data, subgroup_col = subgroup_col,
    estimate_col = estimate_col, se_col = se_col, variance_col = variance_col,
    statistic_col = statistic_col, df_col = df_col, cluster_col = cluster_col,
    estimate_digits_col = estimate_digits_col, se_digits_col = se_digits_col, ...)
  if (!is.null(entries)) {
    d$estimate_precision_source <- entries[[estimate_col]]$source
    d$se_precision_source <- entries[[se_col]]$source
    d$reported_estimate_entry <- entries[[estimate_col]]$entry
    d$reported_se_entry <- entries[[se_col]]$entry
    d$stored_workbook_estimate <- stored_estimate; d$stored_workbook_se <- stored_se
  }
  valid <- d[d$valid, , drop = FALSE]
  if (!nrow(valid)) stop("No valid observations remain; see input estimates, SEs and df.")
  levels <- sort(unique(d$subgroup[!is.na(d$subgroup) & nzchar(d$subgroup)]))
  groups <- c(list(All = rep(TRUE, nrow(d))), lapply(
    levels,
    function(x) !is.na(d$subgroup) & d$subgroup == x
  ))
  if (length(levels)) names(groups)[-1] <- paste0(subgroup_col, ": ", levels)
  modes <- c("unadjusted", "kp_reported_precision",
             if (infer_precision && any(valid$paired_input & valid$precision_source != "reported_se_precision"))
               "kp_inferred_digits_sensitivity")
  results <- list(); coverage <- list(); risks <- list()
  for (group in names(groups)) {
    rows <- groups[[group]]
    base <- d[rows & d$valid, , drop = FALSE]
    coverage[[group]] <- data.frame(
      group = group, n_input = sum(rows), n_valid = nrow(base),
      n_excluded_invalid = sum(rows & !d$valid),
      n_known_se_precision = sum(base$precision_source == "reported_se_precision"),
      n_unknown_pair_precision = sum(base$paired_input & base$precision_source != "reported_se_precision"),
      n_common_reported_decimals = sum(base$common_decimal_precision &
        base$rounding_interval_precision == "reported", na.rm = TRUE),
      n_direct_statistics = sum(!base$paired_input),
      n_retained_kp_reported = sum(base$kp_keep_reported),
      n_retained_kp_sensitivity = sum(base$kp_keep_sensitivity)
    )
    for (mode in modes) {
      keep <- switch(mode, unadjusted = rep(TRUE, nrow(base)),
                     kp_reported_precision = base$kp_keep_reported,
                     kp_inferred_digits_sensitivity = base$kp_keep_sensitivity)
      selected <- base[keep, , drop = FALSE]
      for (h in half_widths) {
        result <- caliper_test(selected$abs_z,
          if (is.null(cluster_col)) NULL else selected$cluster, cutoff, h)
        result$group <- group; result$method <- mode
        result$n_retained <- nrow(selected)
        result$n_omitted_or_unknown <- nrow(base) - nrow(selected)
        result$precision_status <- if (mode == "unadjusted") "not_required" else if (!nrow(selected)) {
          "unavailable_no_reported_precision"
        } else if (mode == "kp_inferred_digits_sensitivity") "assumed_not_verified" else "known_precision_or_direct_statistic"
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
  for (column in c("binom_p_greater_iid", "binom_p_two_sided_iid", "cluster_p_greater", "cluster_p_two_sided")) {
    available <- is.finite(results[[column]])
    adjusted <- rep(NA_real_, nrow(results))
    adjusted[available] <- p.adjust(results[[column]][available], "holm")
    results[[paste0(column, "_holm")]] <- adjusted
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(results, file.path(output_dir, "caliper_tests.csv"), row.names = FALSE)
  write.csv(do.call(rbind, coverage), file.path(output_dir, "precision_coverage.csv"), row.names = FALSE)
  write.csv(do.call(rbind, risks), file.path(output_dir, "rounding_risks.csv"), row.names = FALSE)
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(d, file.path(audit_dir, "observation_audit.rds"))
  pdf(file.path(output_dir, "empirical_discontinuity_diagnostics.pdf"), width = 10, height = 8)
  tryCatch({
    par(mfrow = c(2, 2))
    window <- valid$abs_z[abs(valid$abs_z - cutoff) <= max(half_widths)]
    if (length(window)) {
      hist(window, breaks = seq(cutoff - max(half_widths), cutoff + max(half_widths), length.out = 81),
           main = paste0("Empirical |z| near p = ", signif(2 * pnorm(cutoff, lower.tail = FALSE), 3)),
           xlab = "|z| (equivalent z for t inputs)")
      abline(v = cutoff, col = "red", lwd = 2)
    } else plot.new()
    pooled <- results[results$group == "All" & results$method == "unadjusted", ]
    plot(pooled$half_width, pooled$share_above, ylim = c(0, 1), type = "b",
         xlab = "Caliper half-width", ylab = "Share above cutoff", main = "Unadjusted; cluster 95% CI when available")
    ci <- is.finite(pooled$cluster_ci_lower) & is.finite(pooled$cluster_ci_upper)
    if (any(ci)) arrows(pooled$half_width[ci], pooled$cluster_ci_lower[ci],
                       pooled$half_width[ci], pooled$cluster_ci_upper[ci],
                       angle = 90, code = 3, length = .04)
    abline(h = .5, col = "red", lty = 2)
    pair <- valid$paired_input & abs(valid$estimate) > 0
    if (any(pair)) {
      plot(log10(valid$standard_error[pair]), log10(abs(valid$estimate[pair])), pch = 16,
           cex = .25, col = adjustcolor("black", alpha.f = .15),
           xlab = "log10(SE)", ylab = "log10(|estimate|)", main = "Joint estimate / SE distribution")
      abline(a = log10(cutoff), b = 1, col = "red", lwd = 2)
    } else plot.new()
    sensitivity <- "kp_inferred_digits_sensitivity" %in% modes
    retained <- if (sensitivity) valid$kp_keep_sensitivity else valid$kp_keep_reported
    window <- valid$abs_z[retained & abs(valid$abs_z - cutoff) <= max(half_widths)]
    if (length(window)) {
      hist(window, breaks = seq(cutoff - max(half_widths), cutoff + max(half_widths), length.out = 81),
           main = if (sensitivity) "Omission: inferred-digit SENSITIVITY" else "Omission: workbook/reported precision",
           xlab = "|z| after SE-significand >= 37 rule")
      abline(v = cutoff, col = "red", lwd = 2)
    } else plot.new()
  }, finally = dev.off())
  writeLines(c(
    "Kranz & Puetz (2022), doi:10.1257/aer.20210121; primary omission rule: SE significand >= 37.",
    "Original, empirical observations only; no counterfactual simulations or model outlier screening.",
    paste0("Cutoff = ", format(cutoff, digits = 16), "; two-sided normal p = ", 2 * pnorm(cutoff, lower.tail = FALSE), "."),
    "Window includes both endpoints; exact cutoff ties count as above, following the author replication.",
    "Nominal exact-binomial p-values assume independent observations and probability 0.5 within each window.",
    "Cluster inference is an intercept-share sandwich variance, small-sample G/(G-1), Student t(G-1).",
    "Default clusters are meta-analysis articles (metaID); unrelated articles may still share primary studies.",
    "Holm adjustment is applied across all tested rows separately for each p-value column.",
    "Finite-width symmetry is stronger than density continuity: a rejection is not proof of publication bias.",
    "For MasterData.xlsx, yi and sei are treated as reported entries, as specified by the data owner.",
    "Excel precision uses original string/numeric XML entries or simple explicit number formats; formulas with General format have unknown precision.",
    "Excel General formatting cannot recover trailing zeros already lost when a number was stored.",
    "The paper's Lemma 1 assumes common coefficient/SE decimal places; the audit flags differing precisions and computes individual rounding bounds.",
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
