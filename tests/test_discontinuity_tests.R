source("scripts/discontinuity_tests.R")
stopifnot(identical(reported_decimal_places(c("0.0400", "1.20e-3", "100", "1e3", NA)),
                    c(4, 5, 0, -3, NA_real_)))

c <- qnorm(.975)
z <- c + ifelse(c(0, 1, 1, 0, 0, 1, 1, 1), .05, -.05)
clusters <- c("a", "a", "a", "b", "b", "c", "c", "d")
test <- caliper_test(z, clusters, c, .1)
expected_se <- sqrt(4 / 3 * (.125^2 + (-1.25)^2 + .75^2 + .375^2) / 8^2)
stopifnot(test$n_window == 8, test$n_above_including_ties == 5,
          abs(test$binom_p_greater_iid - binom.test(5, 8, .5, alternative = "greater")$p.value) < 1e-14,
          abs(test$cluster_se - expected_se) < 1e-14,
          abs(test$cluster_p_greater - pt(.125 / expected_se, 3, lower.tail = FALSE)) < 1e-14)
duplicated <- rep(seq_along(z), each = 10)
repeated <- caliper_test(z[duplicated], clusters[duplicated], c, .1)
stopifnot(abs(repeated$cluster_se - test$cluster_se) < 1e-12,
          repeated$binom_p_greater_iid < test$binom_p_greater_iid)
stopifnot(caliper_test(10, "a", c, .1)$n_window == 0,
          is.na(caliper_test(c, "a", c, .1)$cluster_p_greater),
          caliper_test(c, "a", c, .1)$n_at_cutoff == 1)

## Trailing zeros are real reporting information, not floating-point digits.
d <- data.frame(yi = c(.08, .08, .074, .072), sei = c(.04, .04, .037, .036),
                metaID = c("a", "b", "c", "d"), md = c(2, 3, 3, 3), sd = c(2, 3, 3, 3))
prepared <- prepare_discontinuity_data(d, subgroup_col = NULL,
  estimate_digits_col = "md", se_digits_col = "sd")
stopifnot(identical(prepared$se_significand, c(4, 40, 37, 36)),
          identical(prepared$kp_keep_reported, c(FALSE, TRUE, TRUE, FALSE)),
          prepared$rounding_z_lower[3] >= c,
          prepared$rounding_z_lower[4] < c,
          abs(prepared$rounding_z_lower[2] - .0795 / .0405) < 1e-12,
          abs(prepared$rounding_z_upper[2] - .0805 / .0395) < 1e-12)
unknown <- prepare_discontinuity_data(d, subgroup_col = NULL)
stopifnot(!any(unknown$kp_keep_reported),
          !unknown$kp_keep_sensitivity[2],
          identical(inferred_decimal_places(c(.012, .04, 100, 1310)), c(3, 2, 0, 0)))
bad <- d; bad$md[1] <- 1
stopifnot(inherits(try(prepare_discontinuity_data(bad, subgroup_col = NULL,
  estimate_digits_col = "md", se_digits_col = "sd"), silent = TRUE), "try-error"))
bad <- d; bad$sei[1] <- 1e-16
stopifnot(inherits(try(prepare_discontinuity_data(bad, subgroup_col = NULL,
  estimate_digits_col = "md", se_digits_col = "sd"), silent = TRUE), "try-error"))
bad <- d; bad$sei[1:2] <- c(0, NA)
stopifnot(sum(!prepare_discontinuity_data(bad, subgroup_col = NULL)$valid) == 2)

## t statistics need their df; normal calibration must reproduce the same p.
tdata <- data.frame(t_value = qt(.975, c(5, 10, 100)), df = c(5, 10, 100), metaID = 1:3)
converted <- prepare_discontinuity_data(tdata, statistic_col = "t_value", df_col = "df", subgroup_col = NULL)
stopifnot(max(abs(converted$abs_z - c)) < 1e-12,
          max(abs(converted$p_two_sided - .05)) < 1e-12,
          all(converted$kp_keep_reported))
folder <- tempfile("t-discontinuity-")
t_run <- run_discontinuity_tests(tdata, statistic_col = "t_value", df_col = "df",
  cluster_col = NULL, subgroup_col = NULL, half_widths = .1,
  output_dir = folder, audit_dir = file.path(folder, "audit"))
stopifnot(nrow(t_run$tests) == 2L, all(t_run$tests$n_window == 3),
          all(t_run$tests$cluster_status == "no_cluster_column"))
unlink(folder, recursive = TRUE)

## A planted jump must be detected, with article-level inference available.
set.seed(11)
signal <- c(c - runif(80, .001, .099), c + runif(320, .001, .099))
degenerate <- caliper_test(signal, rep(1:80, 5), c, .1)
stopifnot(is.na(degenerate$cluster_p_greater),
          degenerate$cluster_status == "degenerate_cluster_variance")
detected <- caliper_test(signal, rep(1:80, each = 5), c, .1)
stopifnot(detected$share_above == .8, detected$cluster_p_greater < .001)

## End-to-end reported-precision mode; source() itself writes no outputs.
folder <- tempfile("discontinuity-output-")
res <- run_discontinuity_tests(d, output_dir = folder, half_widths = c(.1, .2),
  audit_dir = file.path(folder, "audit"),
  subgroup_col = NULL, estimate_digits_col = "md", se_digits_col = "sd")
stopifnot(nrow(res$tests) == 4L,
          all(file.exists(file.path(folder, c("caliper_tests.csv", "precision_coverage.csv",
                                             "rounding_risks.csv", "METHODS.txt")))))
unlink(folder, recursive = TRUE)

## Excel text and fixed decimal formats must preserve trailing zeros, unlike
## reading the same values as plain R doubles.
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Reported")
openxlsx::writeData(wb, "Reported", data.frame(yi = "0.0800", sei = "0.0400", metaID = 1))
openxlsx::writeData(wb, "Reported", data.frame(yi = .08, sei = .04, metaID = 2),
                    startRow = 3, colNames = FALSE)
openxlsx::addStyle(wb, "Reported", openxlsx::createStyle(numFmt = "0.000"),
                  rows = 3, cols = 1:2, gridExpand = TRUE)
path <- tempfile(fileext = ".xlsx")
openxlsx::saveWorkbook(wb, path)
entries <- excel_reporting_entries(path, c("yi", "sei"), 2L)
stopifnot(identical(entries$sei$decimals, c(4, 3)),
          identical(entries$sei$entry, c("0.0400", "0.040")),
          entries$sei$source[2] == "excel_fixed_decimal_format")
folder <- tempfile("xlsx-discontinuity-")
xlsx <- run_discontinuity_tests(path, output_dir = folder,
  audit_dir = file.path(folder, "audit"), subgroup_col = NULL, half_widths = .1)
stopifnot(all(xlsx$audit$kp_keep_reported),
          identical(xlsx$audit$reported_se_entry, c("0.0400", "0.040")),
          identical(xlsx$audit$se_significand, c(400, 40)))
unlink(c(path, folder), recursive = TRUE)
cat("PASS: exact counts, clustered dependence, rounding omission/bounds, t calibration and planted jump\n")
