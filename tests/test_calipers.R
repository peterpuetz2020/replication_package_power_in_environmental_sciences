source("scripts/functions.R")
c <- qnorm(.975)
z <- c+ifelse(c(0,1,1,0,0,1,1,1),.05,-.05)
clusters <- c("a","a","a","b","b","c","c","d")
counts <- caliper_test(z,clusters,c,.1)
stopifnot(counts$n_window==8,counts$n_above_including_ties==5,
  counts$share_above==5/8,counts$n_clusters==4,
  caliper_test(c(1.76,2.16),cutoff=1.96,half_width=.2)$n_window==2,
  caliper_test(10,"a",c,.1)$n_window==0)

## Duplicating observations within articles preserves bootstrap dependence.
z_boot <- rep(c(c+.05,c-.05),c(30,10))
ids <- rep(1:8,each=5)
set.seed(89); previous_seed <- .Random.seed
b <- cluster_bootstrap_calipers(z_boot,ids,.1,iterations=1000,seed=123)
stopifnot(identical(.Random.seed,previous_seed))
repeated <- cluster_bootstrap_calipers(rep(z_boot,each=3),rep(ids,each=3),.1,
  iterations=1000,seed=123)
stopifnot(identical(b,repeated),
  abs(b$bootstrap_ci_lower-qbinom(.025,8,.75)/8)<=1/8,
  abs(b$bootstrap_ci_upper-qbinom(.975,8,.75)/8)<=1/8)
stopifnot(cluster_bootstrap_calipers(z,NULL,.1)$bootstrap_status=="no_cluster_column",
  cluster_bootstrap_calipers(z,rep("a",8),.1)$bootstrap_status=="fewer_than_two_clusters")

## All finite observed estimates remain eligible, regardless of decimal precision.
d <- data.frame(yi=c(.08,.08,.074,.072),vi=c(.0016,.0016,.001369,.001296),
  sei=999,metaID=1:4,subfield=c("A","B","A","B"))
analysis <- compute_caliper_analysis(d)
stopifnot(analysis$coverage$n_valid==nrow(d),
  identical(analysis$audit$standard_error,sqrt(d$vi)),
  nrow(analysis$tests)==4L,all(analysis$tests$group=="All"),
  all(analysis$tests$bootstrap_iterations==1000),
  !any(grepl("(^p_|_p_|p_value|holm|binom|cluster_ci|precision|kp_)",names(analysis$tests))))
bad <- d; bad$vi[1:3] <- c(0,-1,NA)
stopifnot(sum(!prepare_caliper_data(bad)$valid)==3)
precise <- d; precise$vi[1] <- 1.234567890123456e-10
stopifnot(identical(prepare_caliper_data(precise)$variance,precise$vi))
stopifnot(inherits(try(compute_caliper_analysis(d,half_widths=.3),silent=TRUE),"try-error"))

## Supplied t degrees of freedom calibrate to the same normal threshold.
tdata <- data.frame(t_value=qt(.975,c(5,10,100)),df=c(5,10,100),metaID=1:3)
converted <- compute_caliper_analysis(tdata,statistic_col="t_value",df_col="df")
stopifnot(max(abs(converted$audit$abs_z-c))<1e-12,
  all(converted$tests$n_window==3))

## A planted jump is visible in the cluster-bootstrap confidence interval.
set.seed(11)
signal <- c(c-runif(80,.001,.099),c+runif(320,.001,.099))
detected <- cluster_bootstrap_calipers(signal,rep(1:80,each=5),.1,cutoff=c)
stopifnot(detected$bootstrap_ci_lower>.5)

folder <- tempfile("caliper-derived-")
save_caliper_analysis(analysis,folder)
saved <- openxlsx::read.xlsx(file.path(folder,"caliper_tests.xlsx"),sheet="Results")
stopifnot(isTRUE(all.equal(saved$share_above,analysis$tests$share_above)),
  file.exists(file.path(folder,"caliper_coverage.xlsx")),
  identical(readRDS(file.path(folder,"caliper_analysis.rds"))$settings,analysis$settings))
plots <- make_caliper_plots(analysis)
stopifnot(length(plots)==2L,
  sum(plots[[1]]$data$count)==analysis$tests$n_window[which.max(analysis$tests$half_width)])
unlink(folder,recursive=TRUE)

## Check yi against sei OR vi separately, then use the matching decimal count.
## In particular, a match must not fail because the other input is coarser.
metadata <- function(values,digits) data.frame(entry=as.character(values),decimals=digits,
  source="supplied_original_entries")
paired <- data.frame(yi=rep(.2,4),vi=rep(.01,4),sei=rep(.1,4),metaID=1:4)
entries <- list(yi=metadata(paired$yi,c(4,1,4,4)),
  vi=metadata(paired$vi,c(4,4,2,4)),sei=metadata(paired$sei,c(4,4,4,3)))
previous_scipen <- getOption("scipen")
adjusted <- compute_caliper_analysis(paired,reporting_entries=entries)
negative <- paired; negative$yi <- -paired$yi
negative_audit <- prepare_rounding_calipers(prepare_caliper_data(negative),entries)
stopifnot(identical(adjusted$audit$kp_mu_decimals,c(4,1,4,4)),
  identical(getOption("scipen"),previous_scipen),
  identical(adjusted$audit$kp_sigma_decimals,c(4,4,4,4)),
  identical(adjusted$audit$kp_common_precision,c(TRUE,FALSE,TRUE,TRUE)),
  identical(adjusted$audit$kp_precision_match,c("both","none","sei","vi")),
  identical(adjusted$audit$kp_se_significand,c(1000,1000,1000,1000)),
  identical(adjusted$audit$kp_rule_applied,c(TRUE,FALSE,TRUE,TRUE)),
  all(adjusted$audit$kp_keep),
  adjusted$audit$kp_status[2]=="retained_no_precision_match",
  adjusted$coverage$n_unequal_reporting_precision==1,
  adjusted$coverage$n_matches_sei_only==1,
  adjusted$coverage$n_matches_vi_only==1,
  adjusted$coverage$n_matches_both==1,
  adjusted$coverage$n_omitted_by_precision_rule==0,
  !"kp_reference_z_lower" %in% names(adjusted$audit),
  identical(compute_caliper_analysis(paired,cutoff=2.1,reporting_entries=entries)$audit$kp_keep,
    adjusted$audit$kp_keep),
  nrow(adjusted$tests)==12L,
  all(adjusted$tests$n_retained[adjusted$tests$method=="kranz_puetz"]==4),
  identical(adjusted$audit$variance,paired$vi),
  identical(adjusted$audit$standard_error,sqrt(paired$vi)),
  identical(adjusted$audit$kp_keep,negative_audit$kp_keep),
  identical(adjusted$audit$rounding_abs_z_lower,negative_audit$rounding_abs_z_lower),
  identical(adjusted$audit$rounding_abs_z_upper,negative_audit$rounding_abs_z_upper),
  max(abs(adjusted$audit$rounding_abs_z_lower-
    pmax(0,(paired$yi-.5*10^(-entries$yi$decimals))/
      (sqrt(paired$vi)+.5*10^(-adjusted$audit$kp_sigma_decimals)))))<1e-12)

## A missing nonmatching column must not invalidate a match in the other one.
without_vi <- entries; without_vi$vi$decimals[] <- NA_real_
without_sei <- entries; without_sei$sei$decimals[] <- NA_real_
stopifnot(identical(prepare_rounding_calipers(prepare_caliper_data(paired),without_vi)$kp_rule_applied,
    c(TRUE,FALSE,TRUE,FALSE)),
  identical(prepare_rounding_calipers(prepare_caliper_data(paired),without_sei)$kp_rule_applied,
    c(TRUE,FALSE,FALSE,TRUE)),
  all(prepare_rounding_calipers(prepare_caliper_data(paired),without_vi)$kp_keep),
  all(prepare_rounding_calipers(prepare_caliper_data(paired),without_sei)$kp_keep))

## Recover the published common-precision boundary (36 omitted, 37 retained).
common <- data.frame(yi=c(.072,.074),vi=c(.036,.037)^2,metaID=1:2)
common_entries <- setNames(lapply(c("yi","vi","sei"),function(column)
  metadata(c(0,0),c(3,3))),c("yi","vi","sei"))
common_audit <- prepare_rounding_calipers(prepare_caliper_data(common),common_entries)
stopifnot(identical(common_audit$kp_se_significand,c(36,37)),
  identical(common_audit$kp_keep,c(FALSE,TRUE)),
  identical(common_audit$kp_status,c("omitted_by_precision_rule","retained")))
missing <- entries; missing$yi$decimals[] <- NA_real_
empty <- compute_caliper_analysis(paired,reporting_entries=missing)
stopifnot(all(empty$tests$n_window[empty$tests$method=="kranz_puetz"]==4),
  all(empty$audit$kp_keep),!any(empty$audit$kp_rule_applied),
  all(empty$audit$kp_status=="retained_unknown_precision"),
  all(empty$audit$rounding_filter_keep),
  all(empty$audit$rounding_filter_status=="retained_unknown_precision"),
  empty$coverage$n_unknown_reporting_precision==4,
  identical(empty$tests[empty$tests$method=="empirical",names(analysis$tests)],
    compute_caliper_analysis(paired)$tests))

## Remove only matched failures; the same low SE significand is not an omission
## when coefficient precision does not match. Invalid observations stay excluded.
conditional <- data.frame(yi=c(.072,.074,.072,.074),vi=c(.036,.037,.036,.037)^2,metaID=1:4)
conditional_entries <- list(yi=metadata(conditional$yi,c(3,3,2,2)),
  vi=metadata(conditional$vi,rep(3,4)),sei=metadata(sqrt(conditional$vi),rep(3,4)))
selected <- compute_caliper_analysis(conditional,reporting_entries=conditional_entries)
stopifnot(identical(selected$audit$kp_keep,c(FALSE,TRUE,TRUE,TRUE)),
  identical(selected$audit$kp_rule_applied,c(TRUE,TRUE,FALSE,FALSE)),
  selected$coverage$n_rule_applied==2,
  selected$coverage$n_retained_after_kp_rule==1,
  selected$coverage$n_retained_without_kp_rule==2,
  selected$coverage$n_retained_kranz_puetz==3,
  selected$coverage$n_omitted_by_precision_rule==1,
  all(selected$tests$n_window[selected$tests$method=="kranz_puetz"]==3))
conditional$vi[3] <- 0
stopifnot(!prepare_rounding_calipers(prepare_caliper_data(conditional),conditional_entries)$kp_keep[3])

## The independent third sample uses min(vi,sei) precision even when KP matches
## the more precise input. Exclude only intervals spanning the cutoff, on either
## sign of the coefficient; keep intervals certainly below or certainly above.
ambiguous <- data.frame(yi=c(.196,.18,.22,.196,-.196,0),vi=.01,metaID=1:6)
ambiguity_entries <- list(yi=metadata(ambiguous$yi,c(4,2,2,3,3,0)),
  vi=metadata(ambiguous$vi,c(2,4,4,2,2,4)),
  sei=metadata(sqrt(ambiguous$vi),c(4,3,2,4,4,3)))
sensitivity <- compute_caliper_analysis(ambiguous,reporting_entries=ambiguity_entries)
a <- sensitivity$audit
stopifnot(identical(a$rounding_filter_sigma_decimals,c(2,3,2,2,2,3)),
  a$kp_sigma_decimals[1]==4,a$rounding_filter_sigma_decimals[1]==2,
  identical(a$rounding_filter_straddles_cutoff,c(TRUE,FALSE,FALSE,TRUE,TRUE,TRUE)),
  identical(a$rounding_filter_keep,c(FALSE,TRUE,TRUE,FALSE,FALSE,FALSE)),
  !a$kp_keep[3],a$rounding_filter_keep[3],
  sensitivity$coverage$n_omitted_straddling_cutoff==4,
  sensitivity$coverage$n_retained_rounding_ambiguity==2,
  all(sensitivity$rounding_risks$n_potential_significance_switch[
    sensitivity$rounding_risks$method=="rounding_ambiguity"]==0))
expected_lo <- pmax(0,abs(ambiguous$yi)-.5*10^(-a$reported_yi_decimals))/
  (sqrt(ambiguous$vi)+.5*10^(-a$rounding_filter_sigma_decimals))
expected_hi <- (abs(ambiguous$yi)+.5*10^(-a$reported_yi_decimals))/
  pmax(0,sqrt(ambiguous$vi)-.5*10^(-a$rounding_filter_sigma_decimals))
stopifnot(max(abs(a$rounding_filter_abs_z_lower-expected_lo))<1e-12,
  max(abs(a$rounding_filter_abs_z_upper-expected_hi))<1e-12,
  prepare_rounding_ambiguity_filter(a,cutoff=a$rounding_filter_abs_z_upper[2])$
    rounding_filter_straddles_cutoff[2],
  !prepare_rounding_ambiguity_filter(a,cutoff=a$rounding_filter_abs_z_lower[3])$
    rounding_filter_straddles_cutoff[3])

## Read precision from workbook entries/formatting, including trailing zeros,
## instead of using the digits of computed sqrt(vi) or formatted R doubles.
stopifnot(identical(reported_decimal_places(c("0.0400","1.20e-3","100","1e3",NA)),
  c(4,5,0,-3,NA_real_)))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,"Reported")
openxlsx::writeData(wb,"Reported",data.frame(yi="0.0800",vi="0.001600",sei="0.040"))
openxlsx::writeData(wb,"Reported",data.frame(yi=.08,vi=.0016,sei=.04),startRow=3,colNames=FALSE)
openxlsx::addStyle(wb,"Reported",openxlsx::createStyle(numFmt="0.00000"),rows=3,cols=1:3,gridExpand=TRUE)
path <- tempfile(fileext=".xlsx")
openxlsx::saveWorkbook(wb,path)
xentries <- excel_reporting_entries(path,c("yi","vi","sei"),2L)
stopifnot(identical(xentries$yi$decimals,c(4,5)),
  identical(xentries$vi$decimals,c(6,5)),
  identical(xentries$sei$decimals,c(3,5)),
  xentries$yi$entry[1]=="0.0800",xentries$sei$source[2]=="excel_fixed_decimal_format")
xdata <- as.data.frame(lapply(readxl::read_excel(path),as.numeric)); xdata$metaID <- 1:2
xanalysis <- compute_caliper_analysis(xdata,workbook_path=path)
stopifnot(identical(xanalysis$audit$kp_sigma_decimals,c(3,5)),
  all(xanalysis$audit$kp_keep),
  identical(xanalysis$audit$kp_rule_applied,c(FALSE,TRUE)))
unlink(path)

folder <- tempfile("caliper-adjusted-")
save_caliper_analysis(adjusted,folder)
saved <- openxlsx::read.xlsx(file.path(folder,"caliper_tests.xlsx"),sheet="Results")
plots <- make_caliper_plots(adjusted)
stopifnot(nrow(saved)==12L,length(plots)==6L,
  nrow(openxlsx::read.xlsx(file.path(folder,"caliper_rounding_risks.xlsx"),sheet="Results"))==12L,
  sum(plots[[1]]$data$count)==4L,sum(plots[[3]]$data$count)==4L,
  sum(plots[[5]]$data$count)==adjusted$tests$n_window[
    adjusted$tests$method=="rounding_ambiguity" & adjusted$tests$half_width==.2],
  identical(plots[[2]]$coordinates$limits$y,c(.45,.65)),
  identical(plots[[4]]$coordinates$limits$y,c(.45,.65)),
  identical(plots[[6]]$coordinates$limits$y,c(.45,.65)),
  !any(grepl("p_value|binom|holm|cluster_ci",names(saved))))
unlink(folder,recursive=TRUE)
cat("PASS: conditional KP omissions, retained nonmatches, independent min-precision cutoff filter, original thresholds, bootstrap and six-panel Figure S6.\n")
