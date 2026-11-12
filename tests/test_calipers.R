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
cat("PASS: pooled window counts, unchanged observed data, article bootstrap, t calibration, derived-data handoff and two-panel figure.\n")
