source("scripts/waap_wls.R")
check_reference <- function(d, fit) {
  selected <- if (fit$method == "WAAP") d[fit$adequate, ] else d
  reference <- metafor::rma.mv(yi, vi, mods = ~1, method = "FE", test = "t", data = selected)
  test <- clubSandwich::coef_test(reference, vcov = clubSandwich::vcovCR(
    reference, cluster = droplevels(factor(selected$sID)), type = "CR2"))
  stopifnot(abs(fit$estimate - test$beta[1]) < 1e-10,
    abs(fit$standard_error - test$SE[1]) < 1e-10,
    abs(fit$p_value - test$p_Satt[1]) < 1e-10,
    abs(fit$df - test$df_Satt[1]) < 1e-8)
}
d <- data.frame(yi = c(1,1.1,1.2,.9,.7,1.3,.2,.4),
  vi = c(.01,.01,.04,.04,1,1,4,4), sID = c(1,1,2,2,3,3,4,5))
fit <- fit_waap_wls(d)
stopifnot(fit$method == "WAAP", fit$n_adequate == 4, fit$n_studies_fit == 2)
check_reference(d, fit)
negative <- d; negative$yi <- -negative$yi
neg <- fit_waap_wls(negative)
stopifnot(identical(fit$adequate,neg$adequate),abs(fit$estimate+neg$estimate)<1e-12)
d$yi <- d$yi/100
fallback <- fit_waap_wls(d)
stopifnot(fallback$method == "WLS fallback",fallback$n_adequate == 0)
check_reference(d,fallback)
one <- negative; one$sID[1:4] <- 1
single <- fit_waap_wls(one)
stopifnot(single$method == "WAAP",is.na(single$standard_error),
  single$inference_status == "fewer_than_two_primary_study_clusters")
cat("PASS: WAAP selection, signed effects, fallback and exact agreement with fixed-effect CR2/Satterthwaite inference.\n")
