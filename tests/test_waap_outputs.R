## Integration checks after fitting WAAP and rebuilding Figure S5.
source("scripts/functions.R")
summary <- openxlsx::read.xlsx("data/derived_data/WAAP_WLS_estimates.xlsx",sheet="Estimates")
stopifnot(!anyDuplicated(paste(summary$cID,summary$variant)))
for (variant in c("outliers_removed","all_data")) {
  suffix <- if (variant=="all_data") "_all_data" else ""
  files <- list.files(paste0("data/derived_data/waap_wls",suffix),full.names=TRUE,pattern="\\.rds$")
  stopifnot(length(files)==sum(summary$variant==variant))
  for (file in files) {
    d <- readRDS(file)
    original <- readRDS(file.path(paste0("data/derived_data/multilevel_random",suffix),basename(file)))
    stopifnot(identical(d$yi,original$yi),identical(d$vi,original$vi),
      identical(d$sID,original$sID),identical(d$eID,original$eID),all(d$tau2==0),
      all(d$heterogeneity_source=="not_estimated_zero_assumed"))
    mu <- weighted.mean(d$yi,1/d$vi)
    lambda <- abs(mu)/sqrt(d$vi)
    adequate <- pnorm(-qnorm(.975)-lambda)+pnorm(lambda-qnorm(.975))>=.8
    keep <- if(sum(adequate)>=2) adequate else rep(TRUE,nrow(d))
    stopifnot(all(d$waap_n_adequate==sum(adequate)),all(d$waap_n_fit==sum(keep)),
      max(abs(d$GE-weighted.mean(d$yi[keep],1/d$vi[keep])))<1e-9)
  }
}
rows <- readRDS("data/derived_data/Figure_S5_inputs.rds")$all_esr_results
waap <- subset(rows,estimator=="waap_wls" & measure=="ESR_{0.05}^{sig}")
stopifnot(nrow(waap)==4L,setequal(waap$meta_average_multiplier,c(.25,.5,.75,1)),
  all(waap$heterogeneity_multiplier==0),
  all(rows$heterogeneity_multiplier[rows$estimator=="fixed_effect"]==0),
  setequal(unique(rows$estimator),c("pet_peese","multilevel_random","fixed_effect","waap_wls")))
d <- do.call(rbind,lapply(list.files("data/derived_data/waap_wls",full.names=TRUE,pattern="\\.rds$"),readRDS))
breaks <- c(0,qnorm(c(.9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.01,.001)/2,lower.tail=FALSE),Inf)
for(i in seq_len(nrow(waap))) {
  row <- waap[i,]; filtered <- d; filtered$GE <- filtered$GE*row$meta_average_multiplier
  filtered <- filter_counterfactual_data(filtered,row$heterogeneity_multiplier)
  z <- abs(filtered$yi/sqrt(filtered$vi))
  observed <- vapply(seq_len(13),function(j)sum(z>=breaks[j]&z<=breaks[j+1]),numeric(1))
  stem <- paste0(row$setup_label,"_waap_wls_outliers_removed.rds")
  point <- readRDS(paste0("data/derived_data/p_tab_",stem))$result
  draws <- readRDS(paste0("data/derived_data/p_tab_ci_",stem))$result
  expected_point <- round(sum((observed-point)[11:13])/sum(observed[11:13]),3)
  expected_ci <- round(quantile(sum(observed[11:13])/sum(observed[11:13])-
    rowSums(draws[[3]][,11:13]),c(.025,.975)),3)
  stored_ci <- as.numeric(strsplit(gsub("[\\[\\]]","",row$confidence_interval,perl=TRUE),", ")[[1]])
  stopifnot(nrow(draws[[1]])==1000L,abs(as.numeric(row$estimate)-expected_point)<1e-12,
    max(abs(stored_ci-expected_ci))<1e-12)
}
cat("PASS: shared samples and WAAP fits retained; fixed effects and WAAP restricted to zero heterogeneity; four WAAP points/1000-draw intervals.\n")
