## Optional PET-PEESE model fitting and per-meta-analysis data generation.
## Import the data
meta <- read_excel(here("data","MasterData.xlsx"))
dim(meta)

## Drop trouble making meta-analyses (non-convergence issue)
meta <- meta %>%
  dplyr::filter(cID!="259-1",cID!="259-2",cID!="2204-2",cID!="850-1")

myDat <- list()
jj <- 1
for (i in unique(meta$cID)){
  myDat[[i]] <- meta[which(meta$cID==i), ]
  jj <- jj + 1
}

## check
mss <- NULL
for (i in 1:dim(summary(myDat))[1]) {
  mss[i] <- length(myDat[[i]]$sei)
}
length(mss) #704
sum(mss)    #64711


## -----------------------------------------------------------------------------

## Implementing extended PET-PEESE to estimate genuine effect
## This approach takes into account the non-independence of
## effect sizes & correct for publication bias (Nakagawa et al. 2021)

## Note that influential observations are excluded using standardized residuals (aka internally studentized residuals)
## Trouble makers when fitting PEESE ("Error: Model matrix not of full rank
## (the columns are perfectly correlated). Cannot fit model ("259-1"  "259-2"  "2204-2" "850-1"). For this, these meta-analyses are dropped.

petStud <- petStud.CR2 <- petStud.rob <- petGE1Stud <- pet.pv1.rob.stud <- petpeeseGE2 <- petpeese.pv2.rob <- list()
peeseStud <- peeseStud.CR2 <- peeseStud.rob <- peeseGE1Stud <- peese.pv1.rob.stud <- list()
rstud <- n_outliers <- outliers_dropped <- list()
petMod <- petMod.CR2 <- petMod.rob <- pet.pv1.rob <- petGE1 <- list()
petpeeseGE2 <- petpeese.pv2 <- peeseMod <- peeseMod.CR2 <- peeseMod.rob <- peeseGE1 <- peese.pv1.rob <- list()
pet.tau2.stud <- peese.tau2.stud <- petpeese.tau2 <- list()
pet.tau2.nonstud <- peese.tau2.nonstud <- list()
pet.pv2.rob <- pet.pv2.rob.stud <- peese.pv2.rob.stud <- pet.pv2.rob.nonstud <- peese.pv2.rob.nonstud <- petpeese.pv3 <- list()
pet.isq.stud <- petpeese.isq <- peese.isq.stud <- list()
pet.isq.nonstud <- peese.isq.nonstud <- list()
metaID <- cID <- sID <- eID <- yi <- vi <- list()
etype <- guide <- prere <- subfd <- sdesn <- list()
pkg <- c("metafor","dplyr","clubSandwich","orchaRd", "here")

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
pet_peese_rstud <- foreach(i = 1:length(summary(myDat)[,1]),.packages=pkg,.combine="rbind") %dopar% {

  petMod[[i]] <- rma.mv(yi, vi, data=myDat[[i]], mods = ~1 + sei, method="REML", test="t", random=list(~1 | eID, ~1 | sID), control=list(rel.tol=1e-8))
  rstud[[i]] <- as.data.frame(rstandard.rma.mv(petMod[[i]]))
  n_outliers[[i]] <- length(which(abs(rstud[[i]]$resid) >= 3)) #reference: page 97 of Gareth, Daniela, Trevor, & Robert (2013)

  if(n_outliers[[i]] >= 1){  #when at least one influential observation is identified
    outliers_dropped[[i]] <- myDat[[i]] %>%
      cbind(rstud[[i]]) %>%
      filter(rstud[[i]]$resid < 3)
    petStud[[i]] <- rma.mv(yi, vi, data=outliers_dropped[[i]], mods=~1 + sei, random=list(~1 | eID, ~1 | sID), method="REML", test="t", control=list(rel.tol=1e-8))
    petStud.CR2[[i]] <- vcovCR(petStud[[i]], type="CR2")
    petStud.rob[[i]] <- coef_test(petStud[[i]], vcov=petStud.CR2[[i]])
    pet.pv1.rob.stud[[i]] <- as.numeric(petStud.rob[[i]]$p_Satt[1])
    pet.pv2.rob.stud[[i]] <- as.numeric(petStud.rob[[i]]$p_Satt[2]) #test small-study effect

    if (pet.pv1.rob.stud[[i]] > 0.1) {
      petGE1Stud[[i]] <- round(as.numeric(petStud.rob[[i]]$beta[1]), 3)
      pet.pv1.rob.stud[[i]] <- as.numeric(petStud.rob[[i]]$p_Satt[1])
      pet.pv2.rob.stud[[i]] <- as.numeric(petStud.rob[[i]]$p_Satt[2]) #test small-study effect
      pet.tau2.stud[[i]] <- round(petStud[[i]]$sigma2[2], 3)          #variance of true effect (between-study variance)
      pet.isq.stud[[i]] <- round(as.numeric(i2_ml(petStud[[i]],method="matrix")[1]), 2) #total heterogeneity (I^2)

      petpeeseGE2[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), petGE1Stud[[i]])
      petpeese.pv2[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), pet.pv1.rob.stud[[i]])
      petpeese.pv3[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), pet.pv2.rob.stud[[i]])
      petpeese.tau2[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), pet.tau2.stud[[i]])
      petpeese.isq[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), pet.isq.stud[[i]])
    } else {
      peeseStud[[i]] <-  rma.mv(yi, vi, data=outliers_dropped[[i]], mods=~1 + vi, random=list(~1 | eID, ~1 | sID), method="REML", test="t", control=list(rel.tol=1e-8))
      peeseStud.CR2[[i]] <- vcovCR(peeseStud[[i]], type="CR2")
      peeseStud.rob[[i]] <- coef_test(peeseStud[[i]], vcov=peeseStud.CR2[[i]])
      peeseGE1Stud[[i]] <- round(as.numeric(peeseStud.rob[[i]]$beta[1]), 3)
      peese.pv1.rob.stud[[i]] <- as.numeric(peeseStud.rob[[i]]$p_Satt[1])
      peese.pv2.rob.stud[[i]] <- as.numeric(peeseStud.rob[[i]]$p_Satt[2]) #test small-study effect

      peese.tau2.stud[[i]] <- round(peeseStud[[i]]$sigma2[2], 3)
      peese.isq.stud[[i]] <- round(as.numeric(i2_ml(peeseStud[[i]],method="matrix")[1]), 2)

      petpeeseGE2[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), peeseGE1Stud[[i]])
      petpeese.pv2[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), peese.pv1.rob.stud[[i]])
      petpeese.pv3[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), peese.pv2.rob.stud[[i]])
      petpeese.tau2[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), peese.tau2.stud[[i]])
      petpeese.isq[[i]] <- replicate(n=length(outliers_dropped[[i]]$yi), peese.isq.stud[[i]])
    }

    metaID[[i]] <- outliers_dropped[[i]]$metaID
    cID[[i]] <- outliers_dropped[[i]]$cID
    sID[[i]] <- outliers_dropped[[i]]$sID
    eID[[i]] <- outliers_dropped[[i]]$eID
    yi[[i]] <- outliers_dropped[[i]]$yi
    vi[[i]] <- outliers_dropped[[i]]$vi
    etype[[i]] <- outliers_dropped[[i]]$estype
    guide[[i]] <- outliers_dropped[[i]]$guide
    prere[[i]] <- outliers_dropped[[i]]$prereg
    subfd[[i]] <- outliers_dropped[[i]]$subfield
    sdesn[[i]] <- outliers_dropped[[i]]$sdesign

    df_outliers_dropped <- cbind.data.frame(metaID[[i]],cID[[i]],sID[[i]],eID[[i]],yi[[i]],vi[[i]],petpeeseGE2[[i]],petpeese.tau2[[i]],petpeese.isq[[i]],petpeese.pv2[[i]],petpeese.pv3[[i]],etype[[i]],guide[[i]],prere[[i]],subfd[[i]],sdesn[[i]])
    colnames(df_outliers_dropped) <- c("metaID","cID","sID","eID","yi","vi","GE","tau2","isq","sig_overall","small_study_effect_pval","etype","guide","prere","subfd","sdesn")

    saveRDS(df_outliers_dropped, file=here("results","main","pet_peese_rstandard",paste0("meta_", outliers_dropped[[i]]$cID[1],".rds")))
    return(df_outliers_dropped)

  } else { #when no influential observation is identified
    petMod[[i]] <- rma.mv(yi, vi, data=myDat[[i]],mods = ~1 + sei, random=list(~1 | eID, ~1 | sID), method="REML", test="t", control=list(rel.tol=1e-8))
    petMod.CR2[[i]] <- vcovCR(petMod[[i]], type="CR2")
    petMod.rob[[i]] <- coef_test(petMod[[i]], vcov=petMod.CR2[[i]])
    pet.pv1.rob[[i]] <- as.numeric(petMod.rob[[i]]$p_Satt[1])
    pet.pv2.rob.nonstud[[i]] <- as.numeric(petMod.rob[[i]]$p_Satt[2]) #test small-study effect

    if (pet.pv1.rob[[i]] > 0.1) {
      petGE1[[i]] <- round(as.numeric(petMod.rob[[i]]$beta[1]), 3)
      pet.pv1.rob[[i]] <- as.numeric(petMod.rob[[i]]$p_Satt[1])
      pet.pv2.rob.nonstud[[i]] <- as.numeric(petMod.rob[[i]]$p_Satt[2])
      pet.tau2.nonstud[[i]] <- round(petMod[[i]]$sigma2[2], 3)
      pet.isq.nonstud[[i]] <- round(as.numeric(i2_ml(petMod[[i]],method="matrix")[1]), 2)

      petpeeseGE2[[i]] <- replicate(n=length(myDat[[i]]$yi), petGE1[[i]])
      petpeese.pv2[[i]] <- replicate(n=length(myDat[[i]]$yi), pet.pv1.rob[[i]])
      petpeese.pv3[[i]] <- replicate(n=length(myDat[[i]]$yi), pet.pv2.rob.nonstud[[i]])
      petpeese.tau2[[i]] <- replicate(n=length(myDat[[i]]$yi), pet.tau2.nonstud[[i]])
      petpeese.isq[[i]] <- replicate(n=length(myDat[[i]]$yi), pet.isq.nonstud[[i]])
    } else {
      peeseMod[[i]] <-  rma.mv(yi, vi, data=myDat[[i]], mods = ~1 + vi, random=list(~1 | eID, ~1 | sID), method="REML", test="t", control=list(rel.tol=1e-8))
      peeseMod.CR2[[i]] <- vcovCR(peeseMod[[i]], type="CR2")
      peeseMod.rob[[i]] <- coef_test(peeseMod[[i]], vcov=peeseMod.CR2[[i]])
      peeseGE1[[i]] <- round(as.numeric(peeseMod.rob[[i]]$beta[1]), 3)
      peese.pv1.rob[[i]] <- as.numeric(peeseMod.rob[[i]]$p_Satt[1])
      peese.pv2.rob.nonstud[[i]] <- as.numeric(peeseMod.rob[[i]]$p_Satt[2])
      peese.tau2.nonstud[[i]] <- round(peeseMod[[i]]$sigma2[2], 3)
      peese.isq.nonstud[[i]] <- round(as.numeric(i2_ml(peeseMod[[i]],method="matrix")[1]), 2)

      petpeeseGE2[[i]] <- replicate(n=length(myDat[[i]]$yi), peeseGE1[[i]])
      petpeese.pv2[[i]] <- replicate(n=length(myDat[[i]]$yi), peese.pv1.rob[[i]])
      petpeese.pv3[[i]] <- replicate(n=length(myDat[[i]]$yi), peese.pv2.rob.nonstud[[i]])
      petpeese.tau2[[i]] <- replicate(n=length(myDat[[i]]$yi), peese.tau2.nonstud[[i]])
      petpeese.isq[[i]] <- replicate(n=length(myDat[[i]]$yi), peese.isq.nonstud[[i]])
    }

    metaID[[i]] <- myDat[[i]]$metaID
    cID[[i]] <- myDat[[i]]$cID
    sID[[i]] <- myDat[[i]]$sID
    eID[[i]] <- myDat[[i]]$eID
    yi[[i]] <- myDat[[i]]$yi
    vi[[i]] <- myDat[[i]]$vi
    etype[[i]] <- myDat[[i]]$estype
    guide[[i]] <- myDat[[i]]$guide
    prere[[i]] <- myDat[[i]]$prereg
    subfd[[i]] <- myDat[[i]]$subfield
    sdesn[[i]] <- myDat[[i]]$sdesign

    df <- cbind.data.frame(metaID[[i]],cID[[i]],sID[[i]],eID[[i]],yi[[i]],vi[[i]],petpeeseGE2[[i]],petpeese.tau2[[i]],petpeese.isq[[i]],petpeese.pv2[[i]],petpeese.pv3[[i]],etype[[i]],guide[[i]],prere[[i]],subfd[[i]],sdesn[[i]])
    colnames(df) <- c("metaID","cID","sID","eID","yi","vi","GE","tau2","isq","sig_overall","small_study_effect_pval","etype","guide","prere","subfd","sdesn")
    saveRDS(df, file=here("results","main","pet_peese_rstandard",paste0("meta_", myDat[[i]]$cID[1],".rds")))
    return(df)
  }
}
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) #takes about 1 hrs
