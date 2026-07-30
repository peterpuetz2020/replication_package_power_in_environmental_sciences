## Optional ESR workbook generation for the exploratory regressions.
## -----------------------------------------------------------------
## Estimating the ESR for each meta-analysis for regression analysis
## -----------------------------------------------------------------

## (1). Factual frequency per meta-analysis
kk <- NULL
jj <- NULL
f.per.meta <- list()

for (ii in 1:dim(summary(myDat))[1]) {

  fz <- abs(myDat[[ii]]$yi / myDat[[ii]]$sei)

  temp <- NULL
  for (a in 1:(length(p.grid.tab2)-1)) {
    temp[a] <- length(which(fz >= p.grid.tab2[a] & fz <= p.grid.tab2[a+1]))
  }

  kk[[ii]] <- temp
  jj[[ii]] <- myDat[[ii]]$cID[1]
  f.per.meta[[ii]] <- c(jj[[ii]], kk[[ii]]) #the first element is cluster id
}

#f.per.meta[[1]]

## (2). Counterfactual frequency per meta-analysis
source(here("scripts", "functions.R")) #functions

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cf.per.meta <- cf.disagg(dat=myDat, z.grid=p.grid.tab)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) #less than a minute

#5%
cID <- NULL
esr.all <- NULL
esr.sig <- NULL
esr.all.count <- NULL
esr.sig.count <- NULL
dif <- NULL
tot.all <- NULL #total p-values
tot.sig <- NULL #total significant p-values
matESR <- data.frame(cID=NA,esr.all=NA,esr.sig=NA,esr.all.count=NA,esr.sig.count=NA,dif=NA,tot.all=NA,tot.sig=NA)

for (i in 1:dim(summary(f.per.meta))[1]) {

  cID[i] <- f.per.meta[[i]][1] #cluster id
  esr.all[i] <- (sum(as.numeric(f.per.meta[[i]][12:14])) - sum(as.numeric(cf.per.meta[[i]][11:13]))) / sum(as.numeric(f.per.meta[[i]][2:14])) #as a share of all p-values
  esr.sig[i] <- (sum(as.numeric(f.per.meta[[i]][12:14])) - sum(as.numeric(cf.per.meta[[i]][11:13]))) / sum(as.numeric(f.per.meta[[i]][12:14])) #as a share of only signficant p-values

  # "esr.all.count" = "esr.sig.count" unless their denominator is different
  esr.all.count[i] <- (sum(as.numeric(f.per.meta[[i]][12:14])) - sum(as.numeric(cf.per.meta[[i]][11:13]))) #/ sum(f.per.meta[[i]])
  esr.sig.count[i] <- (sum(as.numeric(f.per.meta[[i]][12:14])) - sum(as.numeric(cf.per.meta[[i]][11:13]))) #/ sum(f.per.meta[[i]][10:13])

  #test
  dif[i] <- (sum(as.numeric(f.per.meta[[i]][12:14])) - sum(as.numeric(cf.per.meta[[i]][11:13])))
  tot.all[i] <- sum(as.numeric(f.per.meta[[i]][2:14]))
  tot.sig[i] <- sum(as.numeric(f.per.meta[[i]][12:14]))

  matESR[i,1] <- cID[i]
  matESR[i,2] <- esr.all[i]
  matESR[i,3] <- esr.sig[i]
  matESR[i,4] <- esr.all.count[i]
  matESR[i,5] <- esr.sig.count[i]
  matESR[i,6] <- dif[i]
  matESR[i,7] <- tot.all[i]
  matESR[i,8] <- tot.sig[i]
}

head(matESR)
tail(matESR)
dim(matESR) #704x8

## NOTE: esr.sig==-Inf if a meta-analysis has no single statistical significant effect sizes.
## Thus, we replace the values of esr.sig & esr.sig.count by zero.
matESR$esr.sig.count[matESR$esr.sig.count<0 & matESR$esr.sig==-Inf] <- 0
matESR$esr.sig[matESR$esr.sig==-Inf] <- 0
write.xlsx(matESR, here("results","main","esr05_pet_peese_rstandard_half meta-average_704.xlsx"), overwrite=TRUE)
#write.xlsx(matESR, here("results","main","esr05_pet_peese_rstandard_half meta-average_0.25xtau2_704.xlsx"), overwrite=TRUE)
#write.xlsx(matESR, here("results","main","esr05_pet_peese_rstandard_full meta-average_704.xlsx"), overwrite=TRUE)
summary(matESR$esr.sig)
