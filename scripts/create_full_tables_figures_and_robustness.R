## Complete legacy manuscript, supplementary, and robustness workflow.
## This script also recreates expensive counterfactual data; source it only when
## those supplied derived data and all legacy outputs need to be regenerated.
## ---

## Load data set generated (rstandard)
pps_rstandard <- list.files(here("results","main","pet_peese_rstandard"),pattern=".rds",full.names=TRUE) %>%
  map_dfr(readRDS)
dim(pps_rstandard) #63956x16

pps_rstandard$sei <- sqrt(pps_rstandard$vi)
pps_rstandard$sse_yn <- ifelse(pps_rstandard$small_study_effect_pval <= 0.05,"yes","no")

## Power analysis (using multilevel PET-PEESE - Table 3)
alpha <- 0.05
q1 <- qnorm(1-alpha/2)  # 1.96
q2 <- qnorm(alpha/2)    # -1.96

pps_rstandard$GE <- meta_average_multiplier * pps_rstandard$GE #default uses half meta-average

pps_rstandard$lambda <- abs(pps_rstandard$GE)/pps_rstandard$sei  #non-centrality parameter
pps_rstandard$power <- 1- pnorm(q1 - pps_rstandard$lambda) + pnorm(q2 - pps_rstandard$lambda) #power of each test
pps_rstandard$yn80 <- ifelse(pps_rstandard$power >= 0.8, "yes","no")
100*round(prop.table(table(pps_rstandard$yn80)), 4) #no=81.71, yes=18.29
xx <- 0.18 #using half meta-average
#xx <- 0.33 #using full meta-average
#xx <- 0.25 #null MAs dropped & half meta-average
#xx <- 0.29 #non-sig. estimates dropped & half meta-average

pps_rstandard$yn20 <- ifelse(pps_rstandard$power <= 0.2,"yes","no")
100*round(prop.table(table(pps_rstandard$yn20)), 4) #no=40.44, yes=59.56

## Summary measures about the power of primary estimates
round(summary(pps_rstandard$power), 2)  #median=13%, mean=32%
a <- round(summary(pps_rstandard$power), 2)[3]
b <- round(summary(pps_rstandard$power), 2)[4]
c <- round(summary(pps_rstandard$power), 2)[2]
d <- round(summary(pps_rstandard$power), 2)[5]

## Summary statistics by subfield - Table 3 (M, N, median, mean, Q25, Q75, sape)
subf_desc <- pps_rstandard %>%
  dplyr::group_by(subfd) %>%
  summarise(M = length(unique(cID)),
            N = length(power),
            median = round(median(power),2),
            mean = round(mean(power),2),
            Q25 = round(quantile(power, probs=0.25),2),
            Q75 = round(quantile(power, probs=0.75),2),
            sape = round(length(which(power >= 0.8))/length(power),2))
print(subf_desc)

## Calculating median of medians by subfield (Table 3)
med_med <- pps_rstandard %>%
  dplyr::group_by(cID) %>%
  summarise(metaID = metaID[1],
            subfd = subfd[1],
            median = median(power))
dim(med_med)
mmedian <- round(summary(med_med$median), 2)[3]

med_med_subf <- med_med %>%
  dplyr::group_by(subfd) %>%
  summarise(mmedian = round(median(median),2))
print(med_med_subf)

## Statistical power in environmental sciences (Table 3)
power.tab3 <- matrix(nrow=8, ncol=8)

# M & N
power.tab3[1,1] <- dim(med_med)[1]; power.tab3[1,2] <- dim(pps_rstandard)[1]
power.tab3[2,1] <- subf_desc$M[1]; power.tab3[2,2] <- subf_desc$N[1]
power.tab3[3,1] <- subf_desc$M[2]; power.tab3[3,2] <- subf_desc$N[2]
power.tab3[4,1] <- subf_desc$M[3]; power.tab3[4,2] <- subf_desc$N[3]
power.tab3[5,1] <- subf_desc$M[4]; power.tab3[5,2] <- subf_desc$N[4]
power.tab3[6,1] <- subf_desc$M[5]; power.tab3[6,2] <- subf_desc$N[5]
power.tab3[7,1] <- subf_desc$M[6]; power.tab3[7,2] <- subf_desc$N[6]
power.tab3[8,1] <- subf_desc$M[7]; power.tab3[8,2] <- subf_desc$N[7]

# Median of medians, raw median and mean
power.tab3[1,3] <- mmedian; power.tab3[1,4] <- a; power.tab3[1,5] <- b
power.tab3[2,3] <- med_med_subf$mmedian[1]; power.tab3[2,4] <- subf_desc$median[1]; power.tab3[2,5] <- subf_desc$mean[1]
power.tab3[3,3] <- med_med_subf$mmedian[2]; power.tab3[3,4] <- subf_desc$median[2]; power.tab3[3,5] <- subf_desc$mean[2]
power.tab3[4,3] <- med_med_subf$mmedian[3]; power.tab3[4,4] <- subf_desc$median[3]; power.tab3[4,5] <- subf_desc$mean[3]
power.tab3[5,3] <- med_med_subf$mmedian[4]; power.tab3[5,4] <- subf_desc$median[4]; power.tab3[5,5] <- subf_desc$mean[4]
power.tab3[6,3] <- med_med_subf$mmedian[5]; power.tab3[6,4] <- subf_desc$median[5]; power.tab3[6,5] <- subf_desc$mean[5]
power.tab3[7,3] <- med_med_subf$mmedian[6]; power.tab3[7,4] <- subf_desc$median[6]; power.tab3[7,5] <- subf_desc$mean[6]
power.tab3[8,3] <- med_med_subf$mmedian[7]; power.tab3[8,4] <- subf_desc$median[7]; power.tab3[8,5] <- subf_desc$mean[7]

# Q25, Q75, & SAPE
power.tab3[1,6] <- c; power.tab3[1,7] <- d; power.tab3[1,8] <- xx
power.tab3[2,6] <- subf_desc$Q25[1]; power.tab3[2,7] <- subf_desc$Q75[1]; power.tab3[2,8] <- subf_desc$sape[1]
power.tab3[3,6] <- subf_desc$Q25[2]; power.tab3[3,7] <- subf_desc$Q75[2]; power.tab3[3,8] <- subf_desc$sape[2]
power.tab3[4,6] <- subf_desc$Q25[3]; power.tab3[4,7] <- subf_desc$Q75[3]; power.tab3[4,8] <- subf_desc$sape[3]
power.tab3[5,6] <- subf_desc$Q25[4]; power.tab3[5,7] <- subf_desc$Q75[4]; power.tab3[5,8] <- subf_desc$sape[4]
power.tab3[6,6] <- subf_desc$Q25[5]; power.tab3[6,7] <- subf_desc$Q75[5]; power.tab3[6,8] <- subf_desc$sape[5]
power.tab3[7,6] <- subf_desc$Q25[6]; power.tab3[7,7] <- subf_desc$Q75[6]; power.tab3[7,8] <- subf_desc$sape[6]
power.tab3[8,6] <- subf_desc$Q25[7]; power.tab3[8,7] <- subf_desc$Q75[7]; power.tab3[8,8] <- subf_desc$sape[7]

colnames(power.tab3) <- c("M","N","mmedian","median","mean","Q25","Q75","SAPE")
rownames(power.tab3) <- c("All meta-analyses","Ecology","Environmental Chemistry","Environmental Engineering","Health, Toxicology and Mutagenesis",
                          "Management, Monitoring, Policy and Law","Nature and Landscape Conservation","Water Science and Technology")
print(power.tab3)
#write.csv(power.tab3,here("results","main","Table_3_legacy_half_meta_average.csv"))
#write.csv(power.tab3,here("results","robustness","Robustness_Table_9_Power_full_meta_average.csv"))
#write.csv(power.tab3,here("results","robustness","Robustness_Table_9_Power_quarter_meta_average.csv"))
#write.csv(power.tab3,here("results","robustness","Robustness_Table_9_Power_null_meta_analyses_dropped.csv"))
write.csv(power.tab3,here("results","robustness","Robustness_Table_9_Power_non_significant_estimates_dropped_half_meta_average.csv"))

## Summary statistics for each meta-analysis (we need for our exploratory regression analyses)
pps_rstandard_median <- pps_rstandard %>%
  dplyr::group_by(cID) %>%
  summarise(metaID=metaID[1],
            median = median(power),
            sape = length(which(power >= 0.8))/length(power),
            nips = length(unique(sID)),
            esty = unique(etype),
            guid = unique(guide),
            prer = unique(prere),
            subf = unique(subfd),
            sdes = unique(sdesn))
dim(pps_rstandard_median)  #704x10

write.xlsx(pps_rstandard_median,here("results","main","Figure_2_data_legacy_half_meta_average.xlsx"),overwrite=T)
#write.xlsx(pps_rstandard_median,here("results","main","Figure_2_data_legacy_full_meta_average.xlsx"),overwrite=T)

length(unique(pps_rstandard_median$metaID)) #202

## Extra summary information at META-level
#>=80% power
pps_rstandard_median$yn80 <- ifelse(pps_rstandard_median$median >=.8, "yes","no")
100*round(prop.table(table(pps_rstandard_median$yn80)), 4) #no=88.21, yes=11.79

#<=20% power
pps_rstandard_median$yn20 <- ifelse(pps_rstandard_median$median <=.2, "yes","no")
100*round(prop.table(table(pps_rstandard_median$yn20)), 4) #no=27.41, yes=72.59

pps_rstandard_median_subf <- pps_rstandard_median %>%
  dplyr::filter(yn80 == "yes")
100*round(prop.table(table(pps_rstandard_median_subf$subf)), 4) #Table S5 in the SI

#Percentage of meta-analyses with 0% share of adequately powered estimates
pps_rstandard_median$yn0 <- ifelse(pps_rstandard_median$sape == 0, "yes","no")
100*round(prop.table(table(pps_rstandard_median$yn0)), 4) #no=40.48, yes=59.52

#Percentage of meta-analyses with >20% share of adequately powered estimates
pps_rstandard_median$sape_yn20 <- ifelse(pps_rstandard_median$sape >.2, "yes","no")
100*round(prop.table(table(pps_rstandard_median$sape_yn20)), 4) #no=81.11, yes=18.89

summary(pps_rstandard_median$median)
pps_rstandard_median$median100 <- round(100*pps_rstandard_median$median,2)
pps_rstandard_median$sape100 <- round(100*pps_rstandard_median$sape,2)


## Figure 2 (median power + share of adequately powered estimates)
med_pwr <- pps_rstandard_median %>%
  ggplot(aes(x=median100, fill=as.factor(yn80))) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins=30, alpha=I(0.6), linewidth=0.1) +
  scale_fill_manual(values=c("brown2", "skyblue2")) +
  xlab("Median statistical power of primary estimates per meta-analysis") +
  ylab("Percentage") +
  ggtitle("(a)") +
  scale_x_continuous(breaks = breaks_width(20),labels=label_percent(scale=1), expand = c(0,0.5)) +
  scale_y_continuous(labels=label_percent(scale=1), expand = c(0,0.5))+
  theme(legend.position="none")+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))
print(med_pwr)

sape <- pps_rstandard_median %>%
  ggplot(aes(x=sape100)) +
  geom_histogram(aes(y = after_stat(count / sum(count) * 100)), bins=30, alpha=I(0.6), linewidth=0.1, fill="skyblue2") +
  xlab("Share of adequately powered primary estimates per meta-analysis") +
  ylab("Percentage") +
  ggtitle("(b)") +
  scale_x_continuous(breaks = breaks_width(20),labels=label_percent(scale=1),expand = c(0, 0.5)) +
  scale_y_continuous(labels=label_percent(scale=1), expand = c(0,0.5))+
  theme(legend.position="none")+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))
print(sape)

pdf(here("results","main","Figure_2_legacy_half_meta_average.pdf"),width=10,height=4)
grid.arrange(med_pwr, sape, ncol = 2)
dev.off()

## Heterogeneity by subfield (for Supplementary Information)
het_subfield <- pps_rstandard %>%
  dplyr::group_by(cID) %>%
  dplyr::summarize(subfd=subfd[1], tau2=tau2[1], isq=isq[1], .groups="drop")

het_distribution <- het_subfield %>%
  dplyr::rename("Subfield"=subfd) %>%
  dplyr::group_by(Subfield) %>%
  dplyr::summarize(M=length(unique(cID)), Median=median(isq), Mean=mean(isq),
                   Q25=quantile(isq,probs=0.25), Q75=quantile(isq,probs=0.75),
                   .groups = "drop") %>%
  dplyr::mutate(dplyr::across(M:Q75, ~ format(round(.x, 2), trim=TRUE)))
print(het_distribution)

## Render with R's native PNG graphics device. Unlike gtsave() for a PNG, this
## does not render an intermediate HTML document and therefore needs no Chrome,
## Chromium, chromote, or other browser automation software.
heterogeneity_table <- gridExtra::tableGrob(
  het_distribution,
  rows=NULL,
  theme=gridExtra::ttheme_minimal(
    base_size=9,
    core=list(fg_params=list(hjust=0.5, x=0.5)),
    colhead=list(fg_params=list(hjust=0.5, x=0.5, fontface="bold"))
  )
)
grDevices::png(
  filename=here("results","robustness","Robustness_Table_1_Heterogeneity_by_subfield.png"),
  width=2400,
  height=700,
  res=200
)
grid::grid.newpage()
grid::grid.draw(heterogeneity_table)
grDevices::dev.off()

## Power analysis (for Supplementary Information)
## NOTE: After subsetting the data, you can go back to Line 207 to calculate power
## for the following two cases.

## Excluding null MAs (Nord et al 2017; Yang et al. 2023)
## Keep these robustness samples in separate objects so the main analysis below
## continues to use the full 704-cluster data set.
pps_rstandard_sig_overall <- pps_rstandard %>% filter(sig_overall < 0.05)
dim(pps_rstandard_sig_overall)

## Excluding statistically non-significant estimates (Lamberink et al. 2018)
pps_rstandard_sig_estimates <- pps_rstandard %>%
  mutate(abs.z = abs(yi/sei)) %>%
  filter(abs.z > 1.96)
dim(pps_rstandard_sig_estimates)

## -----------------------------------------------------------------------------

## Descriptive (using multilevel PET-PEESE + half meta-average - Table 1)
meta <- pps_rstandard
myDat <- list()
jj <- 1
for (i in unique(meta$cID)){
  myDat[[i]] <- meta[which(meta$cID==i), ]
  myDat[[i]]$sape <- length(which(myDat[[i]]$yn80=="yes"))/length(myDat[[i]]$vi)
  jj <- jj + 1
}

## Check
mss <- NULL
for (i in 1:dim(summary(myDat))[1]) {
  mss[i] <- length(myDat[[i]]$sei)
}
length(mss)
sum(mss)

## Experimental vs Observational
e <- 1; o <- 1; m <- 1
i.obs <- NULL; i.exp <- NULL
mydat.exp <- list(); mydat.obs <- list()
for (i in 1:length(summary(myDat)[,1])) {
  if(myDat[[i]]$sdesn[1]=="experimental") {
    mydat.exp[[e]] <- myDat[[i]]
    i.exp[e] <- i
    e <- e+1
  }
  if(myDat[[i]]$sdesn[1]=="observational"|myDat[[i]]$sdesn[1]=="mixed") {
    mydat.obs[[o]] <- myDat[[i]]
    i.obs[o] <- i
    o <- o+1
  }
}
length(mydat.exp) #166
length(mydat.obs) #538

## =0 power vs >0 power
p <- 1; o <- 1
mydat.pow  <- list(); mydat.npow <- list()
i.pow  <- NULL; i.npow <- NULL
for (i in 1:length(summary(myDat)[,1])) {
  if(myDat[[i]]$sape[1] > 0) {
    mydat.pow[[p]] <- myDat[[i]]
    i.pow[p] <- i
    p <- p+1
  }
  if(myDat[[i]]$sape[1]==0) {
    mydat.npow[[o]] <- myDat[[i]]
    i.npow[o] <- i
    o <- o+1
  }
}
length(mydat.pow)  #285
length(mydat.npow) #419

## Guideline (yes vs no)
g <- 1; h <- 1
mydat.gui  <- list(); mydat.ngu <- list()
i.gui  <- NULL; i.ngu <- NULL
for (i in 1:length(summary(myDat)[,1])) {
  if(myDat[[i]]$guide[1]== "yes") {
    mydat.gui[[g]] <- myDat[[i]]
    i.gui[g] <- i
    g <- g+1
  }
  if(myDat[[i]]$guide[1] == "no") {
    mydat.ngu[[h]] <- myDat[[i]]
    i.ngu[h] <- i
    h <- h+1
  }
}
length(mydat.gui) #363
length(mydat.ngu) #341

## Registration (yes vs no)
r <- 1; s <- 1
mydat.reg <- list(); mydat.nre <- list()
i.reg  <- NULL; i.nre <- NULL
for (i in 1:length(summary(myDat)[,1])) {
  if(myDat[[i]]$prere[1]=="yes") {
    mydat.reg[[r]] <- myDat[[i]]
    i.reg[r] <- i
    r <- r+1
  }
  if(myDat[[i]]$prere[1]=="no") {
    mydat.nre[[s]] <- myDat[[i]]
    i.nre[s] <- i
    s <- s+1
  }
}
length(mydat.reg)  #43
length(mydat.nre)  #661

## Table 1
d.tab <- matrix(nrow=9, ncol=9)

## Full
d.tab[1,1] <- length(myDat)
d.tab[1,2] <- sum(mss)
d.tab[1,3] <- round(mean(mss))
d.tab[1,4] <- round(median(mss))
d.tab[1,5] <- min(mss)
d.tab[1,6] <- round(quantile(mss, probs=0.25))
d.tab[1,7] <- round(median(mss))
d.tab[1,8] <- round(quantile(mss, probs=0.75))
d.tab[1,9] <- max(mss)

## Observational & experimental
d.tab[2,1] <- length(mydat.obs);                       d.tab[3,1] <- length(mydat.exp)
d.tab[2,2] <- sum(mss[i.obs]);                         d.tab[3,2] <- sum(mss[i.exp])
d.tab[2,3] <- round(mean(mss[i.obs]));                 d.tab[3,3] <- round(mean(mss[i.exp]))
d.tab[2,4] <- round(median(mss[i.obs]));               d.tab[3,4] <- round(median(mss[i.exp]))
d.tab[2,5] <- min(mss[i.obs]);                         d.tab[3,5] <- min(mss[i.exp])
d.tab[2,6] <- round(quantile(mss[i.obs], probs=0.25)); d.tab[3,6] <- round(quantile(mss[i.exp], probs=0.25))
d.tab[2,7] <- round(median(mss[i.obs]));               d.tab[3,7] <- round(median(mss[i.exp]))
d.tab[2,8] <- round(quantile(mss[i.obs], probs=0.75)); d.tab[3,8] <- round(quantile(mss[i.exp], probs=0.75))
d.tab[2,9] <- max(mss[i.obs]);                         d.tab[3,9] <- max(mss[i.exp])

## >0 power & =0 power
d.tab[4,1] <- length(mydat.pow);                       d.tab[5,1] <- length(mydat.npow)
d.tab[4,2] <- sum(mss[i.pow]);                         d.tab[5,2] <- sum(mss[i.npow])
d.tab[4,3] <- round(mean(mss[i.pow]));                 d.tab[5,3] <- round(mean(mss[i.npow]))
d.tab[4,4] <- round(median(mss[i.pow]));               d.tab[5,4] <- round(median(mss[i.npow]))
d.tab[4,5] <- min(mss[i.pow]);                         d.tab[5,5] <- min(mss[i.npow])
d.tab[4,6] <- round(quantile(mss[i.pow], probs=0.25)); d.tab[5,6] <- round(quantile(mss[i.npow], probs=0.25))
d.tab[4,7] <- round(median(mss[i.pow]));               d.tab[5,7] <- round(median(mss[i.npow]))
d.tab[4,8] <- round(quantile(mss[i.pow], probs=0.75)); d.tab[5,8] <- round(quantile(mss[i.npow], probs=0.75))
d.tab[4,9] <- max(mss[i.pow]);                         d.tab[5,9] <- max(mss[i.npow])

## Guideline (yes vs no)
d.tab[6,1] <- length(mydat.gui);                       d.tab[7,1] <- length(mydat.ngu)
d.tab[6,2] <- sum(mss[i.gui]);                         d.tab[7,2] <- sum(mss[i.ngu])
d.tab[6,3] <- round(mean(mss[i.gui]));                 d.tab[7,3] <- round(mean(mss[i.ngu]))
d.tab[6,4] <- round(median(mss[i.gui]));               d.tab[7,4] <- round(median(mss[i.ngu]))
d.tab[6,5] <- min(mss[i.gui]);                         d.tab[7,5] <- min(mss[i.ngu])
d.tab[6,6] <- round(quantile(mss[i.gui], probs=0.25)); d.tab[7,6] <- round(quantile(mss[i.ngu], probs=0.25))
d.tab[6,7] <- round(median(mss[i.gui]));               d.tab[7,7] <- round(median(mss[i.ngu]))
d.tab[6,8] <- round(quantile(mss[i.gui], probs=0.75)); d.tab[7,8] <- round(quantile(mss[i.ngu], probs=0.75))
d.tab[6,9] <- max(mss[i.gui]);                         d.tab[7,9] <- max(mss[i.ngu])

## Registration (yes vs no)
d.tab[8,1] <- length(mydat.reg);                       d.tab[9,1] <- length(mydat.nre)
d.tab[8,2] <- sum(mss[i.reg]);                         d.tab[9,2] <- sum(mss[i.nre])
d.tab[8,3] <- round(mean(mss[i.reg]));                 d.tab[9,3] <- round(mean(mss[i.nre]))
d.tab[8,4] <- round(median(mss[i.reg]));               d.tab[9,4] <- round(median(mss[i.nre]))
d.tab[8,5] <- min(mss[i.reg]);                         d.tab[9,5] <- min(mss[i.nre])
d.tab[8,6] <- round(quantile(mss[i.reg], probs=0.25)); d.tab[9,6] <- round(quantile(mss[i.nre], probs=0.25))
d.tab[8,7] <- round(median(mss[i.reg]));               d.tab[9,7] <- round(median(mss[i.nre]))
d.tab[8,8] <- round(quantile(mss[i.reg], probs=0.75)); d.tab[9,8] <- round(quantile(mss[i.nre], probs=0.75))
d.tab[8,9] <- max(mss[i.reg]);                         d.tab[9,9] <- max(mss[i.nre])

rownames(d.tab) <- c("All meta-analyses", "Observational","Experimental","SAPE > 0", "SAPE = 0","Yes","No","Yes","No")
colnames(d.tab) <- c("M", "N", "Mean","Median", "Min", "Q25", "Q50", "Q75","Max")
print(d.tab)
write.csv(d.tab,here("results","main","Table_1_legacy_half_meta_average.csv"))


## -----------------------------------------------------------------------------

## Data preparation for subfields
my_dat <- pps_rstandard
my_dat$sdesn_merge <- ifelse(my_dat$sdesn %in% c("observational", "mixed"), "observational","experimental")

myDat <- list()
jj <- 1
for (i in unique(my_dat$cID)){
  myDat[[i]] <- my_dat[which(my_dat$cID==i), ]
  myDat[[i]]$sape <- length(which(myDat[[i]]$yn80=="yes"))/length(myDat[[i]]$vi)
  jj <- jj + 1
}

## Check
mss <- NULL
for (i in 1:dim(summary(myDat))[1]) {
  mss[i] <- length(myDat[[i]]$vi)
}
length(mss) #704
sum(mss)    #63956

p <- 1; o <- 1; q <- 1; r <- 1; s <- 1; t <- 1; v <- 1
mydat.eco <- list(); mydat.enc <- list()
mydat.ene <- list(); mydat.nlc <- list()
mydat.mpl <- list(); mydat.wst <- list()
mydat.htm <- list()

i.eco <- NULL; i.enc <- NULL; i.ene <- NULL
i.mpl <- NULL; i.nlc <- NULL; i.wst <- NULL
i.htm <- NULL

for (i in 1:length(summary(myDat)[,1])) {
  if(myDat[[i]]$subfd[1] == "Ecology") {
    mydat.eco[[p]] <- myDat[[i]]
    i.eco[p] <- i
    p <- p+1
  }
  if(myDat[[i]]$subfd[1] == "Environmental Chemistry") {
    mydat.enc[[o]] <- myDat[[i]]
    i.enc[o] <- i
    o <- o+1
  }
  if(myDat[[i]]$subfd[1] == "Environmental Engineering") {
    mydat.ene[[q]] <- myDat[[i]]
    i.ene[q] <- i
    q <- q+1
  }
  if(myDat[[i]]$subfd[1] == "Nature and Landscape Conservation") {
    mydat.nlc[[r]] <- myDat[[i]]
    i.nlc[r] <- i
    r <- r+1
  }
  if(myDat[[i]]$subfd[1] == "Management, Monitoring, Policy and Law") {
    mydat.mpl[[s]] <- myDat[[i]]
    i.mpl[s] <- i
    s <- s+1
  }
  if(myDat[[i]]$subfd[1] == "Water Science and Technology") {
    mydat.wst[[t]] <- myDat[[i]]
    i.wst[t] <- i
    t <- t+1
  }
  if(myDat[[i]]$subfd[1] == "Health, Toxicology and Mutagenesis") {
    mydat.htm[[v]] <- myDat[[i]]
    i.htm[v] <- i
    v <- v+1
  }
}


## -------------------------
## z-value and p-value grids
## -------------------------

# z.grid is used for visualization and goes from 1 to 10 while
# p.grid is used for exact calculations and goes from z-values of 0 to Inf (absolute values).

## p.grid for main results table
p.grid.tab <- c(-Inf,
                qnorm(0.001/2, lower.tail = TRUE),
                qnorm(0.01/2,  lower.tail = TRUE),
                qnorm(0.05/2,  lower.tail = TRUE),
                qnorm(0.1/2,   lower.tail = TRUE),
                qnorm(0.2/2,   lower.tail = TRUE),
                qnorm(0.3/2,   lower.tail = TRUE),
                qnorm(0.4/2,   lower.tail = TRUE),
                qnorm(0.5/2,   lower.tail = TRUE),
                qnorm(0.6/2,   lower.tail = TRUE),
                qnorm(0.7/2,   lower.tail = TRUE),
                qnorm(0.8/2,   lower.tail = TRUE),
                qnorm(0.9/2,   lower.tail = TRUE),
                0,
                qnorm(0.9/2,   lower.tail = FALSE),
                qnorm(0.8/2,   lower.tail = FALSE),
                qnorm(0.7/2,   lower.tail = FALSE),
                qnorm(0.6/2,   lower.tail = FALSE),
                qnorm(0.5/2,   lower.tail = FALSE),
                qnorm(0.4/2,   lower.tail = FALSE),
                qnorm(0.3/2,   lower.tail = FALSE),
                qnorm(0.2/2,   lower.tail = FALSE),
                qnorm(0.1/2,   lower.tail = FALSE),
                qnorm(0.05/2,  lower.tail = FALSE),
                qnorm(0.01/2,  lower.tail = FALSE),
                qnorm(0.001/2, lower.tail = FALSE),
                Inf)

# We need to calculate for the original z-values also the frequencies
# and we can do this directly for the two-sided test grid
p.grid.tab2 <- c(0,
                 qnorm(0.9/2,   lower.tail = FALSE),
                 qnorm(0.8/2,   lower.tail = FALSE),
                 qnorm(0.7/2,   lower.tail = FALSE),
                 qnorm(0.6/2,   lower.tail = FALSE),
                 qnorm(0.5/2,   lower.tail = FALSE),
                 qnorm(0.4/2,   lower.tail = FALSE),
                 qnorm(0.3/2,   lower.tail = FALSE),
                 qnorm(0.2/2,   lower.tail = FALSE),
                 qnorm(0.1/2,   lower.tail = FALSE),
                 qnorm(0.05/2,  lower.tail = FALSE),
                 qnorm(0.01/2,  lower.tail = FALSE),
                 qnorm(0.001/2, lower.tail = FALSE),
                 Inf)

## p grid for plot -> for exact calculation
p.grid.plot  <- qnorm(seq(1, 0, -0.005), lower.tail = FALSE)
p.grid.plot2 <- p.grid.plot[which(p.grid.plot >= 0)] # this is for actual plotting of two-sided tests

## z.grid for plot
z.grid.plot <- seq(-10.25, 10.25, 0.1025)
z.grid.plot2 <- seq(0, 10.25, 0.1025) # this is for actual plotting of absolute values

## ----------------------------------------------------
## Prepare original data for plotting and p-value table
## ----------------------------------------------------

my_dat <- pps_rstandard
my_dat$GE <- my_dat$GE/2

myDat <- list()
jj <- 1
for (i in unique(my_dat$cID)){
  myDat[[i]] <- my_dat[which(my_dat$cID==i), ]
  jj <- jj + 1
}

## Check
mss <- NULL
for (i in 1:dim(summary(myDat))[1]) {
  mss[i] <- length(myDat[[i]]$vi)
}
length(mss) #704
sum(mss)    #63956

## ---

my_dat <- pps_rstandard
my_dat$sdesn_merge <- ifelse(my_dat$sdesn %in% c("observational", "mixed"), "observational","experimental")

myDat <- list()
jj <- 1
for (i in unique(my_dat$cID)){
  myDat[[i]] <- my_dat[which(my_dat$cID==i), ]
  myDat[[i]]$sape <- length(which(myDat[[i]]$yn80=="yes"))/length(myDat[[i]]$vi)
  jj <- jj + 1
}

## Check
mss <- NULL
for (i in 1:dim(summary(myDat))[1]) {
  mss[i] <- length(myDat[[i]]$vi)
}
length(mss) #704
sum(mss)    #63956

## ---

p <- 1; o <- 1; q <- 1; r <- 1; s <- 1; t <- 1; v <- 1
mydat.eco <- list(); mydat.enc <- list()
mydat.ene <- list(); mydat.nlc <- list()
mydat.mpl <- list(); mydat.wst <- list()
mydat.htm <- list()

i.eco <- NULL; i.enc <- NULL; i.ene <- NULL
i.mpl <- NULL; i.nlc <- NULL; i.wst <- NULL
i.htm <- NULL

for (i in 1:length(summary(myDat)[,1])) {
  if(myDat[[i]]$subfd[1] == "Ecology") {
    mydat.eco[[p]] <- myDat[[i]]
    i.eco[p] <- i
    p <- p+1
  }
  if(myDat[[i]]$subfd[1] == "Environmental Chemistry") {
    mydat.enc[[o]] <- myDat[[i]]
    i.enc[o] <- i
    o <- o+1
  }
  if(myDat[[i]]$subfd[1] == "Environmental Engineering") {
    mydat.ene[[q]] <- myDat[[i]]
    i.ene[q] <- i
    q <- q+1
  }
  if(myDat[[i]]$subfd[1] == "Nature and Landscape Conservation") {
    mydat.nlc[[r]] <- myDat[[i]]
    i.nlc[r] <- i
    r <- r+1
  }
  if(myDat[[i]]$subfd[1] == "Management, Monitoring, Policy and Law") {
    mydat.mpl[[s]] <- myDat[[i]]
    i.mpl[s] <- i
    s <- s+1
  }
  if(myDat[[i]]$subfd[1] == "Water Science and Technology") {
    mydat.wst[[t]] <- myDat[[i]]
    i.wst[t] <- i
    t <- t+1
  }
  if(myDat[[i]]$subfd[1] == "Health, Toxicology and Mutagenesis") {
    mydat.htm[[v]] <- myDat[[i]]
    i.htm[v] <- i
    v <- v+1
  }
}

## ---

fac <- do.call(rbind.data.frame, myDat)
facz <- abs(fac$yi / sqrt(fac$vi))          #factual z
#facz <- abs(fac$yi / sqrt(fac$vi_orig))    #factual z (for 1.5xse)

fac <- do.call(rbind.data.frame, mydat.eco)
facz.eco <- abs(fac$yi / fac$sei)

fac <- do.call(rbind.data.frame, mydat.enc)
facz.enc <- abs(fac$yi / fac$sei)

fac <- do.call(rbind.data.frame, mydat.ene)
facz.ene <- abs(fac$yi / fac$sei)

fac <- do.call(rbind.data.frame, mydat.nlc)
facz.nlc <- abs(fac$yi / fac$sei)

fac <- do.call(rbind.data.frame, mydat.mpl)
facz.mpl <- abs(fac$yi / fac$sei)

fac <- do.call(rbind.data.frame, mydat.wst)
facz.wst <- abs(fac$yi / fac$sei)

fac <- do.call(rbind.data.frame, mydat.htm)
facz.htm <- abs(fac$yi / fac$sei)

## -----------------------------------------------------------------------------

## Full
z.orig <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig[a] <- length(which(facz >= z.grid.plot2[a] & facz <= z.grid.plot2[a+1]))
}

p.orig.tab <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab[a] <- length(which(facz >= p.grid.tab2[a] & facz <= p.grid.tab2[a+1]))
}

p.orig.plot <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot[a] <- length(which(facz >= p.grid.plot2[a] & facz <= p.grid.plot2[a+1]))
}

## Ecology
z.orig.eco <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.eco[a] <- length(which(facz.eco >= z.grid.plot2[a] & facz.eco <= z.grid.plot2[a+1]))
}

p.orig.tab.eco <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.eco[a] <- length(which(facz.eco >= p.grid.tab2[a] & facz.eco <= p.grid.tab2[a+1]))
}

p.orig.plot.eco <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.eco[a] <- length(which(facz.eco >= p.grid.plot2[a] & facz.eco <= p.grid.plot2[a+1]))
}

## Environmental Chemistry
z.orig.enc <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.enc[a] <- length(which(facz.enc >= z.grid.plot2[a] & facz.enc <= z.grid.plot2[a+1]))
}

p.orig.tab.enc <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.enc[a] <- length(which(facz.enc >= p.grid.tab2[a] & facz.enc <= p.grid.tab2[a+1]))
}

p.orig.plot.enc <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.enc[a] <- length(which(facz.enc >= p.grid.plot2[a] & facz.enc <= p.grid.plot2[a+1]))
}

#Environmental Engineering
z.orig.ene <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.ene[a] <- length(which(facz.ene >= z.grid.plot2[a] & facz.ene <= z.grid.plot2[a+1]))
}

p.orig.tab.ene <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.ene[a] <- length(which(facz.ene >= p.grid.tab2[a] & facz.ene <= p.grid.tab2[a+1]))
}

p.orig.plot.ene <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.ene[a] <- length(which(facz.ene >= p.grid.plot2[a] & facz.ene <= p.grid.plot2[a+1]))
}

## Management
z.orig.mpl <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.mpl[a] <- length(which(facz.mpl >= z.grid.plot2[a] & facz.mpl <= z.grid.plot2[a+1]))
}

p.orig.tab.mpl <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.mpl[a] <- length(which(facz.mpl >= p.grid.tab2[a] & facz.mpl <= p.grid.tab2[a+1]))
}

p.orig.plot.mpl <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.mpl[a] <- length(which(facz.mpl >= p.grid.plot2[a] & facz.mpl <= p.grid.plot2[a+1]))
}

## Nature and Landscape Conservation
z.orig.nlc <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.nlc[a] <- length(which(facz.nlc >= z.grid.plot2[a] & facz.nlc <= z.grid.plot2[a+1]))
}

p.orig.tab.nlc <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.nlc[a] <- length(which(facz.nlc >= p.grid.tab2[a] & facz.nlc <= p.grid.tab2[a+1]))
}

p.orig.plot.nlc <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.nlc[a] <- length(which(facz.nlc >= p.grid.plot2[a] & facz.nlc <= p.grid.plot2[a+1]))
}

## Water Science and Technology
z.orig.wst <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.wst[a] <- length(which(facz.wst >= z.grid.plot2[a] & facz.wst <= z.grid.plot2[a+1]))
}

p.orig.tab.wst <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.wst[a] <- length(which(facz.wst >= p.grid.tab2[a] & facz.wst <= p.grid.tab2[a+1]))
}

p.orig.plot.wst <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.wst[a] <- length(which(facz.wst >= p.grid.plot2[a] & facz.wst <= p.grid.plot2[a+1]))
}

## Health, Toxicology and Mutagenesis
z.orig.htm <- NULL
for (a in 1:(length(z.grid.plot2)-1)) {
  #point probability is zero, but I added >= and <= to ensure that 0 and Inf are included
  z.orig.htm[a] <- length(which(facz.htm >= z.grid.plot2[a] & facz.htm <= z.grid.plot2[a+1]))
}

p.orig.tab.htm <- NULL
for (a in 1:(length(p.grid.tab2)-1)) {
  p.orig.tab.htm[a] <- length(which(facz.htm >= p.grid.tab2[a] & facz.htm <= p.grid.tab2[a+1]))
}

p.orig.plot.htm <- NULL
for (a in 1:(length(p.grid.plot2)-1)) {
  p.orig.plot.htm[a] <- length(which(facz.htm >= p.grid.plot2[a] & facz.htm <= p.grid.plot2[a+1]))
}


## ---

## Simulating the counterfactual z-values

cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot <- cf(dat = myDat, z.grid = z.grid.plot, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
saveRDS(z.plot, file=here("results","main","z_plot_pet_peese_rstandard_half.rds"))


## Confidence interval for the counterfactuals
it <- n_iterations                  #bootstrap size
clu <- unique(pps_rstandard$cID)    #cluster id
length(clu)                         #704

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ci <- cf.ci.cluster(dat = myDat, z.grid = z.grid.plot, iters = it, cluster = clu, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) #about 12 hrs with 7 cores
saveRDS(z.plot.ci, file=here("results","main","z_plot_ci_pet_peese_rstandard_half.rds"))


## Load the data
z.plot <- readRDS(here("results","main","z_plot_pet_peese_rstandard_half.rds"))
z.plot.ci <- readRDS(here("results","main","z_plot_ci_pet_peese_rstandard_half.rds"))


## -----------------
## z-plot (Figure 1)
## -----------------

## z-plot Figure
# z.grid
xs <- as.vector(z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2)
N <- sum(p.orig.plot) # Here z.orig goes until 10 and is smaller than the full sample size!
# But the bootstrapping of the confidence interval is based on the entire sample
# and given in shares. So we need to divide here by the full sample to be consistent.

q025 <- as.vector(apply(z.plot.ci[[1]], 2, quantile, na.rm=T, probs=c(0.025)))
q975 <- as.vector(apply(z.plot.ci[[1]], 2, quantile, na.rm=T, probs=c(0.975)))
n.f <- as.vector(z.orig/N)
n.cf <- as.vector(z.plot/N)
datFull <- as.data.frame(cbind(xs, q025, q975, n.f, n.cf))
dim(datFull)

#pdf(here("results","robustness","Robustness_Figure_1_WAAP_z20.pdf"),width=10,height=5)
pdf(here("results","main","Figure_1_legacy_half_meta_average.pdf"),width=10,height=5)
ggplot(datFull) +
  geom_line(aes(xs, q025), color='orange', lty=3) +
  geom_line(aes(xs, n.cf), color='orange', lty=1) +
  geom_point(aes(xs,n.cf), shape=20, fill='orange', color='orange', size=1)+
  geom_line(aes(xs, q975), color='orange', lty=3) +
  geom_line(aes(xs, n.f), color='blue', lty=2) +
  geom_point(aes(xs,n.f), shape=20, fill='blue', color='blue', size=1)+
  coord_cartesian(xlim=c(0,8)) +
  xlab('|z|-value') + ylab('Frequency') +
  geom_vline(xintercept = c(1.64,1.96,2.58), lty=2, color=c(3,2,6), linewidth=0.5)+
  scale_x_continuous(breaks = c(0,1.64,1.96,2.58,4,6,8)) +
  theme(panel.background = element_rect(fill = "gray100"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "gray"))
dev.off()


## -------------------------------
## Table 2 (and tables in the SI)
## -------------------------------

## Functions
source(here("scripts", "functions.R"))
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab <- cf(dat = myDat, z.grid = p.grid.tab, heterogeneity_multiplier = heterogeneity_multiplier) #VERY IMPORTANT: Make sure sure that you're using the appropriate function for each case.
stopCluster(cl)

saveRDS(p.tab, file=here("results","main","p_tab_pps_rstandard_half.rds"))
#saveRDS(p.tab, file=here("results","main","p_tab_pps_rstandard_half_0.25xtau2.rds"))
#saveRDS(p.tab, file=here("results","main","p_tab_pps_rstandard_half_0.5xtau2.rds"))

clu <- unique(pps_rstandard$cID)
length(clu) #704
it <- n_iterations

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ci <- cf.ci.cluster(dat = myDat, z.grid = p.grid.tab, iters = it, cluster = clu, heterogeneity_multiplier = heterogeneity_multiplier) #VERY IMPORTANT: Make sure sure that you're using the appropriate function for each case.
#p.tab.ci <- cf.ci.cluster.se(dat=myDat, z.grid=p.grid.tab, iters=it, cluster=clu) #for 1.5*se
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)

saveRDS(p.tab.ci, file=here("results","main","p_tab_ci_pps_rstandard_half.rds"))
#saveRDS(p.tab.ci, file=here("results","main","p_tab_ci_pps_rstandard_half_0.25xtau2.rds"))
#saveRDS(p.tab.ci, file=here("results","main","p_tab_ci_pps_rstandard_half_0.5xtau2.rds"))

# Load the data
p.tab <- readRDS(here("results","main","p_tab_pps_rstandard_half.rds"))
p.tab.ci <- readRDS(here("results","main","p_tab_ci_pps_rstandard_half.rds"))
#p.tab <- readRDS(here("results","main","p_tab_pps_rstandard_half_0.25xtau2.rds"))
#p.tab.ci <- readRDS(here("results","main","p_tab_ci_pps_rstandard_half_0.25xtau2.rds"))
#p.tab <- readRDS(here("results","main","p_tab_pps_rstandard_half_0.5xtau2.rds"))
#p.tab.ci <- readRDS(here("results","main","p_tab_ci_pps_rstandard_half_0.5xtau2.rds"))

# ---

# Full
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 CI")
N <- sum(p.orig.tab)  #BE CAREFUL HERE! Check whether it's equal to the total!
p.table[,1] <- round((p.orig.tab - p.tab) / N,  3)
q025 <- apply(matrix(p.orig.tab/N, nrow=nrow(p.tab.ci[[1]]), ncol=ncol(p.tab.ci[[1]]), byrow=TRUE)
              - p.tab.ci[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab/N, nrow=nrow(p.tab.ci[[1]]), ncol=ncol(p.tab.ci[[1]]), byrow=TRUE)
              - p.tab.ci[[1]], 2, quantile, probs=c(0.975))
p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab - p.tab)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci[[1]][, 10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[10:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab - p.tab)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci[[1]][, 11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#inflated significance (10%)
shp.10 <- round(sum((p.orig.tab - p.tab)[10:13] / sum(p.orig.tab[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#inflated significance (5%)
shp.5 <- round(sum((p.orig.tab - p.tab)[11:13] / sum(p.orig.tab[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#sample size
p.table <- rbind(p.table, c(length(summary(myDat)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
                       "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
                       "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01",
                       "p < 0.001","ESR_{0.1}^{all}","ESR_{0.05}^{all}","ESR_{0.1}^{sig}","ESR_{0.05}^{sig}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)

write.csv(p.table, here("results","main","Table_2_legacy_half_meta_average.csv"))
#write.csv(p.table, here("results","robustness","Robustness_Table_10_ESR_quarter_heterogeneity.csv"))
#write.csv(p.table, here("results","robustness","Robustness_Table_10_ESR_half_heterogeneity.csv"))

## -------------
## For subfields
## -------------

## Ecology
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.eco <- cf(dat = mydat.eco, z.grid = z.grid.plot, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(z.plot.eco, file=here("results","robustness","pet_peese_rstandard_z_plot.eco.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ci.eco <- cf.ci.cluster(dat = mydat.eco, z.grid = z.grid.plot, iters = it, cluster = clu[i.eco], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time) #about 33 mins
save(z.plot.ci.eco,file=here("results","robustness","pet_peese_rstandard_z_plot_ci.eco.rds"))

## Environmental Chemistry
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.enc <- cf(dat = mydat.enc, z.grid = z.grid.plot, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(z.plot.enc, file=here("results","robustness","pet_peese_rstandard_z_plot.enc.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ci.enc <- cf.ci.cluster(dat = mydat.enc, z.grid = z.grid.plot, iters = it, cluster = clu[i.enc], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)
save(z.plot.ci.enc,file=here("results","robustness","pet_peese_rstandard_z_plot_ci.enc.rds"))

## Environmental Engineering
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ene <- cf(dat = mydat.ene, z.grid = z.grid.plot, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(z.plot.ene, file=here("results","robustness","pet_peese_rstandard_z_plot.ene.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ci.ene <- cf.ci.cluster(dat = mydat.ene, z.grid = z.grid.plot, iters = it, cluster = clu[i.ene], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)
save(z.plot.ci.ene,file=here("results","robustness","pet_peese_rstandard_z_plot_ci.ene.rds"))

## Nature and Landscape Conservation
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.nlc <- cf(dat = mydat.nlc, z.grid = z.grid.plot, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(z.plot.nlc, file=here("results","robustness","pet_peese_rstandard_z_plot.nlc.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ci.nlc <- cf.ci.cluster(dat = mydat.nlc, z.grid = z.grid.plot, iters = it, cluster = clu[i.nlc], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)
save(z.plot.ci.nlc,file=here("results","robustness","pet_peese_rstandard_z_plot_ci.nlc.rds"))

## Management
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.mpl <- cf(dat = mydat.mpl, z.grid = z.grid.plot, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(z.plot.mpl, file=here("results","robustness","pet_peese_rstandard_z_plot.mpl.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ci.mpl <- cf.ci.cluster(dat = mydat.mpl, z.grid = z.grid.plot, iters = it, cluster = clu[i.mpl], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)
save(z.plot.ci.mpl,file=here("results","robustness","pet_peese_rstandard_z_plot_ci.mpl.rds"))

## Water Science and Technology
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.wst <- cf(dat = mydat.wst, z.grid = z.grid.plot, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(z.plot.wst, file=here("results","robustness","pet_peese_rstandard_z_plot.wst.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ci.wst <- cf.ci.cluster(dat = mydat.wst, z.grid = z.grid.plot, iters = it, cluster = clu[i.wst], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)
save(z.plot.ci.wst,file=here("results","robustness","pet_peese_rstandard_z_plot_ci.wst.rds"))

## Health, Toxicology and Mutagenesis
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.htm <- cf(dat = mydat.htm, z.grid = z.grid.plot, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(z.plot.htm, file=here("results","robustness","pet_peese_rstandard_z_plot.htm.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
z.plot.ci.htm <- cf.ci.cluster(dat = mydat.htm, z.grid = z.grid.plot, iters = it, cluster = clu[i.htm], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)
save(z.plot.ci.htm,file=here("results","robustness","pet_peese_rstandard_z_plot_ci.htm.rds"))


## Subfields
load(here("results","robustness","pet_peese_rstandard_z_plot.eco.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot_ci.eco.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot.enc.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot_ci.enc.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot.ene.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot_ci.ene.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot.nlc.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot_ci.nlc.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot.wst.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot_ci.wst.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot.mpl.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot_ci.mpl.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot.htm.rds"))
load(here("results","robustness","pet_peese_rstandard_z_plot_ci.htm.rds"))


## z-plot for subfields
# eco
xs.eco <- as.vector(z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2)
N.eco <- sum(p.orig.plot.eco)
q025.eco <- as.vector(apply(z.plot.ci.eco[[1]], 2, quantile, probs=c(0.025)))
q975.eco <- as.vector(apply(z.plot.ci.eco[[1]], 2, quantile, probs=c(0.975)))
n.f.eco <- as.vector(z.orig.eco/N.eco)
n.cf.eco <- as.vector(z.plot.eco/N.eco)

# enc
xs.enc <- as.vector(z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2)
N.enc <- sum(p.orig.plot.enc)
q025.enc <- as.vector(apply(z.plot.ci.enc[[1]], 2, quantile, probs=c(0.025)))
q975.enc <- as.vector(apply(z.plot.ci.enc[[1]], 2, quantile, probs=c(0.975)))
n.f.enc <- as.vector(z.orig.enc/N.enc)
n.cf.enc <- as.vector(z.plot.enc/N.enc)

# ene
xs.ene <- as.vector(z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2)
N.ene <- sum(p.orig.plot.ene)
q025.ene <- as.vector(apply(z.plot.ci.ene[[1]], 2, quantile, probs=c(0.025)))
q975.ene <- as.vector(apply(z.plot.ci.ene[[1]], 2, quantile, probs=c(0.975)))
n.f.ene <- as.vector(z.orig.ene/N.ene)
n.cf.ene <- as.vector(z.plot.ene/N.ene)

# nlc
xs.nlc <- as.vector(z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2)
N.nlc <- sum(p.orig.plot.nlc)
q025.nlc <- as.vector(apply(z.plot.ci.nlc[[1]], 2, quantile, probs=c(0.025)))
q975.nlc <- as.vector(apply(z.plot.ci.nlc[[1]], 2, quantile, probs=c(0.975)))
n.f.nlc <- as.vector(z.orig.nlc/N.nlc)
n.cf.nlc <- as.vector(z.plot.nlc/N.nlc)

# Management
xs.mpl <- as.vector(z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2)
N.mpl <- sum(p.orig.plot.mpl)
q025.mpl <- as.vector(apply(z.plot.ci.mpl[[1]], 2, quantile, probs=c(0.025)))
q975.mpl <- as.vector(apply(z.plot.ci.mpl[[1]], 2, quantile, probs=c(0.975)))
n.f.mpl <- as.vector(z.orig.mpl/N.mpl)
n.cf.mpl <- as.vector(z.plot.mpl/N.mpl)

# wst
xs.wst <- as.vector(z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2)
N.wst <- sum(p.orig.plot.wst)
q025.wst <- as.vector(apply(z.plot.ci.wst[[1]], 2, quantile, probs=c(0.025)))
q975.wst <- as.vector(apply(z.plot.ci.wst[[1]], 2, quantile, probs=c(0.975)))
n.f.wst <- as.vector(z.orig.wst/N.wst)
n.cf.wst <- as.vector(z.plot.wst/N.wst)

# htm
xs.htm <- as.vector(z.grid.plot2[-length(z.grid.plot2)] + (z.grid.plot2[2]-z.grid.plot2[1])/2)
N.htm <- sum(p.orig.plot.htm)
q025.htm <- as.vector(apply(z.plot.ci.htm[[1]], 2, quantile, probs=c(0.025)))
q975.htm <- as.vector(apply(z.plot.ci.htm[[1]], 2, quantile, probs=c(0.975)))
n.f.htm <- as.vector(z.orig.htm/N.htm)
n.cf.htm <- as.vector(z.plot.htm/N.htm)

datSUBF <- as.data.frame(cbind(xs.eco,q025.eco,q975.eco,n.f.eco,n.cf.eco,
                               xs.enc,q025.enc,q975.enc,n.f.enc,n.cf.enc,
                               xs.ene,q025.ene,q975.ene,n.f.ene,n.cf.ene,
                               xs.nlc,q025.nlc,q975.nlc,n.f.nlc,n.cf.nlc,
                               xs.mpl,q025.mpl,q975.mpl,n.f.mpl,n.cf.mpl,
                               xs.wst,q025.wst,q975.wst,n.f.wst,n.cf.wst,
                               xs.htm,q025.htm,q975.htm,n.f.htm,n.cf.htm))
dim(datSUBF)

eco <- ggplot(datSUBF) +
  geom_line(aes(xs.eco, q025.eco), color='orange', lty=3) +
  geom_line(aes(xs.eco, n.cf.eco), color='orange', lty=1) +
  geom_point(aes(xs.eco, n.cf.eco),shape=20, fill='orange', color='orange', size=1)+
  geom_line(aes(xs.eco, q975.eco), color='orange', lty=3) +
  geom_line(aes(xs.eco, n.f.eco), color='blue', lty=2) +
  geom_point(aes(xs.eco, n.f.eco),shape=20, fill='blue', color='blue', size=1)+
  coord_cartesian(xlim=c(0,8)) + xlab('|z|-value') + ylab('Frequency') +
  ggtitle('Ecology') + easy_center_title()+
  geom_vline(xintercept = c(1.64,1.96,2.58), lty=2, color=c(3,2,6), linewidth=0.5)+
  scale_x_continuous(breaks=c(0,1.64,1.96,2.58,4,6,8)) +
  theme(panel.background = element_rect(fill = "gray100"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

enc <- ggplot(datSUBF) +
  geom_line(aes(xs.enc, q025.enc), color='orange', lty=3) +
  geom_line(aes(xs.enc, n.cf.enc), color='orange', lty=1) +
  geom_point(aes(xs.enc, n.cf.enc),shape=20, fill='orange', color='orange', size=1)+
  geom_line(aes(xs.enc, q975.enc), color='orange', lty=3) +
  geom_line(aes(xs.enc, n.f.enc), color='blue', lty=2) +
  geom_point(aes(xs.enc, n.f.enc),shape=20, fill='blue', color='blue', size=1)+
  coord_cartesian(xlim=c(0,8)) + xlab('|z|-value') + ylab('Frequency') +
  ggtitle('Environmental Chemistry') + easy_center_title()+
  geom_vline(xintercept = c(1.64,1.96,2.58), lty=2, color=c(3,2,6), linewidth=0.5)+
  scale_x_continuous(breaks=c(0,1.64,1.96,2.58,4,6,8)) +
  theme(panel.background = element_rect(fill = "gray100"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

ene <- ggplot(datSUBF) +
  geom_line(aes(xs.ene, q025.ene), color='orange', lty=3) +
  geom_line(aes(xs.ene, n.cf.ene), color='orange', lty=1) +
  geom_point(aes(xs.ene, n.cf.ene),shape=20, fill='orange', color='orange', size=1)+
  geom_line(aes(xs.ene, q975.ene), color='orange', lty=3) +
  geom_line(aes(xs.ene, n.f.ene), color='blue', lty=2) +
  geom_point(aes(xs.ene, n.f.ene),shape=20, fill='blue', color='blue', size=1)+
  coord_cartesian(xlim=c(0,8)) + xlab('|z|-value') + ylab('Frequency') +
  ggtitle('Environmental Engineering') + easy_center_title()+
  geom_vline(xintercept = c(1.64,1.96,2.58), lty=2, color=c(3,2,6), linewidth=0.5)+
  scale_x_continuous(breaks=c(0,1.64,1.96,2.58,4,6,8)) +
  theme(panel.background = element_rect(fill = "gray100"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

nlc <- ggplot(datSUBF) +
  geom_line(aes(xs.nlc, q025.nlc), color='orange', lty=3) +
  geom_line(aes(xs.nlc, n.cf.nlc), color='orange', lty=1) +
  geom_point(aes(xs.nlc, n.cf.nlc),shape=20, fill='orange', color='orange', size=1)+
  geom_line(aes(xs.nlc, q975.nlc), color='orange', lty=3) +
  geom_line(aes(xs.nlc, n.f.nlc), color='blue', lty=2) +
  geom_point(aes(xs.nlc, n.f.nlc),shape=20, fill='blue', color='blue', size=1)+
  coord_cartesian(xlim=c(0,8)) + xlab('|z|-value') + ylab('Frequency') +
  ggtitle('Nature & Landscape Conservation') + easy_center_title()+
  geom_vline(xintercept = c(1.64,1.96,2.58), lty=2, color=c(3,2,6), linewidth=0.5)+
  scale_x_continuous(breaks=c(0,1.64,1.96,2.58,4,6,8)) +
  theme(panel.background = element_rect(fill = "gray100"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

mpl <- ggplot(datSUBF) +
  geom_line(aes(xs.mpl, q025.mpl), color='orange', lty=3) +
  geom_line(aes(xs.mpl, n.cf.mpl), color='orange', lty=1) +
  geom_point(aes(xs.mpl, n.cf.mpl),shape=20, fill='orange', color='orange', size=1)+
  geom_line(aes(xs.mpl, q975.mpl), color='orange', lty=3) +
  geom_line(aes(xs.mpl, n.f.mpl), color='blue', lty=2) +
  geom_point(aes(xs.mpl, n.f.mpl),shape=20, fill='blue', color='blue', size=1)+
  coord_cartesian(xlim=c(0,8)) + xlab('|z|-value') + ylab('Frequency') +
  ggtitle('Management, Monitoring, Policy & Law')+ easy_center_title()+
  geom_vline(xintercept = c(1.64,1.96,2.58), lty=2, color=c(3,2,6), linewidth=0.5)+
  scale_x_continuous(breaks=c(0,1.64,1.96,2.58,4,6,8)) +
  theme(panel.background = element_rect(fill = "gray100"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

wst <- ggplot(datSUBF) +
  geom_line(aes(xs.wst, q025.wst), color='orange', lty=3) +
  geom_line(aes(xs.wst, n.cf.wst), color='orange', lty=1) +
  geom_point(aes(xs.wst, n.cf.wst),shape=20, fill='orange', color='orange', size=1)+
  geom_line(aes(xs.wst, q975.wst), color='orange', lty=3) +
  geom_line(aes(xs.wst, n.f.wst), color='blue', lty=2) +
  geom_point(aes(xs.wst, n.f.wst),shape=20, fill='blue', color='blue', size=1)+
  coord_cartesian(xlim=c(0,8)) + xlab('|z|-value') + ylab('Frequency') +
  ggtitle('Water Science & Technology') + easy_center_title()+
  geom_vline(xintercept = c(1.64,1.96,2.58), lty=2, color=c(3,2,6), linewidth=0.5)+
  scale_x_continuous(breaks=c(0,1.64,1.96,2.58,4,6,8)) +
  theme(panel.background = element_rect(fill = "gray100"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

htm <- ggplot(datSUBF) +
  geom_line(aes(xs.htm, q025.htm), color='orange', lty=3) +
  geom_line(aes(xs.htm, n.cf.htm), color='orange', lty=1) +
  geom_point(aes(xs.htm, n.cf.htm),shape=20, fill='orange', color='orange', size=1)+
  geom_line(aes(xs.htm, q975.htm), color='orange', lty=3) +
  geom_line(aes(xs.htm, n.f.htm), color='blue', lty=2) +
  geom_point(aes(xs.htm, n.f.htm),shape=20, fill='blue', color='blue', size=1)+
  coord_cartesian(xlim=c(0,8)) + xlab('|z|-value') + ylab('Frequency') +
  ggtitle('Health, Toxicology & Mutagenesis') + easy_center_title()+
  geom_vline(xintercept = c(1.64,1.96,2.58), lty=2, color=c(3,2,6), linewidth=0.5)+
  scale_x_continuous(breaks=c(0,1.64,1.96,2.58,4,6,8)) +
  theme(panel.background = element_rect(fill = "gray100"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

pdf(here("results","robustness","Robustness_Figure_1_ESR_by_subfield.pdf"),width=14,height=12)
grid.arrange(eco, enc, ene, nlc, mpl, wst, htm, ncol = 2)
dev.off()


## ESR table for each subfield
# eco
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.eco <- cf(dat = mydat.eco, z.grid = p.grid.tab, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(p.tab.eco, file=here("results","robustness","pet_peese_rstandard_p_tab.eco.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ci.eco <- cf.ci.cluster(dat = mydat.eco, z.grid = p.grid.tab, iters = it, cluster = clu[i.eco], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  #about 44 mins for 1000 iterations.
save(p.tab.ci.eco, file=here("results","robustness","pet_peese_rstandard_p_tab_ci.eco.rds"))

# enc
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.enc <- cf(dat = mydat.enc, z.grid = p.grid.tab, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(p.tab.enc, file=here("results","robustness","pet_peese_rstandard_p_tab.enc.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ci.enc <- cf.ci.cluster(dat = mydat.enc, z.grid = p.grid.tab, iters = it, cluster = clu[i.enc], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  #about 18 mins for 1000 iterations.
save(p.tab.ci.enc, file=here("results","robustness","pet_peese_rstandard_p_tab_ci.enc.rds"))

# ene
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ene <- cf(dat = mydat.ene, z.grid = p.grid.tab, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(p.tab.ene, file=here("results","robustness","pet_peese_rstandard_p_tab.ene.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ci.ene <- cf.ci.cluster(dat = mydat.ene, z.grid = p.grid.tab, iters = it, cluster = clu[i.ene], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  #about 4 mins for 1000 iterations.
save(p.tab.ci.ene, file=here("results","robustness","pet_peese_rstandard_p_tab_ci.ene.rds"))

# nlc
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.nlc <- cf(dat = mydat.nlc, z.grid = p.grid.tab, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(p.tab.nlc, file=here("results","robustness","pet_peese_rstandard_p_tab.nlc.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ci.nlc <- cf.ci.cluster(dat = mydat.nlc, z.grid = p.grid.tab, iters = it, cluster = clu[i.nlc], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  #about 22 mins for 1000 iterations.
save(p.tab.ci.nlc, file=here("results","robustness","pet_peese_rstandard_p_tab_ci.nlc.rds"))

# mpl
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.mpl <- cf(dat = mydat.mpl, z.grid = p.grid.tab, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(p.tab.mpl, file=here("results","robustness","pet_peese_rstandard_p_tab.mpl.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ci.mpl <- cf.ci.cluster(dat = mydat.mpl, z.grid = p.grid.tab, iters = it, cluster = clu[i.mpl], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  #about 2 mins for 1000 iterations.
save(p.tab.ci.mpl, file=here("results","robustness","pet_peese_rstandard_p_tab_ci.mpl.rds"))

# wst
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.wst <- cf(dat = mydat.wst, z.grid = p.grid.tab, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(p.tab.wst, file=here("results","robustness","pet_peese_rstandard_p_tab.wst.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ci.wst <- cf.ci.cluster(dat = mydat.wst, z.grid = p.grid.tab, iters = it, cluster = clu[i.wst], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  #about 1 min for 1000 iterations.
save(p.tab.ci.wst, file=here("results","robustness","pet_peese_rstandard_p_tab_ci.wst.rds"))

# htm
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.htm <- cf(dat = mydat.htm, z.grid = p.grid.tab, heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
save(p.tab.htm, file=here("results","robustness","pet_peese_rstandard_p_tab.htm.rds"))

s.time <- Sys.time()
cl <- makeCluster(n_cores)
registerDoParallel(cl)
p.tab.ci.htm <- cf.ci.cluster(dat = mydat.htm, z.grid = p.grid.tab, iters = it, cluster = clu[i.htm], heterogeneity_multiplier = heterogeneity_multiplier)
stopCluster(cl)
e.time <- Sys.time()
print(e.time - s.time)  #about 12 mins for 1000 iterations.
save(p.tab.ci.htm, file=here("results","robustness","pet_peese_rstandard_p_tab_ci.htm.rds"))

## ---

load(here("results","robustness","pet_peese_rstandard_p_tab.eco.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab_ci.eco.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab.enc.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab_ci.enc.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab.ene.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab_ci.ene.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab.nlc.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab_ci.nlc.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab.mpl.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab_ci.mpl.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab.htm.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab_ci.htm.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab.wst.rds"))
load(here("results","robustness","pet_peese_rstandard_p_tab_ci.wst.rds"))

# eco
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")

N <- sum(p.orig.tab.eco)
p.table[,1] <- round((p.orig.tab.eco - p.tab.eco) / N,  3)

q025 <- apply(matrix(p.orig.tab.eco/N, nrow=nrow(p.tab.ci.eco[[1]]), ncol=ncol(p.tab.ci.eco[[1]]), byrow=TRUE)
              - p.tab.ci.eco[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.eco/N, nrow=nrow(p.tab.ci.eco[[1]]), ncol=ncol(p.tab.ci.eco[[1]]), byrow=TRUE)
              - p.tab.ci.eco[[1]], 2, quantile, probs=c(0.975))

p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.eco - p.tab.eco)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.eco[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.eco/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.eco/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.eco - p.tab.eco)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.eco[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.eco/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.eco/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.eco - p.tab.eco)[10:13] / sum(p.orig.tab.eco[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.eco[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.eco)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.eco/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.eco/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.eco - p.tab.eco)[11:13] / sum(p.orig.tab.eco[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.eco[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.eco)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.eco/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.eco/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(mydat.eco)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
                       "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
                       "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01",
                       "p < 0.001","ESR_{0.1}^{all}","ESR_{0.05}^{all}","ESR_{0.1}^{sig}","ESR_{0.05}^{sig}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)
write.csv(p.table, here("results","robustness","Robustness_Table_2_ESR_Ecology.csv"))

# enc
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")

N <- sum(p.orig.tab.enc)
p.table[,1] <- round((p.orig.tab.enc - p.tab.enc) / N,  3)

q025 <- apply(matrix(p.orig.tab.enc/N, nrow=nrow(p.tab.ci.enc[[1]]), ncol=ncol(p.tab.ci.enc[[1]]), byrow=TRUE)
              - p.tab.ci.enc[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.enc/N, nrow=nrow(p.tab.ci.enc[[1]]), ncol=ncol(p.tab.ci.enc[[1]]), byrow=TRUE)
              - p.tab.ci.enc[[1]], 2, quantile, probs=c(0.975))

p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.enc - p.tab.enc)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.enc[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.enc/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.enc/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.enc - p.tab.enc)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.enc[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.enc/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.enc/N)[11:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.enc - p.tab.enc)[10:13] / sum(p.orig.tab.enc[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.enc[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.enc)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.enc/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.enc/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.enc - p.tab.enc)[11:13] / sum(p.orig.tab.enc[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.enc[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.enc)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.enc/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.enc/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(mydat.enc)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
                       "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
                       "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01",
                       "p < 0.001","ESR_{0.1}^{all}","ESR_{0.05}^{all}","ESR_{0.1}^{sig}","ESR_{0.05}^{sig}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)
write.csv(p.table, here("results","robustness","Robustness_Table_3_ESR_Environmental_Chemistry.csv"))

# ene
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")

N <- sum(p.orig.tab.ene)
p.table[,1] <- round((p.orig.tab.ene - p.tab.ene) / N,  3)

q025 <- apply(matrix(p.orig.tab.ene/N, nrow=nrow(p.tab.ci.ene[[1]]), ncol=ncol(p.tab.ci.ene[[1]]), byrow=TRUE)
              - p.tab.ci.ene[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.ene/N, nrow=nrow(p.tab.ci.ene[[1]]), ncol=ncol(p.tab.ci.ene[[1]]), byrow=TRUE)
              - p.tab.ci.ene[[1]], 2, quantile, probs=c(0.975))

p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.ene - p.tab.ene)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.ene[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.ene/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.ene/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.ene - p.tab.ene)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.ene[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.ene/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.ene/N)[11:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.ene - p.tab.ene)[10:13] / sum(p.orig.tab.ene[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.ene[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.ene)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.ene/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.ene/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.ene - p.tab.ene)[11:13] / sum(p.orig.tab.ene[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.ene[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.ene)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.ene/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.ene/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(mydat.ene)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
                       "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
                       "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01",
                       "p < 0.001","ESR_{0.1}^{all}","ESR_{0.05}^{all}","ESR_{0.1}^{sig}","ESR_{0.05}^{sig}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)
write.csv(p.table, here("results","robustness","Robustness_Table_4_ESR_Environmental_Engineering.csv"))

# nlc
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")

N <- sum(p.orig.tab.nlc)
p.table[,1] <- round((p.orig.tab.nlc - p.tab.nlc) / N,  3)

q025 <- apply(matrix(p.orig.tab.nlc/N, nrow=nrow(p.tab.ci.nlc[[1]]), ncol=ncol(p.tab.ci.nlc[[1]]), byrow=TRUE)
              - p.tab.ci.nlc[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.nlc/N, nrow=nrow(p.tab.ci.nlc[[1]]), ncol=ncol(p.tab.ci.nlc[[1]]), byrow=TRUE)
              - p.tab.ci.nlc[[1]], 2, quantile, probs=c(0.975))

p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.nlc - p.tab.nlc)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.nlc[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.nlc/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.nlc/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.nlc - p.tab.nlc)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.nlc[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.nlc/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.nlc/N)[11:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.nlc - p.tab.nlc)[10:13] / sum(p.orig.tab.nlc[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.nlc[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.nlc)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.nlc/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.nlc/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.nlc - p.tab.nlc)[11:13] / sum(p.orig.tab.nlc[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.nlc[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.nlc)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.nlc/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.nlc/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(mydat.nlc)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
                       "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
                       "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01",
                       "p < 0.001","ESR_{0.1}^{all}","ESR_{0.05}^{all}","ESR_{0.1}^{sig}","ESR_{0.05}^{sig}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)
write.csv(p.table, here("results","robustness","Robustness_Table_5_ESR_Nature_and_Landscape_Conservation.csv"))

# mpl
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")

N <- sum(p.orig.tab.mpl)
p.table[,1] <- round((p.orig.tab.mpl - p.tab.mpl) / N,  3)

q025 <- apply(matrix(p.orig.tab.mpl/N, nrow=nrow(p.tab.ci.mpl[[1]]), ncol=ncol(p.tab.ci.mpl[[1]]), byrow=TRUE)
              - p.tab.ci.mpl[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.mpl/N, nrow=nrow(p.tab.ci.mpl[[1]]), ncol=ncol(p.tab.ci.mpl[[1]]), byrow=TRUE)
              - p.tab.ci.mpl[[1]], 2, quantile, probs=c(0.975))

p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.mpl - p.tab.mpl)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.mpl[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.mpl/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.mpl/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.mpl - p.tab.mpl)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.mpl[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.mpl/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.mpl/N)[11:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.mpl - p.tab.mpl)[10:13] / sum(p.orig.tab.mpl[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.mpl[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.mpl)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.mpl/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.mpl/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.mpl - p.tab.mpl)[11:13] / sum(p.orig.tab.mpl[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.mpl[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.mpl)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.mpl/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.mpl/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(mydat.mpl)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))

rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
                       "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
                       "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01",
                       "p < 0.001","ESR_{0.1}^{all}","ESR_{0.05}^{all}","ESR_{0.1}^{sig}","ESR_{0.05}^{sig}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)
write.csv(p.table, here("results","robustness","Robustness_Table_6_ESR_Management_Monitoring_Policy_and_Law.csv"))

# wst
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")

N <- sum(p.orig.tab.wst)
p.table[,1] <- round((p.orig.tab.wst - p.tab.wst) / N,  3)

q025 <- apply(matrix(p.orig.tab.wst/N, nrow=nrow(p.tab.ci.wst[[1]]), ncol=ncol(p.tab.ci.wst[[1]]), byrow=TRUE)
              - p.tab.ci.wst[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.wst/N, nrow=nrow(p.tab.ci.wst[[1]]), ncol=ncol(p.tab.ci.wst[[1]]), byrow=TRUE)
              - p.tab.ci.wst[[1]], 2, quantile, probs=c(0.975))

p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.wst - p.tab.wst)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.wst[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.wst/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.wst/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.wst - p.tab.wst)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.wst[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.wst/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.wst/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.wst - p.tab.wst)[10:13] / sum(p.orig.tab.wst[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.wst[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.wst)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.wst/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.wst/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)

p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.wst - p.tab.wst)[11:13] / sum(p.orig.tab.wst[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.wst[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.wst)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.wst/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.wst/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(mydat.wst)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
                       "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
                       "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01",
                       "p < 0.001","ESR_{0.1}^{all}","ESR_{0.05}^{all}","ESR_{0.1}^{sig}","ESR_{0.05}^{sig}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)
write.csv(p.table, here("results","robustness","Robustness_Table_7_ESR_Water_Science_and_Technology.csv"))

# htm
p.table <- matrix(ncol=2, nrow=length(p.grid.tab2)-1)
colnames(p.table) <- c("Difference", "0.95 Confidence interval")

N <- sum(p.orig.tab.htm)
p.table[,1] <- round((p.orig.tab.htm - p.tab.htm) / N,  3)

q025 <- apply(matrix(p.orig.tab.htm/N, nrow=nrow(p.tab.ci.htm[[1]]), ncol=ncol(p.tab.ci.htm[[1]]), byrow=TRUE)
              - p.tab.ci.htm[[1]], 2, quantile, probs=c(0.025))
q975 <- apply(matrix(p.orig.tab.htm/N, nrow=nrow(p.tab.ci.htm[[1]]), ncol=ncol(p.tab.ci.htm[[1]]), byrow=TRUE)
              - p.tab.ci.htm[[1]], 2, quantile, probs=c(0.975))

p.table[,2] <- paste("[", round(q025, 3), ", ",  round(q975, 3), "]", sep="")

#Share of inflated p-values in total p-values (10% level)
infp <- round(sum((p.orig.tab.htm - p.tab.htm)[10:13] / N),  3)
infp.bs <- apply(p.tab.ci.htm[[1]][,10:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.htm/N)[10:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.htm/N)[10:13]) - infp.bs, probs=c(0.975)),3)

p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Share of inflated p-values in total p-values (5% level)
infp <- round(sum((p.orig.tab.htm - p.tab.htm)[11:13] / N),  3)
infp.bs <- apply(p.tab.ci.htm[[1]][,11:13], 1, sum)
infp.q025 <- round(quantile(sum((p.orig.tab.htm/N)[11:13]) - infp.bs, probs=c(0.025)),3)
infp.q975 <- round(quantile(sum((p.orig.tab.htm/N)[11:13]) - infp.bs, probs=c(0.975)),3)
p.table <- rbind(p.table, c(infp, paste("[", infp.q025, ", ",  infp.q975, "]", sep="")))

#Inflated significance (10%)
shp.10 <- round(sum((p.orig.tab.htm - p.tab.htm)[10:13] / sum(p.orig.tab.htm[10:13])),  3)
infp.bs.10 <- apply(p.tab.ci.htm[[2]][,10:13], 1, sum)
N.sig.10 <- sum((p.orig.tab.htm)[10:13])
shp.q025.10 <- round(quantile((sum((p.orig.tab.htm/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.025)),3)
shp.q975.10 <- round(quantile((sum((p.orig.tab.htm/N.sig.10)[10:13]) - infp.bs.10) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.10, paste("[", shp.q025.10, ", ",  shp.q975.10, "]", sep="")))

#Inflated significance (5%)
shp.5 <- round(sum((p.orig.tab.htm - p.tab.htm)[11:13] / sum(p.orig.tab.htm[11:13])),  3)
infp.bs.5 <- apply(p.tab.ci.htm[[3]][,11:13], 1, sum) #[[3]] because this is the share in 5%
N.sig.5 <- sum((p.orig.tab.htm)[11:13])
shp.q025.5 <- round(quantile((sum((p.orig.tab.htm/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.025)),3)
shp.q975.5 <- round(quantile((sum((p.orig.tab.htm/N.sig.5)[11:13]) - infp.bs.5) , probs=c(0.975)),3)
p.table <- rbind(p.table, c(shp.5, paste("[", shp.q025.5, ", ",  shp.q975.5, "]", sep="")))

#Sample size
p.table <- rbind(p.table, c(length(summary(mydat.htm)[,1]), 0))
p.table <- rbind(p.table, c(N, 0))
rownames(p.table) <- c("0.9 < p", "0.8 < p < 0.9", "0.7 < p < 0.8", "0.6 < p < 0.7",
                       "0.5 < p < 0.6", "0.4 < p < 0.5", "0.3 < p < 0.4", "0.2 < p < 0.3",
                       "0.1 < p < 0.2", "0.05 < p < 0.1", "0.01 < p < 0.05", "0.001 < p < 0.01",
                       "p < 0.001","ESR_{0.1}^{all}","ESR_{0.05}^{all}","ESR_{0.1}^{sig}","ESR_{0.05}^{sig}",
                       "No. of meta-analysis", "No. of tests")
print(p.table)
write.csv(p.table, here("results","robustness","Robustness_Table_8_ESR_Health_Toxicology_and_Mutagenesis.csv"))
