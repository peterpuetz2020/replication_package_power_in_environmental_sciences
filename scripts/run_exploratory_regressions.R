## Exploratory regression models and diagnostics.
## ----------------------
## Exploratory regression
## ----------------------

## Load the data
first_existing <- function(...) {
  candidates <- c(...)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("None of the expected analysis data files exists: ", paste(candidates, collapse = ", "))
  }
  existing[[1]]
}

esr_path <- first_existing(
  here("results", "main", "esr05_pet_peese_rstandard_half meta-average_704.xlsx"),
  here("data", "esr05_pet_peese_rstandard_half meta-average_704.xlsx")
)
esr <- read_excel(esr_path)
#esr <- read_excel(here("results","main","esr05_pet_peese_rstandard_half meta-average_0.25xtau2_704.xlsx"))
#esr <- read_excel(here("results","main","esr05_pet_peese_rstandard_full meta-average_704.xlsx"))
dim(esr) #704x8

power_path <- first_existing(
  here("results", "main", "median_power_pps_rstandard_half meta-average_704.xlsx"),
  here("data", "median_power_pps_rstandard_half meta-average_704.xlsx")
)
pow <- read.xlsx(power_path)
#pow <- read.xlsx(here("results","main","median_power_pps_rstandard_full meta-average_704.xlsx"))
dim(pow) #704x15

merged <- merge(esr, pow, by="cID")
merged$sape_perc <- 100*merged$sape
merged$med_perc <- 100*merged$median
summary(merged$esr.sig)
summary(merged$esr.all)
length(which(merged$esr.sig < 0)) #74 (half meta-average) #23 (half + 25% tau2)

final <- merged
final$lognps <- log(final$nips)
final$logtotsig <- log(final$tot.sig+0.5)  #add 0.5 b/c log(0) is -Inf.
final$logtotall <- log(final$tot.all)
final$logjif <- log(final$jif_5yr_wos)

## Re-code the categorical variables
final <- final %>%
  dplyr::mutate(subfield = case_when(subf=="Health, Toxicology and Mutagenesis" ~ 0,
                                     subf=="Ecology" ~ 1,
                                     subf=="Environmental Chemistry" ~ 2,
                                     subf=="Environmental Engineering" ~ 3,
                                     subf=="Management, Monitoring, Policy and Law" ~ 4,
                                     subf=="Nature and Landscape Conservation" ~ 5,
                                     subf=="Water Science and Technology" ~ 6),
                design = case_when(sdes=="observational" ~ 0,
                                   sdes=="experimental" ~ 1,
                                   sdes=="mixed" ~ 2),
                #design_merged = case_when(sdes=="observational" | sdes=="mixed" ~ 0,
                #                   sdes=="experimental" ~ 1),
                design_merged = case_when(sdes=="observational" | sdes=="mixed" ~ "no",
                                          sdes=="experimental" ~ "yes"), #for plot
                design_plot = case_when(sdes=="observational" | sdes=="mixed" ~ "Observational",
                                        sdes=="experimental" ~ "Experimental"),
                guidline = case_when(guid=="yes" ~ 1, guid=="no" ~ 0),
                prereg = case_when(prer=="yes" ~ 1, prer=="no" ~ 0),
                metric = case_when(esty=="lnRR" | esty=="log-mean ratio" | esty=="ratio" ~ 0,
                                   esty=="cohen's d" ~ 1,
                                   esty=="hedge's g" ~ 2,
                                   esty=="correlation" ~ 3,
                                   esty=="fisher's z" ~ 4,
                                   esty=="logOR" ~ 5,
                                   esty=="logRR" ~ 6,
                                   esty=="logHR" ~ 7,
                                   esty=="excess risk" ~ 8,
                                   esty=="regression coefficient" ~ 9,
                                   esty=="percentage change" ~ 10,
                                   esty=="mean" | esty=="mean difference" ~ 11))


## Fitting NB regression
final_nb <- final %>% filter(esr.sig.count >= 0)
dim(final_nb) #630x35
final_nb$metric <- as.character(final_nb$metric)

## Negative binomial regression ("logtotsig" is used as an offset variable).
nbMod1 <- glm.nb(esr.sig.count ~ med_perc + as.factor(design_merged) + guid + prer
                 + lognps + logjif + pyear + metric + offset(logtotsig), data=final_nb)
nbMod1.robu <- coeftest(nbMod1, vcov.=vcovCL(nbMod1, cluster=final_nb$metaID))

nbMod2 <- glm.nb(esr.sig.count ~ med_perc + as.factor(design_merged) + guid + prer
                 + lognps + logjif + pyear + metric + subf + offset(logtotsig), data=final_nb)
nbMod2.robu <- coeftest(nbMod2, vcov.=vcovCL(nbMod2, cluster=final_nb$metaID))

## With clustering
stargazer(nbMod1.robu, nbMod2.robu, type="latex",style="demography",
          ci=F,star.cutoffs=c(.1,.05,.01),font.size="small",
          no.space=T, intercept.bottom=F,notes.align="l",
          column.separate=c(0.5, 0.5, 0.5),single.row = T)

## Model evaluation for NB model (based Model 1)
final_nb$resid <- nbMod1$residuals
final_nb$fitted <- nbMod1$fitted.values

## Plots

pdf(here("results","robustness","Robustness_Figure_2_NB_Model_1_continuous_diagnostics.pdf"),width=12,height=12)
par(mfrow=c(4,2))
plot(final_nb$fitted, final_nb$resid,col="skyblue",
     xlab="Fitted values", ylab="Residuals", main="Residuals vs. Fitted")
abline(h=0, col="red", lty="dashed")
qqnorm(final_nb$resid, col="skyblue")
qqline(final_nb$resid)

plot(final_nb$resid ~ final_nb$median, xlab="Fitted", ylab="Residuals",main="Residuals vs. Median Power",col="skyblue")
abline(h=0, lty="dashed", col="red")
plot(final_nb$resid ~ final_nb$lognps, xlab="Fitted", ylab="Residuals",main="Residuals vs. log(NPRS)",col="skyblue")
abline(h=0, lty="dashed", col="red")
plot(final_nb$resid ~ final_nb$logtotall, xlab="Fitted", ylab="Residuals",main="Residuals vs. log(NTPV)",col="skyblue")
abline(h=0, lty="dashed", col="red")
plot(final_nb$resid ~ final_nb$logjif, xlab="Fitted", ylab="Residuals",main="Residuals vs. log(JIF)",col="skyblue")
abline(h=0, lty="dashed", col="red")
plot(final_nb$resid ~ final_nb$pyear, xlab="Fitted", ylab="Residuals",main="Residuals vs. Publication Year",col="skyblue")
abline(h=0, lty="dashed", col="red")
dev.off()
par(mfrow=c(1,1))

## Plots for residuals vs. categorical variables
exp <- ggplot(data=final_nb, aes(x=as.factor(design_merged),y=resid)) +
  geom_boxplot() +
  geom_hline(yintercept=0,color="red") +
  labs(x = "Experimental research design?", y="Residuals") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

guid <- ggplot(data=final_nb, aes(x=as.factor(guid),y=resid)) +
  geom_boxplot() +
  geom_hline(yintercept=0,color="red") +
  labs(x = "Followed reporting guidelines?", y="Residuals") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

prer <- ggplot(data=final_nb, aes(x=as.factor(prer),y=resid)) +
  geom_boxplot() +
  geom_hline(yintercept=0,color="red") +
  labs(x = "Protocol registered?", y="Residuals") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

subf <- ggplot(data=final_nb, aes(x=as.factor(subfield),y=resid)) +
  geom_boxplot() +
  geom_hline(yintercept=0,color="red") +
  labs(x = "Subfield", y="Residuals") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

pdf(here("results","robustness","Robustness_Figure_3_NB_Model_1_categorical_diagnostics.pdf"),width=12,height=6)
grid.arrange(exp, guid, prer, subf, nrow=2, ncol = 2)
dev.off()


## Model evaluation for NB model (based Model 2)
final_nb$resid <- nbMod2$residuals
final_nb$fitted <- nbMod2$fitted.values

## Plots

pdf(here("results","robustness","Robustness_Figure_4_NB_Model_2_continuous_diagnostics.pdf"),width=12,height=12)
par(mfrow=c(4,2))
plot(final_nb$fitted, final_nb$resid,col="skyblue",
     xlab="Fitted values", ylab="Residuals", main="Residuals vs. Fitted")
abline(h=0, col="red", lty="dashed")
qqnorm(final_nb$resid, col="skyblue")
qqline(final_nb$resid)

plot(final_nb$resid ~ final_nb$median, xlab="Fitted", ylab="Residuals",main="Residuals vs. Median Power",col="skyblue")
abline(h=0, lty="dashed", col="red")
plot(final_nb$resid ~ final_nb$lognps, xlab="Fitted", ylab="Residuals",main="Residuals vs. log(NPRS)",col="skyblue")
abline(h=0, lty="dashed", col="red")
plot(final_nb$resid ~ final_nb$logtotall, xlab="Fitted", ylab="Residuals",main="Residuals vs. log(NTPV)",col="skyblue")
abline(h=0, lty="dashed", col="red")
plot(final_nb$resid ~ final_nb$logjif, xlab="Fitted", ylab="Residuals",main="Residuals vs. log(JIF)",col="skyblue")
abline(h=0, lty="dashed", col="red")
plot(final_nb$resid ~ final_nb$pyear, xlab="Fitted", ylab="Residuals",main="Residuals vs. Publication Year",col="skyblue")
abline(h=0, lty="dashed", col="red")
dev.off()
par(mfrow=c(1,1))

## Plots for residuals vs. categorical variables
exp <- ggplot(data=final_nb, aes(x=as.factor(design_merged),y=resid)) +
  geom_boxplot() +
  geom_hline(yintercept=0,color="red") +
  labs(x = "Experimental research design?", y="Residuals") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

guid <- ggplot(data=final_nb, aes(x=as.factor(guid),y=resid)) +
  geom_boxplot() +
  geom_hline(yintercept=0,color="red") +
  labs(x = "Followed reporting guidelines?", y="Residuals") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

prer <- ggplot(data=final_nb, aes(x=as.factor(prer),y=resid)) +
  geom_boxplot() +
  geom_hline(yintercept=0,color="red") +
  labs(x = "Protocol registered?", y="Residuals") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

subf <- ggplot(data=final_nb, aes(x=as.factor(subfield),y=resid)) +
  geom_boxplot() +
  geom_hline(yintercept=0,color="red") +
  labs(x = "Subfield", y="Residuals") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(linewidth = 0.5, color = "gray"))

pdf(here("results","robustness","Robustness_Figure_5_NB_Model_2_categorical_diagnostics.pdf"),width=12,height=6)
grid.arrange(exp, guid, prer, subf, nrow=2, ncol = 2)
dev.off()

## Weighted negative binomial regression
## Define weights
final_nb$wt <- 1 / final_nb$nips

nbMod3 <- glm.nb(esr.sig.count ~ med_perc + as.factor(design_merged) + guid + prer
                 + logjif + pyear + as.factor(metric) + offset(logtotsig),
                 weights = final_nb$wt, data=final_nb)
nbMod3.clu <- coeftest(nbMod3, vcov.=vcovCL(nbMod3, cluster=final_nb$metaID))

nbMod4 <- glm.nb(esr.sig.count ~ med_perc + as.factor(design_merged) + guid + prer
                 + logjif + pyear + as.factor(metric) + subf + offset(logtotsig),
                 weights = final_nb$wt, data=final_nb)
nbMod4.clu <- coeftest(nbMod4, vcov.=vcovCL(nbMod4, cluster=final_nb$metaID))

stargazer(nbMod3.clu, nbMod4.clu, type="text",style="demography",
          ci=F,star.cutoffs=c(.1,.05,.01),font.size="small",
          no.space=T, intercept.bottom=F,notes.align="l",
          column.separate=c(0.5, 0.5, 0.5),single.row = T)

## -----------------------------------------------------------------------------

## Winsorize the dependent variable
esr.sig_win <- Winsorize(final$esr.sig, probs = c(0.05, 0.95))
final_win <- cbind.data.frame(final, esr.sig_win)
summary(final_win$esr.sig_win)
summary(final_win$esr.sig)

## Correlation between median power and research design
cor.test(final$median, final$design_merged) #-0.05442532

## OLS regression (before and after winsorizing dependent variable)
olsMod1 <- lm(esr.sig_win ~ median + as.factor(design_merged) + guid + prer
              + lognps + logjif + pyear + as.factor(metric), data=final_win)

olsMod2 <- lm(esr.sig_win ~ median + as.factor(design_merged) + guid + prer
              + lognps + logjif + pyear + as.factor(metric) + subf, data=final_win)
zz <- round(vif(olsMod2),2) #variance inflation factor
xtable(zz)

## Clustering at "metaID"
olsMod1.clu <- coeftest(olsMod1, vcov.=vcovCL(olsMod1, cluster=final_win$metaID))
olsMod2.clu <- coeftest(olsMod2, vcov.=vcovCL(olsMod2, cluster=final_win$metaID))

## With clustering
stargazer(olsMod1.clu, olsMod2.clu,type="latex",style="demography",
          ci=F,star.cutoffs=c(.1,.05,.01),font.size="small",
          no.space=T, intercept.bottom=F,notes.align="l",
          column.separate=c(0.5, 0.5, 0.5),single.row = T)

olsMod_mixed <- lm(esr.sig_win ~ median + as.factor(sdes) + guid + prer
                   + lognps + logjif + pyear + as.factor(metric) + subf, data=final_win)
olsMod_mixed.clu <- coeftest(olsMod_mixed, vcov.=vcovCL(olsMod_mixed, cluster=final_win$metaID))
stargazer(olsMod_mixed.clu,type="latex",style="demography",
          ci=F,star.cutoffs=c(.1,.05,.01),font.size="small",
          no.space=T, intercept.bottom=F,notes.align="l",
          column.separate=c(0.5, 0.5, 0.5),single.row = T)

## Restrict the sample to where esr.sig >= 0
fin <- final %>% filter(esr.sig >= 0)
olsMod3 <- lm(esr.sig ~ median + as.factor(design_merged) + guid + prer
              + lognps + logjif + pyear + as.factor(metric), data=fin)

olsMod4 <- lm(esr.sig ~ median + as.factor(design_merged) + guid + prer
              + lognps + logjif + pyear + as.factor(metric) + subf, data=fin)

## Clustering at "metaID"
olsMod3.clu <- coeftest(olsMod3, vcov.=vcovCL(olsMod3, cluster=fin$metaID))
olsMod4.clu <- coeftest(olsMod4, vcov.=vcovCL(olsMod4, cluster=fin$metaID))

stargazer(olsMod3.clu, olsMod4.clu,type="latex",style="demography",
          ci=F,star.cutoffs=c(.1,.05,.01),font.size="small",
          no.space=T, intercept.bottom=F,notes.align="l",
          column.separate=c(0.5, 0.5, 0.5),single.row = T)
