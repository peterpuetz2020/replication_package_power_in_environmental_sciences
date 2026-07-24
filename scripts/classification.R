
## -----------------
## classification.R
## -----------------

## The code below was executed using "R version 3.5.2 (2018-12-20) -- Eggshell Igloo 64-bit"

rm(list = ls(all = TRUE))  # clear environment

## Load packages
## First, install "RTools34.exe" (https://cran.r-project.org/bin/windows/Rtools/history.html)
library(here)
library(tm)
library(stringr)
library(dplyr)
library(tidyr)
library(readxl)

## Import meta-articles. Note that 32 articles were dropped as they don't have list of references when imported from Scopus
my_data <- read_excel(here("data", "meta-articles imported from scopus.xlsx"), sheet = "original")
dim(my_data)

## Cleaning references
my_data$references <- gsub("\\.", "", my_data$references, perl = TRUE)
my_data$references <- tolower(my_data$references)
my_data$references <- stripWhitespace(my_data$references)

## Import all journal names (both full & abbreviated names) with their corresponding scimago category
scimago <- read_excel(here("data", "scimago_sjr.xlsx"))

## Cleaning the journal names
scimago$jname <- gsub("\\.", "", scimago$jname)
scimago$jname <- tolower(scimago$jname)
scimago$jname <- stripWhitespace(scimago$jname)

## Cleaning the categories names
scimago$sjr_category <- gsub("[()]", "", scimago$sjr_category)
scimago$sjr_category <- gsub("[[:digit:]]", "", scimago$sjr_category)
scimago$sjr_category <- gsub("Q", "", scimago$sjr_category)
scimago$sjr_category <- stripWhitespace(scimago$sjr_category)

## List of all journal names separated by "|"
pattern <- str_c(scimago$jname, collapse = "|")

## Identifying meta-articles which are not published in environmental science journals
id <- sjr_cat <- jnames <- freq <- nref <- nenvir <- perc <- list()
year <- journal <- references <- list()

ref_env <- data.frame(id = NA, year = NA, journal = NA, sjr_cat = NA, 
                      references = NA, nref = NA, nenvir = NA, perc = NA)

s.time <- Sys.time()
for (i in seq_along(my_data$references)) {
  print(i)
  
  nref[[i]] <- str_count(my_data$references[i], "\\(19[0-9]{2}\\)|\\(20[0-9]{2}\\)")
  jnames[[i]] <- unlist(
    str_extract_all(
      my_data$references[i],
      str_c("(?<=\\(?\\d{4}\\)?\\s)(", pattern, ")(?=,)")
    )
  )
  freq[[i]] <- table(jnames[[i]])
  nenvir[[i]] <- sum(freq[[i]])
  
  ref_env[i, 1] <- my_data$id[i]
  ref_env[i, 2] <- my_data$year[i]
  ref_env[i, 3] <- my_data$journal[i]
  ref_env[i, 4] <- my_data$sjr_category[i]
  ref_env[i, 5] <- my_data$references[i]
  ref_env[i, 6] <- nref[[i]]
  ref_env[i, 7] <- nenvir[[i]]
  ref_env[i, 8] <- ifelse(nref[[i]] > 0, nenvir[[i]] / nref[[i]] * 100, NA)
}
e.time <- Sys.time()
print(e.time - s.time)  # about 15 minutes

length(which(ref_env$nenvir==0))
#ref_env0 <- ref_env %>% filter(nenvir==0)

## Dropping meta-articles that are not published in environmental sciences journal
ref_env1 <- ref_env %>% filter(nenvir != 0)
round(median(ref_env1$perc), digits = 0)

hist(ref_env1$perc, breaks = 60, xlab = "Meta-articles published in environmental sciences journal (in %)",
  main = "",  cex.lab = 0.9,  cex.axis = 0.85,  col = "lightblue")
abline(v = 25, col = "red", lty = "dashed")

## We used the median of the percentages as a cut-off point (25%)
ref_env_final <- ref_env1 %>% filter(perc >= 25)
write.csv(ref_env_final, here("data", "ref_env_percentage.csv"), row.names = FALSE)

## Meta-classification into subfields
refenv <- read.csv(here("data", "ref_env_percentage.csv"), header = TRUE, sep = ",")
dim(refenv)

## Identifying the most frequent subfield for each meta-article
id <- title <- year <- jnames <- sjr_cat <- refs <- list()
subfield <- nref <- nenvir <- perc <- freq_cat1 <- list()
df1 <- df2 <- df3 <- df4 <- list()

subMeta <- data.frame(
  id = NA, title = NA, year = NA, jnames = NA, sjr_cat = NA, refs = NA,
  nref = NA, nenvir = NA, perc = NA, freq_cat1 = NA)

s.time <- Sys.time()
for (i in seq_along(refenv$references)) {
  print(i)
  jnames[[i]] <- unlist(
    str_extract_all(
      refenv$references[i],
      str_c("(?<=\\(?\\d{4}\\)?\\s)(", pattern, ")(?=,)")
    )
  )
  pos <- match(jnames[[i]], scimago$jname)
  subfield[[i]] <- scimago$sjr_category[pos]
  
  df1[[i]] <- structure(
    list(var1 = subfield[[i]]),
    class = "data.frame",
    row.names = c(NA, length(subfield[[i]]))
  )
  
  df2[[i]] <- df1[[i]] %>%
    mutate(rn = row_number()) %>%
    separate_rows(var1, sep = "\\s*;\\s*") %>%
    add_count(rn) %>%
    mutate(n = 1 / n) %>%
    group_by(var1) %>%
    summarise(n = sum(n), .groups = "drop")
  
  df3[[i]] <- df2[[i]] %>%
    group_by(n) %>%
    summarise(var1 = paste(var1, collapse = "|"), .groups = "drop")
  
  df4[[i]] <- df3[[i]][order(df3[[i]]$n, decreasing = TRUE), ]
  freq_cat1[[i]] <- df4[[i]]$var1[1]
  
  subMeta[i, 1] <- refenv$id[i]
  subMeta[i, 2] <- as.character(refenv$title[i])
  subMeta[i, 3] <- refenv$year[i]
  subMeta[i, 4] <- as.character(refenv$journal[i])
  subMeta[i, 5] <- as.character(refenv$sjr_cat[i])
  subMeta[i, 6] <- as.character(refenv$references[i])
  subMeta[i, 7] <- refenv$nref[i]
  subMeta[i, 8] <- refenv$nenvir[i]
  subMeta[i, 9] <- refenv$perc[i]
  subMeta[i, 10] <- freq_cat1[[i]]
}
e.time <- Sys.time()
print(e.time - s.time)

write.csv(subMeta, here("data", "meta-classified into subfields.csv"), row.names = FALSE)

## When I'm unable to resolve the tie issue even after adding the Scimago category to the most frequent categories,
## I resort to manual inspection for few meta-papers (7) and decide the subfield myself.
