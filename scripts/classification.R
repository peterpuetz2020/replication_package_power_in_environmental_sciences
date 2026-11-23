
## -----------------
## classification.R
## -----------------



## Load packages
library(here)
library(tm)
library(stringr)
library(dplyr)
library(tidyr)
library(readxl)

## Ensure generated classification files can be written on a fresh checkout.
dir.create(here("data"), recursive = TRUE, showWarnings = FALSE)

## Import meta-articles. Note that 32 articles were dropped as they don't have list of references when imported from Scopus
meta_articles <- read_excel(here("data", "meta-articles.xlsx"))
dim(meta_articles)

## Import full and abbreviated journal names with their Scimago categories.
journal_categories <- read_excel(here("data", "scimago_sjr.xlsx"))

## Cleaning the journal names
journal_categories$jname <- gsub("\\.", "", journal_categories$jname)
journal_categories$jname <- tolower(journal_categories$jname)
journal_categories$jname <- stripWhitespace(journal_categories$jname)

## Cleaning the categories names
journal_categories$sjr_category <- gsub("[()]", "", journal_categories$sjr_category)
journal_categories$sjr_category <- gsub("[[:digit:]]", "", journal_categories$sjr_category)
journal_categories$sjr_category <- gsub("Q", "", journal_categories$sjr_category)
journal_categories$sjr_category <- stripWhitespace(journal_categories$sjr_category)

## List of all journal names separated by "|"
journal_name_pattern <- str_c(journal_categories$jname, collapse = "|")

## Allocate the per-article intermediate lists once, rather than growing them in
## the loop. Descriptive object names distinguish counts from percentages.
matched_journal_names <- journal_frequencies <- reference_count <- list()
environmental_reference_count <- list()

## Preallocate one row per article. Extending an individual data-frame column with
## article_reference_summary$id[i] fails as soon as i exceeds the data frame's current row count.
article_reference_summary <- data.frame(
  id = meta_articles$id,
  authors = meta_articles$authors,
  title = meta_articles$title,
  year = meta_articles$year,
  journal = meta_articles$journal,
  sjr_cat = meta_articles$sjr_category,
  references = meta_articles$references,
  nref = rep(NA_integer_, nrow(meta_articles)),
  nenvir = rep(NA_integer_, nrow(meta_articles)),
  perc = rep(NA_real_, nrow(meta_articles))
)

start_time <- Sys.time()
for (i in seq_along(meta_articles$references)) {
  print(i)
  
  reference_count[[i]] <- str_count(meta_articles$references[i], "\\(19[0-9]{2}\\)|\\(20[0-9]{2}\\)")
  matched_journal_names[[i]] <- unlist(
    str_extract_all(
      meta_articles$references[i],
      str_c("(?<=\\(?\\d{4}\\)?\\s)(", journal_name_pattern, ")(?=,)")
    )
  )
  journal_frequencies[[i]] <- table(matched_journal_names[[i]])
  environmental_reference_count[[i]] <- sum(journal_frequencies[[i]])
  article_reference_summary$nref[i] <- reference_count[[i]]
  article_reference_summary$nenvir[i] <- environmental_reference_count[[i]]
  article_reference_summary$perc[i] <- ifelse(
    reference_count[[i]] > 0,
    environmental_reference_count[[i]] / reference_count[[i]] * 100,
    NA_real_
  )
}
end_time <- Sys.time()
print(end_time - start_time)  # about 15 minutes

sum(article_reference_summary$nenvir == 0, na.rm = TRUE)

## Dropping meta-articles that are not published in environmental sciences journal
articles_with_environmental_references <- article_reference_summary %>%
  filter(nenvir != 0)
round(median(articles_with_environmental_references$perc), digits = 0)

hist(articles_with_environmental_references$perc, breaks = 60, xlab = "Meta-articles published in environmental sciences journal (in %)",
  main = "",  cex.lab = 0.9,  cex.axis = 0.85,  col = "lightblue")
abline(v = 25, col = "red", lty = "dashed")

## We used the median of the percentages as a cut-off point (25%)
environmental_articles <- articles_with_environmental_references %>%
  filter(perc >= 25)
write.csv(environmental_articles, here("data", "ref_env_percentage.csv"), row.names = FALSE)

## Meta-classification into subfields
environmental_article_references <- read.csv(here("data", "ref_env_percentage.csv"), header = TRUE, sep = ",")

## Older generated files did not contain titles. Add an explicit value for every
## row so this stage also works when it is run on its own with the supplied CSV.
if (!"title" %in% names(environmental_article_references)) {
  environmental_article_references$title <- rep(NA_character_, nrow(environmental_article_references))
}
dim(environmental_article_references)

## Identify the most frequent subfield for each retained meta-article. Each
## journal's weight is divided equally among all categories assigned to it.
matched_journal_names <- matched_subfields <- leading_subfield <- list()
subfield_rows <- weighted_subfield_counts <- tied_subfield_counts <- ranked_subfields <- list()

classified_meta_articles <- data.frame(
  id = NA, title = NA, year = NA, jnames = NA, sjr_cat = NA, refs = NA,
  nref = NA, nenvir = NA, perc = NA, freq_cat1 = NA
)

start_time <- Sys.time()
for (i in seq_along(environmental_article_references$references)) {
  print(i)
  matched_journal_names[[i]] <- unlist(
    str_extract_all(
      environmental_article_references$references[i],
      str_c("(?<=\\(?\\d{4}\\)?\\s)(", journal_name_pattern, ")(?=,)")
    )
  )
  journal_match_positions <- match(
    matched_journal_names[[i]], journal_categories$jname
  )
  matched_subfields[[i]] <-
    journal_categories$sjr_category[journal_match_positions]
  
  subfield_rows[[i]] <- structure(
    list(var1 = matched_subfields[[i]]),
    class = "data.frame",
    row.names = c(NA, length(matched_subfields[[i]]))
  )
  
  weighted_subfield_counts[[i]] <- subfield_rows[[i]] %>%
    mutate(rn = row_number()) %>%
    separate_rows(var1, sep = "\\s*;\\s*") %>%
    add_count(rn) %>%
    mutate(n = 1 / n) %>%
    group_by(var1) %>%
    summarise(n = sum(n), .groups = "drop")
  
  tied_subfield_counts[[i]] <- weighted_subfield_counts[[i]] %>%
    group_by(n) %>%
    summarise(var1 = paste(var1, collapse = "|"), .groups = "drop")
  
  ranked_subfields[[i]] <- tied_subfield_counts[[i]][order(tied_subfield_counts[[i]]$n, decreasing = TRUE), ]
  leading_subfield[[i]] <- ranked_subfields[[i]]$var1[1]
  
  classified_meta_articles[i, 1] <- environmental_article_references$id[i]
  classified_meta_articles[i, 2] <- as.character(environmental_article_references$title[i])
  classified_meta_articles[i, 3] <- environmental_article_references$year[i]
  classified_meta_articles[i, 4] <- as.character(environmental_article_references$journal[i])
  classified_meta_articles[i, 5] <- as.character(environmental_article_references$sjr_cat[i])
  classified_meta_articles[i, 6] <- as.character(environmental_article_references$references[i])
  classified_meta_articles[i, 7] <- environmental_article_references$nref[i]
  classified_meta_articles[i, 8] <- environmental_article_references$nenvir[i]
  classified_meta_articles[i, 9] <- environmental_article_references$perc[i]
  classified_meta_articles[i, 10] <- leading_subfield[[i]]
}
end_time <- Sys.time()
print(end_time - start_time)

write.csv(classified_meta_articles, here("data", "meta-classified into subfields.csv"), row.names = FALSE)
