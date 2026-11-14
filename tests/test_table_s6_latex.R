## Table S6 must carry math labels through direct and manuscript LaTeX writes.
library(dplyr)
library(tidyr)
supplement <- as.list(parse('scripts/create_supplement_tables_and_figures.R'))
main <- as.list(parse('scripts/create_tables_and_figures.R'))
assignment <- function(expressions, name) {
  matches <- Filter(function(expr) is.call(expr) && identical(expr[[1]], as.name('<-')) &&
    identical(expr[[2]], as.name(name)), expressions)
  stopifnot(length(matches) == 1L)
  matches[[1]]
}
write_call <- Filter(function(expr) is.call(expr) &&
  identical(expr[[1]], as.name('write_word_table')) &&
  identical(expr[[2]], as.name('table_s6')), supplement)[[1]]
fallback <- Filter(function(expr) is.call(expr) && identical(expr[[1]], as.name('if')) &&
  any(grepl('write_latex_table', deparse(expr[[2]]), fixed = TRUE)), supplement)[[1]]
labels <- c("Cohen's $d$", 'Correlation coefficient ($r$)', "Fisher's $z$",
            "Hedges's $g$", 'Regression coefficient ($\\beta$)')
root <- tempfile('table-s6-math-')
dir.create(root)
for (mode in c('main', 'standalone')) {
  session <- new.env(parent = .GlobalEnv)
  session$supplement_dir <- file.path(root, mode)
  dir.create(session$supplement_dir)
  if (mode == 'main') {
    eval(assignment(main, 'latex_escape'), session)
    eval(assignment(main, 'write_latex_table'), session)
  } else {
    eval(fallback, session)
  }
  eval(assignment(supplement, 'write_word_table'), session)
  eval(assignment(supplement, 'subfield_levels'), session)
  types <- c("cohen's d", 'correlation', 'excess risk', "fisher's z", "hedge's g",
    'logHR', 'logOR', 'logRR', 'lnRR', 'percentage change', 'mean', 'regression coefficient')
  session$table_s6_all_data <- data.frame(cID = seq_along(types), etype = types,
    subfd = 'Ecology')
  eval(assignment(supplement, 'table_s6'), session)
  stopifnot(all(labels %in% session$table_s6$`Effect size`))
  eval(write_call, session)
  latex <- readLines(file.path(session$supplement_dir, 'Table_S6.tex'))
  stopifnot(all(vapply(labels, function(label) any(grepl(label, latex, fixed = TRUE)), logical(1))),
    !any(grepl('β', latex, fixed = TRUE)))
  direct_path <- file.path(session$supplement_dir, 'direct.tex')
  session$write_latex_table(session$table_s6, direct_path)
  stopifnot(identical(readLines(direct_path), latex))
  unzip(file.path(session$supplement_dir, 'Table_S6.docx'), 'word/document.xml',
    exdir = session$supplement_dir)
  word <- xml2::xml_text(xml2::read_xml(file.path(session$supplement_dir, 'word/document.xml')))
  stopifnot(grepl('Regression coefficient (β)', word, fixed = TRUE),
    !grepl('$', word, fixed = TRUE))
}
unlink(root, recursive = TRUE)
cat('PASS: Table S6 math labels survive direct, main-workflow and standalone LaTeX export; Word labels remain readable.\n')
