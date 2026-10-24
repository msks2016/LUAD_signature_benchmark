# reporting.R - reporting completeness of the reviewed literature (scripts/41_reporting_completeness.R).

# Unweighted Cohen's kappa for two raters' categorical codes; NA when chance agreement is 1
# (both raters used a single, identical category, so kappa is not estimable).
cohen_kappa <- function(a, b) {
  stopifnot(length(a) == length(b))
  lv <- sort(unique(c(a, b)))
  a <- factor(a, levels = lv)
  b <- factor(b, levels = lv)
  po <- mean(a == b)
  pe <- sum(prop.table(table(a)) * prop.table(table(b)))
  if (pe >= 1 - 1e-12) return(NA_real_)
  (po - pe) / (1 - pe)
}

# Agreement of two coders on each item: n, share agreeing, kappa.
coder_agreement <- function(x, y, items, id = "paper_id") {
  m <- merge(x, y, by = id, suffixes = c("_1", "_2"))
  do.call(rbind, lapply(items, function(it) {
    a <- as.character(m[[paste0(it, "_1")]])
    b <- as.character(m[[paste0(it, "_2")]])
    data.frame(item = it, papers = nrow(m), agree = sum(a == b), agreement = mean(a == b),
               kappa = cohen_kappa(a, b), disagreeing_papers = paste(m[[id]][a != b], collapse = ";"),
               stringsAsFactors = FALSE)
  }))
}
