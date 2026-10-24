test_that("all-zero, constant and rarely expressed rows are handled", {
  expr <- rbind(ok   = c(0, 1.2, 3.1, 0.4),
                zero = c(0, 0, 0, 0),
                flat = c(2, 2, 2, 2),
                rare = c(0, 0, 0, 5))
  expect_equal(measurable_rows(expr), c("ok", "rare"))
  expect_equal(measurable_rows(expr, min_expressed_frac = 0.5), "ok")
})

test_that("coverage and pair genes follow the rule", {
  sig_map <- data.frame(signature = "s1", published_symbol = c("A", "B", "C", "D", "E"),
                        coefficient = c(1, -1, 0.5, 0.2, 0.1),
                        ensembl_id = c("E1", "E2", "E3", "E4", NA), stringsAsFactors = FALSE)
  present <- c("E1", "E2", "E3", "E4")
  cov <- coverage_table(sig_map, present, "cohortX")
  expect_equal(cov$n_present, 4)
  expect_equal(cov$frac_present, 0.8)
  expect_equal(cov$missing, "E")
  expect_null(pair_genes(sig_map, "s1", present, min_fraction = 1))
  expect_equal(nrow(pair_genes(sig_map, "s1", present, min_fraction = 0.8)), 4)
})

test_that("pairs skip cohorts the authors used, unless asked for the reference", {
  sig_map <- data.frame(signature = c("s1", "s1", "s2", "s2"),
                        published_symbol = c("A", "B", "C", "D"), coefficient = c(1, 1, 1, 1),
                        ensembl_id = c("E1", "E2", "E3", "E4"), stringsAsFactors = FALSE)
  measurable <- list(X = c("E1", "E2", "E3", "E4"))
  ind <- data.frame(signature = c("s1", "s2"), cohort = "X", used_by_authors = c("no", "yes"),
                    stringsAsFactors = FALSE)
  p <- build_pairs(sig_map, measurable, ind, "X", min_fraction = 1)
  expect_equal(p$pairs$signature, "s1")
  p2 <- build_pairs(sig_map, measurable, ind, "X", min_fraction = 1, include_used = TRUE)
  expect_equal(sort(p2$pairs$signature), c("s1", "s2"))
  ind$used_by_authors[2] <- "to_check"
  expect_error(build_pairs(sig_map, measurable, ind, "X", min_fraction = 1), "yes/no")
})

test_that("an alternative (successor) ID is used only where the original is not measurable", {
  gm <- data.frame(signature = c("S1", "S1", "S2"), published_symbol = c("A", "B", "C"),
                   coefficient = c(0.5, -0.2, 1), ensembl_id = c("E1", "E2", "E3"), stringsAsFactors = FALSE)
  alt <- data.frame(ensembl_id = c("E2", "E3"), alternative_id = c("E20", "E30"), stringsAsFactors = FALSE)
  m <- cohort_gene_map(gm, c("E1", "E20", "E3", "E30"), alt)
  expect_identical(m$ensembl_id, c("E1", "E20", "E3"))                 # E3 measurable: kept
  expect_identical(m$id_used, c("original", "alternative", "original"))
  m2 <- cohort_gene_map(gm, c("E1"), alt)                               # E20 not measurable either
  expect_identical(m2$ensembl_id, c("E1", "E2", "E3"))
  cov <- coverage_table(gm, c("E1", "E20"), "X", alternatives = alt)
  expect_equal(cov$frac_present[cov$signature == "S1"], 1)
  expect_identical(cov$via_alternative_id[cov$signature == "S1"], "B")
  dupl <- data.frame(ensembl_id = "E2", alternative_id = "E1", stringsAsFactors = FALSE)
  expect_error(cohort_gene_map(gm, c("E1"), dupl), "duplicates")
})

test_that("predecessor IDs follow the official history and the decision rule", {
  gm <- data.frame(signature = c("S1", "S1", "S1", "S2", "S2", "S3"),
                   published_symbol = c("A", "B", "C", "D", "E", "F"), coefficient = 1,
                   ensembl_id = c("N1", "N2", "N3", "N4", "N5", "N6"), stringsAsFactors = FALSE)
  alt <- data.frame(ensembl_id = "N5", alternative_id = "A5", stringsAsFactors = FALSE)
  arch <- data.frame(
    id = c("O1", "O2", "M2", "O3", "O4a", "O4b", "O5", "O6", "N1", "N2", "N3", "Z3", "N4", "A5", "N6"),
    current = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, rep(TRUE, 7)),
    replaced_by = c("N1", "M2", "N2", "N3;Z3", "N4", "N4", "A5", "N6", rep("", 7)),
    stringsAsFactors = FALSE)
  w <- history_walk("O2", stats::setNames(arch$current, arch$id), stats::setNames(arch$replaced_by, arch$id))
  expect_identical(w$path, c("O2", "M2", "N2"))
  expect_identical(w$ends, "N2")
  cohort_ids <- c("O1", "O2", "O3", "O4a", "O4b", "O5", "O6")
  meas <- c("O1", "O2", "O3", "O4a", "O4b", "O5")                   # O6 below the presence rule
  h <- predecessor_hits(gm, cohort_ids, meas, arch, alt)
  acc <- h[h$accept, ]
  expect_identical(sort(acc$cohort_id), c("O1", "O2"))              # direct and two-step history
  expect_match(h$note[h$cohort_id == "O3"], "split")
  expect_true(all(grepl("merge", h$note[h$ensembl_id == "N4"])))
  expect_match(h$note[h$cohort_id == "O5"], "alternative")
  expect_match(h$note[h$cohort_id == "O6"], "not measurable")
  h2 <- predecessor_hits(gm, cohort_ids, c(meas, "N1"), arch, alt)   # N1 measurable: no row for it
  expect_false("N1" %in% h2$ensembl_id)
})
