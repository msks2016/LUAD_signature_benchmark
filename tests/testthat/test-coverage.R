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
