test_that("null test returns the observed C-index and a valid p-value", {
  co <- sim_cohort(n = 200, n_genes = 300, beta = 0.8, seed = 11)
  genes <- rownames(co$expr)[1:3]
  res <- expression_matched_null(co$expr, genes, c(1, 1, 1), co$time, co$event,
                                 pool = rownames(co$expr), n_draws = 199, seed = 12)
  expect_equal(res$observed, cindex(co$time, co$event,
                                    score_signature(co$expr, genes, c(1, 1, 1))))
  expect_true(res$p_two_sided > 0 && res$p_two_sided <= 1)
  expect_true(res$p_two_sided < 0.05)   # true signal genes beat random genes
})

test_that("H1 classification follows the registration", {
  expect_equal(classify_h1(c(0.65, 0.35, 0.60), c(0.01, 0.01, 0.20)),
               c("Replicated", "Contradicted", "Inconclusive"))
})

test_that("a pair is scored with its own cohort (guards against bug R1)", {
  a <- sim_cohort(seed = 1)
  genes <- rownames(a$expr)[1:3]
  coefs <- c(0.4, 0.3, 0.2)
  direct <- hr_per_sd(a$time, a$event, score_signature(a$expr, genes, coefs))[["log_hr"]]
  expect_equal(pair_log_hr(a, genes, coefs), direct)
  idx <- sample.int(length(a$time), replace = TRUE)            # resampled patients, as in the bootstrap
  full <- hr_per_sd(a$time[idx], a$event[idx],
                    score_signature(a$expr[, idx, drop = FALSE], genes, coefs))[["log_hr"]]
  expect_equal(pair_log_hr(a, genes, coefs, idx), full, tolerance = 1e-10)   # genes first + coxph.fit: same number
  expect_error(pair_log_hr(a, c(genes[1:2], "ENSG_NOT_THERE"), coefs), "Not in matrix")
})

test_that("the bootstrap CI describes the pooled mean, not the spread of pairs (bug R2)", {
  co <- list(A = sim_cohort(n = 500, beta = 0.8, seed = 3),
             B = sim_cohort(n = 500, beta = 0.0, seed = 4))
  g <- data.frame(ensembl_id = rownames(co$A$expr)[1:3], coefficient = c(1, 1, 1),
                  stringsAsFactors = FALSE)
  pairs <- data.frame(pair_id = c("s1__A", "s1__B"), signature = "s1", cohort = c("A", "B"),
                      stringsAsFactors = FALSE)
  res <- transport_bootstrap(co, pairs, list(s1__A = g, s1__B = g), B = 300, seed = 5, progress = FALSE)
  spread <- diff(range(res$per_pair$log_hr))
  expect_true(res$ci_low < res$pooled && res$pooled < res$ci_high)
  expect_true(res$ci_high - res$ci_low < 0.5 * spread)
  expect_equal(res$n_pairs, 2)
})

test_that("transportability decision rule follows the registration", {
  d <- log(1.25)
  expect_equal(classify_transport(0.05, 0.40, d), "Transportable")
  expect_equal(classify_transport(-0.40, -0.05, d), "Contradicted")
  expect_equal(classify_transport(-0.10, 0.20, d), "Informative null")
  expect_equal(classify_transport(-0.10, 0.50, d), "Inconclusive")
})

test_that("the fast Cox fit gives the same log-HR as coxph(), with and without ties", {
  for (sd in 1:4) {
    a <- sim_cohort(seed = 20 + sd)
    sc <- score_signature(a$expr, rownames(a$expr)[1:4], c(0.5, -0.3, 0.2, 0.1))
    tm <- if (sd > 2) round(a$time) + 1 else a$time            # integer times give ties
    if (sd == 2) tm[2] <- tm[1] * (1 + 1e-12)                  # near-tie, merged by coxph's time fix
    expect_equal(log_hr_per_sd_fast(tm, a$event, sc), hr_per_sd(tm, a$event, sc)[["log_hr"]],
                 tolerance = 1e-10)
  }
})
