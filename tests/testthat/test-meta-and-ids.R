test_that("k = 2 pooling is the inverse-variance mean", {
  y <- log(c(1.516, 1.687))
  se <- c(0.1432, 0.1563)
  r <- pool_log_hr(y, se)
  w <- 1 / se^2
  expect_equal(r$log_hr, sum(w * y) / sum(w), tolerance = 1e-8)
  expect_equal(r$k, 2)
})

test_that("k >= 3 uses REML with Hartung-Knapp (values from the 16 October audit)", {
  y <- log(c(1.516, 1.687, 1.348, 1.278))
  se <- c(0.1432, 0.1563, 0.1371, 0.1513)
  r <- pool_log_hr(y, se, c("CPTAC-3", "GSE31210", "GSE50081", "GSE30219"))
  expect_true(abs(r$hr - 1.442) < 0.003)
  expect_true(abs(r$hr_low - 1.192) < 0.01)
  expect_true(abs(r$hr_high - 1.745) < 0.01)
  loo <- pool_leave_one_out(y, se, c("CPTAC-3", "GSE31210", "GSE50081", "GSE30219"))
  expect_equal(nrow(loo), 4)
})

test_that("Ensembl version suffixes are stripped; PAR_Y ids stay distinct", {
  expect_equal(strip_version(c("ENSG00000228340.6", "ENSG00000233621")),
               c("ENSG00000228340", "ENSG00000233621"))
  expect_equal(strip_version("ENSG00000182378.14_PAR_Y"), "ENSG00000182378.14_PAR_Y")
})

test_that("symbol resolution flags conflicts and misses; the paper's annotation decides", {
  v22 <- data.frame(gene_id = c("E1", "E2", "E3"), gene_name = c("A", "B", "C"),
                    stringsAsFactors = FALSE)
  v36 <- data.frame(gene_id = c("E1", "E9", "E3", "E4"), gene_name = c("A", "B", "C2", "D"),
                    stringsAsFactors = FALSE)
  r <- resolve_symbols(c("A", "B", "C", "Z"), list(v22 = v22, v36 = v36))
  expect_equal(r$status, c("consistent", "conflict", "consistent", "not_found"))
  sig <- data.frame(signature = "s1", published_symbol = c("A", "B", "C", "Z"),
                    coefficient = 1:4, annotation = "v22", stringsAsFactors = FALSE)
  ov <- data.frame(signature = "s1", published_symbol = "Z", ensembl_id = "E7",
                   stringsAsFactors = FALSE)
  m <- choose_ids(sig, r, ov)
  expect_equal(m$ensembl_id, c("E1", "E2", "E3", "E7"))
  expect_false(any(m$needs_review))
})

test_that("probe collapse drops multi-gene probes and prefers _at probes", {
  expr <- rbind(p1_at = c(5, 6, 7), p1_s_at = c(9, 9, 10), p2_at = c(3, 4, 5), px_at = c(8, 8, 9))
  pm <- data.frame(probe_id = c("p1_at", "p1_s_at", "p2_at", "px_at", "px_at"),
                   ensembl_id = c("G1", "G1", "G2", "G2", "G3"), stringsAsFactors = FALSE)
  out <- collapse_probes(expr, pm)
  expect_equal(sort(rownames(out$expr)), c("G1", "G2"))
  expect_equal(out$probe_choice$probe_id[out$probe_choice$ensembl_id == "G1"], "p1_at")
})

test_that("interval variants: modified Hartung-Knapp is never narrower than common effect", {
  # Yao 2021, five independent cohorts (v2 results, 24 Oct 2026): q < 1, tau2 = 0
  y <- c(0.39209, 0.38873, 0.35901, 0.29865, 0.24645)
  se <- c(0.12601, 0.15655, 0.14897, 0.15467, 0.14912)
  v <- pool_ci_variants(y, se, tau2 = 0)
  hk <- v[v$interval == "reml_hk", ]
  ce <- v[v$interval == "common_effect", ]
  mhk <- v[v$interval == "reml_hk_modified", ]
  expect_true(hk$q < 1)
  expect_lt(hk$ci_high - hk$ci_low, ce$ci_high - ce$ci_low)      # the registered HK is narrower here
  expect_gte(mhk$ci_high - mhk$ci_low, ce$ci_high - ce$ci_low)
  expect_equal(ce$log_hr, hk$log_hr)                              # tau2 = 0: same centre
  expect_true(abs(mhk$hr_low - 1.173) < 0.003 && abs(mhk$hr_high - 1.685) < 0.003)
})
