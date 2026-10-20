# Tests for the TPM-composition check and the rescaling sensitivity analysis.

toy_rnaseq <- function() {
  tpm <- matrix(c(300, 100, 600,     # sample 1: 60% of TPM in the "other" RNA
                  450, 450, 100),    # sample 2: 10%
                nrow = 3, dimnames = list(c("G_pc", "G_lnc", "G_other"), c("S1", "S2")))
  list(name = "toy", expr = log2(tpm + 1), meta = list(platform = "RNA-seq"))
}

test_that("the share of TPM outside the kept genes is computed per sample", {
  co <- toy_rnaseq()
  expect_equal(unname(tpm_share_outside(co$expr, c("G_pc", "G_lnc"))), c(0.6, 0.1))
})

test_that("rescaling makes the kept genes sum to 1e6 and keeps every gene", {
  co <- rescale_tpm_cohort(toy_rnaseq(), c("G_pc", "G_lnc"))
  tpm <- 2^co$expr - 1
  expect_equal(unname(colSums(tpm[c("G_pc", "G_lnc"), ])), c(1e6, 1e6))
  expect_equal(unname(tpm["G_other", "S1"]), 600 * 1e6 / 400)
  expect_identical(rownames(co$expr), c("G_pc", "G_lnc", "G_other"))
})

test_that("microarray cohorts are left unchanged", {
  co <- toy_rnaseq()
  co$meta$platform <- "microarray"
  expect_identical(rescale_tpm_cohort(co, "G_pc"), co)
})
