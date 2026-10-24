test_that("Cohen's kappa matches hand calculation and handles single-category ratings", {
  # dual-coded scaling item before discussion: 6 of 10 agree, chance agreement 0.3 -> kappa 0.4286
  p <- c("stated", "not_stated", "not_stated", "stated", "not_stated", "stated", "not_stated", "not_stated", "stated", "stated")
  s <- c("implied", "not_stated", "implied", "stated", "implied", "stated", "not_stated", "implied", "stated", "stated")
  expect_equal(cohen_kappa(p, s), 0.3 / 0.7, tolerance = 1e-12)
  expect_equal(cohen_kappa(c("a", "b", "a"), c("a", "b", "a")), 1)
  expect_true(is.na(cohen_kappa(rep("x", 5), rep("x", 5))))
  ag <- coder_agreement(data.frame(paper_id = 1:3, v = c("y", "n", "n")), data.frame(paper_id = 1:3, v = c("y", "y", "n")), "v")
  expect_equal(ag$agree, 2)
  expect_identical(ag$disagreeing_papers, "2")
})
