test_that("raw-scale score equals a hand calculation", {
  expr <- matrix(c(1, 2, 3,
                   4, 6, 8), nrow = 2, byrow = TRUE,
                 dimnames = list(c("g1", "g2"), c("a", "b", "c")))
  s <- score_signature(expr, c("g1", "g2"), c(0.5, -1), scale = "raw")
  expect_equal(s, c(0.5 * 1 - 4, 0.5 * 2 - 6, 0.5 * 3 - 8))
})

test_that("z-scored score equals a hand calculation", {
  expr <- matrix(c(1, 2, 3,
                   4, 6, 8), nrow = 2, byrow = TRUE,
                 dimnames = list(c("g1", "g2"), c("a", "b", "c")))
  # g1 -> z = (-1, 0, 1); g2 -> z = (-1, 0, 1)
  s <- score_signature(expr, c("g1", "g2"), c(0.5, -1), scale = "zscore")
  expect_equal(s, c(0.5, 0, -0.5))
})

test_that("z-scoring removes each gene's location and scale", {
  set.seed(1)
  expr <- matrix(rnorm(30), nrow = 3, dimnames = list(c("g1", "g2", "g3"), NULL))
  shifted <- expr
  shifted["g2", ] <- 10 + 3 * shifted["g2", ]
  expect_equal(score_signature(expr, c("g1", "g2"), c(1, 2)),
               score_signature(shifted, c("g1", "g2"), c(1, 2)))
})

test_that("absent or constant genes stop scoring", {
  expr <- matrix(c(1, 2, 3,
                   5, 5, 5), nrow = 2, byrow = TRUE,
                 dimnames = list(c("g1", "flat"), NULL))
  expect_error(score_signature(expr, c("g1", "gX"), c(1, 1)), "Not in matrix")
  expect_error(score_signature(expr, c("g1", "flat"), c(1, 1)), "Constant or missing")
})

test_that("duplicated genes are refused", {
  expr <- matrix(rnorm(6), nrow = 2, dimnames = list(c("g1", "g2"), NULL))
  expect_error(score_signature(expr, c("g1", "g1"), c(1, 1)))
})
