test_that("C-index treats a higher score as higher risk", {
  time <- c(1, 2, 3, 4)
  event <- c(1, 1, 1, 1)
  score <- c(4, 3, 2, 1)
  expect_equal(cindex(time, event, score), 1)
  expect_equal(cindex(time, event, -score), 0)
})

test_that("HR per SD recovers a simulated effect", {
  co <- sim_cohort(n = 2000, beta = 0.5, seed = 21)
  s <- score_signature(co$expr, rownames(co$expr)[1:3], c(1, 1, 1))
  est <- hr_per_sd(co$time, co$event, s)
  expect_true(abs(est[["log_hr"]] - 0.5) < 0.1)
  expect_true(est[["lo"]] < est[["hr"]] && est[["hr"]] < est[["hi"]])
})

test_that("a constant score is refused", {
  expect_error(hr_per_sd(c(1, 2, 3), c(1, 0, 1), c(2, 2, 2)), "Constant")
})

test_that("stage text collapses to I / II / III-IV", {
  x <- c("Stage IA1", "IB", "Stage IIB", "stage IIIA", "Stage IV", "IV", "Not Reported", NA)
  expect_equal(as.character(collapse_stage(x)),
               c("I", "I", "II", "III-IV", "III-IV", "III-IV", NA, NA))
})

test_that("H2 returns a likelihood-ratio test and an adjusted HR", {
  co <- sim_cohort(n = 400, beta = 0.6, seed = 31)
  set.seed(32)
  stage <- sample(c("I", "II", "III-IV"), 400, replace = TRUE)
  s <- score_signature(co$expr, rownames(co$expr)[1:3], c(1, 1, 1))
  r <- h2_added_value(co$time, co$event, s, stage)
  expect_equal(r$n, 400)
  expect_true(r$p_lrt < 0.001)
  expect_true(r$hr_adj > 1)
})
