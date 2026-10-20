# Tests for the registered secondary metrics and controls (R/secondary.R, R/added_value.R).

harrell_uncensored <- function(time, score) {   # brute force, no censoring
  o <- outer(time, time, "<")                   # i dies before j
  s <- outer(score, score, ">") + 0.5 * outer(score, score, "==")
  sum(s[o]) / sum(o)
}

test_that("the censoring Kaplan-Meier gives G(u) and its left limit", {
  time <- c(1, 2, 2, 3, 4, 5)
  event <- c(1, 0, 1, 0, 1, 0)
  G <- censoring_km(time, event)
  expect_equal(G(c(0.5, 2, 3, 4.5, 5)), c(1, 0.8, 0.8 * 2 / 3, 0.8 * 2 / 3, 0))
  expect_equal(G(c(2, 5), left = TRUE), c(1, 0.8 * 2 / 3))
})

test_that("time-dependent AUC equals the case-control AUC without censoring", {
  time <- 1:10
  event <- rep(1, 10)
  score <- c(10, 9, 8, 3, 7, 6, 5, 4, 2, 1)
  expect_equal(td_auc(time, event, score, 4), 20 / 24)
  score_tie <- c(10, 9, 8, 7, 7, 6, 5, 4, 2, 1)       # one case-control tie counts one half
  expect_equal(td_auc(time, event, score_tie, 4), 23.5 / 24)
  expect_true(is.na(td_auc(time, event, score, 0.5)))   # no case yet
  expect_true(is.na(td_auc(time, event, score, 10)))    # no control left
})

test_that("with independent censoring the IPCW AUC recovers the uncensored AUC", {
  set.seed(11)
  n <- 4000
  x <- stats::rnorm(n)
  t_true <- stats::rexp(n, rate = 0.05 * exp(0.7 * x))
  c_time <- stats::runif(n, 0, 40)
  time <- pmin(t_true, c_time)
  event <- as.integer(t_true <= c_time)
  truth <- td_auc(t_true, rep(1, n), x, 12)
  est <- td_auc(time, event, x, 12)
  expect_lt(abs(est - truth), 0.02)
  naive_cases_only <- td_auc(time, rep(1, n), x, 12)    # ignoring censoring is biased
  expect_gt(abs(naive_cases_only - truth), abs(est - truth))
})

test_that("positive control: registered rule, monotonicity, and the D1 adaptation", {
  hr <- c(II = 1.6, III = 2.8, IV = 4.0)
  p <- c(II = 0.04, III = 0.001, IV = 0.001)
  r <- positive_control_decision(hr, p, n_iii_iv = 60)
  expect_true(r$pass)
  expect_identical(r$key_level, "III")
  r2 <- positive_control_decision(c(II = 1.6, III = 2.8, IV = 2.0), p, n_iii_iv = 60)
  expect_false(r2$monotonic)
  expect_false(r2$pass)
  r3 <- positive_control_decision(c(II = 1.9, III = 0.2), c(II = 0.01, III = 0.9), n_iii_iv = 1)
  expect_identical(r3$key_level, "II")
  expect_true(r3$pass)                                   # stage III ignored below 10 patients
  r4 <- positive_control_decision(c(II = 1.4), c(II = 0.2), n_iii_iv = 0)
  expect_false(r4$pass)
  r5 <- positive_control_decision(c(II = 1.4), c(II = 0.01), n_iii_iv = 12)
  expect_false(r5$pass)
  expect_match(r5$note, "not estimable")
})

test_that("sex is harmonised with male as the reference", {
  s <- sex_factor(c("M", "Female", "MALE", "", NA, "f"))
  expect_identical(as.character(s), c("male", "female", "male", NA, NA, "female"))
  expect_identical(levels(s)[1:2], c("male", "female"))
})

test_that("Jaccard overlap and clusters (linked signatures share a cluster)", {
  sets <- list(A = c("g1", "g2", "g3"), B = c("g2", "g3", "g4"), C = "g9",
               D = c("g4", "g5", "g6", "g7", "g8"), E = c("g3", "g4", "g10"))
  J <- jaccard_matrix(sets)
  expect_equal(J["A", "B"], 0.5)
  expect_equal(J["B", "D"], 1 / 7)
  expect_equal(unname(diag(J)), rep(1, 5))
  cl <- jaccard_clusters(J, 0.3)
  expect_identical(unname(cl[c("A", "B", "E")]), rep(cl[["A"]], 3))
  expect_true(cl[["C"]] != cl[["A"]] && cl[["D"]] != cl[["A"]] && cl[["C"]] != cl[["D"]])
})

test_that("direction consistency: one observation per cluster per cohort", {
  pairs <- data.frame(signature = c("A", "B", "C", "A", "C", "D"),
                      cohort = c("X", "X", "X", "Y", "Y", "Y"),
                      log_hr = c(0.3, -0.1, 0.2, 0.1, -0.4, 0.5), stringsAsFactors = FALSE)
  clusters <- c(A = 1L, B = 1L, C = 2L, D = 3L)
  dc <- direction_consistency(pairs, clusters)
  expect_equal(nrow(dc$observations), 5)                 # (1,X) (2,X) (1,Y) (2,Y) (3,Y)
  expect_equal(dc$observations$log_hr[dc$observations$cluster == 1 & dc$observations$cohort == "X"], 0.1)
  expect_equal(dc$summary$n, c(5, 6))
  expect_equal(dc$summary$consistent, c(4, 4))
  expect_equal(dc$summary$p_two_sided[1], stats::binom.test(4, 5, 0.5)$p.value)
  expect_error(direction_consistency(pairs, clusters[1:3]), "no cluster")
})

test_that("power simulation: beta gives the target C-index and censoring the share of deaths", {
  expect_equal(c_from_beta(0), 0.5)
  b <- beta_for_c(0.65)
  expect_equal(c_from_beta(b), 0.65, tolerance = 1e-6)
  set.seed(5)
  d <- simulate_ph_cohort(1500, b, 1e-9)                # practically no censoring
  expect_lt(abs(harrell_uncensored(d$time, d$score) - 0.65), 0.02)
  mu <- censoring_rate_for(b, 0.3)
  d2 <- simulate_ph_cohort(20000, b, mu)
  expect_lt(abs(mean(d2$event) - 0.3), 0.015)
  r <- power_simulation(100, 30, 0.65, alphas = c(0.05, 0.001), B = 20, seed = 1,
                        test_fun = function(time, event, score) 0.01)
  expect_equal(r$power, c(1, 0))
  expect_equal(r$replicates, c(20, 20))
})

test_that("the random-gene pool keeps the chosen gene types plus the signature genes", {
  gt <- data.frame(gene_id = c("E1", "E2", "E3", "E4"),
                   gene_type = c("protein_coding", "lncRNA", "snRNA", "lncRNA"), stringsAsFactors = FALSE)
  pool <- null_pool(c("E1", "E2", "E3", "E5"), gt, c("protein_coding", "lncRNA"), "E3")
  expect_setequal(pool, c("E1", "E2", "E3"))
  expect_setequal(null_pool(c("E1", "E3", "E5"), gt, NULL, "E9"), c("E1", "E3", "E5", "E9"))
})

# ---- functions that need the survival package -------------------------------------------------

sim_clinical <- function(n = 500, seed = 3) {
  set.seed(seed)
  stage <- sample(c("I", "II", "III", "IV"), n, replace = TRUE, prob = c(0.5, 0.25, 0.18, 0.07))
  lp <- c(I = 0, II = 0.6, III = 1.1, IV = 1.6)[stage]
  x <- stats::rnorm(n)
  t_event <- stats::rexp(n, rate = 0.02 * exp(lp + 0.5 * x))
  t_cens <- stats::runif(n, 0, 120)
  list(time = pmin(t_event, t_cens), event = as.integer(t_event <= t_cens), stage = stage, x = x,
       sex = sample(c("male", "female"), n, replace = TRUE), age = round(stats::rnorm(n, 65, 9)))
}

test_that("Cox-model C-index equals the score's C-index when the HR is above 1, else 1 - C", {
  skip_if_not_installed("survival")
  d <- sim_clinical()
  expect_equal(cindex_cox(d$time, d$event, d$x), cindex(d$time, d$event, d$x), tolerance = 1e-8)
  expect_equal(cindex_cox(d$time, d$event, -d$x), 1 - cindex(d$time, d$event, -d$x), tolerance = 1e-8)
})

test_that("log-rank, proportional-hazards test and RMST behave on a prognostic score", {
  skip_if_not_installed("survival")
  d <- sim_clinical()
  lr <- logrank_median(d$time, d$event, d$x)
  expect_lt(lr[["p"]], 0.05)
  expect_equal(lr[["n_high"]], sum(d$x > stats::median(d$x)))
  p <- ph_test_p(d$time, d$event, d$x)
  expect_true(p >= 0 && p <= 1)
  rm <- rmst_median_split(d$time, d$event, d$x)
  expect_lt(rm[["diff"]], 0)                               # high score, shorter survival
  expect_true(rm[["lo"]] < rm[["diff"]] && rm[["diff"]] < rm[["hi"]])
})

test_that("stage positive control passes with a real stage effect; D1 rule with few stage III-IV", {
  skip_if_not_installed("survival")
  d <- sim_clinical(n = 800)
  r <- positive_control(d$time, d$event, d$stage)
  expect_true(r$pass)
  expect_match(r$rule, "registered")
  keep <- d$stage %in% c("I", "II") | seq_along(d$stage) %in% which(d$stage == "III")[1:3]
  r2 <- positive_control(d$time[keep], d$event[keep], d$stage[keep])
  expect_match(r2$rule, "adapted")
  expect_equal(r2$n_iii_iv, 3)
})

test_that("H2 uses stage, sex and age, excludes incomplete patients and bootstraps delta C", {
  skip_if_not_installed("survival")
  d <- sim_clinical()
  d$age[1:10] <- NA
  r <- h2_added_value(d$time, d$event, d$x, d$stage, sex = d$sex, age = d$age, B = 30, seed = 1)
  expect_identical(r$covariates, "stage+sex+age")
  expect_equal(r$n, 490)
  expect_lt(r$p_lrt, 0.001)
  expect_true(!is.na(r$delta_c_low) && r$delta_c_low < r$delta_c_high)
  r0 <- h2_added_value(d$time, d$event, d$x, d$stage)
  expect_identical(r0$covariates, "stage")
  expect_true(is.na(r0$delta_c_low))
})

test_that("gene-level Cox signs follow the simulated effect", {
  skip_if_not_installed("survival")
  co <- sim_cohort(n = 400, beta = 0.8, seed = 9)
  gc <- gene_cox(co$expr, rownames(co$expr)[1:3], co$time, co$event)
  expect_true(all(gc$cox_beta > 0))
  expect_identical(gc$ensembl_id, rownames(co$expr)[1:3])
})
