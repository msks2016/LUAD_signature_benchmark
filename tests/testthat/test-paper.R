test_that("manuscript formatting helpers", {
  expect_identical(sig_label(c("Paper62_Yao2021", "Paper9_Zhang2021")), c("Yao 2021", "Zhang 2021"))
  expect_identical(fmt_p(c(0.0004, 0.0056, 0.056, 0.6, 1, NA)), c("<0.001", "0.006", "0.056", "0.60", "1.00", ""))
  expect_identical(p_text(c(0.00002, 0.03)), c("p < 0.001", "p = 0.030"))
  expect_identical(fmt_ci(1.234, 1.1, 1.456), "1.23 (1.10-1.46)")
  expect_identical(fmt_ci(0.018, -0.003, 0.061, 3, " to "), "0.018 (-0.003 to 0.061)")
  expect_identical(fmt_ci(NA, NA, NA), "")
})

test_that("Table 2 age and sex summaries", {
  expect_identical(age_summary(c(50, 60, 70, NA, "unknown")), "60 (55-65)")
  expect_identical(age_summary(c(NA, NA)), "")
  expect_identical(male_summary(c("male", "Female", "M", NA, "not reported")), "2 (67%)")
  expect_identical(male_summary(c(NA, "")), "")
})
