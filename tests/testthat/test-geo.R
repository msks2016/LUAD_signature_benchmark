# Tests for the GEO clinical fields (make_clin_geo, stage_from_tnm) and the CEL download guard.

test_that("stage groups are derived from pathological TNM", {
  t <- c("T1", "T2", "T1", "T3", "T3", "T2", "T4", "T1", "T2", "TX", "T1", "NTL", "TX")
  n <- c("N0", "N0", "N1", "N0", "N1", "N2", "N0", "N3", "N0", "N0", "NX", "NTL", "N2")
  m <- c("M0", "MX", "M0", "M0", "M0", "M0", "M0", "M0", "M1", "M0", "M0", "NTL", "M0")
  expect_identical(as.character(stage_from_tnm(t, n, m)),
                   c("I", "I", "II", "II", "III", "III", "III", "III", "IV",
                     NA, NA, NA, "III"))
})

pd_toy <- data.frame(geo_accession = paste0("GSM", 1:5),
                     "time:ch1" = c("10", "20", "30", "40", "NTL"),
                     "status:ch1" = c("DEAD", "ALIVE", "n/a", "DEAD", "NTL"),
                     "stage:ch1" = c("1A", "2B", "1B", "", NA),
                     "age1:ch1" = c("60", NA, "", "70", NA),
                     "age2:ch1" = c(NA, "55", "65", "99", NA),
                     "sex:ch1" = c("M", "F", "M", "F", NA),
                     "pt:ch1" = c("T1", "T2", "T3", "T4", "NTL"),
                     "pn:ch1" = c("N0", "N1", "N0", "N0", "NTL"),
                     "pm:ch1" = c("M0", "MX", "M0", "M0", "NTL"),
                     "flag:ch1" = c("none", "exclude", "none", "exclude", NA),
                     check.names = FALSE, stringsAsFactors = FALSE)
f_toy <- list(time = list(field = "time:ch1", unit = "months"),
              event = list(field = "status:ch1", dead_values = "DEAD", alive_values = "ALIVE"),
              stage = "stage:ch1",
              stage_recode = list("1A" = "IA", "1B" = "IB", "2A" = "IIA", "2B" = "IIB"),
              age = c("age1:ch1", "age2:ch1"), sex = "sex:ch1",
              sensitivity_exclude = list(field = "flag:ch1", values = "exclude"))

test_that("GEO fields: time, event (unknown values stay missing), recoded stage, combined age", {
  cl <- make_clin_geo(pd_toy, f_toy)
  expect_equal(cl$os_months, c(10, 20, 30, 40, NA))
  expect_equal(cl$os_event, c(1L, 0L, NA, 1L, NA))
  expect_identical(as.character(cl$stage), c("I", "II", "I", NA, NA))
  expect_equal(cl$age, c(60, 55, 65, 70, NA))
  expect_identical(cl$sensitivity_exclude, c(FALSE, TRUE, FALSE, TRUE, FALSE))
})

test_that("GEO stage can come from TNM fields", {
  f <- f_toy
  f$stage <- NULL
  f$stage_recode <- NULL
  f$stage_tnm <- list(tumour = "pt:ch1", node = "pn:ch1", metastasis = "pm:ch1")
  cl <- make_clin_geo(pd_toy, f)
  expect_identical(as.character(cl$stage), c("I", "II", "II", "III", NA))
  expect_identical(cl$stage_raw[1], "T1 N0 M0")
  f$stage_tnm <- list(t = "pt:ch1", "FALSE" = "pn:ch1", m = "pm:ch1")   # what YAML made of {t, n, m}
  expect_error(make_clin_geo(pd_toy, f), "tumour, node and metastasis")
})

test_that("follow-up can be computed from two dates; numeric codes map to stage", {
  pd <- data.frame(geo_accession = paste0("GSM", 1:4),
                   "surgery:ch1" = c("2008-12-18", "2010-04-21", "2009-01-01", NA),
                   "vital:ch1" = c("2013-04-28", "2011-04-21", "2009-01-01", "2013-04-28"),
                   "dead:ch1" = c("0", "1", "n/a", "1"),
                   "stage:ch1" = c("1", "5", "7", "3"), check.names = FALSE, stringsAsFactors = FALSE)
  f <- list(time = list(start_date = "surgery:ch1", end_date = "vital:ch1"),
            event = list(field = "dead:ch1", dead_values = "1", alive_values = "0"),
            stage = "stage:ch1",
            stage_recode = list("1" = "IA", "2" = "IB", "3" = "IIA", "4" = "IIB", "5" = "IIIA",
                                "6" = "IIIB", "7" = "IV"))
  cl <- make_clin_geo(pd, f)
  expect_equal(cl$os_months, c(as.numeric(as.Date("2013-04-28") - as.Date("2008-12-18")), 365, 0, NA) / 30.4375)
  expect_equal(cl$os_event, c(0L, 1L, NA, 1L))
  expect_identical(as.character(cl$stage), c("I", "III", "IV", "II"))
})

test_that("FPKM becomes log2(TPM + 1) with TPM summing to 1e6 per sample", {
  fpkm <- matrix(c(1, 3, 0, 10, 10, 20), nrow = 3, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2")))
  e <- fpkm_to_log2tpm(fpkm)
  expect_equal(unname(colSums(2^e - 1)), c(1e6, 1e6))
  expect_equal(unname(e["g2", "s1"]), log2(0.75e6 + 1))
})

test_that("a missing GEO field stops with its name", {
  f <- f_toy
  f$sex <- "gender:ch1"
  expect_error(make_clin_geo(pd_toy, f), "gender:ch1")
})

test_that("config/cohorts.yaml has every field script 12 needs", {
  skip_if_not_installed("yaml")
  cfg <- yaml::read_yaml(file.path(find_project_root(), "config", "cohorts.yaml"))
  for (g in c("GSE31210", "GSE50081", "GSE30219")) {
    f <- cfg[[g]]
    expect_true(!is.null(f$include$field) && length(f$include$values) > 0)
    expect_true(f$time$unit %in% c("days", "months", "years"))
    expect_true(length(f$event$dead_values) > 0 && length(f$event$alive_values) > 0)
    expect_true(!is.null(f$stage) || !is.null(f$stage_tnm))
    if (!is.null(f$stage_tnm)) {
      expect_true(all(vapply(f$stage_tnm[c("tumour", "node", "metastasis")], is.character, logical(1))))
    }
  }
  expect_identical(unname(unlist(cfg$GSE50081$stage_recode)[c("1A", "2B")]), c("IA", "IIB"))
  g <- cfg$GSE81089
  expect_identical(g$include$values, "2")
  expect_true(is.character(g$time$start_date) && is.character(g$time$end_date))
  expect_identical(unname(unlist(g$stage_recode)[c("1", "7")]), c("IA", "IV"))
  expect_identical(g$event$dead_values, "1")
})

test_that("an interrupted CEL download is never kept; a finished one is reused", {
  skip_if_not(capabilities("libcurl"), "libcurl not available")
  src <- file.path(tempdir(), "geo_src")
  dest <- file.path(tempdir(), "geo_dest")
  unlink(dest, recursive = TRUE)
  dir.create(src, showWarnings = FALSE)
  writeLines("pretend archive", file.path(src, "GSE1_RAW.tar"))
  base <- paste0("file://", normalizePath(src), "/")
  p <- geo_download_raw("GSE1", dest, tries = 1, base = base)
  expect_identical(readLines(p), "pretend archive")
  expect_false(file.exists(paste0(p, ".part")))
  unlink(file.path(src, "GSE1_RAW.tar"))                 # source gone: the kept file is reused
  expect_identical(geo_download_raw("GSE1", dest, tries = 1, base = base), p)
  expect_error(geo_download_raw("GSE2", dest, tries = 1, base = base), "failed")
})

test_that("platform gene IDs are read as text without scientific notation", {
  expect_identical(annotation_ids(c(10000, NA, 1, 100130000)), c("10000", NA, "1", "100130000"))
  expect_identical(annotation_ids(c(" 7157 ", "", "1 /// 2", NA)), c("7157", NA, "1 /// 2", NA))
})

test_that("GSE72094 settings: Entrez probe map and stage codes", {
  skip_if_not_installed("yaml")
  g <- yaml::read_yaml(file.path(find_project_root(), "config", "cohorts.yaml"))$GSE72094
  expect_identical(g$probe_annotation$id_type, "entrez")
  expect_identical(unname(unlist(g$stage_recode)[c("1", "1A", "3", "4")]), c("I", "IA", "III", "IV"))
  pd <- data.frame(geo_accession = paste0("GSM", 1:4),
                   "survival_time_in_days:ch1" = c("365", "NA", "730", "100"),
                   "vital_status:ch1" = c("Dead", "Alive", "NA", "Alive"),
                   "Stage:ch1" = c("1A", "4", "NA", "3"),
                   "age_at_diagnosis:ch1" = c("70", "NA", "66", "81"),
                   "gender:ch1" = c("F", "M", "F", "M"), check.names = FALSE, stringsAsFactors = FALSE)
  cl <- make_clin_geo(pd, g)
  expect_equal(cl$os_months, c(365, NA, 730, 100) / 30.4375)
  expect_equal(cl$os_event, c(1L, 0L, NA, 0L))
  expect_identical(as.character(cl$stage), c("I", "IV", NA, "III"))
  expect_equal(cl$age, c(70, NA, 66, 81))
})
