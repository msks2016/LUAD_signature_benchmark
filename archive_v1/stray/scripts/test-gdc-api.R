# Tests for R/gdc_api.R with small hand-made GDC answers (no internet needed).

file_hit <- function(id, cases) {
  list(id = id, file_id = id, file_name = paste0(id, ".star_gene_counts.tsv"),
       md5sum = "0123", file_size = 4000000L, cases = cases)
}
case_rec <- function(patient, samples) {
  list(case_id = paste0("uuid-", patient), submitter_id = patient,
       primary_site = "Bronchus and lung", disease_type = "Adenomas and Adenocarcinomas",
       samples = samples)
}
fake_files <- list(
  file_hit("f1", list(case_rec("P1", list(list(submitter_id = "S1a", sample_type = "Primary Tumor"))))),
  file_hit("f2", list(case_rec("P1", list(list(submitter_id = "S1b", tissue_type = "Tumor",
                                               tumor_descriptor = "Primary"))))),
  file_hit("f3", list(case_rec("P2", list(list(submitter_id = "S2", tissue_type = "Normal",
                                               tumor_descriptor = "Not Applicable"))))),
  file_hit("f4", list(case_rec("P3", list(list(submitter_id = "S3", tissue_type = "Tumor"))),
                      case_rec("P4", list(list(submitter_id = "S4", tissue_type = "Tumor"))))),
  file_hit("f5", list(case_rec("P5", list(list(submitter_id = "S5", tissue_type = "Tumor",
                                               tumor_descriptor = "Metastatic"))))),
  file_hit("f6", list(case_rec("P6", list(list(submitter_id = "S6"))))),
  file_hit("f7", list(case_rec("P7", list(list(submitter_id = "S7", tissue_type = "Tumor")))))
)

test_that("gdc_value turns missing or empty fields into NA", {
  expect_identical(gdc_value(list(a = "x"), "a"), "x")
  expect_identical(gdc_value(list(a = 12L), "a"), "12")
  expect_identical(gdc_value(list(a = TRUE), "a"), "TRUE")
  expect_identical(gdc_value(list(a = "x"), "b"), NA_character_)
  expect_identical(gdc_value(list(a = NULL), "a"), NA_character_)
  expect_identical(gdc_value(list(), "a"), NA_character_)
  expect_identical(gdc_value(NULL, "a"), NA_character_)
})

test_that("a sample without sample_type no longer stops the file table", {
  ft <- flatten_gdc_files(fake_files)
  expect_equal(nrow(ft), 8)
  expect_identical(ft$patient, c("P1", "P1", "P2", "P3", "P4", "P5", "P6", "P7"))
  expect_true(is.na(ft$sample_type[2]))
  expect_identical(ft$tissue_type[2], "Tumor")
  expect_identical(ft$n_cases, c(1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L))
})

test_that("samples are classified from sample_type or from tissue_type and tumor_descriptor", {
  cls <- classify_gdc_sample(
    sample_type = c("Primary Tumor", NA, NA, NA, NA, "Solid Tissue Normal", "Recurrent Tumor", NA),
    tissue_type = c(NA, "Tumor", "Normal", "Tumor", "Tumor", NA, NA, NA),
    tumor_descriptor = c(NA, "Primary", "Not Applicable", "Metastatic", NA, NA, NA, NA))
  expect_identical(cls, c("primary_tumour", "primary_tumour", "normal", "other_tumour",
                          "primary_tumour", "normal", "other", "other"))
})

more_files <- c(fake_files, list(
  file_hit("f8", list(case_rec("P8", list(list(submitter_id = "S8a", tissue_type = "Tumor",
                                               tumor_descriptor = "Primary"),
                                          list(submitter_id = "S8b", sample_type = "Primary Tumor"))))),
  file_hit("f9", list(case_rec("P9", list(list(submitter_id = "S9a", sample_type = "Primary Tumor"),
                                          list(submitter_id = "S9b", tissue_type = "Normal")))))))

test_that("one primary-tumour file per patient; pooled and mixed files dropped", {
  sel <- choose_rna_files(flatten_gdc_files(more_files))
  expect_identical(sel$chosen$patient, c("P1", "P7", "P8"))
  expect_identical(sel$chosen$file_id, c("f1", "f7", "f8"))   # P1: S1a before S1b
  expect_equal(unname(sel$flow), c(9, 1, 4, 4, 1, 3))
  expect_identical(names(sel$flow), c("files", "pooled_across_patients", "not_primary_tumour",
                                      "primary_tumour_files", "of_which_several_samples",
                                      "patients_with_primary_tumour_file"))
  st <- setNames(sel$files$status, sel$files$file_id)
  expect_identical(unname(st[c("f4", "f8", "f9")]),
                   c("pooled_across_patients", "primary_tumour_several_samples",
                     "not_eligible:normal+primary_tumour"))
})

test_that("field counts show which clinical fields hold values", {
  fc <- gdc_field_counts(list(list(a = 1, b = NULL), list(a = 2, c = "x"), NULL))
  expect_identical(fc$field, c("a", "c"))
  expect_identical(fc$records_with_value, c(2L, 1L))
})

fake_cases <- list(
  list(submitter_id = "C1", case_id = "u1", primary_site = "Bronchus and lung",
       disease_type = "Adenomas and Adenocarcinomas",
       demographic = list(vital_status = "Dead", days_to_death = 300L, age_at_index = 65L,
                          sex_at_birth = "male"),
       diagnoses = list(list(primary_diagnosis = "Adenocarcinoma, NOS",
                             ajcc_pathologic_stage = "Stage IB")),
       follow_ups = list(list(days_to_follow_up = 100L), list(days_to_follow_up = 200L))),
  list(submitter_id = "C2", demographic = list(vital_status = "Alive", gender = "female"),
       diagnoses = list(list(ajcc_pathologic_stage = "Stage IIIA", age_at_diagnosis = 25600L)),
       follow_ups = list(list(days_to_follow_up = 100L), list(days_to_follow_up = 400L),
                         list(timepoint_category = "Follow-up"))),
  list(submitter_id = "C3", demographic = list(vital_status = "Alive", days_to_birth = -20000L),
       diagnoses = list(list(days_to_last_follow_up = 250L, ajcc_pathologic_stage = "Stage IA1"))),
  list(submitter_id = "C4", demographic = list(gender = "female"),
       diagnoses = list(list(diagnosis_is_primary_disease = FALSE, ajcc_pathologic_stage = "Stage IV",
                             primary_diagnosis = "Squamous cell carcinoma, NOS"),
                        list(diagnosis_is_primary_disease = TRUE, ajcc_pathologic_stage = "Stage IIA",
                             primary_diagnosis = "Adenocarcinoma, NOS"))),
  list(submitter_id = "C5")
)

test_that("clinical records give follow-up, stage and sex the way TCGAbiolinks does", {
  cl <- flatten_gdc_cases(fake_cases)
  expect_identical(cl$submitter_id, c("C1", "C2", "C3", "C4", "C5"))
  expect_equal(cl$days_to_death, c(300, NA, NA, NA, NA))
  expect_equal(cl$days_to_last_follow_up, c(200, 400, 250, NA, NA))
  expect_identical(cl$follow_up_source, c("follow_ups", "follow_ups", "diagnoses", NA, NA))
  expect_identical(cl$ajcc_pathologic_stage, c("Stage IB", "Stage IIIA", "Stage IA1", "Stage IIA", NA))
  expect_identical(cl$primary_diagnosis[4], "Adenocarcinoma, NOS")
  expect_identical(cl$sex_at_birth, c("male", "female", NA, "female", NA))
  expect_identical(cl$vital_status, c("Dead", "Alive", "Alive", NA, NA))
  expect_equal(cl$n_diagnoses, c(1, 1, 1, 2, 0))
  expect_equal(cl$age_at_index, c(65, 70, 54, NA, NA))
  expect_identical(cl$age_source, c("age_at_index", "age_at_diagnosis", "days_to_birth", NA, NA))
  mc <- make_clin_gdc(cl)
  expect_equal(mc$os_months, c(300, 400, 250, NA, NA) / 30.4375)
  expect_equal(mc$os_event, c(1L, 0L, 0L, NA, NA))
  expect_identical(as.character(mc$stage), c("I", "III-IV", "I", "II", NA))
})

write_star <- function(path, tpm) {
  writeLines(c("# gene-model: GENCODE v36",
               "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\ttpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded",
               "N_unmapped\t\t\t10\t10\t10\t\t\t",
               "N_multimapping\t\t\t5\t5\t5\t\t\t",
               "N_noFeature\t\t\t7\t7\t7\t\t\t",
               "N_ambiguous\t\t\t3\t3\t3\t\t\t",
               sprintf("ENSG00000000003.15\tTSPAN6\tprotein_coding\t100\t50\t50\t%s\t1\t1", tpm[1]),
               sprintf("ENSG00000002586.20_PAR_Y\tCD99\tprotein_coding\t0\t0\t0\t%s\t0\t0", tpm[2]),
               sprintf("ENSG00000002586.20\tCD99\tprotein_coding\t500\t250\t250\t%s\t2\t2", tpm[3])),
             path)
}

test_that("STAR-Counts files are read into one TPM matrix", {
  d <- file.path(tempdir(), "star_test")
  dir.create(d, showWarnings = FALSE)
  p <- file.path(d, c("a.tsv", "b.tsv"))
  write_star(p[1], c(3, 0, 63))
  write_star(p[2], c(7, 0, 1))
  m <- read_star_tpm(p, c("P1", "P2"))
  expect_identical(dim(m), c(3L, 2L))
  expect_identical(rownames(m), c("ENSG00000000003.15", "ENSG00000002586.20_PAR_Y", "ENSG00000002586.20"))
  expect_equal(unname(m[, "P1"]), c(3, 0, 63))
  expect_equal(unname(m[, "P2"]), c(7, 0, 1))
  expect_identical(attr(m, "gene_model"), "# gene-model: GENCODE v36")
})

tpm_toy <- matrix(c(3, 0, 63, 7, 0, 1, 1, 2, 3), nrow = 3,
                  dimnames = list(c("ENSG00000000003.15", "ENSG00000002586.20_PAR_Y", "ENSG00000002586.20"),
                                  c("C1", "C2", "C3")))
clin_toy <- data.frame(submitter_id = c("C1", "C2", "C3"), vital_status = c("Dead", "Alive", "Alive"),
                       days_to_death = c(300, NA, NA), days_to_last_follow_up = c(200, 400, 250),
                       ajcc_pathologic_stage = c("Stage IB", "Stage IIIA", "Stage IA1"),
                       age_at_index = c(65, 70, 55), sex_at_birth = c("male", "female", "male"),
                       stringsAsFactors = FALSE)

test_that("the CPTAC-3 builder transforms expression like TCGA (log2 TPM + 1, no PAR_Y, no version)", {
  co <- build_star_cohort(tpm_toy, clin_toy, name = "toy")
  expect_identical(rownames(co$expr), c("ENSG00000000003", "ENSG00000002586"))
  expect_equal(unname(co$expr[, "C1"]), log2(c(3, 63) + 1))
  expect_identical(colnames(co$expr), co$clin$patient_id)
  expect_equal(co$event, c(1L, 0L, 0L))
})

test_that("build_star_cohort() and build_gdc_cohort() give identical cohorts from the same data", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")
  tpm_bc <- tpm_toy
  colnames(tpm_bc) <- c("B1", "B2", "B3")
  cd <- S4Vectors::DataFrame(barcode = c("B1", "B2", "B3"), patient = c("C1", "C2", "C3"),
                             sample_type = rep("Primary Tumor", 3), row.names = c("B1", "B2", "B3"))
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(tpm_unstrand = tpm_bc), colData = cd)
  a <- build_gdc_cohort(se, clin_toy, name = "toy")
  b <- build_star_cohort(tpm_toy, clin_toy, name = "toy")
  expect_identical(a$expr, b$expr)
  expect_identical(a$clin, b$clin)
  expect_identical(a$time, b$time)
  expect_identical(a$event, b$event)
})

test_that("GDC filters keep single values as JSON arrays", {
  skip_if_not_installed("jsonlite")
  f <- gdc_filter(list("cases.project.project_id" = "CPTAC-3", "files.data_type" = c("A", "B")))
  expect_identical(
    as.character(jsonlite::toJSON(f, auto_unbox = TRUE)),
    paste0('{"op":"and","content":[',
           '{"op":"in","content":{"field":"cases.project.project_id","value":["CPTAC-3"]}},',
           '{"op":"in","content":{"field":"files.data_type","value":["A","B"]}}]}'))
})

test_that("downloads are checked by md5 and finished files are not fetched again", {
  skip_if_not(capabilities("libcurl"), "libcurl not available")
  src <- file.path(tempdir(), "gdc_src")
  dest <- file.path(tempdir(), "gdc_dest")
  unlink(dest, recursive = TRUE)
  dir.create(src, showWarnings = FALSE)
  writeLines("some content", file.path(src, "id1"))
  writeLines("other content", file.path(src, "id2"))
  md5 <- unname(tools::md5sum(file.path(src, c("id1", "id2"))))
  base <- paste0("file://", normalizePath(src), "/")
  p <- gdc_download(c("id1", "id2"), c("one.tsv", "two.tsv"), md5, dest, tries = 1, base = base)
  expect_true(all(file.exists(p)))
  expect_identical(readLines(p[2]), "other content")
  expect_message(gdc_download(c("id1", "id2"), c("one.tsv", "two.tsv"), md5, dest, tries = 1, base = base),
                 "2 of 2 files already downloaded")
  expect_error(gdc_download("id1", "bad.tsv", "ffff", dest, tries = 1, base = base),
               "Could not download")
})
