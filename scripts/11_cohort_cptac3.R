# 11_cohort_cptac3.R - CPTAC-3 lung adenocarcinoma: GDC STAR-Counts RNA-seq (tpm_unstranded) and
# GDC clinical data, both read straight from the GDC API (functions in R/gdc_api.R).
# TCGAbiolinks::GDCquery() is not used here: on CPTAC-3 it stops with
# "object 'sample_type' not found", because some CPTAC-3 samples are now described only by
# tissue_type and tumor_descriptor. Expression is transformed exactly as for TCGA-LUAD
# (log2(TPM + 1), Ensembl IDs without version, PAR_Y rows dropped); a unit test checks that
# build_star_cohort() and build_gdc_cohort() give identical results from the same data.
# If the script stops (network), run it again: files already downloaded are kept.
# Output to check: results/qc/CPTAC-3_build_summary.txt

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
raw_dir <- "data/raw/GDC/CPTAC-3"
dir.create(file.path(raw_dir, "STAR-Counts"), recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

summary_file <- "results/qc/CPTAC-3_build_summary.txt"
release <- gdc_release()
report <- c("CPTAC-3 build summary",
            paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
            paste("GDC:", release))
add <- function(title, x) {       # add one block to the summary file (rewritten each time)
  report <<- c(report, "", paste("==", title), utils::capture.output(print(x)))
  writeLines(report, summary_file)
}

# 1. Clinical records of every CPTAC-3 case (several cancer types) -> lung adenocarcinoma ----
message("Reading CPTAC-3 clinical records from the GDC ...")
case_hits <- gdc_hits("cases", gdc_filter(list("cases.project.project_id" = "CPTAC-3")),
                      expand = c("demographic", "diagnoses", "follow_ups"))
saveRDS(case_hits, file.path(raw_dir, "cases_api.rds"))
clin_all <- flatten_gdc_cases(case_hits)
lung <- clin_all$primary_site %in% "Bronchus and lung"
clin_raw <- clin_all[lung & clin_all$disease_type %in% "Adenomas and Adenocarcinomas", , drop = FALSE]
utils::write.csv(clin_raw, file.path(raw_dir, "clinical_luad.csv"), row.names = FALSE)

add("CPTAC-3 cases: all / lung (any type) / lung adenocarcinoma",
    c(all = nrow(clin_all), lung = sum(lung), lung_adenocarcinoma = nrow(clin_raw)))
add("Primary site of all CPTAC-3 cases", table(clin_all$primary_site, useNA = "ifany"))
add("Disease type of the lung cases", table(clin_all$disease_type[lung], useNA = "ifany"))
add("Primary diagnosis (lung adenocarcinoma cases)", table(clin_raw$primary_diagnosis, useNA = "ifany"))
add("Diagnoses per case", table(clin_raw$n_diagnoses, useNA = "ifany"))
add("Vital status", table(clin_raw$vital_status, useNA = "ifany"))
add("Where the last follow-up came from", table(clin_raw$follow_up_source, useNA = "ifany"))
add("Dead but no days_to_death", sum(clin_raw$vital_status %in% "Dead" & is.na(clin_raw$days_to_death)))
add("Where the age came from", table(clin_raw$age_source, useNA = "ifany"))
luad_hits <- case_hits[vapply(case_hits, gdc_value, character(1), name = "submitter_id") %in%
                         clin_raw$submitter_id]
add("Fields with a value, lung adenocarcinoma cases: demographic",
    gdc_field_counts(lapply(luad_hits, function(h) h$demographic)))
add("Fields with a value, lung adenocarcinoma cases: diagnoses",
    gdc_field_counts(unlist(lapply(luad_hits, function(h) h$diagnoses), recursive = FALSE)))
add("Fields with a value, lung adenocarcinoma cases: follow_ups",
    gdc_field_counts(unlist(lapply(luad_hits, function(h) h$follow_ups), recursive = FALSE)))
if (!nrow(clin_raw)) stop("No lung adenocarcinoma case found: send ", summary_file)

# 2. CPTAC-3 RNA-seq files (STAR - Counts) -> one primary-tumour file per patient ------------
message("Reading the CPTAC-3 RNA-seq file list from the GDC ...")
file_hits <- gdc_hits("files", gdc_filter(list("cases.project.project_id" = "CPTAC-3",
                                               "files.data_type" = "Gene Expression Quantification",
                                               "files.analysis.workflow_type" = "STAR - Counts",
                                               "files.access" = "open")),
                      expand = c("cases", "cases.samples"))
saveRDS(file_hits, file.path(raw_dir, "files_api.rds"))
ft_all <- flatten_gdc_files(file_hits)
sel <- choose_rna_files(ft_all[ft_all$patient %in% clin_raw$submitter_id, , drop = FALSE])
utils::write.csv(sel$all, file.path(raw_dir, "rnaseq_files_luad.csv"), row.names = FALSE)

by_file_record <- ft_all$primary_site %in% "Bronchus and lung" &
  ft_all$disease_type %in% "Adenomas and Adenocarcinomas"
add("STAR-Counts files: all CPTAC-3 / lung adenocarcinoma (clinical record) / (file record)",
    c(all = length(unique(ft_all$file_id)),
      luad_by_clinical_record = length(unique(sel$all$file_id)),
      luad_by_file_record = length(unique(ft_all$file_id[by_file_record]))))
st <- sel$all
types <- as.data.frame(table(sample_type = st$sample_type, tissue_type = st$tissue_type,
                             tumor_descriptor = st$tumor_descriptor, class = st$class,
                             useNA = "ifany"), responseName = "rows")
add("Sample fields of the lung adenocarcinoma files", types[types$rows > 0, , drop = FALSE])
add("Preservation method", table(st$preservation_method, useNA = "ifany"))
add("Files by status (one row per file)", table(sel$files$status))
add("File selection", sel$flow)
utils::write.csv(sel$files, file.path(raw_dir, "rnaseq_files_status.csv"), row.names = FALSE)
if (!nrow(sel$chosen)) stop("No primary-tumour RNA-seq file found: send ", summary_file)

# 3. Download the chosen files and read TPM --------------------------------------------------
message("Downloading ", nrow(sel$chosen), " STAR-Counts files (one per patient) ...")
paths <- gdc_download(sel$chosen$file_id, sel$chosen$file_name, sel$chosen$md5sum,
                      dir = file.path(raw_dir, "STAR-Counts"))
utils::write.csv(sel$chosen, file.path(raw_dir, "rnaseq_files_chosen.csv"), row.names = FALSE)
message("Reading the files ...")
tpm <- read_star_tpm(paths, sel$chosen$patient)

# 4. Cohort object -------------------------------------------------------------------------
meta <- list(downloaded = format(Sys.Date()), gdc_data_release = release,
             cases_cptac3 = nrow(clin_all), cases_luad = nrow(clin_raw), file_selection = sel$flow,
             patient_rule = paste("one primary-tumour file per patient (one case; every sample a primary",
                                  "tumour): first by sample ID, then file ID"))
co <- build_star_cohort(tpm, clin_raw, name = "CPTAC-3", meta = meta)
check_cohort(co)
saveRDS(co, "data/processed/CPTAC-3.rds")

add("Gene model", co$meta$gene_model)
add("Sample flow (patients)", co$meta$sample_flow)
add("Cohort", c(patients = ncol(co$expr), deaths = sum(co$event),
                median_follow_up_months = round(reverse_km_median(co$time, co$event), 1),
                max_follow_up_months = round(max(co$time), 1), genes = nrow(co$expr)))
add("Stage", table(co$clin$stage, useNA = "ifany"))
add("Sex", table(co$clin$sex, useNA = "ifany"))
add("Age recorded", sum(!is.na(co$clin$age)))
add("Expression range, log2(TPM + 1)", round(range(co$expr), 2))
gt <- read_star_gene_table(paths[1])
top <- names(sort(apply(co$expr, 1, max), decreasing = TRUE))[1:5]
add("Five genes with the highest single values, log2(TPM + 1)",
    data.frame(gene_id = top, gene_name = gt$gene_name[match(top, gt$gene_id)],
               gene_type = gt$gene_type[match(top, gt$gene_id)],
               median = round(apply(co$expr[top, , drop = FALSE], 1, stats::median), 2),
               max = round(apply(co$expr[top, , drop = FALSE], 1, max), 2), row.names = NULL))
cat(report, sep = "\n")
message("\nCPTAC-3: ", ncol(co$expr), " patients, ", sum(co$event), " deaths. Summary saved to ",
        summary_file)
