# 13_cohort_gse81089.R - GSE81089 (Uppsala NSCLC RNA-seq; Djureinovic et al., JCI Insight 2016):
# adenocarcinomas only. Library: Illumina TruSeq RNA v2 (poly(A) selection), HiSeq 2500,
# 2 x 100 bp; TopHat 2.0.8b on GRCh37 with Ensembl 73 genes; FPKM from Cufflinks 2.1.1.
# Expression: FPKM -> TPM (FPKM / sample total x 1e6) -> log2(TPM + 1), Ensembl gene IDs.
# FPKM coded -1 (not quantified) counts as missing; genes with a missing FPKM in any selected
# sample are dropped (counted in the summary and listed in data/processed/).
# Building the cohort is descriptive only. It enters the analyses (cohorts$independent in
# config/analysis.yaml) only after the OSF amendment restores it.
# Output to check: results/qc/GSE81089_build_summary.txt

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
f81 <- yaml::read_yaml("config/cohorts.yaml")$GSE81089
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

summary_file <- "results/qc/GSE81089_build_summary.txt"
report <- c("GSE81089 build summary", paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")))
add <- function(title, x) {
  report <<- c(report, "", paste("==", title), utils::capture.output(print(x)))
  writeLines(report, summary_file)
}
count_table <- function(...) {
  d <- as.data.frame(table(..., useNA = "ifany"), responseName = "n")
  d[d$n > 0, , drop = FALSE]
}

# 1. Clinical data: adenocarcinomas (histology code 2)
pd_all <- geo_pheno("GSE81089")
add("Histology code (1 squamous, 2 adenocarcinoma, 3 large cell/NOS) of all samples",
    table(pd_all[[f81$include$field]], useNA = "ifany"))
pd <- pd_all[tolower(trimws(pd_all[[f81$include$field]])) %in% tolower(f81$include$values), , drop = FALSE]
clin <- make_clin_geo(pd, f81)
add("Event (1 = dead) by raw value", count_table(raw = pd[[f81$event$field]], event = clin$os_event))
add("Stage: raw code -> group", count_table(raw = clin$stage_raw, group = clin$stage))
add("Follow-up in months (summary)", summary(clin$os_months))
add("Latest contact dates (most common)", utils::head(sort(table(pd[[f81$time$end_date]]), decreasing = TRUE), 5))
add("Age recorded / sex", list(age_recorded = sum(!is.na(clin$age)), sex = table(clin$sex, useNA = "ifany")))

# 2. Expression: FPKM file -> log2(TPM + 1), columns renamed to GSM IDs
path <- geo_download_supp("GSE81089", "GSE81089_FPKM_cufflinks.tsv.gz")
raw <- utils::read.delim(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE)
fpkm <- as.matrix(raw[, -1, drop = FALSE])
rownames(fpkm) <- raw[[1]]
# GEO note: columns L608T_2122 and L771T_1 belong to samples L608T and L771T
colnames(fpkm)[colnames(fpkm) == "L608T_2122"] <- "L608T"
colnames(fpkm)[colnames(fpkm) == "L771T_1"] <- "L771T"
if (anyDuplicated(colnames(fpkm))) stop("Duplicated sample columns in the FPKM file")
if (anyDuplicated(rownames(fpkm))) stop("Duplicated gene IDs in the FPKM file")
add("FPKM file: genes x samples", dim(fpkm))
if (!is.numeric(fpkm)) stop("The FPKM file has non-numeric values")
sample_names <- trimws(pd[[f81$sample_name]])
absent <- setdiff(sample_names, colnames(fpkm))
if (length(absent)) stop("No expression column for: ", paste(absent, collapse = ", "))
sel <- fpkm[, sample_names, drop = FALSE]
# The file codes values that could not be quantified as -1 (every negative value is exactly -1):
# treat them as missing. Any other negative value would be an error.
minus_one <- !is.na(sel) & sel == -1
add("FPKM values coded -1 (not quantified) in the adenocarcinomas",
    c(cells = sum(minus_one), genes_with_any = sum(rowSums(minus_one) > 0)))
sel[minus_one] <- NA
if (any(sel < 0, na.rm = TRUE)) stop("Negative FPKM values other than -1: send ", summary_file)
na_by_sample <- colSums(is.na(sel))
add("Missing FPKM values in the adenocarcinomas",
    c(cells = sum(na_by_sample), genes_with_any = sum(rowSums(is.na(sel)) > 0),
      samples_with_any = sum(na_by_sample > 0)))
add("Samples with the most missing FPKM values", utils::head(sort(na_by_sample, decreasing = TRUE), 5))
add("Missing values by gene ID prefix (first 4 characters)",
    table(substr(rownames(sel)[rowSums(is.na(sel)) > 0], 1, 4)))
if (any(na_by_sample > 0.05 * nrow(sel))) {
  stop("More than 5% missing FPKM in some samples: send ", summary_file)
}
# A gene with a missing value in any adenocarcinoma is dropped (it could not count as present
# under the presence rule anyway); TPM is then computed over the remaining genes.
complete <- rowSums(is.na(sel)) == 0
utils::write.csv(data.frame(ensembl_id = rownames(sel)[!complete],
                            n_missing = rowSums(is.na(sel))[!complete]),
                 "data/processed/GSE81089_genes_dropped_missing.csv", row.names = FALSE)
add("Genes dropped for a missing or -1 value in any adenocarcinoma (listed in data/processed/GSE81089_genes_dropped_missing.csv)",
    c(dropped = sum(!complete), kept = sum(complete)))
expr <- fpkm_to_log2tpm(sel[complete, , drop = FALSE])
colnames(expr) <- pd$geo_accession
expr <- expr[grepl("^ENSG", rownames(expr)), , drop = FALSE]
rownames(expr) <- strip_version(rownames(expr))
if (anyDuplicated(rownames(expr))) stop("Duplicated Ensembl IDs after removing versions")

# 3. Cohort object
meta <- list(source = "GEO GSE81089_FPKM_cufflinks.tsv.gz; TPM from FPKM; log2(TPM + 1)",
             gene_model = "Ensembl 73 (GRCh37)", library = "TruSeq RNA v2, poly(A) selection",
             downloaded = format(Sys.Date()))
co <- finalise_cohort(expr, clin, "GSE81089", platform = "RNA-seq", meta = meta)
check_cohort(co)
saveRDS(co, "data/processed/GSE81089.rds")

add("Sample flow (patients)", co$meta$sample_flow)
add("Cohort", c(patients = ncol(co$expr), deaths = sum(co$event),
                median_follow_up_months = round(reverse_km_median(co$time, co$event), 1),
                max_follow_up_months = round(max(co$time), 1), genes = nrow(co$expr)))
add("Stage group", table(co$clin$stage, useNA = "ifany"))
add("Expression range, log2(TPM + 1)", round(range(co$expr), 2))
cat(report, sep = "\n")
message("\nGSE81089: ", ncol(co$expr), " patients, ", sum(co$event), " deaths. Summary saved to ",
        summary_file)
