# 14_cohort_gse72094.R - GSE72094 (Moffitt lung adenocarcinomas; Rosetta/Merck custom Affymetrix
# array GPL15048). Independent cohort in OSF Update 3, analysed for any signature that reaches
# 80% coverage. Expression: the submitter's normalised values from the series matrix (log2 if
# not already on a log scale); probes mapped to Ensembl genes, one probe per gene.
#
# Two runs:
#  1. Without a GSE72094 entry in config/cohorts.yaml the script only writes
#     logs/gse72094_fields.txt (clinical fields with their values, and the platform annotation
#     columns) and stops. Send that file; the config entry is filled from it.
#  2. With the entry, it builds data/processed/GSE72094.rds.
# Output to check: results/qc/GSE72094_build_summary.txt

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
gse <- "GSE72094"
gpl_id <- "GPL15048"
dir.create("logs", showWarnings = FALSE)
dir.create("data/raw/GEO", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

g <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = "data/raw/GEO")[[1]]
pd <- Biobase::pData(g)
gpl <- GEOquery::getGEO(gpl_id, destdir = "data/raw/GEO")
annot <- GEOquery::Table(gpl)
f72 <- yaml::read_yaml("config/cohorts.yaml")[[gse]]

if (is.null(f72)) {
  out <- c(paste(gse, "fields -", format(Sys.time(), "%Y-%m-%d %H:%M")),
           paste("Samples:", nrow(pd), " Probes in the series matrix:", nrow(g)),
           paste("Expression range:", paste(round(range(Biobase::exprs(g), na.rm = TRUE), 2), collapse = " to ")),
           paste("Missing expression values:", sum(is.na(Biobase::exprs(g)))))
  for (fld in grep(":ch1$", names(pd), value = TRUE)) {
    v <- sort(table(pd[[fld]], useNA = "ifany"), decreasing = TRUE)
    out <- c(out, "", paste0(fld, "  (", length(v), " distinct values)"),
             utils::capture.output(print(utils::head(v, 10))))
  }
  out <- c(out, "", paste("Platform", gpl_id, "annotation columns:"), paste(names(annot), collapse = " | "),
           "", "First rows:", utils::capture.output(print(utils::head(annot, 5))),
           "", "Non-empty values per column:",
           utils::capture.output(print(vapply(annot, function(x) sum(!is.na(x) & x != ""), 1L))))
  writeLines(out, "logs/gse72094_fields.txt")
  stop("First run done: send logs/gse72094_fields.txt so the GSE72094 entry of config/cohorts.yaml can be filled.",
       call. = FALSE)
}

summary_file <- "results/qc/GSE72094_build_summary.txt"
report <- c("GSE72094 build summary", paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")))
add <- function(title, x) {
  report <<- c(report, "", paste("==", title), utils::capture.output(print(x)))
  writeLines(report, summary_file)
}
count_table <- function(...) {
  d <- as.data.frame(table(..., useNA = "ifany"), responseName = "n")
  d[d$n > 0, , drop = FALSE]
}

# 1. Clinical data
if (!is.null(f72$include)) {
  keep <- tolower(trimws(pd[[f72$include$field]])) %in% tolower(f72$include$values)
  add("Samples selected / all", c(selected = sum(keep), all = nrow(pd)))
  pd <- pd[keep, , drop = FALSE]
}
clin <- make_clin_geo(pd, f72)
add("Event (1 = dead) by raw value", count_table(raw = pd[[f72$event$field]], event = clin$os_event))
add("Stage: raw value -> group", count_table(raw = clin$stage_raw, group = clin$stage))
add("Follow-up in months (summary)", summary(clin$os_months))
add("Age recorded / sex", list(age_recorded = sum(!is.na(clin$age)), sex = table(clin$sex, useNA = "ifany")))

# 2. Expression and probe -> Ensembl gene map
expr <- Biobase::exprs(g)[, pd$geo_accession, drop = FALSE]
if (max(expr, na.rm = TRUE) > 30) {
  expr <- log2(pmax(expr, 1))
  add("Expression", "values were not on a log scale: log2(max(x, 1)) applied")
}
na_probe <- rowSums(is.na(expr)) > 0
add("Probes with a missing value (dropped)", sum(na_probe))
expr <- expr[!na_probe, , drop = FALSE]

col <- f72$probe_annotation$column
if (is.null(col) || !col %in% names(annot)) stop("probe_annotation$column not found in the ", gpl_id, " table")
ids <- annotation_ids(annot[[col]])
multi <- !is.na(ids) & grepl("///|;|,", ids)
pm <- data.frame(probe_id = as.character(annot$ID), raw_id = ids, stringsAsFactors = FALSE)
pm <- pm[!multi & !is.na(pm$raw_id) & pm$raw_id != "", , drop = FALSE]
add("Probes: annotated with one gene / with several genes (dropped) / none",
    c(one_gene = nrow(pm), several = sum(multi), none = length(ids) - nrow(pm) - sum(multi)))
if (identical(f72$probe_annotation$id_type, "ensembl")) {
  pm$ensembl_id <- strip_version(pm$raw_id)
} else if (identical(f72$probe_annotation$id_type, "entrez")) {
  em <- entrez_ensembl_map()
  em <- unique(em)
  n_ens <- table(em$entrez_id)
  add("Entrez IDs on the array: mapped to one Ensembl gene / to several (dropped) / to none",
      c(one = sum(unique(pm$raw_id) %in% names(n_ens)[n_ens == 1]),
        several = sum(unique(pm$raw_id) %in% names(n_ens)[n_ens > 1]),
        none = sum(!unique(pm$raw_id) %in% em$entrez_id)))
  em <- em[em$entrez_id %in% names(n_ens)[n_ens == 1], , drop = FALSE]   # one gene only
  pm$ensembl_id <- em$ensembl_id[match(pm$raw_id, em$entrez_id)]
} else {
  stop("probe_annotation$id_type must be 'ensembl' or 'entrez' (gene symbols are not used)")
}
pm <- pm[!is.na(pm$ensembl_id), c("probe_id", "ensembl_id")]
col_res <- collapse_probes(expr[intersect(rownames(expr), pm$probe_id), , drop = FALSE], pm)
add("Probes mapped to one Ensembl gene / genes kept (one probe per gene, highest mean)",
    c(probes = sum(rownames(expr) %in% pm$probe_id), genes = nrow(col_res$expr)))

# 3. Cohort object
meta <- list(source = "GEO series matrix (submitter's normalised values)", platform_id = gpl_id,
             probe_annotation = paste(col, f72$probe_annotation$id_type), downloaded = format(Sys.Date()),
             probe_choice = col_res$probe_choice)
co <- finalise_cohort(col_res$expr, clin, gse, platform = gpl_id, meta = meta)
check_cohort(co)
saveRDS(co, "data/processed/GSE72094.rds")
add("Sample flow (patients)", co$meta$sample_flow)
add("Cohort", c(patients = ncol(co$expr), deaths = sum(co$event),
                median_follow_up_months = round(reverse_km_median(co$time, co$event), 1),
                max_follow_up_months = round(max(co$time), 1), genes = nrow(co$expr)))
add("Stage group", table(co$clin$stage, useNA = "ifany"))
add("Expression range", round(range(co$expr), 2))
cat(report, sep = "\n")
message("\nGSE72094: ", ncol(co$expr), " patients, ", sum(co$event), " deaths. Summary saved to ", summary_file)
