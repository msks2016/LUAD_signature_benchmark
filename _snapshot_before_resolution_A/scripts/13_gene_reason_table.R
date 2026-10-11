# ============================================================
# 13_gene_reason_table.R (FIXED)
# Per-gene reason table: why each signature gene is missing
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
})

# --- Load signature genes + mapping ---
sigs <- read.csv("data/processed/signature_coefficients.csv")
sig_genes <- unique(sigs$gene)
gmap <- read.csv("docs/gene_mapping.csv")
paper_to_ens <- setNames(as.character(gmap$paper_ensembl), as.character(gmap$paper_name))

cat("Loaded", length(sig_genes), "unique signature genes\n")
cat("Loaded", nrow(gmap), "gene mappings\n\n")

# --- Load GSE81089 counts ---
counts_gse <- readRDS("data/raw/geo/GSE81089_counts.rds")
gse_ids <- counts_gse[[1]]
gse_ids_base <- sub("\\..*$", "", gse_ids)
cat("GSE81089 has", length(gse_ids), "genes\n\n")

# --- Categorize each signature gene ---
categorize_gene <- function(gene) {
  # 1. Direct symbol match
  if (gene %in% gse_ids) return("present_by_symbol")
  
  # 2. Ensembl base match
  eid <- paper_to_ens[gene]
  if (!is.na(eid) && nchar(eid) > 0) {
    eid_base <- sub("\\..*$", "", eid)
    if (eid_base %in% gse_ids_base) return("present_by_ensembl_base")
  }
  
  # 3. Not found — categorize by name pattern
  if (grepl("^AC[0-9]|^AL[0-9]|^RP11-", gene)) {
    return("clone_based_lncrna_likely_unannotated")
  }
  if (grepl("^LINC", gene)) {
    return("LINC_unannotated_in_cohort")
  }
  if (grepl("-AS1$|-DT$|-AS2$", gene)) {
    return("antisense_lncrna_likely_unannotated")
  }
  return("absent_no_annotation")
}

reason_df <- data.frame(
  cohort = "GSE81089",
  gene = sig_genes,
  reason = sapply(sig_genes, categorize_gene),
  stringsAsFactors = FALSE
)

cat("=== Reason summary ===\n")
print(table(reason_df$reason))

cat("\n=== First 30 missing genes ===\n")
missing_df <- reason_df %>% filter(!grepl("^present_", reason))
print(head(missing_df, 30))

# Save
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(reason_df, "results/tables/gene_reason_table_GSE81089.csv", row.names = FALSE)
cat("\nSaved: results/tables/gene_reason_table_GSE81089.csv\n")
# 🅁 R Console
reason_df <- read.csv("results/tables/gene_reason_table_GSE81089.csv")
cat("=== Summary by reason ===\n")
print(table(reason_df$reason))

cat("\n=== All genes with 'present' status ===\n")
present_genes <- reason_df[grepl("^present", reason_df$reason), ]
print(present_genes)

cat("\n=== Missing genes grouped by category ===\n")
missing_genes <- reason_df[!grepl("^present", reason_df$reason), ]
for (cat in unique(missing_genes$reason)) {
  cat("\n---", cat, "---\n")
  print(missing_genes$gene[missing_genes$reason == cat])
}
reason_df <- read.csv("results/tables/gene_reason_table_GSE81089.csv")
cat("=== Summary by reason ===\n")
print(table(reason_df$reason))

cat("\n=== Present genes ===\n")
present_genes <- reason_df[grepl("^present", reason_df$reason), ]
print(present_genes)

cat("\n=== Missing genes by category ===\n")
missing_genes <- reason_df[!grepl("^present", reason_df$reason), ]
print(table(missing_genes$reason))
reason_df <- read.csv("results/tables/gene_reason_table_GSE81089.csv")
cat("Summary:\n"); print(table(reason_df$reason))
present_genes <- reason_df[grepl("^present", reason_df$reason), ]
cat("\nPresent:", nrow(present_genes), "genes\n")
missing_genes <- reason_df[!grepl("^present", reason_df$reason), ]
cat("\nMissing by category:\n"); print(table(missing_genes$reason))
for (cat_name in unique(missing_genes$reason)) {
  cat("\n---", cat_name, "---\n")
  print(missing_genes$gene[missing_genes$reason == cat_name])
}
getwd()
file.exists("results/tables/gene_reason_table_GSE81089.csv")
list.files("results/tables/", pattern = "gene_reason|GSE81089")
reason_df <- read.csv("results/tables/gene_reason_table_GSE81089.csv"); print(table(reason_df$reason)); missing_genes <- reason_df[!grepl("^present", reason_df$reason), ]; cat("\nMissing count:", nrow(missing_genes), "\n")
# 🅁 R Console — Which 2 genes recover?
present <- reason_df[reason_df$reason == "present_by_ensembl_base", ]
print(present)
# 🅁 R Console
suppressPackageStartupMessages({ library(TCGAbiolinks) })

# CPTAC-3 is a separate GDC project
tryCatch({
  # List available projects to find the right code
  projects <- getGDCprojects()
  cat("=== Available projects with 'CPTAC' in name ===\n")
  print(projects[grepl("CPTAC|LUAD", projects$project_id, ignore.case = TRUE),
                 c("project_id", "name")])
}, error = function(e) cat("Failed:", conditionMessage(e), "\n"))
# 🅁 R Console
suppressPackageStartupMessages({ library(TCGAbiolinks) })

# Check available functions with "project" in name
cat("=== Functions containing 'project' ===\n")
print(grep("project", ls("package:TCGAbiolinks"), value = TRUE, ignore.case = TRUE))
library(TCGAbiolinks)
grep("project", ls("package:TCGAbiolinks"), value = TRUE, ignore.case = TRUE)
tryCatch({
  clin <- GDCquery_clinic(project = "CPTAC-3", type = "clinical")
  cat("Rows:", nrow(clin), "\n")
  print(colnames(clin)[1:15])
  if ("vital_status" %in% colnames(clin)) print(table(clin$vital_status, useNA = "ifany"))
}, error = function(e) cat("Failed:", conditionMessage(e), "\n"))
# 🅁 R Console
proj <- getGDCprojects()
cat("Total projects:", nrow(proj), "\n")
cat("\nProjects with CPTAC or LUAD:\n")
print(proj[grepl("CPTAC|LUAD", proj$project_id, ignore.case = TRUE),
           c("project_id", "name")])
# 🅁 R Console
projects_df <- getGDCprojects()
cat("Total projects:", nrow(projects_df), "\n")
# 🅁 R Console
result <- getGDCprojects()
# 🅁 R Console
class(result)
nrow(result)
# 🅁 R Console
print(result$project_id)
# 🅁 R Console
tryCatch({
  clin_cptac <- GDCquery_clinic(project = "CPTAC-3", type = "clinical")
  cat("Rows:", nrow(clin_cptac), "\n")
  cat("Columns:", ncol(clin_cptac), "\n\n")
  print(colnames(clin_cptac))
}, error = function(e) cat("Failed:", conditionMessage(e), "\n"))
# 🅁 R Console
cat("=== CPTAC-3 ===\n")
cat("Total rows:", nrow(clin_cptac), "\n")
cat("Unique patients:", length(unique(clin_cptac$submitter_id)), "\n\n")

cat("=== Vital status ===\n")
print(table(clin_cptac$vital_status, useNA = "ifany"))

cat("\n=== Death days available ===\n")
cat("Non-NA days_to_death:", sum(!is.na(clin_cptac$days_to_death)), "\n")

cat("\n=== Sex ===\n")
print(table(clin_cptac$sex_at_birth, useNA = "ifany"))

cat("\n=== Age summary ===\n")
print(summary(clin_cptac$age_at_index))
cat("=== CPTAC-3 ===\n")
cat("Total rows:", nrow(clin_cptac), "\n")
cat("Unique patients:", length(unique(clin_cptac$submitter_id)), "\n\n")

cat("=== Vital status ===\n")
print(table(clin_cptac$vital_status, useNA = "ifany"))

cat("\n=== Death days available ===\n")
cat("Non-NA days_to_death:", sum(!is.na(clin_cptac$days_to_death)), "\n")

cat("\n=== Sex ===\n")
print(table(clin_cptac$sex_at_birth, useNA = "ifany"))

cat("\n=== Age summary ===\n")
print(summary(clin_cptac$age_at_index))
# 🅁 R Console — CPTAC-3 structure check
cat("=== CPTAC-3 ===\n")
cat("Total rows:", nrow(clin_cptac), "\n")
cat("Unique patients (submitter_id):", length(unique(clin_cptac$submitter_id)), "\n\n")

cat("=== Vital status ===\n")
print(table(clin_cptac$vital_status, useNA = "ifany"))

cat("\n=== days_to_death non-NA ===\n")
cat(sum(!is.na(clin_cptac$days_to_death)), "\n")

cat("\n=== Sex ===\n")
print(table(clin_cptac$sex_at_birth, useNA = "ifany"))

cat("\n=== Age non-NA ===\n")
cat(sum(!is.na(clin_cptac$age_at_index)), "\n")
cat("Age NAs:", sum(is.na(clin_cptac$age_at_index)), "\n")

cat("\n=== Duplicate submitter_id? ===\n")
cat("Rows per patient (table summary):\n")
print(summary(as.numeric(table(clin_cptac$submitter_id))))
# 🅁 R Console — CPTAC-3 essentials
cat("Total rows:", nrow(clin_cptac), "\n")
cat("Unique patients:", length(unique(clin_cptac$submitter_id)), "\n\n")

cat("=== vital_status ===\n")
print(table(clin_cptac$vital_status, useNA = "ifany"))

cat("\n=== days_to_death non-NA ===\n")
cat(sum(!is.na(clin_cptac$days_to_death)), "\n")

cat("\n=== days_to_last_follow_up present? ===\n")
if ("days_to_last_follow_up" %in% colnames(clin_cptac)) {
  cat("non-NA:", sum(!is.na(clin_cptac$days_to_last_follow_up)), "\n")
} else {
  cat("column not present\n")
}

cat("\n=== sex_at_birth ===\n")
print(table(clin_cptac$sex_at_birth, useNA = "ifany"))

cat("\n=== age_at_index non-NA ===\n")
cat(sum(!is.na(clin_cptac$age_at_index)), "\n")
# 🅁 R Console — CPTAC-3 survival + composition
cat("Total rows:", nrow(clin_cptac), "\n")
cat("Unique patients:", length(unique(clin_cptac$submitter_id)), "\n\n")

cat("=== vital_status ===\n")
print(table(clin_cptac$vital_status, useNA = "ifany"))

cat("\n=== days_to_death non-NA ===\n")
cat(sum(!is.na(clin_cptac$days_to_death)), "\n")

cat("\n=== days_to_last_follow_up non-NA ===\n")
if ("days_to_last_follow_up" %in% colnames(clin_cptac)) {
  cat(sum(!is.na(clin_cptac$days_to_last_follow_up)), "\n")
} else cat("column absent\n")

cat("\n=== sex_at_birth ===\n")
print(table(clin_cptac$sex_at_birth, useNA = "ifany"))

cat("\n=== cancer-type columns ===\n")
print(grep("disease|diagnos|cancer|tumor|site|project",
           colnames(clin_cptac), value = TRUE, ignore.case = TRUE))
# 🅁 R Console — which columns actually vary + LUAD identification
key_cols <- c("primary_diagnosis", "site_of_resection_or_biopsy",
              "last_known_disease_status", "days_to_last_known_disease_status",
              "vital_status", "days_to_death", "days_to_last_follow_up",
              "sex_at_birth", "age_at_index")

for (col in key_cols) {
  if (col %in% colnames(clin_cptac)) {
    x <- clin_cptac[[col]]
    n_nonNA <- sum(!is.na(x))
    n_unique <- length(unique(x[!is.na(x)]))
    cat("\n===", col, "=== non-NA:", n_nonNA, " unique:", n_unique, "\n")
    if (n_unique <= 30) print(table(x, useNA = "ifany"))
    else print(head(sort(table(x), decreasing = TRUE), 15))
  } else {
    cat("\n===", col, "=== COLUMN ABSENT\n")
  }
}
# 🅁 R Console — which columns actually vary + LUAD identification
key_cols <- c("primary_diagnosis", "site_of_resection_or_biopsy",
              "last_known_disease_status", "days_to_last_known_disease_status",
              "vital_status", "days_to_death", "days_to_last_follow_up",
              "sex_at_birth", "age_at_index")

for (col in key_cols) {
  if (col %in% colnames(clin_cptac)) {
    x <- clin_cptac[[col]]
    n_nonNA <- sum(!is.na(x))
    n_unique <- length(unique(x[!is.na(x)]))
    cat("\n===", col, "=== non-NA:", n_nonNA, " unique:", n_unique, "\n")
    if (n_unique <= 30) {
      print(table(x, useNA = "ifany"))
    } else {
      print(head(sort(table(x), decreasing = TRUE), 15))
    }
  } else {
    cat("\n===", col, "=== COLUMN ABSENT\n")
  }
}
key_cols <- c("primary_diagnosis", "site_of_resection_or_biopsy",
              "last_known_disease_status", "days_to_last_known_disease_status",
              "vital_status", "days_to_death", "days_to_last_follow_up",
              "sex_at_birth", "age_at_index")

for (col in key_cols) {
  if (col %in% colnames(clin_cptac)) {
    x <- clin_cptac[[col]]
    n_nonNA <- sum(!is.na(x))
    n_unique <- length(unique(x[!is.na(x)]))
    cat("\n===", col, "=== non-NA:", n_nonNA, " unique:", n_unique, "\n")
    if (n_unique <= 30) {
      print(table(x, useNA = "ifany"))
    } else {
      print(head(sort(table(x), decreasing = TRUE), 15))
    }
  } else {
    cat("\n===", col, "=== COLUMN ABSENT\n")
  }
}
key_cols <- c("primary_diagnosis", "vital_status", "days_to_death",
              "days_to_last_follow_up", "sex_at_birth", "age_at_index")
length(key_cols)
key_cols <- c("primary_diagnosis", "vital_status", "days_to_death",
              "days_to_last_follow_up", "sex_at_birth", "age_at_index")

for (col in key_cols) {
  if (col %in% colnames(clin_cptac)) {
    x <- clin_cptac[[col]]
    n_nonNA <- sum(!is.na(x))
    n_unique <- length(unique(x[!is.na(x)]))
    cat("\n===", col, "=== non-NA:", n_nonNA, " unique:", n_unique, "\n")
    if (n_unique <= 30) {
      print(table(x, useNA = "ifany"))
    } else {
      print(head(sort(table(x), decreasing = TRUE), 15))
    }
  } else {
    cat("\n===", col, "=== COLUMN ABSENT\n")
  }
}
key_cols <- c("primary_diagnosis", "vital_status", "days_to_death",
              "days_to_last_follow_up", "sex_at_birth", "age_at_index")
length(key_cols)
key_cols <- c("primary_diagnosis", "vital_status", "days_to_death",
              "days_to_last_follow_up", "sex_at_birth", "age_at_index")
length(key_cols)
for (col in key_cols) {
  if (col %in% colnames(clin_cptac)) {
    x <- clin_cptac[[col]]
    n_nonNA <- sum(!is.na(x))
    n_unique <- length(unique(x[!is.na(x)]))
    cat("\n===", col, "=== non-NA:", n_nonNA, " unique:", n_unique, "\n")
    if (n_unique <= 30) {
      print(table(x, useNA = "ifany"))
    } else {
      print(head(sort(table(x), decreasing = TRUE), 15))
    }
  } else {
    cat("\n===", col, "=== COLUMN ABSENT\n")
  }
}
cat("\n=== project_id ===\n"); print(table(clin_cptac$project_id, useNA = "ifany"))
cat("\n=== site_of_resection_or_biopsy ===\n"); print(table(clin_cptac$site_of_resection_or_biopsy, useNA = "ifany"))
cat("\n=== disease_type (if present) ===\n")
if ("disease_type" %in% colnames(clin_cptac)) print(table(clin_cptac$disease_type, useNA="ifany")) else cat("absent\n")
cat("=== site_of_resection_or_biopsy ===\n")
print(table(clin_cptac$site_of_resection_or_biopsy, useNA = "ifany"))

cat("\n=== project_id ===\n")
print(table(clin_cptac$project_id, useNA = "ifany"))

cat("\n=== vital_status ===\n")
print(table(clin_cptac$vital_status, useNA = "ifany"))

cat("\n=== days_to_death non-NA ===\n")
cat(sum(!is.na(clin_cptac$days_to_death)), "\n")

cat("\n=== days_to_last_follow_up non-NA ===\n")
cat(sum(!is.na(clin_cptac$days_to_last_follow_up)), "\n")
cat("=== site_of_resection_or_biopsy ===\n")
print(table(clin_cptac$site_of_resection_or_biopsy, useNA = "ifany"))

cat("\n=== project_id ===\n")
print(table(clin_cptac$project_id, useNA = "ifany"))

cat("\n=== vital_status ===\n")
print(table(clin_cptac$vital_status, useNA = "ifany"))

cat("\n=== days_to_death non-NA ===\n")
cat(sum(!is.na(clin_cptac$days_to_death)), "\n")
cat("=== site_of_resection_or_biopsy ===\n")
print(table(clin_cptac$site_of_resection_or_biopsy, useNA = "ifany"))

cat("\n=== tissue_or_organ_of_origin ===\n")
if ("tissue_or_organ_of_origin" %in% colnames(clin_cptac)) {
  print(table(clin_cptac$tissue_or_organ_of_origin, useNA = "ifany"))
} else cat("absent\n")

cat("\n=== primary_site ===\n")
if ("primary_site" %in% colnames(clin_cptac)) {
  print(table(clin_cptac$primary_site, useNA = "ifany"))
} else cat("absent\n")

cat("\n=== project_id ===\n")
print(table(clin_cptac$project_id, useNA = "ifany"))
cat("=== site_of_resection_or_biopsy ===\n")
print(table(clin_cptac$site_of_resection_or_biopsy, useNA = "ifany"))

cat("\n=== tissue_or_organ_of_origin ===\n")
if ("tissue_or_organ_of_origin" %in% colnames(clin_cptac))
  print(table(clin_cptac$tissue_or_organ_of_origin, useNA="ifany")) else cat("absent\n")

cat("\n=== primary_site ===\n")
if ("primary_site" %in% colnames(clin_cptac))
  print(table(clin_cptac$primary_site, useNA="ifany")) else cat("absent\n")

cat("\n=== disease_type ===\n")
if ("disease_type" %in% colnames(clin_cptac))
  print(table(clin_cptac$disease_type, useNA="ifany")) else cat("absent\n")
cat("=== site_of_resection_or_biopsy ===\n")
print(table(clin_cptac$site_of_resection_or_biopsy, useNA = "ifany"))

cat("\n=== tissue_or_organ_of_origin ===\n")
if ("tissue_or_organ_of_origin" %in% colnames(clin_cptac)) print(table(clin_cptac$tissue_or_organ_of_origin, useNA="ifany"))

cat("\n=== primary_site ===\n")
if ("primary_site" %in% colnames(clin_cptac)) print(table(clin_cptac$primary_site, useNA="ifany"))

cat("\n=== disease_type ===\n")
if ("disease_type" %in% colnames(clin_cptac)) print(table(clin_cptac$disease_type, useNA="ifany"))
library(TCGAbiolinks)

cptac_sheet <- GDCquery(
  project       = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

sheet <- getResults(cptac_sheet)
cat("Sample-sheet rows:", nrow(sheet), "\n\n")

# Look at every column that could separate LUAD
for (col in colnames(sheet)) {
  if (grepl("disease|site|tissue|project|cancer|primary|sample_type|tumor",
            col, ignore.case = TRUE)) {
    cat("\n===", col, "===\n")
    print(table(sheet[[col]], useNA = "ifany"))
  }
}
sheet <- getResults(cptac_sheet)
cat("rows:", nrow(sheet), "\n")

cat("\n=== project_id ===\n");    print(table(sheet$project_id, useNA="ifany"))
cat("\n=== sample_type ===\n");   print(table(sheet$sample_type, useNA="ifany"))

for (col in c("disease_type", "primary_site", "tissue_or_organ_of_origin",
              "site_of_resection_or_biopsy")) {
  if (col %in% colnames(sheet)) {
    cat("\n===", col, "===\n")
    print(table(sheet[[col]], useNA="ifany"))
  }
}
library(TCGAbiolinks)
cptac_sheet <- GDCquery(
  project       = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
# Access the raw results, not the pretty helper
sheet_raw <- cptac_sheet$results[[1]]

cat("rows:", nrow(sheet_raw), "\n")
cat("cols:", ncol(sheet_raw), "\n")
print(colnames(sheet_raw))
library(jsonlite)

# GDC API filter for CPTAC-3 RNA-seq STAR-Counts files
filters <- '{"op":"and","content":[
  {"op":"in","content":{"field":"cases.project.project_id","value":["CPTAC-3"]}},
  {"op":"in","content":{"field":"data_category","value":["Transcriptome Profiling"]}},
  {"op":"in","content":{"field":"data_type","value":["Gene Expression Quantification"]}},
  {"op":"in","content":{"field":"analysis.workflow_type","value":["STAR - Counts"]}}
]}'

url <- paste0(
  "https://api.gdc.cancer.gov/files?",
  "filters=", URLencode(filters, reserved = TRUE),
  "&expand=cases,cases.samples",
  "&size=2000",
  "&format=JSON"
)

res <- fromJSON(url, flatten = TRUE)
cat("Total files:", res$data$pagination$total, "\n")
library(jsonlite)
filters <- '{"op":"and","content":[{"op":"in","content":{"field":"cases.project.project_id","value":["CPTAC-3"]}},{"op":"in","content":{"field":"data_category","value":["Transcriptome Profiling"]}},{"op":"in","content":{"field":"data_type","value":["Gene Expression Quantification"]}}]}'
url <- paste0("https://api.gdc.cancer.gov/files?filters=", URLencode(filters, reserved=TRUE), "&expand=cases,cases.samples&size=2000&format=JSON")
nchar(url)
res <- fromJSON(url, flatten=TRUE)
cat("Total files:", res$data$pagination$total, "\n")
files <- res$data$hits
dim(files)
colnames(files)
library(tidyr)
cases_flat <- unnest(files, cols = c(cases), names_sep = "_")
cases_list <- files$cases
first_case <- cases_list[[1]]
print(colnames(first_case))
print(first_case)
dim(cases_flat)
colnames(cases_flat)
print(table(cases_flat$cases_primary_site, useNA="ifany"))
print(table(cases_flat$cases_disease_type, useNA="ifany"))
length(unique(cases_flat$cases_submitter_id))
head(unique(cases_flat$cases_submitter_id))
luad_ids <- unique(cases_flat$cases_submitter_id[cases_flat$cases_primary_site == "Bronchus and lung"])
cat("LUAD patients in CPTAC-3:", length(luad_ids), "\n")
print(table(cases_flat$cases_primary_site, useNA="ifany"))
print(table(cases_flat$cases_primary_site)[
  grep("lung|bronchus|pulmonary", 
       names(table(cases_flat$cases_primary_site)), 
       ignore.case=TRUE)
])
luad_ids <- unique(cases_flat$cases_submitter_id[
  cases_flat$cases_primary_site == "Bronchus and lung" &
    cases_flat$cases_disease_type == "Adenomas and Adenocarcinomas"
])
cat("LUAD patients:", length(luad_ids), "\n")
luad_ids <- unique(cases_flat$cases_submitter_id[
  cases_flat$cases_primary_site == "Bronchus and lung" &
    cases_flat$cases_disease_type == "Adenomas and Adenocarcinomas"
])
cat("LUAD patients (page 1):", length(luad_ids), "\n")
clin_luad <- clin_cptac[clin_cptac$submitter_id %in% luad_ids, ]
cat("Matched in clin_cptac:", nrow(clin_luad), "\n")
cat("Deaths:", sum(clin_luad$vital_status == "Dead", na.rm=TRUE), "\n")
cat("Alive:", sum(clin_luad$vital_status == "Alive", na.rm=TRUE), "\n")
table(clin_luad$vital_status, useNA="ifany")
# 1. Fetch all 3,042 files
url_all <- sub("size=2000", "size=4000", url)
res_all <- fromJSON(url_all, flatten=TRUE)
files_all <- res_all$data$hits
cat("Fetched:", nrow(files_all), "\n")

# 2. Unnest cases
library(tidyr)
cases_all <- unnest(files_all, cols = c(cases), names_sep = "_")

# 3. Filter to LUAD + RNA-seq + primary tumor
luad_manifest <- cases_all[
  cases_all$cases_primary_site == "Bronchus and lung" &
    cases_all$cases_disease_type == "Adenomas and Adenocarcinomas" &
    cases_all$experimental_strategy == "RNA-Seq",
]

cat("LUAD RNA-seq files:", nrow(luad_manifest), "\n")
cat("Unique LUAD patients:", length(unique(luad_manifest$cases_submitter_id)), "\n")

# 4. Save it
write.csv(luad_manifest, "data/raw/cptac3_luad_manifest.csv", row.names=FALSE)
url_all <- sub("size=2000", "size=4000", url)
res_all <- fromJSON(url_all, flatten=TRUE)
files_all <- res_all$data$hits
cat("Fetched:", nrow(files_all), "\n")
library(tidyr)
cases_all <- unnest(files_all, cols = c(cases), names_sep = "_")
cat("cols:", ncol(cases_all), "rows:", nrow(cases_all), "\n")
cat("primary_site present:", "cases_primary_site" %in% colnames(cases_all), "\n")
cat("disease_type present:", "cases_disease_type" %in% colnames(cases_all), "\n")
cat("submitter_id present:", "cases_submitter_id" %in% colnames(cases_all), "\n")
cat("experimental_strategy present:", "experimental_strategy" %in% colnames(cases_all), "\n")
luad_manifest <- cases_all[
  cases_all$cases_primary_site == "Bronchus and lung" &
    cases_all$cases_disease_type == "Adenomas and Adenocarcinomas" &
    cases_all$experimental_strategy == "RNA-Seq",
]
cat("LUAD RNA-seq files:", nrow(luad_manifest), "\n")
cat("Unique LUAD patients:", length(unique(luad_manifest$cases_submitter_id)), "\n")
1 + 1
luad_manifest <- cases_flat[
  cases_flat$cases_primary_site == "Bronchus and lung" &
    cases_flat$cases_disease_type == "Adenomas and Adenocarcinomas" &
    cases_flat$experimental_strategy == "RNA-Seq",
]
nrow(luad_manifest)
length(unique(luad_manifest$cases_submitter_id))
dir.create("data/raw", showWarnings=FALSE, recursive=TRUE)
write.csv(luad_manifest, "data/raw/cptac3_luad_manifest.csv", row.names=FALSE)
luad_manifest$cases <- NULL
write.csv(luad_manifest, "data/raw/cptac3_luad_manifest.csv", row.names=FALSE)
list_cols <- sapply(luad_manifest, is.list)
print(names(luad_manifest)[list_cols])
luad_manifest$cases_samples <- NULL
luad_manifest$acl <- NULL
sapply(luad_manifest, is.list)
write.csv(luad_manifest, "data/raw/cptac3_luad_manifest.csv", row.names=FALSE)
file.exists("data/raw/cptac3_luad_manifest.csv")
saved <- read.csv("data/raw/cptac3_luad_manifest.csv")
dim(saved)
cat("Files:", nrow(luad_manifest), "\n")
cat("Patients:", length(unique(luad_manifest$cases_submitter_id)), "\n")
cat("Has file_id:", "file_id" %in% colnames(luad_manifest), "\n")
cat("Has submitter_id:", "cases_submitter_id" %in% colnames(luad_manifest), "\n")
# =====================================================================
# 08_prepare_CPTAC3.R
# Purpose: Download and harmonize CPTAC-3 LUAD RNA-seq + clinical data
# Output:  data/processed/CPTAC3_LUAD_expr.rds
#          data/processed/CPTAC3_LUAD_clin.rds
# =====================================================================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(edgeR)
  library(dplyr)
  library(tibble)
})

set.seed(42)

# ---- 1. Load manifest + clinical ------------------------------------
manifest <- read.csv("data/raw/cptac3_luad_manifest.csv",
                     stringsAsFactors = FALSE)
cat("Manifest files:", nrow(manifest), "\n")
cat("Unique patients:", length(unique(manifest$cases_submitter_id)), "\n")

# clin_cptac should already exist from GDCquery_clinic("CPTAC-3", "clinical")
# If not, load it:
if (!exists("clin_cptac")) {
  clin_cptac <- GDCquery_clinic("CPTAC-3", type = "clinical")
}

# Restrict clinical to LUAD patients
luad_ids <- unique(manifest$cases_submitter_id)
clin_luad <- clin_cptac[clin_cptac$submitter_id %in% luad_ids, ]
cat("LUAD clinical rows:", nrow(clin_luad), "\n")

# ---- 2. Build survival object --------------------------------------
clin_luad <- clin_luad %>%
  mutate(
    os_time_days = ifelse(vital_status == "Dead",
                          days_to_death,
                          days_to_last_follow_up),
    os_time      = os_time_days / 30.44,       # months
    os_event     = as.integer(vital_status == "Dead")
  )

# Exclude "Not Reported"
clin_luad <- clin_luad[clin_luad$vital_status %in% c("Alive", "Dead"), ]
clin_luad <- clin_luad[!is.na(clin_luad$os_time) & clin_luad$os_time > 0, ]

cat("Evaluable patients:", nrow(clin_luad), "\n")
cat("Events:", sum(clin_luad$os_event), "\n")
print(summary(clin_luad$os_time))

# ---- 3. Download STAR - Counts -------------------------------------
# The trick: bypass GDCquery's strict column checker by using the
# manifest's file_id directly. GDCdownload accepts a query object,
# so we construct a minimal one, then GDCdownload fetches by file UUID.
query <- GDCquery(
  project       = "CPTAC-3",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode       = luad_ids
)

# Restrict to files present in our manifest (drops any normal-tissue files)
query$results[[1]] <- query$results[[1]][
  query$results[[1]]$file_id %in% manifest$file_id, 
]

cat("Files to download:", nrow(query$results[[1]]), "\n")

GDCdownload(query, directory = "data/raw/CPTAC3_LUAD/",
            files.per.chunk = 50)

# ---- 4. Prepare expression matrix ----------------------------------
se <- GDCprepare(query, directory = "data/raw/CPTAC3_LUAD/")
saveRDS(se, "data/raw/CPTAC3_LUAD/SummarizedExperiment.rds")

counts <- assay(se, "unstranded")           # STAR unstranded counts
cat("Raw matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ---- 5. Sample → patient map ---------------------------------------
sample_meta <- as.data.frame(colData(se)) %>%
  select(barcode, patient, sample_type)
sample_meta$submitter_id <- substr(sample_meta$patient, 1, 12)

# Keep primary tumors only, one sample per patient (highest total counts)
sample_meta$total_counts <- colSums(counts[, sample_meta$barcode, drop = FALSE])
sample_meta <- sample_meta %>%
  filter(grepl("Primary Tumor", sample_type, ignore.case = TRUE)) %>%
  arrange(submitter_id, desc(total_counts)) %>%
  distinct(submitter_id, .keep_all = TRUE)

cat("Primary-tumor samples (one per patient):", nrow(sample_meta), "\n")

counts <- counts[, sample_meta$barcode, drop = FALSE]

# ---- 6. Harmonize genes --------------------------------------------
# If you used GENCODE v36 for TCGA, keep it here too. Restrict to
# protein_coding + lncRNA, matching your signature gene space.
rowData_df <- as.data.frame(rowData(se))
keep <- rowData_df$gene_type %in% c("protein_coding", "lncRNA")
counts <- counts[keep, ]
rownames(counts) <- rowData_df$gene_name[keep]

# Collapse duplicate gene symbols: keep highest mean count
counts <- counts[order(-rowMeans(counts)), ]
counts <- counts[!duplicated(rownames(counts)), ]
cat("After gene filter:", nrow(counts), "genes x", ncol(counts), "samples\n")

# ---- 7. Normalize: log2-CPM (matches TCGA pipeline) ----------------
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

saveRDS(logcpm, "data/processed/CPTAC3_LUAD_expr.rds")

# ---- 8. Save clinical -----------------------------------------------
clin_save <- clin_luad %>%
  select(submitter_id, os_time, os_event,
         sex_at_birth, age_at_index,
         ajcc_pathologic_stage) %>%
  filter(submitter_id %in% colnames(logcpm))

saveRDS(clin_save, "data/processed/CPTAC3_LUAD_clin.rds")

cat("\nDone. Files saved:\n")
cat("  data/processed/CPTAC3_LUAD_expr.rds   (",
    nrow(logcpm), "genes x", ncol(logcpm), "patients )\n")
cat("  data/processed/CPTAC3_LUAD_clin.rds   (",
    nrow(clin_save), "patients )\n")
1 + 1
getwd()
file.exists("data/raw/cptac3_luad_manifest.csv")
dir.create("scripts", showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("data/raw/CPTAC3_LUAD", showWarnings = FALSE, recursive = TRUE)
source("scripts/08_prepare_CPTAC3.R")
