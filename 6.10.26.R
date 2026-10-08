Bismillah
install.packages("renv")
renv::init()
install.packages(c("tidyverse", "survival", "survminer", "timeROC", "metafor", "rms", "pROC", "ggplot2"))
Y
install.packages(c("tidyverse", "survival", "survminer", "timeROC", "metafor", "rms", "pROC", "ggplot2"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "GEOquery", "limma", "biomaRt", "ComplexHeatmap"))
renv::snapshot()
R --version
dir.create("data/raw/TCGA-LUAD/GDCdata", recursive = TRUE, showWarnings = FALSE)
query_rna <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
"TCGAbiolinks" %in% rownames(installed.packages())
library(TCGAbiolinks)
library(SummarizedExperiment)
packageVersion("TCGAbiolinks")
exists("GDCquery")
cat("TCGAbiolinks version:", as.character(packageVersion("TCGAbiolinks")), "\n")
dir.create("data/raw/TCGA-LUAD/GDCdata", recursive = TRUE, showWarnings = FALSE)
query_rna <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
cat("Sample types:\n")
print(table(getResults(query_rna)$sample_type))
GDCdownload(query_rna, method = "client",
            directory = "data/raw/TCGA-LUAD/GDCdata")
se_rna <- GDCprepare(query_rna, directory = "data/raw/TCGA-LUAD/GDCdata")
saveRDS(se_rna, "data/raw/TCGA-LUAD/se_rna.rds")
tpm <- assay(se_rna, "tpm_unstrand")
cat("TPM matrix:", dim(tpm), "\n")
clin <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
saveRDS(clin, "data/raw/TCGA-LUAD/clinical.rds")
cat("Clinical table:", dim(clin), "\n")
print(table(clin$vital_status, useNA = "ifany"))
exists("query_rna")
exists("se_rna")
exists("clin")
se_rna <- GDCprepare(query_rna, directory = "data/raw/TCGA-LUAD/GDCdata")
exists("se_rna")
class(se_rna)
dim(se_rna)
exists("se_rna")
saveRDS(se_rna, "data/raw/TCGA-LUAD/se_rna.rds")
tpm <- assay(se_rna, "tpm_unstrand")
dim(tpm)
file.exists("data/raw/TCGA-LUAD/se_rna.rds")
colnames(rowData(se_rna))
head(rowData(se_rna), 3)
clin <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
exists("clin")
dim(clin)
saveRDS(clin, "data/raw/TCGA-LUAD/clinical.rds")
file.exists("data/raw/TCGA-LUAD/clinical.rds")
file.info("data/raw/TCGA-LUAD/clinical.rds")$size
dim(clin)
table(clin$vital_status, useNA = "ifany")
table(rowData(se_rna)$gene_type)
all_genes <- rowData(se_rna)$gene_name
sig_genes <- c("GSEC", "FAM83A-AS1", "AL606489.1", "AL034397.3", "AC010980.2",
               "LINC01843", "MIR193BHG", "AC026355.1", "AC124045.1", "AC091185.1",
               "AC027031.2", "AL021707.2", "AL031667.3", "AL606834.1", "AC025048.4",
               "LINC01671", "LINC01116", "LINC01214", "PTCSC3", "LINC02555",
               "ARHGEF26-AS1", "LINC01137", "C20orf197", "MGC32805", "TMPO-AS1",
               "LINC00324", "AL138778.1", "AP000302.1", "AC007384.1", "AL161431.1",
               "KIAA1671-AS1", "AP000695.2", "AC026355.2", "ITGB1-DT", "LINC01150",
               "AL590226.1", "LINC01133", "PLAC4", "MIR223HG", "CROCCP2",
               "KIAA0125", "LINC00996", "RPARP-AS1", "TBX5-AS1", "LINC00857",
               "RP11-21L23.2", "RP11-295G20.2", "AC018629.1", "AC122134.1",
               "AC119424.1", "AL138789.1", "TRMT2B-AS1", "LINC01480", "FRMD6-AS1",
               "MED4-AS1", "LINC01415", "AC092794.1", "AC069023.1", "AP000695.1",
               "AC091057.1", "HLA-DQB1-AS1", "HSPC324", "AC090559.1", "COLCA1",
               "U91328.1", "AL133445.2", "AL442125.1", "AC091132.2", "AC090948.1",
               "AC020765.2", "CASC8", "LINC00707", "OGFRP1", "GASAL1",
               "AC012085.2", "LINC00941", "AP001453.2", "AC026368.1", "AC236972.3")
present <- sig_genes[sig_genes %in% all_genes]
missing <- sig_genes[!sig_genes %in% all_genes]
cat("Present:", length(present), "/", length(sig_genes), "\n\n")
cat("MISSING GENES:\n")
print(missing)
missing_search <- function(gene, all_genes) {
  cat("\n========== ", gene, " ==========\n", sep = "")
  all_genes <- rowData(se_rna)$gene_name
  cat("\n===== AL034397 =====\n")
  print(grep("^AL034397", all_genes, value = TRUE, ignore.case = TRUE))  
  cat("\n===== AC010980 =====\n")
  print(grep("^AC010980", all_genes, value = TRUE, ignore.case = TRUE))  
  cat("\n===== KIAA0125 =====\n")
  print(grep("^KIAA0125", all_genes, value = TRUE, ignore.case = TRUE))  
  cat("\n===== RP11-21L23 =====\n")
  print(grep("^RP11-21L23", all_genes, value = TRUE, ignore.case = TRUE))  
  cat("\n===== RP11-295G20 =====\n")
  print(grep("^RP11-295G20", all_genes, value = TRUE, ignore.case = TRUE))  
  cat("\n===== HSPC324 =====\n")
  print(grep("^HSPC324", all_genes, value = TRUE, ignore.case = TRUE))  
  all_genes <- rowData(se_rna)$gene_name
  cat("\n===== KIAA0125 aliases =====\n")
  print(grep("FAM30A|MIR325HG|KIAA0125", all_genes, value = TRUE, ignore.case = TRUE))  
  cat("\n===== HSPC324 aliases =====\n")
  print(grep("HSPC324|LINC01187|SPRY4-AS1", all_genes, value = TRUE, ignore.case = TRUE))  
  cat("\n===== RP11-21L23 =====\n")
  print(grep("21L23", all_genes, value = TRUE, ignore.case = TRUE))  
  cat("\n===== RP11-295G20 =====\n")
  print(grep("295G20", all_genes, value = TRUE, ignore.case = TRUE))  
  gene_ids   <- rowData(se_rna)$gene_id
  gene_names <- rowData(se_rna)$gene_name
  ensembl_targets <- c(
    "HSPC324"        = "ENSG00000228401",   # Paper #165, #93
    "KIAA0125"       = "ENSG00000272646",   # Paper #305
    "RP11-21L23.2"   = "ENSG00000272650",   # Paper #305
    "RP11-295G20.2"  = "ENSG00000272651"    # Paper #305
  )  
  for (g in names(ensembl_targets)) {
    eid <- ensembl_targets[[g]]
    cat("\n===== ", g, " (", eid, ") =====\n", sep = "")
    idx <- grep(paste0("^", eid), gene_ids)
    if (length(idx) > 0) {
      cat("Found as gene_name:", paste(gene_names[idx], collapse = ", "), "\n")
    } else {
      cat("Not found in current GENCODE\n")
    }
  }  
  exists("gene_ids")
  exists("gene_names")
  exists("se_rna")
  gene_ids <- rowData(se_rna)$gene_id
  length(gene_ids)  
  length(gene_names)  
  idx <- grep("^ENSG00000228401", gene_ids)
  gene_names[idx]  
  idx <- grep("^ENSG00000228401", gene_ids)  
  gene_names[idx]  
  idx <- grep("^ENSG00000272646", gene_ids); gene_names[idx]
  idx <- which(gene_names == "FAM30A"); gene_ids[idx]
  idx <- which(gene_names == "MIR325HG"); gene_ids[idx]
  idx <- grep("^ENSG00000272650", gene_ids); gene_names[idx]
  idx <- grep("^ENSG00000272651", gene_ids); gene_names[idx]  
  which(gene_names == "LINC01187")
  gene_ids[38510]
  which(gene_names == "SPRY4-AS1")
  gene_ids[29955]
  table(sig_genes %in% gene_names)
  colnames(clin)
  grep("vital|days_to|age|gender|sex|stage|pathologic|follow|death", 
       colnames(clin), value = TRUE)
  "submitter_id" %in% colnames(clin)
  "bcr_patient_barcode" %in% colnames(clin)  
  sigs <- read.csv("data/processed/all_signatures_combined.csv")
  names(sigs)[names(sigs) == "paper_id"] <- "record_no"
  write.csv(sigs, "data/processed/signature_coefficients.csv", row.names = FALSE)
  sigs <- read.csv("data/processed/all_signatures_combined.csv")
  head(sigs)
  nrow(sigs)
  colnames(sigs)
  names(sigs)[names(sigs) == "paper_id"] <- "record_no"
  colnames(sigs)  
  write.csv(sigs, "data/processed/signature_coefficients.csv", row.names = FALSE)  
  file.exists("data/processed/signature_coefficients.csv")  
  surv <- readRDS("data/processed/tcga_luad_survival.rds")
  colnames(surv)
  surv$time_months <- surv$OS_time / 30.4375  
  surv$event <- surv$OS  
  table(surv$event)  
  saveRDS(surv, "data/processed/tcga_luad_survival.rds")  
  file.exists("data/processed/tcga_luad_survival.rds")  
  head(sigs)  
  colnames(sigs)  
  table(surv$event)  
  sigs <- read.csv("data/processed/all_signatures_combined.csv")  
  head(sigs)
  nrow(sigs)
  colnames(sigs)
  names(sigs)[names(sigs) == "paper_id"] <- "record_no"
  colnames(sigs)
  write.csv(sigs, "data/processed/signature_coefficients.csv", row.names = FALSE)  
  file.exists("data/processed/signature_coefficients.csv")  
  surv <- readRDS("data/processed/tcga_luad_survival.rds")
  colnames(surv)
  surv$time_months <- surv$OS_time / 30.4375
  surv$event <- surv$OS
  table(surv$event)
  saveRDS(surv, "data/processed/tcga_luad_survival.rds") 
  sigs <- read.csv("data/processed/all_signatures_combined.csv")
  head(sigs)
  nrow(sigs)
  colnames(sigs)
  names(sigs)[names(sigs) == "paper_id"] <- "record_no"
  colnames(sigs)
  write.csv(sigs, "data/processed/signature_coefficients.csv", row.names = FALSE)  
  file.exists("data/processed/signature_coefficients.csv")  
  expr <- readRDS("data/processed/tcga_luad_expr_clean.rds")
  gn <- rownames(expr)
  expr <- readRDS("data/processed/tcga_luad_expr_clean.rds") 
  dim(expr)  
  gn <- rownames(expr)  
  missing_genes <- c("LINC00941", "AL590226.1", "CROCCP2",
                     "RP11-295G20.2", "AC026355.1", "AC018629.1")  
  for (g in missing_genes) {
    hits <- grep(paste0("^", g), gn, value = TRUE)
    cat(sprintf("%-18s → %s\n", g, 
                ifelse(length(hits) > 0, paste(hits, collapse = ", "), "NOT FOUND")))
  }  
  cat("\n--- Variants ---\n")
  cat("AC026355.x:", grep("^AC026355", gn, value = TRUE), "\n")
  cat("LINC0094x:", grep("^LINC0094", gn, value = TRUE), "\n")
  cat("CROCCPx:", grep("^CROCCP", gn, value = TRUE), "\n")
  cat("AC018629.x:", grep("^AC018629", gn, value = TRUE), "\n")
  cat("RP11-295x:", grep("^RP11-295", gn, value = TRUE), "\n")
  sigs <- read.csv("data/processed/signature_coefficients.csv")
  sigs[sigs$record_no == "Paper151", ]  
  p151 <- c("LINC00941", "AP001453.2", "AC026368.1", "AC236972.3")
  p151 %in% gn  
  expr <- readRDS("data/processed/tcga_luad_expr_clean.rds")
  sigs <- read.csv("data/processed/signature_coefficients.csv")
  gn <- rownames(expr)
  sig <- sigs[sigs$record_no == "Paper151", ]
  sig <- sig[!duplicated(sig$gene), ]
  present <- sig$gene %in% gn
  cat("Signature genes:\n")
  print(sig$gene)
  cat("\nPresent in expr:\n")
  print(present)
  cat("\nSum present:", sum(present), "/", nrow(sig), "\n")
  cat("Missing genes:", paste(sig$gene[!present], collapse = ", "), "\n") 
  expr <- readRDS("data/processed/tcga_luad_expr_clean.rds")
  x <- expr["LINC00941", ]
  cat("LINC00941 stats:\n")
  cat("  min:", min(x, na.rm = TRUE), "\n")
  cat("  max:", max(x, na.rm = TRUE), "\n")
  cat("  sd:", sd(x, na.rm = TRUE), "\n")
  cat("  n unique values:", length(unique(x)), "\n")
  cat("  n zeros:", sum(x == 0, na.rm = TRUE), "\n")
  x_log <- log2(x + 1)
  cat("\nAfter log2(x+1):\n")
  cat("  sd:", sd(x_log, na.rm = TRUE), "\n")
  cat("  n unique values:", length(unique(x_log)), "\n")  
  expr <- readRDS("data/processed/tcga_luad_expr_clean.rds")
  idx <- which(rownames(expr) == "LINC00941")  
  cat("Row index:", idx, "\n")  
  print(expr[idx, 1:10])  
  candidates <- grep("LINC00941", rownames(expr), value = TRUE)  
  cat("\nCandidates in matrix:", candidates, "\n")  
  for (c in candidates) {
    x <- expr[c, ]
    cat(sprintf("%-15s: max=%g, mean=%g\n", c, max(x), mean(x)))
  }  
  str(expr[23089, 1:5])
  class(expr[23089, ])
  head(expr[23089, ], 5)
  se_rna <- readRDS("data/raw/TCGA-LUAD/se_rna.rds")
  gi <- as.data.frame(rowData(se_rna))
  hits <- gi[gi$gene_name == "LINC00941", c("gene_id", "gene_name", "gene_type")]  
  print(hits)  
  tpm <- assay(se_rna, "tpm_unstrand")
  gi <- as.data.frame(rowData(se_rna))
  candidates <- grep("^LINC0094", gi$gene_name, value = TRUE)
  cat("Candidates:", paste(candidates, collapse = ", "), "\n\n") 
  for (c in candidates) {
    idx <- which(gi$gene_name == c)[1]
    x <- tpm[idx, ]
    cat(sprintf("%-15s (%s): max=%g, mean=%g, non-zero=%d/%d\n",
                c, gi$gene_id[idx],
                max(x, na.rm = TRUE), mean(x, na.rm = TRUE),
                sum(x > 0, na.rm = TRUE), length(x)))
  }  
  missing_genes <- c("AL590226.1", "CROCCP2", "RP11-295G20.2",
                     "AC026355.1", "AC018629.1")
  cat("Gene name check:\n")
  for (g in missing_genes) {
    idx <- which(gi$gene_name == g)
    if (length(idx) > 0) {
      x <- tpm[idx[1], ]
      cat(sprintf("%-18s FOUND (%s): max=%g, non-zero=%d\n",
                  g, gi$gene_id[idx[1]],
                  max(x, na.rm = TRUE), sum(x > 0, na.rm = TRUE)))
    } else {
      cat(sprintf("%-18s NOT FOUND\n", g))
    }
  }  
  expr <- readRDS("data/processed/tcga_luad_expr_clean.rds")
  aliases_to_add <- list("HSPC324"       = c("AL590226.1"),
                         "KIAA0125"      = c("AC079766.1"),
                         "RP11-21L23.2"  = c("AC110792.3"),
                         "AL034397.3"    = c("AL034397.1"),
                         "AC010980.2"    = c("AC010980.1"),
                         "AC026355.2"    = c("AC026355.1"),
                         "FAM30A"        = c("KIAA0125"),"MIR325HG"      = c("AC079766.1")
  )  
  if (!"CROCCP2" %in% rownames(expr)) {cat("CROCCP2 not in matrix. Need to add from raw data.\n")
  }  
  added <- 0
  for (current in names(aliases_to_add)) {
    if (current %in% rownames(expr)) {
      for (alias in aliases_to_add[[current]]) {
        if (!(alias %in% rownames(expr))) {
          expr <- rbind(expr, expr[current, , drop = FALSE])
          rownames(expr)[nrow(expr)] <- alias
          added <- added + 1
          cat("Added duplicate row:", current, "→", alias, "\n")
        }
      }
    }
  }  
  cat("\nTotal duplicates added:", added, "\n")
  rm(list = ls())
  gc()  
  message("Loading clean expression...")  
  expr <- readRDS("data/processed/tcga_luad_expr_clean.rds")  
  cat("Before:", nrow(expr), "genes\n")  
  aliases_to_add <- list(
    "HSPC324"       = c("AL590226.1"),
    "KIAA0125"      = c("AC079766.1"),
    "RP11-21L23.2"  = c("AC110792.3"),
    "AL034397.3"    = c("AL034397.1"),
    "AC010980.2"    = c("AC010980.1"),
    "AC026355.2"    = c("AC026355.1"),
    "FAM30A"        = c("KIAA0125"),
    "MIR325HG"      = c("AC079766.1")
  )  
  added <- 0
  for (current in names(aliases_to_add)) {
    if (current %in% rownames(expr)) {
      for (alias in aliases_to_add[[current]]) {
        if (!(alias %in% rownames(expr))) {
          expr <- rbind(expr, expr[current, , drop = FALSE])
          rownames(expr)[nrow(expr)] <- alias
          added <- added + 1
          message("Added: ", current, " → ", alias)
        }
      }
    }
  }  
  