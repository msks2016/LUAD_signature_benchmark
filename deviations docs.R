suppressPackageStartupMessages({ library(tidyverse) })
dir.create("docs", showWarnings = FALSE, recursive = TRUE)

# === 1. Gene mapping ===
gene_mapping <- tribble(
  ~paper_name,     ~paper_ensembl,     ~current_symbol,  ~current_ensembl,     ~gene_type,    ~verified_by,              ~notes,
  "AL034397.3",    "ENSG00000272630",  "AL731563.3",     "ENSG00000272630.1",  "lncRNA",      "Ensembl ID (base match)", "Pre-GENCODE v36 name",
  "AC010980.2",    "ENSG00000272631",  "AC067750.1",     "ENSG00000272631.1",  "lncRNA",      "Ensembl ID (base match)", "Pre-GENCODE v36 name",
  "HSPC324",       "ENSG00000228401",  "AL590226.1",     "ENSG00000228401.5",  "lncRNA",      "Ensembl ID (base match)", "Aliased to paper name",
  "KIAA0125",      "ENSG00000272646",  "AC079766.1",     "ENSG00000272646.1",  "lncRNA",      "Ensembl ID (base match)", "Aliased to paper name",
  "RP11-21L23.2",  "ENSG00000272650",  "AC110792.3",     "ENSG00000272650.1",  "lncRNA",      "Ensembl ID (base match)", "Aliased to paper name",
  "AC018629.1",    "ENSG00000279806",  "AC018629.1",     "ENSG00000279806.1",  "TEC",         "gene_name (TEC)",         "TEC gene_type",
  "LINC00941",     "ENSG00000235884",  "LINC00941",      "ENSG00000235884.4",  "lncRNA",      "gene_name",               "All-zero in GDC STAR counts",
  "RP11-295G20.2", "ENSG00000272651",  NA_character_,    NA_character_,        NA_character_, "not in GENCODE",          "Not in GENCODE v36+"
)
write_csv(gene_mapping, "docs/gene_mapping.csv")
cat("OK: gene_mapping.csv\n")

# === 2. Protocol deviations ===
dev_text <- c(
  "# Protocol Deviations Log",
  "",
  "## Deviation 1: Missing-gene rule relaxation",
  "",
  "Original protocol: 20 percent missing-gene rule.",
  "Deviation applied: 25 percent missing-coefficient-mass rule.",
  "",
  "Reason: Under 20 percent rule, Papers #151 and #305 became not evaluable.",
  "Resolution: Applied both rules. Strict (13/15). Relaxed (15/15).",
  "Both tables retained.",
  "",
  "## Deviation 2: Gene name corrections by Ensembl ID",
  "",
  "Initial gene mappings for Papers #51, #91, #151 used name-based renames.",
  "Corrections:",
  "- AL034397.3 (ENSG00000272630) = AL731563.3, NOT AL034397.1",
  "- AC010980.2 (ENSG00000272631) = AC067750.1, NOT AC010980.1",
  "Built v3 matrix. Rescored affected signatures.",
  "Commit: v2.2 (cd8e30a)"
)
writeLines(dev_text, "docs/PROTOCOL_DEVIATIONS.md")
cat("OK: PROTOCOL_DEVIATIONS.md\n")

# === 3. Cohort overlap ===
cohort_overlap <- data.frame(
  paper = c("Paper51", "Paper53", "Paper62", "Paper91", "Paper141",
            "Paper151", "Paper159", "Paper165", "Paper192", "Paper193",
            "Paper298", "Paper305", "Paper352", "Paper369", "Paper9"),
  development = rep("TCGA-LUAD", 15),
  author_validation = c(
    "TCGA internal", "GSE50081 (partial)", "GSE3141 + GSE37745",
    "None", "TCGA internal", "None", "PDAC exploratory",
    "None", "GSE31210", "TCGA internal", "IMvigor210 (bladder)",
    "GSE30219, GSE31210, GSE37745+GSE50081",
    "TCGA internal", "None", "GSE11969"),
  GSE31210 = c("No","No","No","No","No","No","No","No","YES","No","No","YES","No","No","No"),
  GSE50081 = c("No","YES","No","No","No","No","No","No","No","No","No","YES","No","No","No"),
  GSE30219 = c("No","No","No","No","No","No","No","No","No","No","No","YES","No","No","No"),
  notes = c(
    "TCGA-only", "Skip GSE50081 for Paper53", "None overlap",
    "TCGA-only", "TCGA-only", "GEO attempted; insufficient lncRNA",
    "TCGA-only for LUAD", "TCGA-only", "Skip GSE31210",
    "TCGA-only", "None overlap", "Skip all 3 for Paper305",
    "TCGA-only", "TCGA-only", "No overlap")
)
write.csv(cohort_overlap, "docs/cohort_overlap.csv", row.names = FALSE)
cat("OK: cohort_overlap.csv\n")

# === Verify ===
cat("\n=== Verification ===\n")
for (f in c("docs/gene_mapping.csv",
            "docs/PROTOCOL_DEVIATIONS.md",
            "docs/cohort_overlap.csv")) {
  cat(sprintf("  %-35s exists: %s\n", f, file.exists(f)))
}
line1 <- "# Protocol Deviations Log"
line2 <- ""
line3 <- "## Deviation 1: Missing-gene rule relaxation"
line4 <- ""

line5 <- "Original protocol: 20 percent missing-gene rule."
line6 <- "Deviation applied: 25 percent missing-coefficient-mass rule."
line7 <- ""
line8 <- "Reason: Papers #151 and #305 became not evaluable under 20 percent rule."
line9 <- "Resolution: Applied both rules. Strict (13/15). Relaxed (15/15)."
line10 <- ""
line11 <- "## Deviation 2: Gene name corrections by Ensembl ID"
line12 <- ""
line13 <- "Paper AL034397.3 (ENSG00000272630) = current AL731563.3, NOT AL034397.1"
line14 <- "Paper AC010980.2 (ENSG00000272631) = current AC067750.1, NOT AC010980.1"
line15 <- "Built v3 matrix with correct mappings. Rescored affected signatures."
line16 <- "Commit: v2.2 (cd8e30a)"
all_lines <- c(line1, line2, line3, line4, line5, line6, line7, line8,
               line9, line10, line11, line12, line13, line14, line15, line16)
writeLines(all_lines, "docs/PROTOCOL_DEVIATIONS.md")
file.exists("docs/PROTOCOL_DEVIATIONS.md")
papers <- c("Paper51","Paper53","Paper62","Paper91","Paper141",
            "Paper151","Paper159","Paper165","Paper192","Paper193",
            "Paper298","Paper305","Paper352","Paper369","Paper9")
GSE31210 <- c("No","No","No","No","No","No","No","No","YES","No","No","YES","No","No","No")
GSE50081 <- c("No","YES","No","No","No","No","No","No","No","No","No","YES","No","No","No")
GSE30219 <- c("No","No","No","No","No","No","No","No","No","No","No","YES","No","No","No")
cohort_overlap <- data.frame(
  paper = papers,
  development = rep("TCGA-LUAD", 15),
  GSE31210 = GSE31210,
  GSE50081 = GSE50081,
  GSE30219 = GSE30219,
  stringsAsFactors = FALSE
)
write.csv(cohort_overlap, "docs/cohort_overlap.csv", row.names = FALSE)
file.exists("docs/cohort_overlap.csv")
for (f in c("docs/gene_mapping.csv",
            "docs/PROTOCOL_DEVIATIONS.md",
            "docs/cohort_overlap.csv")) {
  cat(sprintf("  %-35s exists: %s\n", f, file.exists(f)))
}
