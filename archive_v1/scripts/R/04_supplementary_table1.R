# ============================================================
# 04_supplementary_table1.R - Supplementary Table 1
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(writexl)
})

all_sigs <- read_csv("data/processed/all_signatures_combined.csv",
                     show_col_types = FALSE)

metadata <- tribble(
  ~paper_id, ~first_author, ~year, ~journal, ~doi, ~biological_category, ~data_source, ~expression_scale, ~annotation_version, ~sample_size,
  "Paper51",  "Song J",   2021, "Bioengineered",                 "10.1080/21655979.2021.1972078", "Pyroptosis-related",              "TCGA-LUAD",              "FPKM",         "Not specified",      504,
  "Paper193", "Zheng Z",  2021, "Frontiers in Genetics",         "10.3389/fgene.2021.690509",     "Ferroptosis-related",             "TCGA-LUAD",              "Not specified","GENCODE GRCh38.p12", 443,
  "Paper53",  "Yang L",   2022, "Frontiers in Oncology",         "10.3389/fonc.2022.780631",      "Genomic instability-derived",     "TCGA-LUAD + GSE50081",   "FPKM",         "Not specified",      497,
  "Paper62",  "Yao J",    2021, "Cancer Cell International",     "10.1186/s12935-021-02027-2",    "Ferroptosis and iron-metabolism", "TCGA-LUAD + GEO",        "Not specified","GENCODE",            477,
  "Paper298", "Li L",     2023, "Scientific Reports",            "10.1038/s41598-023-29684-9",    "Cuproptosis-related",             "TCGA-LUAD + IMvigor210", "FPKM",         "Not specified",      507,
  "Paper159", "He C",     2021, "Translational Cancer Research", "10.21037/tcr-20-2827",          "Immune-associated",               "TCGA-LUAD",              "FPKM",         "GENCODE",            499,
  "Paper192", "Geng W",   2021, "Front Cell Dev Biol",           "10.3389/fcell.2021.657667",     "Genome instability-related",      "TCGA-LUAD + GSE31210",   "FPKM",         "GENCODE release 37", 457,
  "Paper305", "Yu L",     2024, "Scientific Reports",            "10.1038/s41598-024-62836-z",    "Endoplasmic reticulum stress",    "TCGA-LUAD + 4 GEO",      "FPKM",         "Not specified",      515,
  "Paper9",   "Zhang Z",  2021, "Open Medicine",                 "10.1515/med-2021-0276",         "Prognostic lncRNA signature",     "TCGA-LUAD + GSE11969",   "RPKM",         "GENCODE v7",         466,
  "Paper369", "Wu J",     2022, "Frontiers in Genetics",         "10.3389/fgene.2022.966896",     "Necroptosis-related",             "TCGA-LUAD",              "Not specified","GRCh38",             490,
  "Paper165", "Jin D",    2020, "BioMed Research International", "10.1155/2020/7929132",          "Immune-related",                  "TCGA-LUAD",              "Not specified","GENCODE",            437,
  "Paper352", "Zhang G",  2025, "Int J General Medicine",        "Not provided",                  "m6A-related",                     "TCGA-LUAD",              "Not specified","Not specified",      526,
  "Paper141", "Yang Z",   2023, "Genetics Research",             "Not provided",                  "Disulfidptosis-related",          "TCGA-LUAD",              "STAR-Counts",  "GRCh38.gtf",         500,
  "Paper91",  "Lin C",    2024, "Scientific Reports",            "10.1038/s41598-024-65327-3",    "Genomic instability-associated",  "TCGA-LUAD",              "Not specified","GENCODE GRCh37",     490,
  "Paper151", "Diao X",   2022, "Frontiers in Genetics",         "10.3389/fgene.2022.833362",     "Necroptosis-related",             "TCGA-LUAD",              "Not specified","GENCODE",            535
)

table1 <- all_sigs %>%
  group_by(paper_id) %>%
  summarise(
    n_lncrnas = n(),
    lncrna_names = paste(gene, collapse = "; "),
    coefficients = paste(coefficient, collapse = "; "),
    direction = paste(ifelse(coefficient > 0, "Risk", "Protective"), collapse = "; "),
    .groups = "drop"
  ) %>%
  left_join(metadata, by = "paper_id") %>%
  select(
    paper_id, first_author, year, journal, doi,
    biological_category, data_source,
    expression_scale, annotation_version, sample_size,
    n_lncrnas, lncrna_names, coefficients, direction
  )

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write_csv(table1, "results/tables/Supplementary_Table_1.csv")
write_xlsx(table1, "results/tables/Supplementary_Table_1.xlsx")

cat("========================================\n")
cat("SUPPLEMENTARY TABLE 1 SAVED\n")
cat("========================================\n")
cat("  - results/tables/Supplementary_Table_1.csv\n")
cat("  - results/tables/Supplementary_Table_1.xlsx\n")
cat("\nRows:", nrow(table1), "\n\n")
print(table1 %>% select(paper_id, first_author, year, n_lncrnas))
