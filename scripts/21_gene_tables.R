# 21_gene_tables.R - Step 3a: gene tables for the annotations the papers used.
# GENCODE v22 = the GDC gene model behind HTSeq-FPKM (most 2018-2021 TCGA papers).
# GENCODE v36 = the GDC STAR-Counts model (read from any STAR file downloaded in Step 2).
#
# Download gencode.v22.annotation.gtf.gz once (GENCODE human release 22, or the GDC
# reference-files page) into data/reference/.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
dir.create("data/reference", recursive = TRUE, showWarnings = FALSE)

v22 <- read_gencode_genes("data/reference/gencode.v22.annotation.gtf.gz")
star_files <- list.files("data/raw/GDC", pattern = "star_gene_counts\\.tsv$",
                         recursive = TRUE, full.names = TRUE)
if (!length(star_files)) stop("Run scripts 10-11 first: v36 is read from a STAR-Counts file")
v36 <- read_star_gene_table(star_files[1])

tables <- list(v22 = v22, v36 = v36)
# Add another release only when a paper states it, for example:
# tables$v37 <- read_gencode_genes("data/reference/gencode.v37.annotation.gtf.gz")
saveRDS(tables, "data/reference/gene_tables.rds")
message("Gene tables: ", paste(names(tables), vapply(tables, nrow, 1L), sep = " = ", collapse = "; "))
