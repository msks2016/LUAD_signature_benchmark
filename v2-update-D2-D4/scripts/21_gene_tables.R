# 21_gene_tables.R - Step 3a: gene tables for the annotations the papers used.
# GENCODE v22 = the GDC gene model behind HTSeq-FPKM (most 2018-2021 TCGA papers).
# GENCODE v36 = the GDC STAR-Counts model (read from any STAR file downloaded in Step 2).
#
# GENCODE v22 is downloaded once from the EBI GENCODE archive into data/reference/ (about
# 40 MB; not committed to git). If the download stops, run the script again.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
dir.create("data/reference", recursive = TRUE, showWarnings = FALSE)

gtf22 <- "data/reference/gencode.v22.annotation.gtf.gz"
if (!file.exists(gtf22)) {
  url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.annotation.gtf.gz"
  message("Downloading GENCODE v22 (about 40 MB) ...")
  part <- paste0(gtf22, ".part")
  old_timeout <- options(timeout = 3600)
  ok <- tryCatch(utils::download.file(url, part, mode = "wb") == 0, error = function(e) FALSE)
  options(old_timeout)
  if (!ok) {
    unlink(part)
    stop("GENCODE v22 download failed: run the script again")
  }
  file.rename(part, gtf22)
}
v22 <- read_gencode_genes(gtf22)
if (nrow(v22) < 58000) stop("GENCODE v22 looks incomplete (", nrow(v22), " genes): delete ", gtf22,
                            " and run the script again")
star_files <- list.files("data/raw/GDC", pattern = "star_gene_counts\\.tsv$",
                         recursive = TRUE, full.names = TRUE)
if (!length(star_files)) stop("Run scripts 10-11 first: v36 is read from a STAR-Counts file")
v36 <- read_star_gene_table(star_files[1])

tables <- list(v22 = v22, v36 = v36)
# Add another release only when a paper states it, for example:
# tables$v37 <- read_gencode_genes("data/reference/gencode.v37.annotation.gtf.gz")
saveRDS(tables, "data/reference/gene_tables.rds")
message("Gene tables: ", paste(names(tables), vapply(tables, nrow, 1L), sep = " = ", collapse = "; "))
