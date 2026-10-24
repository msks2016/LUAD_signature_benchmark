# ============================================================
# 00_setup.R - Environment setup
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

required_packages <- c(
  "tidyverse", "data.table",
  "survival", "survminer",
  "glmnet", "timeROC", "survivalROC",
  "ggplot2", "patchwork", "cowplot", "RColorBrewer",
  "pheatmap", "forestplot", "ggpubr", "ggrepel",
  "pROC", "boot", "broom", "janitor",
  "readxl", "writexl", "openxlsx",
  "here", "fs", "glue", "sessioninfo"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg)
  }
}

cat("=== Checking R packages ===\n")
invisible(lapply(required_packages, install_if_missing))

cat("\n=== Setup complete ===\n")
cat("R version:", R.version.string, "\n")
cat("Packages checked:", length(required_packages), "\n")

dir.create("logs", showWarnings = FALSE, recursive = TRUE)
sink("logs/session_info.txt")
sessionInfo()
sink()
