# 00_setup.R - run from the PROJECT ROOT (open the RStudio project first).
# Creates the project library with renv and records exact package versions in renv.lock.
# Safe to run again: if R restarts after the renv step, just source this file once more.
# Commit renv.lock to git.

if (!file.exists("config/analysis.yaml")) {
  stop("Run this from the project root (the folder that contains config/, R/ and tests/). ",
       "Open the RStudio project, or setwd() to that folder first.")
}

if (!file.exists("renv/activate.R")) {
  if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
  message("Setting up renv. If R restarts, run source('scripts/00_setup.R') again.")
  renv::init(bare = TRUE, restart = FALSE)
}

renv::install(c("BiocManager", "yaml", "survival", "metafor", "testthat", "clubSandwich",
                "data.table", "ggplot2", "rmarkdown"), prompt = FALSE)
renv::install(c("bioc::TCGAbiolinks", "bioc::SummarizedExperiment", "bioc::GEOquery",
                "bioc::Biobase", "bioc::affy", "bioc::hgu133plus2cdf", "bioc::hgu133plus2.db",
                "bioc::AnnotationDbi", "bioc::biomaRt", "bioc::rtracklayer", "bioc::S4Vectors"),
              prompt = FALSE)
renv::snapshot(prompt = FALSE)

# Then check the functions before touching any data:
#   library(testthat); test_dir("tests/testthat")
