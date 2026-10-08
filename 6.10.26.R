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
