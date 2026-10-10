# scripts/12_paired_TCGA_CPTAC3.R
# Paired analysis: TCGA-LUAD apparent performance vs CPTAC-3 independent validation

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# --- 1. Load both tables --------------------------------------------
tcga <- read.csv("results/tables/tcga_apparent_performance.csv",
                 stringsAsFactors = FALSE)
cptac <- read.csv("results/tables/CPTAC3_signature_null_matched.csv",
                  stringsAsFactors = FALSE)

# --- 2. Harmonize signature IDs -------------------------------------
# TCGA uses "Paper62"; CPTAC-3 uses "Paper62_Yao2021"
tcga$paper_no   <- tcga$record_no                                   # already "Paper62"
cptac$paper_no  <- sub("_.*$", "", cptac$signature)                 # strip "_Yao2021"

# --- 3. Join --------------------------------------------------------
paired <- merge(
  tcga[, c("paper_no", "C_index", "HR_per_SD", "logrank_p", "status")],
  cptac[, c("paper_no", "C_observed", "C_lo", "C_hi",
            "null_2.5", "null_97.5", "p_holm", "class")],
  by = "paper_no"
)
colnames(paired)[colnames(paired) == "C_index"]    <- "C_TCGA"
colnames(paired)[colnames(paired) == "C_observed"] <- "C_CPTAC3"

cat("Paired signatures:", nrow(paired), "\n")
print(paired)

# --- 4. Save --------------------------------------------------------
write.csv(paired, "results/tables/paired_TCGA_CPTAC3.csv", row.names = FALSE)

# --- 5. Plot --------------------------------------------------------
paired$class_color <- ifelse(paired$class == "Contradicted",
                             "Contradicted",
                             "Inconclusive")

p <- ggplot(paired, aes(x = C_TCGA, y = C_CPTAC3)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey70") +
  geom_point(aes(colour = class_color, size = n_genes)) +
  geom_errorbar(aes(ymin = C_lo, ymax = C_hi), width = 0.005, alpha = 0.5) +
  scale_colour_manual(values = c("Contradicted" = "#D7263D",
                                 "Inconclusive" = "#3A7CA5")) +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    x = "TCGA-LUAD apparent C-index",
    y = "CPTAC-3 independent C-index",
    title = "Published LUAD lncRNA signatures fail independent validation",
    subtitle = paste0("13 signatures with paired C-indexes; ",
                      "all drop by 0.11 to 0.34"),
    colour = "CPTAC-3 classification",
    size   = "Signature size"
  ) +
  theme_minimal(base_size = 12) +
  coord_cartesian(xlim = c(0.5, 0.75), ylim = c(0.2, 0.65))

ggsave("results/figures/Figure2_TCGA_vs_CPTAC3.pdf", p,
       width = 7, height = 5.5)
ggsave("results/figures/Figure2_TCGA_vs_CPTAC3.png", p,
       width = 7, height = 5.5, dpi = 300)

cat("\nSaved:\n")
cat("  results/tables/paired_TCGA_CPTAC3.csv\n")
cat("  results/figures/Figure2_TCGA_vs_CPTAC3.pdf\n")
cat("  results/figures/Figure2_TCGA_vs_CPTAC3.png\n")