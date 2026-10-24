# scripts/17_H3_cross_cohort.R
# Cross-cohort synthesis for H3 (generalization).
# - Computes HR per SD per signature per cohort
# - Pools HRs with metafor::rma for signatures with >=3 cohorts
# - Produces a paired HR panel for 2-cohort signatures
# - Forest plot: Figure 3.

suppressPackageStartupMessages({
  library(survival)
  library(metafor)
  library(ggplot2)
})

set.seed(2026)

# ============================================================
# 1. TCGA HRs (already computed)
# ============================================================
tcga <- read.csv("results/tables/tcga_apparent_performance.csv",
                 stringsAsFactors = FALSE)
tcga$signature <- tcga$record_no          # "Paper62"
tcga <- tcga[tcga$status == "evaluable", ]
tcga_hr <- data.frame(
  signature = tcga$signature,
  cohort    = "TCGA-LUAD",
  platform  = "RNA-seq",
  n         = NA,
  events    = NA,
  HR        = tcga$HR_per_SD,
  HR_lo     = tcga$HR_low,
  HR_hi     = tcga$HR_high,
  C_index   = tcga$C_index,
  stringsAsFactors = FALSE
)
cat("TCGA rows:", nrow(tcga_hr), "\n")

# ============================================================
# 2. CPTAC-3 HRs (compute from scores)
# ============================================================
cptac_scores <- read.csv("results/tables/CPTAC3_LUAD_signature_scores.csv",
                         stringsAsFactors = FALSE)

compute_hr <- function(d, cohort, platform) {
  d <- d[!is.na(d$os_time) & d$os_time > 0 & !is.na(d$os_event), ]
  d$risk_z <- as.numeric(scale(d$risk_score))
  fit <- coxph(Surv(os_time, os_event) ~ risk_z, data = d)
  hr  <- exp(coef(fit)); ci <- exp(confint(fit))
  data.frame(
    signature = unique(d$signature),
    cohort    = cohort,
    platform  = platform,
    n         = nrow(d),
    events    = sum(d$os_event),
    HR        = round(hr, 3),
    HR_lo     = round(ci[1], 3),
    HR_hi     = round(ci[2], 3),
    C_index   = round(concordance(fit)$concordance, 3),
    stringsAsFactors = FALSE
  )
}

cptac_hr <- do.call(rbind, lapply(unique(cptac_scores$signature), function(s) {
  d <- cptac_scores[cptac_scores$signature == s, ]
  if (nrow(d) < 50) return(NULL)
  compute_hr(d, "CPTAC-3", "RNA-seq")
}))

# Convert long-form signature names to short form for join
shorten_sig <- function(x) sub("_.*$", "", x)
cptac_hr$signature <- shorten_sig(cptac_hr$signature)

cat("CPTAC-3 rows:", nrow(cptac_hr), "\n")

# ============================================================
# 3. GSE31210 HRs (score + compute)
# ============================================================
expr_31210 <- readRDS("data/processed/GSE31210_expr_matched.rds")
clin_31210 <- readRDS("data/processed/GSE31210_clin_matched.rds")
stopifnot(all(colnames(expr_31210) == clin_31210$sample_id))

zexpr_31210 <- t(scale(t(expr_31210)))
zexpr_31210 <- zexpr_31210[complete.cases(zexpr_31210), ]

# Only the two signatures evaluable in GSE31210
gse_sig_files <- c("signatures/Paper62_Yao2021.tsv",
                   "signatures/Paper305_Yu2024.tsv")

gse_hr <- do.call(rbind, lapply(gse_sig_files, function(f) {
  sig <- read.delim(f, stringsAsFactors = FALSE)
  sig <- sig[sig$gene %in% rownames(zexpr_31210), ]
  if (nrow(sig) == 0) return(NULL)
  
  sc  <- as.numeric(crossprod(sig$coefficient,
                              zexpr_31210[sig$gene, , drop = FALSE]))
  d <- data.frame(
    signature  = sub("_.*$", "", tools::file_path_sans_ext(basename(f))),
    risk_score = sc,
    os_time    = clin_31210$os_time,
    os_event   = clin_31210$os_event,
    stringsAsFactors = FALSE
  )
  compute_hr(d, "GSE31210", "Microarray")
}))

cat("GSE31210 rows:", nrow(gse_hr), "\n\n")

# ============================================================
# 4. Combine into cross-cohort table
# ============================================================
all_hr <- rbind(tcga_hr, cptac_hr, gse_hr)

# Report coverage per signature
coverage_table <- table(all_hr$signature)
cat("Cohorts per signature:\n")
print(coverage_table)

write.csv(all_hr, "results/tables/H3_cross_cohort_HRs.csv",
          row.names = FALSE)

# ============================================================
# 5. Forest plot for Yao 2021 (3 cohorts)
# ============================================================
yao <- all_hr[all_hr$signature == "Paper62", ]

if (nrow(yao) >= 2) {
  yao$se <- (log(yao$HR_hi) - log(yao$HR_lo)) / (2 * 1.96)
  
  res <- rma(yi = log(yao$HR), sei = yao$se,
             method = "REML", slab = yao$cohort)
  
  cat("\n=== Yao 2021 meta-analysis ===\n")
  print(res)
  
  # Forest plot
  pdf("results/figures/Figure3_Yao2021_forest.pdf", width = 8, height = 4.5)
  forest(res,
         atransf = exp,
         header = c("Cohort", "HR per SD [95% CI]"),
         xlab = "Hazard ratio per SD of risk score",
         mlab = paste0("Random-effects model (I\u00B2 = ",
                       round(res$I2, 1), "%)"),
         col = "#3A7CA5",
         border = "#3A7CA5")
  dev.off()
  
  # PNG version too
  png("results/figures/Figure3_Yao2021_forest.png",
      width = 8, height = 4.5, units = "in", res = 300)
  forest(res,
         atransf = exp,
         header = c("Cohort", "HR per SD [95% CI]"),
         xlab = "Hazard ratio per SD of risk score",
         mlab = paste0("Random-effects model (I\u00B2 = ",
                       round(res$I2, 1), "%)"),
         col = "#3A7CA5",
         border = "#3A7CA5")
  dev.off()
  
  cat("Saved Figure3_Yao2021_forest.pdf and .png\n")
}

# ============================================================
# 6. Paired HR panel for 2-cohort signatures
# ============================================================
paired_sigs <- names(coverage_table[coverage_table == 2])
paired_df <- all_hr[all_hr$signature %in% paired_sigs, ]

if (nrow(paired_df) > 0) {
  paired_df$cohort <- factor(paired_df$cohort,
                             levels = c("TCGA-LUAD", "CPTAC-3", "GSE31210"))
  
  p <- ggplot(paired_df, aes(x = cohort, y = HR, colour = cohort)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = HR_lo, ymax = HR_hi), width = 0.15) +
    facet_wrap(~ signature, scales = "free_y") +
    scale_colour_manual(values = c("TCGA-LUAD" = "#3A7CA5",
                                   "CPTAC-3" = "#D7263D",
                                   "GSE31210" = "#F4A259")) +
    labs(
      x = NULL, y = "HR per SD (95% CI)",
      title = "Cross-cohort HR per SD by signature",
      subtitle = "Paired comparison of apparent vs independent validation"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  ggsave("results/figures/Figure3b_paired_HRs.pdf", p,
         width = 10, height = 6)
  ggsave("results/figures/Figure3b_paired_HRs.png", p,
         width = 10, height = 6, dpi = 300)
  cat("Saved Figure3b_paired_HRs.pdf and .png\n")
}

# ============================================================
# 7. Direction consistency summary
# ============================================================
all_hr$dir_pub <- "increases risk"     # all published signatures are risk-increasing
all_hr$dir_obs <- ifelse(all_hr$HR > 1, "increases risk",
                         ifelse(all_hr$HR < 1, "decreases risk", "null"))
all_hr$consistent <- all_hr$dir_obs == all_hr$dir_pub

cat("\n=== Direction consistency ===\n")
print(table(all_hr$cohort, all_hr$consistent))

write.csv(all_hr,
          "results/tables/H3_cross_cohort_HRs_with_direction.csv",
          row.names = FALSE)

cat("\nDone. Outputs:\n")
cat("  results/tables/H3_cross_cohort_HRs.csv\n")
cat("  results/tables/H3_cross_cohort_HRs_with_direction.csv\n")
cat("  results/figures/Figure3_Yao2021_forest.pdf\n")
cat("  results/figures/Figure3b_paired_HRs.pdf\n")