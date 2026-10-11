# ============================================================
# 09_stage_lrt.R
# Likelihood Ratio Test: stage-only vs stage + signature
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

message("Loading data...")
expr <- readRDS("data/processed/tcga_luad_expr_clean_v3.rds")
surv <- readRDS("data/processed/tcga_luad_survival.rds")
sigs <- read.csv("data/processed/signature_coefficients.csv")

# --- Filter + align ---
expr <- expr[, substr(colnames(expr), 14, 15) == "01", drop = FALSE]
colnames(expr) <- substr(colnames(expr), 1, 12)
expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]

common <- intersect(colnames(expr), surv$patient)
expr <- expr[, common, drop = FALSE]
surv <- surv[match(common, surv$patient), ]

# --- Log2 + z-score ---
lexpr <- log2(expr + 1)
sd_vec <- apply(lexpr, 1, sd, na.rm = TRUE)
lexpr  <- lexpr[sd_vec > 0, , drop = FALSE]
zexpr  <- t(scale(t(lexpr)))

# --- Simplify stage to I-IV ---
simplify_stage <- function(x) {
  x <- as.character(x)
  x[grepl("^Stage I$|^Stage IA$|^Stage IB$", x)] <- "I"
  x[grepl("^Stage II$|^Stage IIA$|^Stage IIB$", x)] <- "II"
  x[grepl("^Stage III$|^Stage IIIA$|^Stage IIIB$", x)] <- "III"
  x[grepl("^Stage IV$", x)] <- "IV"
  x[!x %in% c("I","II","III","IV")] <- NA
  factor(x, levels = c("I","II","III","IV"))
}

surv$stage_simple <- simplify_stage(surv$ajcc_pathologic_stage)

cat("Stage distribution (non-NA):\n")
print(table(surv$stage_simple, useNA = "ifany"))

# Keep patients with stage info
surv_stage <- surv %>% filter(!is.na(stage_simple))
cat("\nPatients with stage:", nrow(surv_stage), "\n")

# Align expression to stage-available patients
zexpr_stage <- zexpr[, colnames(zexpr) %in% surv_stage$patient, drop = FALSE]
surv_stage <- surv_stage %>% filter(patient %in% colnames(zexpr_stage))
zexpr_stage <- zexpr_stage[, surv_stage$patient, drop = FALSE]

cat("Final cohort:", ncol(zexpr_stage), "patients,", sum(surv_stage$event), "events\n\n")

# --- LRT per signature ---
results <- list()

for (rid in unique(sigs$record_no)) {
  sig <- sigs %>% filter(record_no == !!rid)
  sig <- sig[!duplicated(sig$gene), ]
  
  present <- sig$gene %in% rownames(zexpr_stage)
  if (mean(!present) > 0.25) {
    cat(sprintf("%-10s | n=%2d | SKIPPED (missing > 25%%)\n", rid, nrow(sig)))
    next
  }
  
  s <- sig[present, ]
  rs <- as.numeric(crossprod(zexpr_stage[s$gene, , drop = FALSE], s$coefficient))
  z_rs <- as.numeric(scale(rs))
  
  d <- data.frame(
    time = surv_stage$time_months,
    event = surv_stage$event,
    stage = surv_stage$stage_simple,
    z_rs = z_rs
  )
  
  fit_a <- coxph(Surv(time, event) ~ stage, data = d)
  fit_b <- coxph(Surv(time, event) ~ stage + z_rs, data = d)
  
  lrt <- anova(fit_a, fit_b)
  lrt_p <- lrt$`Pr(>|Chi|)`[2]
  lrt_chi <- lrt$Chisq[2]
  
  cox_b <- summary(fit_b)
  hr_rs <- exp(cox_b$coef["z_rs", "coef"])
  hr_low <- exp(cox_b$conf.int["z_rs", "lower .95"])
  hr_high <- exp(cox_b$conf.int["z_rs", "upper .95"])
  p_rs <- cox_b$coef["z_rs", "Pr(>|z|)"]
  
  cidx_a <- fit_a$concordance[1]
  cidx_b <- fit_b$concordance[1]
  
  cat(sprintf("%-10s | n=%2d | C_stage=%.3f | C_stage+sig=%.3f | HR_adj=%.3f (%.3f-%.3f) | LRT chi2=%.2f p=%.4f\n",
              rid, nrow(sig), cidx_a, cidx_b, hr_rs, hr_low, hr_high, lrt_chi, lrt_p))
  
  results[[rid]] <- data.frame(
    record_no = rid,
    n_genes = nrow(sig),
    n_present = sum(present),
    C_index_stage_only = cidx_a,
    C_index_stage_plus_sig = cidx_b,
    delta_C = cidx_b - cidx_a,
    HR_adjusted = hr_rs,
    HR_low = hr_low,
    HR_high = hr_high,
    p_risk_score_in_joint = p_rs,
    LRT_chi2 = lrt_chi,
    LRT_p = lrt_p
  )
}

res_df <- do.call(rbind, results)
res_df <- res_df[order(res_df$LRT_p), ]

write.csv(res_df, "results/tables/stage_adjusted_LRT.csv", row.names = FALSE)

cat("\n=========================================\n")
cat("STAGE-ADJUSTED LRT\n")
cat("=========================================\n")
print(res_df, row.names = FALSE)

cat("\nSaved: results/tables/stage_adjusted_LRT.csv\n")
