# scripts/10_metrics_CPTAC3.R
# Phase 5 metrics for CPTAC-3: C-index, time-dependent AUC, log-rank, HR
# Output: results/tables/CPTAC3_LUAD_metrics.csv
#         results/figures/CPTAC3_KM_<signature>.pdf

suppressPackageStartupMessages({
  library(survival)
  library(timeROC)
  library(survminer)
  library(dplyr)
})

set.seed(42)

scores <- read.csv("results/tables/CPTAC3_LUAD_signature_scores.csv",
                   stringsAsFactors = FALSE)
cat("Loaded:", nrow(scores), "rows,",
    length(unique(scores$signature)), "signatures\n")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

metrics_list <- list()
sig_names <- unique(scores$signature)

for (sig in sig_names) {
  cat("\n=== ", sig, " ===\n")
  
  d <- scores[scores$signature == sig, ]
  d <- d[!is.na(d$os_time) & d$os_time > 0, ]
  d$risk_z <- as.numeric(scale(d$risk_score))
  
  n         <- nrow(d)
  n_events  <- sum(d$os_event)
  cat("  N =", n, " Events =", n_events, "\n")
  
  # ---- C-index (bootstrap 95% CI) ----
  cfit <- concordance(Surv(os_time, os_event) ~ risk_z, data = d)
  cidx <- cfit$concordance
  
  boot_c <- replicate(1000, {
    idx <- sample(n, replace = TRUE)
    tryCatch({
      cf <- concordance(Surv(os_time, os_event) ~ risk_z, data = d[idx, ])
      cf$concordance
    }, error = function(e) NA_real_)
  })
  c_ci <- quantile(boot_c, c(0.025, 0.975), na.rm = TRUE)
  
  # ---- Time-dependent AUC ----
  auc_1 <- NA; auc_3 <- NA
  auc_1 <- tryCatch({
    t <- timeROC(T = d$os_time, delta = d$os_event,
                 marker = d$risk_z, cause = 1,
                 times = 12, iid = FALSE)
    as.numeric(t$AUC[2])
  }, error = function(e) NA_real_)
  
  auc_3 <- tryCatch({
    t <- timeROC(T = d$os_time, delta = d$os_event,
                 marker = d$risk_z, cause = 1,
                 times = 36, iid = FALSE)
    as.numeric(t$AUC[2])
  }, error = function(e) NA_real_)
  
  # ---- Log-rank (median split) ----
  d$group <- ifelse(d$risk_score > median(d$risk_score), "High", "Low")
  sd_fit  <- survdiff(Surv(os_time, os_event) ~ group, data = d)
  lr_p    <- 1 - pchisq(sd_fit$chisq, df = 1)
  
  # ---- Cox univariable (per SD) ----
  cox1 <- coxph(Surv(os_time, os_event) ~ risk_z, data = d)
  hr_sd <- exp(coef(cox1))
  hr_sd_ci <- exp(confint(cox1))
  
  # ---- Cox univariable (high vs low) ----
  cox2 <- coxph(Surv(os_time, os_event) ~ group, data = d)
  hr_grp <- exp(coef(cox2))
  hr_grp_ci <- exp(confint(cox2))
  
  # ---- Cox adjusted for stage ----
  hr_adj <- NA; hr_adj_lo <- NA; hr_adj_hi <- NA; adj_p <- NA
  if ("ajcc_pathologic_stage" %in% colnames(d) && sum(!is.na(d$ajcc_pathologic_stage)) > 100) {
    d$stage_simple <- sub("Stage (I{1,3}|IV).*", "\\1", d$ajcc_pathologic_stage)
    d$stage_simple[!d$stage_simple %in% c("I","II","III","IV")] <- NA
    cox3 <- tryCatch(
      coxph(Surv(os_time, os_event) ~ risk_z + stage_simple, data = d),
      error = function(e) NULL
    )
    if (!is.null(cox3)) {
      hr_adj    <- exp(coef(cox3))["risk_z"]
      hr_adj_ci <- exp(confint(cox3))["risk_z", ]
      hr_adj_lo <- hr_adj_ci[1]
      hr_adj_hi <- hr_adj_ci[2]
      adj_p     <- summary(cox3)$coefficients["risk_z", "Pr(>|z|)"]
    }
  }
  
  # ---- PH assumption ----
  ph_p <- tryCatch({
    z <- cox.zph(cox1)
    z$table["risk_z", "p"]
  }, error = function(e) NA_real_)
  
  # ---- Save metrics row ----
  metrics_list[[sig]] <- data.frame(
    signature    = sig,
    n_patients   = n,
    n_events     = n_events,
    C_index      = round(cidx, 3),
    C_index_lo   = round(c_ci[1], 3),
    C_index_hi   = round(c_ci[2], 3),
    AUC_1yr      = round(auc_1, 3),
    AUC_3yr      = round(auc_3, 3),
    logrank_p    = signif(lr_p, 3),
    HR_per_SD    = round(hr_sd, 3),
    HR_per_SD_lo = round(hr_sd_ci[1], 3),
    HR_per_SD_hi = round(hr_sd_ci[2], 3),
    HR_high_vs_low    = round(hr_grp, 3),
    HR_high_vs_low_lo = round(hr_grp_ci[1], 3),
    HR_high_vs_low_hi = round(hr_grp_ci[2], 3),
    HR_adj_stage    = round(hr_adj, 3),
    HR_adj_stage_lo = round(hr_adj_lo, 3),
    HR_adj_stage_hi = round(hr_adj_hi, 3),
    adj_p        = signif(adj_p, 3),
    ph_p         = signif(ph_p, 3),
    stringsAsFactors = FALSE
  )
  
  # ---- Kaplan-Meier curve ----
  km_fit <- survfit(Surv(os_time, os_event) ~ group, data = d)
  pdf(file.path("results/figures", paste0("CPTAC3_KM_", sig, ".pdf")),
      width = 6, height = 5)
  print(ggsurvplot(km_fit, data = d,
                   pval = TRUE, risk.table = TRUE,
                   xlab = "Time (months)", ylab = "Overall survival",
                   title = sig, conf.int = FALSE))
  dev.off()
}

metrics <- do.call(rbind, metrics_list)

# ---- Multiple testing correction ----
metrics$logrank_p_BH <- signif(p.adjust(metrics$logrank_p, method = "BH"), 3)

write.csv(metrics, "results/tables/CPTAC3_LUAD_metrics.csv", row.names = FALSE)

cat("\n===================\n")
cat("Metrics saved:", nrow(metrics), "signatures\n")
print(metrics[, c("signature", "n_patients", "n_events",
                  "C_index", "C_index_lo", "C_index_hi",
                  "HR_per_SD", "logrank_p", "logrank_p_BH")])

sig_names_to_test <- c(
  "Paper369_Wu2022",    # 0.335 / 0.665
  "Paper298_Li2023",    # 0.372 / 0.628
  "Paper53_Yang2022"    # 0.385 / 0.615
)
for (s in sig_names_to_test) {
  sig <- read.delim(paste0("signatures/", s, ".tsv"), stringsAsFactors = FALSE)
  sig <- sig[sig$gene %in% rownames(expr), ]
  if (nrow(sig) == 0) { cat(s, " — 0 genes\n"); next }
  k <- nrow(sig)
  coefs <- sig$coefficient  
  nd <- replicate(200, {
    rg <- sample(valid_genes, k)
    sub <- expr[rg, , drop = FALSE]
    sub_z <- t(scale(t(sub)))
    sc <- as.numeric(crossprod(coefs, sub_z))
    df <- data.frame(score = sc, os_time = clin$os_time, os_event = clin$os_event)
    concordance(Surv(os_time, os_event) ~ score, data = df)$concordance
  })  
  cat("\n===", s, "===\n")
  cat("  Median:", round(median(nd), 3), 
      " 2.5%:", round(quantile(nd, 0.025), 3),
      " 97.5%:", round(quantile(nd, 0.975), 3), "\n")
  cat("  Null range:", round(min(nd), 3), "–", round(max(nd), 3), "\n")
}  
