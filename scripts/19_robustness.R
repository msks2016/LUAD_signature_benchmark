# scripts/19_robustness.R
# Phase 8 robustness checks:
#   1. Sensitivity to missing-gene threshold (10%, 20%, 30%)
#   2. Stage-restricted metrics (exclude stage IV; stage I-II only)
#   3. Leave-one-cohort-out meta-analysis for Yao 2021
#   4. Early-death exclusion (remove deaths within 3 months)
# Outputs: results/tables/robustness_*.csv

suppressPackageStartupMessages({
  library(survival)
  library(metafor)
})

set.seed(2026)

# ============================================================
# PART 1 — Missing-gene threshold sensitivity
# ============================================================
cat("========== PART 1: Missing-gene threshold ==========\n\n")

sig_files <- list.files("signatures", pattern = "\\.tsv$", full.names = TRUE)
expr_cptac <- readRDS("data/processed/CPTAC3_LUAD_expr_matched.rds")

# Which signatures are evaluable at each threshold
thresholds <- c(10, 20, 30)

coverage_test <- do.call(rbind, lapply(sig_files, function(f) {
  s <- read.delim(f, stringsAsFactors = FALSE)
  sig_name <- tools::file_path_sans_ext(basename(f))
  pct_present <- 100 * sum(s$gene %in% rownames(expr_cptac)) / nrow(s)
  pct_missing <- 100 - pct_present
  data.frame(
    signature     = sig_name,
    n_genes       = nrow(s),
    pct_present   = round(pct_present, 1),
    pct_missing   = round(pct_missing, 1),
    eval_10pct    = pct_missing <= 10,
    eval_20pct    = pct_missing <= 20,
    eval_30pct    = pct_missing <= 30,
    stringsAsFactors = FALSE
  )
}))

cat("Evaluable signatures by threshold:\n")
cat("  <=10% missing:", sum(coverage_test$eval_10pct), "\n")
cat("  <=20% missing:", sum(coverage_test$eval_20pct), " (primary rule)\n")
cat("  <=30% missing:", sum(coverage_test$eval_30pct), "\n\n")

print(coverage_test)
write.csv(coverage_test, "results/tables/robustness_missing_gene_thresholds.csv",
          row.names = FALSE)

# ============================================================
# PART 2 — Stage-restricted metrics (CPTAC-3)
# ============================================================
cat("\n========== PART 2: Stage-restricted metrics ==========\n\n")

scores <- read.csv("results/tables/CPTAC3_LUAD_signature_scores.csv",
                   stringsAsFactors = FALSE)

scores$stage_simple <- sub("Stage (I{1,3}|IV).*", "\\1",
                           as.character(scores$ajcc_pathologic_stage))
scores$stage_simple[!scores$stage_simple %in% c("I","II","III","IV")] <- NA

compute_metrics <- function(d, label) {
  d <- d[!is.na(d$os_time) & d$os_time > 0 & !is.na(d$os_event), ]
  if (nrow(d) < 30 || sum(d$os_event) < 10) return(NULL)
  d$risk_z <- as.numeric(scale(d$risk_score))
  fit <- coxph(Surv(os_time, os_event) ~ risk_z, data = d)
  hr  <- exp(coef(fit))
  ci  <- exp(confint(fit))
  c_cox <- concordance(fit)$concordance
  c_raw <- concordance(Surv(os_time, os_event) ~ risk_score, data = d, reverse = TRUE)$concordance
  data.frame(
    subset = label,
    n = nrow(d),
    events = sum(d$os_event),
    C_raw = round(c_raw, 3),
    C_cox = round(c_cox, 3),
    HR_SD = round(hr, 3),
    HR_lo = round(ci[1], 3),
    HR_hi = round(ci[2], 3),
    stringsAsFactors = FALSE
  )
}

stage_results <- do.call(rbind, lapply(unique(scores$signature), function(s) {
  d <- scores[scores$signature == s, ]
  r_all <- compute_metrics(d, "all")
  r_noIV <- compute_metrics(d[!is.na(d$stage_simple) & d$stage_simple != "IV", ],
                            "exclude_stage_IV")
  r_I_II <- compute_metrics(d[!is.na(d$stage_simple) & d$stage_simple %in% c("I","II"), ],
                            "stage_I_II_only")
  out <- rbind(r_all, r_noIV, r_I_II)
  if (!is.null(out)) out$signature <- s
  out
}))

print(head(stage_results, 12))
write.csv(stage_results, "results/tables/robustness_stage_restricted.csv",
          row.names = FALSE)

# ============================================================
# PART 3 — Leave-one-cohort-out meta-analysis (Yao 2021)
# ============================================================
cat("\n========== PART 3: Leave-one-cohort-out Yao meta-analysis ==========\n\n")

yao_hr <- data.frame(
  cohort = c("TCGA-LUAD", "CPTAC-3", "GSE31210"),
  HR     = c(1.800583, 1.578, 1.817),
  HR_lo  = c(1.557857, 1.194, 1.323),
  HR_hi  = c(2.081128, 2.085, 2.494),
  stringsAsFactors = FALSE
)
yao_hr$se <- (log(yao_hr$HR_hi) - log(yao_hr$HR_lo)) / (2 * 1.96)

cat("Full meta-analysis (all 3 cohorts):\n")
res_full <- rma(yi = log(yao_hr$HR), sei = yao_hr$se, method = "REML")
cat("  Pooled HR:", round(exp(res_full$b), 3),
    " (95% CI", round(exp(res_full$ci.lb), 3), "-",
    round(exp(res_full$ci.ub), 3), ")\n")
cat("  I^2:", round(res_full$I2, 1), "%\n\n")

loo_results <- do.call(rbind, lapply(seq_len(nrow(yao_hr)), function(i) {
  sub <- yao_hr[-i, ]
  r <- rma(yi = log(sub$HR), sei = sub$se, method = "REML")
  data.frame(
    left_out = yao_hr$cohort[i],
    k = nrow(sub),
    HR_pooled = round(exp(r$b), 3),
    HR_lo = round(exp(r$ci.lb), 3),
    HR_hi = round(exp(r$ci.ub), 3),
    I2 = round(r$I2, 1),
    stringsAsFactors = FALSE
  )
}))

print(loo_results)
write.csv(loo_results, "results/tables/robustness_leave_one_cohort_out.csv",
          row.names = FALSE)

# ============================================================
# PART 4 — Early-death exclusion (CPTAC-3)
# ============================================================
cat("\n========== PART 4: Early-death exclusion ==========\n\n")

early_death_results <- do.call(rbind, lapply(unique(scores$signature), function(s) {
  d <- scores[scores$signature == s, ]
  r_all <- compute_metrics(d, "all")
  r_no_3mo <- compute_metrics(d[!(d$os_event == 1 & d$os_time < 3), ],
                              "no_death_within_3mo")
  out <- rbind(r_all, r_no_3mo)
  if (!is.null(out)) out$signature <- s
  out
}))

print(head(early_death_results, 8))
write.csv(early_death_results, "results/tables/robustness_early_death.csv",
          row.names = FALSE)

cat("\n========== All robustness checks complete ==========\n")
cat("Saved:\n")
cat("  results/tables/robustness_missing_gene_thresholds.csv\n")
cat("  results/tables/robustness_stage_restricted.csv\n")
cat("  results/tables/robustness_leave_one_cohort_out.csv\n")
cat("  results/tables/robustness_early_death.csv\n")
