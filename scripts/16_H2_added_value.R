# scripts/16_H2_added_value.R
# H2: added value beyond stage in CPTAC-3.
# Uses columns already in scores (os_time, os_event, ajcc_pathologic_stage).

suppressPackageStartupMessages(library(survival))
set.seed(2026)

scores <- read.csv("results/tables/CPTAC3_LUAD_signature_scores.csv",
                   stringsAsFactors = FALSE)

# Build stage_simple
scores$stage_simple <- sub("Stage (I{1,3}|IV).*", "\\1",
                           as.character(scores$ajcc_pathologic_stage))
scores$stage_simple[!scores$stage_simple %in% c("I","II","III","IV")] <- NA
scores$stage_simple <- factor(scores$stage_simple, levels = c("I","II","III","IV"))

cat("Rows:", nrow(scores), "\n")
cat("Signatures:", length(unique(scores$signature)), "\n")
cat("Stage distribution:\n")
print(table(scores$stage_simple, useNA = "ifany"))
cat("\n")

# Subset to valid rows
scores_h2 <- scores[!is.na(scores$stage_simple) &
                      !is.na(scores$os_time) &
                      scores$os_time > 0 &
                      !is.na(scores$os_event), ]

cat("Rows for H2:", nrow(scores_h2), "\n")
cat("Events:", sum(scores_h2$os_event), "\n\n")

# Stage-only baseline
stage_df  <- scores_h2[scores_h2$signature == unique(scores_h2$signature)[1], ]
stage_fit <- coxph(Surv(os_time, os_event) ~ stage_simple, data = stage_df)
stage_c   <- concordance(stage_fit)$concordance
cat("Stage-only C-index:", round(stage_c, 3), "\n\n")

# Per-signature
h2_results <- do.call(rbind, lapply(unique(scores_h2$signature), function(s) {
  d <- scores_h2[scores_h2$signature == s, ]
  if (nrow(d) < 50 || sum(d$os_event) < 10) return(NULL)
  
  d$risk_z <- as.numeric(scale(d$risk_score))
  
  m1 <- coxph(Surv(os_time, os_event) ~ stage_simple, data = d)
  m2 <- coxph(Surv(os_time, os_event) ~ stage_simple + risk_z, data = d)
  
  lrt   <- anova(m1, m2, test = "LRT")
  p_lrt <- lrt$`Pr(>|Chi|)`[2]
  
  c1 <- concordance(m1)$concordance
  c2 <- concordance(m2)$concordance
  
  hr_adj    <- exp(coef(m2))["risk_z"]
  hr_adj_ci <- exp(confint(m2))["risk_z", ]
  p_adj     <- summary(m2)$coefficients["risk_z", "Pr(>|z|)"]
  
  data.frame(
    signature     = s,
    n             = nrow(d),
    events        = sum(d$os_event),
    C_stage_only  = round(c1, 3),
    C_full        = round(c2, 3),
    delta_C       = round(c2 - c1, 3),
    LRT_chi2      = round(lrt$Chisq[2], 2),
    LRT_p         = p_lrt,
    HR_sig_adj    = round(hr_adj, 3),
    HR_sig_adj_lo = round(hr_adj_ci[1], 3),
    HR_sig_adj_hi = round(hr_adj_ci[2], 3),
    p_sig_adj     = p_adj,
    stringsAsFactors = FALSE
  )
}))

if (is.null(h2_results) || nrow(h2_results) == 0) {
  stop("No signatures passed the filter.")
}

h2_results$LRT_p_BH <- signif(p.adjust(h2_results$LRT_p,     method = "BH"), 3)
h2_results$p_sig_BH <- signif(p.adjust(h2_results$p_sig_adj, method = "BH"), 3)

h2_results$H2_result <- ifelse(
  h2_results$LRT_p_BH < 0.05 & h2_results$delta_C >  0.01, "Adds value",
  ifelse(h2_results$LRT_p_BH < 0.05 & h2_results$delta_C < -0.01, "Reduces (overfit)",
         "No added value"))

cat("=== H2: added value beyond stage (CPTAC-3) ===\n\n")
print(h2_results)

write.csv(h2_results, "results/tables/CPTAC3_H2_added_value.csv",
          row.names = FALSE)

cat("\n=== Summary ===\n")
print(table(h2_results$H2_result))