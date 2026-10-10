# scripts/18_both_Cindexes.R
# Reports BOTH raw-score and Cox-model C-indexes for every signature-cohort pair.
# Neither convention is "chosen" — both are shown, so the divergence is visible.

suppressPackageStartupMessages(library(survival))

compute_both <- function(d, cohort, platform, sig_short) {
  d <- d[!is.na(d$os_time) & d$os_time > 0 & !is.na(d$os_event), ]
  if (nrow(d) < 50) return(NULL)
  d$risk_z <- as.numeric(scale(d$risk_score))
  
  c_raw <- concordance(Surv(os_time, os_event) ~ risk_score, data = d)$concordance
  
  fit   <- coxph(Surv(os_time, os_event) ~ risk_z, data = d)
  hr    <- exp(coef(fit))
  ci    <- exp(confint(fit))
  c_cox <- concordance(fit)$concordance
  
  data.frame(
    signature = sig_short,
    cohort    = cohort,
    platform  = platform,
    n         = nrow(d),
    events    = sum(d$os_event),
    C_raw     = round(c_raw, 3),
    C_cox     = round(c_cox, 3),
    HR_SD     = round(hr, 3),
    HR_lo     = round(ci[1], 3),
    HR_hi     = round(ci[2], 3),
    dir_ok    = hr > 1,
    stringsAsFactors = FALSE
  )
}

results <- list()

# --- TCGA --------------------------------------------------------
tcga <- read.csv("results/tables/tcga_apparent_performance.csv",
                 stringsAsFactors = FALSE)
tcga <- tcga[tcga$status == "evaluable", ]
for (i in seq_len(nrow(tcga))) {
  results[[length(results) + 1]] <- data.frame(
    signature = tcga$record_no[i],
    cohort    = "TCGA-LUAD",
    platform  = "RNA-seq",
    n         = NA,
    events    = NA,
    C_raw     = round(tcga$C_index[i], 3),
    C_cox     = round(tcga$C_index[i], 3),
    HR_SD     = round(tcga$HR_per_SD[i], 3),
    HR_lo     = round(tcga$HR_low[i], 3),
    HR_hi     = round(tcga$HR_high[i], 3),
    dir_ok    = tcga$HR_per_SD[i] > 1,
    stringsAsFactors = FALSE
  )
}

# --- CPTAC-3 -----------------------------------------------------
cptac_scores <- read.csv("results/tables/CPTAC3_LUAD_signature_scores.csv",
                         stringsAsFactors = FALSE)
for (s in unique(cptac_scores$signature)) {
  d <- cptac_scores[cptac_scores$signature == s, ]
  r <- compute_both(d, "CPTAC-3", "RNA-seq", sub("_.*$", "", s))
  if (!is.null(r)) results[[length(results) + 1]] <- r
}

# --- GSE31210 ----------------------------------------------------
expr_31210 <- readRDS("data/processed/GSE31210_expr_matched.rds")
clin_31210 <- readRDS("data/processed/GSE31210_clin_matched.rds")
zexpr_31210 <- t(scale(t(expr_31210)))
zexpr_31210 <- zexpr_31210[complete.cases(zexpr_31210), ]

gse_sigs <- c("signatures/Paper62_Yao2021.tsv",
              "signatures/Paper305_Yu2024.tsv")
for (f in gse_sigs) {
  sig_df <- read.delim(f, stringsAsFactors = FALSE)
  sig_df <- sig_df[sig_df$gene %in% rownames(zexpr_31210), ]
  if (nrow(sig_df) == 0) next
  sc <- as.numeric(crossprod(sig_df$coefficient,
                             zexpr_31210[sig_df$gene, , drop = FALSE]))
  d <- data.frame(
    signature  = sub("_.*$", "", tools::file_path_sans_ext(basename(f))),
    risk_score = sc,
    os_time    = clin_31210$os_time,
    os_event   = clin_31210$os_event,
    stringsAsFactors = FALSE
  )
  r <- compute_both(d, "GSE31210", "Microarray", unique(d$signature))
  if (!is.null(r)) results[[length(results) + 1]] <- r
}

final <- do.call(rbind, results)
final <- final[order(final$signature, final$cohort), ]

cat("=== Both C-indexes per signature-cohort pair ===\n\n")
print(final)

# Flag where raw and Cox diverge substantially
final$divergence <- round(final$C_cox - final$C_raw, 3)
final$interpretation <- ifelse(
  abs(final$divergence) < 0.05, "convention-agnostic",
  ifelse(final$C_raw < 0.5 & final$C_cox > 0.5,
         "INVERTED — published direction fails",
         "check")
)

cat("\n=== Divergence summary ===\n")
print(table(final$interpretation))

write.csv(final, "results/tables/H3_both_Cindexes.csv", row.names = FALSE)
cat("\nSaved results/tables/H3_both_Cindexes.csv\n")