# scripts/23_transportability.R
# Calibration-slope transportability analysis + divergence null test,
# adapted from supervisor's template to the actual project objects.

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
})

set.seed(2026)
B_null <- 200
B_opt  <- 200
B_boot <- 1000

# --- 1. Load project data ------------------------------------------
expr_cptac <- readRDS("data/processed/CPTAC3_LUAD_expr_matched.rds")
clin_cptac <- readRDS("data/processed/CPTAC3_LUAD_clin_matched.rds")

expr_gse <- readRDS("data/processed/GSE31210_expr_matched.rds")
clin_gse <- readRDS("data/processed/GSE31210_clin_matched.rds")

# Combine into a named list matching supervisor's structure
cohorts <- list(
  CPTAC3 = list(
    expr  = expr_cptac,
    time  = as.numeric(clin_cptac$os_time),
    event = as.integer(clin_cptac$os_event)
  ),
  GSE31210 = list(
    expr  = expr_gse,
    time  = as.numeric(clin_gse$os_time),
    event = as.integer(clin_gse$os_event)
  )
)

# (z-scoring disabled; raw expression used)
for (cn in names(cohorts)) {
  e <- cohorts[[cn]]$expr
  # z-scoring disabled; raw expression used
  # cohorts[[cn]]$z <- cohorts[[cn]]$expr[complete.cases(cohorts[[cn]]$z), ]
}

# --- 2. Load signatures and build independent-pairs table ----------
sig_files <- list.files("signatures", pattern = "\\.tsv$", full.names = TRUE)

sigs <- do.call(rbind, lapply(sig_files, function(f) {
  s <- read.delim(f, stringsAsFactors = FALSE)
  s$record_no <- sub("_.*$", "", tools::file_path_sans_ext(basename(f)))
  s
}))

# Define independent pairs: only cohorts that did not develop the signature
# Convention: signatures developed on TCGA are excluded from TCGA pairs.
# Since all 15 signatures are TCGA-developed (per screening), only CPTAC3 and
# GSE31210 are independent pairs.

# Coverage rule: 100% of genes present (primary analysis)
pairs_list <- list()
for (cn in c("CPTAC3", "GSE31210")) {
  expr_cohort <- cohorts[[cn]]$expr
  for (sig_id in unique(sigs$record_no)) {
    s <- sigs[sigs$record_no == sig_id, ]
    n_present <- sum(s$gene %in% rownames(expr_cohort))
    if (n_present / nrow(s) >= 1.0) {
      pairs_list[[length(pairs_list) + 1]] <- data.frame(
        record_no = sig_id, cohort = cn,
        n_genes = nrow(s), n_present = n_present,
        stringsAsFactors = FALSE
      )
    }
  }
}
pairs <- do.call(rbind, pairs_list)
cat("Independent pairs:", nrow(pairs), "\n")
print(pairs)

# --- 3. Helper functions -------------------------------------------
cidx <- function(score, time, event, i = seq_along(time)) {
  unname(concordance(Surv(time[i], event[i]) ~ score[i], reverse = TRUE)$concordance)
}

lp_published <- function(sig_genes, sig_coefs, expr) {
  keep <- sig_genes %in% rownames(expr_cohort)
  as.numeric(crossprod(sig_coefs[keep], expr[sig_genes[keep], , drop = FALSE]))
}

calib <- function(lp, time, event) {
  if (sd(lp, na.rm = TRUE) < 1e-6) return(c(beta = NA_real_, se = NA_real_))
  s <- as.numeric(scale(lp))
  f <- coxph(Surv(time, event) ~ s)
  c(beta = unname(coef(f)), se = sqrt(vcov(f)[1, 1]))
}

ridge_lp <- function(x, t, e) {
  f <- cv.glmnet(x, Surv(t, e), family = 'cox', alpha = 0, nfolds = 5)
  as.numeric(predict(f, x, s = 'lambda.min', type = 'link'))
}

refit_c <- function(sig_genes, co) {
  g <- intersect(sig_genes, rownames(co$expr))
  x <- t(co$expr[g, , drop = FALSE])
  t <- co$time; e <- co$event
  if (ncol(x) < 2 || length(unique(e)) < 2) {
    return(c(refit_apparent = NA, optimism = NA, refit_corrected = NA))
  }
  app <- tryCatch(cidx(ridge_lp(x, t, e), t, e), error = function(e) NA_real_)
  opt <- replicate(B_opt, tryCatch({
    i <- sample(length(t), replace = TRUE)
    f <- cv.glmnet(x[i, ], Surv(t[i], e[i]), family = "cox",
                   alpha = 0, nfolds = 5)
    cidx(as.numeric(predict(f, x[i, ], s = "lambda.min")), t[i], e[i]) -
      cidx(as.numeric(predict(f, x,     s = "lambda.min")), t, e)
  }, error = function(err) NA_real_))
  c(refit_apparent    = app,
    optimism          = mean(opt, na.rm = TRUE),
    refit_corrected   = app - mean(opt, na.rm = TRUE))
}

# --- 4. Per-pair table: raw C, slope, refit C ----------------------
est <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(k) {
  cn   <- pairs$cohort[k]
  co   <- cohorts[[cn]]
  s    <- sigs[sigs$record_no == pairs$record_no[k], ]
  lp   <- lp_published(s$gene, s$coefficient, co$expr)
  raw_c <- cidx(lp, co$time, co$event)
  sl   <- calib(lp, co$time, co$event)
  rc   <- refit_c(s$gene, co)
  data.frame(
    record_no = pairs$record_no[k],
    cohort    = cn,
    n_genes   = nrow(s),
    raw_C     = round(raw_c, 3),
    slope     = round(sl["beta"], 3),
    slope_se  = round(sl["se"], 3),
    HR_per_SD = round(exp(sl["beta"]), 3),
    refit_C   = round(rc["refit_apparent"], 3),
    refit_C_corrected = round(rc["refit_corrected"], 3),
    stringsAsFactors = FALSE
  )
}))

cat("\n=== Per-pair results ===\n")
print(est)

# --- 5. Pooled slope with cluster bootstrap -------------------------
pool_once <- function() {
  idx <- lapply(cohorts, function(co) sample(length(co$time), replace = TRUE))
  vapply(seq_len(nrow(pairs)), function(k) {
    cn <- pairs$cohort[k]
    co <- cohorts[[cn]]
    i  <- idx[[cn]]
    s  <- sigs[sigs$record_no == pairs$record_no[k], ]
    tryCatch(
      calib(lp_published(s$gene, s$coefficient, co$expr[, i, drop = FALSE]),
            co$time[i], co$event[i])[["beta"]],
      error = function(e) NA_real_
    )
  }, numeric(1))
}

boot <- replicate(B_boot, pool_once())
pooled <- rowMeans(boot, na.rm = TRUE)

cat("\n=== Pooled calibration slope ===\n")
cat("Mean slope:", round(mean(est$slope, na.rm = TRUE), 3), "\n")
cat("95% CI:", round(quantile(pooled, c(0.025, 0.975), na.rm = TRUE), 3), "\n")

# Classify by pre-specified rules
delta <- log(1.25)
lo <- quantile(pooled, 0.025, na.rm = TRUE)
hi <- quantile(pooled, 0.975, na.rm = TRUE)
class <- if (lo > 0) "Transportable" else if (hi < 0) "Contradicted" else if (hi < delta) "Informative null" else "Inconclusive"

# --- 6. Divergence null test ----------------------------------------
# For each pair: does the divergence rate (raw C < 0.5 AND Cox C > 0.5)
# exceed what expression-matched random gene sets produce?

# Expression bins per cohort for random-gene matching
bins <- lapply(cohorts, function(co) {
  mu <- rowMeans(co$expr, na.rm = TRUE)
  setNames(cut(mu, quantile(mu, 0:10 / 10), include.lowest = TRUE, labels = FALSE),
           names(mu))
})

div_once <- function(genes, coefs, co) {
  keep <- genes %in% rownames(co$expr)
  if (sum(keep) < 2) return(NA_real_)
  genes <- genes[keep]; coefs <- coefs[keep]
  x <- t(co$expr[genes, , drop = FALSE])
  raw <- cidx(as.numeric(x %*% coefs), co$time, co$event)
  app <- tryCatch(cidx(ridge_lp(x, co$time, co$event), co$time, co$event),
                  error = function(e) NA_real_)
  as.numeric(raw < 0.5 && !is.na(app) && app > 0.5)
}

obs_div  <- numeric(nrow(pairs))
null_div <- numeric(nrow(pairs))

for (k in seq_len(nrow(pairs))) {
  cn <- pairs$cohort[k]
  co <- cohorts[[cn]]
  s  <- sigs[sigs$record_no == pairs$record_no[k], ]
  keep <- s$gene %in% rownames(co$expr)
  s <- s[keep, ]
  obs_div[k] <- div_once(s$gene, s$coefficient, co)
  
  # Build null via expression-matched gene swaps
  pool <- setdiff(names(bins[[cn]]), s$gene)
  b <- bins[[cn]]
  pool_by_bin <- split(pool, b[pool])
  sig_bins <- b[s$gene]
  
  null_div[k] <- mean(replicate(B_null, {
    g <- vapply(as.character(sig_bins),
                function(bin) sample(pool_by_bin[[bin]], 1), character(1))
    div_once(g, s$coefficient, co)
  }, simplify = TRUE), na.rm = TRUE)
}

cat("\n=== Divergence null test ===\n")
cat("Observed divergences:", sum(obs_div, na.rm = TRUE), "\n")
cat("Expected under null:  ", round(sum(null_div, na.rm = TRUE), 1), "\n")
cat("Observed proportion:  ", round(mean(obs_div, na.rm = TRUE), 3), "\n")
cat("Null proportion:      ", round(mean(null_div, na.rm = TRUE), 3), "\n")

# --- 7. Save outputs -------------------------------------------------
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
est$obs_div       <- obs_div
est$null_div_prob <- null_div
write.csv(est, "results/tables/transportability_by_pair.csv", row.names = FALSE)
cat("\nSaved results/tables/transportability_by_pair.csv\n")

