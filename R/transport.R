# transport.R - aggregate transportability: mean log-HR per SD over independent pairs
# with a patient-level cluster bootstrap (fixes bugs R1 and R2 of the old pipeline).

# log-HR per SD of the published score for one pair, optionally on resampled patients.
# The pair's own cohort object is passed in, so no global variable can leak in (bug R1).
# The signature genes are taken before the resampled patients: the same numbers as subsetting
# patients first, without copying the whole expression matrix in every bootstrap replicate
# (which made scripts 33 and 34 run for hours).
pair_log_hr <- function(cohort, genes, coefs, idx = NULL, scale = "zscore") {
  if (is.null(idx)) idx <- seq_along(cohort$time)
  absent <- setdiff(genes, rownames(cohort$expr))
  if (length(absent)) stop("Not in matrix: ", paste(absent, collapse = ", "))
  s <- score_signature(cohort$expr[genes, idx, drop = FALSE], genes, coefs, scale)
  log_hr_per_sd_fast(cohort$time[idx], cohort$event[idx], s)
}

# Decision rule written in the v0.3 deviations log on 11 Oct 2026, never filed on OSF before
# Update 3: an analysis added after registration (reported as not pre-registered). delta = log(1.25).
classify_transport <- function(lo, hi, delta = log(1.25)) {
  if (lo > 0) return("Transportable")
  if (hi < 0) return("Contradicted")
  if (hi < delta) return("Informative null")
  "Inconclusive"
}

# cohorts: named list of cohort objects (expr, time, event)
# pairs:   data.frame with pair_id, signature, cohort
# genes_by_pair: named list (by pair_id) of data.frames with ensembl_id and coefficient
transport_bootstrap <- function(cohorts, pairs, genes_by_pair, B = 2000, scale = "zscore",
                                delta = log(1.25), seed = NULL, progress = TRUE) {
  if (!is.null(seed)) set.seed(seed)
  K <- nrow(pairs)
  one_pair <- function(k, idx = NULL) {
    g <- genes_by_pair[[pairs$pair_id[k]]]
    pair_log_hr(cohorts[[pairs$cohort[k]]], g$ensembl_id, g$coefficient, idx, scale)
  }
  estimate <- vapply(seq_len(K), function(k) one_pair(k), numeric(1))

  # One resample of patients per cohort, shared by every pair in that cohort,
  # so the correlation between pairs that share patients is kept.
  boot <- matrix(NA_real_, nrow = B, ncol = K, dimnames = list(NULL, pairs$pair_id))
  t0 <- Sys.time()
  for (b in seq_len(B)) {
    idx <- lapply(cohorts, function(co) sample.int(length(co$time), replace = TRUE))
    for (k in seq_len(K)) {
      boot[b, k] <- tryCatch(one_pair(k, idx[[pairs$cohort[k]]]), error = function(e) NA_real_)
    }
    if (progress && (b %% 100 == 0 || b == B)) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      message(sprintf("  bootstrap %d of %d (%d pairs): %.1f min so far, about %.1f min left",
                      b, B, K, el, el / b * (B - b)))
    }
  }
  # Rows are replicates: rowMeans gives ONE pooled mean per replicate (bug R2 averaged pairs).
  pooled_boot <- rowMeans(boot)
  ci <- unname(stats::quantile(pooled_boot, c(0.025, 0.975), na.rm = TRUE))
  list(per_pair = data.frame(pairs, log_hr = estimate, stringsAsFactors = FALSE),
       pooled = mean(estimate), ci_low = ci[1], ci_high = ci[2],
       verdict = classify_transport(ci[1], ci[2], delta),
       n_pairs = K, n_cohorts = length(unique(pairs$cohort)),
       failed_replicates = sum(is.na(pooled_boot)), boot = boot)
}

# Optional cross-check: multilevel meta-analysis with cluster-robust (CR2) inference.
# per_pair needs columns log_hr, se, cohort, pair_id.
robust_crosscheck <- function(per_pair) {
  if (!requireNamespace("clubSandwich", quietly = TRUE)) {
    message("Install clubSandwich for the robust cross-check")
    return(NULL)
  }
  fit <- metafor::rma.mv(yi = log_hr, V = se^2, random = ~ 1 | cohort / pair_id,
                         data = per_pair)
  ci <- as.data.frame(clubSandwich::conf_int(fit, vcov = "CR2", cluster = per_pair$cohort))
  df_col <- grep("^df", names(ci), value = TRUE)[1]
  data.frame(log_hr = ci$beta, ci_low = ci$CI_L, ci_high = ci$CI_U, df = ci[[df_col]])
}
