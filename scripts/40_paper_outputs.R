# 40_paper_outputs.R - manuscript tables and figures, built only from the v2 outputs of scripts
# 01-37 (results/tables, results/numbers.csv, data/signatures, config) and the cohort objects.
# Nothing is recomputed except descriptive cohort statistics (median follow-up, stage counts)
# and the Kaplan-Meier curves. Run after scripts 29-37.
# Output: results/paper/tables/*.csv, results/paper/figures/*.pdf|png, results/paper/README.md

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
out_t <- "results/paper/tables"
out_f <- "results/paper/figures"
dir.create(out_t, recursive = TRUE, showWarnings = FALSE)
dir.create(out_f, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) utils::read.csv(file.path("results/tables", f), stringsAsFactors = FALSE, check.names = FALSE)
wr <- function(x, f) utils::write.csv(x, file.path(out_t, f), row.names = FALSE, na = "")

sig <- utils::read.csv("data/signatures/signature_master.csv", stringsAsFactors = FALSE)
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
ovr <- utils::read.csv("config/gene_overrides.csv", stringsAsFactors = FALSE)
alts <- read_gene_alternatives()
nums <- utils::read.csv("results/numbers.csv", stringsAsFactors = FALSE)
num <- function(k) nums$value[match(k, nums$key)]
cov <- rd("coverage_long.csv")
pc <- rd("positive_control.csv")
h1 <- rd("h1_null_test.csv")
h2 <- rd("h2_added_value.csv")
h3 <- rd("h3_per_pair.csv")
pooled <- rd("h3_pooled_exploratory.csv")
loo <- rd("h3_leave_one_out.csv")
psens <- rd("h3_pooled_sensitivity.csv")
tr <- rd("transportability_pooled.csv")
rc <- if (file.exists("results/tables/transportability_robust_crosscheck.csv")) rd("transportability_robust_crosscheck.csv") else NULL
str_ <- rd("sensitivity_transportability.csv")
sh1 <- rd("sensitivity_h1.csv")
dcs <- rd("direction_consistency.csv")
pw <- rd("power_simulation.csv")
cohort_order <- c(cfg$cohorts$reference, cfg$cohorts$independent)
sig_ids <- sort(unique(sig$signature))
sig_ids <- sig_ids[order(sig_label(sig_ids))]
cols <- cohort_colours()

# ---- Table: signatures -------------------------------------------------------------------
indep_pairs <- h1[, c("signature", "cohort")]
t_sig <- do.call(rbind, lapply(sig_ids, function(s) {
  d <- sig[sig$signature == s, ]
  used <- ind$cohort[ind$signature == s & ind$used_by_authors == "yes"]
  data.frame(Signature = sig_label(s), ID = s, DOI = d$doi[1], Genes = nrow(d),
             `Risk coefficients (>0)` = sum(d$coefficient > 0),
             `Protective coefficients (<0)` = sum(d$coefficient < 0),
             `Gene-name release` = paste("GENCODE", unique(d$annotation), collapse = "/"),
             `Genes with a documented manual ID` = sum(ovr$signature == s),
             `Cohorts used by the authors` = paste(used, collapse = ", "),
             `Independent cohorts evaluable` = paste(indep_pairs$cohort[indep_pairs$signature == s], collapse = ", "),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
wr(t_sig, "table_signatures.csv")

# ---- Table: cohorts ----------------------------------------------------------------------
cohorts <- load_cohorts(cohort_order)
t_coh <- do.call(rbind, lapply(cohort_order, function(cn) {
  co <- cohorts[[cn]]
  st <- collapse_stage(co$clin$stage)
  cv <- cov[cov$cohort == cn, ]
  p <- pc[pc$cohort == cn, ]
  data.frame(Cohort = cn, Role = if (cn %in% cfg$cohorts$reference) "reference (apparent performance)" else "independent",
             Platform = paste(c(co$meta$platform, co$meta$source), collapse = "; "),
             Patients = length(co$time), Deaths = sum(co$event),
             `Median follow-up, months (reverse KM)` = round(reverse_km_median(co$time, co$event), 1),
             `Stage I` = sum(st %in% "I"), `Stage II` = sum(st %in% "II"), `Stage III` = sum(st %in% "III"),
             `Stage IV` = sum(st %in% "IV"), `Stage missing` = sum(is.na(st)),
             `Stage I-II, %` = round(100 * mean(st[!is.na(st)] %in% c("I", "II")), 1),
             `Positive control` = if (nrow(p)) paste0(ifelse(p$pass, "pass", "fail"), " (",
                                                      ifelse(grepl("^adapted", p$rule), "HR II vs I", "HR III vs I"), ")") else "",
             `Signatures with >= 80% of genes` = sum(cv$frac_present >= cfg$coverage$primary_min_fraction - 1e-9),
             `Signatures with 100% of genes` = sum(cv$frac_present >= 1 - 1e-9),
             `Independent evaluable pairs` = sum(h1$cohort == cn),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
wr(t_coh, "table_cohorts.csv")

# ---- Table: the independent signature x cohort pairs (H1, HR, H2) --------------------------
m <- merge(h1, h3[, c("pair_id", "hr", "hr_low", "hr_high")], by = "pair_id")
m <- merge(m, h2[, c("pair_id", "delta_c", "delta_c_low", "delta_c_high", "p_lrt", "p_holm", "h2_result")],
           by = "pair_id", all.x = TRUE, suffixes = c("", "_h2"))
m <- m[order(match(m$cohort, cohort_order), sig_label(m$signature)), ]
t_pairs <- data.frame(Signature = sig_label(m$signature), Cohort = m$cohort,
                      `Genes present` = paste0(m$n_present, "/", m$n_genes), Patients = m$patients, Deaths = m$deaths,
                      `C-index (95% CI)` = fmt_ci(m$c_index, m$c_low, m$c_high),
                      `Random-gene null, 2.5-97.5%` = sprintf("%.2f-%.2f", m$null_q025, m$null_q975),
                      `H1 p` = fmt_p(m$p_two_sided), `H1 p (Holm)` = fmt_p(m$p_holm), `H1` = m$classification,
                      `HR per SD (95% CI)` = fmt_ci(m$hr, m$hr_low, m$hr_high),
                      `Delta C-index (95% CI)` = fmt_ci(m$delta_c, m$delta_c_low, m$delta_c_high, 3, " to "),
                      `H2 p (Holm)` = fmt_p(m$p_holm_h2), `H2` = m$h2_result,
                      check.names = FALSE, stringsAsFactors = FALSE)
wr(t_pairs, "table_independent_pairs.csv")

# ---- Table: aggregate analyses ------------------------------------------------------------
yao <- pooled[1, ]
v <- function(type) psens[psens$cohorts == "all" & psens$interval == type & psens$signature == yao$signature, ]
t_agg <- rbind(
  data.frame(Analysis = sprintf("Aggregate transportability: mean log-HR per SD over %d independent pairs in %d cohorts, patient-level cluster bootstrap (%d resamples); not pre-registered",
                                tr$n_pairs, tr$n_cohorts, tr$resamples),
             Estimate = sprintf("HR %s; log-HR %s", fmt_ci(exp(tr$pooled_log_hr), exp(tr$ci_low), exp(tr$ci_high)),
                                fmt_ci(tr$pooled_log_hr, tr$ci_low, tr$ci_high, 3, " to ")),
             Result = tr$verdict),
  if (!is.null(rc)) data.frame(Analysis = "Cross-check: multilevel meta-analysis (pairs within cohorts), CR2 cluster-robust CI",
                               Estimate = sprintf("HR %s; df %.1f", fmt_ci(exp(rc$log_hr), exp(rc$ci_low), exp(rc$ci_high)), rc$df),
                               Result = ifelse(rc$ci_low > 0, "CI above HR 1", "CI includes HR 1")),
  data.frame(Analysis = paste("Direction consistency (Update 1):", dcs$level),
             Estimate = sprintf("%d of %d in the published direction (%.0f%%)", dcs$consistent, dcs$n, 100 * dcs$proportion),
             Result = p_text(dcs$p_two_sided, "binomial p")),
  data.frame(Analysis = sprintf("%s pooled over %d independent cohorts (%s; %d deaths), REML + Hartung-Knapp (registered, exploratory)",
                                sig_label(yao$signature), yao$k, gsub(";", ", ", yao$cohorts), yao$deaths),
             Estimate = sprintf("HR %s; I2 %.0f%%; tau2 %.3f", fmt_ci(yao$hr, yao$hr_low, yao$hr_high), yao$i2, yao$tau2),
             Result = ""),
  data.frame(Analysis = paste(sig_label(yao$signature), "pooled: modified Hartung-Knapp (q* = max(1, q)); not pre-registered"),
             Estimate = sprintf("HR %s; q = %.2f", fmt_ci(v("reml_hk_modified")$hr, v("reml_hk_modified")$hr_low, v("reml_hk_modified")$hr_high), v("reml_hk_modified")$q),
             Result = ""),
  data.frame(Analysis = paste(sig_label(yao$signature), "pooled: common effect"),
             Estimate = paste("HR", fmt_ci(v("common_effect")$hr, v("common_effect")$hr_low, v("common_effect")$hr_high)),
             Result = ""),
  data.frame(Analysis = paste(sig_label(yao$signature), "leave-one-cohort-out (REML + Hartung-Knapp)"),
             Estimate = sprintf("HR %.2f to %.2f; lowest CI bound %.2f", min(loo$hr), max(loo$hr), min(loo$hr_low)),
             Result = ""),
  data.frame(Analysis = "Gene level: share of (gene, independent cohort) pairs whose univariable Cox sign matches the published coefficient",
             Estimate = num("gene_sign_agreement_independent"), Result = "")
)
wr(t_agg, "table_aggregate.csv")

# ---- Table: sensitivity analyses ---------------------------------------------------------
t_sens_tr <- data.frame(Setting = str_$setting, Pairs = str_$n_pairs,
                        `HR per SD (95% CI)` = fmt_ci(exp(str_$pooled_log_hr), exp(str_$ci_low), exp(str_$ci_high)),
                        `log-HR (95% CI)` = fmt_ci(str_$pooled_log_hr, str_$ci_low, str_$ci_high, 3, " to "),
                        Verdict = str_$verdict, check.names = FALSE, stringsAsFactors = FALSE)
h1_sum <- do.call(rbind, lapply(unique(sh1$setting), function(s) {
  d <- sh1[sh1$setting == s, ]
  data.frame(Setting = s, Pairs = nrow(d), Replicated = sum(d$classification == "Replicated"),
             Contradicted = sum(d$classification == "Contradicted"), Inconclusive = sum(d$classification == "Inconclusive"),
             `Replicated pairs` = paste(paste(sig_label(d$signature), d$cohort)[d$classification == "Replicated"], collapse = "; "),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
h1_sum <- rbind(data.frame(Setting = "primary (registered)", Pairs = nrow(h1), Replicated = sum(h1$classification == "Replicated"),
                           Contradicted = sum(h1$classification == "Contradicted"), Inconclusive = sum(h1$classification == "Inconclusive"),
                           `Replicated pairs` = paste(paste(sig_label(h1$signature), h1$cohort)[h1$classification == "Replicated"], collapse = "; "),
                           check.names = FALSE, stringsAsFactors = FALSE), h1_sum)
wr(t_sens_tr, "table_sensitivity_transportability.csv")
wr(h1_sum, "table_sensitivity_h1.csv")

# ---- Supplementary tables ------------------------------------------------------------------
wr(data.frame(Cohort = pc$cohort, Role = pc$role, Patients = pc$patients, Deaths = pc$deaths, Stages = pc$n_by_stage,
              `Stage III-IV patients` = pc$n_iii_iv, Rule = pc$rule,
              `HR II vs I` = round(pc$hr_ii, 2), `p II` = fmt_p(pc$p_ii), `HR III vs I` = round(pc$hr_iii, 2),
              `p III` = fmt_p(pc$p_iii), `HR IV vs I` = round(pc$hr_iv, 2), `p IV` = fmt_p(pc$p_iv),
              `Stage C-index` = round(pc$c_stage, 3), Monotonic = pc$monotonic, Pass = pc$pass, Note = pc$note,
              check.names = FALSE), "supp_positive_control.csv")
cw <- stats::reshape(cov[, c("signature", "cohort", "n_present", "n_genes")], idvar = "signature", timevar = "cohort",
                     direction = "wide")
cov_wide <- data.frame(Signature = sig_label(cw$signature), Genes = cw[[paste0("n_genes.", cohort_order[1])]])
for (cn in cohort_order) cov_wide[[cn]] <- paste0(cw[[paste0("n_present.", cn)]], "/", cw[[paste0("n_genes.", cn)]])
cov_wide <- cov_wide[order(cov_wide$Signature), ]
wr(cov_wide, "supp_coverage_matrix.csv")
wr(data.frame(Cohort = cov$cohort, Signature = sig_label(cov$signature), Present = cov$n_present, Genes = cov$n_genes,
              `Fraction present` = round(cov$frac_present, 3), `Missing genes` = cov$missing,
              `Present through a successor ID` = cov$via_alternative_id, check.names = FALSE), "supp_coverage_long.csv")
if (file.exists("results/tables/presence_check.csv")) {
  prc <- rd("presence_check.csv")
  prc <- prc[prc$reason != "present", ]
  wr(data.frame(Cohort = prc$cohort, Signature = sig_label(prc$signature), Gene = prc$published_symbol,
                `Ensembl ID` = prc$ensembl_id, Reason = prc$reason, `In Ensembl 114` = prc$in_ensembl114,
                `Replaced by (Ensembl history)` = prc$replaced_by, check.names = FALSE),
     "supp_missing_genes_reasons.csv")
}
gid <- merge(gm[, intersect(c("signature", "published_symbol", "coefficient", "annotation", "ensembl_id", "id_source", "status"), names(gm))],
             alts[, c("ensembl_id", "alternative_id")], by = "ensembl_id", all.x = TRUE)
gid$override_evidence <- ovr$evidence[match(paste(gid$signature, gid$published_symbol), paste(ovr$signature, ovr$published_symbol))]
gid$successor_evidence <- alts$evidence[match(paste(gid$ensembl_id, gid$alternative_id), paste(alts$ensembl_id, alts$alternative_id))]
gid <- gid[order(sig_label(gid$signature), gid$published_symbol), ]
gid$signature <- sig_label(gid$signature)
gid <- gid[, c(intersect(c("signature", "published_symbol", "coefficient", "annotation", "ensembl_id", "id_source", "status"), names(gid)),
               "alternative_id", "override_evidence", "successor_evidence")]
wr(gid, "supp_gene_identity.csv")
wr(data.frame(Signature = sig_label(h1$signature), Cohort = h1$cohort, Patients = h1$patients, Deaths = h1$deaths,
              `C-index` = round(h1$c_index, 3), `C-index 95% CI` = sprintf("%.3f-%.3f", h1$c_low, h1$c_high),
              `Cox-model C-index` = round(h1$c_index_cox, 3), `Null pool genes` = h1$null_pool_genes,
              `Null 2.5%` = round(h1$null_q025, 3), `Null 97.5%` = round(h1$null_q975, 3), p = fmt_p(h1$p_two_sided),
              `p (Holm)` = fmt_p(h1$p_holm), `Holm family` = h1$holm_family_size, Classification = h1$classification,
              `Cohort status` = h1$cohort_status, check.names = FALSE), "supp_h1_null_test.csv")
wr(data.frame(Signature = sig_label(h2$signature), Cohort = h2$cohort, Patients = h2$n, Deaths = h2$events, Covariates = h2$covariates,
              `C base` = round(h2$c_base, 3), `C with score` = round(h2$c_full, 3),
              `Delta C (95% CI)` = fmt_ci(h2$delta_c, h2$delta_c_low, h2$delta_c_high, 3, " to "),
              `LRT chi2` = round(h2$lrt_chisq, 2), `p (LRT)` = fmt_p(h2$p_lrt), `p (Holm)` = fmt_p(h2$p_holm),
              `Adjusted HR per SD (95% CI)` = fmt_ci(h2$hr_adj, h2$hr_lo, h2$hr_hi), Result = h2$h2_result,
              check.names = FALSE), "supp_h2_added_value.csv")
h3o <- h3[order(h3$role != "reference", match(h3$cohort, cohort_order), sig_label(h3$signature)), ]
wr(data.frame(Signature = sig_label(h3o$signature), Cohort = h3o$cohort, Role = h3o$role, Patients = h3o$patients, Deaths = h3o$deaths,
              `Genes present` = paste0(h3o$n_present, "/", h3o$n_genes),
              `HR per SD (95% CI)` = fmt_ci(h3o$hr, h3o$hr_low, h3o$hr_high),
              `C-index (95% CI)` = fmt_ci(h3o$c_index, h3o$c_low, h3o$c_high, 3), `Cox-model C-index` = round(h3o$c_index_cox, 3),
              `Uno C (60 months)` = round(h3o$c_uno, 3),
              `td-AUC 12 m (95% CI)` = fmt_ci(h3o$td_auc_12m, h3o$td_auc_low_12m, h3o$td_auc_high_12m, 3),
              `td-AUC 36 m (95% CI)` = fmt_ci(h3o$td_auc_36m, h3o$td_auc_low_36m, h3o$td_auc_high_36m, 3),
              `Log-rank p (median split)` = fmt_p(h3o$logrank_p), `PH test p` = fmt_p(h3o$ph_p),
              `RMST tau, months` = round(h3o$rmst_tau, 1),
              `RMST difference high - low, months (95% CI)` = fmt_ci(h3o$rmst_diff, h3o$rmst_low, h3o$rmst_high, 1, " to "),
              check.names = FALSE), "supp_per_pair_metrics.csv")
wr(rbind(data.frame(signature = sig_label(psens$signature), cohorts = psens$cohorts, interval = psens$interval, k = psens$k,
                    tau2 = round(psens$tau2, 4), q = round(psens$q, 3), hr = round(psens$hr, 3),
                    hr_low = round(psens$hr_low, 3), hr_high = round(psens$hr_high, 3))),
   "supp_yao_pooling_intervals.csv")
wr(data.frame(Cohort = pw$cohort, `Evaluable signatures` = pw$evaluable_signatures, Patients = pw$n, Deaths = pw$deaths,
              `True C-index` = pw$target_c, alpha = signif(pw$alpha, 3), Power = pw$power, Replicates = pw$replicates,
              check.names = FALSE), "supp_power.csv")
wr(data.frame(Setting = sh1$setting, Signature = sig_label(sh1$signature), Cohort = sh1$cohort, Patients = sh1$patients,
              Deaths = sh1$deaths, `C-index` = round(sh1$c_index, 3), p = fmt_p(sh1$p_two_sided),
              `p (Holm)` = fmt_p(sh1$p_holm), Classification = sh1$classification, check.names = FALSE),
   "supp_sensitivity_h1_detail.csv")
for (f in c("direction_consistency_observations.csv", "gene_concordance_by_pair.csv", "gene_concordance_by_gene.csv",
            "jaccard_overlap.csv", "signature_clusters.csv", "extraction_agreement.csv")) {
  if (file.exists(file.path("results/tables", f))) wr(rd(f), paste0("supp_", f))
}
wr(ind, "supp_independence.csv")

# ---- Figure: coverage heatmap -------------------------------------------------------------
draw_coverage <- function() {
  ns <- length(sig_ids); nc <- length(cohort_order)
  op <- graphics::par(mar = c(6, 9, 2, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0.5, nc + 0.5), ylim = c(0.5, ns + 0.5), axes = FALSE, xlab = "", ylab = "")
  pal <- grDevices::colorRampPalette(c("#FFFFFF", "#9ECAE1", "#08519C"))(101)
  for (i in seq_len(ns)) for (j in seq_len(nc)) {
    d <- cov[cov$signature == sig_ids[i] & cov$cohort == cohort_order[j], ]
    if (!nrow(d)) next
    y <- ns + 1 - i
    used <- any(ind$used_by_authors[ind$signature == sig_ids[i] & ind$cohort == cohort_order[j]] == "yes")
    graphics::rect(j - 0.5, y - 0.5, j + 0.5, y + 0.5, col = pal[round(100 * d$frac_present) + 1], border = "grey80")
    ev <- d$frac_present >= cfg$coverage$primary_min_fraction - 1e-9
    if (ev) graphics::rect(j - 0.45, y - 0.45, j + 0.45, y + 0.45, border = "black", lwd = 1.6)
    graphics::text(j, y, paste0(d$n_present, "/", d$n_genes, if (used) "*" else ""), cex = 0.62,
                   col = if (d$frac_present > 0.6) "white" else "black")
  }
  graphics::axis(1, at = seq_len(nc), labels = cohort_order, las = 2, tick = FALSE, cex.axis = 0.75)
  graphics::axis(2, at = rev(seq_len(ns)), labels = sig_label(sig_ids), las = 1, tick = FALSE, cex.axis = 0.75)
  graphics::mtext("Signature genes measurable in each cohort", side = 3, line = 0.5, cex = 0.85, font = 2)
  graphics::mtext("Black frame: >= 80% of genes (evaluable). *: cohort used by the authors (not independent).",
                  side = 1, line = 4.8, cex = 0.62)
}
save_figure(file.path(out_f, "fig_coverage_heatmap"), draw_coverage, width = 7, height = 6)

# ---- Figure: C-index in TCGA-LUAD vs independent cohorts ---------------------------------
draw_cpair <- function() {
  ref <- h3[h3$role == "reference", c("signature", "c_index")]
  d <- merge(h3[h3$role == "independent", c("signature", "cohort", "c_index", "c_low", "c_high")], ref,
             by = "signature", suffixes = c("", "_tcga"))
  op <- graphics::par(mar = c(4.5, 4.5, 1.5, 1)); on.exit(graphics::par(op))
  lim <- c(0.35, 0.8)
  graphics::plot(NA, xlim = c(0.58, 0.72), ylim = lim, xlab = "C-index in TCGA-LUAD (apparent performance)",
                 ylab = "C-index in the independent cohort (95% CI)")
  graphics::abline(0, 1, lty = 3, col = "grey50"); graphics::abline(h = 0.5, lty = 2, col = "grey50")
  jit <- (match(d$cohort, cohort_order) - 4) * 0.0012
  graphics::segments(d$c_index_tcga + jit, d$c_low, d$c_index_tcga + jit, d$c_high, col = cols[d$cohort], lwd = 0.8)
  graphics::points(d$c_index_tcga + jit, d$c_index, pch = 19, col = cols[d$cohort], cex = 0.9)
  present <- cohort_order[cohort_order %in% d$cohort]
  graphics::legend("topleft", legend = c(present, "equal to TCGA-LUAD", "no discrimination (0.5)"),
                   col = c(cols[present], "grey50", "grey50"), pch = c(rep(19, length(present)), NA, NA),
                   lty = c(rep(NA, length(present)), 3, 2), bty = "n", cex = 0.72)
}
save_figure(file.path(out_f, "fig_cindex_tcga_vs_independent"), draw_cpair, width = 6, height = 5.5)

# ---- Figure: HR per SD by signature and cohort --------------------------------------------
rows <- do.call(rbind, lapply(sig_ids, function(s) {
  d <- h3[h3$signature == s, ]
  d <- d[order(match(d$cohort, cohort_order)), ]
  rbind(data.frame(label = sig_label(s), est = NA, lo = NA, hi = NA, col = "black", bold = TRUE, right = "", pch = 15),
        data.frame(label = paste0(d$cohort, ifelse(d$role == "reference", " (ref.)", "")), est = d$hr, lo = d$hr_low,
                   hi = d$hr_high, col = unname(cols[d$cohort]), bold = FALSE,
                   right = fmt_ci(d$hr, d$hr_low, d$hr_high), pch = ifelse(d$role == "reference", 1, 15)))
}))
save_figure(file.path(out_f, "fig_hr_by_signature"), function() {
  forest_plot(rows, xlim = c(0.6, 3), right_title = "HR (95% CI)")
}, width = 7.5, height = 12)

# ---- Figure: Yao 2021 forest -------------------------------------------------------------
yd <- h3[h3$signature == yao$signature, ]
yd <- yd[order(yd$role != "reference", match(yd$cohort, cohort_order)), ]
yrows <- rbind(
  data.frame(label = paste0(yd$cohort, ifelse(yd$role == "reference", " (development; not pooled)", "")),
             est = yd$hr, lo = yd$hr_low, hi = yd$hr_high, col = unname(cols[yd$cohort]),
             pch = ifelse(yd$role == "reference", 1, 15), cex = 1, bold = FALSE,
             right = sprintf("%s  %d deaths", fmt_ci(yd$hr, yd$hr_low, yd$hr_high), yd$deaths)),
  data.frame(label = "Pooled, REML + HK (registered)", est = yao$hr, lo = yao$hr_low, hi = yao$hr_high, col = "black",
             pch = 18, cex = 1.8, bold = TRUE, right = sprintf("%s  %d deaths", fmt_ci(yao$hr, yao$hr_low, yao$hr_high), yao$deaths)),
  data.frame(label = "Pooled, modified HK", est = v("reml_hk_modified")$hr, lo = v("reml_hk_modified")$hr_low,
             hi = v("reml_hk_modified")$hr_high, col = "grey30", pch = 18, cex = 1.8, bold = FALSE,
             right = fmt_ci(v("reml_hk_modified")$hr, v("reml_hk_modified")$hr_low, v("reml_hk_modified")$hr_high)))
save_figure(file.path(out_f, "fig_yao2021_forest"), function() {
  forest_plot(yrows, xlim = c(0.8, 2.5), right_title = "HR per SD (95% CI)",
              main = sprintf("%s: published score in independent cohorts", sig_label(yao$signature)))
}, width = 8, height = 4)

# ---- Figure: aggregate transportability across settings ----------------------------------
trows <- data.frame(label = str_$setting, est = exp(str_$pooled_log_hr), lo = exp(str_$ci_low), hi = exp(str_$ci_high),
                    col = ifelse(grepl("^primary", str_$setting), "black", "grey35"), pch = 15, cex = 1,
                    bold = grepl("^primary", str_$setting),
                    right = sprintf("%s  %d pairs", fmt_ci(exp(str_$pooled_log_hr), exp(str_$ci_low), exp(str_$ci_high)), str_$n_pairs))
short <- c("RNA-seq TPM rescaled to protein-coding + lncRNA genes" = "TPM rescaled (protein-coding + lncRNA)",
           "GSE31210 without the 22 flagged cases" = "GSE31210 without 22 flagged cases",
           "deaths within 3 months excluded" = "deaths within 3 months excluded")
trows$label <- ifelse(trows$label %in% names(short), short[trows$label], trows$label)
trows$label[grepl("^primary", trows$label)] <- sprintf("primary (B = %d)", tr$resamples)
trows$est[grepl("^primary", str_$setting)] <- exp(tr$pooled_log_hr)
trows$lo[grepl("^primary", str_$setting)] <- exp(tr$ci_low)
trows$hi[grepl("^primary", str_$setting)] <- exp(tr$ci_high)
trows$right[grepl("^primary", str_$setting)] <- sprintf("%s  %d pairs", fmt_ci(exp(tr$pooled_log_hr), exp(tr$ci_low), exp(tr$ci_high)), tr$n_pairs)
if (!is.null(rc)) trows <- rbind(trows, data.frame(label = "cross-check: multilevel model, CR2", est = exp(rc$log_hr),
                                                   lo = exp(rc$ci_low), hi = exp(rc$ci_high), col = "grey35", pch = 5,
                                                   cex = 1, bold = FALSE, right = fmt_ci(exp(rc$log_hr), exp(rc$ci_low), exp(rc$ci_high))))
save_figure(file.path(out_f, "fig_transportability_settings"), function() {
  forest_plot(trows, xlab = "Mean HR per SD over independent pairs (95% CI)", xlim = c(0.9, 1.7),
              extra_lines = exp(cfg$transportability$delta), extra_labels = "HR 1.25 (delta)",
              right_title = "HR (95% CI)", left = 15)
}, width = 8, height = 4.5)

# ---- Supplementary figure: observed C-index against the random-gene null -----------------
draw_null <- function() {
  d <- h1[order(match(h1$cohort, cohort_order), sig_label(h1$signature)), ]
  n <- nrow(d); y <- rev(seq_len(n))
  op <- graphics::par(mar = c(4.5, 12, 1.5, 1)); on.exit(graphics::par(op))
  graphics::plot(NA, xlim = c(0.3, 0.8), ylim = c(0.5, n + 0.5), yaxt = "n", ylab = "",
                 xlab = "C-index (95% CI); grey band: random-gene null, 2.5-97.5%", bty = "n")
  graphics::abline(v = 0.5, lty = 2, col = "grey40")
  graphics::rect(d$null_q025, y - 0.3, d$null_q975, y + 0.3, col = "grey88", border = NA)
  graphics::segments(d$c_low, y, d$c_high, y, col = cols[d$cohort])
  graphics::points(d$c_index, y, pch = 19, col = cols[d$cohort], cex = 0.8)
  graphics::mtext(paste(sig_label(d$signature), "-", d$cohort), side = 2, at = y, las = 1, line = 0.5, cex = 0.65)
}
save_figure(file.path(out_f, "figS_h1_null"), draw_null, width = 7, height = 8)

# ---- Supplementary figure: Kaplan-Meier curves, independent pairs -------------------------
P <- build_pairs(gm, readRDS("data/processed/measurable_ids.rds"), ind, cfg$cohorts$independent,
                 cfg$coverage$primary_min_fraction)
grDevices::pdf(file.path(out_f, "figS_km_independent_pairs.pdf"), width = 8.5, height = 11)
graphics::par(mfrow = c(4, 3), mar = c(4, 4, 2.5, 1))
for (k in seq_len(nrow(P$pairs))) {
  p <- P$pairs[k, ]; g <- P$genes[[p$pair_id]]; co <- cohorts[[p$cohort]]
  s <- score_signature(co$expr, g$ensembl_id, g$coefficient, cfg$scoring$scale)
  km_median_split(co$time, co$event, s, main = paste(sig_label(p$signature), "-", p$cohort))
}
grDevices::dev.off()

# ---- README with the suggested numbering --------------------------------------------------
writeLines(c(
  "# Manuscript tables and figures (v2)", "",
  "Built by scripts/40_paper_outputs.R from the v2 outputs only. Suggested numbering:", "",
  "Main tables: 1 table_signatures; 2 table_cohorts; 3 table_independent_pairs; 4 table_aggregate;",
  "5 table_sensitivity_transportability + table_sensitivity_h1.",
  "Main figures: 1 PRISMA flow (literature search; not built here); 2 fig_coverage_heatmap;",
  "3 fig_cindex_tcga_vs_independent; 4 fig_hr_by_signature; 5 fig_yao2021_forest; 6 fig_transportability_settings.",
  "Supplementary: figS_h1_null, figS_km_independent_pairs, and every supp_*.csv table.", "",
  paste("Built:", format(Sys.time(), "%Y-%m-%d %H:%M"))), "results/paper/README.md")
message("Paper tables: ", length(list.files(out_t)), " files; figures: ", length(list.files(out_f)), " files (results/paper)")
