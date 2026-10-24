# supplement.R - Supplementary Tables S1-S19 as one Excel workbook (script 42).
# Only labels change (column names, signature names, codes); values are copied from the tables
# that scripts 40 and 41 write to results/paper/tables.

# Signature IDs to labels ("Paper62_Yao2021" -> "Yao 2021"); other values unchanged.
sig_or_value <- function(x) ifelse(grepl("^Paper[0-9]+_", x), sig_label(x), x)

# Rename columns by name; stops if a column is missing (the table layout has changed).
rename_cols <- function(d, map) {
  miss <- setdiff(names(map), names(d))
  if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "))
  names(d)[match(names(map), names(d))] <- unname(map)
  d
}

# Replace codes by labels; values without a label are kept.
recode <- function(x, map) {
  i <- match(x, names(map))
  ifelse(is.na(i), x, unname(map[i]))
}

# Read a table as text, then turn columns into numbers only where no formatting is lost
# ("0.56", "33850", "-0.0328" become numbers; "1.00", "<0.001", "0.47-0.64" stay text) and
# TRUE/FALSE columns into yes/no.
read_supp_csv <- function(path) {
  if (!file.exists(path)) stop("Missing ", path, ": run scripts 40 and 41 first")
  d <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
  for (j in seq_along(d)) {
    v <- d[[j]]
    v[is.na(v)] <- ""
    full <- v[v != ""]
    if (!length(full)) { d[[j]] <- v; next }
    if (all(full %in% c("TRUE", "FALSE"))) {
      d[[j]] <- ifelse(v == "TRUE", "yes", ifelse(v == "FALSE", "no", ""))
    } else if (all(grepl("^-?[0-9]+$", full) | grepl("^-?[0-9]*\\.[0-9]*[1-9]$", full))) {
      d[[j]] <- ifelse(v == "", NA_real_, suppressWarnings(as.numeric(v)))
    } else {
      d[[j]] <- v
    }
  }
  d
}

# Titles, notes and the tidied data of Supplementary Tables S1-S19, in manuscript order.
supp_workbook_tables <- function(td = "results/paper/tables") {
  rd <- function(f) read_supp_csv(file.path(td, f))
  sig_cols <- function(d, cols) { for (cc in cols) d[[cc]] <- sig_or_value(d[[cc]]); d }
  dash <- "–"
  tab <- function(id, file, title, note, fix = identity) {
    list(id = id, file = file, title = title, note = note, data = fix(rd(file)))
  }
  list(
    tab("S1", "supp_positive_control.csv", "Stage as a positive control in each cohort",
        paste0("Cox model of overall survival on stage (I", dash, "IV); patients with known stage. HR, hazard ratio against stage I. ",
               "Rule: the registered criterion, HR(III vs I) > 1 with p < 0.05 and hazard ratios increasing with stage; in cohorts with fewer than ",
               "10 stage III", dash, "IV patients, HR(II vs I) > 1 with p < 0.05 (Update 3). Stage C-index: C-index of stage alone. ",
               "A cohort that fails is analysed as exploratory.")),
    tab("S2", "supp_coverage_matrix.csv", "Signature genes measurable in each cohort",
        paste0("Measurable genes / signature genes. A gene is measurable if it has no missing values, is not constant and, for RNA-seq, ",
               "has TPM > 0 in at least 10% of patients. A signature is evaluable in a cohort if at least 80% of its genes are measurable.")),
    tab("S3", "supp_coverage_long.csv", "Coverage by signature and cohort, with the missing genes",
        paste0("One row per signature and cohort. Present through a successor ID: genes measured under the official successor of a retired ",
               "or merged Ensembl identifier, used only where the original identifier is not measurable in that cohort.")),
    tab("S4", "supp_gene_identity.csv", "Ensembl identity of every published gene name",
        paste0("One row per signature gene. Annotation release of the paper: the GENCODE release whose gene names the paper followed. ",
               "Source of the ID: the release used for the mapping, or a documented manual mapping (config/gene_overrides.csv) with its ",
               "evidence. Successor Ensembl ID: official successor of a retired or merged identifier (Ensembl ID history), used in a cohort ",
               "only where the original identifier is not measurable."),
        function(d) {
          d$id_source <- recode(d$id_source, c(v36 = "GENCODE v36", v22 = "GENCODE v22", manual_override = "manual mapping"))
          d$annotation <- recode(d$annotation, c(v36 = "GENCODE v36", v22 = "GENCODE v22"))
          rename_cols(d, c(signature = "Signature", published_symbol = "Published gene name", coefficient = "Published coefficient",
                           annotation = "Annotation release of the paper", ensembl_id = "Ensembl gene ID", id_source = "Source of the ID",
                           status = "Status", alternative_id = "Successor Ensembl ID", override_evidence = "Evidence for the manual mapping",
                           successor_evidence = "Evidence for the successor ID"))
        }),
    tab("S5", "supp_missing_genes_reasons.csv", "Signature genes not measurable in each cohort, with the reason",
        paste0("One row per missing gene and cohort. In Ensembl 114: whether the identifier exists in Ensembl release 114. ",
               "Replaced by: the current identifier according to the Ensembl ID history.")),
    tab("S6", "supp_h1_null_test.csv", "H1: published scores compared with expression-matched random gene sets (independent pairs)",
        paste0("C-index: Harrell's C of the published score (higher score = higher risk), 95% CI from 1,000 bootstrap resamples; ",
               "Cox-model C-index for comparison. Null pool genes: measurable protein-coding and lncRNA genes from which random genes were drawn. ",
               "Null 2.5% and 97.5%: percentiles of the C-index of 2,000 random gene sets, each gene drawn from the same decile of mean ",
               "expression and scored with the published coefficients. p: two-sided; Holm correction across the signatures evaluable in the ",
               "cohort (Holm family = their number). Replicated: adjusted p < 0.05 and C > 0.5; Contradicted: adjusted p < 0.05 and C < 0.5; ",
               "otherwise Inconclusive. Cohort status: valid if the cohort passed the stage positive control (Supplementary Table S1).")),
    tab("S7", "supp_h2_added_value.csv", "H2: added value of the published score beyond stage, sex and age (independent pairs)",
        paste0("Cox models with stage (four levels), sex and age, without and with the standardised score; patients with complete covariates. ",
               "C base and C with score: C-indexes of the two models; Delta C with a 1,000-resample bootstrap 95% CI. LRT: likelihood-ratio ",
               "test (1 df), Holm-corrected within cohort. Adjusted HR per SD: HR of the score in the full model. Adds value: Holm-adjusted p < 0.05.")),
    tab("S8", "supp_per_pair_metrics.csv", "Performance in TCGA-LUAD and in every independent evaluable pair",
        paste0("Reference: TCGA-LUAD, apparent performance (every signature was developed in it). HR per SD of the standardised published ",
               "score. C-index: Harrell's C of the published score with a bootstrap 95% CI. Uno C: Uno's C-index truncated at 60 months. ",
               "td-AUC: time-dependent AUC (inverse-probability-of-censoring weighting). Log-rank test at the median score. PH test: ",
               "Schoenfeld-residual test of proportional hazards. RMST: restricted mean survival time up to tau, patients above minus below ",
               "the median score.")),
    tab("S9", "supp_power.csv", "Power simulation at each independent cohort's numbers of patients and deaths",
        paste0("Proportion of simulated cohorts (Replicates) in which C = 0.5 was rejected, for a score with the true C-index shown; ",
               "exponential censoring matched to the observed proportion of deaths. alpha 0.05: nominal level; smaller alpha: the strictest ",
               "Holm step, 0.05 divided by the number of signatures evaluable in the cohort.")),
    tab("S10", "supp_yao_pooling_intervals.csv", "Yao 2021: pooled HR per SD by interval method, with leave-one-cohort-out estimates",
        paste0("Random-effects meta-analysis (REML) of the HR per SD over the independent cohorts in which Yao 2021 is evaluable ",
               "(Cohorts: all, or all but the one named). Interval methods: REML with the Hartung", dash, "Knapp adjustment (registered, primary); ",
               "REML with a Wald interval; modified Hartung", dash, "Knapp, q* = max(1, q) (post hoc; Update 4, B3); common effect. ",
               "q = Q/(k ", "−", " 1). Pooled estimates are exploratory because every signature was developed in TCGA-LUAD."),
        function(d) {
          d$signature <- sig_or_value(d$signature)
          d$interval <- recode(d$interval, c(reml_hk = paste0("REML, Hartung", dash, "Knapp (registered)"),
                                             reml_hk_modified = paste0("REML, modified Hartung", dash, "Knapp (post hoc)"),
                                             reml_wald = "REML, Wald", common_effect = "Common effect"))
          rename_cols(d, c(signature = "Signature", cohorts = "Cohorts", interval = "Interval method", tau2 = "tau²",
                           q = "q", hr = "HR per SD", hr_low = "95% CI, lower", hr_high = "95% CI, upper"))
        }),
    tab("S11", "supp_sensitivity_h1_detail.csv", "H1 sensitivity analyses, by pair",
        paste0("Settings: random genes drawn from lncRNA genes only or from all measurable genes; stage I", dash, "II only; deaths within ",
               "3 months excluded; GSE31210 without the 22 cases flagged by the submitters. Measurable genes were recomputed in every subset. ",
               "Columns as in Supplementary Table S6.")),
    tab("S12", "supp_direction_consistency_observations.csv", "Direction consistency: one observation per signature cluster and independent cohort",
        paste0("Mean log-HR per SD of the cluster's signatures in the cohort (clusters: Supplementary Table S16). An observation is ",
               "consistent with the published direction if it is above 0; two-sided binomial test against 50%."),
        function(d) {
          d$log_hr <- round(d$log_hr, 3)
          rename_cols(d, c(cluster = "Signature cluster", cohort = "Cohort", log_hr = "Mean log-HR per SD",
                           n_signatures = "Signatures of the cluster evaluable in the cohort"))
        }),
    tab("S13", "supp_gene_concordance_by_pair.csv", "Gene-level concordance by signature and cohort",
        paste0("For every measurable signature gene, the sign of its univariable Cox coefficient in the cohort compared with the sign of its ",
               "published coefficient. Cohort used by the authors: Supplementary Table S18. Independent evaluable pair: the authors did not use ",
               "the cohort and the signature is evaluable there."),
        function(d) {
          d$signature <- sig_or_value(d$signature)
          d$agreement_rate <- round(d$agreement_rate, 3)
          rename_cols(d, c(signature = "Signature", cohort = "Cohort", used_by_authors = "Cohort used by the authors",
                           signature_evaluable_independent = "Independent evaluable pair", genes = "Measurable signature genes",
                           agree = "Genes with the published sign", agreement_rate = "Proportion with the published sign"))
        }),
    tab("S14", "supp_gene_concordance_by_gene.csv", "Gene-level concordance by gene, over the independent cohorts",
        "Counts over the independent cohorts, not used by the signature's authors, in which the gene is measurable.",
        function(d) {
          d$signature <- sig_or_value(d$signature)
          d$agreement_rate <- round(d$agreement_rate, 3)
          rename_cols(d, c(signature = "Signature", published_symbol = "Published gene name", ensembl_id = "Ensembl gene ID",
                           cohorts = "Independent cohorts with the gene measurable", agree = "Cohorts with the published sign",
                           agreement_rate = "Proportion with the published sign"))
        }),
    tab("S15", "supp_jaccard_overlap.csv", "Overlap of the signatures' gene sets (Jaccard index)",
        "Jaccard index of the mapped Ensembl gene sets of each pair of signatures.",
        function(d) {
          d$signature <- sig_or_value(d$signature)
          names(d) <- sig_or_value(names(d))
          rename_cols(d, c(signature = "Signature"))
        }),
    tab("S16", "supp_signature_clusters.csv", "Signature clusters used in the direction-consistency test",
        "Signatures linked by a Jaccard index above 0.3, directly or through other signatures, form one cluster.",
        function(d) rename_cols(sig_cols(d, "signature"), c(signature = "Signature", cluster = "Cluster"))),
    tab("S17", "supp_extraction_agreement.csv", "Agreement of the two independent coefficient extractions",
        "Coefficient used: the full-precision printed value. Rounding only: one extractor had rounded the printed value.",
        function(d) rename_cols(sig_cols(d, "signature"),
                                c(signature = "Signature", published_symbol = "Published gene name", coef_extractor1 = "Coefficient, extractor 1",
                                  coef_extractor2 = "Coefficient, extractor 2", coefficient = "Coefficient used", status = "Agreement"))),
    tab("S18", "supp_independence.csv", "Independence of every signature-cohort pair",
        paste0("Recorded from the full texts by two investigators. A pair is independent only if the authors did not use the cohort for ",
               "development, feature selection or validation."),
        function(d) rename_cols(sig_cols(d, "signature"),
                                c(signature = "Signature", cohort = "Cohort", used_by_authors = "Used by the authors", role = "Use",
                                  evidence = "Evidence"))),
    tab("S19", "supp_reporting_completeness_76.csv", "Reporting completeness of the 76 reviewed papers",
        paste0("Record: screening record number; bibliographic details from PubMed. Full coefficient table: numeric coefficient of every gene ",
               "(full) or not (none). Risk cutoff: method only, or value stated. Expression scale: stated (named unit), implied (RNA-seq ",
               "mentioned without a unit) or not stated. Text-table conflict: the text (including the risk-score formula) and the tables ",
               "disagree on a coefficient's sign, a coefficient's value beyond rounding of the same value, or the number of genes or samples. ",
               "Second coding: dual-coded papers (all items) or the blinded recheck of the text-table item. Notes: reason for every change ",
               "after re-reading."),
        function(d) {
          d$cutoff <- recode(d$cutoff, c(method_only = "method only", explicit_value = "value stated"))
          d$scaling <- recode(d$scaling, c(not_stated = "not stated"))
          rename_cols(d, c(paper_id = "Record", first_author = "First author", year = "Year", title = "Title", journal = "Journal",
                           volume = "Volume (issue)", pages = "Pages", pmid = "PMID", doi = "DOI",
                           benchmark_signature = "Signature in the benchmark", coef_table = "Full coefficient table",
                           cutoff = "Risk cutoff", scaling = "Expression scale", external_validation = "Independent external validation",
                           text_table_conflict = "Text-table conflict", notes = "Notes", second_coding = "Second coding"))
        })
  )
}

# Column widths (characters) from the longest header word and value, between 8 and 60.
supp_col_widths <- function(d) {
  head <- vapply(strsplit(names(d), " "), function(w) max(nchar(w)), 0)
  vals <- vapply(d, function(x) {
    x <- as.character(x); x <- x[!is.na(x)]
    if (length(x)) max(nchar(x)) else 0
  }, 0)
  pmin(60, pmax(8, pmax(head, pmin(vals, 60)) + 2))
}

# Write the workbook: a contents sheet, then one sheet per table with its title (row 1), note
# (row 2) and the table (from row 4). `ox` returns openxlsx functions; tests pass a stub.
write_supp_workbook <- function(tabs, file, heading, ox = function(f) getExportedValue("openxlsx", f)) {
  wb <- ox("createWorkbook")(creator = "")
  bold <- ox("createStyle")(fontSize = 12, textDecoration = "bold")
  wrap <- ox("createStyle")(wrapText = TRUE, valign = "top")
  head <- ox("createStyle")(textDecoration = "bold", fgFill = "#E7E6E6", border = "Bottom", wrapText = TRUE, valign = "top")
  toc <- data.frame(Table = vapply(tabs, function(t) t$id, ""),
                    Title = vapply(tabs, function(t) t$title, ""),
                    Rows = vapply(tabs, function(t) nrow(t$data), 0L), check.names = FALSE, stringsAsFactors = FALSE)
  ox("addWorksheet")(wb, sheetName = "Contents")
  ox("writeData")(wb, sheet = "Contents", x = heading, startCol = 1, startRow = 1)
  ox("addStyle")(wb, sheet = "Contents", style = bold, rows = 1, cols = 1)
  ox("writeData")(wb, sheet = "Contents", x = toc, startCol = 1, startRow = 3, headerStyle = head)
  ox("setColWidths")(wb, sheet = "Contents", cols = 1:3, widths = c(10, 100, 8))
  for (t in tabs) {
    s <- t$id
    d <- t$data
    w <- supp_col_widths(d)
    span <- seq_len(max(length(w), 6))
    note_width <- sum(c(w, rep(8.43, 6))[span])
    ox("addWorksheet")(wb, sheetName = s)
    ox("writeData")(wb, sheet = s, x = paste0("Supplementary Table ", s, ". ", t$title), startCol = 1, startRow = 1)
    ox("addStyle")(wb, sheet = s, style = bold, rows = 1, cols = 1)
    ox("writeData")(wb, sheet = s, x = t$note, startCol = 1, startRow = 2)
    ox("mergeCells")(wb, sheet = s, cols = span, rows = 2)
    ox("addStyle")(wb, sheet = s, style = wrap, rows = 2, cols = 1)
    ox("setRowHeights")(wb, sheet = s, rows = 2, heights = 15 * ceiling(1.1 * nchar(t$note) / note_width) + 4)
    ox("writeData")(wb, sheet = s, x = d, startCol = 1, startRow = 4, headerStyle = head)
    ox("setColWidths")(wb, sheet = s, cols = seq_along(w), widths = w)
    long <- which(w >= 40)
    if (length(long) && nrow(d))
      ox("addStyle")(wb, sheet = s, style = wrap, rows = 4 + seq_len(nrow(d)), cols = long, gridExpand = TRUE)
    ox("freezePane")(wb, sheet = s, firstActiveRow = 5)
  }
  ox("saveWorkbook")(wb, file = file, overwrite = TRUE)
  invisible(file)
}
