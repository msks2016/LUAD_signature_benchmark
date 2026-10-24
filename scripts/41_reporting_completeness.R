# 41_reporting_completeness.R - reporting completeness of the 76 reviewed papers (registered,
# Update 2; correction in Update 3). Recomputes every reported number from the coding files:
#   final coding (both coders' resolved codes) .......... docs/coding/reporting_completeness_76papers_2026-10-12_v2.csv
#   primary coding before discussion (dual-coded 10) .... docs/coding/reporting_completeness_76papers_2026-10-12.csv
#   second coder, before discussion ..................... docs/coding/dual_coding_moinul_original.csv
#   blinded recheck of the text-table item (other 66) ... docs/coding/recheck_form_moinul_filled.csv + recheck_key.csv
#   bibliographic record of each paper (PubMed) ....... docs/coding/records_76_pubmed.csv
# Output: results/tables/reporting_completeness_items.csv, reporting_completeness_agreement.csv,
#         reporting_completeness_recheck.csv, results/paper/tables/supp_reporting_completeness_76.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cd <- "docs/coding"
rd <- function(f) utils::read.csv(file.path(cd, f), stringsAsFactors = FALSE, colClasses = "character")
final <- rd("reporting_completeness_76papers_2026-10-12_v2.csv")
pre <- rd("reporting_completeness_76papers_2026-10-12.csv")
sec <- rd("dual_coding_moinul_original.csv")
key <- rd("recheck_key.csv")
rc <- rd("recheck_form_moinul_filled.csv")
items <- c("coef_table", "cutoff", "scaling", "external_validation", "text_table_conflict")
if (nrow(final) != 76 || anyDuplicated(final$paper_id)) stop("The final coding must hold 76 distinct papers")

# 1. Item counts in the final coding
n <- nrow(final)
cnt <- function(item, value) sum(final[[item]] == value)
it <- data.frame(
  item = c("full coefficient table", "numeric cutoff stated", "scale stated (named unit)", "scale implied (RNA-seq, no unit)",
           "scale not stated", "independent external validation", "text-table conflict"),
  papers = c(cnt("coef_table", "full"), cnt("cutoff", "explicit_value"), cnt("scaling", "stated"), cnt("scaling", "implied"),
             cnt("scaling", "not_stated"), cnt("external_validation", "yes"), cnt("text_table_conflict", "yes")),
  of = n, stringsAsFactors = FALSE)
it$percent <- round(100 * it$papers / it$of, 1)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(it, "results/tables/reporting_completeness_items.csv", row.names = FALSE)
keys <- c("rc_coef_full", "rc_cutoff_value", "rc_scale_stated", "rc_scale_implied", "rc_scale_not_stated",
          "rc_external_validation", "rc_text_table_conflict")
for (i in seq_len(nrow(it))) record_number(keys[i], sprintf("%d/%d (%.1f%%)", it$papers[i], it$of[i], it$percent[i]))

# 2. Dual coding of 10 papers, before discussion (primary pre-discussion codes vs second coder)
dual_ids <- sec$paper_id
ag <- coder_agreement(pre[pre$paper_id %in% dual_ids, c("paper_id", items)], sec[, c("paper_id", items)], items)
ag$stage <- "10 dual-coded papers, before discussion"

# 3. Blinded recheck of the text-table item in the other 66 papers: primary code before resolution
#    (the 2026-10-12 file) vs the second coder; then the resolved (final) code.
rc$paper_id <- key$real_id[match(rc$blind_id, key$blind_id)]
if (anyNA(rc$paper_id)) stop("Every blind_id in the recheck form needs a real_id in recheck_key.csv")
rcx <- merge(pre[, c("paper_id", "text_table_conflict")], rc[, c("paper_id", "text_vs_table_conflict", "specific_disagreement", "notes")],
             by = "paper_id")
names(rcx)[names(rcx) == "text_table_conflict"] <- "primary_before"
names(rcx)[names(rcx) == "text_vs_table_conflict"] <- "recheck"
rcx$final <- final$text_table_conflict[match(rcx$paper_id, final$paper_id)]
rcx$final_note <- final$notes[match(rcx$paper_id, final$paper_id)]
ag2 <- coder_agreement(data.frame(paper_id = rcx$paper_id, text_table_conflict = rcx$primary_before),
                       data.frame(paper_id = rcx$paper_id, text_table_conflict = rcx$recheck), "text_table_conflict")
ag2$stage <- "66 papers, blinded recheck (primary code before resolution vs second coder)"
utils::write.csv(rbind(ag, ag2), "results/tables/reporting_completeness_agreement.csv", row.names = FALSE)
changed <- rcx[rcx$primary_before == "yes" | rcx$recheck == "yes" | rcx$final == "yes", ]
utils::write.csv(changed, "results/tables/reporting_completeness_recheck.csv", row.names = FALSE, na = "")
record_number("rc_text_table_recoded", sum(rcx$primary_before != rcx$final),
              note = "papers whose text-table code changed after the recheck and re-reading")
record_number("rc_recheck_agreement", sprintf("%d/%d", ag2$agree, ag2$papers))
record_number("rc_recheck_kappa", fmt_num(ag2$kappa, 2))
record_number("rc_dual_scaling_kappa", fmt_num(ag$kappa[ag$item == "scaling"], 2))

# 4. Supplementary table S19: one row per paper, with the bibliographic record checked in PubMed
#    (docs/coding/records_76_pubmed.csv; paper_id = screening record number)
bib <- rd("records_76_pubmed.csv")
if (!setequal(bib$paper_id, final$paper_id)) stop("records_76_pubmed.csv must list the same 76 paper_id values as the coding")
sup <- final[, c("paper_id", items, "notes")]
sup$second_coding <- ifelse(sup$paper_id %in% dual_ids, "dual-coded (all items)",
                            ifelse(sup$paper_id %in% rcx$paper_id, paste("recheck of text-table item:", rcx$recheck[match(sup$paper_id, rcx$paper_id)]), ""))
sup <- merge(bib[, c("paper_id", "first_author", "year", "title", "journal", "volume", "pages", "pmid", "doi", "benchmark_signature")],
             sup, by = "paper_id")
sup$benchmark_signature <- ifelse(sup$benchmark_signature == "", "", sig_label(sup$benchmark_signature))
sup <- sup[order(as.integer(sup$paper_id)), ]
record_number("rc_records_pubmed", sum(bib$pmid != ""), note = "of 76 records matched to a PubMed record")
record_number("rc_records_bib_corrected", sum(bib$verification_note != ""),
              note = "extraction records whose title, author, journal or DOI differed from PubMed (coding unaffected)")
dir.create("results/paper/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(sup, "results/paper/tables/supp_reporting_completeness_76.csv", row.names = FALSE, na = "")
print(it, row.names = FALSE)
print(rbind(ag, ag2)[, c("stage", "item", "papers", "agree", "kappa", "disagreeing_papers")], row.names = FALSE)
print(changed[, c("paper_id", "primary_before", "recheck", "final")], row.names = FALSE)
