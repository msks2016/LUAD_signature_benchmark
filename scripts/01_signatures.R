# 01_signatures.R - Step 1 output check.
# Input: data/signatures/signature_master.csv (copy of the TEMPLATE, filled after blind
# double extraction). Stops until both extractors agree on every coefficient.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

sig <- utils::read.csv("data/signatures/signature_master.csv", stringsAsFactors = FALSE)
required <- c("signature", "doi", "published_symbol", "coefficient", "annotation",
              "scale_stated", "source_location", "coef_extractor1", "coef_extractor2")
missing_cols <- setdiff(required, names(sig))
if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

if (anyNA(sig$coef_extractor2)) stop("coef_extractor2 is not filled for every gene")
sig$notes[is.na(sig$notes)] <- ""
n_doi <- tapply(sig$doi, sig$signature, function(x) length(unique(x)))
if (any(n_doi != 1)) stop("Each signature needs one DOI; check: ", paste(names(n_doi)[n_doi != 1], collapse = ", "))

# Agreement of the two blind extractions. "rounding" = one extractor rounded the value printed
# in the paper (same sign, difference below half a unit of its last decimal); any other
# difference stops the script until it is resolved against the PDF.
decimals_of <- function(x) {
  s <- sub("0+$", "", formatC(abs(x), format = "f", digits = 12))
  nchar(sub("^[0-9]*\\.?", "", s))
}
e1 <- sig$coef_extractor1
e2 <- sig$coef_extractor2
same <- abs(e1 - e2) <= 1e-12
rounding <- !same & sign(e1) == sign(e2) & abs(e1 - e2) <= 0.5 * 10^(-decimals_of(e1)) + 1e-12
real <- !same & !rounding
if (any(real)) {
  print(sig[real, c("signature", "published_symbol", "coef_extractor1", "coef_extractor2")])
  stop("Resolve these disagreements against the PDF, then fill `coefficient` with the agreed value.")
}
if (anyNA(sig$coefficient) || any(sig$coefficient == 0)) stop("`coefficient` must hold the agreed, non-zero value")
if (any(!same & !(abs(sig$coefficient - e1) < 1e-12 | abs(sig$coefficient - e2) < 1e-12))) {
  stop("Where the extractors differ, `coefficient` must be one of the two extracted values")
}
if (any(!same & !nzchar(trimws(sig$notes)))) stop("Document every rounding difference in `notes`")
agreement <- data.frame(sig[, c("signature", "published_symbol", "coef_extractor1", "coef_extractor2", "coefficient")],
                        status = ifelse(same, "identical", "rounding only"), stringsAsFactors = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(agreement, "results/tables/extraction_agreement.csv", row.names = FALSE)
record_number("extraction_identical", sum(same))
record_number("extraction_rounding_only", sum(rounding))
record_number("extraction_disagreements", sum(real))
if (anyDuplicated(paste(sig$signature, sig$published_symbol))) stop("A gene is listed twice in one signature")
if (any(grepl("TO_FILL", sig$annotation))) stop("Fill `annotation` (the gene table the paper used, e.g. v22 or v36)")

table1 <- do.call(rbind, lapply(split(sig, sig$signature), function(d) {
  data.frame(signature = d$signature[1], doi = d$doi[1], n_genes = nrow(d),
             n_risk = sum(d$coefficient > 0), n_protective = sum(d$coefficient < 0),
             annotation = d$annotation[1], scale_stated = d$scale_stated[1],
             stringsAsFactors = FALSE)
}))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(table1, "results/tables/table1_signatures.csv", row.names = FALSE)

record_number("n_signatures", nrow(table1))
record_number("n_signature_genes", nrow(sig))
record_number("signature_size_range", paste(range(table1$n_genes), collapse = "-"))
record_number("signature_size_mean", fmt_num(mean(table1$n_genes), 1))
message("Signatures checked: ", nrow(table1), " signatures, ", nrow(sig), " genes")
