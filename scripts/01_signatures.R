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
disagree <- sig[abs(sig$coef_extractor1 - sig$coef_extractor2) > 1e-12, ]
if (nrow(disagree)) {
  print(disagree[, c("signature", "published_symbol", "coef_extractor1", "coef_extractor2")])
  stop("Resolve these disagreements against the PDF, then fill `coefficient` with the agreed value.")
}
if (anyNA(sig$coefficient) || any(sig$coefficient == 0)) stop("`coefficient` must hold the agreed, non-zero value")
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
