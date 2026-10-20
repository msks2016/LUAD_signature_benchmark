# load_data.R - building cohort objects from raw public files (used by scripts 10-13).
# Every cohort ends as list(name, expr = log2 matrix [Ensembl ID x patient], clin, time, event, meta).

# ---- GDC (TCGA-LUAD, CPTAC-3) --------------------------------------------------------

# GDC clinical table -> standard clinical columns.
make_clin_gdc <- function(clin_raw) {
  vs <- clin_raw$vital_status
  dead <- vs == "Dead"
  os_days <- ifelse(dead, clin_raw$days_to_death, clin_raw$days_to_last_follow_up)
  sex <- if (!is.null(clin_raw$sex_at_birth)) clin_raw$sex_at_birth else clin_raw$gender
  data.frame(patient_id = clin_raw$submitter_id,
             os_months = suppressWarnings(as.numeric(os_days)) / 30.4375,
             os_event = ifelse(vs %in% c("Dead", "Alive"), as.integer(dead), NA_integer_),
             stage_raw = clin_raw$ajcc_pathologic_stage,
             stage = collapse_stage(clin_raw$ajcc_pathologic_stage),
             age = suppressWarnings(as.numeric(clin_raw$age_at_index)),
             sex = if (is.null(sex)) NA_character_ else sex,
             stringsAsFactors = FALSE)
}

# SummarizedExperiment from TCGAbiolinks::GDCprepare() -> cohort object.
# Primary tumours only; one sample per patient = first barcode alphabetically (outcome-blind).
build_gdc_cohort <- function(se, clin_raw, name) {
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  n_files <- nrow(cd)
  cd <- cd[cd$sample_type == "Primary Tumor", , drop = FALSE]
  n_tumour <- nrow(cd)
  cd <- cd[order(cd$patient, cd$barcode), , drop = FALSE]
  cd <- cd[!duplicated(cd$patient), , drop = FALSE]
  tpm <- SummarizedExperiment::assay(se, "tpm_unstrand")[, cd$barcode, drop = FALSE]
  ids <- rownames(tpm)
  keep <- grepl("^ENSG", ids) & !grepl("_PAR_Y$", ids)
  expr <- log2(as.matrix(tpm[keep, , drop = FALSE]) + 1)
  rownames(expr) <- strip_version(ids[keep])
  colnames(expr) <- cd$patient
  meta <- list(source = "GDC STAR-Counts, tpm_unstranded, log2(TPM + 1)",
               gene_model = "GENCODE v36",
               files = n_files, primary_tumour_samples = n_tumour, patients_with_rna = nrow(cd))
  finalise_cohort(expr, make_clin_gdc(clin_raw), name, platform = "RNA-seq", meta = meta)
}

# Share of TPM outside the genes in keep_ids, per sample (RNA-seq; expr = log2(TPM + 1)).
# Abundant short or mitochondrial RNAs (e.g. RN7SL1 in total-RNA libraries) take TPM away from
# every other gene by a sample-specific amount.
tpm_share_outside <- function(expr, keep_ids) {
  tpm <- 2^expr - 1
  1 - colSums(tpm[rownames(tpm) %in% keep_ids, , drop = FALSE]) / colSums(tpm)
}

# Sensitivity analysis: rescale TPM so that the genes in keep_ids (protein-coding and lncRNA)
# sum to 10^6 in every sample, then back to log2(TPM + 1). Every gene is kept; only the
# per-sample scaling changes. Microarray cohorts are returned unchanged.
rescale_tpm_cohort <- function(co, keep_ids) {
  if (!identical(co$meta$platform, "RNA-seq")) return(co)
  tpm <- 2^co$expr - 1
  in_keep <- rownames(tpm) %in% keep_ids
  if (!any(in_keep)) stop(co$name, ": none of keep_ids is in the expression matrix")
  factor <- 1e6 / colSums(tpm[in_keep, , drop = FALSE])
  co$expr <- log2(sweep(tpm, 2, factor, "*") + 1)
  co$meta$tpm_rescaled <- "TPM rescaled so that protein-coding and lncRNA genes sum to 1e6"
  co
}

# ---- GEO (GSE31210, GSE50081, GSE30219, GSE81089) ------------------------------------

geo_pheno <- function(gse, destdir = "data/raw/GEO") {
  dir.create(destdir, recursive = TRUE, showWarnings = FALSE)
  g <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = destdir)
  Biobase::pData(g[[1]])
}

# Print every characteristic field with its most common values (use it to fill config/cohorts.yaml).
show_geo_fields <- function(gse) {
  pd <- geo_pheno(gse)
  for (f in grep(":ch1$", names(pd), value = TRUE)) {
    v <- sort(table(pd[[f]], useNA = "ifany"), decreasing = TRUE)
    cat("\n", f, "  (", length(v), " distinct values)\n", sep = "")
    print(utils::head(v, 8))
  }
  invisible(pd)
}

# Stage group (I / II / III-IV) from pathological T, N and M categories, for series that give
# TNM but no stage. AJCC 6th-edition grouping, because such series do not split T into a/b:
#   I = T1-2 N0 M0;  II = T1-2 N1 M0 or T3 N0 M0;  III-IV = T3 N1-2, T4, N2-3 or M1.
# MX or a missing M counts as M0 (resected series); TX or NX gives NA unless the stage is
# III-IV anyway (for example N2 with TX).
stage_from_tnm <- function(t, n, m) {
  cat_num <- function(x, letter) {
    x <- toupper(trimws(as.character(x)))
    suppressWarnings(as.integer(ifelse(grepl(paste0("^", letter, "[0-4]"), x),
                                       substr(x, 2, 2), NA_character_)))
  }
  tt <- cat_num(t, "T")
  nn <- cat_num(n, "N")
  mm <- cat_num(m, "M")
  m0 <- is.na(mm) | mm == 0
  out <- rep(NA_character_, length(tt))
  out[m0 & tt %in% c(1, 2) & nn %in% 0] <- "I"
  out[m0 & ((tt %in% c(1, 2) & nn %in% 1) | (tt %in% 3 & nn %in% 0))] <- "II"
  out[mm %in% 1 | nn %in% c(2, 3) | tt %in% 4 | (tt %in% 3 & nn %in% 1)] <- "III-IV"
  factor(out, levels = c("I", "II", "III-IV"))
}

# GEO phenotype table + one entry of config/cohorts.yaml -> standard clinical columns.
# Entry fields: time {field, unit} or {start_date, end_date}, event {field, dead_values,
# alive_values}, stage (field),
# stage_recode (raw value -> stage, e.g. "1A" -> "IA"), stage_tnm {tumour, node, metastasis}
# (field names, used instead of stage),
# age (one field, or several combined: first non-empty value), sex, sensitivity_exclude.
# An event value that is neither dead nor alive becomes NA (not "alive").
make_clin_geo <- function(pd, f) {
  field <- function(cols) {
    if (is.null(cols)) return(rep(NA_character_, nrow(pd)))
    absent <- setdiff(cols, names(pd))
    if (length(absent)) stop("Field not found in GEO table: ", paste(absent, collapse = ", "))
    out <- rep(NA_character_, nrow(pd))
    for (col in cols) {
      v <- trimws(as.character(pd[[col]]))
      fill <- is.na(out) & !is.na(v) & v != ""
      out[fill] <- v[fill]
    }
    out
  }
  if (!is.null(f$time$start_date)) {          # follow-up = end date - start date (YYYY-MM-DD)
    d0 <- as.Date(field(f$time$start_date), format = "%Y-%m-%d")
    d1 <- as.Date(field(f$time$end_date), format = "%Y-%m-%d")
    os_months <- as.numeric(d1 - d0) / 30.4375
  } else {
    per_month <- c(days = 1 / 30.4375, months = 1, years = 12)[[f$time$unit]]
    os_months <- suppressWarnings(as.numeric(field(f$time$field))) * per_month
  }
  ev <- tolower(field(f$event$field))
  dead <- ev %in% tolower(f$event$dead_values)
  known <- if (is.null(f$event$alive_values)) !is.na(ev) & ev != "" else
    dead | ev %in% tolower(f$event$alive_values)
  if (!is.null(f$stage_tnm)) {
    # keys are tumour / node / metastasis: YAML reads an unquoted key "n" as FALSE
    parts <- c("tumour", "node", "metastasis")
    if (!all(parts %in% names(f$stage_tnm))) {
      stop("stage_tnm needs the keys tumour, node and metastasis (found: ",
           paste(names(f$stage_tnm), collapse = ", "), ")")
    }
    tnm <- lapply(f$stage_tnm[parts], field)
    stage_raw <- paste(tnm$tumour, tnm$node, tnm$metastasis)
    stage <- stage_from_tnm(tnm$tumour, tnm$node, tnm$metastasis)
  } else {
    stage_raw <- field(f$stage)
    st <- stage_raw
    if (!is.null(f$stage_recode)) st <- unname(unlist(f$stage_recode)[st])
    stage <- collapse_stage(st)
  }
  out <- data.frame(patient_id = pd$geo_accession,
                    os_months = os_months,
                    os_event = ifelse(known, as.integer(dead), NA_integer_),
                    stage_raw = stage_raw,
                    stage = stage,
                    age = suppressWarnings(as.numeric(field(f$age))),
                    sex = field(f$sex),
                    stringsAsFactors = FALSE)
  if (!is.null(f$sensitivity_exclude)) {
    flag <- tolower(field(f$sensitivity_exclude$field))
    out$sensitivity_exclude <- flag %in% tolower(f$sensitivity_exclude$values)
  }
  out
}

# Size of a remote file from its HTTP headers (NA when the server does not say).
remote_size <- function(url) {
  h <- tryCatch(curlGetHeaders(url), error = function(e) character(0))
  cl <- grep("^content-length:", h, ignore.case = TRUE, value = TRUE)
  if (!length(cl)) return(NA_real_)
  as.numeric(trimws(sub("^[^:]*:", "", cl[length(cl)])))
}

# Download one GEO supplementary file once. It is written under a temporary name and kept only
# when its size matches the server's, so an interrupted download is never mistaken for a
# finished one. Run again after a failure: the download simply starts over.
geo_download_supp <- function(gse, fname, destdir = "data/raw/GEO", tries = 3,
                              base = GEOquery::getGEOSuppFileURL(gse)) {
  dir <- file.path(destdir, gse)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, fname)
  if (file.exists(path)) return(path)
  url <- paste0(base, fname)
  expected <- remote_size(url)
  part <- paste0(path, ".part")
  old <- options(timeout = max(10800, getOption("timeout")))
  on.exit(options(old))
  for (t in seq_len(tries)) {
    message(gse, ": downloading ", fname,
            if (!is.na(expected)) sprintf(" (%.0f MB)", expected / 1e6) else "", " ...")
    status <- tryCatch(utils::download.file(url, part, mode = "wb", quiet = FALSE),
                       error = function(e) 1L, warning = function(w) 1L)
    if (isTRUE(status == 0) && file.exists(part) &&
        (is.na(expected) || file.size(part) == expected)) {
      file.rename(part, path)
      return(path)
    }
    unlink(part)
    if (t < tries) Sys.sleep(30 * t)
  }
  stop(gse, ": download of ", fname, " failed ", tries, " times. Run the script again.")
}

# The CEL archive of a series (GSExxxx_RAW.tar).
geo_download_raw <- function(gse, destdir = "data/raw/GEO", tries = 3,
                             base = GEOquery::getGEOSuppFileURL(gse)) {
  geo_download_supp(gse, paste0(gse, "_RAW.tar"), destdir, tries, base)
}

# FPKM matrix (genes x samples) -> log2(TPM + 1). TPM = FPKM / sum of FPKM in the sample x 1e6.
fpkm_to_log2tpm <- function(fpkm) {
  stopifnot(is.matrix(fpkm), is.numeric(fpkm), !anyNA(fpkm), all(fpkm >= 0))
  log2(sweep(fpkm, 2, colSums(fpkm), "/") * 1e6 + 1)
}

# RMA (background correction, quantile normalisation, log2) on the CEL files of the chosen
# samples only, so each cohort is normalised on its own analysis set. Stops if any selected
# sample has no CEL file (for example after an interrupted extraction).
rma_from_cel <- function(gse, gsm_ids, destdir = "data/raw/GEO") {
  tar <- geo_download_raw(gse, destdir)
  cel_dir <- file.path(destdir, gse, "CEL")
  find_cels <- function() {
    cels <- list.files(cel_dir, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE, full.names = TRUE)
    cels[sub("^(GSM[0-9]+).*$", "\\1", basename(cels)) %in% gsm_ids]
  }
  cels <- find_cels()
  if (length(cels) < length(gsm_ids)) {
    utils::untar(tar, exdir = cel_dir)
    cels <- find_cels()
  }
  gsm <- sub("^(GSM[0-9]+).*$", "\\1", basename(cels))
  absent <- setdiff(gsm_ids, gsm)
  if (length(absent)) stop(gse, ": no CEL file for ", length(absent), " selected samples, e.g. ", absent[1])
  if (anyDuplicated(gsm)) stop(gse, ": several CEL files for one sample")
  eset <- affy::justRMA(filenames = cels)
  expr <- Biobase::exprs(eset)
  colnames(expr) <- gsm
  expr
}

# GPL570 probe -> Ensembl gene map from a pinned Ensembl release, saved once and reused.
# BioMart can be slow or drop connections, so every request is tried up to three times.
gpl570_probe_map <- function(probes, ensembl_version = 114,
                             file = "data/reference/gpl570_probe_map_ensembl.csv") {
  if (file.exists(file)) return(utils::read.csv(file, stringsAsFactors = FALSE))
  retry <- function(what, expr_fun) {
    for (t in 1:3) {
      r <- tryCatch(expr_fun(), error = function(e) e)
      if (!inherits(r, "error")) return(r)
      message("  Ensembl BioMart (", what, "): ", conditionMessage(r), " - trying again")
      Sys.sleep(20 * t)
    }
    stop("Ensembl BioMart did not answer (", what, "). Run the script again later.")
  }
  mart <- retry("connect", function() {
    biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = ensembl_version)
  })
  chunks <- split(probes, ceiling(seq_along(probes) / 5000))
  pm <- do.call(rbind, lapply(seq_along(chunks), function(k) {
    retry(paste("chunk", k, "of", length(chunks)), function() {
      biomaRt::getBM(attributes = c("affy_hg_u133_plus_2", "ensembl_gene_id"),
                     filters = "affy_hg_u133_plus_2", values = chunks[[k]], mart = mart)
    })
  }))
  names(pm) <- c("probe_id", "ensembl_id")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(pm, file, row.names = FALSE)
  pm
}

# Keep a subset of patients in a cohort object (sensitivity analyses).
subset_cohort <- function(co, keep) {
  co$expr <- co$expr[, keep, drop = FALSE]
  co$clin <- co$clin[keep, , drop = FALSE]
  co$time <- co$time[keep]
  co$event <- co$event[keep]
  co
}
