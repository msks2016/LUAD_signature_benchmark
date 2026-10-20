# gdc_api.R - reading GDC (Genomic Data Commons) records and files straight from the GDC API.
# Used for CPTAC-3 (scripts/11_cohort_cptac3.R). TCGAbiolinks::GDCquery() stops on CPTAC-3 with
# "object 'sample_type' not found": the GDC now describes some samples only with tissue_type,
# tumor_descriptor, specimen_type and preservation_method, and GDCquery() expects sample_type
# on every sample. Here every field the GDC leaves out simply becomes NA.

# Filter in the GDC API syntax: every condition must hold; each element is field = value(s).
gdc_filter <- function(conditions) {
  stopifnot(is.list(conditions), length(conditions) > 0, !is.null(names(conditions)))
  list(op = "and",
       content = unname(lapply(names(conditions), function(f) {
         # I() keeps a single value as a JSON array, as the GDC API expects
         list(op = "in", content = list(field = f, value = I(as.character(conditions[[f]]))))
       })))
}

# One POST request to the GDC API, tried again after network failures.
gdc_request <- function(endpoint, body, tries = 4, base = "https://api.gdc.cancer.gov/") {
  json <- as.character(jsonlite::toJSON(body, auto_unbox = TRUE))
  url <- paste0(base, endpoint)
  for (i in seq_len(tries)) {
    r <- tryCatch(httr::POST(url, body = json, httr::content_type_json(), httr::timeout(600)),
                  error = function(e) e)
    if (inherits(r, "error")) {
      why <- conditionMessage(r)
    } else {
      code <- httr::status_code(r)
      txt <- httr::content(r, as = "text", encoding = "UTF-8")
      if (code == 200) return(jsonlite::fromJSON(txt, simplifyVector = FALSE))
      if (code >= 400 && code < 500) {
        stop("GDC API refused the request to /", endpoint, " (HTTP ", code, "): ", substr(txt, 1, 300))
      }
      why <- paste("HTTP status", code)
    }
    if (i == tries) {
      stop("GDC API request to /", endpoint, " failed ", tries, " times (last: ", why, "). ",
           "Check the internet connection and run the script again.")
    }
    message("  GDC API: ", why, " - trying again in ", 15 * i, " s")
    Sys.sleep(15 * i)
  }
}

# All records of `endpoint` ("cases" or "files") that match `filters`, as parsed JSON lists.
gdc_hits <- function(endpoint, filters, expand = NULL, size = 10000) {
  body <- list(filters = filters, format = "JSON", size = as.character(size))
  if (length(expand)) body$expand <- paste(expand, collapse = ",")
  res <- gdc_request(endpoint, body)
  hits <- res$data$hits
  total <- res$data$pagination$total
  if (is.null(hits)) stop("Unexpected answer from the GDC API (no data$hits)")
  if (!is.null(total) && total > length(hits)) {
    stop("GDC returned ", length(hits), " of ", total, " records: raise `size`")
  }
  hits
}

# GDC data release on the day of the download (for the methods section).
gdc_release <- function() {
  r <- tryCatch(httr::GET("https://api.gdc.cancer.gov/status", httr::timeout(120)),
                error = function(e) NULL)
  if (is.null(r) || httr::status_code(r) != 200) return(NA_character_)
  s <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))
  if (is.null(s$data_release)) NA_character_ else as.character(s$data_release)
}

# First value of field `name` in one parsed JSON record; NA when the GDC left it out.
gdc_value <- function(x, name) {
  v <- if (is.list(x)) x[[name]] else NULL
  if (is.null(v) || !length(v)) return(NA_character_)
  v <- v[[1]]
  if (is.null(v) || is.list(v)) NA_character_ else as.character(v)
}

# /files records (expanded with cases and cases.samples) -> one row per file x case x sample.
flatten_gdc_files <- function(hits) {
  sample_fields <- c("sample_type", "tissue_type", "tumor_descriptor", "specimen_type",
                     "preservation_method")
  rows <- list()
  for (h in hits) {
    file_id <- gdc_value(h, "file_id")
    if (is.na(file_id)) file_id <- gdc_value(h, "id")
    file_part <- c(file_id = file_id, file_name = gdc_value(h, "file_name"),
                   md5sum = gdc_value(h, "md5sum"), file_size = gdc_value(h, "file_size"),
                   n_cases = as.character(length(h$cases)))
    cases <- if (length(h$cases)) h$cases else list(list())
    for (cs in cases) {
      case_part <- c(patient = gdc_value(cs, "submitter_id"), case_id = gdc_value(cs, "case_id"),
                     primary_site = gdc_value(cs, "primary_site"),
                     disease_type = gdc_value(cs, "disease_type"),
                     n_samples = as.character(length(cs$samples)))
      samples <- if (length(cs$samples)) cs$samples else list(list())
      for (s in samples) {
        rows[[length(rows) + 1]] <- c(file_part, case_part,
                                      sample_submitter_id = gdc_value(s, "submitter_id"),
                                      vapply(sample_fields, function(f) gdc_value(s, f), character(1)))
      }
    }
  }
  if (!length(rows)) stop("The GDC returned no files")
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  out$file_size <- suppressWarnings(as.numeric(out$file_size))
  out$n_cases <- as.integer(out$n_cases)
  out$n_samples <- as.integer(out$n_samples)
  out
}

# Sample class from the GDC sample fields. Newer records use tissue_type (Tumor / Normal) and
# tumor_descriptor (Primary / Metastatic / Recurrence ...); older ones only sample_type
# ("Primary Tumor", "Solid Tissue Normal", ...). A tumour with no descriptor counts as primary
# (CPTAC-3 collected treatment-naive primary tumours); metastatic or recurrent ones never do.
classify_gdc_sample <- function(sample_type, tissue_type, tumor_descriptor) {
  norm <- function(x) {
    x <- tolower(trimws(as.character(x)))
    x[is.na(x)] <- ""
    x
  }
  st <- norm(sample_type)
  tt <- norm(tissue_type)
  td <- norm(tumor_descriptor)
  normal <- tt == "normal" | grepl("normal", st)
  tumour <- !normal & (tt == "tumor" | st %in% c("primary tumor", "tumor"))
  later <- td %in% c("metastatic", "recurrence", "progression") |
    grepl("metasta|recurr|new primary", st)
  out <- rep("other", length(st))
  out[normal] <- "normal"
  out[tumour & !later] <- "primary_tumour"
  out[tumour & later] <- "other_tumour"
  out
}

# One RNA-seq file per patient. A file is eligible when it comes from one case (no pooling across
# patients) and every sample behind it is a primary tumour. CPTAC-3 often makes one RNA aliquot
# from several pieces of the same tumour, so an eligible file can list several samples.
# If a patient has several eligible files, the first by sample ID and then file ID is kept
# (outcome-blind, like the first-barcode rule used for TCGA-LUAD).
choose_rna_files <- function(ft) {
  ft$class <- classify_gdc_sample(ft$sample_type, ft$tissue_type, ft$tumor_descriptor)
  files <- do.call(rbind, lapply(split(ft, ft$file_id), function(d) {
    ids <- sort(unique(stats::na.omit(d$sample_submitter_id)))
    data.frame(file_id = d$file_id[1], file_name = d$file_name[1], md5sum = d$md5sum[1],
               patient = d$patient[1], n_cases = max(d$n_cases, length(unique(d$patient))),
               n_samples = nrow(d), first_sample = if (length(ids)) ids[1] else NA_character_,
               samples = paste(ids, collapse = ";"),
               classes = paste(sort(unique(d$class)), collapse = "+"), stringsAsFactors = FALSE)
  }))
  rownames(files) <- NULL
  pooled <- files$n_cases > 1
  files$status <- ifelse(pooled, "pooled_across_patients",
                         ifelse(files$classes != "primary_tumour", paste0("not_eligible:", files$classes),
                                ifelse(files$n_samples > 1, "primary_tumour_several_samples",
                                       "primary_tumour")))
  tum <- files[files$status %in% c("primary_tumour", "primary_tumour_several_samples"), , drop = FALSE]
  tum <- tum[order(tum$patient, tum$first_sample, tum$file_id), , drop = FALSE]
  chosen <- tum[!duplicated(tum$patient), , drop = FALSE]
  flow <- c(files = nrow(files),
            pooled_across_patients = sum(pooled),
            not_primary_tumour = sum(!pooled & files$classes != "primary_tumour"),
            primary_tumour_files = nrow(tum),
            of_which_several_samples = sum(tum$n_samples > 1),
            patients_with_primary_tumour_file = nrow(chosen))
  list(chosen = chosen, flow = flow, files = files, all = ft)
}

# Number of records with a value for each field (to see which clinical fields a project uses).
gdc_field_counts <- function(records) {
  records <- Filter(function(r) is.list(r) && length(r) > 0 && !is.null(names(r)), records)
  if (!length(records)) return(data.frame(field = character(0), records_with_value = integer(0)))
  nm <- unlist(lapply(records, function(r) {
    names(r)[vapply(r, function(v) length(v) > 0 && !is.list(v), logical(1))]
  }))
  tab <- sort(table(nm), decreasing = TRUE)
  data.frame(field = names(tab), records_with_value = as.integer(tab), stringsAsFactors = FALSE)
}

# /cases records (expanded with demographic, diagnoses and follow_ups) -> one row per case.
# Column names follow TCGAbiolinks::GDCquery_clinic(), so make_clin_gdc() reads both alike.
# Last follow-up = the latest follow_ups$days_to_follow_up (as TCGAbiolinks does for TCGA);
# diagnoses$days_to_last_follow_up (older GDC model) only when no follow-up record has a value.
# Diagnosis fields come from the diagnoses flagged as the primary disease when there are flags.
# age_at_index holds the age in completed years, from the field named in age_source.
flatten_gdc_cases <- function(hits) {
  num <- function(x) suppressWarnings(as.numeric(x))
  rows <- lapply(hits, function(h) {
    dem <- h$demographic
    if (length(dem) && is.null(names(dem))) dem <- dem[[1]]
    dx <- if (is.null(h$diagnoses)) list() else h$diagnoses
    fu <- if (is.null(h$follow_ups)) list() else h$follow_ups
    is_primary <- vapply(dx, function(d) {
      identical(toupper(gdc_value(d, "diagnosis_is_primary_disease")), "TRUE")
    }, logical(1))
    dx_use <- if (any(is_primary)) dx[is_primary] else dx
    dx_first <- function(field) {
      v <- vapply(dx_use, gdc_value, character(1), name = field)
      v <- v[!is.na(v)]
      if (length(v)) v[[1]] else NA_character_
    }
    fu_days <- num(vapply(fu, gdc_value, character(1), name = "days_to_follow_up"))
    dx_days <- num(vapply(dx, gdc_value, character(1), name = "days_to_last_follow_up"))
    if (any(!is.na(fu_days))) {
      last_fu <- max(fu_days, na.rm = TRUE)
      fu_source <- "follow_ups"
    } else if (any(!is.na(dx_days))) {
      last_fu <- max(dx_days, na.rm = TRUE)
      fu_source <- "diagnoses"
    } else {
      last_fu <- NA_real_
      fu_source <- NA_character_
    }
    vital <- gdc_value(dem, "vital_status")
    if (is.na(vital)) vital <- dx_first("vital_status")
    death <- num(gdc_value(dem, "days_to_death"))
    if (is.na(death)) death <- num(dx_first("days_to_death"))
    sex <- gdc_value(dem, "sex_at_birth")
    if (is.na(sex)) sex <- gdc_value(dem, "gender")
    # Age in completed years: age_at_index, else age_at_diagnosis (days), else days_to_birth
    age <- num(gdc_value(dem, "age_at_index"))
    age_source <- if (is.na(age)) NA_character_ else "age_at_index"
    if (is.na(age) && !is.na(num(dx_first("age_at_diagnosis")))) {
      age <- floor(num(dx_first("age_at_diagnosis")) / 365.25)
      age_source <- "age_at_diagnosis"
    }
    if (is.na(age) && !is.na(num(gdc_value(dem, "days_to_birth")))) {
      age <- floor(-num(gdc_value(dem, "days_to_birth")) / 365.25)
      age_source <- "days_to_birth"
    }
    c(submitter_id = gdc_value(h, "submitter_id"),
      case_id = gdc_value(h, "case_id"),
      primary_site = gdc_value(h, "primary_site"),
      disease_type = gdc_value(h, "disease_type"),
      primary_diagnosis = dx_first("primary_diagnosis"),
      ajcc_pathologic_stage = dx_first("ajcc_pathologic_stage"),
      vital_status = vital,
      days_to_death = as.character(death),
      days_to_last_follow_up = as.character(last_fu),
      follow_up_source = fu_source,
      age_at_index = as.character(age),
      age_source = age_source,
      sex_at_birth = sex,
      n_diagnoses = as.character(length(dx)),
      n_follow_ups = as.character(length(fu)))
  })
  if (!length(rows)) stop("The GDC returned no cases")
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  for (v in c("days_to_death", "days_to_last_follow_up", "age_at_index")) out[[v]] <- num(out[[v]])
  for (v in c("n_diagnoses", "n_follow_ups")) out[[v]] <- as.integer(out[[v]])
  out
}

# Download GDC files into dir/<file_id>/<file_name>, checking every file's md5 checksum.
# Files already there with the right checksum are skipped, so a stopped run can be resumed.
gdc_download <- function(file_id, file_name, md5sum, dir, tries = 4,
                         base = "https://api.gdc.cancer.gov/data/") {
  stopifnot(length(file_id) == length(file_name), length(file_id) == length(md5sum))
  old <- options(timeout = max(900, getOption("timeout")))
  on.exit(options(old))
  paths <- file.path(dir, file_id, file_name)
  md5_ok <- function(p, m) {
    file.exists(p) && (is.na(m) || identical(unname(tools::md5sum(p)), m))
  }
  todo <- which(!mapply(md5_ok, paths, md5sum, USE.NAMES = FALSE))
  message("  ", length(paths) - length(todo), " of ", length(paths), " files already downloaded")
  for (k in seq_along(todo)) {
    i <- todo[k]
    dir.create(dirname(paths[i]), recursive = TRUE, showWarnings = FALSE)
    ok <- FALSE
    for (t in seq_len(tries)) {
      status <- tryCatch(utils::download.file(paste0(base, file_id[i]), paths[i],
                                              mode = "wb", quiet = TRUE),
                         error = function(e) 1L, warning = function(w) 1L)
      if (isTRUE(status == 0) && md5_ok(paths[i], md5sum[i])) {
        ok <- TRUE
        break
      }
      if (t < tries) Sys.sleep(10 * t)
    }
    if (!ok) {
      stop("Could not download ", file_name[i], " after ", tries, " tries. ",
           "Run the script again: files already downloaded are kept.")
    }
    if (k %% 10 == 0 || k == length(todo)) message("  downloaded ", k, " of ", length(todo))
  }
  paths
}

# tpm_unstranded of GDC STAR-Counts files -> matrix (gene_id x sample). Every file must list
# the same genes in the same order; the gene-model line of each file is kept as an attribute.
read_star_tpm <- function(paths, sample_names) {
  stopifnot(length(paths) > 0, length(paths) == length(sample_names), !anyDuplicated(sample_names))
  read_one <- function(p) {
    top <- readLines(p, n = 10)
    hdr <- strsplit(top[!startsWith(top, "#")][1], "\t", fixed = TRUE)[[1]]
    if (!all(c("gene_id", "tpm_unstranded") %in% hdr)) {
      stop(p, ": columns gene_id and tpm_unstranded not found")
    }
    cc <- ifelse(hdr == "gene_id", "character",
                 ifelse(hdr == "tpm_unstranded", "numeric", "NULL"))
    d <- utils::read.delim(p, comment.char = "#", colClasses = cc, quote = "")
    d <- d[grepl("^ENSG", d$gene_id), , drop = FALSE]
    list(gene_id = d$gene_id, tpm = d$tpm_unstranded,
         model = if (startsWith(top[1], "#")) top[1] else NA_character_)
  }
  first <- read_one(paths[1])
  tpm <- matrix(NA_real_, nrow = length(first$gene_id), ncol = length(paths),
                dimnames = list(first$gene_id, sample_names))
  tpm[, 1] <- first$tpm
  models <- first$model
  for (i in seq_along(paths)[-1]) {
    d <- read_one(paths[i])
    if (!identical(d$gene_id, first$gene_id)) stop(paths[i], ": genes differ from the first file")
    tpm[, i] <- d$tpm
    models <- c(models, d$model)
    if (i %% 25 == 0) message("  read ", i, " of ", length(paths))
  }
  attr(tpm, "gene_model") <- unique(models)
  tpm
}

# TPM matrix (gene_id x patient) + clinical table -> cohort object, with exactly the gene filter
# and transform of build_gdc_cohort(): ENSG rows without PAR_Y, log2(TPM + 1), Ensembl IDs
# without version.
build_star_cohort <- function(tpm, clin_raw, name, meta = list()) {
  stopifnot(is.matrix(tpm), !is.null(colnames(tpm)), !anyDuplicated(colnames(tpm)))
  model <- attr(tpm, "gene_model")
  ids <- rownames(tpm)
  keep <- grepl("^ENSG", ids) & !grepl("_PAR_Y$", ids)
  expr <- log2(tpm[keep, , drop = FALSE] + 1)
  rownames(expr) <- strip_version(ids[keep])
  meta <- c(list(source = "GDC STAR-Counts, tpm_unstranded, log2(TPM + 1)",
                 gene_model = if (length(model)) paste(sub("^#[[:space:]]*gene-model:[[:space:]]*", "",
                                                           model), collapse = "; ") else NA_character_),
            meta)
  finalise_cohort(expr, make_clin_gdc(clin_raw), name, platform = "RNA-seq", meta = meta)
}
