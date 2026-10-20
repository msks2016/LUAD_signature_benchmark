# utils.R - small helpers used across the project.
# Files in R/ only define functions: nothing runs when they are sourced.

load_config <- function(path = "config/analysis.yaml") {
  yaml::read_yaml(path)
}

# Load every function file in R/ (each script starts with this line).
load_functions <- function(dir = "R") {
  for (f in list.files(dir, pattern = "\\.R$", full.names = TRUE)) source(f)
  invisible(TRUE)
}

# Record one number for the manuscript in results/numbers.csv.
# Re-running a script replaces its own keys, so the file always holds the latest run.
record_number <- function(key, value, note = "", file = "results/numbers.csv") {
  stopifnot(length(key) == 1, length(value) == 1)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  row <- data.frame(key = key, value = as.character(value), note = note,
                    written = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                    stringsAsFactors = FALSE)
  if (file.exists(file)) {
    old <- utils::read.csv(file, stringsAsFactors = FALSE, colClasses = "character")
    row <- rbind(old[old$key != key, , drop = FALSE], row)
  }
  utils::write.csv(row, file, row.names = FALSE)
  invisible(value)
}

fmt_num <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

fmt_hr <- function(hr, lo, hi, digits = 2) {
  sprintf("%s (95%% CI %s-%s)", fmt_num(hr, digits), fmt_num(lo, digits), fmt_num(hi, digits))
}

# Collapse AJCC stage text ("Stage IA1", "IIB", "stage IIIA", "IV") to I / II / III-IV.
collapse_stage <- function(x) {
  s <- toupper(trimws(as.character(x)))
  s <- sub("^STAGE[[:space:]]*", "", s)
  out <- rep(NA_character_, length(s))
  out[grepl("^IV", s)] <- "III-IV"
  out[is.na(out) & grepl("^III", s)] <- "III-IV"
  out[is.na(out) & grepl("^II", s)] <- "II"
  out[is.na(out) & grepl("^I", s)] <- "I"
  factor(out, levels = c("I", "II", "III-IV"))
}

# Median follow-up by the reverse Kaplan-Meier method (censoring treated as the event).
reverse_km_median <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, 1 - event) ~ 1)
  unname(summary(fit)$table["median"])
}

# Read the cohorts named in `cohort_names` from data/processed/<name>.rds.
load_cohorts <- function(cohort_names, dir = "data/processed") {
  out <- lapply(cohort_names, function(n) readRDS(file.path(dir, paste0(n, ".rds"))))
  stats::setNames(out, cohort_names)
}

# Checks every cohort object must pass (run in scripts/19_cohort_qc.R).
check_cohort <- function(co, max_log2 = 20, min_max_log2 = 6) {
  if (!is.matrix(co$expr)) stop(co$name, ": expr must be a numeric matrix")
  if (!identical(colnames(co$expr), co$clin$patient_id)) stop(co$name, ": expr columns and clin rows differ")
  if (anyDuplicated(co$clin$patient_id)) stop(co$name, ": duplicated patients")
  if (any(co$time <= 0)) stop(co$name, ": follow-up must be > 0")
  if (!all(co$event %in% c(0L, 1L))) stop(co$name, ": event must be 0/1")
  if (anyNA(co$expr)) stop(co$name, ": missing expression values")
  rng <- range(co$expr)
  if (rng[2] > max_log2) stop(co$name, ": values up to ", round(rng[2]), " - not on a log2 scale")
  if (rng[2] < min_max_log2) stop(co$name, ": values only up to ", round(rng[2], 1), " - looks double-logged")
  invisible(TRUE)
}

# Turn an aligned expression matrix and clinical table into the standard cohort object.
finalise_cohort <- function(expr, clin, name, platform, meta = list()) {
  clin <- clin[!is.na(clin$patient_id) & !duplicated(clin$patient_id), , drop = FALSE]
  common <- intersect(colnames(expr), clin$patient_id)
  flow <- c(samples_with_expression = ncol(expr), with_clinical_record = length(common))
  clin <- clin[match(common, clin$patient_id), , drop = FALSE]
  ok <- !is.na(clin$os_months) & clin$os_months > 0 & !is.na(clin$os_event)
  flow["with_valid_survival"] <- sum(ok)
  clin <- clin[ok, , drop = FALSE]
  expr <- as.matrix(expr[, clin$patient_id, drop = FALSE])
  meta$platform <- platform
  meta$sample_flow <- flow
  meta$built <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  list(name = name, expr = expr, clin = clin,
       time = clin$os_months, event = as.integer(clin$os_event), meta = meta)
}
