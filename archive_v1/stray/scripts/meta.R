# meta.R - pooling one signature's log-HRs across independent cohorts (H3, exploratory).

pool_log_hr <- function(log_hr, se, labels = NULL) {
  if (is.null(labels)) labels <- paste0("study", seq_along(log_hr))   # metafor 5 refuses slab = NULL
  k <- length(log_hr)
  if (k < 2) stop("Need at least two cohorts")
  if (k == 2) {
    # Heterogeneity cannot be estimated from two studies: common-effect, descriptive only.
    fit <- metafor::rma(yi = log_hr, sei = se, method = "EE", slab = labels)
    model <- "common-effect (k = 2, descriptive)"
  } else {
    fit <- metafor::rma(yi = log_hr, sei = se, method = "REML", test = "knha", slab = labels)
    model <- "REML + Hartung-Knapp"
  }
  b <- as.numeric(fit$b)
  data.frame(k = k, model = model, log_hr = b, ci_low = fit$ci.lb, ci_high = fit$ci.ub,
             hr = exp(b), hr_low = exp(fit$ci.lb), hr_high = exp(fit$ci.ub),
             i2 = if (is.null(fit$I2)) NA_real_ else fit$I2,
             tau2 = if (is.null(fit$tau2)) NA_real_ else fit$tau2,
             stringsAsFactors = FALSE)
}

pool_leave_one_out <- function(log_hr, se, labels) {
  out <- lapply(seq_along(log_hr), function(i) {
    r <- pool_log_hr(log_hr[-i], se[-i], labels[-i])
    r$left_out <- labels[i]
    r
  })
  do.call(rbind, out)
}
