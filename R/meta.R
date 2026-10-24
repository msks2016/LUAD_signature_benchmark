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

# Confidence-interval variants for one pooled estimate (sensitivity to the registered
# REML + Hartung-Knapp interval). With little observed heterogeneity (q = Q/(k-1) < 1) the
# Hartung-Knapp interval is narrower than the common-effect interval; the modified version
# (q* = max(1, q); Knapp & Hartung 2003, Roever et al. 2015) never is. tau2: the REML estimate
# (computed with metafor when not given).
pool_ci_variants <- function(log_hr, se, tau2 = NULL, level = 0.95) {
  k <- length(log_hr)
  if (k < 3) stop("Interval variants are reported for k >= 3 only")
  if (is.null(tau2)) tau2 <- metafor::rma(yi = log_hr, sei = se, method = "REML")$tau2
  w <- 1 / (se^2 + tau2)
  mu <- sum(w * log_hr) / sum(w)
  q <- sum(w * (log_hr - mu)^2) / (k - 1)
  zc <- stats::qnorm(1 - (1 - level) / 2)
  tc <- stats::qt(1 - (1 - level) / 2, df = k - 1)
  w0 <- 1 / se^2
  mu0 <- sum(w0 * log_hr) / sum(w0)
  half <- c(common_effect = zc * sqrt(1 / sum(w0)), reml_wald = zc * sqrt(1 / sum(w)),
            reml_hk = tc * sqrt(q / sum(w)), reml_hk_modified = tc * sqrt(max(1, q) / sum(w)))
  centre <- c(common_effect = mu0, reml_wald = mu, reml_hk = mu, reml_hk_modified = mu)
  data.frame(interval = names(half), k = k, tau2 = tau2, q = q, log_hr = unname(centre),
             ci_low = unname(centre - half), ci_high = unname(centre + half),
             hr = exp(unname(centre)), hr_low = exp(unname(centre - half)),
             hr_high = exp(unname(centre + half)), stringsAsFactors = FALSE)
}
