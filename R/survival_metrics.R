# survival_metrics.R - discrimination and effect size for a risk score
# (higher score = higher risk).

# Harrell's C-index.
cindex <- function(time, event, score) {
  unname(survival::concordance(survival::Surv(time, event) ~ score, reverse = TRUE)$concordance)
}

# Uno's C-index truncated at tau months (sensitivity for unequal follow-up).
cindex_uno <- function(time, event, score, tau = 60) {
  unname(survival::concordance(survival::Surv(time, event) ~ score, reverse = TRUE,
                               timewt = "n/G2", ymax = tau)$concordance)
}

# Hazard ratio per 1 SD of the score from a univariable Cox model.
hr_per_sd <- function(time, event, score) {
  if (stats::sd(score) < 1e-8) stop("Constant risk score: check the presence rule")
  s <- (score - mean(score)) / stats::sd(score)
  fit <- survival::coxph(survival::Surv(time, event) ~ s)
  b <- unname(stats::coef(fit))
  se <- sqrt(unname(stats::vcov(fit)[1, 1]))
  c(log_hr = b, se = se, hr = exp(b), lo = exp(b - 1.959964 * se), hi = exp(b + 1.959964 * se))
}

# Percentile bootstrap 95% CI for the C-index of a fixed score.
boot_cindex_ci <- function(time, event, score, B = 1000, seed = NULL, min_events = 5) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(time)
  v <- vapply(seq_len(B), function(b) {
    i <- sample.int(n, replace = TRUE)
    if (sum(event[i]) < min_events) return(NA_real_)
    cindex(time[i], event[i], score[i])
  }, numeric(1))
  unname(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))
}
