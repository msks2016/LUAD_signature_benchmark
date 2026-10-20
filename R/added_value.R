# added_value.R - H2: does the score add prognostic information beyond stage, sex and age?
# Registered model (Update 3): Model A = stage (four levels) + sex + age; Model B = Model A +
# scale(risk score); likelihood-ratio test on 1 df. Sex is used where recorded; age where
# recorded for >= 80% of patients (the caller passes age = NULL otherwise). Patients missing
# any covariate in the model are excluded from H2 only. B > 0 adds a percentile bootstrap
# 95% CI for the C-index difference (patients resampled, both models refitted).

h2_added_value <- function(time, event, score, stage, sex = NULL, age = NULL, B = 0, seed = NULL) {
  d <- data.frame(time = time, event = event, s = as.numeric(scale(score)),
                  stage = factor(as.character(stage), levels = c("I", "II", "III", "IV")))
  if (!is.null(sex)) d$sex <- sex_factor(sex)
  if (!is.null(age)) d$age <- suppressWarnings(as.numeric(age))
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d$stage <- droplevels(d$stage)
  if (!is.null(sex)) d$sex <- droplevels(d$sex)
  base_terms <- c(if (nlevels(d$stage) > 1) "stage",
                  if (!is.null(sex) && nlevels(d$sex) > 1) "sex",
                  if (!is.null(age)) "age")
  if (!length(base_terms)) stop("No clinical covariate varies in this cohort")
  f0 <- stats::as.formula(paste("survival::Surv(time, event) ~", paste(base_terms, collapse = " + ")))
  f1 <- stats::update(f0, . ~ . + s)
  fit_both <- function(dd) {
    suppressWarnings(list(m0 = survival::coxph(f0, data = dd), m1 = survival::coxph(f1, data = dd)))
  }
  m <- fit_both(d)
  lrt <- 2 * (m$m1$loglik[2] - m$m0$loglik[2])
  b <- unname(stats::coef(m$m1)["s"])
  se <- sqrt(stats::vcov(m$m1)["s", "s"])
  c0 <- unname(survival::concordance(m$m0)$concordance)
  c1 <- unname(survival::concordance(m$m1)$concordance)
  ci <- c(NA_real_, NA_real_)
  if (B > 0) {
    if (!is.null(seed)) set.seed(seed)
    v <- vapply(seq_len(B), function(i) {
      dd <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]
      r <- tryCatch(fit_both(dd), error = function(e) NULL)
      if (is.null(r)) return(NA_real_)
      unname(survival::concordance(r$m1)$concordance - survival::concordance(r$m0)$concordance)
    }, numeric(1))
    ci <- unname(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))
  }
  data.frame(n = nrow(d), events = sum(d$event), covariates = paste(base_terms, collapse = "+"),
             c_base = c0, c_full = c1, delta_c = c1 - c0, delta_c_low = ci[1], delta_c_high = ci[2],
             lrt_chisq = lrt, p_lrt = stats::pchisq(lrt, df = 1, lower.tail = FALSE),
             hr_adj = exp(b), hr_lo = exp(b - 1.959964 * se), hr_hi = exp(b + 1.959964 * se),
             stringsAsFactors = FALSE)
}
