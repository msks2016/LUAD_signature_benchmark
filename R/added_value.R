# added_value.R - H2: does the score add prognostic information beyond stage (and age)?

h2_added_value <- function(time, event, score, stage, age = NULL) {
  d <- data.frame(time = time, event = event, s = as.numeric(scale(score)), stage = stage)
  if (!is.null(age)) d$age <- age
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d$stage <- droplevels(as.factor(d$stage))
  base_terms <- c(if (nlevels(d$stage) > 1) "stage", if (!is.null(age)) "age")
  if (!length(base_terms)) stop("No clinical covariate varies in this cohort")
  f0 <- stats::as.formula(paste("survival::Surv(time, event) ~", paste(base_terms, collapse = " + ")))
  f1 <- stats::update(f0, . ~ . + s)
  m0 <- survival::coxph(f0, data = d)
  m1 <- survival::coxph(f1, data = d)
  lrt <- 2 * (m1$loglik[2] - m0$loglik[2])
  b <- unname(stats::coef(m1)["s"])
  se <- sqrt(stats::vcov(m1)["s", "s"])
  c0 <- unname(survival::concordance(m0)$concordance)
  c1 <- unname(survival::concordance(m1)$concordance)
  data.frame(n = nrow(d), events = sum(d$event), covariates = paste(base_terms, collapse = "+"),
             c_base = c0, c_full = c1, delta_c = c1 - c0, lrt_chisq = lrt,
             p_lrt = stats::pchisq(lrt, df = 1, lower.tail = FALSE),
             hr_adj = exp(b), hr_lo = exp(b - 1.959964 * se), hr_hi = exp(b + 1.959964 * se),
             stringsAsFactors = FALSE)
}
