# secondary.R - registered secondary metrics and controls (OSF registration and Update 3):
# Cox-model C-index (paired with the C-index of the published score), log-rank test at the
# median score, time-dependent AUC, proportional-hazards check with RMST, the stage positive
# control, gene-level sign concordance, signature overlap (Jaccard), direction consistency,
# and the power simulation. Functions only; the scripts call them.

# ---- Per-pair metrics --------------------------------------------------------------------

# C-index of the fitted Cox model Surv ~ scale(score). Reported beside the C-index of the
# published score (cindex()), which is the primary C-index (Update 3, clarification 10).
cindex_cox <- function(time, event, score) {
  s <- (score - mean(score)) / stats::sd(score)
  fit <- survival::coxph(survival::Surv(time, event) ~ s)
  unname(survival::concordance(fit)$concordance)
}

# Log-rank test comparing patients above the within-cohort median score with the rest.
logrank_median <- function(time, event, score) {
  high <- score > stats::median(score)
  if (length(unique(high)) < 2) return(c(chisq = NA_real_, p = NA_real_, n_high = sum(high)))
  fit <- survival::survdiff(survival::Surv(time, event) ~ high)
  chisq <- unname(fit$chisq)
  c(chisq = chisq, p = stats::pchisq(chisq, df = 1, lower.tail = FALSE), n_high = sum(high))
}

# Kaplan-Meier estimate of the censoring distribution G. Returns a function of u giving G(u),
# or the left limit G(u-) with left = TRUE. Deaths at a censoring time stay in the risk set.
censoring_km <- function(time, event) {
  cens_times <- sort(unique(time[event == 0]))
  surv <- numeric(length(cens_times))
  s <- 1
  for (k in seq_along(cens_times)) {
    u <- cens_times[k]
    s <- s * (1 - sum(time == u & event == 0) / sum(time >= u))
    surv[k] <- s
  }
  function(u, left = FALSE) {
    idx <- findInterval(u, cens_times, left.open = left)
    out <- rep(1, length(u))
    out[idx > 0] <- surv[idx[idx > 0]]
    out
  }
}

# Cumulative/dynamic time-dependent AUC at time t with inverse-probability-of-censoring
# weights (Uno et al. 2007; the "marginal" weighting of timeROC). Cases: died by t (weight
# 1 / G(T-)); controls: still followed after t (common weight, which cancels). Ties count 1/2.
td_auc <- function(time, event, score, t) {
  case <- time <= t & event == 1
  ctrl <- time > t
  if (!any(case) || !any(ctrl)) return(NA_real_)
  G <- censoring_km(time, event)
  g <- G(time[case], left = TRUE)
  if (any(g <= 0)) return(NA_real_)
  w <- 1 / g
  m_ctrl <- sort(score[ctrl])
  below <- findInterval(score[case], m_ctrl, left.open = TRUE)   # controls strictly lower
  ties <- findInterval(score[case], m_ctrl) - below
  sum(w * (below + 0.5 * ties)) / (sum(w) * length(m_ctrl))
}

# Percentile bootstrap 95% CI of td_auc (patients resampled).
boot_td_auc_ci <- function(time, event, score, t, B = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(time)
  v <- vapply(seq_len(B), function(b) {
    i <- sample.int(n, replace = TRUE)
    td_auc(time[i], event[i], score[i], t)
  }, numeric(1))
  if (all(is.na(v))) return(c(NA_real_, NA_real_))
  unname(stats::quantile(v, c(0.025, 0.975), na.rm = TRUE))
}

# Proportional-hazards test (cox.zph) of the Cox model Surv ~ scale(score).
ph_test_p <- function(time, event, score) {
  s <- (score - mean(score)) / stats::sd(score)
  fit <- survival::coxph(survival::Surv(time, event) ~ s)
  unname(survival::cox.zph(fit)$table["s", "p"])
}

# Restricted mean survival time, patients above vs at-or-below the median score. tau defaults
# to the smaller of the two groups' largest follow-up times (the usual default, as in survRM2).
# A negative difference means shorter survival in the high-score group (published direction).
rmst_median_split <- function(time, event, score, tau = NULL) {
  high <- score > stats::median(score)
  if (length(unique(high)) < 2) stop("The median split leaves one group empty")
  if (is.null(tau)) tau <- min(max(time[high]), max(time[!high]))
  one <- function(sel) {
    fit <- survival::survfit(survival::Surv(time[sel], event[sel]) ~ 1)
    tab <- summary(fit, rmean = tau)$table
    c(rmst = unname(tab[["rmean"]]), se = unname(tab[["se(rmean)"]]))
  }
  h <- one(high)
  l <- one(!high)
  d <- h[["rmst"]] - l[["rmst"]]
  se <- sqrt(h[["se"]]^2 + l[["se"]]^2)
  c(tau = tau, rmst_high = h[["rmst"]], rmst_low = l[["rmst"]], diff = d,
    lo = d - 1.959964 * se, hi = d + 1.959964 * se, p = 2 * stats::pnorm(-abs(d / se)))
}

# ---- Clinical covariates -----------------------------------------------------------------

# Sex as a factor with male as the reference level (values differ between cohorts).
sex_factor <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s[s %in% c("m", "male")] <- "male"
  s[s %in% c("f", "female")] <- "female"
  s[is.na(s) | s %in% c("", "na", "n/a", "not reported", "unknown")] <- NA
  other <- setdiff(sort(unique(stats::na.omit(s))), c("male", "female"))
  factor(s, levels = c("male", "female", other))
}

# ---- Stage positive control (registered; adapted in OSF Update 3, decision D1) -------------

# hr, p: named by stage level ("II", "III", "IV"), each versus stage I.
# Registered rule: HR(III vs I) > 1 with p < alpha and a monotonic HR across stages.
# Update 3 (D1): in a cohort with fewer than min_iii_iv patients in stage III or IV, the
# registered criterion cannot be evaluated and HR(II vs I) > 1 with p < alpha is used instead.
positive_control_decision <- function(hr, p, n_iii_iv, min_iii_iv = 10, alpha = 0.05) {
  present <- intersect(c("II", "III", "IV"), names(hr)[!is.na(hr)])
  if (n_iii_iv >= min_iii_iv) {
    rule <- "registered: HR(III vs I) > 1, p < 0.05, monotonic HR across stages"
    key <- "III"
    chain_levels <- present
  } else {
    rule <- "adapted (Update 3, D1): fewer than 10 stage III-IV patients; HR(II vs I) > 1, p < 0.05"
    key <- "II"
    chain_levels <- intersect("II", present)
  }
  if (!key %in% present) {
    return(list(rule = rule, key_level = key, monotonic = NA, pass = FALSE,
                note = paste0("HR(", key, " vs I) not estimable")))
  }
  monotonic <- all(diff(c(1, unname(hr[chain_levels]))) >= 0)
  pass <- isTRUE(hr[[key]] > 1 && p[[key]] < alpha && monotonic)
  list(rule = rule, key_level = key, monotonic = monotonic, pass = pass, note = "")
}

# Cox model on stage (factor, stage I = reference) in one cohort, and the gate decision.
positive_control <- function(time, event, stage, min_iii_iv = 10, alpha = 0.05) {
  d <- data.frame(time = time, event = event,
                  stage = factor(as.character(stage), levels = c("I", "II", "III", "IV")))
  d <- d[!is.na(d$stage), , drop = FALSE]
  d$stage <- droplevels(d$stage)
  if (!"I" %in% levels(d$stage)) stop("No stage I patients: the reference level is missing")
  warn <- character(0)
  fit <- withCallingHandlers(
    survival::coxph(survival::Surv(time, event) ~ stage, data = d),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  b <- stats::coef(fit)
  se <- sqrt(diag(stats::vcov(fit)))
  lev <- sub("^stage", "", names(b))
  hr <- stats::setNames(exp(b), lev)
  p <- stats::setNames(2 * stats::pnorm(-abs(b / se)), lev)
  n_by <- table(d$stage)
  n_iii_iv <- sum(d$stage %in% c("III", "IV"))
  dec <- positive_control_decision(hr, p, n_iii_iv, min_iii_iv, alpha)
  get <- function(x, l) if (l %in% names(x)) unname(x[[l]]) else NA_real_
  data.frame(patients = nrow(d), deaths = sum(d$event),
             n_by_stage = paste(names(n_by), as.integer(n_by), sep = ":", collapse = " "),
             n_iii_iv = n_iii_iv, rule = dec$rule,
             hr_ii = get(hr, "II"), p_ii = get(p, "II"), hr_iii = get(hr, "III"), p_iii = get(p, "III"),
             hr_iv = get(hr, "IV"), p_iv = get(p, "IV"),
             c_stage = unname(survival::concordance(fit)$concordance),
             monotonic = dec$monotonic, pass = dec$pass,
             note = paste(c(if (nzchar(dec$note)) dec$note, unique(warn)), collapse = "; "),
             stringsAsFactors = FALSE)
}

# ---- Random-gene pool of the H1 null test (Update 3, decision D3) -------------------------

# Measurable genes of the given GENCODE v36 types, plus the signature genes themselves (they
# are needed to place each signature gene in its expression decile but are never drawn).
null_pool <- function(measurable_ids, gene_table, types, signature_ids) {
  if (is.null(types)) return(union(measurable_ids, signature_ids))
  typed <- gene_table$gene_id[gene_table$gene_type %in% types]
  union(intersect(measurable_ids, typed), signature_ids)
}

# ---- Gene-level sign concordance -------------------------------------------------------------

# Univariable Cox coefficient of each gene's z-scored expression.
gene_cox <- function(expr, genes, time, event) {
  out <- lapply(genes, function(g) {
    z <- as.numeric(scale(expr[g, ]))
    fit <- survival::coxph(survival::Surv(time, event) ~ z)
    b <- unname(stats::coef(fit))
    se <- sqrt(stats::vcov(fit)[1, 1])
    data.frame(ensembl_id = g, cox_beta = b, cox_p = 2 * stats::pnorm(-abs(b / se)),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# ---- Signature overlap and direction consistency -----------------------------------------

# Jaccard index between every pair of gene sets (named list of character vectors).
jaccard_matrix <- function(sets) {
  n <- names(sets)
  m <- matrix(NA_real_, length(n), length(n), dimnames = list(n, n))
  for (i in seq_along(n)) {
    for (j in seq_along(n)) {
      m[i, j] <- length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]]))
    }
  }
  m
}

# Clusters = connected groups of signatures linked by Jaccard > threshold (registered: 0.3).
jaccard_clusters <- function(J, threshold = 0.3) {
  n <- nrow(J)
  cl <- rep(NA_integer_, n)
  k <- 0L
  for (i in seq_len(n)) {
    if (!is.na(cl[i])) next
    k <- k + 1L
    cl[i] <- k
    queue <- i
    while (length(queue)) {
      v <- queue[1]
      queue <- queue[-1]
      nb <- which(J[v, ] > threshold & is.na(cl))
      cl[nb] <- k
      queue <- c(queue, nb)
    }
  }
  stats::setNames(cl, rownames(J))
}

# Direction consistency (secondary exploratory outcome, Update 1). pairs: one row per
# independent evaluable signature x cohort pair with log_hr (HR per SD of the published score).
# Consistent = log_hr > 0 (higher published score, higher hazard). One observation per cluster
# per cohort: the mean log_hr of the cluster's signatures in that cohort. Two-sided binomial
# sign test against 50%, at the nominal (pair) and the effective (cluster x cohort) level.
direction_consistency <- function(pairs, clusters) {
  if (!nrow(pairs)) stop("No independent pairs")
  if (!all(pairs$signature %in% names(clusters))) stop("Some signatures have no cluster")
  pairs$cluster <- unname(clusters[pairs$signature])
  obs <- stats::aggregate(log_hr ~ cluster + cohort, data = pairs, FUN = mean)
  obs$n_signatures <- stats::aggregate(log_hr ~ cluster + cohort, data = pairs, FUN = length)$log_hr
  test <- function(x, level) {
    k <- sum(x > 0)
    bt <- stats::binom.test(k, length(x), 0.5)
    data.frame(level = level, n = length(x), consistent = k, proportion = k / length(x),
               p_two_sided = bt$p.value, stringsAsFactors = FALSE)
  }
  list(observations = obs,
       summary = rbind(test(obs$log_hr, "effective: one per cluster per cohort"),
                       test(pairs$log_hr, "nominal: every signature x cohort pair")))
}

# ---- Power simulation -------------------------------------------------------------------------

# True (uncensored) C-index when the hazard is proportional to exp(beta * x), x ~ N(0, 1):
# C = 2 E[1(D > 0) plogis(beta D)], D = x1 - x2 ~ N(0, 2).
c_from_beta <- function(beta) {
  if (beta == 0) return(0.5)
  f <- function(d) stats::plogis(beta * d) * stats::dnorm(d, sd = sqrt(2))
  2 * stats::integrate(f, 0, Inf)$value
}

beta_for_c <- function(target_c) {
  if (target_c <= 0.5 || target_c >= 1) stop("target_c must be between 0.5 and 1")
  stats::uniroot(function(b) c_from_beta(b) - target_c, c(1e-6, 20), tol = 1e-10)$root
}

# Exponential censoring rate giving the wanted share of deaths (event hazard exp(beta * x)).
# P(death | x) = exp(beta x) / (exp(beta x) + mu) = plogis(beta x - log mu), written with plogis
# so that large x cannot overflow exp() (which stopped integrate() for C = 0.70).
censoring_rate_for <- function(beta, event_fraction) {
  if (event_fraction <= 0 || event_fraction >= 1) stop("event_fraction must be between 0 and 1")
  ef <- function(log_mu) {
    stats::integrate(function(x) stats::plogis(beta * x - log_mu) * stats::dnorm(x),
                     -Inf, Inf)$value
  }
  exp(stats::uniroot(function(l) ef(l) - event_fraction, c(-15, 15), tol = 1e-10)$root)
}

simulate_ph_cohort <- function(n, beta, censoring_rate) {
  x <- stats::rnorm(n)
  t_event <- stats::rexp(n, rate = exp(beta * x))
  t_cens <- stats::rexp(n, rate = censoring_rate)
  list(time = pmin(t_event, t_cens), event = as.integer(t_event <= t_cens), score = x)
}

# Two-sided test of C = 0.5 from the variance survival::concordance() reports.
cindex_test_p <- function(time, event, score) {
  cc <- survival::concordance(survival::Surv(time, event) ~ score, reverse = TRUE)
  z <- (cc$concordance - 0.5) / sqrt(cc$var)
  unname(2 * stats::pnorm(-abs(z)))
}

# Share of simulated cohorts (n patients, the cohort's share of deaths) in which a score with
# true C-index target_c is detected at each alpha. The H1 test compares with random genes,
# whose null spread is wider than chance alone, so this power is an upper bound for H1.
power_simulation <- function(n, deaths, target_c, alphas = 0.05, B = 1000, seed = NULL,
                             test_fun = cindex_test_p) {
  if (!is.null(seed)) set.seed(seed)
  beta <- beta_for_c(target_c)
  mu <- censoring_rate_for(beta, deaths / n)
  p <- vapply(seq_len(B), function(b) {
    d <- simulate_ph_cohort(n, beta, mu)
    if (sum(d$event) < 2) return(NA_real_)
    test_fun(d$time, d$event, d$score)
  }, numeric(1))
  data.frame(n = n, deaths = deaths, target_c = target_c, beta = beta, alpha = alphas,
             power = vapply(alphas, function(a) mean(p < a, na.rm = TRUE), numeric(1)),
             replicates = sum(!is.na(p)), stringsAsFactors = FALSE)
}

# ---- Helpers for the analysis scripts -------------------------------------------------------

# Positive-control status of each cohort (written by scripts/29_positive_control.R).
cohort_gate <- function(file = "results/tables/positive_control.csv") {
  if (!file.exists(file)) stop("Run scripts/29_positive_control.R first")
  pc <- utils::read.csv(file, stringsAsFactors = FALSE)
  stats::setNames(ifelse(pc$pass %in% TRUE, "valid", "exploratory: stage control failed"), pc$cohort)
}

# Whether age enters H2 in a cohort (registered: recorded for >= min_recorded of patients).
use_age_in_h2 <- function(age, min_recorded = 0.8) mean(!is.na(age)) >= min_recorded
