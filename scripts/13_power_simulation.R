# scripts/13_power_simulation.R
# Simulates statistical power to detect a truly prognostic signature
# at n=150 patients / 38 events (CPTAC-3 sample size).
# Purpose: support the interpretation of "Inconclusive" signatures as
# underpowered rather than non-prognostic.

suppressPackageStartupMessages(library(survival))
set.seed(2026)

sim_cohort <- function(n = 150, events = 38, true_beta = 0.4) {
  score   <- rnorm(n)
  rate    <- exp(true_beta * score) * 0.02
  t_event <- rexp(n, rate = rate)
  t_cens  <- runif(n, 0, 60)
  time    <- pmin(t_event, t_cens)
  event   <- as.integer(t_event <= t_cens)
  
  # force exactly `events` events
  if (sum(event) < events) {
    need <- events - sum(event)
    idx  <- order(time * (!event))[1:need]
    event[idx] <- 1L
  } else if (sum(event) > events) {
    excess <- sum(event) - events
    idx    <- order(-time * event)[1:excess]
    event[idx] <- 0L
  }
  data.frame(score = score, time = time, event = event)
}

beta_to_c <- function(b) pnorm(b / sqrt(2 + b^2))

run_power <- function(true_beta, B = 500) {
  p_vals <- replicate(B, {
    d   <- sim_cohort(true_beta = true_beta)
    fit <- coxph(Surv(time, event) ~ score, data = d)
    summary(fit)$coefficients["score", "Pr(>|z|)"]
  })
  mean(p_vals < 0.05)
}

cat("Power simulation (n=150, events=38, B=500):\n\n")
for (b in c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) {
  p <- run_power(b)
  cat(sprintf("  true C-index=%.2f  (beta=%.2f)  power=%.2f\n",
              beta_to_c(b), b, p))
}