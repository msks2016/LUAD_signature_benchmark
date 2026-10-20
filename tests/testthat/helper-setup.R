# testthat runs helper-*.R files first: load the project functions and a small simulator.

find_project_root <- function(start = getwd()) {
  d <- normalizePath(start)
  for (i in 1:6) {
    if (file.exists(file.path(d, "config", "analysis.yaml"))) return(d)
    d <- dirname(d)
  }
  stop("Project root (folder containing config/analysis.yaml) not found")
}

for (f in list.files(file.path(find_project_root(), "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# A cohort with known truth: the first `signal` genes drive survival with log-HR `beta`
# per SD of their standardised sum.
sim_cohort <- function(n = 300, n_genes = 200, beta = 0.5, signal = 1:3, seed = 1) {
  set.seed(seed)
  expr <- matrix(stats::rnorm(n_genes * n, mean = 6, sd = 1), nrow = n_genes,
                 dimnames = list(sprintf("ENSG%011d", seq_len(n_genes)), sprintf("P%04d", seq_len(n))))
  z <- t(scale(t(expr[signal, , drop = FALSE])))
  u <- colSums(z) / sqrt(length(signal))
  t_event <- stats::rexp(n, rate = 0.02 * exp(beta * u))
  t_cens <- stats::runif(n, 0, 100)
  list(name = paste0("sim", seed), expr = expr, time = pmin(t_event, t_cens),
       event = as.integer(t_event <= t_cens))
}
