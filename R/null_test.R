# null_test.R - H1: is a signature better than random genes of similar expression?

# Each signature gene is replaced by a random measurable gene from the same decile of mean
# expression (drawn without replacement within a draw); the published coefficients are
# applied unchanged and scored exactly like the real signature.
expression_matched_null <- function(expr, genes, coefs, time, event, pool,
                                    scale = "zscore", n_draws = 2000, n_bins = 10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!all(genes %in% pool)) stop("Every signature gene must be in the pool of measurable genes")
  mu <- rowMeans(expr[pool, , drop = FALSE])
  breaks <- unique(stats::quantile(mu, probs = seq(0, 1, length.out = n_bins + 1)))
  bin <- cut(mu, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  names(bin) <- pool
  others <- setdiff(pool, genes)
  candidates <- split(others, bin[others])
  gene_bin <- bin[genes]

  draw <- function() {
    g <- character(length(genes))
    for (b in unique(gene_bin)) {
      idx <- which(gene_bin == b)
      cand <- candidates[[as.character(b)]]
      g[idx] <- cand[sample.int(length(cand), length(idx))]
    }
    g
  }

  observed <- cindex(time, event, score_signature(expr, genes, coefs, scale))
  null <- vapply(seq_len(n_draws), function(i) {
    cindex(time, event, score_signature(expr, draw(), coefs, scale))
  }, numeric(1))
  p <- (1 + sum(abs(null - 0.5) >= abs(observed - 0.5))) / (n_draws + 1)
  list(observed = observed, p_two_sided = p,
       null_q025 = unname(stats::quantile(null, 0.025)),
       null_q975 = unname(stats::quantile(null, 0.975)),
       null = null)
}

# Registered three-way classification (S21, excerpt 1).
classify_h1 <- function(c_index, p_adjusted, alpha = 0.05) {
  ifelse(p_adjusted < alpha & c_index > 0.5, "Replicated",
         ifelse(p_adjusted < alpha & c_index < 0.5, "Contradicted", "Inconclusive"))
}
