# scoring.R - the only place a risk score is computed.
# Observed scores, null draws, bootstrap resamples and sensitivity analyses all call
# score_signature(), so the scale cannot drift between scripts.

# Per-gene z-scores within one cohort (rows = genes, columns = patients).
zscore_rows <- function(x) {
  mu <- rowMeans(x)
  sdv <- apply(x, 1, stats::sd)
  z <- (x - mu) / sdv
  z[sdv < 1e-8, ] <- NA_real_
  z
}

# Registered formula: risk = sum(beta_i * z(expr_i)), z computed within the cohort.
# scale = "raw" applies the coefficients to log2 expression (sensitivity analysis only).
score_signature <- function(expr, genes, coefs, scale = c("zscore", "raw")) {
  scale <- match.arg(scale)
  if (!is.matrix(expr)) expr <- as.matrix(expr)
  stopifnot(length(genes) == length(coefs), length(genes) > 0, !anyDuplicated(genes))
  absent <- setdiff(genes, rownames(expr))
  if (length(absent)) stop("Not in matrix: ", paste(absent, collapse = ", "))
  x <- expr[genes, , drop = FALSE]
  if (scale == "zscore") x <- zscore_rows(x)
  if (anyNA(x)) stop("Constant or missing values in a scored gene: check the presence rule")
  as.numeric(crossprod(coefs, x))
}
