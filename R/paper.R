# paper.R - helpers for the manuscript tables and figures (scripts/40_paper_outputs.R).
# Base graphics only, so the figures need no extra package.

# "Paper62_Yao2021" -> "Yao 2021"
sig_label <- function(id) {
  x <- sub("^Paper[0-9]+_", "", id)
  sub("^([A-Za-z]+)([0-9]{4})$", "\\1 \\2", x)
}

# Estimate with 95% CI as text: "0.65 (0.58-0.73)"; sep = " to " for quantities that can be
# negative (delta C-index, log-HR, RMST difference).
fmt_ci <- function(est, lo, hi, digits = 2, sep = "-") {
  f <- function(x) formatC(x, format = "f", digits = digits)
  ifelse(is.na(est), "", paste0(f(est), " (", f(lo), sep, f(hi), ")"))
}

# p-values as text: "<0.001"; 3 decimals below 0.1 (e.g. 0.056); 2 decimals otherwise.
fmt_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001",
         ifelse(p < 0.1, formatC(p, format = "f", digits = 3), formatC(p, format = "f", digits = 2))))
}

# "p = 0.012" or "p < 0.001"
p_text <- function(p, prefix = "p") {
  x <- fmt_p(p)
  ifelse(substr(x, 1, 1) == "<", paste(prefix, "<", substring(x, 2)), paste(prefix, "=", x))
}

# Colour-blind-safe cohort colours (Okabe-Ito); TCGA-LUAD (reference) in grey.
cohort_colours <- function() {
  c("TCGA-LUAD" = "#7F7F7F", "CPTAC-3" = "#0072B2", "GSE31210" = "#E69F00",
    "GSE50081" = "#009E73", "GSE30219" = "#CC79A7", "GSE81089" = "#D55E00", "GSE72094" = "#56B4E9")
}

# Write one figure to PDF and PNG with the same drawing function.
save_figure <- function(path_stem, draw, width = 7, height = 5) {
  dir.create(dirname(path_stem), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(paste0(path_stem, ".pdf"), width = width, height = height)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(paste0(path_stem, ".png"), width = width, height = height, units = "in", res = 300)
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(path_stem)
}

# Forest plot (base graphics). rows: data.frame with label, est, lo, hi (ratio scale), and
# optional col, pch, cex, bold, right (text printed at the right margin). log axis.
forest_plot <- function(rows, xlab = "Hazard ratio per SD (95% CI)", ref = 1, xlim = NULL,
                        main = NULL, extra_lines = NULL, extra_labels = NULL, right_title = NULL,
                        left = 12) {
  n <- nrow(rows)
  y <- rev(seq_len(n))
  if (is.null(rows$col)) rows$col <- "black"
  if (is.null(rows$pch)) rows$pch <- 15
  if (is.null(rows$cex)) rows$cex <- 1
  if (is.null(rows$bold)) rows$bold <- FALSE
  ok <- !is.na(rows$est)
  if (is.null(xlim)) xlim <- range(c(rows$lo[ok], rows$hi[ok], ref)) * c(0.9, 1.1)
  op <- graphics::par(mar = c(4, left, if (is.null(main)) 1.5 else 3, if (is.null(rows$right)) 1 else 11))
  on.exit(graphics::par(op))
  graphics::plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5), log = "x", yaxt = "n", xlab = xlab,
                 ylab = "", main = main, bty = "n")
  graphics::abline(v = ref, lty = 2, col = "grey40")
  for (i in seq_along(extra_lines)) {
    graphics::abline(v = extra_lines[i], lty = 3, col = "grey60")
    if (!is.null(extra_labels)) graphics::text(extra_lines[i], n + 0.5, extra_labels[i], cex = 0.65,
                                               col = "grey40", pos = 4, xpd = TRUE)
  }
  graphics::segments(rows$lo[ok], y[ok], rows$hi[ok], y[ok], col = rows$col[ok], lwd = 1.5)
  graphics::points(rows$est[ok], y[ok], pch = rows$pch[ok], col = rows$col[ok], cex = rows$cex[ok])
  graphics::axis(2, at = y, labels = FALSE, tick = FALSE)
  graphics::mtext(rows$label, side = 2, at = y, las = 1, line = 0.5, cex = 0.75,
                  font = ifelse(rows$bold, 2, 1))
  if (!is.null(rows$right)) {
    graphics::mtext(rows$right, side = 4, at = y, las = 1, line = 0.5, cex = 0.7,
                    font = ifelse(rows$bold, 2, 1))
    if (!is.null(right_title)) graphics::mtext(right_title, side = 4, at = n + 0.9, las = 1,
                                               line = 0.5, cex = 0.7, font = 2)
  }
  invisible(y)
}

# Kaplan-Meier curves for one score split at its median (needs survival).
km_median_split <- function(time, event, score, main = "", col = c("#0072B2", "#D55E00")) {
  high <- factor(ifelse(score > stats::median(score), "High risk", "Low risk"),
                 levels = c("Low risk", "High risk"))
  fit <- survival::survfit(survival::Surv(time, event) ~ high)
  lr <- survival::survdiff(survival::Surv(time, event) ~ high)
  p <- stats::pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  graphics::plot(fit, col = col, lwd = 1.5, xlab = "Months", ylab = "Overall survival",
                 main = main, cex.main = 0.85, mark.time = TRUE, cex = 0.5)
  graphics::legend("bottomleft", bty = "n", cex = 0.7, lwd = 1.5, col = col,
                   legend = sprintf("%s (n = %d)", levels(high), as.integer(table(high))))
  graphics::legend("topright", bty = "n", cex = 0.7, legend = paste("log-rank p =", fmt_p(p)))
  invisible(p)
}
