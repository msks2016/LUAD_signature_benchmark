# 43_graphical_abstract.R - graphical abstract for the journal submission. Drawn with base R
# graphics only from the v2 result tables (no image generator): left, gene coverage (as Figure 2);
# centre, aggregate transportability across settings (as Figure 6); right, the H1 result.
# Labels and numbers are computed from the tables, never typed.
# Size: 7 x 2.8 inches at 400 dpi (2800 x 1120 pixels; Elsevier minimum 1328 x 531, same 2.5:1 ratio),
# with large text so that it stays readable when shown 13 cm wide.
# Output: results/paper/figures/graphical_abstract.pdf, .png and .tif

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
rd <- function(f) {
  p <- file.path("results/tables", f)
  if (!file.exists(p)) stop("Missing ", p, ": run scripts 23, 30, 33 and 34 first")
  utils::read.csv(p, stringsAsFactors = FALSE)
}
cov <- rd("coverage_long.csv")
h1 <- rd("h1_null_test.csv")
tr <- rd("transportability_pooled.csv")
str_ <- rd("sensitivity_transportability.csv")
rc <- if (file.exists("results/tables/transportability_robust_crosscheck.csv")) rd("transportability_robust_crosscheck.csv") else NULL

cohorts <- c(cfg$cohorts$reference, cfg$cohorts$independent)
sigs <- sort(unique(cov$signature))
n_beat <- sum(h1$classification == "Replicated")
hr <- exp(tr$pooled_log_hr); hr_lo <- exp(tr$ci_low); hr_hi <- exp(tr$ci_high)
label_left <- sprintf("%d signatures, %d cohorts", length(sigs), length(cohorts))
label_mid <- sprintf("Shared signal HR %.2f", hr)
label_right <- if (n_beat == 0) "None beats random gene sets" else sprintf("%d of %d beat random gene sets", n_beat, nrow(h1))

navy <- "#1F3A5F"; blue <- "#2E6DA4"; grey <- "grey45"

draw_ga <- function() {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::layout(matrix(1:3, nrow = 1), widths = c(1.05, 1.45, 0.95))
  graphics::par(family = "sans", oma = c(0, 0, 0, 0))

  # Left: share of each signature's genes measurable in each cohort (frame: >= 80%, evaluable)
  graphics::par(mar = c(4.6, 2.2, 4.0, 0.6))
  short <- sub("^GSE", "", sub("TCGA-LUAD", "TCGA", cohorts))
  m <- matrix(NA_real_, length(sigs), length(cohorts))
  for (i in seq_along(sigs)) for (j in seq_along(cohorts)) {
    v <- cov$frac_present[cov$signature == sigs[i] & cov$cohort == cohorts[j]]
    if (length(v)) m[i, j] <- v[1]
  }
  pal <- grDevices::colorRampPalette(c("#EEF2F6", "#C9DAEA", blue))(101)
  graphics::plot(NA, xlim = c(0.5, length(cohorts) + 0.5), ylim = c(0.5, length(sigs) + 0.5), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_along(sigs)) for (j in seq_along(cohorts)) {
    y <- length(sigs) + 1 - i
    v <- m[i, j]
    graphics::rect(j - 0.46, y - 0.42, j + 0.46, y + 0.42, border = NA,
                   col = if (is.na(v)) "grey95" else pal[round(100 * v) + 1])
    if (!is.na(v) && v >= cfg$coverage$primary_min_fraction - 1e-9)
      graphics::rect(j - 0.46, y - 0.42, j + 0.46, y + 0.42, border = navy, lwd = 1.6)
  }
  graphics::axis(1, at = seq_along(cohorts), labels = short, las = 2, tick = FALSE, line = -0.6, cex.axis = 1.05)
  graphics::mtext("signatures", side = 2, line = 0.6, cex = 0.85, col = grey)
  graphics::mtext(label_left, side = 3, line = 1.6, font = 2, cex = 1.15, col = navy)
  graphics::mtext("framed: >= 80% of genes measured", side = 3, line = 0.4, cex = 0.78, col = grey)

  # Centre: mean HR per SD over the independent pairs, primary and sensitivity settings
  graphics::par(mar = c(4.6, 10.2, 4.0, 1.0))
  lab <- c("primary (repeated with B = 1000)" = "primary, repeated", "raw log2 scale" = "raw log2 scale",
           "100% coverage" = "100% coverage", "70% coverage" = "70% coverage",
           "RNA-seq TPM rescaled to protein-coding + lncRNA genes" = "TPM rescaled", "stage I-II only" = "stage I-II only",
           "deaths within 3 months excluded" = "early deaths excluded", "GSE31210 without the 22 flagged cases" = "GSE31210 minus 22")
  s <- str_[!grepl("^primary", str_$setting), ]
  rows <- data.frame(label = c(sprintf("Primary (%d pairs)", tr$n_pairs), unname(ifelse(s$setting %in% names(lab), lab[s$setting], s$setting))),
                     est = c(hr, exp(s$pooled_log_hr)), lo = c(hr_lo, exp(s$ci_low)), hi = c(hr_hi, exp(s$ci_high)),
                     col = c(navy, rep(grey, nrow(s))), pch = c(15, rep(15, nrow(s))), cex = c(2.2, rep(1.3, nrow(s))),
                     stringsAsFactors = FALSE)
  if (!is.null(rc)) rows <- rbind(rows, data.frame(label = "multilevel cross-check", est = exp(rc$log_hr), lo = exp(rc$ci_low),
                                                    hi = exp(rc$ci_high), col = grey, pch = 5, cex = 1.3))
  n <- nrow(rows); y <- rev(seq_len(n))
  graphics::plot(NA, xlim = c(0.9, 1.7), ylim = c(0.5, n + 0.5), log = "x", axes = FALSE, xlab = "", ylab = "")
  graphics::abline(v = 1, lty = 2, col = "grey55")
  graphics::abline(v = exp(cfg$transportability$delta), lty = 3, col = "grey55")
  graphics::segments(rows$lo, y, rows$hi, y, col = rows$col, lwd = ifelse(seq_len(n) == 1, 3, 1.6))
  graphics::points(rows$est, y, pch = rows$pch, col = rows$col, cex = rows$cex)
  graphics::axis(1, at = c(1, 1.25, 1.5), labels = c("1", "1.25", "1.5"), cex.axis = 1.05)
  graphics::mtext("mean HR per SD of the published scores (95% CI)", side = 1, line = 2.6, cex = 0.8, col = grey)
  graphics::mtext(rows$label, side = 2, at = y, las = 1, line = 0.4, cex = ifelse(seq_len(n) == 1, 0.92, 0.78),
                  font = ifelse(seq_len(n) == 1, 2, 1), col = ifelse(seq_len(n) == 1, navy, "grey25"))
  graphics::mtext(label_mid, side = 3, line = 1.6, font = 2, cex = 1.15, col = navy)
  graphics::mtext(sprintf("95%% CI %.2f-%.2f%s", hr_lo, hr_hi, if (all(str_$verdict == tr$verdict)) "; same in every setting" else ""),
                  side = 3, line = 0.4, cex = 0.78, col = grey)

  # Right: H1 against expression-matched random gene sets
  graphics::par(mar = c(4.6, 0.6, 4.0, 0.6))
  graphics::plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  graphics::rect(0.04, 0.12, 0.96, 0.88, border = navy, lwd = 1.6, col = "#F3F6FA")
  graphics::text(0.5, 0.68, gsub(" random", "\nrandom", label_right), font = 2, cex = 1.45, col = navy)
  graphics::text(0.5, 0.38, sprintf("%d of %d independent pairs\nafter Holm correction", n_beat, nrow(h1)), cex = 1.05, col = "grey20")
  graphics::text(0.5, 0.20, "public cohorts: too few lncRNAs\nmeasured, too few deaths", cex = 0.85, col = grey)
}

out <- "results/paper/figures/graphical_abstract"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
w <- 7; h <- 2.8; dpi <- 400
grDevices::pdf(paste0(out, ".pdf"), width = w, height = h); draw_ga(); grDevices::dev.off()
grDevices::png(paste0(out, ".png"), width = w, height = h, units = "in", res = dpi); draw_ga(); grDevices::dev.off()
grDevices::tiff(paste0(out, ".tif"), width = w, height = h, units = "in", res = dpi, compression = "lzw"); draw_ga(); grDevices::dev.off()
message("Graphical abstract: ", out, ".pdf/.png/.tif (", w * dpi, " x ", h * dpi, " pixels)")
message("Labels: ", label_left, " | ", label_mid, " | ", label_right)
