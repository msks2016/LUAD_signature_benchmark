# scripts/21_prisma_diagram.R
# PRISMA-style flow diagram for signature selection.

# --- Fill these with your real screening counts ---
initial_records               <- 404
after_dup_removal             <- 404     # adjust if duplicates existed
titles_abstracts_screened     <- 404
full_text_screened            <- 90      # ← change to your real number
included_signatures           <- 15
excluded_after_title_abstract <- titles_abstracts_screened - full_text_screened
excluded_after_full_text      <- full_text_screened - included_signatures

prisma_counts <- data.frame(
  step = c("Records identified",
           "After duplicates removed",
           "Title/abstract screening",
           "Full-text screening",
           "Included in benchmark"),
  n = c(initial_records, after_dup_removal,
        titles_abstracts_screened, full_text_screened, included_signatures)
)
write.csv(prisma_counts, "results/tables/prisma_counts.csv", row.names = FALSE)
cat("=== PRISMA counts ===\n"); print(prisma_counts)

# --- Diagram (base R) ---
draw_prisma <- function() {
  par(mar = c(1, 1, 1, 1))
  plot(0:100, 0:100, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 100), ylim = c(0, 100))
  
  box <- function(x, y, w, h, label, cex = 0.85, col = "white") {
    rect(x-w/2, y-h/2, x+w/2, y+h/2, col = col, border = "black", lwd = 1.2)
    text(x, y, label, cex = cex)
  }
  
  left_x  <- 33
  right_x <- 78
  
  box(left_x, 92, 48, 10,
      paste0("Records identified via PubMed\n(n = ", initial_records, ")"))
  box(left_x, 77, 48, 9,
      paste0("Records after duplicates removed\n(n = ", after_dup_removal, ")"))
  box(left_x, 62, 48, 9,
      paste0("Records screened on title and abstract\n(n = ", titles_abstracts_screened, ")"))
  box(right_x, 62, 30, 10,
      paste0("Excluded (n = ", excluded_after_title_abstract, "):\n",
             "reviews, no coefficients,\nnon-LUAD, non-human"),
      cex = 0.7, col = "#F5F5F5")
  box(left_x, 44, 48, 9,
      paste0("Full-text articles assessed\n(n = ", full_text_screened, ")"))
  box(right_x, 44, 30, 11,
      paste0("Excluded (n = ", excluded_after_full_text, "):\n",
             "no risk formula,\nno gene list,\nimmune-only models,\nmixed cancers"),
      cex = 0.7, col = "#F5F5F5")
  box(left_x, 25, 48, 10,
      paste0("Signatures included in benchmark\n(n = ", included_signatures, ")"),
      col = "#E6F1EA")
  
  arrows(left_x, 87,   left_x, 81.5, length = 0.1, lwd = 1.2)
  arrows(left_x, 72.5, left_x, 66.5, length = 0.1, lwd = 1.2)
  arrows(left_x, 57.5, left_x, 48.5, length = 0.1, lwd = 1.2)
  arrows(left_x, 39.5, left_x, 30,   length = 0.1, lwd = 1.2)
  arrows(left_x + 24, 62, right_x - 15, 62, length = 0.1, lwd = 1.2)
  arrows(left_x + 24, 44, right_x - 15, 44, length = 0.1, lwd = 1.2)
}

pdf("results/figures/Figure1_PRISMA.pdf", width = 7.5, height = 8)
draw_prisma()
dev.off()

png("results/figures/Figure1_PRISMA.png", width = 7.5, height = 8,
    units = "in", res = 300)
draw_prisma()
dev.off()

cat("\nSaved Figure1_PRISMA.pdf and .png\n")