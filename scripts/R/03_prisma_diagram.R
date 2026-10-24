# ============================================================
# 03_prisma_diagram.R - PRISMA flow diagram
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

prisma_data <- data.frame(
  stage = c(
    "Records identified from databases\n(n = ~500)",
    "Records after duplicates removed\n(n = ~400)",
    "Records screened by title/abstract\n(n = ~400)",
    "Full-text articles assessed\n(n = ~100)",
    "Studies included in qualitative synthesis\n(n = ~75)",
    "Signatures with complete formulas\n(n = 55)",
    "Top signatures selected for analysis\n(n = 15)"
  ),
  y = 7:1,
  type = c("id", "id", "sc", "sc", "in", "in", "in")
)

p <- ggplot(prisma_data, aes(x = 1, y = y)) +
  geom_label(aes(label = stage, fill = type),
             size = 3.5,
             label.padding = grid::unit(0.4, "lines"),
             color = "black",
             fontface = "bold") +
  scale_fill_manual(values = c(
    "id" = "#AED6F1",
    "sc" = "#A9DFBF",
    "in" = "#F9E79F"
  )) +
  geom_segment(aes(x = 1, xend = 1, y = y - 0.18, yend = y - 0.82),
               arrow = grid::arrow(length = grid::unit(0.25, "cm")),
               color = "gray40") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  labs(title = "PRISMA Flow Diagram\nLUAD lncRNA Prognostic Signatures")

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

ggsave("results/figures/PRISMA_flow_diagram.pdf", p, width = 9, height = 11)
ggsave("results/figures/PRISMA_flow_diagram.png", p, width = 9, height = 11, dpi = 300)

cat("========================================\n")
cat("PRISMA DIAGRAM SAVED\n")
cat("========================================\n")
cat("  - results/figures/PRISMA_flow_diagram.pdf\n")
cat("  - results/figures/PRISMA_flow_diagram.png\n")
