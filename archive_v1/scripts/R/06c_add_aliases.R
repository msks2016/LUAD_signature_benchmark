message("Loading clean expression...")
expr <- readRDS("data/processed/tcga_luad_expr_clean.rds")
cat("Before:", nrow(expr), "genes\n")
aliases_to_add <- list(
  "HSPC324"       = c("AL590226.1"),
  "KIAA0125"      = c("AC079766.1"),
  "RP11-21L23.2"  = c("AC110792.3"),
  "AL034397.3"    = c("AL034397.1"),
  "AC010980.2"    = c("AC010980.1"),
  "AC026355.2"    = c("AC026355.1"),
  "FAM30A"        = c("KIAA0125"),
  "MIR325HG"      = c("AC079766.1")
)
added <- 0
for (current in names(aliases_to_add)) {
  if (current %in% rownames(expr)) {
    for (alias in aliases_to_add[[current]]) {
      if (!(alias %in% rownames(expr))) {
        expr <- rbind(expr, expr[current, , drop = FALSE])
        rownames(expr)[nrow(expr)] <- alias
        added <- added + 1
        message("Added: ", current, " -> ", alias)
      }
    }
  }
}
print("cleared")
