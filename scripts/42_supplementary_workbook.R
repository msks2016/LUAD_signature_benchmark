# 42_supplementary_workbook.R - Supplementary Tables S1-S19 in one Excel workbook: a contents
# sheet, then one sheet per table with its title and note. Values are copied from the tables in
# results/paper/tables (scripts 40 and 41); only labels change (column names, signature names,
# codes). Needs the openxlsx package: install.packages("openxlsx"), then renv::snapshot().
# Output: results/paper/Supplementary_Tables_S1-S19.xlsx

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
if (!requireNamespace("openxlsx", quietly = TRUE)) stop('Install openxlsx first: install.packages("openxlsx")')

tabs <- supp_workbook_tables("results/paper/tables")
out <- "results/paper/Supplementary_Tables_S1-S19.xlsx"
write_supp_workbook(tabs, out, heading = paste0(
  "Supplementary Tables S1–S19. Benchmark of published lncRNA prognostic signatures for lung adenocarcinoma"))
print(data.frame(Table = vapply(tabs, function(t) t$id, ""), Rows = vapply(tabs, function(t) nrow(t$data), 0L),
                 Columns = vapply(tabs, function(t) ncol(t$data), 0L), Source = vapply(tabs, function(t) t$file, "")),
      row.names = FALSE)
message("Written ", out)
