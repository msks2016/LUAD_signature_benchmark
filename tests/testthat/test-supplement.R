test_that("supplementary table helpers keep values and change labels only", {
  expect_identical(sig_or_value(c("Paper62_Yao2021", "GSE31210", "")), c("Yao 2021", "GSE31210", ""))
  expect_identical(recode(c("a", "b", "c"), c(a = "A", c = "C")), c("A", "b", "C"))
  d <- data.frame(x = 1, y = 2)
  expect_identical(names(rename_cols(d, c(y = "Y"))), c("x", "Y"))
  expect_error(rename_cols(d, c(z = "Z")), "Columns not found: z")

  f <- tempfile(fileext = ".csv")
  writeLines(c('"num","fmt","p","flag","mixed","empty"',
               '0.56,"1.00","<0.001",TRUE,"0.47-0.64",""',
               '-0.0328,"0.50","0.32",FALSE,"1",""',
               '33850,"2.10","0.05",,"x",""'), f)
  r <- read_supp_csv(f)
  expect_identical(r$num, c(0.56, -0.0328, 33850))     # plain numbers become numbers
  expect_identical(r$fmt, c("1.00", "0.50", "2.10"))    # formatted numbers stay text
  expect_identical(r$p, c("<0.001", "0.32", "0.05"))
  expect_identical(r$flag, c("yes", "no", ""))
  expect_identical(r$mixed, c("0.47-0.64", "1", "x"))
  expect_identical(r$empty, c("", "", ""))
  expect_error(read_supp_csv(tempfile()), "run scripts 40 and 41 first")

  expect_identical(supp_col_widths(data.frame(`A long header` = "x", b = strrep("y", 100), check.names = FALSE)), c(8, 60))
})

test_that("the workbook writer makes a contents sheet and one sheet per table", {
  calls <- list()
  rec <- function(name, f) { force(f); function(...) { calls[[length(calls) + 1]] <<- list(name = name, args = list(...)); f(...) } }
  # Stubs with the argument names of openxlsx 4.2, so a wrong argument name fails here.
  stub <- list(
    createWorkbook = function(creator = "", title = NULL, subject = NULL, category = NULL) "wb",
    createStyle = function(fontName = NULL, fontSize = NULL, fontColour = NULL, numFmt = "GENERAL", border = NULL,
                           borderColour = NULL, borderStyle = NULL, bgFill = NULL, fgFill = NULL, halign = NULL,
                           valign = NULL, textDecoration = NULL, wrapText = FALSE, textRotation = NULL, indent = NULL,
                           locked = NULL, hidden = NULL) "style",
    addWorksheet = function(wb, sheetName, gridLines = TRUE, tabColour = NULL, zoom = 100, header = NULL, footer = NULL) NULL,
    writeData = function(wb, sheet, x, startCol = 1, startRow = 1, array = FALSE, xy = NULL, colNames = TRUE,
                         rowNames = FALSE, headerStyle = NULL, borders = "none", borderColour = NULL, borderStyle = NULL,
                         withFilter = FALSE, keepNA = FALSE, na.string = NULL, name = NULL, sep = ", ") NULL,
    addStyle = function(wb, sheet, style, rows, cols, gridExpand = FALSE, stack = FALSE) NULL,
    mergeCells = function(wb, sheet, cols, rows) NULL,
    setRowHeights = function(wb, sheet, rows, heights) NULL,
    setColWidths = function(wb, sheet, cols, widths = 8.43, hidden = rep(FALSE, length(cols)), ignoreMergedCells = FALSE) {
      stopifnot(length(widths) == 1 || length(widths) == length(cols)); NULL
    },
    freezePane = function(wb, sheet, firstActiveRow = NULL, firstActiveCol = NULL, firstRow = FALSE, firstCol = FALSE) NULL,
    saveWorkbook = function(wb, file, overwrite = FALSE, returnValue = FALSE) NULL)
  ox <- function(f) rec(f, stub[[f]])
  tabs <- list(list(id = "S1", title = "First", note = "A note.", data = data.frame(a = 1:2, b = c("x", "y"))),
               list(id = "S2", title = "Second", note = strrep("Long note. ", 40), data = data.frame(c = "z")))
  write_supp_workbook(tabs, "out.xlsx", heading = "Heading", ox = ox)
  sheets <- vapply(Filter(function(x) x$name == "addWorksheet", calls), function(x) x$args$sheetName, "")
  expect_identical(sheets, c("Contents", "S1", "S2"))
  saved <- Filter(function(x) x$name == "saveWorkbook", calls)
  expect_length(saved, 1)
  expect_identical(saved[[1]]$args$file, "out.xlsx")
  titles <- vapply(Filter(function(x) x$name == "writeData" && is.character(x$args$x) && x$args$startRow == 1, calls),
                   function(x) x$args$x, "")
  expect_true("Supplementary Table S2. Second" %in% titles)
})
