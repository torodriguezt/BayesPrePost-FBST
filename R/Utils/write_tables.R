# Export numerical results as CSV and plain LaTeX tables.

escape_latex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)

  for (s in c("&", "%", "$", "#", "_")) x <- gsub(s, paste0("\\", s), x, fixed = TRUE)
  x
}

write_result_table <- function(data, name) {
  path <- file.path(tables_directory(), name)
  write.csv(data, paste0(path, ".csv"), row.names = FALSE, na = "NA")
  formatted <- lapply(
    data,
    function(x) {
      if (is.numeric(x)) {
        ifelse(
          is.na(x),
          "NA",
          format(x, digits = 10, trim = TRUE, scientific = NA)
        )
      } else {
        as.character(x)
      }
    }
  )
  values <- do.call(cbind, lapply(formatted, escape_latex))
  rows <- apply(values, 1, function(x) {
    paste0(paste(x, collapse = " & "), " \\\\")
  })

  writeLines(
    c(
      paste0("% Article label: tab:", name),
      paste0(
        "\\begin{tabular}{",
        paste(rep("l", ncol(data)), collapse = ""),
        "}"
      ),
      "\\hline",
      paste0(paste(escape_latex(names(data)), collapse = " & "), " \\\\"),
      "\\hline",
      rows,
      "\\hline",
      "\\end{tabular}"
    ),
    paste0(path, ".tex")
  )

  invisible(paste0(path, ".csv"))
}
