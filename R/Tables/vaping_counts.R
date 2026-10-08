# Summarize linked transition counts and full-sample margins.

summarize_vaping_counts <- function() {
  data <- read_vaping_counts()
  columns <- c(
    "item",
    "sample",
    "n1",
    "n2",
    "x1",
    "x2",
    "n00",
    "n01",
    "n10",
    "n11",
    "psi"
  )

  write_result_table(data[columns], "vaping_data")
}
