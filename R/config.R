# Paths relative to the repository root.

PATHS <- list(
  vaping_data = "R/Data/vaping_mccauley2023.csv",
  results = "results",
  figures = "figures"
)

results_directory <- function() {
  dir.create(PATHS$results, recursive = TRUE, showWarnings = FALSE)
  PATHS$results
}

tables_directory <- function() {
  path <- file.path(results_directory(), "tables")
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

vaping_data_path <- function() {
  PATHS$vaping_data
}

figures_directory <- function() {
  dir.create(PATHS$figures, recursive = TRUE, showWarnings = FALSE)
  PATHS$figures
}
