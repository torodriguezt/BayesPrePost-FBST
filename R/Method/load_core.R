# Load the numerical routines for posterior integration and FBST e-values.
CORE_VERSION <- "r-2026-09-24-1"
.numerical_core <- new.env(parent = baseenv())
.quadrature_cache <- new.env(parent = emptyenv())

.cpp_file <- "R/Method/numerical_core.cpp"
.core_modules <- file.path("R/Method", c(
  "posterior_kernel.R",
  "fbst_evalue.R",
  "prior_predictive.R",
  "adaptive_cutoff.R",
  "posterior_summary.R"
))

# Include the numerical source files in enumeration checkpoints.
.source_hash <- unname(tools::md5sum(c(
  "R/Method/load_core.R", .cpp_file, .core_modules,
  "R/Method/study_settings.R"
)))

Rcpp::sourceCpp(.cpp_file, env = .numerical_core, showOutput = FALSE)

for (.module in .core_modules) {
  source(.module, local = TRUE)
}

rm(.module)
