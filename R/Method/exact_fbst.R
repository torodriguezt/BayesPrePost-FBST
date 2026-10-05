# Load the deterministic FBST core. Prior order: (lambda0, lambda1, lambda2).
# Numerical definitions are grouped by responsibility; the public API is unchanged.
FBST_CORE_VERSION <- "r-2026-09-24-1"
.fbst_env <- new.env(parent = baseenv())
.fbst_axis_cache <- new.env(parent = emptyenv())
.fbst_sources <- vapply(sys.frames(), function(frame) {
  if (is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)[1L]
}, character(1))
.fbst_sources <- .fbst_sources[!is.na(.fbst_sources) & basename(.fbst_sources) == "exact_fbst.R"]
.fbst_path <- normalizePath(
  if (length(.fbst_sources)) tail(.fbst_sources, 1L) else "R/Method/exact_fbst.R",
  mustWork = TRUE
)
.fbst_core_dir <- dirname(.fbst_path)
.fbst_cpp <- file.path(.fbst_core_dir, "exact_fbst_core.cpp")
.fbst_core_modules <- file.path(.fbst_core_dir, c(
  "posterior_kernel.R", "evidence.R", "predictive.R",
  "decision_rules.R", "posterior_summary.R"
))
.fbst_core_files <- c(
  .fbst_path, .fbst_cpp, .fbst_core_modules,
  file.path(.fbst_core_dir, "priors.R")
)
if (any(!file.exists(.fbst_core_files))) {
  stop(
    "Required numerical companions not found: ",
    paste(.fbst_core_files[!file.exists(.fbst_core_files)], collapse = ", ")
  )
}
# Checkpoints must account for every numerical module after the split.
.fbst_source_hash <- unname(tools::md5sum(.fbst_core_files))
if (requireNamespace("Rcpp", quietly = TRUE)) {
  Rcpp::sourceCpp(.fbst_cpp, env = .fbst_env, rebuild = FALSE, showOutput = FALSE)
}
for (.fbst_module in .fbst_core_modules) {
  source(.fbst_module, local = TRUE)
}
rm(.fbst_module)
