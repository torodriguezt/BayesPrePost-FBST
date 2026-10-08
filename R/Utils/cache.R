# Reuse calculations with the same parameters and numerical code.

object_hash <- function(x) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(x, path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}

.numerical_code_snapshot <- local({
  files <- sort(list.files("R/Method", pattern = "\\.(R|cpp)$", full.names = TRUE))
  tools::md5sum(files)
})

numerical_code_hashes <- function() {
  .numerical_code_snapshot
}

cached_calculation <- function(kind, parameters, compute) {
  signature <- list(
    version = ANALYSIS_VERSION,
    code = numerical_code_hashes(),
    parameters = parameters
  )
  directory <- file.path(results_directory(), "cache", kind)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, paste0(object_hash(signature), ".rds"))

  if (file.exists(path)) {
    return(readRDS(path)$value)
  }

  started <- proc.time()[[3]]
  value <- compute()
  cached <- list(
    signature = signature,
    value = value,
    elapsed_seconds = proc.time()[[3]] - started,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(cached, temporary)

  if (!file.rename(temporary, path)) {
    stop("Could not save cache: ", path)
  }

  value
}
