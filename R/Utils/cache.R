# Content-addressed calculation caches in results/cache.
# The signature includes the numerical core (R/Method), so changing it creates
# new caches, while reorganising studies, tables or figures keeps them valid.

fbst_hash <- function(x) {
  p <- tempfile(fileext = ".rds"); on.exit(unlink(p))
  saveRDS(x, p, version = 3, compress = FALSE); unname(tools::md5sum(p))
}

.fbst_code_snapshot <- local({
  files <- sort(list.files("R/Method", pattern = "\\.(R|cpp)$", full.names = TRUE))
  tools::md5sum(files)
})
fbst_code_hashes <- function() .fbst_code_snapshot

fbst_cache_compute <- function(kind, parameters, compute) {
  signature <- list(version = FBST_VERSION, code = fbst_code_hashes(), parameters = parameters)
  d <- file.path(fbst_output_dir(), "cache", kind)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  p <- file.path(d, paste0(fbst_hash(signature), ".rds"))
  if (file.exists(p)) {
    obj <- readRDS(p)
    if (identical(obj$signature, signature)) return(obj$value)
    stop("Cache signature mismatch: ", p)
  }
  t <- proc.time()[[3]]
  value <- if (is.function(compute)) compute() else force(compute)
  obj <- list(
    signature = signature, value = value, elapsed_seconds = proc.time()[[3]] - t,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  tmp <- paste0(p, ".tmp-", Sys.getpid()); saveRDS(obj, tmp)
  if (!file.rename(tmp, p)) stop("Could not atomically save cache: ", p)
  value
}
