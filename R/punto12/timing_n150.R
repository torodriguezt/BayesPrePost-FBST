# Task 5(e): wall time of the complete calibration for n1=n2=150 under the KL prior,
# from an empty cache: all e-values of Omega (grid enumeration + adaptive refinement
# near competitive cutoffs), m (predictive_A), m_H (predictive_H) and k*.
args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(FBST_PROFILE = "full", FBST_OUTPUT_ROOT = args[1])
t_source <- system.time(source("R/pipeline_helpers.R"))
n <- 150L
t0 <- proc.time()
d <- fbst_get_design(n, n, FBST_PRIORS$KL)
el <- proc.time() - t0
cache <- file.path(fbst_output_dir(), "cache")
part <- function(kind) {
  f <- list.files(file.path(cache, kind), full.names = TRUE)
  if (!length(f)) return(NA_real_)
  sum(vapply(f, function(p) readRDS(p)$elapsed_seconds, 0))
}
diag <- attr(d$EV, "diagnostics")
out <- data.frame(n1 = n, n2 = n, prior = "KL", points = length(d$EV),
  wall_seconds = unname(el["elapsed"]), cpu_user_seconds = unname(el["user.self"]),
  evidence_grid_seconds = part("evidence"), predictive_A_seconds = part("predictive_A"),
  predictive_H_seconds = part("predictive_H"), refinement_and_kstar_seconds = part("calibration_refined"),
  refined_evalues = diag$refinement_count, refinement_passes = diag$refinement_passes,
  kstar = d$calibration$kstar, alpha = d$calibration$alpha, beta = d$calibration$beta,
  sum_H = sum(d$pH), sum_A = sum(d$pA), rcpp_compile_and_source_seconds = unname(t_source["elapsed"]),
  started = format(Sys.time() - el["elapsed"], tz = "UTC", usetz = TRUE))
print(out, digits = 8)
saveRDS(out, file.path(args[1], "timing_n150.rds"))
