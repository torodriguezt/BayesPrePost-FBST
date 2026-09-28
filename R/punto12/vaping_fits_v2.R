# Posterior summaries of delta for the 18 vaping analyses (3 items x 2 samples x 3
# priors) at D = 1, and the calibration of the six KL designs used by the error-curve
# figures. Same functions as run_application_vaping() / run_vaping_figures(); only
# D = 1 is computed and only the KL designs are calibrated. Writes to a scratch root.
# A KL calibration already cached by the full run of 2026-09-27 is reused only if the
# md5 of every file of the numerical path is identical to the current one; otherwise
# the design is enumerated and calibrated from scratch with fbst_get_design().
args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(FBST_PROFILE = "full", FBST_OUTPUT_ROOT = args[1], FBST_VAPING_DATA = args[2])
source("R/pipeline_helpers.R"); source("R/manuscript_applications.R"); source("R/application_vaping.R")
numeric_files <- c("R/exact_fbst.R", "R/exact_fbst_core.cpp", "R/pipeline_helpers.R", "R/manuscript_applications.R")
current <- tools::md5sum(numeric_files)
cached_design <- function(n1, n2, prior) {
  d <- file.path(fbst_output_dir(), "cache", "calibration_refined")
  for (f in list.files(d, full.names = TRUE)) {
    o <- readRDS(f); e <- o$signature$parameters$evidence
    if (e$n1 == n1 && e$n2 == n2 && e$reference == "flat" && e$family == "olkin_liu" &&
        isTRUE(all.equal(unname(e$prior), unname(prior), tolerance = 0)) &&
        identical(unname(o$signature$code[numeric_files]), unname(current)))
      return(list(value = o$value, source = paste("cache", basename(f), o$created_utc)))
  }
  NULL
}
dat <- fbst_vaping_data(); fbst_write_table(dat, "vaping_transition")
priors <- fbst_vaping_priors(); rows <- list(); designs <- list()
for (p in names(priors)) for (i in seq_len(nrow(dat))) for (analysis in c("complete", "observed")) {
  g <- dat[i, ]; t0 <- proc.time()[[3]]
  r <- fbst_vaping_fit(g, priors[[p]], p, 1, analysis, calibrate = FALSE)
  if (p == "KL") {
    args_d <- fbst_app_args(g, priors[[p]], 1, analysis)
    hit <- cached_design(args_d$n1, args_d$n2, priors[[p]])
    t1 <- proc.time()[[3]]
    design <- if (is.null(hit)) fbst_get_design(args_d$n1, args_d$n2, priors[[p]]) else
      local({
        pH <- predictive_H(args_d$n1, args_d$n2, priors[[p]]); pA <- predictive_A(args_d$n1, args_d$n2, priors[[p]])
        # The cached rule must be reproduced by the current predictives.
        chk <- adaptive_cutoff(hit$value$EV, pH, pA, tol = FBST_NUMERICS$near_risk)
        stopifnot(chk$kstar == hit$value$calibration$kstar, abs(chk$risk - hit$value$calibration$risk) < 1e-12)
        list(EV = hit$value$EV, calibration = hit$value$calibration, pH = pH, pA = pA)
      })
    cal <- design$calibration
    designs[[paste(g$group, r$sample)]] <- list(n1 = args_d$n1, n2 = args_d$n2, EV = design$EV,
      pH = design$pH, pA = if (!is.null(design$pA)) design$pA else NULL, calibration = cal,
      source = if (is.null(hit)) sprintf("computed now (%.0f s)", proc.time()[[3]] - t1) else hit$source)
    r$k_star <- cal$kstar; r$alpha_star_original <- cal$alpha; r$beta_star_original <- cal$beta
    r$risk_original <- cal$risk; r$reject <- r$ev <= cal$kstar; r$calibration_status <- "enumerated"
    saveRDS(designs, file.path(fbst_output_dir(), "vaping_KL_designs.rds"))
  }
  rows[[length(rows) + 1L]] <- r
  message(sprintf("%s / %s / %s: %.1f s", g$group, r$sample, p, proc.time()[[3]] - t0))
  fbst_write_table(do.call(rbind, rows), "vaping_results")
}
