# Reproduce the tables and figures of the article. From the repository root:
#   Rscript main.R                 every block, in the order of the article
#   Rscript main.R real_data       one or more blocks: simulation, sensitivity, real_data
# Tables are written to results/ (CSV and TeX) and figures to figures/.
# The full run takes many hours. Calibrated designs are cached in results/cache,
# so an interrupted run resumes where it stopped. See R/README.md.

source("R/load.R")

steps <- list(
  simulation = list(
    power = fbst_run_power, # Table 2 (tab:design)
    association = fbst_run_association, # Table 3 (tab:coverage)
    estimation = fbst_run_estimation # Appendix: tab:estimation, fig:est_prior_post
  ),
  sensitivity = list(
    prior_sensitivity = fbst_run_prior_sensitivity, # Table 4 (tab:sens); SI tab:sens_summary, tab:S_sens
    boundaries = fbst_run_boundaries, # Table 5 (tab:boundary)
    independent_priors = fbst_run_independent_priors, # Table 6 (tab:indep_decisions); SI tab:indep
    reference_density = fbst_run_references, # SI tab:reference
    unequal_samples = fbst_run_unequal, # SI tab:unequal
    null_prior = fbst_run_null_prior # Appendix / SI tab:arcmarg
  ),
  real_data = list(
    fits = fbst_run_vaping, # Table 7 (tab:vaping_data); inputs of the next steps
    tables = fbst_vaping_tables, # Table 8 (tab:vaping_ni); SI tab:S_vaping_inf, tab:S_vaping_conf, tab:vaping_freq
    figures = fbst_vaping_figures, # fig:vaping_errors, fig:vaping_forest
    design_effect = fbst_run_design_effect, # Section 4, clustering
    missing_data = fbst_run_missing_data, # Section 4, unlinked surveys (eq. mar-vape)
    recoding = fbst_run_recoding, # SI tab:S_recoding
    reconstruction = fbst_run_reconstruction # SI tab:S_reconstruction
  )
)

blocks <- commandArgs(trailingOnly = TRUE)
if (!length(blocks)) blocks <- names(steps)
unknown <- setdiff(blocks, names(steps))
if (length(unknown)) {
  stop("Unknown block: ", paste(unknown, collapse = ", "), ". Use: ", paste(names(steps), collapse = ", "))
}

for (block in blocks) {
  for (step in names(steps[[block]])) {
    message("\n== ", block, " / ", step, " ==")
    t0 <- proc.time()[[3L]]
    steps[[block]][[step]]()
    message(sprintf("%s / %s: %.1f min", block, step, (proc.time()[[3L]] - t0) / 60))
  }
}
writeLines(capture.output(sessionInfo()), file.path(fbst_output_dir(), "sessionInfo.txt"))
