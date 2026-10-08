# Compute operating characteristics, application results and figures.
# Run from the repository root: Rscript main.R [one or more steps below]
source("R/load.R")

steps <- list(
  design = compute_power,
  boundary = compute_typeI_error,
  reference = compare_reference_densities,
  coverage = compute_interval_coverage,
  sensitivity = compare_prior_sensitivity,
  independent = compare_prior_dependence,

  vaping_data = summarize_vaping_counts,
  vaping = analyze_vaping,
  vaping_freq = compare_frequentist_tests,

  arcmarg = compare_null_priors,
  unequal = compute_overlap_error,
  estimation = compute_estimation_error,
  reconstruction = summarize_reconstructions,
  recoding = compare_outcome_coding,

  diagrams = plot_model_diagrams,
  intervals = plot_posterior_intervals,
  errors = plot_error_curves,
  densities = plot_prior_posterior
)

selected <- commandArgs(trailingOnly = TRUE)

if (!length(selected)) {
  selected <- names(steps)
}

figure_steps <- c("diagrams", "intervals", "errors", "densities")
selected <- unlist(lapply(selected, function(step) {
  if (step == "figures") {
    figure_steps
  } else if (step == "tables") {
    setdiff(names(steps), figure_steps)
  } else {
    step
  }
}), use.names = FALSE)
unknown <- setdiff(selected, names(steps))

if (length(unknown)) {
  stop(
    "Unknown step: ", paste(unknown, collapse = ", "),
    ". Available: ", paste(names(steps), collapse = ", ")
  )
}

for (step in unique(selected)) {
  message("\nReproduce: ", step)
  steps[[step]]()
}

writeLines(
  capture.output(sessionInfo()),
  file.path(results_directory(), "sessionInfo.txt")
)
