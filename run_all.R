# Reproduces every table and figure of the paper. Run from the repository
# root with: Rscript run_all.R
# Outputs go to Figures/ and output/. Each step can also be run on its own.
# Total runtime: roughly 2 hours.

steps <- c(
  "R/02_fit_priors_kl.R",
  "R/03_calibrate_weights.R",
  "R/04_simulation_estimation_figures.R",
  "R/05_simulation_prior_sensitivity.R",
  "R/06_application_tvsfp.R",
  "R/07_application_figures.R"
)

t0 <- Sys.time()
for (s in steps) {
  cat(sprintf("\n========== %s ==========\n", s))
  t1 <- Sys.time()
  source(s, echo = FALSE)
  cat(sprintf("[%s finished in %.1f min | %.1f min total]\n",
              s,
              as.numeric(difftime(Sys.time(), t1, units = "mins")),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat(sprintf("\nAll steps completed in %.1f min.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
