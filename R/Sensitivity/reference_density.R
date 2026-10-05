# SI tab:reference: flat reference against the prior reference r_ref = f, with the
# cutoff calibrated separately for each reference.
# Output: results/reference_sensitivity.csv.

# KL-optimal, informative, conflict and the five weak priors (N0 = 2).
fbst_reference_priors <- function() {
  priors <- FBST_PRIORS[c("KL", "informative", "conflict")]
  for (mu in FBST_SCENARIOS$weak_mu) priors[[paste0("weak_mu", mu)]] <- fbst_prior(mu, 2)
  priors
}

fbst_run_references <- function() {
  priors <- fbst_reference_priors()
  sizes <- FBST_SCENARIOS$boundary_sizes
  scenarios <- data.frame(
    theta1 = c(
      .001, .01, .05, .4, .95, .99, .999,
      .001, .05, .4, .6, .95, .999
    ),
    theta2 = c(
      .001, .01, .05, .4, .95, .99, .999,
      .01, .01, .5, .4, .99, .99
    )
  )
  rows <- list()
  for (i in seq_len(nrow(sizes))) {
    for (name in names(priors)) {
      n1 <- sizes[i, 1]
      n2 <- sizes[i, 2]
      for (reference in c("flat", "prior")) {
        design <- fbst_get_design(n1, n2, priors[[name]], reference = reference)
        reject <- design$EV <= design$calibration$kstar
        for (j in seq_len(nrow(scenarios))) {
          s <- scenarios[j, ]
          rows[[length(rows) + 1L]] <- data.frame(
            n1 = n1, n2 = n2, prior = name,
            reference = reference, null_prior = "marginal", overlap = 0, psi = 1,
            s, as.list(fbst_design_fields(design)),
            rate = sum(fbst_sampling_mass(n1, n2, s$theta1, s$theta2, overlap = 0) * reject)
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "reference_sensitivity")
  invisible(out)
}
