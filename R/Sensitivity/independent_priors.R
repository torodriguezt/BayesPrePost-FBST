# Table 6 (tab:indep_decisions) and SI tab:indep: separately calibrated procedures
# under Olkin-Liu and independent beta priors with matching marginals.
# Outputs: results/independent_prior_by_n.csv, independent_evalue_differences.csv,
# prior_probability_near_null.csv.

fbst_run_independent_priors <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  rows <- list()
  prior_rows <- list()
  for (name in c("KL", "informative")) {
    prior <- FBST_PRIORS[[name]]
    for (family in c("olkin_liu", "independent")) {
      p_near <- delta_cdf(.05, 0, 0, 0, 0, prior, family = family) -
        delta_cdf(-.05, 0, 0, 0, 0, prior, family = family)
      prior_rows[[length(prior_rows) + 1L]] <- data.frame(
        prior = name,
        family = family, delta_band = .05, probability = p_near
      )
    }
    for (n in sort(unique(c(FBST_SCENARIOS$independent_n, 50L, 200L)))) {
      biv <- fbst_get_design(n, n, prior, family = "olkin_liu")
      ind <- fbst_get_design(n, n, prior, family = "independent")
      rb <- biv$EV <= biv$calibration$kstar
      ri <- ind$EV <= ind$calibration$kstar
      for (delta in FBST_SCENARIOS$deltas) {
        for (psi in c(.2, 1, 5)) {
          mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta, psi)
          rows[[length(rows) + 1L]] <- data.frame(
            prior = name, n = n,
            theta1 = theta1, delta = delta, psi = psi, overlap = n,
            k_biv = biv$calibration$kstar, alpha_biv = biv$calibration$alpha,
            beta_biv = biv$calibration$beta, k_ind = ind$calibration$kstar,
            alpha_ind = ind$calibration$alpha, beta_ind = ind$calibration$beta,
            num_status_biv = fbst_design_fields(biv)$num_status,
            num_status_ind = fbst_design_fields(ind)$num_status,
            calibration_status_biv = fbst_design_fields(biv)$calibration_status,
            calibration_status_ind = fbst_design_fields(ind)$calibration_status,
            evidence_max_error_biv = fbst_design_fields(biv)$evidence_max_error,
            evidence_max_error_ind = fbst_design_fields(ind)$evidence_max_error,
            evidence_unresolved_biv = fbst_design_fields(biv)$evidence_unresolved,
            evidence_unresolved_ind = fbst_design_fields(ind)$evidence_unresolved,
            rejection_biv = sum(mass * rb), rejection_ind = sum(mass * ri),
            disagreement = sum(mass * (rb != ri)), mean_abs_ev_diff = sum(mass * abs(biv$EV - ind$EV))
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "independent_prior_by_n")
  fbst_write_table(do.call(rbind, prior_rows), "prior_probability_near_null")
  fbst_write_table(subset(out, delta == .2 & psi == 1), "independent_evalue_differences")
  invisible(out)
}
