# Section 4 (clustering): daily use and addiction with the margins as observed,
# recalibrating k* at the effective sample sizes n/D, rounded (e.g. 400/271 at D = 1.5).
# The fixed-cutoff comparison is vaping_design_effect.csv, written by tables.R.
# Output: results/design_effect_recalibrated.csv.

fbst_run_design_effect <- function() {
  dat <- fbst_vaping_data()
  g <- dat[dat$group == "daily_use_addiction", , drop = FALSE]
  stopifnot(nrow(g) == 1L, g$full_available)
  prior <- FBST_PRIORS$KL
  rows <- list()
  for (D in FBST_SCENARIOS$design_effects) {
    # The sample space needs integer sizes: the effective sizes are rounded.
    n1 <- round(g$n1 / D)
    n2 <- round(g$n2 / D)
    message("Design effect ", D, ": calibration at ", n1, "/", n2)
    d <- fbst_get_design(n1, n2, prior)
    r <- fbst_vaping_fit(g, prior, "KL", D, "observed", calibrate = FALSE)
    cal <- d$calibration
    rows[[length(rows) + 1L]] <- data.frame(
      item = g$group, sample = "full", D = D,
      n1_effective = g$n1 / D, n2_effective = g$n2 / D,
      n1_calibration = n1, n2_calibration = n2,
      ev = r$ev, P_increase = r$prob_gt,
      delta_mean = r$delta_mean, delta_lo = r$delta_lo, delta_hi = r$delta_hi,
      k_star = cal$kstar, k_near_lo = cal$k_range[1], k_near_hi = cal$k_range[2],
      alpha_star = cal$alpha, beta_star = cal$beta, reject = r$ev <= cal$kstar
    )
  }
  out <- do.call(rbind, rows)
  out$k_star_D1 <- out$k_star[out$D == 1]
  out$reject_with_D1_cutoff <- out$ev <= out$k_star_D1
  fbst_write_table(out, "design_effect_recalibrated")
  invisible(out)
}
