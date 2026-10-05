# Table 2 (tab:design): operating characteristics under the KL-optimal prior and
# stage independence, against the exact z-test at the procedure's level and at 0.05.
# Also the asymptotic null rejection under association, eq. (eq:assocnull), Section 3.
# Outputs: results/power_exact.csv, results/association_check.csv.

fbst_run_power <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  rows <- assoc <- list()
  for (n in FBST_SCENARIOS$main_n) {
    message("Power: calibration for n = ", n)
    d <- fbst_get_design(n, n, FBST_PRIORS$KL)
    reject <- d$EV <= d$calibration$kstar
    bayes <- fbst_bayes_rule(d$pH, d$pA)
    summ <- fbst_design_summary(d)
    null <- fbst_sampling_mass(n, n, theta1, theta1)
    level <- sum(null * reject)
    zp <- fbst_z_pvalues(n, n)
    for (delta in FBST_SCENARIOS$deltas) {
      mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta)
      rate <- sum(mass * reject)
      z_matched <- sum(mass[zp <= level])
      rows[[length(rows) + 1L]] <- data.frame(
        n = n, theta1 = theta1, delta = delta, summ,
        rate = rate, matched_level = level,
        z_exact_matched = z_matched, z_exact_05 = sum(mass[zp <= .05]),
        z_exact_null_at_matched = sum(null[zp <= level]),
        z_normal_matched = fbst_z_power_normal(n, n, theta1, theta1 + delta, level),
        z_normal_05 = fbst_z_power_normal(n, n, theta1, theta1 + delta, .05),
        difference_exact = rate - z_matched,
        rate_bayes_rule = sum(mass * bayes$reject)
      )
    }
    # Asymptotic null rejection under association against the exact value with
    # every subject paired.
    if (n %in% FBST_SCENARIOS$paired_n) {
      for (psi in c(.2, 1, 3, 5)) {
        p <- fbst_cell_probs(theta1, theta1, psi)
        s2w <- 2 * theta1 * (1 - theta1)
        s2 <- s2w - 2 * (p[["p11"]] - theta1^2)
        assoc[[length(assoc) + 1L]] <- data.frame(
          n = n, psi = psi,
          k_star = d$calibration$kstar, variance_ratio = s2w / s2,
          approximation = stats::pchisq(-2 * log(d$calibration$kstar) * s2w / s2,
            df = 1,
            lower.tail = FALSE
          ),
          exact = sum(fbst_sampling_mass(n, n, theta1, theta1, psi) * reject)
        )
      }
    }
    fbst_write_table(do.call(rbind, rows), "power_exact")
    if (length(assoc)) fbst_write_table(do.call(rbind, assoc), "association_check")
  }
  invisible(do.call(rbind, rows))
}
