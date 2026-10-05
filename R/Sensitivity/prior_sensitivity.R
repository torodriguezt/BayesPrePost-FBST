# Table 4 (tab:sens) and SI tab:sens_summary, tab:S_sens: rejection rates and the
# difference with the z-test at a matched level, by prior category.
# Categories: W weak (N0 = 2), N strong near the baseline and F strong far from it
# (N0 = 50); five values of mu0 each.
# Outputs: results/sensitivity_individual.csv, sensitivity_average.csv, table_S1.csv.

fbst_run_prior_sensitivity <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  priors <- rbind(
    data.frame(category = "W", N0 = 2, mu0 = FBST_SCENARIOS$weak_mu),
    data.frame(category = "N", N0 = 50, mu0 = FBST_SCENARIOS$informative_mu),
    data.frame(category = "F", N0 = 50, mu0 = FBST_SCENARIOS$conflict_mu)
  )
  rows <- list()
  for (n in FBST_SCENARIOS$sensitivity_n) {
    zp <- fbst_z_pvalues(n, n)
    null <- fbst_sampling_mass(n, n, theta1, theta1)
    for (i in seq_len(nrow(priors))) {
      pr <- priors[i, ]
      prior <- fbst_prior(pr$mu0, pr$N0)
      message("Sensitivity: n = ", n, ", ", pr$category, ", mu0 = ", pr$mu0)
      d <- fbst_get_design(n, n, prior)
      reject <- d$EV <= d$calibration$kstar
      summ <- fbst_design_summary(d)
      level <- sum(null * reject)
      for (delta in FBST_SCENARIOS$deltas) {
        mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta)
        rate <- sum(mass * reject)
        z_exact <- sum(mass[zp <= level])
        z_normal <- fbst_z_power_normal(n, n, theta1, theta1 + delta, level)
        rows[[length(rows) + 1L]] <- data.frame(
          category = pr$category, mu0 = pr$mu0,
          N0 = pr$N0, n = n, delta = delta, summ, rate = rate, matched_level = level,
          z_exact_matched = z_exact, difference_exact = rate - z_exact,
          z_normal_matched = z_normal, difference_normal = rate - z_normal
        )
      }
    }
    fbst_write_table(do.call(rbind, rows), "sensitivity_individual")
  }
  out <- do.call(rbind, rows)
  fbst_write_table(aggregate(cbind(alpha_star, rate, difference_exact, difference_normal, excess) ~
    category + n + delta, out, mean), "sensitivity_average")
  at10 <- out[abs(out$delta - .10) < 1e-9, ]
  s1 <- aggregate(difference_exact ~ category + n, at10, mean)
  s1$max_abs_individual <- aggregate(
    difference_exact ~ category + n, at10,
    function(x) max(abs(x))
  )$difference_exact
  fbst_write_table(s1, "table_S1")
  invisible(out)
}
