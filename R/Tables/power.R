# Compute rejection rates and matched-level power across sample sizes.

compute_power <- function() {
  theta <- SCENARIOS$baseline

  rows <- lapply(SCENARIOS$main_n, function(n) {
    message("Design: n = ", n)

    design <- calibrate_test(n, n, PRIORS$KL)
    calibration <- design$calibration
    reject <- design$EV <= calibration$kstar

    rates <- vapply(SCENARIOS$deltas, function(delta) {
      sum(sampling_distribution(n, n, theta, theta + delta) * reject)
    }, numeric(1))

    z <- z_test_pvalue_grid(n, n)
    alternative <- sampling_distribution(n, n, theta, theta + .10)

    data.frame(
      n = n,
      k_star = calibration$kstar,
      k_near_lo = calibration$k_range[1],
      k_near_hi = calibration$k_range[2],
      alpha_star = calibration$alpha,
      rejection_0 = rates[1],
      rejection_005 = rates[2],
      rejection_010 = rates[3],
      rejection_020 = rates[4],
      z_matched = sum(alternative[z <= rates[1]]),
      z_005 = sum(alternative[z <= .05])
    )
  })

  write_result_table(do.call(rbind, rows), "design")
}
