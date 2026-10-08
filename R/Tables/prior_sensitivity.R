# Compute prior sensitivity of type-I error and power.

compare_prior_sensitivity <- function() {
  theta <- SCENARIOS$baseline
  priors <- rbind(
    data.frame(
      category = "W",
      N0 = 2,
      mu0 = SCENARIOS$weak_mu
    ),

    data.frame(
      category = "I",
      N0 = 50,
      mu0 = SCENARIOS$informative_mu
    ),

    data.frame(
      category = "C",
      N0 = 50,
      mu0 = SCENARIOS$conflict_mu
    )
  )
  rows <- list()

  for (n in SCENARIOS$sensitivity_n) {
    z <- z_test_pvalue_grid(n, n)
    null <- sampling_distribution(n, n, theta, theta)

    for (i in seq_len(nrow(priors))) {
      specification <- priors[i, ]
      prior <- beta_prior(specification$mu0, specification$N0)
      message(
        "Sensitivity: n = ", n, ", ", specification$category,
        ", mu0 = ", specification$mu0
      )
      design <- calibrate_test(n, n, prior)
      reject <- design$EV <= design$calibration$kstar
      level <- sum(null * reject)

      for (delta in SCENARIOS$deltas) {
        mass <- sampling_distribution(n, n, theta, theta + delta)
        rate <- sum(mass * reject)
        rows[[length(rows) + 1L]] <- data.frame(
          category = specification$category,
          n = n,
          delta = delta,
          alpha_star = design$calibration$alpha,
          rejection = rate,
          difference_exact = rate - sum(mass[z <= level]),
          difference_normal = rate - z_test_power(n, n, theta, theta + delta, level)
        )
      }
    }
  }

  results <- do.call(rbind, rows)
  at10 <- results[results$delta == .10, ]

  write_result_table(
    aggregate(rejection ~ n + category + delta, results, mean),
    "sens_summary"
  )

  write_result_table(aggregate(
    cbind(alpha_star, difference_normal) ~ n + category,
    at10, mean
  ), "sens_alpha")
  exact <- aggregate(difference_exact ~ n + category, at10, mean)
  maximum <- aggregate(difference_exact ~ n + category, at10, function(x) {
    max(abs(x))
  })
  names(maximum)[3] <- "max_abs_individual"

  write_result_table(merge(exact, maximum, by = c("n", "category")), "S_sens")

  invisible(results)
}
