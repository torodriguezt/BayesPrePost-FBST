# Compute interval coverage and power under within-subject association.

# Cache one row of equal-tailed 95% intervals at a time.
posterior_interval_grid <- function(n, prior) {
  rows <- lapply(0:n, function(x1) {
    cached_calculation("interval_row", list(
      n1 = n,
      n2 = n,
      x1 = x1,
      prior = prior,
      family = "olkin_liu",
      interval = "equal_tailed_95"
    ), function() {
      t(vapply(0:n, function(x2) {
        summary <- posterior_summary(x1, n, x2, n, prior)
        c(lo = summary$lo, hi = summary$hi)
      }, c(lo = 0, hi = 0)))
    })
  })
  list(
    lo = do.call(rbind, lapply(rows, function(x) {
      x[, "lo"]
    })),
    hi = do.call(rbind, lapply(rows, function(x) {
      x[, "hi"]
    }))
  )
}

compute_interval_coverage <- function() {
  theta <- SCENARIOS$baseline
  rows <- list()

  for (n in SCENARIOS$paired_n) {
    message("Coverage: n = ", n)
    intervals <- posterior_interval_grid(n, PRIORS$KL)
    design <- calibrate_test(n, n, PRIORS$KL)
    reject <- design$EV <= design$calibration$kstar

    for (psi in SCENARIOS$association_psi) {
      coverage <- vapply(
        SCENARIOS$coverage_deltas,
        function(delta) {
          mass <- sampling_distribution(n, n, theta, theta + delta, psi)
          wald <- wald_interval_properties(
            discordant_distribution(n, theta, theta + delta, psi),
            n,
            delta
          )
          c(
            proposed = sum(mass * (intervals$lo <= delta & intervals$hi >= delta)),
            wald_coverage = unname(wald["coverage"]),
            wald_width = unname(wald["width"])
          )
        },
        c(
          proposed = 0,
          wald_coverage = 0,
          wald_width = 0
        )
      )
      level <- sum(sampling_distribution(n, n, theta, theta, psi) * reject)
      discordants <- discordant_distribution(n, theta, theta + .10, psi)
      rows[[length(rows) + 1L]] <- data.frame(
        n = n,
        psi = psi,
        coverage = mean(coverage["proposed", ]),
        wald_coverage = mean(coverage["wald_coverage", ]),
        wald_width = mean(coverage["wald_width", ]),
        type_I = level,
        power = sum(sampling_distribution(n, n, theta, theta + .10, psi) * reject),
        mcnemar_matched = mcnemar_rejection_rate(discordants, level),
        mcnemar_at_alpha = mcnemar_rejection_rate(discordants, design$calibration$alpha)
      )
    }
  }

  write_result_table(do.call(rbind, rows), "coverage")
}
