# Compare e-value distributions under restricted and marginal null priors.

weighted_ecdf <- function(evidence, mass) {
  order <- order(as.vector(evidence))
  values <- as.vector(evidence)[order]
  cumulative <- cumsum(as.vector(mass)[order])
  last <- !duplicated(values, fromLast = TRUE)

  data.frame(
    k = values[last],
    alpha = cumulative[last]
  )
}

compare_null_priors <- function() {
  rows <- list()

  for (name in c("informative", "conflict", "KL")) {
    prior <- PRIORS[[name]]

    # The arc-length restriction is improper for the KL-optimal prior.
    null <- if (name == "KL") {
      "lor"
    } else {
      "arc"
    }

    for (n in SCENARIOS$null_comparison_n) {
      message("Null prior: ", name, ", n = ", n)
      marginal <- calibrate_test(n, n, prior)
      restricted <- calibrate_test(n, n, prior, null = null)

      # Use the same evidence in both CDFs, taking the more precise refinement.
      evidence <- marginal$EV
      better <- attr(restricted$EV, "error_estimate") < attr(evidence, "error_estimate")
      evidence[better] <- restricted$EV[better]
      marginal_cdf <- weighted_ecdf(evidence, marginal$pH)
      restricted_cdf <- weighted_ecdf(evidence, restricted$pH)
      difference <- abs(restricted_cdf$alpha - marginal_cdf$alpha)

      for (k in c(.1, .3, .5, .7)) {
        rows[[length(rows) + 1L]] <- data.frame(
          prior = name,
          n = n,
          restriction = null,
          k = k,
          alpha_restricted = sum(restricted$pH[evidence <= k]),
          alpha_marginal = sum(marginal$pH[evidence <= k]),
          sup_below_1 = max(c(0, difference[marginal_cdf$k < 1])),
          sup_at_most_03 = max(c(0, difference[marginal_cdf$k <= .3]))
        )
      }
    }
  }

  write_result_table(do.call(rbind, rows), "arcmarg")
}
