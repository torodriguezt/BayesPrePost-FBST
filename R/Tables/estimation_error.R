# Compute bias and RMSE of the posterior mean.

compute_estimation_error <- function() {
  scenarios <- list(
    KL = list(
      prior = PRIORS$KL,
      theta = .1
    ),
    centred = list(
      prior = PRIORS$estimation,
      theta = .5
    ),
    conflict = list(
      prior = PRIORS$estimation,
      theta = .1
    )
  )
  rows <- list()

  for (name in names(scenarios)) {
    scenario <- scenarios[[name]]

    for (n in SCENARIOS$appendix_n) {
      message("Estimation: ", name, ", n = ", n)
      values <- cached_calculation("estimation_posterior_mean", list(
        n = n,
        prior = scenario$prior,
        algorithm = "series-normaliser-ratio-v1"
      ), function() {
        counts <- expand.grid(x1 = 0:n, x2 = 0:n)
        counts$mean <- mapply(posterior_mean_theta1, counts$x1, n, counts$x2, n,
          MoreArgs = list(prior = scenario$prior)
        )
        counts
      })
      theta <- scenario$theta
      mass <- dbinom(values$x1, n, theta) * dbinom(values$x2, n, theta)
      rows[[length(rows) + 1L]] <- data.frame(
        scenario = name,
        n = n,
        bias = sum(mass * (values$mean - theta)),
        rmse = sqrt(sum(mass * (values$mean - theta)^2)),
        rmse_sample_proportion = sqrt(theta * (1 - theta) / n)
      )
    }
  }

  write_result_table(do.call(rbind, rows), "estimation")
}
