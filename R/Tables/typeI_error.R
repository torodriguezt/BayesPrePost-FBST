# Compute type-I error across null proportions and prior specifications.

# Refine the theta grid around each local maximum.
typeI_error_grid <- function(reject, n1, n2) {
  theta <- sort(unique(c(
    SCENARIOS$boundary_theta,
    seq(0, 1, length.out = 201), 10^seq(-5, -1, length.out = 31),
    1 - 10^seq(-5, -1, length.out = 31)
  )))

  evaluate <- function(t) {
    sum(sampling_distribution(n1, n2, t, t) * reject)
  }

  rates <- vapply(theta, evaluate, numeric(1))
  peaks <- which(rates >= c(-Inf, head(rates, -1)) &
    rates >= c(tail(rates, -1), -Inf))
  refinement <- unlist(lapply(peaks, function(i) {
    seq(theta[max(1, i - 1)], theta[min(length(theta), i + 1)], length.out = 21)
  }))
  theta <- sort(unique(c(theta, refinement)))

  data.frame(
    theta = theta,
    rate = vapply(theta, evaluate, numeric(1))
  )
}

format_intervals <- function(theta, ok) {
  if (!any(ok)) {
    return("")
  }

  runs <- rle(ok)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  paste(sprintf(
    "[%.3f, %.3f]", theta[starts[runs$values]],
    theta[ends[runs$values]]
  ), collapse = " U ")
}

compute_typeI_error <- function() {
  sizes <- SCENARIOS$boundary_sizes
  rows <- list()

  for (name in c("KL", "informative", "conflict")) {
    prior <- PRIORS[[name]]

    for (i in seq_len(nrow(sizes))) {
      n1 <- sizes[i, 1]
      n2 <- sizes[i, 2]
      message("Boundary: ", name, ", ", n1, "/", n2)
      design <- calibrate_test(n1, n2, prior)
      reject <- design$EV <= design$calibration$kstar
      grid <- typeI_error_grid(reject, n1, n2)
      total <- outer(0:n1, 0:n2, "+")
      singular <- prior[2] + prior[3] + total - 2 < 0 |
        prior[1] + n1 + n2 - total - 2 < 0
      best <- which.max(grid$rate)
      rows[[length(rows) + 1L]] <- data.frame(
        prior = name,
        n1 = n1,
        n2 = n2,
        alpha_star = design$calibration$alpha,
        rejection_at_040 = sum(sampling_distribution(n1, n2, .4, .4) * reject),
        max_rejection = grid$rate[best],
        theta_at_max = grid$theta[best],
        theta_rejection_ge_0999 = format_intervals(grid$theta, grid$rate >= .999),
        singular_mass = sum(design$pH[singular])
      )
    }
  }

  write_result_table(do.call(rbind, rows), "boundary")
}
