# Compare type-I error and power under flat and prior references.

compare_reference_densities <- function() {
  rows <- list()

  for (name in c("KL", "informative", "conflict")) {
    sizes <- SCENARIOS$boundary_sizes

    if (name != "KL") {
      sizes <- sizes[-1, , drop = FALSE]
    }

    for (i in seq_len(nrow(sizes))) {
      n1 <- sizes[i, 1]
      n2 <- sizes[i, 2]

      for (reference in c("flat", "prior")) {
        message("Reference: ", name, ", ", n1, "/", n2, ", ", reference)
        design <- calibrate_test(n1, n2, PRIORS[[name]], reference = reference)
        calibration <- design$calibration
        reject <- design$EV <= calibration$kstar
        grid <- typeI_error_grid(reject, n1, n2)
        rows[[length(rows) + 1L]] <- data.frame(
          prior = name,
          n1 = n1,
          n2 = n2,
          reference = reference,
          k_star = calibration$kstar,
          alpha_star = calibration$alpha,
          beta_star = calibration$beta,
          null_at_040 = sum(sampling_distribution(n1, n2, .4, .4) * reject),
          null_max = max(grid$rate),
          power = sum(sampling_distribution(n1, n2, .4, .5) * reject)
        )
      }
    }
  }

  write_result_table(do.call(rbind, rows), "reference")
}
