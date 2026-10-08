# Compute type-I error under unequal sample sizes and partial overlap.

compute_overlap_error <- function() {
  sizes <- rbind(
    SCENARIOS$boundary_sizes,
    cbind(150, c(75, 105, 120, 135, 150))
  )
  theta <- SCENARIOS$baseline
  rows <- list()

  for (i in seq_len(nrow(sizes))) {
    n1 <- sizes[i, 1]
    n2 <- sizes[i, 2]
    design <- calibrate_test(n1, n2, PRIORS$KL)
    reject <- design$EV <= design$calibration$kstar

    for (overlap in c(0, floor(min(n1, n2) / 2), min(n1, n2))) {
      for (psi in SCENARIOS$overlap_psi) {
        rows[[length(rows) + 1L]] <- data.frame(
          n1 = n1,
          n2 = n2,
          overlap = overlap,
          psi = psi,
          alpha_star = design$calibration$alpha,
          rejection = sum(sampling_distribution(n1, n2, theta, theta, psi, overlap) * reject)
        )
      }
    }
  }

  write_result_table(do.call(rbind, rows), "unequal")
}
