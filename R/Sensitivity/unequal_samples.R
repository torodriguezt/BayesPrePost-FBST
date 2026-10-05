# SI tab:unequal: null rejection and power with unequal stage sizes and partial
# overlap between the stages, against McNemar's test and the z-test at 0.05.
# Output: results/unequal_samples.csv.

fbst_run_unequal <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  sizes <- unique(rbind(
    FBST_SCENARIOS$boundary_sizes,
    cbind(150, 150 - round(150 * c(0, .1, .2, .3, .5)))
  ))
  rows <- list()
  for (i in seq_len(nrow(sizes))) {
    n1 <- sizes[i, 1]
    n2 <- sizes[i, 2]
    design <- fbst_get_design(n1, n2, FBST_PRIORS$KL)
    reject <- design$EV <= design$calibration$kstar
    for (overlap in unique(c(0L, floor(min(n1, n2) / 2), min(n1, n2)))) {
      for (psi in fbst_association_grid()) {
        for (delta in c(-.2, -.1, 0, .1, .2)) {
          mass <- fbst_sampling_mass(n1, n2, theta1, theta1 + delta, psi, overlap)
          disc <- fbst_discordants(overlap, theta1, theta1 + delta, psi)
          rows[[length(rows) + 1L]] <- data.frame(
            n1 = n1, n2 = n2,
            overlap = overlap, unpaired_pre = n1 - overlap, unpaired_post = n2 - overlap,
            theta1 = theta1, delta = delta, psi = psi,
            as.list(fbst_design_fields(design)), rate_fbst = sum(mass * reject),
            mcnemar_nominal = .05,
            rate_mcnemar = if (overlap > 0) fbst_mcnemar_rate(disc, .05) else NA_real_,
            z_finite_05 = sum(mass[fbst_z_pvalues(n1, n2) <= .05])
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "unequal_samples")
  invisible(out)
}
