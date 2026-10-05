# Table 3 (tab:coverage): credible-interval coverage and tests under within-subject
# association, with every subject paired and the KL-optimal prior.
# Outputs: results/delta_coverage.csv (averaged over Delta), delta_coverage_by_delta.csv
# and paired_vs_mcnemar.csv (procedure against McNemar's test).

# Equal-tailed 95% intervals of delta for every (x1, x2) of the sample space.
fbst_interval_grid <- function(n1, n2, prior, family = "olkin_liu") {
  rows <- lapply(0:n1, function(x1) {
    fbst_cache_compute("interval_row", list(
      n1 = n1, n2 = n2, x1 = x1,
      prior = prior, family = family, interval = "equal_tailed_95"
    ), function() {
      t(vapply(0:n2, function(x2) {
        s <- posterior_summary(x1, n1, x2, n2, prior, family = family)
        c(lo = s$lo, hi = s$hi, mean = s$mean)
      }, c(lo = 0, hi = 0, mean = 0)))
    })
  })
  list(
    lo = do.call(rbind, lapply(rows, function(r) r[, "lo"])),
    hi = do.call(rbind, lapply(rows, function(r) r[, "hi"])),
    mean = do.call(rbind, lapply(rows, function(r) r[, "mean"]))
  )
}

fbst_run_coverage <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  rows <- list()
  for (n in FBST_SCENARIOS$paired_n) {
    intervals <- fbst_interval_grid(n, n, FBST_PRIORS$KL)
    design <- fbst_get_design(n, n, FBST_PRIORS$KL)
    for (psi in fbst_association_grid()) {
      for (delta in FBST_SCENARIOS$coverage_deltas) {
        mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta, psi)
        wald <- fbst_wald_properties(fbst_discordants(n, theta1, theta1 + delta, psi), n, delta)
        rows[[length(rows) + 1L]] <- data.frame(
          n = n, theta1 = theta1, delta = delta,
          psi = psi, overlap = n, as.list(fbst_design_fields(design)),
          coverage = sum(mass * (intervals$lo <= delta & intervals$hi >= delta)),
          width = sum(mass * (intervals$hi - intervals$lo)),
          wald_coverage = wald["coverage"], wald_width = wald["width"],
          rejection = sum(mass * (design$EV <= design$calibration$kstar))
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  average <- aggregate(cbind(coverage, width, wald_coverage, wald_width) ~ n + psi, out, mean)
  fbst_write_table(out, "delta_coverage_by_delta")
  fbst_write_table(average, "delta_coverage")
  invisible(out)
}

fbst_run_paired <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  rows <- list()
  for (n in FBST_SCENARIOS$paired_n) {
    design <- fbst_get_design(n, n, FBST_PRIORS$KL)
    reject <- design$EV <= design$calibration$kstar
    for (psi in fbst_association_grid()) {
      nominal <- sum(fbst_sampling_mass(n, n, theta1, theta1, psi) * reject)
      for (delta in FBST_SCENARIOS$deltas) {
        p <- fbst_cell_probs(theta1, theta1 + delta, psi)
        disc <- fbst_discordants(n, theta1, theta1 + delta, psi)
        mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta, psi)
        rows[[length(rows) + 1L]] <- data.frame(
          n = n, theta1 = theta1,
          delta = delta, psi = psi, overlap = n, as.list(p),
          as.list(fbst_design_fields(design)), rate_fbst = sum(mass * reject),
          mcnemar_matched_nominal = nominal,
          mcnemar_matched = fbst_mcnemar_rate(disc, nominal),
          mcnemar_alpha_nominal = design$calibration$alpha,
          mcnemar_at_alpha = fbst_mcnemar_rate(disc, design$calibration$alpha),
          mcnemar_05_nominal = .05, mcnemar_at_05 = fbst_mcnemar_rate(disc, .05)
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "paired_vs_mcnemar")
  invisible(out)
}

fbst_run_association <- function() {
  fbst_run_coverage()
  fbst_run_paired()
}
