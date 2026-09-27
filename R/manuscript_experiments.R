# Deterministic experiments for the manuscript.
# Source this file to define functions; the numbered scripts select experiments.
# Full profiles preserve manuscript sample sizes. Pilot n = 6/8 exercises the
# same algorithms and must never be reported as a manuscript reproduction.
if (!exists("fbst_get_design", mode = "function")) source("R/pipeline_helpers.R")

fbst_experiment_sizes <- function(full) {
  if (fbst_profile() == "pilot") c(6L, 8L) else full
}

fbst_mu_prior <- function(mu0, N0) {
  fbst_prior(mu0, N0)
}

fbst_association_grid <- function() sort(unique(c(FBST_SCENARIOS$coverage_psi,
                              FBST_SCENARIOS$test_psi, FBST_SCENARIOS$extra_psi)))

# Cells are ordered 00,01,10,11. Endpoint margins are allowed; association is
# not identified there, but the unique compatible cell distribution is valid.
fbst_cell_probs <- function(theta1, theta2, psi = 1) {
  stopifnot(length(theta1) == 1L, length(theta2) == 1L,
            is.finite(theta1), is.finite(theta2), is.finite(psi),
            theta1 >= 0, theta1 <= 1, theta2 >= 0, theta2 <= 1, psi > 0)
  lo <- max(0, theta1 + theta2 - 1)
  hi <- min(theta1, theta2)
  p11 <- if (hi == lo) lo else if (psi == 1) theta1 * theta2 else {
    uniroot(function(z) z * (1 - theta1 - theta2 + z) -
              psi * (theta1 - z) * (theta2 - z), c(lo, hi),
            tol = 1e-14)$root
  }
  p <- c(p00 = 1 - theta1 - theta2 + p11, p01 = theta2 - p11,
         p10 = theta1 - p11, p11 = p11)
  stopifnot(min(p) >= -1e-12, abs(sum(p) - 1) < 1e-12,
            abs(p[3] + p[4] - theta1) < 1e-12,
            abs(p[2] + p[4] - theta2) < 1e-12)
  # Only roundoff at exactly degenerate margins can produce a negative zero.
  p[p < 0] <- 0
  p
}

# Exact distribution of the margins: the overlap subjects contribute iid
# bivariate Bernoulli pairs; the remaining n1-overlap and n2-overlap subjects
# contribute independent Bernoulli observations at their respective stages.
# Matrix rows index X1 and columns index X2, both beginning at zero.
fbst_sampling_mass <- function(n1, n2, theta1, theta2, psi = 1,
                               overlap = min(n1, n2)) {
  stopifnot(n1 >= 0, n2 >= 0, n1 == as.integer(n1), n2 == as.integer(n2),
            overlap >= 0, overlap <= min(n1, n2), overlap == as.integer(overlap))
  p <- fbst_cell_probs(theta1, theta2, psi)
  if (psi == 1 || overlap == 0) {
    return(outer(dbinom(0:n1, n1, theta1), dbinom(0:n2, n2, theta2)))
  }
  mass <- matrix(1, 1L, 1L)
  if (overlap > 0) for (j in seq_len(overlap)) {
    next_mass <- matrix(0, j + 1L, j + 1L)
    i <- seq_len(j)
    next_mass[i, i] <- next_mass[i, i] + p[1] * mass
    next_mass[i, i + 1L] <- next_mass[i, i + 1L] + p[2] * mass
    next_mass[i + 1L, i] <- next_mass[i + 1L, i] + p[3] * mass
    next_mass[i + 1L, i + 1L] <- next_mass[i + 1L, i + 1L] + p[4] * mass
    mass <- next_mass
  }
  if (n1 > overlap) {
    new_mass <- matrix(0, n1 + 1L, ncol(mass))
    for (i in 0:(n1 - overlap)) {
      idx <- seq_len(nrow(mass)) + i
      new_mass[idx, ] <- new_mass[idx, ] + dbinom(i, n1 - overlap, theta1) * mass
    }
    mass <- new_mass
  }
  if (n2 > overlap) {
    new_mass <- matrix(0, nrow(mass), n2 + 1L)
    for (i in 0:(n2 - overlap)) {
      idx <- seq_len(ncol(mass)) + i
      new_mass[, idx] <- new_mass[, idx] + dbinom(i, n2 - overlap, theta2) * mass
    }
    mass <- new_mass
  }
  if (abs(sum(mass) - 1) > 5e-11) stop("Sampling distribution lost probability mass")
  mass
}

# Collapsing the concordant cells gives a trinomial distribution. Conditional
# on D discordant pairs, N01 is binomial(D, p01/(p01+p10)). This enumerates all
# discordant configurations without enumerating unnecessary concordant splits.
fbst_discordants <- function(n, theta1, theta2, psi) {
  stopifnot(n >= 0, n == as.integer(n))
  p <- fbst_cell_probs(theta1, theta2, psi)
  pd <- unname(p[2] + p[3])
  out <- lapply(0:n, function(d) {
    b <- 0:d
    w <- dbinom(d, n, pd) * dbinom(b, d, if (pd > 0) p[2] / pd else 0)
    statistic <- if (d == 0) 0 else (2 * b - d)^2 / d
    data.frame(n01 = b, n10 = d - b, mass = w, Q = statistic,
               pvalue = if (d == 0) 1 else pchisq(statistic, 1, lower.tail = FALSE))
  })
  out <- do.call(rbind, out)
  stopifnot(abs(sum(out$mass) - 1) < 5e-11)
  out
}

fbst_mcnemar_rate <- function(discordants, level) {
  stopifnot(is.finite(level), level >= 0, level <= 1)
  sum(discordants$mass[discordants$n01 + discordants$n10 > 0 &
                         discordants$pvalue <= level])
}

fbst_wald_properties <- function(discordants, n, delta) {
  if (n == 0) return(c(coverage = NA_real_, width = NA_real_))
  estimate <- (discordants$n01 - discordants$n10) / n
  half <- qnorm(.975) / n * sqrt(pmax(0, discordants$n01 +
                          discordants$n10 - n * estimate^2))
  c(coverage = sum(discordants$mass * (estimate - half <= delta &
                                      estimate + half >= delta)),
    width = sum(discordants$mass * 2 * half))
}

# Two-sided pooled z test. At all failures/all successes, the variance and
# observed difference are both zero; define p=1 and do not reject.
fbst_z_pvalues <- function(n1, n2) {
  x1 <- row(matrix(0, n1 + 1L, n2 + 1L)) - 1L
  x2 <- col(x1) - 1L
  pooled <- (x1 + x2) / (n1 + n2)
  se <- sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))
  z <- matrix(0, n1 + 1L, n2 + 1L)
  ok <- se > 0
  z[ok] <- abs(x2[ok] / n2 - x1[ok] / n1) / se[ok]
  2 * pnorm(-z)
}

# Retain the manuscript's normal approximation (unpooled alternative SE).
# A common nominal level does not imply equal finite-sample type-I errors.
fbst_z_power_normal <- function(n1, n2, theta1, theta2, level) {
  stopifnot(level >= 0, level <= 1)
  if (level == 0) return(0)
  if (level == 1) return(1)
  se <- sqrt(theta1 * (1 - theta1) / n1 + theta2 * (1 - theta2) / n2)
  if (se == 0) return(as.numeric(theta1 != theta2))
  d <- (theta2 - theta1) / se
  z <- qnorm(1 - level / 2)
  pnorm(d - z) + pnorm(-d - z)
}

fbst_design_fields <- function(design) {
  diagnostics <- attr(design$EV, "diagnostics")
  errors <- attr(design$EV, "error_estimate")
  unresolved <- attr(design$EV, "unresolved")
  count <- if (is.null(unresolved)) NA_integer_ else sum(unresolved)
  max_error <- if (is.null(errors)) NA_real_ else max(errors)
  status <- if (is.na(count)) "numerical_diagnostics_unavailable" else
    if (count > 0) "estimated_grid_requires_review" else "quadrature_tolerance_met"
  calibration_status <- design$calibration$numeric_status
  if (is.null(calibration_status)) calibration_status <- "grid_calibration_only"
  list(k_star = design$calibration$kstar, alpha_star = design$calibration$alpha,
    beta_star = design$calibration$beta, risk = design$calibration$risk,
    num_status = status, calibration_status = calibration_status,
    evidence_max_error = max_error, evidence_unresolved = count)
}

fbst_experiment_plot <- function(data, x, y, group, name, xlab = x, ylab = y) {
  if (!nrow(data)) return(invisible(NULL))
  dir <- fbst_figure_dir()
  draw <- function() {
    groups <- unique(as.character(data[[group]]))
    colors <- grDevices::hcl.colors(max(3, length(groups)), "Dark 3")
    plot(range(data[[x]], finite = TRUE), range(data[[y]], finite = TRUE),
         type = "n", xlab = xlab, ylab = ylab)
    for (i in seq_along(groups)) {
      d <- data[as.character(data[[group]]) == groups[i], ]
      d <- d[order(d[[x]]), ]
      lines(d[[x]], d[[y]], col = colors[i], type = "b", pch = i)
    }
    legend("bottomright", groups, col = colors[seq_along(groups)],
           pch = seq_along(groups), lty = 1, cex = .8, bty = "n")
  }
  grDevices::pdf(file.path(dir, paste0(name, ".pdf")), width = 7, height = 4.5)
  tryCatch(draw(), finally = grDevices::dev.off())
  grDevices::png(file.path(dir, paste0(name, ".png")), width = 1400, height = 900, res = 200)
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(NULL)
}

fbst_run_power <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  rows <- list()
  for (n in fbst_experiment_sizes(FBST_SCENARIOS$main_n)) {
    design <- fbst_get_design(n, n, FBST_PRIORS$KL)
    reject <- design$EV <= design$calibration$kstar
    null <- fbst_sampling_mass(n, n, theta1, theta1)
    level <- sum(null * reject)
    zp <- fbst_z_pvalues(n, n)
    for (delta in FBST_SCENARIOS$deltas) {
      mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta)
      rows[[length(rows) + 1L]] <- data.frame(n = n, theta1 = theta1,
        theta2 = theta1 + delta, delta = delta, psi = 1, overlap = n,
        as.list(fbst_design_fields(design)), rate = sum(mass * reject),
        matched_nominal = level,
        z_normal_matched = fbst_z_power_normal(n, n, theta1, theta1 + delta, level),
        z_normal_05 = fbst_z_power_normal(n, n, theta1, theta1 + delta, .05),
        z_finite_matched = sum(mass[zp <= level]), z_finite_05 = sum(mass[zp <= .05]))
    }
  }
  out <- do.call(rbind, rows)
  out$difference_normal_matched <- out$rate - out$z_normal_matched
  fbst_write_table(out, "power_vs_n")
  fbst_experiment_plot(out, "n", "rate", "delta", "power_vs_n", "n", "Rejection probability")
  invisible(out)
}

fbst_run_sensitivity <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  priors <- rbind(data.frame(category = "weak", N0 = 2, mu0 = FBST_SCENARIOS$weak_mu),
    data.frame(category = "informative", N0 = 50, mu0 = FBST_SCENARIOS$informative_mu),
    data.frame(category = "conflict", N0 = 50, mu0 = FBST_SCENARIOS$conflict_mu))
  if (fbst_profile() == "pilot") priors <- priors[c(1, 8, 12), ]
  rows <- list()
  for (n in fbst_experiment_sizes(FBST_SCENARIOS$sensitivity_n)) {
    zp <- fbst_z_pvalues(n, n)
    for (i in seq_len(nrow(priors))) {
      pr <- priors[i, ]
      prior <- fbst_mu_prior(pr$mu0, pr$N0)
      design <- fbst_get_design(n, n, prior)
      reject <- design$EV <= design$calibration$kstar
      level <- sum(fbst_sampling_mass(n, n, theta1, theta1) * reject)
      for (delta in FBST_SCENARIOS$deltas) {
        mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta)
        rate <- sum(mass * reject)
        z_normal <- fbst_z_power_normal(n, n, theta1, theta1 + delta, level)
        rows[[length(rows) + 1L]] <- data.frame(pr, n = n, delta = delta,
          theta1 = theta1, as.list(fbst_design_fields(design)),
          lambda0 = prior[1], lambda1 = prior[2], lambda2 = prior[3],
          rate = rate, matched_nominal = level, z_normal_matched = z_normal,
          z_finite_matched = sum(mass[zp <= level]), difference = rate - z_normal)
      }
    }
  }
  out <- do.call(rbind, rows)
  average <- aggregate(cbind(rate, alpha_star, beta_star, risk, difference,
      z_normal_matched, z_finite_matched) ~ category + n + delta, out, mean)
  fbst_write_table(out, "prior_sensitivity_individual")
  fbst_write_table(average, "prior_sensitivity")
  fbst_experiment_plot(subset(average, delta == .1), "n", "rate", "category",
                       "prior_sensitivity", "n", "Rejection probability at Delta = 0.10")
  invisible(out)
}

fbst_run_paired <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  rows <- list()
  for (n in fbst_experiment_sizes(FBST_SCENARIOS$paired_n)) {
    design <- fbst_get_design(n, n, FBST_PRIORS$KL)
    reject <- design$EV <= design$calibration$kstar
    for (psi in fbst_association_grid()) {
      nominal <- sum(fbst_sampling_mass(n, n, theta1, theta1, psi) * reject)
      for (delta in FBST_SCENARIOS$deltas) {
        p <- fbst_cell_probs(theta1, theta1 + delta, psi)
        disc <- fbst_discordants(n, theta1, theta1 + delta, psi)
        mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta, psi)
        rows[[length(rows) + 1L]] <- data.frame(n = n, theta1 = theta1,
          delta = delta, psi = psi, overlap = n, as.list(p),
          as.list(fbst_design_fields(design)), rate_fbst = sum(mass * reject),
          mcnemar_matched_nominal = nominal,
          mcnemar_matched = fbst_mcnemar_rate(disc, nominal),
          mcnemar_alpha_nominal = design$calibration$alpha,
          mcnemar_at_alpha = fbst_mcnemar_rate(disc, design$calibration$alpha),
          mcnemar_05_nominal = .05, mcnemar_at_05 = fbst_mcnemar_rate(disc, .05))
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "paired_vs_mcnemar")
  fbst_experiment_plot(subset(out, delta == .1), "n", "rate_fbst", "psi",
                       "paired_vs_mcnemar", "n", "Rejection probability at Delta = 0.10")
  invisible(out)
}

fbst_interval_grid <- function(n1, n2, prior, family = "olkin_liu") {
  rows <- lapply(0:n1, function(x1) {
    fbst_cache_compute("interval_row", list(n1 = n1, n2 = n2, x1 = x1,
      prior = prior, family = family, interval = "equal_tailed_95"), function() {
      t(vapply(0:n2, function(x2) {
        s <- posterior_summary(x1, n1, x2, n2, prior, family = family)
        c(lo = s$lo, hi = s$hi, mean = s$mean)
      }, c(lo = 0, hi = 0, mean = 0)))
    })
  })
  list(lo = do.call(rbind, lapply(rows, function(r) r[, "lo"])),
       hi = do.call(rbind, lapply(rows, function(r) r[, "hi"])),
       mean = do.call(rbind, lapply(rows, function(r) r[, "mean"])))
}

fbst_run_coverage <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  rows <- list()
  for (n in fbst_experiment_sizes(FBST_SCENARIOS$paired_n)) {
    intervals <- fbst_interval_grid(n, n, FBST_PRIORS$KL)
    design <- fbst_get_design(n, n, FBST_PRIORS$KL)
    for (psi in fbst_association_grid()) for (delta in FBST_SCENARIOS$coverage_deltas) {
      mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta, psi)
      wald <- fbst_wald_properties(fbst_discordants(n, theta1, theta1 + delta, psi), n, delta)
      rows[[length(rows) + 1L]] <- data.frame(n = n, theta1 = theta1, delta = delta,
        psi = psi, overlap = n, as.list(fbst_design_fields(design)),
        coverage = sum(mass * (intervals$lo <= delta & intervals$hi >= delta)),
        width = sum(mass * (intervals$hi - intervals$lo)),
        wald_coverage = wald["coverage"], wald_width = wald["width"],
        rejection = sum(mass * (design$EV <= design$calibration$kstar)))
    }
  }
  out <- do.call(rbind, rows)
  average <- aggregate(cbind(coverage, width, wald_coverage, wald_width) ~ n + psi, out, mean)
  fbst_write_table(out, "delta_coverage_by_delta")
  fbst_write_table(average, "delta_coverage")
  fbst_experiment_plot(average, "n", "coverage", "psi", "delta_coverage", "n", "Coverage")
  invisible(out)
}

# This benchmark is inexpensive and can be reproduced independently of the
# posterior interval grid, which is the costly part of the coverage study.
fbst_run_wald_benchmark <- function() {
  rows <- list()
  theta1 <- FBST_SCENARIOS$baseline
  for (n in fbst_experiment_sizes(FBST_SCENARIOS$paired_n)) {
    for (psi in fbst_association_grid()) for (delta in FBST_SCENARIOS$coverage_deltas) {
      properties <- fbst_wald_properties(fbst_discordants(n, theta1, theta1 + delta, psi), n, delta)
      rows[[length(rows) + 1L]] <- data.frame(n = n, theta1 = theta1,
        delta = delta, psi = psi, coverage = properties["coverage"], width = properties["width"])
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "paired_wald_by_delta")
  fbst_write_table(aggregate(cbind(coverage, width) ~ n + psi, out, mean), "paired_wald_average")
  invisible(out)
}

fbst_boundary_priors <- function() {
  priors <- FBST_PRIORS[c("KL", "informative", "conflict")]
  for (mu in FBST_SCENARIOS$weak_mu) priors[[paste0("weak_mu", mu)]] <- fbst_mu_prior(mu, 2)
  if (fbst_profile() == "pilot") priors <- priors[c("KL", "weak_mu0.1")]
  priors
}

fbst_boundary_sizes <- function() {
  if (fbst_profile() == "pilot") matrix(c(6, 6, 6, 8, 8, 6), ncol = 2, byrow = TRUE)
  else FBST_SCENARIOS$boundary_sizes
}

fbst_null_grid <- function(reject, n1, n2) {
  theta <- sort(unique(c(FBST_SCENARIOS$boundary_theta,
      seq(0, 1, length.out = 201), 10^seq(-5, -1, length.out = 31),
      1 - 10^seq(-5, -1, length.out = 31))))
  evaluate <- function(t) sum(fbst_sampling_mass(n1, n2, t, t, overlap = 0) * reject)
  rates <- vapply(theta, evaluate, 0)
  # Refine around every local grid maximum; report a grid maximum only.
  peaks <- which(rates >= c(-Inf, head(rates, -1L)) & rates >= c(tail(rates, -1L), -Inf))
  refinement <- unique(unlist(lapply(peaks, function(i)
    seq(theta[max(1, i - 1L)], theta[min(length(theta), i + 1L)], length.out = 21))))
  theta <- sort(unique(c(theta, refinement)))
  data.frame(theta = theta, rate = vapply(theta, evaluate, 0))
}

fbst_run_boundaries <- function() {
  priors <- fbst_boundary_priors()
  sizes <- fbst_boundary_sizes()
  null_rows <- alt_rows <- maxima <- singular <- list()
  alternatives <- data.frame(theta1 = c(.001, .01, .05, .4, .4, .6, .95, .99, .999),
    theta2 = c(.01, .05, .01, .5, .6, .4, .99, .95, .99))
  for (i in seq_len(nrow(sizes))) for (name in names(priors)) {
    n1 <- sizes[i, 1]; n2 <- sizes[i, 2]; prior <- priors[[name]]
    design <- fbst_get_design(n1, n2, prior)
    reject <- design$EV <= design$calibration$kstar
    grid <- fbst_null_grid(reject, n1, n2)
    base <- data.frame(n1 = n1, n2 = n2, prior = name,
                       as.list(fbst_design_fields(design)))
    null_rows[[length(null_rows) + 1L]] <- cbind(base, grid)
    best <- grid[grid$rate == max(grid$rate), , drop = FALSE]
    maxima[[length(maxima) + 1L]] <- cbind(base, best, n_grid = nrow(grid),
                                            scope = "maximum_on_refined_grid")
    S <- outer(0:n1, 0:n2, "+")
    singular_set <- prior[2] + prior[3] + S - 2 < 0 |
                   prior[1] + n1 + n2 - S - 2 < 0
    singular[[length(singular) + 1L]] <- data.frame(n1 = n1, n2 = n2, prior = name,
      n_singular = sum(singular_set), mass_null = sum(design$pH[singular_set]),
      mass_alternative = sum(design$pA[singular_set]),
      all_ev_one = all(design$EV[singular_set] == 1))
    for (j in seq_len(nrow(alternatives))) {
      t1 <- alternatives$theta1[j]; t2 <- alternatives$theta2[j]
      mass <- fbst_sampling_mass(n1, n2, t1, t2, overlap = 0)
      alt_rows[[length(alt_rows) + 1L]] <- cbind(base, alternatives[j, ],
        delta = t2 - t1, psi = 1, overlap = 0, rate = sum(mass * reject),
        z_finite_05 = sum(mass[fbst_z_pvalues(n1, n2) <= .05]))
    }
  }
  # The manuscript's two singular masses do not require e-value enumeration.
  a <- unname(FBST_PRIORS$KL[1])
  checks <- do.call(rbind, lapply(c(50L, 150L), function(n) data.frame(n = n,
    mass_null_formula = 2 * exp(lbeta(a, a + 2 * n) - lbeta(a, a)) +
      2 * n * exp(lbeta(a + 2 * n - 1, a + 1) - lbeta(a, a)),
    mass_null_enumerated = {
      p <- predictive_H(n, n, FBST_PRIORS$KL)
      sum(p[cbind(c(1, n + 1, n + 1, n), c(1, n + 1, n, n + 1))])
    })))
  stopifnot(max(abs(checks$mass_null_formula - checks$mass_null_enumerated)) < 1e-10)
  fbst_write_table(do.call(rbind, null_rows), "boundary_null_grid")
  fbst_write_table(do.call(rbind, maxima), "boundary_grid_maxima")
  fbst_write_table(do.call(rbind, alt_rows), "boundary_alternatives")
  fbst_write_table(do.call(rbind, singular), "boundary_singular_mass")
  fbst_write_table(checks, "kl_singular_mass")
  invisible(checks)
}

fbst_run_references <- function() {
  priors <- fbst_boundary_priors()
  sizes <- fbst_boundary_sizes()
  scenarios <- data.frame(theta1 = c(.001, .01, .05, .4, .95, .99, .999,
                                      .001, .05, .4, .6, .95, .999),
                          theta2 = c(.001, .01, .05, .4, .95, .99, .999,
                                      .01, .01, .5, .4, .99, .99))
  rows <- list()
  for (i in seq_len(nrow(sizes))) for (name in names(priors)) {
    n1 <- sizes[i, 1]; n2 <- sizes[i, 2]
    for (reference in c("flat", "prior")) {
      design <- fbst_get_design(n1, n2, priors[[name]], reference = reference)
      reject <- design$EV <= design$calibration$kstar
      for (j in seq_len(nrow(scenarios))) {
        s <- scenarios[j, ]
        rows[[length(rows) + 1L]] <- data.frame(n1 = n1, n2 = n2, prior = name,
          reference = reference, null_prior = "marginal", overlap = 0, psi = 1,
          s, as.list(fbst_design_fields(design)),
          rate = sum(fbst_sampling_mass(n1, n2, s$theta1, s$theta2, overlap = 0) * reject))
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "reference_sensitivity")
  invisible(out)
}

fbst_run_unequal <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  sizes <- fbst_boundary_sizes()
  if (fbst_profile() != "pilot") sizes <- unique(rbind(sizes,
    cbind(150, 150 - round(150 * c(0, .1, .2, .3, .5)))))
  rows <- list()
  for (i in seq_len(nrow(sizes))) {
    n1 <- sizes[i, 1]; n2 <- sizes[i, 2]
    design <- fbst_get_design(n1, n2, FBST_PRIORS$KL)
    reject <- design$EV <= design$calibration$kstar
    for (overlap in unique(c(0L, floor(min(n1, n2) / 2), min(n1, n2)))) {
      for (psi in fbst_association_grid()) for (delta in c(-.2, -.1, 0, .1, .2)) {
        mass <- fbst_sampling_mass(n1, n2, theta1, theta1 + delta, psi, overlap)
        disc <- fbst_discordants(overlap, theta1, theta1 + delta, psi)
        rows[[length(rows) + 1L]] <- data.frame(n1 = n1, n2 = n2,
          overlap = overlap, unpaired_pre = n1 - overlap, unpaired_post = n2 - overlap,
          theta1 = theta1, delta = delta, psi = psi,
          as.list(fbst_design_fields(design)), rate_fbst = sum(mass * reject),
          mcnemar_nominal = .05,
          rate_mcnemar = if (overlap > 0) fbst_mcnemar_rate(disc, .05) else NA_real_,
          z_finite_05 = sum(mass[fbst_z_pvalues(n1, n2) <= .05]))
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "unequal_samples")
  invisible(out)
}

# Restriction on the difference: beta(a1+a2-1,a0-1) times (1+t)^-lambda.
# Restriction on log odds: beta(a1+a2,a0) times the same smooth factor.
# Beta-quantile integration handles integrable endpoint singularities without
# replacing the mathematical domain by an arbitrary truncated interval.
fbst_null_tail <- function(prior, threshold, null = c("marginal", "arc", "lor")) {
  null <- match.arg(null)
  prior <- unname(prior)
  if (null == "marginal") return(pbeta(threshold, prior[2], prior[1], lower.tail = FALSE))
  A <- prior[2] + prior[3] - as.numeric(null == "arc")
  B <- prior[1] - as.numeric(null == "arc")
  if (A <= 0 || B <= 0) return(NA_real_)
  lambda <- sum(prior)
  shift <- -lambda * log1p(A / (A + B))
  density <- function(u) exp(-lambda * log1p(qbeta(u, A, B, lower.tail = FALSE)) - shift)
  # A beta-weighted Gaussian rule integrates the smooth tilt directly. The
  # full quantile transform is less regular at zero for concentrated priors.
  denominator <- NA_real_
  relative_error <- Inf
  for (G in c(32L, 64L, 128L, 256L)) {
    quadrature <- statmod::gauss.quad.prob(G, dist = "beta", alpha = A, beta = B)
    value <- sum(quadrature$weights * exp(-lambda * log1p(quadrature$nodes) - shift))
    if (is.finite(denominator)) relative_error <- abs(value / denominator - 1)
    if (relative_error < 1e-10) break
    denominator <- value
  }
  if (relative_error >= 1e-10) stop("Null-prior quadrature did not converge")
  denominator <- value
  upper <- pbeta(threshold, A, B, lower.tail = FALSE)
  if (upper == 0) return(0)
  # Scale the tail to [0,1] before integration so even a tiny numerator is
  # evaluated to relative precision instead of passing an absolute tolerance.
  upper * integrate(function(v) density(upper * v), 0, 1, rel.tol = 1e-9,
                    abs.tol = 0, subdivisions = 1000L)$value / denominator
}

fbst_run_restricted_prior <- function() {
  priors <- c(FBST_PRIORS[c("KL", "informative", "conflict")],
    list(illustrative = c(a0 = 1.1, a1 = 1, a2 = 1),
         finite_boundary = c(a0 = 2, a1 = 1, a2 = 1)))
  rows <- list()
  for (name in names(priors)) for (null in c("marginal", "arc", "lor")) {
    prior <- priors[[name]]
    proper <- null != "arc" || (prior[1] > 1 && sum(prior[2:3]) > 1)
    for (threshold in c(.9, .99)) {
      rows[[length(rows) + 1L]] <- data.frame(prior = name, null = null,
        lambda0 = prior[1], lambda1 = prior[2], lambda2 = prior[3],
        proper = proper, threshold = threshold,
        tail_probability = if (proper) fbst_null_tail(prior, threshold, null) else NA_real_)
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "restricted_prior")
  invisible(out)
}

fbst_step_cdf <- function(EV, mass) {
  o <- order(as.vector(EV))
  ev <- as.vector(EV)[o]
  total <- cumsum(as.vector(mass)[o])
  last <- !duplicated(ev, fromLast = TRUE)
  data.frame(k = ev[last], alpha = total[last])
}

fbst_run_null_priors <- function() {
  rows <- curves <- list()
  for (name in c("KL", "informative", "conflict")) {
    prior <- FBST_PRIORS[[name]]
    nulls <- if (prior[1] > 1 && sum(prior[2:3]) > 1) c("arc", "lor") else "lor"
    for (n in fbst_experiment_sizes(FBST_SCENARIOS$null_comparison_n)) {
      design <- fbst_get_design(n, n, prior)
      for (null in nulls) {
        restricted_design <- fbst_get_design(n, n, prior, null = null)
        restricted <- restricted_design$pH
        # Refinement can select different cells under the two null priors.
        # Both CDFs must use the same evidence table: retain the more precise
        # value at each cell from the union of their numerical refinements.
        EV <- design$EV
        errors <- attr(EV, "error_estimate")
        other_errors <- attr(restricted_design$EV, "error_estimate")
        better <- other_errors < errors
        EV[better] <- restricted_design$EV[better]
        errors[better] <- other_errors[better]
        marginal <- fbst_step_cdf(EV, design$pH)
        alt <- fbst_step_cdf(EV, restricted)
        stopifnot(identical(alt$k, marginal$k))
        diff <- abs(alt$alpha - marginal$alpha)
        local <- marginal$k >= .1 & marginal$k <= .7
        # Include the plateau spanning k=.1, whose jump can precede .1.
        idx <- findInterval(.1, marginal$k)
        if (idx > 0) local[idx] <- TRUE
        curves[[length(curves) + 1L]] <- data.frame(prior = name, n = n,
          null = null, k = marginal$k, alpha_marginal = marginal$alpha,
          alpha_restricted = alt$alpha)
        cal <- adaptive_cutoff(EV, restricted, design$pA)
        marginal_cal <- adaptive_cutoff(EV, design$pH, design$pA)
        for (k in c(.1, .3, .5, .7)) {
          rows[[length(rows) + 1L]] <- data.frame(prior = name, n = n, null = null, k = k,
            alpha_marginal = sum(design$pH[EV <= k]),
            alpha_restricted = sum(restricted[EV <= k]),
            sup_all_jumps = max(diff), sup_01_07 = max(c(0, diff[local])),
            k_star_marginal = marginal_cal$kstar, k_star_restricted = cal$kstar,
            evidence_max_error = max(errors),
            evidence_unresolved = sum(attr(design$EV, "unresolved") &
                                        attr(restricted_design$EV, "unresolved")))
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "arc_vs_marginal")
  # Complete curves can be very large; retain all jump points in CSV/RDS.
  full_curves <- do.call(rbind, curves)
  write.csv(full_curves, file.path(fbst_output_dir(), "null_prior_alpha_curves.csv"), row.names = FALSE)
  saveRDS(full_curves, file.path(fbst_output_dir(), "null_prior_alpha_curves.rds"))
  invisible(out)
}

fbst_run_independent_priors <- function() {
  theta1 <- FBST_SCENARIOS$baseline
  rows <- list()
  prior_rows <- list()
  for (name in c("KL", "informative")) {
    prior <- FBST_PRIORS[[name]]
    for (family in c("olkin_liu", "independent")) {
      p_near <- delta_cdf(.05, 0, 0, 0, 0, prior, family = family) -
                delta_cdf(-.05, 0, 0, 0, 0, prior, family = family)
      prior_rows[[length(prior_rows) + 1L]] <- data.frame(prior = name,
         family = family, delta_band = .05, probability = p_near)
    }
    for (n in fbst_experiment_sizes(sort(unique(c(FBST_SCENARIOS$independent_n, 50L, 200L))))) {
      biv <- fbst_get_design(n, n, prior, family = "olkin_liu")
      ind <- fbst_get_design(n, n, prior, family = "independent")
      rb <- biv$EV <= biv$calibration$kstar
      ri <- ind$EV <= ind$calibration$kstar
      for (delta in FBST_SCENARIOS$deltas) for (psi in c(.2, 1, 5)) {
        mass <- fbst_sampling_mass(n, n, theta1, theta1 + delta, psi)
        rows[[length(rows) + 1L]] <- data.frame(prior = name, n = n,
          theta1 = theta1, delta = delta, psi = psi, overlap = n,
          k_biv = biv$calibration$kstar, alpha_biv = biv$calibration$alpha,
          beta_biv = biv$calibration$beta, k_ind = ind$calibration$kstar,
          alpha_ind = ind$calibration$alpha, beta_ind = ind$calibration$beta,
          num_status_biv = fbst_design_fields(biv)$num_status,
          num_status_ind = fbst_design_fields(ind)$num_status,
          calibration_status_biv = fbst_design_fields(biv)$calibration_status,
          calibration_status_ind = fbst_design_fields(ind)$calibration_status,
          evidence_max_error_biv = fbst_design_fields(biv)$evidence_max_error,
          evidence_max_error_ind = fbst_design_fields(ind)$evidence_max_error,
          evidence_unresolved_biv = fbst_design_fields(biv)$evidence_unresolved,
          evidence_unresolved_ind = fbst_design_fields(ind)$evidence_unresolved,
          rejection_biv = sum(mass * rb), rejection_ind = sum(mass * ri),
          disagreement = sum(mass * (rb != ri)), mean_abs_ev_diff = sum(mass * abs(biv$EV - ind$EV)))
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "independent_prior_by_n")
  fbst_write_table(do.call(rbind, prior_rows), "prior_probability_near_null")
  fbst_write_table(subset(out, delta == .2 & psi == 1), "independent_evalue_differences")
  invisible(out)
}

fbst_run_application_size <- function() {
  # These are published aggregate counts, not invented individual records.
  groups <- data.frame(group = c("CC-TV", "CC-NoTV", "NoCC-TV", "Control"),
                        n = c(383L, 380L, 416L, 421L), x1 = c(118L, 128L, 145L, 159L))
  if (fbst_profile() == "pilot") {
    groups$n <- 6L
    groups$x1 <- 2L
    groups$group <- paste0("PILOT-", groups$group)
  }
  rows <- list()
  for (i in seq_len(nrow(groups))) for (name in c("KL", "informative", "conflict")) {
    g <- groups[i, ]
    design <- fbst_get_design(g$n, g$n, FBST_PRIORS[[name]])
    theta <- g$x1 / g$n
    rows[[length(rows) + 1L]] <- data.frame(group = g$group, prior = name,
      n = g$n, theta_null = theta, psi = 1, as.list(fbst_design_fields(design)),
      pointwise_type_I = sum(fbst_sampling_mass(g$n, g$n, theta, theta) *
                             (design$EV <= design$calibration$kstar)))
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "application_size")
  invisible(out)
}

# Independent, inexpensive checks of the sampling layer used by every study.
fbst_validate_sampling <- function() {
  for (psi in c(.2, 1, 3, 25)) for (overlap in 0:3) {
    p <- fbst_cell_probs(.3, .6, psi)
    q <- fbst_sampling_mass(3, 5, .3, .6, psi, overlap)
    stopifnot(max(abs(rowSums(q) - dbinom(0:3, 3, .3))) < 1e-12,
              max(abs(colSums(q) - dbinom(0:5, 5, .6))) < 1e-12)
    cov <- sum(q * outer(0:3 - 3 * .3, 0:5 - 5 * .6))
    stopifnot(abs(cov - overlap * (p[4] - .3 * .6)) < 1e-12)
  }
  p <- fbst_cell_probs(.3, .6, 5)
  brute <- matrix(0, 4, 4)
  for (n00 in 0:3) for (n01 in 0:(3 - n00)) for (n10 in 0:(3 - n00 - n01)) {
    n11 <- 3 - n00 - n01 - n10
    brute[n10 + n11 + 1, n01 + n11 + 1] <- brute[n10 + n11 + 1, n01 + n11 + 1] +
      dmultinom(c(n00, n01, n10, n11), 3, p)
  }
  stopifnot(max(abs(brute - fbst_sampling_mass(3, 3, .3, .6, 5))) < 1e-12,
            fbst_mcnemar_rate(fbst_discordants(0, .3, .6, 5), 1) == 0,
            fbst_z_pvalues(3, 5)[1, 1] == 1,
            fbst_z_pvalues(3, 5)[4, 6] == 1)
  invisible(TRUE)
}
