# Calibrate the rejection cutoff using null and alternative predictive errors.

weighted_cdf <- function(values, weights, cutoffs) {
  order <- order(values)
  cumulative <- c(0, cumsum(as.numeric(weights)[order]))
  cumulative[findInterval(cutoffs, as.numeric(values)[order]) + 1L]
}

# Find e-values whose integration error could change the optimal cutoff.
competitive_cells <- function(EV, error, upper_bound, pH, pA, symmetric) {
  low <- pmax(EV - error, 0)
  high <- pmin(EV + error, 1)
  low[upper_bound] <- 0
  high[upper_bound] <- EV[upper_bound]

  cutoffs <- sort(unique(c(as.numeric(EV), as.numeric(low), as.numeric(high), 0, 1)))
  lower_risk <- weighted_cdf(high, pH, cutoffs) + sum(pA) - weighted_cdf(low, pA, cutoffs)
  upper_risk <- weighted_cdf(low, pH, cutoffs) + sum(pA) - weighted_cdf(high, pA, cutoffs)
  plausible <- cutoffs[lower_risk <= min(upper_risk) + NUMERICS$near_risk]

  cells <- which(
    error > 1e-7 & low <= max(plausible) & high >= min(plausible),
    arr.ind = TRUE
  )

  if (symmetric) {
    cells <- cells[cells[, 1] <= cells[, 2], , drop = FALSE]
  }

  cells
}

refine_calibration <- function(EV, pH, pA, n1, n2, prior, reference, family) {
  error <- attr(EV, "error_estimate")
  unresolved <- attr(EV, "unresolved")
  upper_bound <- attr(EV, "upper_bound")
  symmetric <- n1 == n2 && prior[2] == prior[3]

  if (any(!is.finite(EV)) || any(!is.finite(error)) || any(error < 0)) {
    stop("Evidence integration did not produce finite values and error estimates")
  }

  # Keep roundoff corrections within the integration error estimate.
  outside <- EV < 0 | EV > 1
  clipped <- pmin(pmax(EV[outside], 0), 1)
  error[outside] <- pmax(error[outside], abs(EV[outside] - clipped))
  EV[outside] <- clipped

  max_attempts <- 4L
  attempts <- matrix(0L, n1 + 1, n2 + 1)

  repeat {
    cells <- competitive_cells(EV, error, upper_bound, pH, pA, symmetric)

    if (nrow(cells) == 0) break

    if (any(attempts[cells] >= max_attempts)) {
      stop("Calibration did not converge after four precision levels")
    }

    message("Refining ", nrow(cells), " e-values near the cutoff (", n1, ", ", n2, ")")

    for (i in seq_len(nrow(cells))) {
      r <- cells[i, 1]
      c <- cells[i, 2]
      attempt <- attempts[r, c] + 1L
      args <- list(
        x1 = r - 1,
        n1 = n1,
        x2 = c - 1,
        n2 = n2,
        prior = unname(prior),
        reference = reference,
        family = family,
        method = "adaptive",
        rel_tol = 1e-8 / 10^(attempt - 1L),
        abs_tol = 1e-10 / 100^(attempt - 1L)
      )
      value <- cached_calculation(
        "adaptive_evidence", args,
        function() do.call(fbst_evalue, args)
      )

      EV[r, c] <- as.numeric(value)
      error[r, c] <- attr(value, "error_estimate")
      unresolved[r, c] <- FALSE
      upper_bound[r, c] <- FALSE
      attempts[r, c] <- attempt

      if (symmetric) {
        EV[c, r] <- EV[r, c]
        error[c, r] <- error[r, c]
        unresolved[c, r] <- FALSE
        upper_bound[c, r] <- FALSE
        attempts[c, r] <- attempt
      }
    }
  }

  attr(EV, "error_estimate") <- error
  attr(EV, "unresolved") <- unresolved
  attr(EV, "upper_bound") <- upper_bound
  attr(EV, "diagnostics") <- NULL

  list(
    EV = EV,
    calibration = adaptive_cutoff(EV, pH, pA, tol = NUMERICS$near_risk)
  )
}

calibrate_test <- function(
    n1, n2, prior, reference = "flat", family = "olkin_liu",
    null = "marginal", ...) {
  settings <- modifyList(NUMERICS$enum, list(...))
  design_args <- list(n1 = n1, n2 = n2, prior = unname(prior), family = family)
  evidence_args <- c(design_args, list(reference = reference), settings)

  EV <- cached_calculation("evidence", evidence_args, function() {
    checkpoint <- file.path(results_directory(), "cache", paste0(
      "partial-evidence-",
      object_hash(list(arguments = evidence_args, code = numerical_code_hashes())),
      ".rds"
    ))
    do.call(enumerate_fbst_evalues, c(evidence_args, list(checkpoint_file = checkpoint)))
  })

  alternative_args <- c(design_args, NUMERICS$predictive)
  pA <- cached_calculation("alternative_predictive", alternative_args, function() {
    do.call(alternative_predictive, alternative_args)
  })

  null_args <- list(n1 = n1, n2 = n2, prior = unname(prior), null = null)
  pH <- cached_calculation("null_predictive", null_args, function() {
    do.call(null_predictive, null_args)
  })

  result <- cached_calculation(
    "calibration_refined",
    list(evidence = evidence_args, null = null_args, near_risk = NUMERICS$near_risk),
    function() refine_calibration(EV, pH, pA, n1, n2, prior, reference, family)
  )

  list(
    n1 = n1,
    n2 = n2,
    EV = result$EV,
    pH = pH,
    pA = pA,
    calibration = result$calibration
  )
}
