# Cached design construction and adaptive calibration refinement.
# Depends on the numerical core and the cache/table infrastructure.

fbst_call <- function(fun, args) {
  f <- get(fun, mode = "function")
  allowed <- names(formals(f))
  if (!"..." %in% allowed) {
    bad <- setdiff(names(args), allowed)
    if (length(bad)) stop(fun, " does not accept: ", paste(bad, collapse = ", "))
  }
  do.call(f, args)
}
fbst_check_mass <- function(p, label, tol = FBST_NUMERICS$predictive_mass_tol) {
  if (any(!is.finite(p)) || any(p < 0)) stop(label, ": invalid probability masses")
  err <- abs(sum(p) - 1)
  if (err > tol) {
    stop(sprintf(
      "%s sum %.12g differs from one by %.3g (tolerance %.3g); no renormalisation applied",
      label, sum(p), err, tol
    ))
  }
  invisible(err)
}
fbst_refine_design <- function(EV, pH, pA, n1, n2, prior, reference, family) {
  err <- attr(EV, "error_estimate")
  if (is.null(err)) stop("Enumeration must carry numerical error diagnostics")
  initial <- attr(EV, "diagnostics")
  initial_error <- max(err)
  unresolved <- attr(EV, "unresolved")
  upper <- attr(EV, "upper_bound")
  if (is.null(upper)) upper <- matrix(FALSE, n1 + 1, n2 + 1)
  if (any(!is.finite(EV)) || any(!is.finite(err)) || any(err < 0)) {
    stop("Evidence enumeration has non-finite values or error diagnostics")
  }
  # Grid integration can overshoot the probability range by a few ulps. Project
  # such estimates to [0,1] and retain the projection distance in the error
  # diagnostic, so calibration can use only mathematically admissible values.
  outside <- EV < 0 | EV > 1
  clipped_count <- sum(outside)
  clipped_distance <- if (clipped_count) abs(EV[outside] - pmin(pmax(EV[outside], 0), 1)) else numeric()
  if (clipped_count) {
    err[outside] <- pmax(err[outside], clipped_distance)
    EV[outside] <- pmin(pmax(EV[outside], 0), 1)
    message(
      "Projected ", clipped_count, " grid e-values onto [0,1]; maximum correction ",
      format(max(clipped_distance), scientific = TRUE), " retained in error diagnostics"
    )
  }
  cdf <- function(x, w, k) {
    o <- order(x)
    cumulative <- c(0, cumsum(as.numeric(w)[o]))
    cumulative[findInterval(k, as.numeric(x)[o]) + 1L]
  }
  candidates <- function(EV, err, upper) {
    low <- pmax(EV - err, 0)
    high <- pmin(EV + err, 1)
    low[upper] <- 0
    high[upper] <- EV[upper]
    ks <- sort(unique(c(as.numeric(EV), as.numeric(low), as.numeric(high), 0, 1)))
    lower <- cdf(high, pH, ks) + sum(pA) - cdf(low, pA, ks)
    higher <- cdf(low, pH, ks) + sum(pA) - cdf(high, pA, ks)
    plausible <- ks[lower <= min(higher) + FBST_NUMERICS$near_risk]
    cells <- which(err > 1e-7 & low <= max(plausible) & high >= min(plausible), arr.ind = TRUE)
    if (n1 == n2 && prior[2] == prior[3] && nrow(cells)) {
      cells <- cells[cells[, 1] <= cells[, 2], , drop = FALSE]
    }
    cells
  }
  refined <- 0L
  max_attempts <- 4L
  attempts <- matrix(0L, n1 + 1, n2 + 1)
  original_EV <- EV
  trace <- list()
  symmetric <- n1 == n2 && prior[2] == prior[3]
  # Each pass consumes at least one previously unused attempt. This bound also
  # allows new competitive cells to enter after earlier values are updated.
  max_passes <- max_attempts * length(EV)
  describe_cells <- function(cells) {
    paste(vapply(seq_len(nrow(cells)), function(i) {
      r <- cells[i, 1]
      c <- cells[i, 2]
      sprintf(
        "(n1,n2,x1,x2)=(%d,%d,%d,%d), attempts=%d, estimated error=%.3g",
        n1, n2, r - 1, c - 1, attempts[r, c], err[r, c]
      )
    }, character(1)), collapse = "; ")
  }
  for (pass in seq_len(max_passes)) {
    cells <- candidates(EV, err, upper)
    if (!nrow(cells)) break
    exhausted <- cells[attempts[cells] >= max_attempts, , drop = FALSE]
    if (nrow(exhausted)) {
      stop(
        "Adaptive calibration refinement remains incomplete after ",
        max_attempts, " precision levels for ", nrow(exhausted), " competitive points: ",
        describe_cells(exhausted), ". Evaluations at each tolerance are cached; ",
        "estimated errors are not certified bounds."
      )
    }
    message(
      "Refining ", nrow(cells), " e-values near competitive cutoffs (", n1, ",", n2,
      ", pass ", pass, "; per-cell precision levels ",
      paste(sort(unique(attempts[cells] + 1L)), collapse = ","), "/", max_attempts, ")"
    )
    for (i in seq_len(nrow(cells))) {
      r <- cells[i, 1]
      c <- cells[i, 2]
      attempt <- attempts[r, c] + 1L
      previous <- EV[r, c]
      # The first call retains the original cache key. Persistent cells receive
      # genuinely tighter controls and therefore a distinct cached evaluation.
      args <- list(
        x1 = r - 1, n1 = n1, x2 = c - 1, n2 = n2, prior = unname(prior),
        reference = reference, family = family, method = "adaptive",
        rel_tol = 1e-8 / 10^(attempt - 1L), abs_tol = 1e-10 / 100^(attempt - 1L)
      )
      z <- fbst_cache_compute("adaptive_evidence", args, function() {
        tryCatch(fbst_call("evalue", args), error = function(e) {
          stop(
            "Adaptive evidence failed for (n1,n2,x1,x2)=(",
            n1, ",", n2, ",", r - 1, ",", c - 1, "), attempt ", attempt,
            " (rel_tol=", args$rel_tol, ", abs_tol=", args$abs_tol, "): ", conditionMessage(e)
          )
        })
      })
      value <- as.numeric(z)
      estimated_error <- attr(z, "error_estimate")
      if (length(value) != 1L || !is.finite(value) || value < 0 || value > 1 ||
        length(estimated_error) != 1L || !is.finite(estimated_error) || estimated_error < 0) {
        stop(
          "Invalid adaptive evidence or error estimate for (n1,n2,x1,x2)=(",
          n1, ",", n2, ",", r - 1, ",", c - 1, "), attempt ", attempt
        )
      }
      attempts[r, c] <- attempt
      trace[[length(trace) + 1L]] <- data.frame(
        pass = pass, x1 = r - 1, x2 = c - 1,
        attempt = attempt, rel_tol = args$rel_tol, abs_tol = args$abs_tol,
        evalue = value, error_estimate = estimated_error,
        change_from_grid = abs(value - original_EV[r, c]),
        change_between_precisions = if (attempt > 1L) abs(value - previous) else NA_real_
      )
      EV[r, c] <- value
      err[r, c] <- estimated_error
      unresolved[r, c] <- FALSE
      upper[r, c] <- FALSE
      if (symmetric) {
        EV[c, r] <- EV[r, c]
        err[c, r] <- err[r, c]
        unresolved[c, r] <- FALSE
        upper[c, r] <- FALSE
        attempts[c, r] <- attempt
      }
      refined <- refined + 1L
    }
  }
  remaining <- candidates(EV, err, upper)
  if (nrow(remaining)) {
    stop(
      "Adaptive calibration refinement reached its finite pass limit: ",
      describe_cells(remaining), ". Evaluations at each tolerance are cached for diagnosis."
    )
  }
  trace <- if (length(trace)) {
    do.call(rbind, trace)
  } else {
    data.frame(
      pass = integer(),
      x1 = integer(), x2 = integer(), attempt = integer(), rel_tol = numeric(), abs_tol = numeric(),
      evalue = numeric(), error_estimate = numeric(), change_from_grid = numeric(),
      change_between_precisions = numeric()
    )
  }
  attr(EV, "error_estimate") <- err
  attr(EV, "unresolved") <- unresolved
  attr(EV, "upper_bound") <- upper
  attr(EV, "refinement_attempts") <- attempts
  attr(EV, "refinement_trace") <- trace
  attr(EV, "diagnostics") <- modifyList(initial, list(
    initial_max_error = initial_error,
    initial_grid_values_outside_probability_range = clipped_count,
    max_grid_range_projection = if (clipped_count) max(clipped_distance) else 0,
    max_error = max(err), unresolved_count = sum(unresolved),
    refinement_count = refined, remaining_grid_unresolved = sum(unresolved),
    refinement_max_attempts = max(attempts), refinement_attempt_limit = max_attempts,
    refinement_passes = if (nrow(trace)) max(trace$pass) else 0L,
    max_change_between_precisions = if (any(is.finite(trace$change_between_precisions))) {
      max(trace$change_between_precisions, na.rm = TRUE)
    } else {
      NA_real_
    },
    error_estimates_are_certified_bounds = FALSE,
    refinement = "adaptive near competitive cutoffs; persistent cells use tighter tolerances; remaining evidence retains grid diagnostics"
  ))
  cal <- adaptive_cutoff(EV, pH, pA, tol = FBST_NUMERICS$near_risk)
  cal$numeric_status <- "adaptive_refinement_of_competitive_cutoffs"
  list(EV = EV, calibration = cal)
}
fbst_get_design <- function(n1, n2, prior, reference = "flat", family = "olkin_liu", null = "marginal", ...) {
  opts <- modifyList(FBST_NUMERICS$enum, list(...))
  base <- list(n1 = n1, n2 = n2, prior = unname(prior), family = family)
  evargs <- c(base, list(reference = reference), opts)
  EV <- fbst_cache_compute("evidence", evargs, function() {
    checkpoint <- file.path(fbst_output_dir(), "cache", paste0(
      "partial-evidence-",
      fbst_hash(list(arguments = evargs, code = fbst_code_hashes())), ".rds"
    ))
    fbst_call("enumerate_evalues", c(evargs, list(checkpoint_file = checkpoint)))
  })
  paargs <- c(base, FBST_NUMERICS$predictive)
  pA <- fbst_cache_compute("predictive_A", paargs, function() fbst_call("predictive_A", paargs))
  phargs <- list(n1 = n1, n2 = n2, prior = unname(prior), null = null)
  pH <- fbst_cache_compute("predictive_H", phargs, function() fbst_call("predictive_H", phargs))
  fbst_check_mass(pA, "Alternative predictive")
  fbst_check_mass(pH, "Null predictive")
  result <- fbst_cache_compute(
    "calibration_refined", list(evidence = evargs, null = phargs, near_risk = FBST_NUMERICS$near_risk),
    function() fbst_refine_design(EV, pH, pA, n1, n2, prior, reference, family)
  )
  EV <- result$EV
  calibration <- result$calibration
  list(
    n1 = n1, n2 = n2, EV = EV, pH = pH, pA = pA, calibration = calibration,
    parameters = c(base, list(reference = reference, null = null)),
    numeric_controls = opts, predictive_sums = c(H = sum(pH), A = sum(pA))
  )
}

# Calibration fields and numerical status of a design, for the result tables.
fbst_design_fields <- function(design) {
  errors <- attr(design$EV, "error_estimate")
  unresolved <- attr(design$EV, "unresolved")
  count <- if (is.null(unresolved)) NA_integer_ else sum(unresolved)
  max_error <- if (is.null(errors)) NA_real_ else max(errors)
  status <- if (is.na(count)) {
    "numerical_diagnostics_unavailable"
  } else if (count > 0) {
    "estimated_grid_requires_review"
  } else {
    "quadrature_tolerance_met"
  }
  calibration_status <- design$calibration$numeric_status
  if (is.null(calibration_status)) calibration_status <- "grid_calibration_only"
  list(
    k_star = design$calibration$kstar, alpha_star = design$calibration$alpha,
    beta_star = design$calibration$beta, risk = design$calibration$risk,
    num_status = status, calibration_status = calibration_status,
    evidence_max_error = max_error, evidence_unresolved = count
  )
}

# Level of the asymptotically equivalent Wald test, eq. (eq:zlevel).
fbst_alpha_z <- function(k) stats::pchisq(-2 * log(k), df = 1, lower.tail = FALSE)

# Best test for a*alpha + b*beta over ALL tests (Pericchi and Pereira, 2016): reject H
# when a m_H(x) < b m(x). Its risk is sum(min(a m_H, b m)), eq. (eq:bayesmin).
fbst_bayes_rule <- function(pH, pA, a = 1, b = 1) {
  reject <- a * pH < b * pA
  list(
    reject = reject, alpha = sum(pH[reject]), beta = sum(pA[!reject]),
    risk = sum(pmin(a * pH, b * pA))
  )
}

# Cutoff, near-optimal range, attained errors and excess over the Bayes rule.
fbst_design_summary <- function(d) {
  cal <- d$calibration
  bayes <- fbst_bayes_rule(d$pH, d$pA)
  ev_rule <- d$EV <= cal$kstar
  data.frame(
    k_star = cal$kstar, k_near_lo = cal$k_range[1], k_near_hi = cal$k_range[2],
    alpha_star = cal$alpha, beta_star = cal$beta, total_star = cal$risk,
    alpha_z = fbst_alpha_z(cal$kstar),
    R_min = bayes$risk, alpha_bayes = bayes$alpha, beta_bayes = bayes$beta,
    excess = cal$risk - bayes$risk,
    # Probabilities, under m_H and under m, that the two rules decide differently.
    disagree_H = sum(d$pH[ev_rule != bayes$reject]),
    disagree_A = sum(d$pA[ev_rule != bayes$reject])
  )
}
