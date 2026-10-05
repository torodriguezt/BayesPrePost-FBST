# Appendix / SI tab:arcmarg: prior under H restricted to the null line against the
# marginal definition (priorH), and the tail probabilities of the restricted prior
# quoted in the appendix (e.g. 0.61 for t > 0.9 with lambda = (1.1, 1, 1)).
# Outputs: results/restricted_prior.csv, results/arc_vs_marginal.csv.

# Restriction on the difference: beta(a1+a2-1,a0-1) times (1+t)^-lambda.
# Restriction on log odds: beta(a1+a2,a0) times the same smooth factor.
# Beta-quantile integration handles integrable endpoint singularities without
# replacing the mathematical domain by an arbitrary truncated interval.
fbst_null_tail <- function(prior, threshold, null = c("marginal", "arc", "lor")) {
  null <- match.arg(null)
  prior <- unname(prior)
  if (null == "marginal") {
    return(pbeta(threshold, prior[2], prior[1], lower.tail = FALSE))
  }
  A <- prior[2] + prior[3] - as.numeric(null == "arc")
  B <- prior[1] - as.numeric(null == "arc")
  if (A <= 0 || B <= 0) {
    return(NA_real_)
  }
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
  if (upper == 0) {
    return(0)
  }
  # Scale the tail to [0,1] before integration so even a tiny numerator is
  # evaluated to relative precision instead of passing an absolute tolerance.
  upper * integrate(function(v) density(upper * v), 0, 1,
    rel.tol = 1e-9,
    abs.tol = 0, subdivisions = 1000L
  )$value / denominator
}

fbst_run_restricted_prior <- function() {
  priors <- c(
    FBST_PRIORS[c("KL", "informative", "conflict")],
    list(
      illustrative = c(a0 = 1.1, a1 = 1, a2 = 1),
      finite_boundary = c(a0 = 2, a1 = 1, a2 = 1)
    )
  )
  rows <- list()
  for (name in names(priors)) {
    for (null in c("marginal", "arc", "lor")) {
      prior <- priors[[name]]
      proper <- null != "arc" || (prior[1] > 1 && sum(prior[2:3]) > 1)
      for (threshold in c(.9, .99)) {
        rows[[length(rows) + 1L]] <- data.frame(
          prior = name, null = null,
          lambda0 = prior[1], lambda1 = prior[2], lambda2 = prior[3],
          proper = proper, threshold = threshold,
          tail_probability = if (proper) fbst_null_tail(prior, threshold, null) else NA_real_
        )
      }
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
  rows <- list()
  for (name in c("KL", "informative", "conflict")) {
    prior <- FBST_PRIORS[[name]]
    nulls <- if (prior[1] > 1 && sum(prior[2:3]) > 1) c("arc", "lor") else "lor"
    for (n in FBST_SCENARIOS$null_comparison_n) {
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
        cal <- adaptive_cutoff(EV, restricted, design$pA)
        marginal_cal <- adaptive_cutoff(EV, design$pH, design$pA)
        for (k in c(.1, .3, .5, .7)) {
          rows[[length(rows) + 1L]] <- data.frame(
            prior = name, n = n, null = null, k = k,
            alpha_marginal = sum(design$pH[EV <= k]),
            alpha_restricted = sum(restricted[EV <= k]),
            sup_all_jumps = max(diff), sup_01_07 = max(c(0, diff[local])),
            k_star_marginal = marginal_cal$kstar, k_star_restricted = cal$kstar,
            evidence_max_error = max(errors),
            evidence_unresolved = sum(attr(design$EV, "unresolved") &
              attr(restricted_design$EV, "unresolved"))
          )
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  fbst_write_table(out, "arc_vs_marginal")
  invisible(out)
}

fbst_run_null_prior <- function() {
  fbst_run_restricted_prior()
  fbst_run_null_priors()
}
