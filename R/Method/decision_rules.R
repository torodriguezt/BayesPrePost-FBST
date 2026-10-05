# Adaptive cutoff: the k minimising the averaged errors alpha(k) + beta(k).

adaptive_cutoff <- function(EV, pH, pA, tol = .002) {
  if (!identical(dim(EV), dim(pH)) || !identical(dim(EV), dim(pA)) || anyNA(EV) || any(!is.finite(EV)) ||
    any(EV < 0 | EV > 1)) {
    stop("Incompatible/invalid evidence and predictive matrices")
  }
  if (abs(sum(pH) - 1) > 1e-7 || abs(sum(pA) - 1) > 1e-7) {
    stop("Predictives must sum to one before calibration")
  }
  ev <- as.numeric(EV)
  o <- order(ev)
  v <- ev[o]
  ends <- which(c(diff(v) != 0, TRUE))
  ks <- v[ends]
  alpha <- cumsum(as.numeric(pH)[o])[ends]
  beta <- sum(pA) - cumsum(as.numeric(pA)[o])[ends]
  beta <- pmax(0, beta)
  if (min(ks) > 0) {
    ks <- c(0, ks)
    alpha <- c(0, alpha)
    beta <- c(sum(pA), beta)
  }
  # At k=1 the last region is closed on the right; every other interval is [lo,hi).
  hi <- c(ks[-1], 1)
  risk <- alpha + beta
  rmin <- min(risk)
  j <- which.min(risk)
  intervals <- data.frame(
    lower = ks,
    upper = hi,
    lower_closed = TRUE,
    upper_closed = c(rep(FALSE, length(ks) - 1L), TRUE),
    alpha = alpha,
    beta = beta,
    risk = risk
  )
  mini <- which(risk == rmin)
  near <- which(risk <= rmin + tol)
  list(
    kstar = ks[j], alpha = alpha[j], beta = beta[j], risk = rmin,
    tie_break = "smallest attainable cutoff in [0,1] attaining the numerically smallest risk; numerical near-ties listed separately",
    ks = ks,
    alpha_k = alpha, beta_k = beta, risk_k = risk,
    intervals = intervals,
    minimizers = intervals[mini, , drop = FALSE],
    numerical_ties = intervals[abs(risk - rmin) <= 1e-12, , drop = FALSE],
    near_optimal = intervals[near, , drop = FALSE],
    k_range = range(c(intervals$lower[near], intervals$upper[near])), alpha_range = range(alpha[near]),
    no_rejection = list(k = -Inf, alpha = 0, beta = sum(pA), risk = sum(pA))
  )
}
