# Posterior parameters, null supremum and quadrature primitives.

.fbst_par <- function(x1, n1, x2, n2, prior, family = "olkin_liu", power = 1) {
  family <- match.arg(family, c("olkin_liu", "independent"))
  if (length(prior) != 3 || any(!is.finite(prior)) ||
    any(prior <= 0)) {
    stop("Three positive prior parameters required")
  }
  if (any(!is.finite(c(x1, n1, x2, n2, power))) || power <= 0 || min(x1, n1 - x1, x2, n2 -
    x2) < 0) {
    stop("Invalid counts or likelihood power")
  }
  a0 <- prior[1]
  a1 <- prior[2]
  a2 <- prior[3]
  c(
    A1 = a1 + power * x1, B1 = a0 + if (family == "olkin_liu") a2 + power * (n1 - x1) else power * (n1 - x1),
    A2 = a2 + power * x2, B2 = a0 + if (family == "olkin_liu") a1 + power * (n2 - x2) else power * (n2 - x2),
    lambda = if (family == "olkin_liu") sum(prior) else 0
  )
}
.fbst_xlog <- function(a, x) if (a == 0) rep(0, length(x)) else a * log(x)
.fbst_logkernel <- function(x, y, p) {
  .fbst_xlog(p[1] - 1, x) + .fbst_xlog(p[2] - 1, 1 - x) +
    .fbst_xlog(p[3] - 1, y) + .fbst_xlog(p[4] - 1, 1 - y) - .fbst_xlog(p[5], 1 - x * y)
}
.fbst_roots <- function(a, b, c) {
  if (abs(a) < 1e-14 * max(1, abs(b), abs(c))) {
    return(if (abs(b) < 1e-14) numeric() else -c / b)
  }
  d <- b * b - 4 * a * c
  if (d < 0) {
    return(numeric())
  }
  q <- -.5 * (b + if (b >= 0) sqrt(d) else -sqrt(d))
  if (q == 0) {
    return(-b / (2 * a))
  }
  unique(c(q / a, c / q))
}
# Global supremum includes finite boundary limits and all stationary roots.
null_supremum <- function(x1, n1, x2, n2, prior, reference = "flat", family = "olkin_liu", power = 1) {
  reference <- match.arg(reference, c("flat", "prior"))
  p <- .fbst_par(x1, n1, x2, n2, prior, family, power)
  r <- if (reference == "prior") {
    c(
      power * x1, power * (n1 - x1), power * x2, power * (n2 - x2),
      0
    )
  } else {
    c(p[1:4] - 1, p[5])
  }
  # Direct expressions prevent cancellation from turning an exact zero exponent
  # (for example a0=1, F=3, power=1/3) into a spurious negative exponent.
  S <- x1 + x2
  F <- n1 + n2 - S
  N <- n1 + n2
  lam <- r[5]
  if (reference == "prior") {
    u <- power * S
    v <- power * F
    quadratic <- -power * N
  } else {
    u <- prior[2] + prior[3] + power * S - 2
    if (family == "olkin_liu") {
      v <- prior[1] + power * F - 2
      quadratic <- 4 - power * N
    } else {
      v <- 2 * prior[1] + power * F - 2
      quadratic <- 4 - (2 * prior[1] + prior[2] + prior[3]) - power * N
    }
  }
  if (u < 0 || v < 0) {
    return(list(log_sup = Inf, t = NA_real_, singular = TRUE, u = u, v = v))
  }
  candidates <- c(0, 1, .fbst_roots(quadratic, -v - lam, u))
  candidates <- sort(unique(candidates[candidates >= 0 & candidates <= 1]))
  vals <- .fbst_xlog(u, candidates) + .fbst_xlog(v, 1 - candidates) - lam * log1p(candidates)
  j <- which.max(vals)
  list(
    log_sup = unname(vals[j]), t = candidates[j], singular = FALSE, u = unname(u), v = unname(v),
    candidates = candidates, log_values = vals
  )
}
.fbst_gauss <- function(A, B, G) {
  if (!requireNamespace("statmod",
    quietly = TRUE
  )) {
    stop("Package 'statmod' is required: install.packages('statmod')")
  }
  key <- paste(sprintf("%.17g", c(A, B, G)), collapse = ":")
  if (exists(key, envir = .fbst_axis_cache, inherits = FALSE)) {
    return(get(key, envir = .fbst_axis_cache, inherits = FALSE))
  }
  q <- statmod::gauss.quad.prob(G, dist = "beta", alpha = A, beta = B)
  if (any(q$nodes <= 0 | q$nodes >= 1) || any(q$weights < 0)) stop("Invalid Gaussian quadrature nodes")
  assign(key, q, envir = .fbst_axis_cache)
  q
}
.fbst_grid <- function(p, G) {
  q1 <- .fbst_gauss(p[1], p[2], G)
  q2 <- .fbst_gauss(p[3], p[4], G)
  L <- -p[5] * log1p(-outer(q1$nodes, q2$nodes))
  shift <- max(L)
  W <- exp(L - shift) * outer(q1$weights, q2$weights)
  z <- sum(W)
  list(x = q1$nodes, y = q2$nodes, w = W / z, logZ = lbeta(p[1], p[2]) + lbeta(p[3], p[4]) + shift + log(z))
}
.fbst_normalizer <- function(p, rel_tol = 1e-8, max_G = 1024) {
  if (p[5] == 0) {
    return(list(logZ = lbeta(p[1], p[2]) + lbeta(p[3], p[4]), log_lower = lbeta(p[1], p[2]) + lbeta(
      p[3],
      p[4]
    ) - 1e-12, error = 0, G = 0))
  }
  if (exists("fbst_normalizer_series", envir = .fbst_env, inherits = FALSE)) {
    z <- .fbst_env$fbst_normalizer_series(unname(p), rel_tol / 10, 1048576L)
    if (isTRUE(z$unresolved)) stop("Posterior series normalizer unresolved: ", z$error)
    return(z)
  }
  old <- NA_real_
  for (G in unique(pmin(max_G, 2^(5:11)))) {
    g <- .fbst_grid(p, G)
    err <- if (is.na(old)) Inf else abs(expm1(g$logZ - old))
    if (err < rel_tol) {
      return(list(
        logZ = g$logZ, log_lower = lbeta(p[1], p[2]) + lbeta(p[3], p[4]) - 1e-12, error = err,
        G = G
      ))
    }
    old <- g$logZ
    if (G >= max_G) break
  }
  stop(sprintf("Posterior normalization did not converge: relative difference %.3g at G=%d", err, G))
}
.fbst_integrate <- function(f, lo = 0, hi = 1, rel_tol = 1e-8, abs_tol = 1e-11, .depth = 0L) {
  if (lo >= hi) {
    return(0)
  }
  z <- integrate(f, lo, hi, rel.tol = rel_tol, abs.tol = abs_tol, subdivisions = 600L, stop.on.error = FALSE)
  if (z$message != "OK") {
    if (.depth < 4L) {
      mid <- (lo + hi) / 2
      return(.fbst_integrate(f, lo, mid, rel_tol, abs_tol / 2, .depth + 1L) + .fbst_integrate(
        f, mid, hi,
        rel_tol, abs_tol / 2, .depth + 1L
      ))
    }
    stop("Adaptive integration failed after interval subdivision: ", z$message)
  }
  z$value
}
