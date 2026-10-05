# Prior predictive probabilities under the alternative and null.

# Positive quadrature over the exact independent representation:
# X~Beta(a1,a0), V~Beta(a2,a0+a1), Y=V/(1-X+X*V).
predictive_A <- function(
    n1, n2, prior, K = NULL, G = 64, max_G = 1024, abs_tol = 1e-8,
    family = "olkin_liu", strict = TRUE, method = c("series", "quadrature"), max_terms = 1048576) {
  .fbst_par(0, n1, 0, n2, prior, family)
  a0 <- prior[1]
  a1 <- prior[2]
  a2 <- prior[3]
  if (family == "independent") {
    v1 <- exp(lchoose(n1, 0:n1) + lbeta(a1 + 0:n1, a0 + n1 - 0:n1) - lbeta(a1, a0))
    v2 <- exp(lchoose(n2, 0:n2) + lbeta(a2 + 0:n2, a0 + n2 - 0:n2) - lbeta(a2, a0))
    ans <- outer(v1, v2)
    attr(ans, "normalization_error") <- abs(sum(ans) - 1)
    attr(ans, "error_estimate") <- 0
    return(ans)
  }
  method <- match.arg(method)
  if (method == "series" && exists("fbst_predictive_series", envir = .fbst_env, inherits = FALSE)) {
    z <- .fbst_env$fbst_predictive_series(n1, n2, as.numeric(prior), min(1e-10, abs_tol / 100), max_terms)
    ans <- z$probability
    ne <- abs(sum(ans) - 1)
    err <- sum(z$error)
    if (any(!is.finite(ans)) || any(ans < 0) || ne > abs_tol || z$unresolved_count > 0) {
      if (strict) stop(sprintf("Predictive series unresolved: normalization %.3g, L1 error estimate %.3g, %d cells", ne, err, z$unresolved_count))
    }
    attr(ans, "normalization_error") <- ne
    attr(ans, "error_estimate") <- err
    attr(ans, "diagnostics") <- list(
      method = "positive_3F2_series_Richardson6",
      max_terms = z$max_terms_used, L1_error_estimate = err, unresolved = z$unresolved_count > 0
    )
    return(ans)
  }
  old <- NULL
  for (g in unique(pmin(max_G, G * 2^(0:10)))) {
    q1 <- .fbst_gauss(a1, a0, g)
    q2 <- .fbst_gauss(a2, a0 + a1, g)
    B1 <- vapply(q1$nodes, function(x) dbinom(0:n1, n1, x), numeric(n1 + 1))
    conditional <- vapply(q1$nodes, function(x) {
      y <- q2$nodes / (1 - x + x * q2$nodes)
      as.numeric(vapply(0:n2, function(j) sum(q2$weights * dbinom(j, n2, y)), numeric(1)))
    }, numeric(n2 + 1))
    ans <- (B1 * rep(q1$weights, each = n1 + 1)) %*% t(conditional)
    normalization_error <- abs(sum(ans) - 1)
    err <- if (is.null(old)) Inf else sum(abs(ans - old))
    if (err < abs_tol) break
    old <- ans
    if (g >= max_G) break
  }
  if (normalization_error > 1e-9) stop("Alternative predictive failed normalization: ", normalization_error)
  if (strict &&
    err >= abs_tol) {
    stop(sprintf(
      "Alternative predictive quadrature unresolved: L1 change %.3g at G=%d",
      err, g
    ))
  }
  attr(ans, "normalization_error") <- normalization_error
  attr(ans, "error_estimate") <- err
  attr(ans, "diagnostics") <- list(
    method = "independent_prior_representation", G = g,
    L1_difference = err, unresolved = err >= abs_tol
  )
  ans
}
predictive_H <- function(n1, n2, prior, null = "marginal", G = 128, rel_tol = 1e-10) {
  null <- match.arg(null, c("marginal", "arc", "lor"))
  a0 <- prior[1]
  a1 <- prior[2]
  a2 <- prior[3]
  if (null == "marginal" && abs(a1 - a2) > 1e-12) stop("Main marginal null prior requires lambda1=lambda2")
  S <- outer(0:n1, 0:n2, "+")
  N <- n1 + n2
  LC <- outer(lchoose(n1, 0:n1), lchoose(n2, 0:n2), "+")
  if (null == "marginal") {
    ans <- exp(LC + lbeta(a1 + S, a0 + N - S) - lbeta(a1, a0))
  } else {
    ca <- a1 + a2 - if (null == "arc") 1 else 0
    cb <- a0 - if (null == "arc") 1 else 0
    if (ca <= 0 ||
      cb <= 0) {
      stop("Improper restriction on difference: requires lambda0>1 and lambda1+lambda2>1")
    }
    logI <- function(a, b) {
      q <- .fbst_gauss(a, b, G)
      lbeta(a, b) + log(sum(q$weights * exp(-sum(prior) * log1p(q$nodes))))
    }
    vals <- vapply(0:N, function(s) logI(ca + s, cb + N - s), numeric(1)) - logI(ca, cb)
    ans <- exp(LC + vals[S + 1])
  }
  err <- abs(sum(ans) - 1)
  if (err > max(1e-9, rel_tol * 10)) stop("Null predictive normalization failed: ", err)
  attr(ans, "normalization_error") <- err
  ans
}
