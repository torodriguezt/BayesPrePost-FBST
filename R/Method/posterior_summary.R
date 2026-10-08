# Posterior distribution of delta = theta2 - theta1 and posterior summaries.

# CDF of delta by direct integration, using the same Lebesgue reference throughout.
delta_cdf <- function(
    d,
    x1,
    n1,
    x2,
    n2,
    prior,
    family = "olkin_liu",
    power = 1,
    rel_tol = 1e-7,
    abs_tol = 1e-9,
    max_G = 1024) {
  p <- posterior_parameters(x1, n1, x2, n2, prior, family, power)
  norm <- posterior_normalizer(p, rel_tol / 10)

  one <- function(v) {
    if (v <= -1) {
      return(0)
    }

    if (v >= 1) {
      return(1)
    }

    f <- function(xx) {
      vapply(xx, function(x) {
        h <- min(1, x + v)

        if (h <= 0) {
          return(0)
        }

        integrate_posterior_slice(x, 0, h, p, norm$logZ, rel_tol / 5, abs_tol / 10)
      }, numeric(1))
    }

    split <- sort(unique(c(max(0, -v), min(1, 1 - v), 1)))
    split <- split[split >= max(0, -v)]
    ans <- sum(vapply(seq_len(length(split) - 1L), function(j) {
      adaptive_integral(
        f, split[j], split[j + 1],
        rel_tol, abs_tol
      )
    }, numeric(1)))

    if (ans < -abs_tol || ans > 1 + 1e-6) {
      stop("Delta CDF violates probability bounds")
    }

    min(1, max(0, ans))
  }

  vapply(d, one, numeric(1))
}

posterior_summary <- function(
    x1,
    n1,
    x2,
    n2,
    prior,
    G = 128,
    family = "olkin_liu",
    power = 1,
    method = c("adaptive", "grid"),
    rel_tol = 1e-7,
    abs_tol = 1e-9,
    max_G = 1024) {
  method <- match.arg(method)
  p <- posterior_parameters(x1, n1, x2, n2, prior, family, power)
  q <- posterior_grid(p, G)
  d <- outer(q$x, q$y, function(x, y) {
    y - x
  })
  mu <- sum(q$w * d)
  q2 <- posterior_grid(p, min(max_G, 2 * G))
  mu2 <- sum(q2$w * outer(q2$x, q2$y, function(x, y) {
    y - x
  }))
  norm <- posterior_normalizer(p, rel_tol / 10)
  pp <- p
  pp[1] <- pp[1] + 1
  m1 <- exp(posterior_normalizer(pp, rel_tol / 10)$logZ - norm$logZ)
  pp <- p
  pp[3] <- pp[3] + 1
  m2 <- exp(posterior_normalizer(pp, rel_tol / 10)$logZ - norm$logZ)
  mu2 <- m2 - m1
  mean_error <- norm$error * 3

  if (method == "grid") {
    describe <- function(q) {
      dd <- outer(q$x, q$y, function(x, y) {
        y - x
      })
      o <- order(dd)
      cum <- cumsum(q$w[o])
      c(
        lo = dd[o][which(cum >= .025)[1]],
        hi = dd[o][which(cum >= .975)[1]],
        P_increase = sum(q$w * (dd > 0)) + .5 * sum(q$w * (dd == 0)),
        P_decrease = sum(q$w * (dd < 0)) +
          .5 * sum(q$w * (dd == 0))
      )
    }

    coarse <- describe(q)
    fine <- describe(q2)
    component_errors <- c(mean = max(abs(mu - mu2), mean_error), abs(fine - coarse))

    return(list(
      mean = mu2,
      lo = unname(fine["lo"]),
      hi = unname(fine["hi"]),
      P_increase = unname(fine["P_increase"]),
      P_decrease = unname(fine["P_decrease"]),
      mean_error_estimate = unname(component_errors["mean"]),
      error_estimate = max(component_errors),
      component_error_estimates = component_errors,
      diagnostics = list(
        not_certified = TRUE,
        coarse_G = G,
        fine_G = min(max_G, 2 * G),
        error_method = "successive quadratures; not a mathematical error bound"
      ),
      method = "grid"
    ))
  }

  cdf <- function(d) {
    delta_cdf(
      d,
      x1,
      n1,
      x2,
      n2,
      prior,
      family,
      power,
      rel_tol,
      abs_tol,
      max_G
    )
  }

  qq <- function(v) {
    uniroot(
      function(d) {
        cdf(d) - v
      },
      c(-1, 1),
      tol = max(1e-8, abs_tol)
    )$root
  }

  list(
    mean = mu2,
    lo = qq(.025),
    hi = qq(.975),
    P_increase = delta_cdf(
      0, x2, n2, x1, n1, prior[c(1, 3, 2)],
      family, power, rel_tol, abs_tol, max_G
    ),
    P_decrease = cdf(0),
    error_estimate = max(mean_error, rel_tol),
    method = "adaptive"
  )
}

# Posterior mean of theta1: E(theta1 | x) = Z(A1+1,B1,A2,B2,lambda) / Z(A1,B1,A2,B2,lambda).
posterior_mean_theta1 <- function(x1, n1, x2, n2, prior) {
  p <- posterior_parameters(x1, n1, x2, n2, prior)
  exp(posterior_normalizer(p + c(1, 0, 0, 0, 0))$logZ - posterior_normalizer(p)$logZ)
}
