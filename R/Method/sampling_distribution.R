# Exact sampling distributions for paired and partially paired counts.
# psi is the within-subject odds ratio.

# Cells are ordered 00,01,10,11. Endpoint margins are allowed; association is
# not identified there, but the unique compatible cell distribution is valid.
transition_probabilities <- function(theta1, theta2, psi = 1) {
  stopifnot(
    length(theta1) == 1L, length(theta2) == 1L,
    is.finite(theta1), is.finite(theta2), is.finite(psi),
    theta1 >= 0, theta1 <= 1, theta2 >= 0, theta2 <= 1, psi > 0
  )
  lo <- max(0, theta1 + theta2 - 1)
  hi <- min(theta1, theta2)
  p11 <- if (hi == lo) {
    lo
  } else if (psi == 1) {
    theta1 * theta2
  } else {
    uniroot(
      function(z) {
        z * (1 - theta1 - theta2 + z) -
          psi * (theta1 - z) * (theta2 - z)
      }, c(lo, hi),
      tol = 1e-14
    )$root
  }

  p <- c(
    p00 = 1 - theta1 - theta2 + p11,
    p01 = theta2 - p11,
    p10 = theta1 - p11,
    p11 = p11
  )
  # Only roundoff at exactly degenerate margins can produce a negative zero.
  p[p < 0] <- 0
  p
}

# Exact distribution of the margins: the overlap subjects contribute iid
# bivariate Bernoulli pairs; the remaining n1-overlap and n2-overlap subjects
# contribute independent Bernoulli observations at their respective stages.
# Matrix rows index X1 and columns index X2, both beginning at zero.
sampling_distribution <- function(
    n1,
    n2,
    theta1,
    theta2,
    psi = 1,
    overlap = min(n1, n2)) {
  stopifnot(
    n1 >= 0, n2 >= 0, n1 == as.integer(n1), n2 == as.integer(n2),
    overlap >= 0, overlap <= min(n1, n2), overlap == as.integer(overlap)
  )
  p <- transition_probabilities(theta1, theta2, psi)

  if (psi == 1 || overlap == 0) {
    return(outer(dbinom(0:n1, n1, theta1), dbinom(0:n2, n2, theta2)))
  }

  mass <- matrix(1, 1L, 1L)

  for (j in seq_len(overlap)) {
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

  if (abs(sum(mass) - 1) > 5e-11) {
    stop("Sampling distribution lost probability mass")
  }

  mass
}

# Collapsing the concordant cells gives a trinomial distribution. Conditional
# on D discordant pairs, N01 is binomial(D, p01/(p01+p10)). This enumerates all
# discordant configurations without enumerating unnecessary concordant splits.
discordant_distribution <- function(n, theta1, theta2, psi) {
  stopifnot(n >= 0, n == as.integer(n))
  p <- transition_probabilities(theta1, theta2, psi)
  pd <- unname(p[2] + p[3])
  out <- lapply(0:n, function(d) {
    b <- 0:d
    w <- dbinom(d, n, pd) * dbinom(b, d, if (pd > 0) {
      p[2] / pd
    } else {
      0
    })
    statistic <- if (d == 0) {
      0
    } else {
      (2 * b - d)^2 / d
    }

    data.frame(
      n01 = b,
      n10 = d - b,
      mass = w,
      Q = statistic,
      pvalue = if (d == 0) {
        1
      } else {
        pchisq(statistic, 1, lower.tail = FALSE)
      }
    )
  })
  out <- do.call(rbind, out)
  stopifnot(abs(sum(out$mass) - 1) < 5e-11)
  out
}
