# Frequentist comparators: McNemar, paired Wald interval and two-proportion z-test.

mcnemar_rejection_rate <- function(discordants, level) {
  sum(discordants$mass[discordants$n01 + discordants$n10 > 0 &
    discordants$pvalue <= level])
}

wald_interval_properties <- function(discordants, n, delta) {
  if (n == 0) {
    return(c(coverage = NA_real_, width = NA_real_))
  }

  estimate <- (discordants$n01 - discordants$n10) / n
  half <- qnorm(.975) / n * sqrt(pmax(0, discordants$n01 +
    discordants$n10 - n * estimate^2))
  c(
    coverage = sum(discordants$mass * (estimate - half <= delta &
      estimate + half >= delta)),
    width = sum(discordants$mass * 2 * half)
  )
}

# Two-sided pooled z test. At all failures/all successes, the variance and
# observed difference are both zero; define p=1 and do not reject.
z_test_pvalue_grid <- function(n1, n2) {
  x1 <- row(matrix(0, n1 + 1L, n2 + 1L)) - 1L
  x2 <- col(x1) - 1L
  pooled <- (x1 + x2) / (n1 + n2)
  se <- sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))
  z <- matrix(0, n1 + 1L, n2 + 1L)
  ok <- se > 0
  z[ok] <- abs(x2[ok] / n2 - x1[ok] / n1) / se[ok]
  2 * pnorm(-z)
}

# Normal approximation to the z-test power (unpooled alternative SE).
# A common nominal level does not imply equal finite-sample type-I errors.
z_test_power <- function(n1, n2, theta1, theta2, level) {
  if (level == 0) {
    return(0)
  }

  if (level == 1) {
    return(1)
  }

  se <- sqrt(theta1 * (1 - theta1) / n1 + theta2 * (1 - theta2) / n2)

  if (se == 0) {
    return(as.numeric(theta1 != theta2))
  }

  d <- (theta2 - theta1) / se
  z <- qnorm(1 - level / 2)
  pnorm(d - z) + pnorm(-d - z)
}

mcnemar_pvalue <- function(n01, n10) {
  if (anyNA(c(n01, n10))) {
    return(NA_real_)
  }

  if (n01 + n10 == 0) {
    return(1)
  }

  stats::pchisq((n01 - n10)^2 / (n01 + n10), df = 1, lower.tail = FALSE)
}

z_test_pvalue <- function(x1, n1, x2, n2) {
  pooled <- (x1 + x2) / (n1 + n2)
  se <- sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))

  if (se == 0) {
    return(1)
  }

  2 * stats::pnorm(-abs((x2 / n2 - x1 / n1) / se))
}
