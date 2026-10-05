# Complement integration, FBST e-values and evidence enumeration.

# Conditional odds transform resolves the joint (1,1) corner without truncating.
.fbst_slice_integral <- function(x, lo, hi, p, shift, rel_tol, abs_tol) {
  if (lo >= hi) {
    return(0)
  }
  if (p[5] == 0) {
    return(exp((p[1] - 1) * log(x) + (p[2] - 1) * log1p(-x) + lbeta(p[3], p[4]) - shift) *
      if (lo == 0) {
        pbeta(hi, p[3], p[4])
      } else if (hi == 1) {
        pbeta(lo, p[3], p[4],
          lower.tail = FALSE
        )
      } else if (pbeta(lo, p[3], p[4]) > .5) {
        (pbeta(lo, p[3], p[4],
          lower.tail = FALSE
        ) - pbeta(hi, p[3], p[4], lower.tail = FALSE))
      } else {
        (pbeta(hi, p[3], p[4]) -
          pbeta(lo, p[3], p[4]))
      })
  }
  if (p[4] - p[5] > 0) {
    # If the conditional limit at x=1 remains integrable, original coordinates
    # avoid the narrow artificial peak introduced by the odds transformation.
    cc <- (p[1] - 1) * log(x) + (p[2] - 1) * log1p(-x) - shift
    direct_piece <- function(l, h) {
      if (p[3] < 1) {
        return(.fbst_integrate(function(z) {
          y <- h * z^(1 / p[3])
          exp(cc + p[3] * log(h) - log(p[3]) + (p[4] - 1) * log1p(-y) - p[5] * log1p(-x * y))
        }, (l / h)^p[3], 1, rel_tol, abs_tol))
      }
      b <- min(p[4], p[4] - p[5])
      if (h == 1 && b < 1) {
        return(.fbst_integrate(function(z) {
          t <- (1 - l) * z^(1 / b)
          exp(cc + p[4] * log1p(-l) - log(b) + (p[4] / b - 1) * log(z) + (p[3] - 1) * log1p(-t) - p[5] *
            log((1 - x) + x * t))
        }, 0, 1, rel_tol, abs_tol))
      }
      .fbst_integrate(function(y) {
        exp(cc + (p[3] - 1) * log(y) + (p[4] - 1) * log1p(-y) - p[5] *
          log1p(-x * y))
      }, l, h, rel_tol, abs_tol)
    }
    return(if (lo < .5 && hi == 1) direct_piece(lo, .5) + direct_piece(.5, 1) else direct_piece(lo, hi))
  }
  trans <- function(y) if (y == 1) 1 else if (y == 0) 0 else (1 - x) * y / (1 - x * y)
  vl <- trans(lo)
  vh <- trans(hi)
  cst <- (p[1] - 1) * log(x) + (p[2] + p[4] - p[5] - 1) * log1p(-x) - shift
  logrest <- function(v) cst + (p[5] - p[3] - p[4]) * log1p(-x + x * v)
  piece <- function(vl, vh) {
    if (p[3] < 1) {
      # v=vh*z^(1/A2) cancels the integrable power singularity exactly.
      return(.fbst_integrate(function(z) {
        v <- vh * z^(1 / p[3])
        exp(logrest(v) + (p[4] - 1) * log1p(-v) + p[3] * log(vh) - log(p[3]))
      }, (vl / vh)^p[3], 1, rel_tol, abs_tol))
    }
    if (vh == 1 && p[4] < 1) {
      return(.fbst_integrate(function(z) {
        w <- (1 - vl) * z^(1 / p[4])
        v <- 1 - w
        exp(logrest(v) + (p[3] - 1) * log(v) + p[4] * log1p(-vl) - log(p[4]))
      }, 0, 1, rel_tol, abs_tol))
    }
    .fbst_integrate(
      function(v) exp(logrest(v) + (p[3] - 1) * log(v) + (p[4] - 1) * log1p(-v)), vl, vh,
      rel_tol, abs_tol
    )
  }
  tryCatch(if (vl == 0 && vh == 1 && p[3] < 1 && p[4] < 1) piece(0, .5) + piece(.5, 1) else piece(vl, vh),
    error = function(e) {
      stop(
        conditionMessage(e), "; slice x=", format(x, digits = 17), ", y interval=",
        paste(format(c(lo, hi), digits = 17), collapse = ",")
      )
    }
  )
}
# Integrate the complement directly. Each slice is split at all stationary points.
.fbst_complement_intervals <- function(x, r, target) {
  cst <- .fbst_xlog(r[1], x) + .fbst_xlog(r[2], 1 - x)
  f <- function(y) cst + .fbst_xlog(r[3], y) + .fbst_xlog(r[4], 1 - y) - r[5] * log1p(-x * y) - target
  st <- .fbst_roots(x * (r[3] + r[4] - r[5]), -r[3] * (1 + x) - r[4] + r[5] * x, r[3])
  pts <- sort(unique(c(0, st[st > 0 & st < 1], 1)))
  roots <- numeric()
  for (j in seq_len(length(pts) - 1L)) {
    l <- pts[j]
    h <- pts[j + 1L]
    fl <- f(l)
    fh <- f(h)
    if (is.finite(fl) && fl == 0) roots <- c(roots, l)
    if (is.finite(fh) && fh == 0) roots <- c(roots, h)
    if (!is.nan(fl) && !is.nan(fh) && sign(fl) * sign(fh) < 0) {
      rr <- if (l == 0 && r[3] != 0) {
        lf <- f(exp(-744))
        if (sign(lf) * sign(fh) >= 0) {
          0
        } else {
          exp(uniroot(function(z) f(exp(z)), c(-744, log(h)),
            tol = 1e-12
          )$root)
        }
      } else {
        uniroot(f, c(l, h), tol = 1e-13)$root
      }
      roots <- c(roots, rr)
    }
  }
  b <- sort(unique(c(0, roots, 1)))
  ans <- matrix(numeric(), ncol = 2)
  for (j in seq_len(length(b) - 1L)) if (f(mean(b[j:(j + 1L)])) <= 0) ans <- rbind(ans, b[j:(j + 1L)])
  ans
}
evalue <- function(
    x1, n1, x2, n2, prior, G = 64, W = NULL, reference = "flat", family = "olkin_liu",
    power = 1,
    method = c("adaptive", "grid"), rel_tol = 1e-8, abs_tol = 1e-10, max_G = 1024,
    strict = TRUE) {
  method <- match.arg(method)
  reference <- match.arg(reference, c("flat", "prior"))
  p <- .fbst_par(x1, n1, x2, n2, prior, family, power)
  if (method == "adaptive" && p[1] / (p[1] + p[2]) > p[3] / (p[3] + p[4])) {
    # Exchange parameter coordinates AND their prior shapes; this is a
    # measure-preserving symmetry, not a success/failure complement.
    return(evalue(x2, n2, x1, n1, prior[c(1, 3, 2)],
      G = G, reference = reference, family = family,
      power = power, method = method, rel_tol = rel_tol, abs_tol = abs_tol, max_G = max_G, strict = strict
    ))
  }
  ns <- null_supremum(x1, n1, x2, n2, prior, reference, family, power)
  if (ns$singular) {
    return(structure(1, error_estimate = 0, log_evalue = 0, diagnostics = list(
      method = "analytic",
      singular = TRUE
    )))
  }
  r <- if (reference == "prior") {
    c(
      power * x1, power * (n1 - x1), power * x2, power * (n2 - x2),
      0
    )
  } else {
    c(p[1:4] - 1, p[5])
  }
  if (method == "grid") {
    old <- NA_real_
    for (g in unique(pmin(max_G, G * 2^(0:10)))) {
      q <- .fbst_grid(p, g)
      if (exists("fbst_grid_ev", envir = .fbst_env, inherits = FALSE)) {
        val <- .fbst_env$fbst_grid_ev(q$x, q$y, q$w, p, r, ns$log_sup)
      } else {
        L <- outer(q$x, q$y, Vectorize(function(x, y) .fbst_logkernel(x, y, c(r[1:4] + 1, r[5]))))
        val <- sum(q$w * (L <= ns$log_sup))
      }
      err <- if (is.na(old)) Inf else abs(val - old)
      bound <- Inf
      if (val < abs_tol && reference == "flat") {
        norm_bound <- .fbst_normalizer(p, rel_tol = rel_tol, max_G = max_G)
        # A finite partial sum of the positive series gives a lower normalizer;
        # hence exp(sup)/Z_lower bounds the complement (the square has area one).
        logbound <- ns$log_sup - norm_bound$log_lower
        bound <- max(exp(logbound), .Machine$double.xmin)
      }
      isbound <- val < abs_tol && is.finite(bound) && bound < abs_tol
      if (isbound) {
        val <- bound
        err <- bound
        break
      }
      if (err < abs_tol && val > 0) break
      old <- val
      if (g >= max_G) break
    }
    unresolved <- (!is.finite(err) || err > abs_tol || val == 0)
    if (strict && unresolved) stop(sprintf("E-value grid unresolved: difference %.3g at G=%d", err, g))
    return(structure(val, error_estimate = err, log_evalue = log(val), diagnostics = list(
      method = "grid",
      G = g, unresolved = unresolved, upper_bound = isbound,
      log_upper_bound = if (isbound) logbound else NA_real_
    )))
  }
  norm <- .fbst_normalizer(p, rel_tol = rel_tol / 5, max_G = max_G)
  # Flat-reference complement is bounded by exp(log_sup), retaining tiny probabilities.
  shift <- if (reference == "flat") ns$log_sup else norm$logZ
  fun <- function(xx) {
    vapply(xx, function(x) {
      intervals <- .fbst_complement_intervals(x, r, ns$log_sup)
      if (!nrow(intervals)) {
        return(0)
      }
      sum(vapply(seq_len(nrow(intervals)), function(j) {
        .fbst_slice_integral(
          x, intervals[j, 1],
          intervals[j, 2], p, shift, rel_tol, abs_tol / 10
        )
      }, numeric(1)))
    }, numeric(1))
  }
  z <- .fbst_integrate(fun, rel_tol = rel_tol, abs_tol = abs_tol)
  if (z <= 0) stop("Complement integral is numerically zero; increase numerical precision")
  logev <- log(z) + shift - norm$logZ
  if (logev < log(.Machine$double.xmin)) stop(sprintf("e-value below double range; log(e-value)=%.16g. Use log scale / arbitrary precision", logev))
  ev <- exp(logev)
  if (ev > 1 + max(1e-7, rel_tol * 10)) stop("e-value exceeds one: integration inconsistency")
  structure(min(ev, 1),
    error_estimate = ev * (norm$error + rel_tol) + abs_tol * exp(shift - norm$logZ), log_evalue = logev,
    diagnostics = list(
      method = "adaptive", normalizer_G = norm$G, normalizer_error = norm$error,
      singular = FALSE
    )
  )
}
enumerate_evalues <- function(n1, n2, prior, G = 64, reference = "flat", family = "olkin_liu", power = 1,
                              method = "grid", max_G = 256, abs_tol = 1e-5, rel_tol = 1e-8,
                              strict = FALSE, progress = FALSE, checkpoint_file = NULL) {
  stopifnot(n1 == as.integer(n1), n2 == as.integer(n2), n1 >= 0, n2 >= 0)
  started <- proc.time()[3]
  EV <- matrix(NA_real_, n1 + 1, n2 + 1)
  errors <- EV
  unresolved <- matrix(FALSE, n1 + 1, n2 + 1)
  upper_bound <- unresolved
  symmetric <- n1 == n2 && identical(unname(prior[2]), unname(prior[3]))
  signature <- list(
    version = FBST_CORE_VERSION, core_hash = .fbst_source_hash, n1 = n1, n2 = n2,
    prior = prior, G = G, reference = reference, family = family, power = power, method = method,
    max_G = max_G, abs_tol = abs_tol, rel_tol = rel_tol
  )
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    saved <- readRDS(checkpoint_file)
    if (identical(saved$signature, signature)) {
      EV <- saved$EV
      errors <- saved$errors
      unresolved <- saved$unresolved
      upper_bound <- saved$upper_bound
    }
  }
  for (i in 0:n1) {
    for (j in 0:n2) {
      if (!is.na(EV[i + 1, j + 1])) next
      if (symmetric && j < i) {
        EV[i + 1, j + 1] <- EV[j + 1, i + 1]
        errors[i + 1, j + 1] <- errors[j + 1, i + 1]
        unresolved[i + 1, j + 1] <- unresolved[j + 1, i + 1]
        upper_bound[i + 1, j + 1] <- upper_bound[j + 1, i + 1]
        next
      }
      z <- evalue(i, n1, j, n2, prior,
        G = G, reference = reference, family = family, power = power, method = method,
        max_G = max_G, abs_tol = abs_tol, rel_tol = rel_tol, strict = strict
      )
      EV[i + 1, j + 1] <- as.numeric(z)
      errors[i + 1, j + 1] <- attr(z, "error_estimate")
      unresolved[i + 1, j + 1] <- isTRUE(attr(z, "diagnostics")$unresolved)
      upper_bound[i + 1, j + 1] <- isTRUE(attr(z, "diagnostics")$upper_bound)
    }
    if (!is.null(checkpoint_file)) {
      dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
      tmp <- paste0(checkpoint_file, ".tmp")
      saveRDS(list(
        signature = signature, EV = EV, errors = errors, unresolved = unresolved,
        upper_bound = upper_bound
      ), tmp)
      if (!file.rename(tmp, checkpoint_file)) stop("Cannot save enumeration checkpoint")
    }
    if (progress && i %% 10 == 0) message("e-values: ", i, "/", n1)
  }
  attr(EV, "error_estimate") <- errors
  attr(EV, "unresolved") <- unresolved
  attr(EV, "upper_bound") <- upper_bound
  attr(EV, "diagnostics") <- list(
    method = method, G = G, max_G = max_G, max_error = max(errors),
    unresolved_count = sum(unresolved), elapsed = unname(proc.time()[3] - started)
  )
  EV
}
