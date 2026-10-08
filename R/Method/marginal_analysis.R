# Posterior summaries and adaptive decisions for one row of marginal counts.
analyze_margins <- function(row, prior, calibrate = TRUE) {
  args <- list(
    x1 = row$x1,
    n1 = row$n1,
    x2 = row$x2,
    n2 = row$n2,
    prior = prior
  )
  fit <- cached_calculation(
    "application_posterior",
    args,
    function() {
      list(
        evalue = do.call(fbst_evalue, args),
        summary = do.call(posterior_summary, args)
      )
    }
  )
  calibration <- list(kstar = NA_real_, alpha = NA_real_, beta = NA_real_)

  if (calibrate) {
    calibration <- calibrate_test(row$n1, row$n2, prior)$calibration
    distance <- abs(as.numeric(fit$evalue) - calibration$kstar)

    if (distance <= attr(fit$evalue, "error_estimate")) {
      fit$evalue <- do.call(fbst_evalue, c(args, list(rel_tol = 1e-10, abs_tol = 1e-12)))
    }
  }

  summary <- fit$summary

  data.frame(
    item = row$item,
    sample = row$sample,
    n1 = row$n1,
    n2 = row$n2,
    delta_mean = summary$mean,
    delta_lo = summary$lo,
    delta_hi = summary$hi,
    P_increase = summary$P_increase,
    P_decrease = summary$P_decrease,
    evalue = as.numeric(fit$evalue),
    k_star = calibration$kstar,
    alpha_star = calibration$alpha,
    beta_star = calibration$beta,
    reject = as.numeric(fit$evalue) <= calibration$kstar
  )
}
