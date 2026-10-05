# Fit the vaping items (McCauley et al., 2023) under the three priors, for the
# linked students and for all students observed at each stage, with design
# effects D = 1, 1.5, 2, 3 (likelihood raised to 1/D, cutoff of D = 1 kept fixed).
# Outputs: results/vaping_transition.csv (Table 7, tab:vaping_data) and
# results/vaping_results.csv, read by tables.R, figures.R and the other studies.

fbst_app_mcnemar <- function(n01, n10) {
  # Always asymptotic, without continuity correction. No discordants => no rejection.
  if (n01 + n10 == 0) return(1)
  stats::pchisq((n01 - n10)^2 / (n01 + n10), df = 1, lower.tail = FALSE)
}

fbst_app_z <- function(x1, n1, x2, n2) {
  pooled <- (x1 + x2) / (n1 + n2)
  se <- sqrt(pooled * (1 - pooled) * (1 / n1 + 1 / n2))
  if (se == 0) return(1)
  2 * stats::pnorm(-abs((x2 / n2 - x1 / n1) / se))
}

# "complete": linked students; "observed": margins of all students at each stage.
fbst_app_args <- function(g, prior, D = 1, analysis = "observed") {
  if (analysis == "complete") {
    list(x1 = g$x1_cc, n1 = g$n_cc, x2 = g$x2_cc, n2 = g$n_cc, prior = prior, power = 1 / D)
  } else list(x1 = g$x1, n1 = g$n1, x2 = g$x2, n2 = g$n2, prior = prior, power = 1 / D)
}

fbst_app_fit <- function(g, prior, prior_key, D = 1, analysis = "observed", calibrate = TRUE) {
  args <- fbst_app_args(g, prior, D, analysis)
  t0 <- proc.time()[[3L]]
  fit <- fbst_cache_compute("application_posterior", args, function() {
    summary <- do.call(posterior_summary, args)
    if (is.null(summary$P_decrease)) summary$P_decrease <- do.call(delta_cdf, c(list(d = 0), args))
    list(evalue = do.call(evalue, args), summary = summary)
  })
  # Only the ORIGINAL integer-sample experiment has a sampling calibration.
  design <- if (calibrate) fbst_get_design(args$n1, args$n2, prior) else NULL
  cal <- if (is.null(design)) NULL else design$calibration
  scalar <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)[1L]
  s <- fit$summary
  ev_error <- scalar(attr(fit$evalue, "error_estimate"))
  if (!is.null(cal) && abs(as.numeric(fit$evalue) - cal$kstar) <= ev_error) {
    fit$evalue <- do.call(evalue, c(args, list(rel_tol = 1e-10, abs_tol = 1e-12)))
    ev_error <- scalar(attr(fit$evalue, "error_estimate"))
  }
  data.frame(
    study = g$study, group = g$group, analysis = analysis,
    prior_key = prior_key, lambda0 = unname(prior[1]), lambda1 = unname(prior[2]), lambda2 = unname(prior[3]),
    n1 = args$n1, n2 = args$n2, x1 = args$x1, x2 = args$x2,
    D = D, likelihood_power = 1 / D, reference = "flat", null = "marginal",
    ev = as.numeric(fit$evalue), ev_error_estimate = ev_error,
    posterior_error_estimate = scalar(s$error_estimate),
    delta_mean = scalar(s$mean), delta_lo = scalar(s$lo), delta_hi = scalar(s$hi),
    prob_gt = scalar(s$P_increase), prob_lt = scalar(s$P_decrease),
    k_star = scalar(cal$kstar), alpha_star_original = scalar(cal$alpha),
    beta_star_original = scalar(cal$beta), risk_original = scalar(cal$risk),
    reject = if (is.null(cal)) NA else as.numeric(fit$evalue) <= cal$kstar,
    decision_clear_of_ev_error = if (is.null(cal)) NA else abs(as.numeric(fit$evalue) - cal$kstar) > ev_error,
    calibration_status = if (is.null(cal)) "not_run" else if (D == 1) "enumerated" else "original_cutoff_fixed",
    calibration_numeric_status = if (is.null(cal)) "not_run" else cal$numeric_status,
    calibration_grid_unresolved = if (is.null(design)) NA_integer_ else sum(attr(design$EV, "unresolved")),
    sampling_errors_apply_to_D = D == 1 && !is.null(cal),
    p_z = fbst_app_z(args$x1, args$n1, args$x2, args$n2),
    p_mcn = if (g$n1 == g$n_cc && g$n2 == g$n_cc || analysis == "complete")
      fbst_app_mcnemar(g$n01, g$n10) else NA_real_,
    frequentist_level = 0.05, elapsed_seconds = proc.time()[[3L]] - t0,
    data_source = g$data_source
  )
}

fbst_vaping_fit <- function(g, prior, prior_key, D, analysis, calibrate = TRUE) {
  h <- g
  # fbst_app_fit evaluates McNemar for complete pairs; with unknown cells it receives
  # a placeholder and the result is reset to NA below.
  if (!isTRUE(g$cells_available)) h$n01 <- h$n10 <- 0
  r <- fbst_app_fit(h, prior, prior_key, D, analysis = analysis, calibrate = calibrate)
  if (!isTRUE(g$cells_available) || analysis != "complete") r$p_mcn <- NA_real_
  r$sample <- if (analysis == "complete") "linked" else "full"
  r$description <- g$description
  r
}

# The six KL-optimal designs of the application (three items x two samples).
fbst_vaping_designs <- function(dat = fbst_vaping_data()) {
  designs <- list()
  for (i in seq_len(nrow(dat))) {
    for (analysis in c("complete", "observed")) {
      g <- dat[i, ]
      if (analysis == "complete" && !g$linked_available) next
      if (analysis == "observed" && !g$full_available) next
      a <- fbst_app_args(g, FBST_PRIORS$KL, 1, analysis)
      sample <- if (analysis == "complete") "linked" else "full"
      designs[[paste(g$group, sample)]] <- fbst_get_design(a$n1, a$n2, FBST_PRIORS$KL)
    }
  }
  designs
}

fbst_run_vaping <- function() {
  dat <- fbst_vaping_data()
  fbst_write_table(dat, "vaping_transition")
  priors <- FBST_PRIORS[c("KL", "informative", "conflict")]
  rows <- list()
  for (p in names(priors)) for (D in FBST_SCENARIOS$design_effects) {
    for (i in seq_len(nrow(dat))) for (analysis in c("complete", "observed")) {
      g <- dat[i, ]
      if (analysis == "complete" && !g$linked_available) next
      if (analysis == "observed" && !g$full_available) next
      message(
        "Vaping ", g$group, "; sample=", if (analysis == "complete") "linked" else "full",
        "; prior=", p, "; D=", D
      )
      rows[[length(rows) + 1L]] <- fbst_vaping_fit(g, priors[[p]], p, D, analysis)
      fbst_write_table(do.call(rbind, rows), "vaping_results")
    }
  }
  # No between-item contrasts: the items are answered by the same students, so their
  # posteriors are not independent, and the study has no control group.
  invisible(do.call(rbind, rows))
}
