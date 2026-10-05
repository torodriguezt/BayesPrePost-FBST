# Vaping application tables, read from results/vaping_results.csv (run fits.R first).
#   app_vaping_KL.csv           Table 8 (tab:vaping_ni)
#   app_vaping_informative.csv  SI tab:S_vaping_inf
#   app_vaping_conflict.csv     SI tab:S_vaping_conf
#   vaping_frequentist.csv      SI tab:vaping_freq (McNemar and z-test at 0.05)
#   vaping_design_effect.csv    Section 4: design effects with the original cutoffs
#   vaping_calibration.csv      Section 4: k*, near-optimal ranges and excess over
#                               the Bayes rule for the six KL-optimal designs

fbst_vaping_tables <- function() {
  res <- fbst_read_table("vaping_results")
  if (is.null(res)) stop("vaping_results.csv not found; run fbst_run_vaping() first.")
  dat <- fbst_read_table("vaping_transition")
  main <- res[res$D == 1, , drop = FALSE]
  cols <- c(
    "group", "sample", "n1", "n2", "x1", "x2", "delta_mean", "delta_lo", "delta_hi",
    "prob_gt", "ev", "k_star", "alpha_star_original", "beta_star_original", "reject",
    "calibration_status"
  )
  for (p in unique(main$prior_key)) {
    fbst_write_table(main[main$prior_key == p, cols, drop = FALSE], paste0("app_vaping_", p))
  }
  # Asymptotic McNemar without continuity correction on the linked pairs and
  # two-sample z-test on the margins, at level 0.05.
  base <- main[main$prior_key == "KL", , drop = FALSE]
  freq <- base[, c("group", "sample", "n1", "n2", "x1", "x2", "p_z", "p_mcn", "reject")]
  freq$reject_z <- freq$p_z <= 0.05
  freq$reject_mcnemar <- ifelse(is.na(freq$p_mcn), NA, freq$p_mcn <= 0.05)
  freq$chi2_cc_computed <- dat$chi2_cc[match(freq$group, dat$group)]
  freq$chi2_cc_published <- dat$chi2_published[match(freq$group, dat$group)]
  freq$chi2_cc_computed[freq$sample != "linked"] <- NA
  freq$chi2_cc_published[freq$sample != "linked"] <- NA
  names(freq)[names(freq) == "reject"] <- "reject_fbst_KL"
  fbst_write_table(freq, "vaping_frequentist")
  # Design-effect sensitivity with the ORIGINAL cutoffs.
  de <- res[
    res$prior_key == "KL",
    c(
      "group", "sample", "D", "likelihood_power", "ev", "k_star", "reject",
      "delta_mean", "delta_lo", "delta_hi", "prob_gt"
    )
  ]
  fbst_write_table(de[order(de$group, de$sample, de$D), ], "vaping_design_effect")
  # Calibration of the six KL-optimal designs.
  designs <- fbst_vaping_designs()
  calibration <- do.call(rbind, lapply(names(designs), function(key) {
    d <- designs[[key]]
    parts <- strsplit(key, " ", fixed = TRUE)[[1]]
    cbind(
      data.frame(item = parts[1], sample = parts[2], n1 = d$n1, n2 = d$n2),
      fbst_design_summary(d)
    )
  }))
  fbst_write_table(calibration, "vaping_calibration")
  invisible(main)
}
