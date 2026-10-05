# Section 4 (unlinked surveys): estimator under missingness at random given the
# pretest answer, eq. (mar-vape): baseline distribution of every student seen at
# pretest combined with the transition probabilities of the linked students.
# Assumes the linked students are a subset of the pretest sample.
# Output: results/vaping_missingness_MAR.csv.

fbst_run_missing_data <- function(dat = fbst_vaping_data()) {
  d <- dat[dat$linked_available & dat$full_available & dat$cells_available, , drop = FALSE]
  if (!nrow(d)) {
    message("The MAR estimator needs full-sample margins and linked transition cells; skipped.")
    return(invisible(NULL))
  }
  out <- d
  out$loss_fraction_pre <- (d$n1 - d$n_cc) / d$n1
  out$baseline_rate_lost <- ifelse(d$n1 > d$n_cc, (d$x1 - d$x1_cc) / (d$n1 - d$n_cc), NA_real_)
  out$baseline_rate_linked <- d$x1_cc / d$n_cc
  out$post_rate_post_only <- ifelse(d$n2 > d$n_cc, (d$x2 - d$x2_cc) / (d$n2 - d$n_cc), NA_real_)
  out$post_rate_linked <- d$x2_cc / d$n_cc
  out$p_post_given_pre0 <- d$n01 / (d$n00 + d$n01)
  out$p_post_given_pre1 <- d$n11 / (d$n10 + d$n11)
  out$baseline_rate_observed <- d$x1 / d$n1
  out$post_rate_MAR <- with(out, p_post_given_pre0 * (1 - baseline_rate_observed) +
    p_post_given_pre1 * baseline_rate_observed)
  out$delta_MAR <- out$post_rate_MAR - out$baseline_rate_observed
  out$delta_observed <- d$x2 / d$n2 - d$x1 / d$n1
  out$delta_linked <- (d$x2_cc - d$x1_cc) / d$n_cc
  out$p_mcn_linked <- mapply(fbst_app_mcnemar, d$n01, d$n10)
  out$p_z_observed <- mapply(fbst_app_z, d$x1, d$n1, d$x2, d$n2)
  out$p_z_linked <- mapply(fbst_app_z, d$x1_cc, d$n_cc, d$x2_cc, d$n_cc)
  out$MAR_assumption <- "Y2 independent of posttest observation given Y1; posttest-only students are not used"
  fbst_write_table(out, "vaping_missingness_MAR")
  invisible(out)
}
