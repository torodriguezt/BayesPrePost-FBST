# SI tab:S_recoding: reverse the outcome coding (incorrect answer = 1) and compare the
# e-value and the probability of an increase with the original analysis.
# The KL-optimal and informative priors are unchanged by the recoding, so the
# original cutoffs apply; the conflict prior would need new calibrations.
# Output: results/recoding.csv.

fbst_recode <- function(g) {
  # Incorrect answer coded as 1: every success count becomes its complement.
  h <- g
  h$x1 <- g$n1 - g$x1
  h$x2 <- g$n2 - g$x2
  h$x1_cc <- g$n_cc - g$x1_cc
  h$x2_cc <- g$n_cc - g$x2_cc
  h$n00 <- g$n11
  h$n11 <- g$n00
  h$n01 <- g$n10
  h$n10 <- g$n01
  h
}

fbst_run_recoding <- function(priors = c("KL", "informative")) {
  dat <- fbst_vaping_data()
  rows <- list()
  for (p in priors) {
    for (i in seq_len(nrow(dat))) {
      for (analysis in c("complete", "observed")) {
        g <- dat[i, ]
        if (analysis == "complete" && !g$linked_available) next
        if (analysis == "observed" && !g$full_available) next
        sample <- if (analysis == "complete") "linked" else "full"
        prior <- FBST_PRIORS[[p]]
        # Same beliefs about the recoded probabilities: mu0 -> 1 - mu0 with N0 fixed, i.e.
        # lambda0 and lambda1 = lambda2 swap. Unchanged when mu0 = 0.5 (KL, informative).
        prior_rec <- setNames(unname(prior[c(2, 1, 1)]), c("a0", "a1", "a2"))
        message("Recoding: ", g$group, " / ", sample, " / ", p)
        r0 <- fbst_vaping_fit(g, prior, p, 1, analysis, calibrate = FALSE)
        r1 <- fbst_vaping_fit(fbst_recode(g), prior_rec, p, 1, analysis, calibrate = FALSE)
        k0 <- fbst_get_design(r0$n1, r0$n2, prior)$calibration$kstar
        k1 <- if (identical(unname(prior_rec), unname(prior))) {
          k0
        } else {
          fbst_get_design(r1$n1, r1$n2, prior_rec)$calibration$kstar
        }
        rows[[length(rows) + 1L]] <- data.frame(
          prior = p, item = g$group, sample = sample,
          n1 = r0$n1, n2 = r0$n2, x1 = r0$x1, x2 = r0$x2,
          k_star_original = k0, k_star_reversed = k1,
          ev_original = r0$ev, ev_reversed = r1$ev,
          ev_change = r1$ev - r0$ev,
          # For correct answers P(theta2 > theta1) is P(theta2' < theta1') after recoding.
          P_increase_original = r0$prob_gt, P_increase_reversed = r1$prob_lt,
          P_change = r1$prob_lt - r0$prob_gt,
          delta_original = r0$delta_mean, delta_reversed = -r1$delta_mean,
          lo_original = r0$delta_lo, hi_original = r0$delta_hi,
          lo_reversed = -r1$delta_hi, hi_reversed = -r1$delta_lo,
          reject_original = r0$ev <= k0, reject_reversed = r1$ev <= k1
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$decision_changed <- out$reject_original != out$reject_reversed
  fbst_write_table(out, "recoding")
  invisible(out)
}
