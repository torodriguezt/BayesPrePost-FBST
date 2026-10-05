# SI tab:S_reconstruction: full-sample counts allowed by the published summaries.
# The counts must reproduce the published rounded percentages and corrected
# two-proportion p-values, include the linked respondents, and fit within the
# 600 pretest and 410 posttest survey totals (McCauley et al., 2023).
# calibrate_all = TRUE also computes the exact cutoff of every admissible size (slow).
# Outputs: results/reconstruction_solutions.csv, reconstruction_summary.csv.

fbst_run_reconstruction <- function(calibrate_all = FALSE) {
  stage_totals <- c(600, 410) # completed pretest and posttest surveys (McCauley et al., 2023)
  items <- c("daily_use_addiction", "addiction_definition")
  raw <- utils::read.csv(fbst_vaping_path(),
    stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA")
  )
  dat <- fbst_vaping_data()

  # Integers x with |100 x/n - pct| <= 0.05, i.e. rounding to the published percentage.
  counts_for <- function(n, pct) {
    lo <- max(0, ceiling(n * (pct - 0.05) / 100 - 1e-9))
    hi <- min(n, floor(n * (pct + 0.05) / 100 + 1e-9))
    if (hi < lo) integer() else seq.int(lo, hi)
  }
  p_cc <- function(x1, n1, x2, n2) {
    suppressWarnings(stats::prop.test(c(x1, x2), c(n1, n2), correct = TRUE)$p.value)
  }

  admissible <- function(item) {
    full <- raw[raw$sample == "full" & raw$item == item, , drop = FALSE]
    g <- dat[dat$group == item, , drop = FALSE]
    stopifnot(nrow(full) == 1L, nrow(g) == 1L, isTRUE(g$linked_available))
    p_txt <- trimws(as.character(full$p_published))
    p_pub <- suppressWarnings(as.numeric(p_txt))
    if (!is.finite(p_pub)) stop("Published p-value of ", item, " is not a number: ", p_txt)
    decimals <- if (grepl(".", p_txt, fixed = TRUE)) nchar(sub("^[^.]*\\.", "", p_txt)) else 0L
    tol_p <- 0.5 * 10^(-decimals)
    rows <- list()
    for (n1 in seq.int(g$n_cc, stage_totals[1])) {
      x1s <- counts_for(n1, full$pct_pre)
      x1s <- x1s[x1s >= g$x1_cc & n1 - x1s >= g$n_cc - g$x1_cc]
      if (!length(x1s)) next
      for (n2 in seq.int(g$n_cc, stage_totals[2])) {
        x2s <- counts_for(n2, full$pct_post)
        x2s <- x2s[x2s >= g$x2_cc & n2 - x2s >= g$n_cc - g$x2_cc]
        for (x1 in x1s) {
          for (x2 in x2s) {
            p <- p_cc(x1, n1, x2, n2)
            if (abs(p - p_pub) <= tol_p + 1e-12) {
              rows[[length(rows) + 1L]] <- data.frame(
                item = item, n1 = n1, x1 = x1, n2 = n2, x2 = x2,
                item_nonresponse = (stage_totals[1] - n1) + (stage_totals[2] - n2),
                p_cc = p, p_published = p_pub, p_uncorrected = fbst_app_z(x1, n1, x2, n2),
                pct_pre = 100 * x1 / n1, pct_post = 100 * x2 / n2, delta_hat = x2 / n2 - x1 / n1
              )
            }
          }
        }
      }
    }
    if (!length(rows)) stop("No admissible reconstruction for ", item)
    out <- do.call(rbind, rows)
    out <- out[order(out$item_nonresponse, -out$n1), ]
    used <- out$n1 == g$n1 & out$x1 == g$x1 & out$n2 == g$n2 & out$x2 == g$x2
    if (!any(used)) warning("The reconstruction used in the article is not admissible for ", item)
    out$used_in_article <- used
    out
  }

  # Every admissible (n1, n2) has n1 <= 600 and n2 <= 410; since k* decreases with the stage
  # sizes in every calibrated design, k*(600, 410) bounds every admissible cutoff from below.
  k_lower <- fbst_vaping_designs(dat)[["nicotine_delivery_form full"]]$calibration$kstar
  k_upper <- 0.255 # largest cutoff of Table 2 (n = 50); admissible sizes are far larger

  solutions <- list()
  for (item in items) {
    message("Enumerating admissible reconstructions: ", item)
    s <- admissible(item)
    message(nrow(s), " admissible solutions; computing the KL-optimal e-values and posteriors")
    g <- dat[dat$group == item, , drop = FALSE]
    res <- lapply(seq_len(nrow(s)), function(i) {
      h <- g
      h$n1 <- s$n1[i]
      h$x1 <- s$x1[i]
      h$n2 <- s$n2[i]
      h$x2 <- s$x2[i]
      r <- fbst_vaping_fit(h, FBST_PRIORS$KL, "KL", 1, "observed", calibrate = FALSE)
      if (i %% 10 == 0) message("  ", i, " / ", nrow(s))
      data.frame(
        ev = r$ev, P_increase = r$prob_gt, delta_mean = r$delta_mean,
        delta_lo = r$delta_lo, delta_hi = r$delta_hi
      )
    })
    s <- cbind(s, do.call(rbind, res))
    s$k_star <- NA_real_
    if (calibrate_all) {
      for (key in unique(paste(s$n1, s$n2))) {
        nn <- as.integer(strsplit(key, " ")[[1]])
        message("  calibrating k* at ", nn[1], "/", nn[2])
        k <- fbst_get_design(nn[1], nn[2], FBST_PRIORS$KL)$calibration$kstar
        s$k_star[paste(s$n1, s$n2) == key] <- k
      }
    }
    s$reject_at_lower_bound <- s$ev <= k_lower
    s$retain_above_upper_bound <- s$ev > k_upper
    s$reject_exact <- if (calibrate_all) s$ev <= s$k_star else NA
    solutions[[item]] <- s
    fbst_write_table(do.call(rbind, solutions), "reconstruction_solutions")
  }
  all_s <- do.call(rbind, solutions)
  fbst_write_table(all_s, "reconstruction_solutions")

  summary <- do.call(rbind, lapply(solutions, function(s) {
    data.frame(
      item = s$item[1],
      solutions = nrow(s), used_solution_admissible = any(s$used_in_article),
      n1_range = paste(range(s$n1), collapse = "-"), n2_range = paste(range(s$n2), collapse = "-"),
      delta_hat_min = min(s$delta_hat), delta_hat_max = max(s$delta_hat),
      ev_min = min(s$ev), ev_max = max(s$ev), P_min = min(s$P_increase), P_max = max(s$P_increase),
      k_lower_bound = k_lower, all_ev_below_k_lower = all(s$ev <= k_lower),
      all_ev_above_0.255 = all(s$ev > k_upper),
      k_star_min = if (calibrate_all) min(s$k_star) else NA_real_,
      k_star_max = if (calibrate_all) max(s$k_star) else NA_real_
    )
  }))
  fbst_write_table(summary, "reconstruction_summary")
  invisible(summary)
}
