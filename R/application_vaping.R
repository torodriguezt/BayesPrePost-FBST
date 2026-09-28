# Vaping-knowledge application (McCauley, Baiocchi, Cruse and Halpern-Felsher, 2023).
# Same analysis as the TVSFP and onychomycosis applications: FBST at the adaptive
# cutoff, posterior of delta, P(theta2 > theta1 | X), frequentist comparators,
# design-effect sensitivity, and the incomplete-follow-up analysis.
# Sourcing defines functions only.
#
# Input: a CSV of published summary counts, one row per (sample, item).
#   sample = "linked": students linked pre and post. Uses n1 = n2, x1, x2 and, when
#            available, the transition cells n00, n01, n10, n11.
#   sample = "full":   every student observed at each stage. Uses n1, n2, x1, x2
#            (if x1/x2 are empty they are rounded from pct_pre/pct_post).
#            Rows with empty margins are skipped as pending.
# Items: nicotine_delivery_form, daily_use_addiction, addiction_definition.
# The path defaults to vaping_mccauley2023.csv; override it with FBST_VAPING_DATA.
# Priors default to KL, informative and conflict; restrict them with
# FBST_VAPING_PRIORS, e.g. FBST_VAPING_PRIORS=KL.
if (!exists("fbst_app_fit", mode = "function")) source("R/manuscript_applications.R")

fbst_vaping_path <- function() Sys.getenv("FBST_VAPING_DATA", "vaping_mccauley2023.csv")

fbst_vaping_priors <- function() {
  keys <- trimws(strsplit(Sys.getenv("FBST_VAPING_PRIORS", "KL,informative,conflict"), ",")[[1]])
  available <- fbst_application_priors()
  bad <- setdiff(keys, names(available))
  if (length(bad)) stop("Unknown prior in FBST_VAPING_PRIORS: ", paste(bad, collapse = ", "))
  available[keys]
}

fbst_vaping_design_effects <- function() c(1, 1.5, 2, 3)

fbst_vaping_data <- function(path = fbst_vaping_path()) {
  if (!file.exists(path)) stop("Vaping data not found: ", path, " (set FBST_VAPING_DATA).")
  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                         na.strings = c("", "NA"))
  miss <- setdiff(c("sample", "item", "x1", "x2"), names(raw))
  if (length(miss)) stop("Vaping data lacks columns: ", paste(miss, collapse = ", "))
  # Sample sizes: n1/n2 per row (current layout) or N / n_pre / n_post (earlier layout).
  if (!"n1" %in% names(raw)) raw$n1 <- if ("n_pre" %in% names(raw)) raw$n_pre else NA
  if (!"n2" %in% names(raw)) raw$n2 <- if ("n_post" %in% names(raw)) raw$n_post else NA
  if ("N" %in% names(raw)) {
    raw$n1 <- ifelse(is.na(raw$n1), raw$N, raw$n1); raw$n2 <- ifelse(is.na(raw$n2), raw$N, raw$n2)
  }
  # Published McNemar statistic: its own column, or "chi2=..." inside test_published.
  if (!"chi2_mcnemar_published" %in% names(raw) && "test_published" %in% names(raw)) {
    hit <- regmatches(raw$test_published, regexpr("chi2 *= *[0-9.]+", raw$test_published))
    raw$chi2_mcnemar_published <- NA
    raw$chi2_mcnemar_published[grepl("chi2 *= *[0-9.]+", raw$test_published)] <- sub("chi2 *= *", "", hit)
  }
  optional <- c("description", "pct_pre", "pct_post", "n00", "n01", "n10", "n11",
                "psi", "chi2_mcnemar_published")
  for (v in setdiff(optional, names(raw))) raw[[v]] <- NA
  raw$description <- ifelse(is.na(raw$description), raw$item, raw$description)
  num <- function(x) suppressWarnings(as.numeric(x))
  one <- function(d, v) if (nrow(d)) num(d[[v]][1L]) else NA_real_
  rows <- lapply(unique(raw$item), function(it) {
    L <- raw[raw$sample == "linked" & raw$item == it, , drop = FALSE]
    C <- raw[raw$sample == "full" & raw$item == it, , drop = FALSE]
    if (nrow(L) > 1L || nrow(C) > 1L) stop("Duplicated item/sample rows for item ", it)
    # Linked (complete-pairs) sample.
    n_cc <- one(L, "n1"); x1_cc <- one(L, "x1"); x2_cc <- one(L, "x2")
    if (nrow(L) && is.finite(one(L, "n2")) && one(L, "n2") != n_cc)
      stop("Item ", it, ": the linked sample must have n1 = n2.")
    linked <- all(is.finite(c(n_cc, x1_cc, x2_cc)))
    cells <- c(n00 = one(L, "n00"), n01 = one(L, "n01"), n10 = one(L, "n10"), n11 = one(L, "n11"))
    has_cells <- linked && all(is.finite(cells))
    if (has_cells && (sum(cells) != n_cc || cells[["n10"]] + cells[["n11"]] != x1_cc ||
                      cells[["n01"]] + cells[["n11"]] != x2_cc))
      stop("Transition cells of item ", it, " do not reproduce N, x1 and x2.")
    # Full sample: margins as observed, n_pre != n_post allowed.
    n1 <- one(C, "n1"); n2 <- one(C, "n2"); x1 <- one(C, "x1"); x2 <- one(C, "x2")
    source_x <- "counts"
    if (all(is.finite(c(n1, n2))) && !all(is.finite(c(x1, x2)))) {
      p1 <- one(C, "pct_pre"); p2 <- one(C, "pct_post")
      if (all(is.finite(c(p1, p2)))) {
        x1 <- round(p1 * n1 / 100); x2 <- round(p2 * n2 / 100); source_x <- "rounded_from_percent"
      }
    }
    full <- all(is.finite(c(n1, n2, x1, x2)))
    # The linked students are a subset of those seen at each stage.
    if (full && linked && (n_cc > n1 || n_cc > n2 || x1_cc > x1 || x2_cc > x2 ||
                           n_cc - x1_cc > n1 - x1 || n_cc - x2_cc > n2 - x2))
      stop("Item ", it, ": the linked sample is not contained in the full sample.")
    if (!linked && !full) {
      message("Vaping item ", it, ": no complete margins in either sample; skipped.")
      return(NULL)
    }
    if (!full) { n1 <- n2 <- x1 <- x2 <- NA_real_; source_x <- NA_character_ }
    if (!has_cells) cells[] <- NA_real_
    data.frame(study = "vaping", group = it,
      description = if (nrow(L)) L$description[1L] else C$description[1L],
      n1 = n1, x1 = x1, n2 = n2, x2 = x2, full_available = full, full_margins_source = source_x,
      n_cc = n_cc, x1_cc = x1_cc, x2_cc = x2_cc, linked_available = linked,
      n00 = cells[["n00"]], n01 = cells[["n01"]], n10 = cells[["n10"]], n11 = cells[["n11"]],
      cells_available = has_cells,
      baseline_only = if (full && linked) n1 - n_cc else NA_real_,
      final_only = if (full && linked) n2 - n_cc else NA_real_,
      baseline_success_lost = if (full && linked) x1 - x1_cc else NA_real_,
      psi = if (has_cells && cells[["n01"]] * cells[["n10"]] > 0)
        cells[["n00"]] * cells[["n11"]] / (cells[["n01"]] * cells[["n10"]]) else NA_real_,
      psi_published = one(L, "psi"),
      chi2_cc = if (has_cells && cells[["n01"]] + cells[["n10"]] > 0)
        (abs(cells[["n01"]] - cells[["n10"]]) - 1)^2 / (cells[["n01"]] + cells[["n10"]]) else NA_real_,
      chi2_published = one(L, "chi2_mcnemar_published"),
      pct_pre_linked = one(L, "pct_pre"), pct_post_linked = one(L, "pct_post"),
      pct_pre_full = one(C, "pct_pre"), pct_post_full = one(C, "pct_post"),
      data_source = basename(path), ALA_version = NA_character_, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out)) stop("No analysable item in ", path)
  out
}

fbst_audit_vaping_data <- function(dat = fbst_vaping_data()) {
  # Published summaries are audit targets only; inference uses the counts.
  fbst_write_table(dat, "vaping_transition")
  add <- function(item, sample, check, computed, reported, tolerance)
    if (is.finite(computed) && is.finite(reported))
      data.frame(item = item, sample = sample, check = check, computed = computed,
                 reported = reported, difference = computed - reported,
                 agrees = abs(computed - reported) <= tolerance)
  rows <- list()
  for (i in seq_len(nrow(dat))) {
    d <- dat[i, ]
    rows <- c(rows, list(
      add(d$group, "linked", "pct_pre", 100 * d$x1_cc / d$n_cc, d$pct_pre_linked, 0.05),
      add(d$group, "linked", "pct_post", 100 * d$x2_cc / d$n_cc, d$pct_post_linked, 0.05),
      add(d$group, "linked", "psi", d$psi, d$psi_published, 0.006),
      # The published McNemar statistic includes the continuity correction.
      add(d$group, "linked", "chi2_mcnemar_continuity_corrected", d$chi2_cc, d$chi2_published, 0.006),
      add(d$group, "full", "pct_pre", 100 * d$x1 / d$n1, d$pct_pre_full, 0.05),
      add(d$group, "full", "pct_post", 100 * d$x2 / d$n2, d$pct_post_full, 0.05)))
  }
  rows <- Filter(Negate(is.null), rows)
  out <- if (length(rows)) do.call(rbind, rows) else
    data.frame(item = character(), sample = character(), check = character())
  fbst_write_table(out, "vaping_data_audit")
  fbst_vaping_missingness(dat)
  if (nrow(out) && any(!out$agrees)) warning("Vaping data audit: some published summaries disagree; see vaping_data_audit.csv")
  invisible(out)
}

fbst_vaping_fit <- function(g, prior, prior_key, D, analysis, calibrate) {
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

run_application_vaping <- function(calibrate = fbst_profile() == "full") {
  dat <- fbst_vaping_data(); fbst_write_table(dat, "vaping_transition")
  priors <- fbst_vaping_priors(); rows <- list()
  for (p in names(priors)) for (D in fbst_vaping_design_effects()) {
    for (i in seq_len(nrow(dat))) for (analysis in c("complete", "observed")) {
      g <- dat[i, ]
      if (analysis == "complete" && !g$linked_available) next
      if (analysis == "observed" && !g$full_available) next
      message("Vaping ", g$group, "; sample=", if (analysis == "complete") "linked" else "full",
              "; prior=", p, "; D=", D)
      rows[[length(rows) + 1L]] <- fbst_vaping_fit(g, priors[[p]], p, D, analysis, calibrate)
      fbst_write_table(do.call(rbind, rows), "vaping_results")
    }
  }
  # No between-item contrasts: the items are answered by the same students, so their
  # posteriors are not independent, and the study has no control group.
  result <- do.call(rbind, rows)
  saveRDS(result, file.path(fbst_output_dir(), "vaping_results.rds"))
  result
}

fbst_vaping_missingness <- function(dat = fbst_vaping_data()) {
  # Analogue of the onychomycosis MAR estimator: baseline distribution of every
  # student seen at pretest combined with the transition probabilities of the
  # linked students. Assumes the linked students are a subset of the pretest sample.
  d <- dat[dat$linked_available & dat$full_available & dat$cells_available, , drop = FALSE]
  if (!nrow(d)) {
    message("Vaping missingness analysis needs full-sample margins and linked transition cells; skipped.")
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
  out$frequentist_level <- 0.05
  out$MAR_assumption <- "Y2 independent of posttest observation given Y1; posttest-only students are not used"
  fbst_write_table(out, "vaping_missingness_MAR")
  invisible(out)
}

run_vaping_tables <- function() {
  res <- fbst_app_read("vaping_results")
  if (is.null(res)) stop("vaping_results.csv not found; run the vaping results task first.")
  dat <- fbst_app_read("vaping_transition")
  main <- res[res$D == 1, , drop = FALSE]
  cols <- c("group", "sample", "n1", "n2", "x1", "x2", "delta_mean", "delta_lo", "delta_hi",
            "prob_gt", "ev", "k_star", "alpha_star_original", "beta_star_original", "reject",
            "calibration_status")
  # Analogues of tab:thks_ni / tab:thks_inf / tab:thks_conf.
  for (p in unique(main$prior_key))
    fbst_write_table(main[main$prior_key == p, cols, drop = FALSE], paste0("app_vaping_", p))
  # Analogue of tab:thks_mcnemar: asymptotic McNemar without continuity correction on
  # the linked pairs and two-sample z-test on the margins, at level 0.05.
  base <- main[main$prior_key == unique(main$prior_key)[1L], , drop = FALSE]
  freq <- base[, c("group", "sample", "n1", "n2", "x1", "x2", "p_z", "p_mcn", "reject")]
  freq$reject_z <- freq$p_z <= 0.05
  freq$reject_mcnemar <- ifelse(is.na(freq$p_mcn), NA, freq$p_mcn <= 0.05)
  freq$chi2_cc_computed <- dat$chi2_cc[match(freq$group, dat$group)]
  freq$chi2_cc_published <- dat$chi2_published[match(freq$group, dat$group)]
  freq$chi2_cc_computed[freq$sample != "linked"] <- NA
  freq$chi2_cc_published[freq$sample != "linked"] <- NA
  names(freq)[names(freq) == "reject"] <- paste0("reject_fbst_", unique(main$prior_key)[1L])
  fbst_write_table(freq, "vaping_frequentist")
  # Design-effect sensitivity with the ORIGINAL cutoffs (TVSFP design-effect paragraph).
  de <- res[res$prior_key == unique(res$prior_key)[1L],
            c("group", "sample", "D", "likelihood_power", "ev", "k_star", "reject",
              "delta_mean", "delta_lo", "delta_hi", "prob_gt")]
  fbst_write_table(de[order(de$group, de$sample, de$D), ], "vaping_design_effect")
  fbst_vaping_validation()
  invisible(main)
}

fbst_vaping_validation <- function() {
  fits <- fbst_app_read("vaping_results"); dat <- fbst_app_read("vaping_transition")
  if (is.null(fits) || is.null(dat)) return(invisible(NULL))
  expected <- length(unique(fits$prior_key)) * length(unique(fits$D)) *
    (sum(dat$linked_available) + sum(dat$full_available))
  base <- fits[fits$D == 1, c("group", "sample", "prior_key", "k_star")]
  other <- fits[fits$D != 1, c("group", "sample", "prior_key", "k_star")]
  key <- function(x) paste(x$group, x$sample, x$prior_key)
  kdiff <- if (nrow(other) && any(is.finite(other$k_star)))
    max(abs(other$k_star - base$k_star[match(key(other), key(base))]), na.rm = TRUE) else 0
  out <- data.frame(
    check = c("vaping posterior row coverage", "vaping directional integrals sum to one",
              "vaping design effects keep the D=1 cutoff"),
    error = c(abs(nrow(fits) - expected), max(abs(fits$prob_gt + fits$prob_lt - 1)), kdiff),
    tolerance = c(0, 2e-6, 0))
  out$passed <- is.finite(out$error) & out$error <= out$tolerance
  fbst_write_table(out, "vaping_validation")
  if (any(!out$passed)) stop("Vaping validation failed; see vaping_validation.csv")
  out
}

fbst_vaping_error_curve <- function(design) {
  # Same construction as the TVSFP error-curve figure.
  ev <- as.numeric(design$EV); ord <- order(ev); ev <- ev[ord]
  alpha <- cumsum(as.numeric(design$pH)[ord]); beta <- 1 - cumsum(as.numeric(design$pA)[ord])
  keep <- !duplicated(ev, fromLast = TRUE)
  data.frame(k = c(-1, ev[keep]), alpha = c(0, alpha[keep]), beta = c(1, beta[keep]))
}

run_vaping_figures <- function() {
  dat <- fbst_vaping_data(); prior <- FBST_PRIORS$KL
  cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  item_col <- setNames(cols[(seq_len(nrow(dat)) - 1L) %% length(cols) + 1L], dat$group)
  samples <- list(linked = function(g) if (g$linked_available)
                    list(x1 = g$x1_cc, n1 = g$n_cc, x2 = g$x2_cc, n2 = g$n_cc),
                  full = function(g) if (g$full_available)
                    list(x1 = g$x1, n1 = g$n1, x2 = g$x2, n2 = g$n2))
  samples <- Filter(function(f) any(vapply(seq_len(nrow(dat)), function(i) !is.null(f(dat[i, ])), TRUE)), samples)
  # Display labels; the data keep the item codes as identifiers.
  item_labels <- c(nicotine_delivery_form = "Nicotine delivery form",
                   daily_use_addiction = "Daily use and addiction",
                   addiction_definition = "Definition of addiction")
  item_label <- function(x) ifelse(x %in% names(item_labels), item_labels[x],
    tools::toTitleCase(gsub("_", " ", x)))
  sample_labels <- c(linked = "linked sample", full = "full sample")
  prior_labels <- c(KL = "KL-optimal prior", informative = "informative prior", conflict = "conflicting prior")
  status <- list(); paths <- character()
  # 1. Joint posteriors under the KL-optimal prior (analogue of app_posterior_contours).
  paths <- c(paths, fbst_app_plot_files("vaping_posterior_contours", function() {
    graphics::par(mfrow = c(1, length(samples)), mar = c(4.3, 4.3, 2, 1))
    for (s in names(samples)) {
      graphics::plot(NA_real_, NA_real_, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
        xlab = expression(theta[1]), ylab = expression(theta[2]), main = sub("^(.)", "\\U\\1", sample_labels[[s]], perl = TRUE))
      graphics::abline(a = 0, b = 1, lty = 2, col = "grey40")
      shown <- character()
      for (i in seq_len(nrow(dat))) {
        m <- samples[[s]](dat[i, ]); if (is.null(m)) next
        gr <- fbst_app_contour_grid(m, prior, size = 801L)
        graphics::contour(gr$x, gr$x, gr$z, levels = gr$levels, col = item_col[dat$group[i]],
          lty = c(2, 1), lwd = c(1, 1.7), drawlabels = FALSE, add = TRUE)
        graphics::points(m$x1 / m$n1, m$x2 / m$n2, pch = 19, cex = .6, col = item_col[dat$group[i]])
        shown <- c(shown, dat$group[i])
      }
      graphics::legend("bottomright", legend = item_label(shown), col = item_col[shown], lty = 1, lwd = 2, bty = "n", cex = .75)
      graphics::legend("topleft", legend = c("50% HPD", "95% HPD"), lty = c(1, 2), bty = "n", cex = .75)
    }
  }, width = 5 * length(samples), height = 5.2))
  status[[length(status) + 1L]] <- data.frame(figure = "vaping_posterior_contours", status = "generated", reason = "")
  # 2. Averaged error curves under the KL prior (analogue of app_error_curves).
  res <- fbst_app_read("vaping_results")
  kl <- if (is.null(res)) NULL else res[res$D == 1 & res$prior_key == "KL", , drop = FALSE]
  if (!is.null(kl) && nrow(kl) && all(kl$calibration_status == "enumerated")) {
    panels <- list()
    for (s in names(samples)) for (i in seq_len(nrow(dat))) {
      m <- samples[[s]](dat[i, ]); if (is.null(m)) next
      panels[[length(panels) + 1L]] <- list(title = paste0(item_label(dat$group[i]), " (", sample_labels[[s]], ")"),
        design = fbst_get_design(m$n1, m$n2, prior))
    }
    paths <- c(paths, fbst_app_plot_files("vaping_error_curves", function() {
      graphics::par(mfrow = c(ceiling(length(panels) / 3), min(3, length(panels))), mar = c(4, 4, 2, 1))
      for (pn in panels) {
        cr <- fbst_vaping_error_curve(pn$design); risk <- cr$alpha + cr$beta
        graphics::plot(cr$k, risk, type = "n", xlim = c(0, 1), ylim = c(0, 1),
          xlab = "k", ylab = "Prior-averaged error", main = pn$title, cex.main = .85)
        # Shade every near-optimal decision interval separately, preserving gaps.
        for (j in which(risk <= min(risk) + 0.002)) {
          right <- if (j < nrow(cr)) cr$k[j + 1L] else 1
          if (right >= 0) graphics::rect(max(0, cr$k[j]), 0, right, 1,
            col = grDevices::adjustcolor("grey60", alpha.f = .25), border = NA)
        }
        graphics::lines(cr$k, cr$alpha, type = "s", col = "#0072B2")
        graphics::lines(cr$k, cr$beta, type = "s", col = "#D55E00")
        graphics::lines(cr$k, risk, type = "s", lty = 2)
        graphics::abline(v = pn$design$calibration$kstar, lty = 3)
        graphics::legend("top", legend = expression(alpha, beta, alpha + beta),
          col = c("#0072B2", "#D55E00", "black"), lty = c(1, 1, 2), bty = "n", cex = .7)
      }
    }, width = 3.6 * min(3, length(panels)), height = 3.4 * ceiling(length(panels) / 3)))
    status[[length(status) + 1L]] <- data.frame(figure = "vaping_error_curves", status = "generated", reason = "")
  } else status[[length(status) + 1L]] <- data.frame(figure = "vaping_error_curves", status = "not_run",
    reason = "Requires the KL prior and the full-profile calibration (FBST_PROFILE=full)")
  # 3. Posterior of delta per item, sample and prior (takes the place of the TVSFP
  #    contrast figure, which has no counterpart here).
  if (!is.null(res)) {
    ct <- res[res$D == 1, , drop = FALSE]
    ct <- ct[order(ct$group, ct$sample, ct$prior_key), , drop = FALSE]
    paths <- c(paths, fbst_app_plot_files("vaping_delta_intervals", function() {
      y <- rev(seq_len(nrow(ct))); graphics::par(mar = c(4, 20, 2, 1))
      graphics::plot(ct$delta_mean, y, xlim = range(c(0, ct$delta_lo, ct$delta_hi)), yaxt = "n",
        ylab = "", xlab = expression(delta == theta[2] - theta[1]), pch = 19,
        col = item_col[ct$group], main = "Posterior mean and 95% credible interval (D = 1)", cex.main = .9)
      graphics::segments(ct$delta_lo, y, ct$delta_hi, y, col = item_col[ct$group])
      graphics::abline(v = 0, lty = 2, col = "grey50")
      graphics::axis(2, at = y, labels = paste0(item_label(ct$group), " (", sample_labels[ct$sample], ", ",
          ifelse(ct$prior_key %in% names(prior_labels), prior_labels[ct$prior_key], ct$prior_key), ")"),
        las = 1, cex.axis = .7)
    }, height = max(4, nrow(ct) * .32)))
    status[[length(status) + 1L]] <- data.frame(figure = "vaping_delta_intervals", status = "generated", reason = "")
  }
  fbst_write_table(do.call(rbind, status), "vaping_figure_status")
  invisible(paths)
}
