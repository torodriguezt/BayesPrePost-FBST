# Applications for the deterministic R reproduction. Sourcing defines functions only.
if (!exists("fbst_get_design", mode = "function")) source("R/pipeline_helpers.R")

fbst_application_priors <- function() FBST_PRIORS[c("KL", "informative", "conflict")]

fbst_pair_counts <- function(pre, post) {
  complete <- !is.na(pre) & !is.na(post)
  if (any(!is.na(pre) & !pre %in% 0:1) || any(!is.na(post) & !post %in% 0:1))
    stop("Outcomes must be binary or missing.")
  c(n1 = sum(!is.na(pre)), x1 = sum(pre, na.rm = TRUE),
    n2 = sum(!is.na(post)), x2 = sum(post, na.rm = TRUE),
    n_cc = sum(complete), x1_cc = sum(pre[complete]), x2_cc = sum(post[complete]),
    n00 = sum(pre[complete] == 0 & post[complete] == 0),
    n01 = sum(pre[complete] == 0 & post[complete] == 1),
    n10 = sum(pre[complete] == 1 & post[complete] == 0),
    n11 = sum(pre[complete] == 1 & post[complete] == 1),
    baseline_only = sum(!is.na(pre) & is.na(post)),
    final_only = sum(is.na(pre) & !is.na(post)),
    baseline_success_lost = sum(pre[!is.na(pre) & is.na(post)]))
}

fbst_application_data <- function(study = c("tvsfp", "toenail")) {
  study <- match.arg(study)
  if (!requireNamespace("ALA", quietly = TRUE))
    stop("Install ALA with install.packages('ALA'). No individual records are fabricated.")
  env <- new.env(parent = emptyenv())
  utils::data(list = study, package = "ALA", envir = env)
  dat <- env[[study]]
  stage <- if (study == "tvsfp") as.character(dat$stage) else as.character(dat$week)
  keep <- stage %in% if (study == "tvsfp") c("pre", "post") else c("0", "48")
  dat <- dat[keep, , drop = FALSE]; stage <- stage[keep]
  if (anyDuplicated(paste(dat$id, stage, sep = ":"))) stop("Duplicate subject/stage records: ", study)
  if (study == "tvsfp") {
    groups <- data.frame(group = c("CC - TV", "CC - No TV", "No CC - TV", "No CC - No TV"),
                         sb = c("yes", "yes", "no", "no"), tv = c("yes", "no", "yes", "no"))
    stable <- split(dat[c("school", "class", "school.based", "tv.based")], dat$id)
    if (!all(vapply(stable, function(x) all(vapply(x, function(y) length(unique(y)) == 1L, TRUE)), TRUE)))
      stop("TVSFP has inconsistent subject membership.")
    rows <- lapply(seq_len(nrow(groups)), function(i) {
      d <- dat[dat$school.based == groups$sb[i] & dat$tv.based == groups$tv[i], , drop = FALSE]
      ids <- unique(as.character(d$id)); pre <- d[d$stage == "pre", ]; post <- d[d$stage == "post", ]
      counts <- fbst_pair_counts(as.numeric(pre$THKS[match(ids, pre$id)] >= 3),
                                as.numeric(post$THKS[match(ids, post$id)] >= 3))
      data.frame(study = study, group = groups$group[i], as.list(counts),
                 schools = length(unique(d$school)), classrooms = length(unique(d$class)))
    })
  } else {
    rows <- lapply(c("A", "B"), function(arm) {
      d <- dat[dat$treatment == arm, , drop = FALSE]; ids <- unique(as.character(d$id))
      pre <- d[d$week == 0, ]; post <- d[d$week == 48, ]
      counts <- fbst_pair_counts(pre$onycholysis[match(ids, pre$id)], post$onycholysis[match(ids, post$id)])
      data.frame(study = study, group = arm, as.list(counts))
    })
  }
  out <- do.call(rbind, rows)
  out$psi <- with(out, ifelse(n01 * n10 > 0, n00 * n11 / (n01 * n10), NA_real_))
  out$data_source <- paste0("ALA::", study)
  out$ALA_version <- as.character(utils::packageVersion("ALA"))
  out
}

fbst_audit_application_data <- function(tvsfp = fbst_application_data("tvsfp"),
                                        toenail = fbst_application_data("toenail")) {
  # Reported counts are audit targets only. All analyses use the original records.
  expected <- rbind(
    data.frame(study = "tvsfp", group = c("CC - TV", "CC - No TV", "No CC - TV", "No CC - No TV"),
      n1 = c(383, 380, 416, 421), n2 = c(383, 380, 416, 421),
      x1 = c(118, 128, 145, 159), x2 = c(231, 240, 201, 175),
      n00 = c(113, 113, 160, 178), n01 = c(152, 139, 111, 84),
      n10 = c(39, 27, 55, 68), n11 = c(79, 101, 90, 91)),
    data.frame(study = "toenail", group = c("A", "B"), n1 = c(146, 148), n2 = c(133, 131),
      x1 = c(54, 55), x2 = c(14, 6), n00 = c(77, 79), n01 = c(5, 3),
      n10 = c(42, 46), n11 = c(9, 3)))
  actual <- rbind(tvsfp[names(expected)], toenail[names(expected)])
  out <- do.call(rbind, lapply(seq_len(nrow(expected)), function(i) {
    fields <- setdiff(names(expected), c("study", "group"))
    j <- which(actual$study == expected$study[i] & actual$group == expected$group[i])
    data.frame(study = expected$study[i], group = expected$group[i], quantity = fields,
      manuscript = as.numeric(expected[i, fields]), verified = as.numeric(actual[j, fields]))
  }))
  out$difference <- out$verified - out$manuscript; out$agrees <- out$difference == 0
  fbst_write_table(out, "application_data_audit")
  out
}

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

fbst_app_args <- function(g, prior, D = 1, analysis = "observed") {
  if (analysis == "complete") {
    list(x1 = g$x1_cc, n1 = g$n_cc, x2 = g$x2_cc, n2 = g$n_cc, prior = prior, power = 1 / D)
  } else list(x1 = g$x1, n1 = g$n1, x2 = g$x2, n2 = g$n2, prior = prior, power = 1 / D)
}

fbst_app_fit <- function(g, prior, prior_key, D = 1, analysis = "observed",
                         calibrate = fbst_profile() == "full") {
  args <- fbst_app_args(g, prior, D, analysis); t0 <- proc.time()[[3L]]
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
  data.frame(study = g$study, group = g$group, analysis = analysis,
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
    data_source = g$data_source, ALA_version = g$ALA_version)
}

fbst_app_density <- function(args, grid_size = 4001L) {
  fbst_cache_compute("application_delta_density", c(args, list(grid_size = grid_size)), function() {
    do.call(delta_density, c(args, list(grid = seq(-1, 1, length.out = grid_size))))
  })
}

fbst_app_contrast <- function(g1, g2, prior, prior_key, D = 1, analysis = "observed") {
  a <- fbst_app_args(g1, prior, D, analysis); b <- fbst_app_args(g2, prior, D, analysis)
  # Explicit order: delta(group 1) minus delta(group 2); independent group posteriors.
  fit <- fbst_cache_compute("application_contrast", list(first = a, second = b, grid_size = 4001L), function() {
    h1 <- fbst_app_density(a); h2 <- fbst_app_density(b)
    fine <- contrast(h1, h2)
    # Nested-grid refinement quantifies convolution/discretisation sensitivity.
    coarse <- contrast(h1[seq.int(1L, nrow(h1), 2L), ], h2[seq.int(1L, nrow(h2), 2L), ])
    fields <- c("mean", "lo", "hi", "P_positive")
    fine$refinement_max_change <- max(abs(unlist(fine[fields]) - unlist(coarse[fields])))
    fine$P_negative <- with(fine$distribution, sum(mass[delta < 0]) + 0.5 * sum(mass[delta == 0]))
    fine
  })
  data.frame(study = g1$study, comparison = paste(g1$group, "minus", g2$group),
    analysis = analysis, prior_key = prior_key, D = D, likelihood_power = 1 / D,
    mean = fit$mean, lo = fit$lo, hi = fit$hi,
    P_positive = fit$P_positive, P_negative = fit$P_negative,
    grid_size = 4001L, refinement_max_change = fit$refinement_max_change,
    method = "deterministic_delta_convolution", stringsAsFactors = FALSE)
}

fbst_toenail_missingness <- function(dat = fbst_application_data("toenail")) {
  out <- dat; out$loss_fraction <- out$baseline_only / out$n1
  out$baseline_rate_lost <- out$baseline_success_lost / out$baseline_only
  out$baseline_rate_complete <- out$x1_cc / out$n_cc
  out$p_post_given_pre0 <- out$n01 / (out$n00 + out$n01)
  out$p_post_given_pre1 <- out$n11 / (out$n10 + out$n11)
  out$baseline_rate_observed <- out$x1 / out$n1
  out$post_rate_MAR <- with(out, p_post_given_pre0 * (1 - baseline_rate_observed) +
                                  p_post_given_pre1 * baseline_rate_observed)
  out$delta_MAR <- out$post_rate_MAR - out$baseline_rate_observed
  out$delta_observed <- out$x2 / out$n2 - out$x1 / out$n1
  out$delta_complete <- (out$x2_cc - out$x1_cc) / out$n_cc
  out$p_mcn_complete <- mapply(fbst_app_mcnemar, out$n01, out$n10)
  out$p_z_observed <- mapply(fbst_app_z, out$x1, out$n1, out$x2, out$n2)
  out$p_z_complete <- mapply(fbst_app_z, out$x1_cc, out$n_cc, out$x2_cc, out$n_cc)
  out$frequentist_level <- 0.05
  out$MAR_assumption <- "Y2 independent of follow-up observation conditional on baseline Y1"
  fbst_write_table(out, "toenail_missingness_MAR")
  out
}

run_application_tvsfp <- function(calibrate = fbst_profile() == "full") {
  dat <- fbst_application_data("tvsfp"); fbst_write_table(dat, "tvsfp_transition")
  fbst_audit_application_data(tvsfp = dat)
  priors <- fbst_application_priors(); rows <- list(); ctr <- list()
  for (p in names(priors)) for (D in c(1, 1.5, 2, 3)) {
    for (i in seq_len(nrow(dat))) {
      message("TVSFP ", dat$group[i], "; prior=", p, "; D=", D)
      rows[[length(rows) + 1L]] <- fbst_app_fit(dat[i, ], priors[[p]], p, D, calibrate = calibrate)
      fbst_write_table(do.call(rbind, rows), "tvsfp_results")
    }
    for (i in seq_len(nrow(dat) - 1L)) {
      ctr[[length(ctr) + 1L]] <- fbst_app_contrast(dat[i, ], dat[nrow(dat), ], priors[[p]], p, D)
      fbst_write_table(do.call(rbind, ctr), "tvsfp_contrasts")
    }
  }
  result <- list(results = do.call(rbind, rows), contrasts = do.call(rbind, ctr))
  saveRDS(result, file.path(fbst_output_dir(), "tvsfp_results.rds"))
  result
}

run_application_toenail <- function(calibrate = fbst_profile() == "full") {
  dat <- fbst_application_data("toenail"); fbst_write_table(dat, "toenail_transition")
  fbst_toenail_missingness(dat); fbst_audit_application_data(toenail = dat)
  priors <- fbst_application_priors(); rows <- list(); ctr <- list()
  for (p in names(priors)) for (analysis in c("observed", "complete")) {
    for (i in seq_len(nrow(dat))) {
      message("Toenail ", dat$group[i], "; prior=", p, "; analysis=", analysis)
      rows[[length(rows) + 1L]] <- fbst_app_fit(dat[i, ], priors[[p]], p, analysis = analysis, calibrate = calibrate)
      fbst_write_table(do.call(rbind, rows), "toenail_results")
    }
    ctr[[length(ctr) + 1L]] <- fbst_app_contrast(dat[2, ], dat[1, ], priors[[p]], p, analysis = analysis)
    fbst_write_table(do.call(rbind, ctr), "toenail_contrasts")
  }
  result <- list(results = do.call(rbind, rows), contrasts = do.call(rbind, ctr))
  saveRDS(result, file.path(fbst_output_dir(), "toenail_results.rds"))
  result
}

fbst_app_read <- function(name) {
  path <- file.path(fbst_output_dir(), paste0(name, ".csv"))
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

run_application_tables <- function() {
  # CSV values remain unrounded. TeX is regenerated from these computed results.
  names <- c("tvsfp_transition", "tvsfp_results", "tvsfp_contrasts", "toenail_transition",
             "toenail_results", "toenail_contrasts", "toenail_missingness_MAR", "application_data_audit")
  status <- do.call(rbind, lapply(names, function(name) {
    dat <- fbst_app_read(name)
    if (!is.null(dat)) fbst_write_table(dat, name)
    data.frame(table = name, status = if (is.null(dat)) "not_run" else "generated_from_results")
  }))
  for (study in c("tvsfp", "toenail")) {
    results <- fbst_app_read(paste0(study, "_results")); if (is.null(results)) next
    main <- results[results$D == 1, , drop = FALSE]
    for (p in unique(main$prior_key)) {
      small <- main[main$prior_key == p,
        c("group", "analysis", "n1", "n2", "x1", "x2", "delta_mean", "delta_lo", "delta_hi",
          "prob_gt", "ev", "k_star", "alpha_star_original", "beta_star_original", "reject", "calibration_status")]
      fbst_write_table(small, paste0("app_", study, "_", p))
    }
  }
  fbst_write_table(status, "application_table_status")
  fbst_application_manuscript_changes()
  fbst_application_validation()
  invisible(status)
}

fbst_app_plot_files <- function(name, draw, width = 8, height = 6) {
  directory <- fbst_figure_dir(); dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- character()
  for (ext in c("pdf", "svg", "png")) {
    path <- file.path(directory, paste0(name, ".", ext))
    if (ext == "pdf") grDevices::pdf(path, width = width, height = height)
    else if (ext == "svg") grDevices::svg(path, width = width, height = height)
    else grDevices::png(path, width = width, height = height, units = "in", res = 180)
    tryCatch(draw(), finally = grDevices::dev.off()); paths <- c(paths, path)
  }
  paths
}

fbst_app_contour_grid <- function(g, prior, size = 321L) {
  # Unit-square midpoint quadrature for displaying HPD contours only.
  # This display grid never enters estimates, e-values or decisions.
  xs <- (seq_len(size) - 0.5) / size
  a0 <- prior[1]; a1 <- prior[2]; a2 <- prior[3]
  l1 <- (a1 + g$x1 - 1) * log(xs) + (a0 + a2 + g$n1 - g$x1 - 1) * log1p(-xs)
  l2 <- (a2 + g$x2 - 1) * log(xs) + (a0 + a1 + g$n2 - g$x2 - 1) * log1p(-xs)
  logz <- outer(l1, l2, "+") - sum(prior) * log1p(-outer(xs, xs))
  z <- exp(logz - max(logz)); zs <- sort(as.numeric(z), decreasing = TRUE)
  cum <- cumsum(zs) / sum(zs)
  levels <- vapply(c(0.95, 0.50), function(p) zs[which(cum >= p)[1L]], 0.0)
  list(x = xs, z = z, levels = levels)
}

run_application_figures <- function() {
  dat <- fbst_application_data("tvsfp"); paths <- character(); status <- list()
  cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  grids <- lapply(seq_len(nrow(dat)), function(i) fbst_app_contour_grid(dat[i, ], FBST_PRIORS$KL))
  paths <- c(paths, fbst_app_plot_files("app_posterior_contours", function() {
    graphics::par(mar = c(4.3, 4.3, 1, 1))
    graphics::plot(NA_real_, NA_real_, xlim = c(0.20, 0.72), ylim = c(0.25, 0.75),
      xlab = expression(theta[1]), ylab = expression(theta[2]), asp = 1)
    graphics::abline(a = 0, b = 1, lty = 2, col = "grey40")
    for (i in seq_len(nrow(dat))) {
      gr <- grids[[i]]
      graphics::contour(gr$x, gr$x, gr$z, levels = gr$levels, col = cols[i],
        lty = c(2, 1), lwd = c(1, 1.7), drawlabels = FALSE, add = TRUE)
      graphics::points(dat$x1[i] / dat$n1[i], dat$x2[i] / dat$n2[i], pch = 19, col = cols[i])
    }
    graphics::legend("bottomright", legend = dat$group, col = cols, lty = 1, bty = "n", cex = .82)
    graphics::legend("topleft", legend = c("50% HPD", "95% HPD"), lty = c(1, 2), bty = "n", cex = .8)
  }))
  status[[1L]] <- data.frame(figure = "app_posterior_contours", status = "generated", reason = "")
  results <- fbst_app_read("tvsfp_results")
  available <- !is.null(results) && any(results$D == 1 & results$prior_key == "KL") &&
    all(results$calibration_status[results$D == 1 & results$prior_key == "KL"] == "enumerated")
  if (available) {
    designs <- lapply(seq_len(nrow(dat)), function(i) fbst_get_design(dat$n1[i], dat$n2[i], FBST_PRIORS$KL))
    curves <- lapply(designs, function(d) {
      ev <- as.numeric(d$EV); ord <- order(ev); ev <- ev[ord]
      alpha <- cumsum(as.numeric(d$pH)[ord]); beta <- 1 - cumsum(as.numeric(d$pA)[ord])
      keep <- !duplicated(ev, fromLast = TRUE)
      data.frame(k = c(-1, ev[keep]), alpha = c(0, alpha[keep]), beta = c(1, beta[keep]))
    })
    paths <- c(paths, fbst_app_plot_files("app_error_curves", function() {
      graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
      for (i in seq_along(curves)) {
        cr <- curves[[i]]; risk <- cr$alpha + cr$beta
        graphics::plot(cr$k, risk, type = "n", xlim = c(0, 1), ylim = c(0, 1),
          xlab = "k", ylab = "Prior-averaged error", main = dat$group[i])
        near <- which(risk <= min(risk) + 0.002)
        # Shade EVERY qualifying decision interval separately, preserving gaps.
        for (j in near) {
          right <- if (j < nrow(cr)) cr$k[j + 1L] else 1
          if (right >= 0) graphics::rect(max(0, cr$k[j]), 0, right, 1,
            col = grDevices::adjustcolor("grey60", alpha.f = .25), border = NA)
        }
        graphics::lines(cr$k, cr$alpha, type = "s", col = "#0072B2")
        graphics::lines(cr$k, cr$beta, type = "s", col = "#D55E00")
        graphics::lines(cr$k, risk, type = "s", lty = 2)
        graphics::abline(v = designs[[i]]$calibration$kstar, lty = 3)
        graphics::legend("topright", legend = c("alpha", "beta", "alpha + beta"),
          col = c("#0072B2", "#D55E00", "black"), lty = c(1, 1, 2), bty = "n", cex = .7)
      }
    }))
    status[[2L]] <- data.frame(figure = "app_error_curves", status = "generated", reason = "")
  } else status[[2L]] <- data.frame(figure = "app_error_curves", status = "not_run",
    reason = "Full TVSFP enumeration required; run R/06_application_tvsfp.R with FBST_PROFILE=full")
  for (study in c("tvsfp", "toenail")) {
    ct <- fbst_app_read(paste0(study, "_contrasts")); if (is.null(ct)) next
    ct <- ct[ct$prior_key == "KL", , drop = FALSE]
    paths <- c(paths, fbst_app_plot_files(paste0(study, "_contrasts"), function() {
      y <- seq_len(nrow(ct)); graphics::par(mar = c(4, 11, 2, 1))
      labels <- if (study == "tvsfp") paste0(sub(" minus.*$", "", ct$comparison), " (D = ", ct$D, ")")
        else paste0("B minus A (", ct$analysis, ")")
      graphics::plot(ct$mean, y, xlim = range(c(ct$lo, ct$hi)), yaxt = "n", ylab = "",
        xlab = if (study == "tvsfp") expression(delta[group] - delta[control]) else expression(delta[B] - delta[A]),
        main = if (study == "tvsfp") "TVSFP: change relative to control" else "Onychomycosis: treatment contrast",
        pch = 19)
      graphics::segments(ct$lo, y, ct$hi, y); graphics::abline(v = 0, lty = 2, col = "grey50")
      graphics::axis(2, at = y, labels = labels, las = 1, cex.axis = .8)
    }, height = max(4, nrow(ct) * .30)))
    status[[length(status) + 1L]] <- data.frame(figure = paste0(study, "_contrasts"), status = "generated", reason = "")
  }
  fbst_write_table(do.call(rbind, status), "application_figure_status")
  invisible(paths)
}

fbst_application_manuscript_changes <- function() {
  # Audit values below are taken from the manuscript; they never enter inference.
  tv <- fbst_app_read("tvsfp_results"); tn <- fbst_app_read("toenail_results")
  ct <- fbst_app_read("tvsfp_contrasts"); cn <- fbst_app_read("toenail_contrasts")
  rows <- list()
  compare <- function(location, quantity, old, verified, digits, source) {
    if (!length(verified) || is.na(verified)) return(invisible(NULL))
    if (sprintf(paste0("%.", digits, "f"), old) == sprintf(paste0("%.", digits, "f"), verified))
      return(invisible(NULL))
    rows[[length(rows) + 1L]] <<- data.frame(location = location, quantity = quantity,
      old = as.character(old), new = format(verified, digits = 12),
      reason = "Deterministic R calculation changes the displayed rounding", source = source)
  }
  if (!is.null(tv)) {
    groups <- c("CC - TV", "CC - No TV", "No CC - TV", "No CC - No TV")
    values <- list(
      KL = list(delta_mean = c(.293,.292,.134,.038), delta_lo = c(.225,.224,.067,-.028),
        delta_hi = c(.360,.359,.199,.103), k_star = c(.106,.106,.099,.099),
        alpha_star_original = c(.033,.033,.030,.030), beta_star_original = c(.189,.189,.184,.183)),
      informative = list(delta_mean = c(.237,.235,.110,.031), delta_lo = c(.175,.173,.050,-.029),
        delta_hi = c(.298,.297,.171,.091), k_star = c(.436,.434,.425,.428),
        alpha_star_original = c(.155,.154,.152,.154), beta_star_original = c(.482,.485,.472,.468)),
      conflict = list(delta_mean = c(.252,.250,.117,.033), delta_lo = c(.189,.187,.055,-.028),
        delta_hi = c(.314,.313,.177,.094), k_star = c(.350,.350,.331,.331),
        alpha_star_original = c(.123,.123,.115,.115), beta_star_original = c(.418,.420,.413,.411)))
    tags <- c(KL = "tab:thks_ni", informative = "tab:thks_inf", conflict = "tab:thks_conf")
    for (p in names(values)) for (field in names(values[[p]])) for (i in seq_along(groups)) {
      r <- tv[tv$prior_key == p & tv$D == 1 & tv$group == groups[i], ]
      if (nrow(r)) compare(paste(tags[p], groups[i]), field, values[[p]][[field]][i], r[[field]], 3, "tvsfp_results.csv")
    }
    for (j in seq_along(c(1.5,2,3))) {
      D <- c(1.5,2,3)[j]; r <- tv[tv$prior_key == "KL" & tv$D == D & tv$group == "No CC - TV", ]
      if (nrow(r)) compare(paste("TVSFP design-effect paragraph; No CC - TV; D", D), "ev", c(.006,.020,.074)[j], r$ev, 3, "tvsfp_results.csv")
    }
  }
  if (!is.null(tn)) {
    values <- list(delta_mean = c(-.260,-.320,-.273,-.322), delta_lo = c(-.353,-.406,-.370,-.412),
      delta_hi = c(-.166,-.235,-.176,-.233), k_star = c(.161,.171,.162,.162),
      alpha_star_original = c(.053,.057,.054,.053), beta_star_original = c(.267,.264,.272,.274))
    r <- tn[tn$prior_key == "KL", ]; id <- match(paste(c("A","B","A","B"),c("observed","observed","complete","complete")), paste(r$group,r$analysis))
    r <- r[id, ]
    for (field in names(values)) for (i in seq_len(nrow(r)))
      compare(paste("tab:toenail_ni", r$group[i], r$analysis[i]), field, values[[field]][i], r[[field]][i], 3, "toenail_results.csv")
  }
  if (!is.null(ct)) {
    r <- ct[ct$prior_key == "KL", ]
    for (D in c(1,1.5,2,3)) {
      j <- match(D,c(1,1.5,2,3)); x <- r[r$D == D & r$comparison == "No CC - TV minus No CC - No TV", ]
      if (nrow(x)) compare(paste("tab:between; No CC - TV; D",D), "P_positive", c(.978,.950,.923,.877)[j],x$P_positive,3,"tvsfp_contrasts.csv")
    }
    base <- r[r$D == 1, ]; groups <- c("CC - TV", "CC - No TV", "No CC - TV")
    for (i in seq_along(groups)) {
      x <- base[base$comparison == paste(groups[i],"minus No CC - No TV"), ]
      if (nrow(x)) {
        compare(paste("tab:between",groups[i]), "lo",c(.161,.160,.003)[i],x$lo,3,"tvsfp_contrasts.csv")
        compare(paste("tab:between",groups[i]), "hi",c(.349,.349,.189)[i],x$hi,3,"tvsfp_contrasts.csv")
      }
    }
  }
  if (!is.null(cn)) for (analysis in c("observed","complete")) {
    i <- match(analysis,c("observed","complete")); r <- cn[cn$prior_key=="KL" & cn$analysis==analysis, ]
    if (nrow(r)) {
      compare(paste("Toenail treatment-contrast paragraph",analysis),"lo",c(-.187,-.182)[i],r$lo,3,"toenail_contrasts.csv")
      compare(paste("Toenail treatment-contrast paragraph",analysis),"hi",c(.066,.083)[i],r$hi,3,"toenail_contrasts.csv")
    }
  }
  rows[[length(rows)+1L]] <- data.frame(location="TVSFP between-condition methods paragraph", quantity="contrast method",
    old="Subtract independent posterior draws",new="Deterministic integration of delta densities and FFT convolution; nested-grid precision check",
    reason="Computational implementation now uses no Monte Carlo draws",source="tvsfp_contrasts.csv; toenail_contrasts.csv")
  out <- do.call(rbind, rows)
  fbst_write_table(out,"application_manuscript_changes")
  out
}

fbst_application_validation <- function() {
  rows <- list()
  add <- function(check, value, tolerance) {
    rows[[length(rows) + 1L]] <<- data.frame(check = check, error = value,
      tolerance = tolerance, passed = is.finite(value) && value <= tolerance)
  }
  for (study in c("tvsfp", "toenail")) {
    fits <- fbst_app_read(paste0(study, "_results"))
    ctr <- fbst_app_read(paste0(study, "_contrasts"))
    if (is.null(fits) || is.null(ctr)) next
    add(paste(study, "posterior row coverage"), abs(nrow(fits) - if (study == "tvsfp") 48L else 12L), 0)
    add(paste(study, "contrast row coverage"), abs(nrow(ctr) - if (study == "tvsfp") 36L else 6L), 0)
    add(paste(study, "convolution refinement"), max(ctr$refinement_max_change), 2e-5)
    mean_errors <- vapply(seq_len(nrow(ctr)), function(i) {
      r <- ctr[i, ]; g <- strsplit(r$comparison, " minus ", fixed = TRUE)[[1]]
      first <- fits[fits$prior_key == r$prior_key & fits$D == r$D & fits$analysis == r$analysis & fits$group == g[1], ]
      second <- fits[fits$prior_key == r$prior_key & fits$D == r$D & fits$analysis == r$analysis & fits$group == g[2], ]
      if (nrow(first) != 1 || nrow(second) != 1) return(Inf)
      abs(r$mean - (first$delta_mean - second$delta_mean))
    }, 0.0)
    add(paste(study, "convolution mean vs independent posterior means"), max(mean_errors), 2e-7)
    if (all(c("prob_gt", "prob_lt") %in% names(fits)))
      add(paste(study, "direct directional integrals sum to one"), max(abs(fits$prob_gt + fits$prob_lt - 1)), 2e-6)
  }
  if (!length(rows)) return(invisible(NULL))
  out <- do.call(rbind, rows)
  fbst_write_table(out, "application_validation")
  if(any(!out$passed)) stop("Application validation failed; see application_validation.csv")
  out
}
