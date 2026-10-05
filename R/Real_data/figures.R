# Vaping application figures, read from results/vaping_results.csv (run fits.R first).
#   figures/vaping_error_curves.*     Figure fig:vaping_errors: averaged errors under the
#                                     KL-optimal prior, linked (top) and full samples (bottom)
#   figures/vaping_delta_intervals.*  Figure fig:vaping_forest: posterior mean and 95%
#                                     interval of delta per item, sample and prior (D = 1)

fbst_vaping_figures <- function() {
  dat <- fbst_vaping_data()
  res <- fbst_read_table("vaping_results")
  if (is.null(res)) stop("vaping_results.csv not found; run fbst_run_vaping() first.")
  cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
  item_col <- setNames(cols[(seq_len(nrow(dat)) - 1L) %% length(cols) + 1L], dat$group)
  # Display labels; the data keep the item codes as identifiers.
  item_labels <- c(
    nicotine_delivery_form = "Nicotine delivery form",
    daily_use_addiction = "Daily use and addiction",
    addiction_definition = "Definition of addiction"
  )
  item_label <- function(x) ifelse(x %in% names(item_labels), item_labels[x],
    tools::toTitleCase(gsub("_", " ", x))
  )
  sample_labels <- c(linked = "linked sample", full = "full sample")
  prior_labels <- c(KL = "KL-optimal prior", informative = "informative prior", conflict = "conflict prior")

  # Averaged error curves of the six KL-optimal designs.
  designs <- fbst_vaping_designs(dat)
  panels <- list()
  for (s in names(sample_labels)) for (g in dat$group) {
    d <- designs[[paste(g, s)]]
    if (is.null(d)) next
    panels[[length(panels) + 1L]] <- list(
      title = paste0(item_label(g), " (", sample_labels[[s]], ")"), design = d
    )
  }
  fbst_save_figure("vaping_error_curves", function() {
    graphics::par(mfrow = c(ceiling(length(panels) / 3), min(3, length(panels))), mar = c(4, 4, 2, 1))
    for (pn in panels) {
      cr <- fbst_error_curve(pn$design)
      risk <- cr$alpha + cr$beta
      graphics::plot(cr$k, risk,
        type = "n", xlim = c(0, 1), ylim = c(0, 1),
        xlab = "k", ylab = "Prior-averaged error", main = pn$title, cex.main = .85
      )
      # Shade every near-optimal decision interval separately, preserving gaps.
      for (j in which(risk <= min(risk) + 0.002)) {
        right <- if (j < nrow(cr)) cr$k[j + 1L] else 1
        if (right >= 0) graphics::rect(max(0, cr$k[j]), 0, right, 1,
          col = grDevices::adjustcolor("grey60", alpha.f = .25), border = NA
        )
      }
      graphics::lines(cr$k, cr$alpha, type = "s", col = "#0072B2")
      graphics::lines(cr$k, cr$beta, type = "s", col = "#D55E00")
      graphics::lines(cr$k, risk, type = "s", lty = 2)
      graphics::abline(v = pn$design$calibration$kstar, lty = 3)
      graphics::legend("top",
        legend = expression(alpha, beta, alpha + beta),
        col = c("#0072B2", "#D55E00", "black"), lty = c(1, 1, 2), bty = "n", cex = .7
      )
    }
  }, width = 3.6 * min(3, length(panels)), height = 3.4 * ceiling(length(panels) / 3))

  # Posterior of delta per item, sample and prior.
  ct <- res[res$D == 1, , drop = FALSE]
  ct <- ct[order(ct$group, ct$sample, ct$prior_key), , drop = FALSE]
  fbst_save_figure("vaping_delta_intervals", function() {
    y <- rev(seq_len(nrow(ct)))
    graphics::par(mar = c(4, 20, 2, 1))
    graphics::plot(ct$delta_mean, y,
      xlim = range(c(0, ct$delta_lo, ct$delta_hi)), yaxt = "n",
      ylab = "", xlab = expression(delta == theta[2] - theta[1]), pch = 19,
      col = item_col[ct$group], main = "Posterior mean and 95% credible interval (D = 1)", cex.main = .9
    )
    graphics::segments(ct$delta_lo, y, ct$delta_hi, y, col = item_col[ct$group])
    graphics::abline(v = 0, lty = 2, col = "grey50")
    graphics::axis(2,
      at = y, labels = paste0(
        item_label(ct$group), " (", sample_labels[ct$sample], ", ",
        ifelse(ct$prior_key %in% names(prior_labels), prior_labels[ct$prior_key], ct$prior_key), ")"
      ),
      las = 1, cex.axis = .7
    )
  }, height = max(4, nrow(ct) * .32))
  invisible(NULL)
}
