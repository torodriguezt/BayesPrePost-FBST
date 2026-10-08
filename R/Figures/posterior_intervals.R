# Plot posterior means and equal-tailed credible intervals for the change.

plot_posterior_intervals <- function() {
  data <- read_vaping_counts()
  fits <- lapply(names(PRIOR_LABELS), function(name) {
    rows <- lapply(seq_len(nrow(data)), function(i) {
      fit <- analyze_margins(data[i, ], PRIORS[[name]], calibrate = FALSE)
      fit$prior <- name
      fit
    })
    do.call(rbind, rows)
  })
  fits <- do.call(rbind, fits)
  save_figure(
    "vaping_delta_intervals",
    function() {
      graphics::par(mfrow = c(1, 3), mar = c(4, 6, 3, 1), oma = c(2, 0, 0, 0))

      for (item in names(ITEM_LABELS)) {
        panel <- fits[fits$item == item, ]
        y <- match(panel$sample, c("full", "linked")) +
          c(
            KL = .16,
            informative = 0,
            conflict = -.16
          )[panel$prior]
        graphics::plot(panel$delta_mean, y,
          xlim = range(c(0, panel$delta_lo, panel$delta_hi)), ylim = c(.6, 2.4),
          yaxt = "n", ylab = "", xlab = expression(delta == theta[2] - theta[1]),
          main = ITEM_LABELS[[item]], cex.main = .9,
          pch = 19, col = PRIOR_COLORS[panel$prior]
        )
        graphics::abline(v = 0, lty = 2, col = "grey60")
        graphics::segments(
          panel$delta_lo,
          y,
          panel$delta_hi,
          y,
          col = PRIOR_COLORS[panel$prior],
          lwd = 2
        )
        graphics::axis(
          2,
          at = c(1, 2),
          labels = c("All margins", "Linked"),
          las = 1,
          cex.axis = .9
        )
        graphics::legend("bottom",
          legend = PRIOR_LABELS, col = PRIOR_COLORS,
          pch = 19, bty = "n", cex = .75
        )
      }

      graphics::mtext(
        "Posterior mean and 95% credible interval",
        side = 1,
        outer = TRUE,
        cex = .9
      )
    },
    width = 11,
    height = 4
  )
}
