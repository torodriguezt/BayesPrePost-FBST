# Plot type-I error, type-II error and total error across cutoffs.

plot_error_curves <- function(data = read_vaping_counts()) {
  data <- data[order(
    match(data$sample, c("linked", "full")),
    match(data$item, names(ITEM_LABELS))
  ), ]
  designs <- lapply(seq_len(nrow(data)), function(i) {
    calibrate_test(data$n1[i], data$n2[i], PRIORS$KL)
  })
  save_figure("vaping_error_curves", function() {
    graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

    for (i in seq_len(nrow(data))) {
      curve <- error_curve(designs[[i]])
      risk <- curve$alpha + curve$beta
      graphics::plot(curve$k, risk,
        type = "n", xlim = c(0, 1), ylim = c(0, 1),
        xlab = "Cutoff k", ylab = "Prior-averaged error",
        main = paste(ITEM_LABELS[[data$item[i]]],
          SAMPLE_LABELS[[data$sample[i]]],
          sep = "\n"
        ), cex.main = .85
      )

      # Shade near-optimal intervals separately, preserving gaps.
      for (j in which(risk <= min(risk) + NUMERICS$near_risk)) {
        right <- if (j < nrow(curve)) {
          curve$k[j + 1L]
        } else {
          1
        }

        if (right >= 0) {
          graphics::rect(max(0, curve$k[j]), 0, right, 1,
            col = "grey90", border = NA
          )
        }
      }

      graphics::lines(curve$k, curve$alpha, type = "s", col = "#0072B2")
      graphics::lines(curve$k, curve$beta, type = "s", col = "#D55E00")
      graphics::lines(curve$k, risk, type = "s", lty = 2)
      graphics::abline(v = designs[[i]]$calibration$kstar, lty = 3)
      graphics::legend("top",
        legend = expression(alpha, beta, alpha + beta),
        col = c("#0072B2", "#D55E00", "black"),
        lty = c(1, 1, 2), bty = "n", cex = .75
      )
    }
  }, width = 10.8, height = 6.8)
}
