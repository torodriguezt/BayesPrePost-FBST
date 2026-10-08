# Draw the sampling structure and restriction of the prior to the null.

plot_sampling_diagram <- function() {
  save_figure("diagrama_muestras", function() {
    graphics::par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))

    for (paired in c(FALSE, TRUE)) {
      graphics::plot.new()
      graphics::plot.window(xlim = c(-3, 3), ylim = c(-4.8, .7))
      graphics::text(0, 0, expression(f(theta[1], theta[2])), cex = 1.3)
      graphics::text(c(-2, 2), -2, expression(theta[1], theta[2]), cex = 1.3)
      graphics::text(c(-2, 2), -4, expression(X[1], X[2]), cex = 1.3)
      graphics::arrows(c(-.4, .4), -.35, c(-1.8, 1.8), -1.65, length = .1)
      graphics::segments(-1.6, -2, 1.6, -2)
      graphics::arrows(c(-2, 2), -2.35, c(-2, 2), -3.65, length = .1)

      if (paired) {
        graphics::arrows(-1.6, -4, 1.6, -4, code = 3, lty = 2, length = .1)
        graphics::text(0, -4.5, expression(psi), cex = 1.3)
      }

      graphics::title(main = if (paired) {
        "(b) Subjects observed at both stages"
      } else {
        "(a) Composite likelihood and prior"
      }, cex.main = .95)
    }
  }, width = 9, height = 4)
}

# A proper informative prior makes the schematic restriction integrable.
plot_prior_restriction <- function() {
  parameters <- posterior_parameters(0, 0, 0, 0, PRIORS$informative)
  log_normalizer <- posterior_normalizer(parameters)$logZ
  axis <- seq(.01, .99, length.out = 80)

  density <- function(x, y) {
    exp(log_posterior_kernel(x, y, parameters) - log_normalizer)
  }

  surface <- outer(axis, axis, density)
  section <- density(axis, axis)
  save_figure(
    "line_integral_corregido",
    function() {
      graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
      projection <- graphics::persp(axis, axis, surface,
        theta = 40, phi = 25,
        col = "#D7B34B", border = "grey60", lwd = .3, ticktype = "simple",
        xlab = "", ylab = "", zlab = "density", main = "Joint prior"
      )
      graphics::text(grDevices::trans3d(.5, -.2, 0, projection), expression(theta[1]))
      graphics::text(grDevices::trans3d(1.2, .5, 0, projection), expression(theta[2]))
      graphics::lines(
        grDevices::trans3d(axis, axis, section, projection),
        col = "#D55E00",
        lwd = 2
      )
      graphics::plot(axis, section,
        type = "l", col = "#D55E00", lwd = 2,
        xlab = expression(t), ylab = expression(f(t, t)),
        main = "Restriction to the null line"
      )
    },
    width = 9,
    height = 4
  )
}

plot_model_diagrams <- function() {
  plot_sampling_diagram()
  plot_prior_restriction()
}
