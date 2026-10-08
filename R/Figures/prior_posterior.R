# Plot prior and posterior density surfaces.

plot_prior_posterior <- function(n = 20L, grid_size = 60L) {
  scenarios <- list(
    KL = list(
      prior = PRIORS$KL,
      theta = .1,
      cap = 4
    ),
    centred = list(
      prior = PRIORS$estimation,
      theta = .5,
      cap = Inf
    ),
    conflict = list(
      prior = PRIORS$estimation,
      theta = .1,
      cap = Inf
    )
  )
  axis <- (seq_len(grid_size) - .5) / grid_size

  density <- function(parameters) {
    exp(outer(axis, axis, function(x, y) {
      log_posterior_kernel(x, y, parameters)
    }) -
      posterior_normalizer(parameters)$logZ)
  }

  save_figure(
    "est_prior_posterior",
    function() {
      graphics::par(mfcol = c(2, 3), mar = c(1, 1.5, 2.4, .5))

      for (name in names(scenarios)) {
        scenario <- scenarios[[name]]
        prior <- scenario$prior
        x <- round(n * scenario$theta)
        surfaces <- list(
          prior = pmin(density(posterior_parameters(0, 0, 0, 0, prior)), scenario$cap),
          posterior = density(posterior_parameters(x, n, x, n, prior))
        )

        for (kind in names(surfaces)) {
          projection <- graphics::persp(axis, axis, surfaces[[kind]],
            theta = 40, phi = 25, expand = .75, col = "#D7B34B",
            border = "grey30", lwd = .3, ticktype = "detailed", nticks = 4,
            cex.axis = .65, xlab = "", ylab = "", zlab = "density",
            zlim = c(0, max(surfaces[[kind]]))
          )
          graphics::text(
            grDevices::trans3d(.5, -.28, 0, projection),
            expression(theta[1])
          )
          graphics::text(
            grDevices::trans3d(1.28, .5, 0, projection),
            expression(theta[2])
          )
          title <- if (kind == "prior") {
            paste("Prior:", name)
          } else {
            bquote(list(n == .(n), x[1] == .(x), x[2] == .(x)))
          }

          graphics::title(main = title, cex.main = .95)
        }
      }
    },
    width = 9,
    height = 6
  )
}
