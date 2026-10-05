# Table 5 (tab:boundary): pointwise null rejection pi(theta, theta) near the
# boundaries under stage independence, for the KL-optimal, informative and conflict priors.
# Outputs: results/boundary_summary.csv, boundary_null_curves.csv,
# figures/boundary_null_rejection.*

# Null rejection on a grid of theta, refined around every local maximum.
fbst_null_grid <- function(reject, n1, n2) {
  theta <- sort(unique(c(
    FBST_SCENARIOS$boundary_theta,
    seq(0, 1, length.out = 201), 10^seq(-5, -1, length.out = 31),
    1 - 10^seq(-5, -1, length.out = 31)
  )))
  evaluate <- function(t) sum(fbst_sampling_mass(n1, n2, t, t, overlap = 0) * reject)
  rates <- vapply(theta, evaluate, 0)
  # Refine around every local grid maximum; report a grid maximum only.
  peaks <- which(rates >= c(-Inf, head(rates, -1L)) & rates >= c(tail(rates, -1L), -Inf))
  refinement <- unique(unlist(lapply(peaks, function(i) {
    seq(theta[max(1, i - 1L)], theta[min(length(theta), i + 1L)], length.out = 21)
  })))
  theta <- sort(unique(c(theta, refinement)))
  data.frame(theta = theta, rate = vapply(theta, evaluate, 0))
}

fbst_intervals_text <- function(theta, ok) {
  if (!any(ok)) {
    return("")
  }
  r <- rle(ok)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  paste(sprintf("[%.3f, %.3f]", theta[starts[r$values]], theta[ends[r$values]]), collapse = " U ")
}

fbst_boundary_figure <- function(curves) {
  # Okabe-Ito hues in the order of the article's other figures; the line type is a
  # second channel, so the priors stay distinguishable without colour.
  cols <- c(KL = "#0072B2", informative = "#D55E00", conflict = "#009E73")
  ltys <- c(KL = 1, informative = 2, conflict = 4)
  labels <- c(
    KL = "KL-optimal prior", informative = "Informative prior",
    conflict = "Conflict prior"
  )
  sizes <- unique(curves[c("n1", "n2")])
  draw <- function() {
    graphics::par(
      mfrow = c(1, nrow(sizes)), mar = c(4.2, 4.2, 2, .8),
      mgp = c(2.5, .7, 0), las = 1
    )
    for (i in seq_len(nrow(sizes))) {
      n1 <- sizes$n1[i]
      n2 <- sizes$n2[i]
      graphics::plot(NA_real_, NA_real_,
        xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
        xlab = expression("Common success probability " * theta),
        ylab = if (i == 1L) "Null rejection probability" else "",
        main = bquote(n[1] == .(n1) * "," ~ n[2] == .(n2)), cex.main = 1, font.main = 1
      )
      graphics::abline(h = seq(.2, .8, .2), v = seq(.2, .8, .2), col = "grey90", lwd = .7)
      for (p in names(cols)) {
        cv <- curves[curves$prior == p & curves$n1 == n1 & curves$n2 == n2, , drop = FALSE]
        if (!nrow(cv)) next
        cv <- cv[order(cv$theta), ]
        graphics::lines(c(0, 1), rep(cv$alpha_star[1], 2),
          col = grDevices::adjustcolor(cols[[p]], alpha.f = .6), lwd = .9, lty = 3
        )
        graphics::lines(cv$theta, cv$rate, col = cols[[p]], lty = ltys[[p]], lwd = 1.8)
      }
      graphics::box(col = "grey35")
      if (i == 1L) {
        graphics::legend("topright",
          legend = c(labels, expression("prior-averaged " * alpha^"*")),
          col = c(cols, "grey45"), lty = c(ltys, 3), lwd = c(1.8, 1.8, 1.8, .9),
          bg = "white", box.col = "grey80", cex = .8, seg.len = 2.8
        )
      }
    }
  }
  fbst_save_figure("boundary_null_rejection", draw, width = 9.5, height = 3.5)
}

fbst_run_boundaries <- function() {
  sizes <- matrix(c(20, 20, 50, 50, 20, 100), ncol = 2, byrow = TRUE)
  priors <- FBST_PRIORS[c("KL", "informative", "conflict")]
  curves <- summaries <- list()
  for (i in seq_len(nrow(sizes))) {
    for (p in names(priors)) {
      n1 <- sizes[i, 1]
      n2 <- sizes[i, 2]
      prior <- priors[[p]]
      message("Boundary: ", n1, "/", n2, " / ", p)
      d <- fbst_get_design(n1, n2, prior)
      cal <- d$calibration
      reject <- d$EV <= cal$kstar
      grid <- fbst_null_grid(reject, n1, n2)
      grid <- grid[order(grid$theta), ]
      S <- outer(0:n1, 0:n2, "+")
      singular <- prior[2] + prior[3] + S - 2 < 0 | prior[1] + n1 + n2 - S - 2 < 0
      best <- which.max(grid$rate)
      summaries[[length(summaries) + 1L]] <- data.frame(
        prior = p, n1 = n1, n2 = n2,
        k_star = cal$kstar, alpha_star = cal$alpha, beta_star = cal$beta,
        rate_at_040 = sum(fbst_sampling_mass(n1, n2, .4, .4, overlap = 0) * reject),
        max_rate = grid$rate[best], theta_at_max = grid$theta[best],
        theta_rate_ge_0999 = fbst_intervals_text(grid$theta, grid$rate >= .999),
        singular_mass_null = sum(d$pH[singular]), n_grid = nrow(grid)
      )
      curves[[length(curves) + 1L]] <- data.frame(
        prior = p, n1 = n1, n2 = n2,
        alpha_star = cal$alpha, theta = grid$theta, rate = grid$rate
      )
    }
  }
  summary <- do.call(rbind, summaries)
  curves <- do.call(rbind, curves)
  fbst_write_table(summary, "boundary_summary")
  fbst_write_table(curves, "boundary_null_curves")
  fbst_boundary_figure(curves)
  invisible(summary)
}
