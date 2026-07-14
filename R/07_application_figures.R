# TVSFP figures: posterior density surfaces and averaged error curves with
# the adaptive cutoff k* marked, for the four treatment conditions. The seed
# matches R/06_application_tvsfp.R so both report the same k*.
# Runtime: about 10 min.

library(ALA)
library(dplyr)
library(Rcpp)
library(ggplot2)
library(tidyr)
library(lattice)

sourceCpp("src/BivBetaBinom.cpp")
source("R/01_priors_config.R")
dir.create("Figures", showWarnings = FALSE)

a0 <- prior_NI["a0"]; a1 <- prior_NI["a1"]; a2 <- prior_NI["a2"]
M_POST <- 2000
SEED   <- 7

# Error curves with the k* location marked.
plot_error_curves <- function(ev_H, ev_A, title, save_path,
                              k_grid = seq(0, 1, length.out = 401)) {
  curves <- error_curves(ev_H, ev_A, k_grid, 1, 1)
  opt <- find_kstar(ev_H, ev_A, 1, 1)
  long <- pivot_longer(curves, c("alpha", "beta", "sum"),
                       names_to = "metric", values_to = "value")
  long$metric <- factor(long$metric,
                        levels = c("alpha", "beta", "sum"),
                        labels = c("alpha", "beta", "alpha+beta"))
  p <- ggplot(long, aes(k, value, color = metric, linetype = metric)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = opt$k_star, linetype = "dashed",
               color = "grey40") +
    annotate("text", x = opt$k_star, y = 0.05,
             label = sprintf("k* = %.3f", opt$k_star),
             hjust = -0.1, size = 3.5) +
    scale_color_manual(values = c("steelblue", "tomato", "black")) +
    scale_linetype_manual(values = c("solid", "solid", "dashed")) +
    labs(x = "k", y = "Averaged error", title = title,
         color = NULL, linetype = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave(save_path, p, width = 6, height = 4, dpi = 150)
  cat("  ->", save_path, sprintf("(k* = %.4f)\n", opt$k_star))
  invisible(opt)
}

# Complete pre/post pairs for one condition.
prep_group <- function(d, sb, tv) {
  wide <- d %>%
    mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
    filter(school.based == sb, tv.based == tv) %>%
    select(id, stage, binTHKS) %>%
    pivot_wider(names_from = stage, values_from = binTHKS) %>%
    filter(!is.na(pre), !is.na(post))
  list(X = c(sum(wide$pre), sum(wide$post)),
       n = c(nrow(wide), nrow(wide)))
}

groups <- list(
  yy = list(g = prep_group(tvsfp, "yes", "yes"), lab = "CC + TV"),
  yn = list(g = prep_group(tvsfp, "yes", "no"),  lab = "CC, no TV"),
  ny = list(g = prep_group(tvsfp, "no",  "yes"), lab = "no CC, TV"),
  nn = list(g = prep_group(tvsfp, "no",  "no"),  lab = "no CC, no TV")
)

cat("\nPosterior-based error curves (TVSFP)\n")
for (k in names(groups)) {
  g <- groups[[k]]$g; lab <- groups[[k]]$lab
  n1 <- g$n[1]; n2 <- g$n[2]; x1 <- g$X[1]; x2 <- g$X[2]
  cat(sprintf("%s: n1=%d n2=%d x1=%d x2=%d (M = %d)\n",
              lab, n1, n2, x1, x2, M_POST))
  set.seed(SEED)
  ev_H <- simulate_evs_H_post(n1, n2, x1, x2, a0, a1, a2, M_POST, 401)
  ev_A <- simulate_evs_A_post(n1, n2, x1, x2, a0, a1, a2, M_POST, 401)
  plot_error_curves(
    ev_H, ev_A,
    title = sprintf("Posterior-based: %s (n1=%d, n2=%d)", lab, n1, n2),
    save_path = sprintf("Figures/error_curves_post_%s.png", k))
}

cat("\nBivariate posterior surfaces (TVSFP)\n")
for (k in names(groups)) {
  g <- groups[[k]]$g; lab <- groups[[k]]$lab
  n1 <- g$n[1]; n2 <- g$n[2]; x1 <- g$X[1]; x2 <- g$X[2]
  consts <- bb_constants(n1, n2, x1, x2, a0, a1, a2)
  xs <- seq(0.01, 0.99, length.out = 80)
  ys <- seq(0.01, 0.99, length.out = 80)
  z  <- densBB_grid(xs, ys, consts)
  df <- expand.grid(theta1 = xs, theta2 = ys)
  df$density <- as.vector(z)
  p <- wireframe(
    density ~ theta1 * theta2, data = df,
    xlab = expression(theta[1]),
    ylab = expression(theta[2]),
    zlab = NULL,
    main = sprintf("%s (n1=%d, n2=%d, x1=%d, x2=%d)",
                   lab, n1, n2, x1, x2),
    scales = list(arrows = FALSE, cex = 0.7),
    drape = TRUE,
    col.regions = "gold",
    colorkey = FALSE,
    screen = list(z = -30, x = -60),
    par.settings = list(axis.line = list(col = "transparent"))
  )
  png(sprintf("Figures/posterior_surface_%s.png", k),
      width = 700, height = 600, res = 100)
  print(p)
  dev.off()
  cat(sprintf("  -> Figures/posterior_surface_%s.png\n", k))
}

cat("\nDone. Application figures written to Figures/\n")
