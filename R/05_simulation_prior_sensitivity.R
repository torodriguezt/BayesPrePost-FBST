# Prior sensitivity of the FBST decision on simulated data.
# For each prior category (non-informative, informative, conflict) and each
# combination of prior mean, sample size and effect size: calibrate a* from
# the prior, compute the e-value and k*, and record the decision.
# Writes the sim_decision and sim_ev heatmaps plus a summary CSV.
# Runtime: about 35 min.

library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcpp)

sourceCpp("src/BivBetaBinom.cpp")
dir.create("output",  showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

set.seed(7)
SEED  <- 7
M_cal <- 3000     # draws for the prior-based a* calibration
M_dec <- 1000     # draws for the posterior-based decision
NGP   <- 401      # quadrature grid for ev
NGS   <- 201      # grid for posterior grid-sampling (ev_A)
EPS   <- 1e-6
TARGET_ALPHA <- 0.05
A_GRID <- c(1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 7, 10, 15, 20)

theta1  <- 0.40
deltas  <- c(0.05, 0.10, 0.20)
n_grid  <- c(30, 50, 75, 100, 150)

categories <- list(
  noninf = list(label = "Non-informative", N0 = 2,
               mu0_grid = c(0.1, 0.3, 0.5, 0.7, 0.9)),
  inf    = list(label = "Informative", N0 = 50,
               mu0_grid = c(0.30, 0.40, 0.45, 0.50, 0.60)),
  conf   = list(label = "Conflict", N0 = 50,
               mu0_grid = c(0.05, 0.10, 0.15, 0.20, 0.25))
)

xss <- seq(EPS, 1 - EPS, length.out = NGS)

# e-values of datasets simulated from the posterior (grid sampling,
# stable under strong prior-data conflict).
ev_A_grid <- function(consts, n, a0, a1, a2, M) {
  Z <- densBB_grid(xss, xss, consts)
  p <- as.numeric(Z) / sum(Z)
  idx <- sample.int(length(p), M, replace = TRUE, prob = p)
  ii  <- ((idx - 1) %% NGS) + 1
  jj  <- ((idx - 1) %/% NGS) + 1
  th1 <- xss[ii]; th2 <- xss[jj]
  x1s <- rbinom(M, n, th1); x2s <- rbinom(M, n, th2)
  vapply(seq_len(M), function(m)
    ev_quad_from_data(n, n, x1s[m], x2s[m], a0, a1, a2), 0.0)
}

# Smallest weight a with alpha_phi below the target, from the prior alone.
calibrate_a_prior_based <- function(n, a0, a1, a2, seed = SEED) {
  set.seed(seed)
  evH <- simulate_evs_H(n, n, a0, a1, a2, M_cal, NGP)
  evA <- simulate_evs_A(n, n, a0, a1, a2, M_cal, NGP)
  alpha_of_a <- vapply(A_GRID, function(a)
    mean(evH <= find_kstar(evH, evA, a, 1)$k_star), 0.0)
  below <- which(alpha_of_a <= TARGET_ALPHA)
  if (!length(below)) return(max(A_GRID))
  if (below[1] == 1)  return(A_GRID[1])
  j <- below[1]
  approx(x = c(alpha_of_a[j - 1], alpha_of_a[j]),
         y = c(A_GRID[j - 1], A_GRID[j]), xout = TARGET_ALPHA)$y
}

t0 <- Sys.time()
all_rows <- list()

for (cat_key in names(categories)) {
  cat_info <- categories[[cat_key]]
  N0 <- cat_info$N0
  cat(sprintf("\n== Category: %s (N0 = %d) ==\n", cat_info$label, N0))

  # Calibrate a* once per (mu0, n) cell; it does not depend on delta.
  cal_grid <- expand.grid(mu0 = cat_info$mu0_grid, n = n_grid)
  cal_grid$a1 <- cal_grid$mu0 * N0
  cal_grid$a2 <- cal_grid$mu0 * N0
  cal_grid$a0 <- (1 - cal_grid$mu0) * N0
  cal_grid$a_star <- NA_real_
  cat("  calibrating a* (25 cells)")
  for (i in seq_len(nrow(cal_grid))) {
    r <- cal_grid[i, ]
    cal_grid$a_star[i] <- calibrate_a_prior_based(r$n, r$a0, r$a1, r$a2)
    cat(".")
  }
  cat(sprintf(" done (%.1f min)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  # Decision for each cell x delta, using the cell's calibrated a*.
  cat("  decision grid (25 cells x 3 deltas)")
  for (i in seq_len(nrow(cal_grid))) {
    r <- cal_grid[i, ]
    for (dl in deltas) {
      x1 <- round(r$n * theta1); x2 <- round(r$n * (theta1 + dl))
      consts <- bb_constants(r$n, r$n, x1, x2, r$a0, r$a1, r$a2)
      ev <- ev_quad(consts, find_sup_H(consts)$sup_H, NGP, EPS)
      set.seed(SEED)
      evH <- simulate_evs_H_post(r$n, r$n, x1, x2, r$a0, r$a1, r$a2, M_dec, NGP)
      evA <- ev_A_grid(consts, r$n, r$a0, r$a1, r$a2, M_dec)
      ks  <- find_kstar(evH, evA, r$a_star, 1)$k_star
      all_rows[[length(all_rows) + 1]] <- data.frame(
        category = cat_info$label, N0 = N0, mu0 = r$mu0, n = r$n,
        delta = dl, x1 = x1, x2 = x2, a_star = r$a_star,
        ev = ev, kstar = ks, reject = ev <= ks)
    }
    cat(".")
  }
  cat(sprintf(" done (%.1f min)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

res <- bind_rows(all_rows)
write.csv(res, "output/sim_decision_by_category.csv", row.names = FALSE)
cat(sprintf("\n-> output/sim_decision_by_category.csv (%.1f min total)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

res$delta_lab <- factor(sprintf("delta = %.2f", res$delta),
                        levels = sprintf("delta = %.2f", deltas))
res$n_f   <- factor(res$n)
res$mu0_f <- factor(res$mu0)

big <- theme_bw(base_size = 18) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.position = "top")

file_map     <- c("Non-informative" = "sim_decision_noninf",
                  "Informative"     = "sim_decision_inf",
                  "Conflict"        = "sim_decision_conf")
file_map_ev  <- c("Non-informative" = "sim_ev_noninf",
                  "Informative"     = "sim_ev_inf",
                  "Conflict"        = "sim_ev_conf")
n0_of        <- c("Non-informative" = 2,
                  "Informative"     = 50,
                  "Conflict"        = 50)

for (lab in names(file_map)) {
  d <- res[res$category == lab, ]

  p_dec <- ggplot(d, aes(n_f, mu0_f, fill = reject)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(reject, "R", "-")), size = 5.5) +
    facet_grid(delta_lab ~ .) +
    scale_fill_manual(values = c("FALSE" = "grey85", "TRUE" = "tomato"),
                      labels = c("Do not reject", "Reject H"), name = NULL) +
    labs(x = "n (sample size)", y = expression(mu[0]~"(prior mean)"),
         title = sprintf("FBST decision, %s prior (N0=%d)", lab, n0_of[lab])) +
    big
  fn_dec <- sprintf("Figures/%s.png", file_map[lab])
  ggsave(fn_dec, p_dec, width = 8, height = 9, dpi = 150)
  cat("  ->", fn_dec, "\n")

  p_ev <- ggplot(d, aes(n_f, mu0_f, fill = ev)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", ev)), size = 4.2) +
    facet_grid(delta_lab ~ .) +
    scale_fill_gradient(low = "firebrick", high = "steelblue", limits = c(0, 1)) +
    labs(x = "n (sample size)", y = expression(mu[0]~"(prior mean)"),
         fill = "ev(H;X)",
         title = sprintf("e-value, %s prior (N0=%d)", lab, n0_of[lab])) +
    big + theme(legend.position = "right")
  fn_ev <- sprintf("Figures/%s.png", file_map_ev[lab])
  ggsave(fn_ev, p_ev, width = 8, height = 9, dpi = 150)
  cat("  ->", fn_ev, "\n")
}

cat("\nRejection counts by category and effect size:\n")
summ <- res %>% group_by(category, delta) %>%
  summarise(n_reject = sum(reject), n_total = n(), .groups = "drop")
print(as.data.frame(summ))
