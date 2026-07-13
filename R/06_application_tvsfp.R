# TVSFP application (paper Tables 2-4 and the McNemar table).
# For each treatment condition and each of the three priors:
#   1) calibrate the loss weight a* prior-based (simulate_evs_H / simulate_evs_A,
#      using only n and the prior, never the observed counts);
#   2) compute the observed e-value by 2D quadrature;
#   3) obtain the adaptive cutoff k* posterior-based (simulate_evs_*_post with
#      the observed counts) minimising a* alpha + beta, and decide ev <= k*;
#   4) compute P(theta1 <= theta2 | X) from SIR posterior draws.
# McNemar's test on the discordant pairs is reported as a frequentist benchmark.
#
# Data: tvsfp from the ALA package; groups defined by treatment condition
# (school.based x tv.based), complete pre/post pairs only, THKS >= 3 as success.
# Outputs: output/thks_fbst_{ni,inf,conf}.tex, output/thks_mcnemar.tex
# Runtime: ~15-25 min (12 prior-based calibrations + 12 posterior runs).

library(ALA)
library(dplyr)
library(tidyr)
library(Rcpp)

sourceCpp("src/BivBetaBinom.cpp")
source("R/01_priors_config.R")
dir.create("output", showWarnings = FALSE)

M_POST <- 2000    # posterior-based draws for k* and P(H_np)
M_CAL  <- 3000    # prior-based draws for the a* calibration
SEED   <- 7       # shared with R/07_application_figures.R so k* matches
TARGET <- 0.05
A_GRID <- c(1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 7, 10, 15, 20)
set.seed(SEED)

priors <- list(
  list(key = "non_informative", label = "Non-informative",
       a0 = prior_NI["a0"], a1 = prior_NI["a1"], a2 = prior_NI["a2"]),
  list(key = "informative", label = "Informative",
       a0 = prior_INF["a0"], a1 = prior_INF["a1"], a2 = prior_INF["a2"]),
  list(key = "conflict", label = "Conflict",
       a0 = prior_CONF["a0"], a1 = prior_CONF["a1"], a2 = prior_CONF["a2"])
)

# Complete pre/post pairs for one treatment condition, with the 2x2 cell
# counts needed by McNemar's test.
extract_pairs <- function(d, label, sb, tv) {
  wide <- d %>%
    mutate(binTHKS = ifelse(THKS >= 3, 1, 0)) %>%
    filter(school.based == sb, tv.based == tv) %>%
    select(id, stage, binTHKS) %>%
    pivot_wider(names_from = stage, values_from = binTHKS) %>%
    filter(!is.na(pre), !is.na(post))

  list(
    label = label,
    n     = nrow(wide),
    x1    = sum(wide$pre),
    x2    = sum(wide$post),
    n00   = sum(wide$pre == 0 & wide$post == 0),
    n01   = sum(wide$pre == 0 & wide$post == 1),
    n10   = sum(wide$pre == 1 & wide$post == 0),
    n11   = sum(wide$pre == 1 & wide$post == 1)
  )
}

groups <- list(
  extract_pairs(tvsfp, "CC + TV",      sb = "yes", tv = "yes"),
  extract_pairs(tvsfp, "CC, no TV",    sb = "yes", tv = "no"),
  extract_pairs(tvsfp, "no CC, TV",    sb = "no",  tv = "yes"),
  extract_pairs(tvsfp, "no CC, no TV", sb = "no",  tv = "no")
)

# Exact binomial version for few discordant pairs, chi-square otherwise.
mcnemar_pval <- function(n01, n10) {
  bc <- n01 + n10
  if (bc == 0L) return(1.0)
  if (bc < 25L) {
    2 * min(pbinom(n01, bc, 0.5), pbinom(n10, bc, 0.5))
  } else {
    pchisq((n01 - n10)^2 / bc, df = 1, lower.tail = FALSE)
  }
}

# Prior-based a*: smallest weight with alpha_phi <= TARGET, interpolated
# between the bracketing points of A_GRID. Depends only on n and the prior.
calibrate_a <- function(n, a0, a1, a2, seed = SEED) {
  set.seed(seed)
  evH <- simulate_evs_H(n, n, a0, a1, a2, M_CAL, 401)
  evA <- simulate_evs_A(n, n, a0, a1, a2, M_CAL, 401)
  alpha_of_a <- vapply(A_GRID, function(a)
    mean(evH <= find_kstar(evH, evA, a, 1)$k_star), 0.0)
  below <- which(alpha_of_a <= TARGET)
  if (!length(below)) return(max(A_GRID))
  if (below[1] == 1)  return(A_GRID[1])
  j <- below[1]
  approx(x = c(alpha_of_a[j - 1], alpha_of_a[j]),
         y = c(A_GRID[j - 1], A_GRID[j]), xout = TARGET)$y
}

# FBST with calibrated weight for one group under one prior. The seed is
# reset before the posterior-based pair so the k* here matches the error
# curves drawn by R/07_application_figures.R.
run_fbst <- function(g, prior) {
  a0 <- prior$a0; a1 <- prior$a1; a2 <- prior$a2

  a_star <- calibrate_a(g$n, a0, a1, a2)
  ev_obs <- ev_quad_from_data(g$n, g$n, g$x1, g$x2, a0, a1, a2)

  set.seed(SEED)
  ev_H <- simulate_evs_H_post(g$n, g$n, g$x1, g$x2, a0, a1, a2, M_POST, 401)
  ev_A <- simulate_evs_A_post(g$n, g$n, g$x1, g$x2, a0, a1, a2, M_POST, 401)
  opt  <- find_kstar(ev_H, ev_A, a_star, 1)

  samp      <- sample_posterior(M_POST, g$n, g$n, g$x1, g$x2, a0, a1, a2)
  prob_H_np <- mean(samp[, 1] <= samp[, 2])

  list(
    a_star = round(a_star, 3),
    ev_obs = round(ev_obs, 4),
    k_star = round(opt$k_star, 4),
    reject = ev_obs <= opt$k_star,
    prob_H = round(prob_H_np, 4)
  )
}

cat("\nTVSFP application: FBST (3 priors, calibrated a*) + McNemar\n")

rows <- list()
for (g in groups) {
  pval    <- mcnemar_pval(g$n01, g$n10)
  mcn_dec <- ifelse(pval < 0.05, "Reject", "No reject")

  cat(sprintf("\n-- %s (n=%d, x1=%d, x2=%d | n01=%d, n10=%d) --\n",
              g$label, g$n, g$x1, g$x2, g$n01, g$n10))
  cat(sprintf("  McNemar: p = %.4f -> %s\n", pval, mcn_dec))

  row <- list(group = g$label, n = g$n, x1 = g$x1, x2 = g$x2,
              n01 = g$n01, n10 = g$n10,
              p_mcn = round(pval, 4), mcn_dec = mcn_dec)

  for (pr in priors) {
    res <- run_fbst(g, pr)
    dec <- ifelse(res$reject, "Reject", "No reject")
    cat(sprintf("  FBST %-16s P(H_np)=%.4f | ev=%.4f a*=%.3f k*=%.4f -> %s\n",
                paste0("[", pr$label, "]:"),
                res$prob_H, res$ev_obs, res$a_star, res$k_star, dec))
    row[[paste0("prob_H_", pr$key)]] <- res$prob_H
    row[[paste0("ev_",     pr$key)]] <- res$ev_obs
    row[[paste0("as_",     pr$key)]] <- res$a_star
    row[[paste0("ks_",     pr$key)]] <- res$k_star
    row[[paste0("dec_",    pr$key)]] <- dec
  }

  rows[[length(rows) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE)
}

tab <- do.call(rbind, rows)
cat("\nSummary:\n")
print(tab, row.names = FALSE)

dec_sym <- function(dec) ifelse(dec == "Reject", "Reject $H$", "Do not reject")

# LaTeX table for one prior, including the calibrated a* column.
write_fbst_table <- function(tab, key, prior_label, prior_spec,
                              label, filename) {
  ph  <- paste0("prob_H_", key)
  evc <- paste0("ev_", key)
  asc <- paste0("as_", key)
  ksc <- paste0("ks_", key)
  dcc <- paste0("dec_", key)

  tex <- c(
    "\\begin{table}[!h]",
    "\\centering",
    paste0("\\caption{TVSFP: FBST results under the ", prior_label,
           " prior (", prior_spec, "), posterior-based formulation, $M=",
           M_POST, "$. The weight $a^{*}$ ($b=1$) is calibrated prior-based ",
           "so that $\\alpha_{\\varphi}\\approx", TARGET, "$. ",
           "$P(H_{\\mathrm{np}}\\mid\\mathbf{X}) = ",
           "P(\\theta_1 \\leq \\theta_2 \\mid \\mathbf{X})$.}"),
    paste0("\\label{", label, "}"),
    "\\begin{tabular}{lcccccccc}",
    "\\toprule",
    paste0("\\textbf{Group} & $n$ & $x_1$ & $x_2$ & ",
           "$P(H_{\\mathrm{np}}\\mid\\mathbf{X})$ & ",
           "$ev(\\mathbf{H};\\mathbf{X})$ & $a^{*}$ & $k^*$ & Decision \\\\"),
    "\\midrule"
  )
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    tex <- c(tex, sprintf(
      "%s & %d & %d & %d & %.4f & %.4f & %.3f & %.4f & %s \\\\",
      r$group, r$n, r$x1, r$x2,
      r[[ph]], r[[evc]], r[[asc]], r[[ksc]], dec_sym(r[[dcc]])
    ))
  }
  tex <- c(tex, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, filename)
  cat(sprintf("-> %s\n", filename))
}

write_fbst_table(tab,
  key         = "non_informative",
  prior_label = "non-informative KL-optimal",
  prior_spec  = "$\\alpha_0=0.76,\\,\\alpha_1=0.76,\\,\\alpha_2=0.76$",
  label       = "tab:thks_ni",
  filename    = "output/thks_fbst_ni.tex")

write_fbst_table(tab,
  key         = "informative",
  prior_label = "informative",
  prior_spec  = "$\\alpha_0=24.99,\\,\\alpha_1=24.80,\\,\\alpha_2=24.88$",
  label       = "tab:thks_inf",
  filename    = "output/thks_fbst_inf.tex")

write_fbst_table(tab,
  key         = "conflict",
  prior_label = "informative (conflict)",
  prior_spec  = "$\\alpha_0=45.87,\\,\\alpha_1=5.11,\\,\\alpha_2=5.07$",
  label       = "tab:thks_conf",
  filename    = "output/thks_fbst_conf.tex")

tex_mcn <- c(
  "\\begin{table}[!h]",
  "\\centering",
  paste0("\\caption{TVSFP: McNemar's test for each treatment group. ",
         "$n_{01}$: pre$=0\\to$post$=1$ (improvement); ",
         "$n_{10}$: pre$=1\\to$post$=0$ (deterioration).}"),
  "\\label{tab:thks_mcnemar}",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "\\textbf{Group} & $n$ & $n_{01}$ & $n_{10}$ & $p$-value & Decision \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  tex_mcn <- c(tex_mcn, sprintf(
    "%s & %d & %d & %d & %.4f & %s \\\\",
    r$group, r$n, r$n01, r$n10, r$p_mcn, dec_sym(r$mcn_dec)
  ))
}
tex_mcn <- c(tex_mcn, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex_mcn, "output/thks_mcnemar.tex")
cat("-> output/thks_mcnemar.tex\n")
