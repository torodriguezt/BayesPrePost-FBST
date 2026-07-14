# Fits the three Olkin-Liu priors by KL minimisation: non-informative
# (closest to the uniform), informative (mean 0.5) and conflict (mean 0.1).
# The fitted values are cached in R/01_priors_config.R, so this only needs
# to be re-run if the targets change. Writes output/kl_priors.rds.
# Runtime: about 5 min.

library(Rcpp)
sourceCpp("src/BivBetaBinom.cpp")

# Olkin-Liu log density in closed form.
log_p_olkin_liu <- function(y1, y2, a0, a1, a2) {
  a <- a0 + a1 + a2
  lgamma(a) - lgamma(a0) - lgamma(a1) - lgamma(a2) +
    (a1 - 1) * log(y1) + (a2 - 1) * log(y2) +
    (a0 + a2 - 1) * log1p(-y1) + (a0 + a1 - 1) * log1p(-y2) -
    a * log1p(-y1 * y2)
}

# Monte Carlo estimate of D_KL(uniform || p_alpha).
kl_U_to_p <- function(par, theta_unif) {
  a <- exp(par)  # keeps alpha positive
  -mean(log_p_olkin_liu(theta_unif[, 1], theta_unif[, 2], a[1], a[2], a[3]))
}

# Monte Carlo estimate of D_KL(p_alpha || target).
kl_p_to_target <- function(par, log_target_fun, M, eps = 1e-10) {
  a <- exp(par)
  th <- sample_prior(M, a[1], a[2], a[3])
  y1 <- pmin(pmax(th[, 1], eps), 1 - eps)
  y2 <- pmin(pmax(th[, 2], eps), 1 - eps)
  mean(log_p_olkin_liu(y1, y2, a[1], a[2], a[3]) - log_target_fun(y1, y2))
}

fit_noninformative_kl <- function(M = 100000, seed = 42, init = c(1, 1, 1)) {
  set.seed(seed)
  theta_unif <- cbind(runif(M), runif(M))
  res <- optim(log(init), kl_U_to_p, theta_unif = theta_unif,
               method = "Nelder-Mead",
               control = list(reltol = 1e-10, maxit = 5000))
  alpha_hat <- exp(res$par)
  list(alpha = alpha_hat,
       kl    = res$value,
       N0    = sum(alpha_hat),
       conv  = res$convergence,
       M     = M)
}

make_log_target_indep_beta <- function(mu = c(0.5, 0.5), N = 50) {
  a_t <- N * mu
  b_t <- N * (1 - mu)
  function(y1, y2) {
    dbeta(y1, a_t[1], b_t[1], log = TRUE) +
      dbeta(y2, a_t[2], b_t[2], log = TRUE)
  }
}

fit_informative_kl <- function(mu = c(0.5, 0.5), N_target = 50,
                               M = 5000, n_replicates = 10,
                               seed = 42, init = c(10, 10, 10)) {
  log_target <- make_log_target_indep_beta(mu, N_target)
  # replicate average smooths Monte Carlo noise in the objective
  obj <- function(par) {
    set.seed(seed)
    mean(replicate(n_replicates, kl_p_to_target(par, log_target, M)))
  }
  res <- optim(log(init), obj, method = "Nelder-Mead",
               control = list(reltol = 1e-6, maxit = 3000))
  alpha_hat <- exp(res$par)
  list(alpha   = alpha_hat,
       kl      = res$value,
       N0      = sum(alpha_hat),
       conv    = res$convergence,
       target  = list(mu = mu, N = N_target),
       M = M, n_replicates = n_replicates)
}

cat("Fitting Olkin-Liu priors by KL minimisation\n\n")
dir.create("output", showWarnings = FALSE)

cat("[1] Non-informative prior: min D_KL(U || p_alpha)\n")
fit_NI <- fit_noninformative_kl(M = 100000, seed = 42)
cat(sprintf("    alpha = (%.6f, %.6f, %.6f) | ESS = %.3f | KL = %.6f\n\n",
            fit_NI$alpha[1], fit_NI$alpha[2], fit_NI$alpha[3],
            fit_NI$N0, fit_NI$kl))

cat("[2] Informative prior: target mu = (0.5, 0.5), N = 50\n")
fit_INF_sim <- fit_informative_kl(mu = c(0.5, 0.5), N_target = 50,
                                  M = 5000, n_replicates = 10)
cat(sprintf("    alpha = (%.4f, %.4f, %.4f) | ESS = %.3f | KL = %.6f\n\n",
            fit_INF_sim$alpha[1], fit_INF_sim$alpha[2], fit_INF_sim$alpha[3],
            fit_INF_sim$N0, fit_INF_sim$kl))

cat("[3] Conflict prior: target mu = (0.1, 0.1), N = 50\n")
fit_INF_conf <- fit_informative_kl(mu = c(0.1, 0.1), N_target = 50,
                                   M = 5000, n_replicates = 10,
                                   init = c(40, 5, 5))
cat(sprintf("    alpha = (%.4f, %.4f, %.4f) | ESS = %.3f | KL = %.6f\n\n",
            fit_INF_conf$alpha[1], fit_INF_conf$alpha[2], fit_INF_conf$alpha[3],
            fit_INF_conf$N0, fit_INF_conf$kl))

resumen <- data.frame(
  prior = c("Non-informative (KL to uniform)",
            "Informative symmetric (mu=0.5, N=50)",
            "Informative conflict (mu=0.1, N=50)"),
  alpha0 = c(fit_NI$alpha[1], fit_INF_sim$alpha[1], fit_INF_conf$alpha[1]),
  alpha1 = c(fit_NI$alpha[2], fit_INF_sim$alpha[2], fit_INF_conf$alpha[2]),
  alpha2 = c(fit_NI$alpha[3], fit_INF_sim$alpha[3], fit_INF_conf$alpha[3]),
  ESS_N0 = c(fit_NI$N0, fit_INF_sim$N0, fit_INF_conf$N0),
  KL_min = c(fit_NI$kl, fit_INF_sim$kl, fit_INF_conf$kl)
)
cat("Summary:\n")
print(resumen, row.names = FALSE)

saveRDS(list(noninformative       = fit_NI,
             informative_sym      = fit_INF_sim,
             informative_conflict = fit_INF_conf,
             resumen              = resumen),
        "output/kl_priors.rds")
cat("\n-> output/kl_priors.rds\n")
