# Single source of truth for the Olkin-Liu bivariate beta hyperparameters
# (alpha0, alpha1, alpha2) used throughout the project. Every script should
# load them with: source("R/01_priors_config.R")
#
# The three sets were obtained by Kullback-Leibler minimisation in
# R/02_fit_priors_kl.R (seed = 42, M = 1e5 / 5e3):
#   prior_NI   : argmin D_KL( U[0,1]^2 || p_alpha )              (non-informative)
#   prior_INF  : argmin D_KL( p_alpha || Beta(N*0.5, N*0.5)^2 ), N = 50
#   prior_CONF : argmin D_KL( p_alpha || Beta(N*0.1, N*0.9)^2 ), N = 50

# Non-informative prior (KL-closest to the bivariate uniform), ESS ~ 2.3
prior_NI <- c(a0 = 0.760595, a1 = 0.762204, a2 = 0.758178)

# Symmetric informative prior (mean 0.5, target ESS = 50)
prior_INF <- c(a0 = 24.9938, a1 = 24.8026, a2 = 24.8843)

# Informative prior in conflict with the data (mean 0.1, target ESS = 50);
# E[theta_j] = a_j / (a_j + a0) ~ 0.10, against THKS data (theta ~ 0.3-0.6)
prior_CONF <- c(a0 = 45.8654, a1 = 5.1072, a2 = 5.0664)

priors_summary <- data.frame(
  prior  = c("Non-informative (KL to uniform)",
             "Informative symmetric (mu=0.5, N=50)",
             "Conflict (mu=0.1, N=50)"),
  alpha0 = c(prior_NI["a0"],  prior_INF["a0"],  prior_CONF["a0"]),
  alpha1 = c(prior_NI["a1"],  prior_INF["a1"],  prior_CONF["a1"]),
  alpha2 = c(prior_NI["a2"],  prior_INF["a2"],  prior_CONF["a2"]),
  ESS_N0 = c(sum(prior_NI),   sum(prior_INF),   sum(prior_CONF)),
  row.names = NULL
)
