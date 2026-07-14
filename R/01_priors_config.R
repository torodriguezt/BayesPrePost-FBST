# Hyperparameters of the three Olkin-Liu priors used across the project,
# fitted by KL minimisation in R/02_fit_priors_kl.R.

# Non-informative, closest to the bivariate uniform (ESS 2.3)
prior_NI <- c(a0 = 0.760595, a1 = 0.762204, a2 = 0.758178)

# Informative, mean 0.5, ESS 50
prior_INF <- c(a0 = 24.9938, a1 = 24.8026, a2 = 24.8843)

# Informative in conflict with the THKS data, mean 0.1, ESS 50
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
