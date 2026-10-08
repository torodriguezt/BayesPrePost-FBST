# Prior definitions, fixed study scenarios and numerical controls.
# Shared by every study; edit here to change the article's settings.

# Manuscript specification. Tuples are (lambda0, lambda1, lambda2).
# N0 is the marginal concentration, not the sum of all three parameters.
ANALYSIS_VERSION <- "2.0.0"
KL_SHAPE <- uniroot(
  function(a) {
    digamma(3 * a) - digamma(a) - pi^2 / 6
  },
  c(0.05, 10),
  tol = 1e-14
)$root

beta_prior <- function(mu0, N0) {
  stopifnot(length(mu0) == 1L, mu0 > 0, mu0 < 1, length(N0) == 1L, N0 > 0)
  c(
    a0 = (1 - mu0) * N0,
    a1 = mu0 * N0,
    a2 = mu0 * N0
  )
}

PRIORS <- list(
  KL = setNames(rep(KL_SHAPE, 3), c("a0", "a1", "a2")),
  informative = beta_prior(.5, 50),
  conflict = beta_prior(.1, 50),
  estimation = setNames(rep(10, 3), c("a0", "a1", "a2"))
)
SCENARIOS <- list(
  baseline = .40,
  deltas = c(0, .05, .10, .20),
  main_n = c(50L, 100L, 150L, 250L, 400L, 600L),
  paired_n = c(50L, 100L, 150L),
  coverage_deltas = c(0, .10, .20),
  association_psi = c(.2, 1, 3, 5),
  overlap_psi = c(.2, 1, 3, 5, 25),
  sensitivity_n = c(30L, 50L, 75L, 100L, 150L),
  weak_mu = c(.1, .3, .5, .7, .9),
  informative_mu = c(.30, .40, .45, .50, .60),
  conflict_mu = c(.05, .10, .15, .20, .25),
  boundary_sizes = matrix(c(20, 20, 50, 50, 20, 100), ncol = 2, byrow = TRUE),
  boundary_theta = c(.001, .01, .05, .5, .95, .99, .999),
  appendix_n = c(5L, 10L, 20L, 40L, 60L),
  null_comparison_n = c(30L, 50L, 100L, 200L, 400L),
  independent_n = c(30L, 100L, 400L)
)
NUMERICS <- list(
  enum = list(
    method = "grid",
    G = 64L,
    max_G = 256L,
    abs_tol = 1e-5
  ),
  predictive = list(),
  near_risk = .002,
  predictive_mass_tol = 1e-7
)
