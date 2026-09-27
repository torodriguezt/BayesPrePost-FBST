# Manuscript specification. Tuples are (lambda0, lambda1, lambda2).
# N0 is the marginal concentration, not the sum of all three parameters.
FBST_VERSION <- "2.0.0"
FBST_KL_A <- uniroot(function(a) digamma(3*a)-digamma(a)-pi^2/6,
                     c(0.05,10), tol=1e-14)$root
fbst_prior <- function(mu0,N0) {
  stopifnot(length(mu0)==1L,mu0>0,mu0<1,length(N0)==1L,N0>0)
  c(a0=(1-mu0)*N0,a1=mu0*N0,a2=mu0*N0)
}
FBST_PRIORS <- list(KL=setNames(rep(FBST_KL_A,3),c("a0","a1","a2")),
                    informative=fbst_prior(.5,50),conflict=fbst_prior(.1,50),
                    estimation=setNames(rep(10,3),c("a0","a1","a2")))
prior_NI <- FBST_PRIORS$KL
prior_INF <- FBST_PRIORS$informative
prior_CONF <- FBST_PRIORS$conflict
priors_summary <- do.call(rbind,lapply(names(FBST_PRIORS),function(nm) {
  p <- FBST_PRIORS[[nm]]
  data.frame(prior=nm,alpha0=unname(p[1]),alpha1=unname(p[2]),alpha2=unname(p[3]),
             mu0=unname(p[2]/(p[1]+p[2])),ESS_N0=unname(p[1]+p[2]),lambda=sum(p))
}))
FBST_SCENARIOS <- list(
  baseline=.40,deltas=c(0,.05,.10,.20),main_n=c(50L,100L,150L,250L,400L,600L),
  paired_n=c(50L,100L,150L),coverage_deltas=c(0,.10,.20),
  coverage_psi=c(1,5,25),test_psi=c(1,3,5),extra_psi=.2,
  sensitivity_n=c(30L,50L,75L,100L,150L),
  weak_mu=c(.1,.3,.5,.7,.9),informative_mu=c(.30,.40,.45,.50,.60),
  conflict_mu=c(.05,.10,.15,.20,.25),
  boundary_sizes=matrix(c(20,20,50,50,20,100,100,20),ncol=2,byrow=TRUE),
  boundary_theta=c(.001,.01,.05,.5,.95,.99,.999),
  design_effects=c(1,1.5,2,3),appendix_n=c(5L,10L,20L,40L,60L),
  null_comparison_n=c(30L,50L,100L,200L,400L),independent_n=c(30L,100L,400L))
FBST_NUMERICS <- list(enum=list(method="grid",G=64L,max_G=256L,abs_tol=1e-5),
                      predictive=list(),near_risk=.002,predictive_mass_tol=1e-7,seed=42L)
