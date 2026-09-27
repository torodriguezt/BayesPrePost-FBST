# Validate the analytic KL optimum. Informative priors are specified directly.
source("R/pipeline_helpers.R")
fbst_validate_priors <- function() {
  a <- kl_vague_a()
  p <- priors_summary
  p$kl_uniform_to_prior <- vapply(FBST_PRIORS,function(x) do.call(kl_divergence,as.list(unname(x))),0.0)
  stopifnot(abs(a-FBST_KL_A)<1e-12,
            abs(digamma(3*a)-digamma(a)-pi^2/6)<1e-12,
            identical(unname(prior_INF),c(25,25,25)),
            identical(unname(prior_CONF),c(45,5,5)))
  fbst_write_table(p,"kl_priors")
  invisible(p)
}
if (sys.nframe()==0L) fbst_validate_priors()
