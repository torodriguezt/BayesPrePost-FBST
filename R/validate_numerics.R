# Numerical contracts and an execution pilot. No manuscript values are used in estimators.
source("R/pipeline_helpers.R")
source("R/02_fit_priors_kl.R")
fbst_validate_core <- function() {
  rows <- list()
  check <- function(name,observed,expected,tolerance,kind="numeric") {
    passed <- is.finite(observed) && abs(observed-expected)<=tolerance
    rows[[length(rows)+1L]] <<- data.frame(check=name,observed=observed,
      expected=expected,tolerance=tolerance,absolute_error=abs(observed-expected),passed=passed,kind=kind)
    invisible(passed)
  }
  fbst_validate_priors()
  a <- kl_vague_a()
  check("KL_root",a,0.7658512497788699578,2e-13)
  check("KL_divergence",kl_divergence(a,a,a),0.18575491070711565,2e-13)
  for (n in c(50,150)) {
    mass <- exp(lbeta(a,a+2*n)-lbeta(a,a))*2+
      2*n*exp(lbeta(a+2*n-1,a+1)-lbeta(a,a))
    check(paste0("singular_null_mass_n",n),mass,
          if(n==50)0.05964316530021684 else 0.02578703668055189,1e-12)
  }
  check("KL_all_failures",evalue(0,6,0,6,prior_NI),1,0)
  check("KL_single_failure_unequal",evalue(6,6,8,9,prior_NI),1,0)
  check("weak_mu_0.1_one_success",evalue(0,4,1,9,fbst_prior(.1,2)),1,0)
  left <- null_supremum(0,0,0,0,c(2,1,1))
  check("zero_endpoint_exponents_finite",as.numeric(!left$singular),1,0)
  check("zero_endpoint_exponents_max_at_zero",left$t,0,0)
  right <- null_supremum(4,4,4,4,c(2,1,1))
  check("finite_upper_endpoint_location",right$t,1,0)
  check("finite_upper_endpoint_log_limit",right$log_sup,-4*log(2),1e-13)
  check("independent_distinct_corner_condition",
        as.numeric(null_supremum(6,6,8,9,prior_NI,family="independent")$singular),0,0)
  check("tempered_singularity",evalue(6,6,2,6,prior_NI,power=.25),1,0)
  fractional_zero <- null_supremum(3,4,3,5,c(1,1,1),power=1/3)
  check("tempered_exact_zero_not_negative",as.numeric(fractional_zero$singular),0,0)
  check("tempered_exact_zero_endpoint_exponent",fractional_zero$v,0,0)
  check("flat_independent_analytic_evidence",evalue(0,1,1,1,c(1,1,1),family="independent"),
        (1+2*log(4))/16,2e-7)
  ev_flat <- evalue(0,1,1,1,c(1,1,1),family="independent")
  ev_prior <- evalue(0,1,1,1,c(1,1,1),family="independent",reference="prior")
  check("same_reference_for_uniform_independent_prior",ev_flat,ev_prior,2e-7)
  for (nm in c("KL","informative","conflict")) {
    p <- FBST_PRIORS[[nm]]
    A <- predictive_A(6,9,p)
    H <- predictive_H(6,9,p)
    check(paste0("predictive_A_sum_",nm),sum(A),1,1e-7)
    check(paste0("predictive_H_sum_",nm),sum(H),1,1e-10)
    target <- exp(lchoose(6,0:6)+lbeta(p[2]+0:6,p[1]+6-0:6)-lbeta(p[2],p[1]))
    check(paste0("predictive_A_marginal_",nm),max(abs(rowSums(A)-target)),0,1e-7)
  }
  # Finite decision spaces: exact ties and the no-rejection interval.
  EV <- matrix(c(.1,.1,.7,1),2)
  H <- matrix(c(.05,.05,.6,.3),2); A <- matrix(c(.3,.3,.2,.2),2)
  opt <- adaptive_cutoff(EV,H,A)
  check("calibration_exact_tie_alpha",opt$alpha,.1,1e-14)
  check("calibration_exact_tie_beta",opt$beta,.4,1e-14)
  same <- adaptive_cutoff(EV,H,H)
  check("calibration_all_rules_same_risk",min(same$alpha_k+same$beta_k),1,1e-14)
  # Published rounded values are validation targets only, never computational inputs.
  cases <- list(control=c(159,421,175,421),tv=c(145,416,201,416),toenail_A=c(54,146,14,133))
  targets <- c(control=.529071,tv=.000410017,toenail_A=7.39e-7)
  tolerances <- c(control=2e-6,tv=2e-8,toenail_A=1e-9)
  for (nm in names(cases)) {
    v <- cases[[nm]]
    ev <- evalue(v[1],v[2],v[3],v[4],prior_NI)
    check(paste0("application_evidence_",nm),ev,targets[[nm]],tolerances[[nm]],"manuscript_anchor")
    if (nm %in% c("control","toenail_A")) {
      s <- posterior_summary(v[1],v[2],v[3],v[4],prior_NI)
      if(nm=="control") check("control_direction",s$P_increase,.869672,2e-6,"manuscript_anchor")
      else {
        check("toenail_A_delta_mean",s$mean,-.259936,2e-6,"manuscript_anchor")
        check("toenail_A_delta_q025",s$lo,-.352937,3e-6,"manuscript_anchor")
        check("toenail_A_delta_q975",s$hi,-.165843,3e-6,"manuscript_anchor")
      }
    }
  }
  out <- do.call(rbind,rows)
  fbst_write_table(out,"validation_checks")
  writeLines(c("# Numerical validation",paste("Version:",FBST_VERSION),
      paste("Passed:",sum(out$passed),"/",nrow(out)),
      "Targets marked manuscript_anchor check the supplied numbers independently; differences are reported.",
      "A passed check is not certification of every enumeration or manuscript table."),
      file.path(fbst_output_dir(),"validation_report.md"))
  if(any(!out$passed & out$kind=="numeric")) stop("Numerical checks failed; see validation_checks.csv. Results were retained.")
  if(any(!out$passed & out$kind=="manuscript_anchor"))
    warning("Some manuscript anchors differ: investigate validation_checks.csv; targets are not imposed on the computation.")
  invisible(out)
}
fbst_run_pilot <- function() {
  sizes <- list(c(6,6),c(6,9),c(20,20))
  rows <- lapply(sizes,function(ns) {
    start <- proc.time()[[3]]
    d <- fbst_get_design(ns[1],ns[2],prior_NI)
    elapsed <- proc.time()[[3]]-start
    data.frame(n1=ns[1],n2=ns[2],points=length(d$EV),seconds=elapsed,
               object_bytes=as.numeric(object.size(d)),k=d$calibration$kstar,
               alpha=d$calibration$alpha,beta=d$calibration$beta,
               sum_H=sum(d$pH),sum_A=sum(d$pA))
  })
  out <- do.call(rbind,rows)
  fbst_write_table(out,"runtime_pilot")
  invisible(out)
}
if(sys.nframe()==0L) { fbst_record_environment(); fbst_validate_core(); fbst_run_pilot() }
