# Appendix estimator study: exact sampling distribution of the manuscript's
# 60x60 grid estimators. Sampling error is eliminated; grid approximation remains.
source("R/pipeline_helpers.R")
fbst_estimation_grid <- function(n,x1,x2,prior,G=60L) {
  xs <- seq(.001,.999,length.out=G)
  a0 <- prior[1]; a1 <- prior[2]; a2 <- prior[3]
  p1 <- (a1+x1-1)*log(xs)+(a2+a0+n-x1-1)*log1p(-xs)
  p2 <- (a2+x2-1)*log(xs)+(a1+a0+n-x2-1)*log1p(-xs)
  L <- outer(p1,p2,"+")-sum(prior)*log1p(-outer(xs,xs))
  z <- exp(L-max(L)); marg <- rowSums(z)
  mean <- sum(xs*marg)/sum(marg)
  candidates <- which(marg>=max(marg)*(1-1e-12))
  mode_grid <- xs[candidates[1]]
  # Marginal endpoint exponents from the hypergeometric marginal density.
  left <- a1+x1-1
  right <- min(a2+a0+n-x1-1,a0+2*n-x1-x2-1)
  log_boundary <- abs(n-x2-a2)<1e-14 && abs(right)<1e-14
  status <- if(left<0 && (right<0 || log_boundary)) "both_endpoints_unbounded" else
    if(left<0) "unbounded_at_zero" else if(right<0 || log_boundary) "unbounded_at_one" else
    if(length(candidates)>1) "multiple_grid_maxima" else
    if(candidates %in% c(1,G)) "grid_endpoint_not_interior_mode" else "interior_grid_maximum"
  c(mean=mean,grid_mode=mode_grid,
    boundary_mode=if(status=="unbounded_at_zero")0 else if(status=="unbounded_at_one")1 else mode_grid,
    status_code=match(status,c("both_endpoints_unbounded","unbounded_at_zero","unbounded_at_one",
                              "multiple_grid_maxima","grid_endpoint_not_interior_mode","interior_grid_maximum")))
}
fbst_weighted_quantile <- function(x,w,p) {
  o <- order(x); x[o][which(cumsum(w[o])/sum(w)>=p)[1]]
}
fbst_run_estimation <- function() {
  scenarios <- list(KL=list(prior=FBST_PRIORS$KL,theta=.1),
                    informative=list(prior=FBST_PRIORS$estimation,theta=.5),
                    conflict=list(prior=FBST_PRIORS$estimation,theta=.1))
  ns <- if(fbst_profile()=="pilot")c(2L,5L,10L) else FBST_SCENARIOS$appendix_n
  rows <- list()
  for(nm in names(scenarios)) {
    sc <- scenarios[[nm]]
    for(n in ns) {
      message("Appendix grid estimator: ",nm,", n=",n)
      values <- fbst_cache_compute("estimation_grid",list(n=n,prior=sc$prior,G=60L,
              algorithm="60x60-marginals-v2"),function() {
        xy <- expand.grid(x1=0:n,x2=0:n)
        v <- t(vapply(seq_len(nrow(xy)),function(i)
          fbst_estimation_grid(n,xy$x1[i],xy$x2[i],sc$prior),numeric(4)))
        cbind(xy,as.data.frame(v))
      })
      w <- dbinom(values$x1,n,sc$theta)*dbinom(values$x2,n,sc$theta)
      for(estimator in c("mean","grid_mode","boundary_mode")) {
        v <- values[[estimator]]
        rows[[length(rows)+1L]] <- data.frame(scenario=nm,n=n,theta=sc$theta,estimator=estimator,
          expectation=sum(w*v),bias=sum(w*(v-sc$theta)),rmse=sqrt(sum(w*(v-sc$theta)^2)),
          q025=fbst_weighted_quantile(v,w,.025),q975=fbst_weighted_quantile(v,w,.975),
          probability_unbounded_at_zero=sum(w[values$status_code==2]),
          probability_unbounded_at_one=sum(w[values$status_code==3]),
          sampling_method="finite_enumeration",estimator_method="60x60_grid",
          mc_standard_error=0)
      }
    }
  }
  out <- do.call(rbind,rows); fbst_write_table(out,"estimation_exact_sampling")
  # Refinement diagnostics for the actual grid approximation, not for sampling error.
  diagnostic <- list()
  for(nm in names(scenarios)) for(n in c(2,5,20,60)) for(x in unique(c(0,round(n*scenarios[[nm]]$theta),n))) {
    p <- scenarios[[nm]]$prior
    g60 <- fbst_estimation_grid(n,x,x,p,60); g120 <- fbst_estimation_grid(n,x,x,p,120)
    parameters <- .fbst_par(x,n,x,n,p)
    denominator <- .fbst_normalizer(parameters)$logZ
    numerator <- .fbst_normalizer(parameters+c(1,0,0,0,0))$logZ
    continuous_mean <- exp(numerator-denominator)
    diagnostic[[length(diagnostic)+1L]] <- data.frame(scenario=nm,n=n,x1=x,x2=x,
      mean_grid60=g60["mean"],mean_grid120=g120["mean"],difference=g120["mean"]-g60["mean"],
      continuous_mean=continuous_mean,grid60_error=g60["mean"]-continuous_mean,
      mode_grid60=g60["grid_mode"],mode_grid120=g120["grid_mode"],status_code=g60["status_code"])
  }
  fbst_write_table(do.call(rbind,diagnostic),"estimation_grid_refinement")
  figure <- function(draw,name) {
    for(ext in c("pdf","png")) {
      path <- file.path(fbst_figure_dir(),paste0(name,".",ext))
      if(ext=="pdf")pdf(path,width=10,height=4) else png(path,width=1800,height=720,res=180)
      tryCatch(draw(),finally=dev.off())
    }
  }
  names_est <- c(KL="est_fig4_mode_mean_noninformative",informative="est_fig6_mode_mean_informative",
                 conflict="est_fig8_mode_mean_conflict")
  for(nm in names(scenarios)) {
    figure(function() {
      par(mfrow=c(1,2),mar=c(4,4,3,1))
      for(est in c("grid_mode","mean")) {
        d <- out[out$scenario==nm & out$estimator==est,]
        plot(d$n,d$expectation,type="n",ylim=range(d$q025,d$q975,d$theta),
             xlab="n",ylab=if(est=="mean")"Posterior grid mean" else "Marginal grid mode",
             main=nm)
        polygon(c(d$n,rev(d$n)),c(d$q025,rev(d$q975)),col="#C8E8F2",border=NA)
        lines(d$n,d$expectation,col="#19768E",lwd=2);abline(h=d$theta[1],lty=2)
      }
    },names_est[[nm]])
  }
  density_names <- c(KL="est_fig3_prior_post_noninformative",informative="est_fig5_prior_post_informative",
                     conflict="est_fig7_prior_post_conflict")
  for(nm in names(scenarios)) {
    sc <- scenarios[[nm]]
    figure(function() {
      par(mfrow=c(1,2),mar=c(3,3,3,1))
      xs <- seq(.001,.999,length.out=80)
      for(n in c(0L,20L)) {
        x <- round(n*sc$theta);p <- sc$prior
        lp <- (p[2]+x-1)*log(xs)+(p[1]+p[3]+n-x-1)*log1p(-xs)
        L <- outer(lp,lp,"+")-sum(p)*log1p(-outer(xs,xs))
        z <- exp(L-max(L));z <- z/(sum(z)*diff(xs)[1]^2)
        persp(xs,xs,z,theta=40,phi=25,col="#D7B34B",border=NA,ticktype="simple",
              xlab="theta1",ylab="theta2",zlab="Grid density",
              main=if(n==0)"Prior (grid display)" else paste("Posterior n=20, x1=x2=",x))
      }
    },density_names[[nm]])
  }
  writeLines(c("# Appendix methodology changes",
    "The sampling distribution is enumerated, replacing 100 simulation replicates. No sampling seeds or Monte Carlo bands are involved.",
    "The retained estimator is explicitly the manuscript's 60x60 grid mean / marginal grid mode on [0.001,0.999].",
    "It remains a numerical approximation; see estimation_grid_refinement.csv before interpreting it as the continuous posterior mean.",
    "The q025/q975 bands describe the estimator's sampling distribution, not posterior uncertainty or the MC error of its average.",
    "boundary_mode replaces grid endpoints by 0 or 1 only when marginal density is analytically unbounded there.",
    "Interior grid maxima are not certified continuous modes; status_code identifies this limitation.",
    "Density figures are finite-grid displays; singular prior corners are not represented as finite exact density maxima."),
    file.path(fbst_output_dir(),"appendix_method_changes.md"))
  invisible(out)
}
if(sys.nframe()==0L)fbst_run_estimation()
