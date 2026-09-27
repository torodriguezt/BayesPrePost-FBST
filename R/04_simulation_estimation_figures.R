# Appendix: the posterior mean as a point estimator. Posterior means are ratios of
# series normalisers and the sampling distribution is enumerated, so the bias/RMSE
# table carries no simulation error. The prior/posterior figure is a display grid only.
if (!exists("fbst_app_contour_grid", mode = "function")) source("R/manuscript_applications.R")
fbst_estimation_scenarios <- function() list(
  KL=list(prior=FBST_PRIORS$KL,theta=.1,label="KL-optimal",display_cap=4),
  centred=list(prior=FBST_PRIORS$estimation,theta=.5,label="centred",display_cap=Inf),
  conflict=list(prior=FBST_PRIORS$estimation,theta=.1,label="conflict",display_cap=Inf))
fbst_posterior_mean1 <- function(x1,n1,x2,n2,prior) {
  # E(theta1 | x) = Z(A1+1,B1,A2,B2,lambda) / Z(A1,B1,A2,B2,lambda).
  p <- .fbst_par(x1,n1,x2,n2,prior)
  exp(.fbst_normalizer(p+c(1,0,0,0,0))$logZ-.fbst_normalizer(p)$logZ)
}
fbst_run_estimation <- function() {
  sc <- fbst_estimation_scenarios()
  ns <- if(fbst_profile()=="pilot")c(5L,10L) else FBST_SCENARIOS$appendix_n
  rows <- list()
  for(nm in names(sc)) for(n in ns) {
    message("Appendix posterior mean: ",nm,", n=",n)
    values <- fbst_cache_compute("estimation_posterior_mean",list(n=n,prior=sc[[nm]]$prior,
              algorithm="series-normaliser-ratio-v1"),function() {
      xy <- expand.grid(x1=0:n,x2=0:n)
      cbind(xy,mean=mapply(fbst_posterior_mean1,xy$x1,n,xy$x2,n,MoreArgs=list(prior=sc[[nm]]$prior)))
    })
    th <- sc[[nm]]$theta
    w <- dbinom(values$x1,n,th)*dbinom(values$x2,n,th)
    fbst_check_mass(w,paste("estimation",nm,n))
    rows[[length(rows)+1L]] <- data.frame(scenario=nm,n=n,theta=th,
      expectation=sum(w*values$mean),bias=sum(w*(values$mean-th)),
      rmse=sqrt(sum(w*(values$mean-th)^2)),rmse_sample_proportion=sqrt(th*(1-th)/n),
      probability_no_successes=(1-th)^n,sampling_method="finite_enumeration",mc_standard_error=0)
  }
  out <- do.call(rbind,rows); fbst_write_table(out,"estimation_posterior_mean")
  fbst_estimation_figure()
  invisible(out)
}
fbst_estimation_figure <- function(n=20L,G=60L) {
  # Exactly normalised densities on a midpoint display grid. The KL-optimal prior is
  # unbounded as theta_j -> 0 and at the vertex (1,1), so its surface is capped at display_cap.
  sc <- fbst_estimation_scenarios()
  xs <- (seq_len(G)-.5)/G
  density <- function(p) exp(outer(xs,xs,function(a,b) .fbst_logkernel(a,b,p))-.fbst_normalizer(p)$logZ)
  fbst_app_plot_files("est_prior_posterior",function() {
    graphics::par(mfcol=c(2,3),mar=c(1,1.5,2.4,.5))
    for(nm in names(sc)) {
      p <- sc[[nm]]$prior; th <- sc[[nm]]$theta; x <- round(n*th)
      surfaces <- list(prior=pmin(density(.fbst_par(0,0,0,0,p)),sc[[nm]]$display_cap),
                       posterior=density(.fbst_par(x,n,x,n,p)))
      for(k in names(surfaces)) {
        pm <- graphics::persp(xs,xs,surfaces[[k]],theta=40,phi=25,expand=.75,col="#D7B34B",
          border="grey30",lwd=.3,ticktype="detailed",nticks=4,cex.axis=.65,
          xlab="",ylab="",zlab="density",zlim=c(0,max(surfaces[[k]])))
        graphics::text(grDevices::trans3d(.5,-.28,0,pm),expression(theta[1]))
        graphics::text(grDevices::trans3d(1.28,.5,0,pm),expression(theta[2]))
        graphics::title(main=if(k=="prior")bquote(.(paste0("Prior: ",sc[[nm]]$label," ("))*theta==.(th)*")") else
          bquote(list("Posterior",n==.(n),x[1]==.(x),x[2]==.(x))),cex.main=.95)
      }
    }
  },width=9,height=6)
}
if(sys.nframe()==0L)fbst_run_estimation()
