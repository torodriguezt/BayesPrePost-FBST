# Deterministic FBST computations. Prior order: (lambda0, lambda1, lambda2).
# Version changes intentionally invalidate numerical caches.
FBST_CORE_VERSION <- "r-2026-09-24-1"
.fbst_env <- new.env(parent=baseenv())
.fbst_sources <- vapply(sys.frames(),function(frame) {
  if(is.null(frame$ofile)) NA_character_ else as.character(frame$ofile)[1L]
},character(1))
.fbst_sources <- .fbst_sources[!is.na(.fbst_sources) & basename(.fbst_sources)=="exact_fbst.R"]
.fbst_path <- normalizePath(if(length(.fbst_sources))tail(.fbst_sources,1L) else "R/exact_fbst.R",
                           mustWork=TRUE)
.fbst_axis_cache <- new.env(parent=emptyenv())
.fbst_cpp <- file.path(dirname(.fbst_path), "exact_fbst_core.cpp")
if(!file.exists(.fbst_cpp))stop("Required numerical companion not found: ",.fbst_cpp)
.fbst_source_hash <- unname(tools::md5sum(c(.fbst_path,.fbst_cpp)))
if (requireNamespace("Rcpp", quietly=TRUE) && file.exists(.fbst_cpp)) {
  Rcpp::sourceCpp(.fbst_cpp, env=.fbst_env, rebuild=FALSE, showOutput=FALSE)
}

kl_vague_a <- function() uniroot(function(a) digamma(3*a)-digamma(a)-pi^2/6,c(.05,10),tol=1e-14)$root
kl_divergence <- function(a0,a1,a2) lgamma(a0)+lgamma(a1)+lgamma(a2)-lgamma(a0+a1+a2)+(a0+a1+a2)*pi^2/6-4
prior_mu_N0 <- function(mu0,N0) {
  stopifnot(length(mu0)==1L,mu0>0,mu0<1,N0>0)
  c(lambda0=(1-mu0)*N0,lambda1=mu0*N0,lambda2=mu0*N0)
}
.fbst_par <- function(x1,n1,x2,n2,prior,family="olkin_liu",power=1) {
  family <- match.arg(family,c("olkin_liu","independent"))
  if(length(prior)!=3 || any(!is.finite(prior)) || any(prior<=0)) stop("Three positive prior parameters required")
  if(any(!is.finite(c(x1,n1,x2,n2,power))) || power<=0 || min(x1,n1-x1,x2,n2-x2)<0) stop("Invalid counts or likelihood power")
  a0<-prior[1];a1<-prior[2];a2<-prior[3]
  c(A1=a1+power*x1,B1=a0+if(family=="olkin_liu")a2+power*(n1-x1) else power*(n1-x1),
    A2=a2+power*x2,B2=a0+if(family=="olkin_liu")a1+power*(n2-x2) else power*(n2-x2),
    lambda=if(family=="olkin_liu")sum(prior) else 0)
}
.fbst_xlog <- function(a,x) if(a==0) rep(0,length(x)) else a*log(x)
.fbst_logkernel <- function(x,y,p) .fbst_xlog(p[1]-1,x)+.fbst_xlog(p[2]-1,1-x)+.fbst_xlog(p[3]-1,y)+.fbst_xlog(p[4]-1,1-y)-.fbst_xlog(p[5],1-x*y)
.fbst_roots <- function(a,b,c) {
  if(abs(a)<1e-14*max(1,abs(b),abs(c))) return(if(abs(b)<1e-14) numeric() else -c/b)
  d<-b*b-4*a*c
  if(d<0)return(numeric())
  q<--.5*(b+if(b>=0)sqrt(d) else -sqrt(d))
  if(q==0)return(-b/(2*a))
  unique(c(q/a,c/q))
}
# Global supremum includes finite boundary limits and all stationary roots.
null_supremum <- function(x1,n1,x2,n2,prior,reference="flat",family="olkin_liu",power=1) {
  reference<-match.arg(reference,c("flat","prior"));p<-.fbst_par(x1,n1,x2,n2,prior,family,power)
  r<-if(reference=="prior") c(power*x1,power*(n1-x1),power*x2,power*(n2-x2),0) else c(p[1:4]-1,p[5])
  # Direct expressions prevent cancellation from turning an exact zero exponent
  # (for example a0=1, F=3, power=1/3) into a spurious negative exponent.
  S<-x1+x2;F<-n1+n2-S;N<-n1+n2;lam<-r[5]
  if(reference=="prior") {
    u<-power*S;v<-power*F;quadratic<--power*N
  } else {
    u<-prior[2]+prior[3]+power*S-2
    if(family=="olkin_liu") {
      v<-prior[1]+power*F-2;quadratic<-4-power*N
    } else {
      v<-2*prior[1]+power*F-2
      quadratic<-4-(2*prior[1]+prior[2]+prior[3])-power*N
    }
  }
  if(u<0 || v<0)return(list(log_sup=Inf,t=NA_real_,singular=TRUE,u=u,v=v))
  candidates<-c(0,1,.fbst_roots(quadratic,-v-lam,u))
  candidates<-sort(unique(candidates[candidates>=0 & candidates<=1]))
  vals<-.fbst_xlog(u,candidates)+.fbst_xlog(v,1-candidates)-lam*log1p(candidates)
  j<-which.max(vals)
  list(log_sup=unname(vals[j]),t=candidates[j],singular=FALSE,u=unname(u),v=unname(v),candidates=candidates,log_values=vals)
}
.fbst_gauss <- function(A,B,G) {
  if(!requireNamespace("statmod",quietly=TRUE))stop("Package 'statmod' is required: install.packages('statmod')")
  key<-paste(sprintf("%.17g",c(A,B,G)),collapse=":")
  if(exists(key,envir=.fbst_axis_cache,inherits=FALSE))return(get(key,envir=.fbst_axis_cache,inherits=FALSE))
  q<-statmod::gauss.quad.prob(G,dist="beta",alpha=A,beta=B)
  if(any(q$nodes<=0|q$nodes>=1)||any(q$weights<0))stop("Invalid Gaussian quadrature nodes")
  assign(key,q,envir=.fbst_axis_cache)
  q
}
.fbst_grid <- function(p,G) {
  q1<-.fbst_gauss(p[1],p[2],G);q2<-.fbst_gauss(p[3],p[4],G)
  L<--p[5]*log1p(-outer(q1$nodes,q2$nodes)); shift<-max(L)
  W<-exp(L-shift)*outer(q1$weights,q2$weights)
  z<-sum(W)
  list(x=q1$nodes,y=q2$nodes,w=W/z,logZ=lbeta(p[1],p[2])+lbeta(p[3],p[4])+shift+log(z))
}
.fbst_normalizer <- function(p,rel_tol=1e-8,max_G=1024) {
  if(p[5]==0)return(list(logZ=lbeta(p[1],p[2])+lbeta(p[3],p[4]),log_lower=lbeta(p[1],p[2])+lbeta(p[3],p[4])-1e-12,error=0,G=0))
  if(exists("fbst_normalizer_series",envir=.fbst_env,inherits=FALSE)) {
    z<-.fbst_env$fbst_normalizer_series(unname(p),rel_tol/10,1048576L)
    if(isTRUE(z$unresolved))stop("Posterior series normalizer unresolved: ",z$error)
    return(z)
  }
  old<-NA_real_
  for(G in unique(pmin(max_G,2^(5:11)))) {
    g<-.fbst_grid(p,G);err<-if(is.na(old))Inf else abs(expm1(g$logZ-old))
    if(err<rel_tol)return(list(logZ=g$logZ,log_lower=lbeta(p[1],p[2])+lbeta(p[3],p[4])-1e-12,error=err,G=G))
    old<-g$logZ
    if(G>=max_G)break
  }
  stop(sprintf("Posterior normalization did not converge: relative difference %.3g at G=%d",err,G))
}
.fbst_integrate <- function(f,lo=0,hi=1,rel_tol=1e-8,abs_tol=1e-11,.depth=0L) {
  if(lo>=hi)return(0)
  z<-integrate(f,lo,hi,rel.tol=rel_tol,abs.tol=abs_tol,subdivisions=600L,stop.on.error=FALSE)
  if(z$message!="OK") {
    if(.depth<4L) {
      mid<-(lo+hi)/2
      return(.fbst_integrate(f,lo,mid,rel_tol,abs_tol/2,.depth+1L)+.fbst_integrate(f,mid,hi,rel_tol,abs_tol/2,.depth+1L))
    }
    stop("Adaptive integration failed after interval subdivision: ",z$message)
  }
  z$value
}
# Conditional odds transform resolves the joint (1,1) corner without truncating.
.fbst_slice_integral <- function(x,lo,hi,p,shift,rel_tol,abs_tol) {
  if(lo>=hi)return(0)
  if(p[5]==0)return(exp((p[1]-1)*log(x)+(p[2]-1)*log1p(-x)+lbeta(p[3],p[4])-shift)*
    if(lo==0) pbeta(hi,p[3],p[4]) else if(hi==1)pbeta(lo,p[3],p[4],lower.tail=FALSE) else if(pbeta(lo,p[3],p[4])>.5) (pbeta(lo,p[3],p[4],lower.tail=FALSE)-pbeta(hi,p[3],p[4],lower.tail=FALSE)) else (pbeta(hi,p[3],p[4])-pbeta(lo,p[3],p[4])))
  if(p[4]-p[5]>0) {
    # If the conditional limit at x=1 remains integrable, original coordinates
    # avoid the narrow artificial peak introduced by the odds transformation.
    cc<-(p[1]-1)*log(x)+(p[2]-1)*log1p(-x)-shift
    direct_piece<-function(l,h) {
      if(p[3]<1) return(.fbst_integrate(function(z) {
        y<-h*z^(1/p[3]);exp(cc+p[3]*log(h)-log(p[3])+(p[4]-1)*log1p(-y)-p[5]*log1p(-x*y))
      },(l/h)^p[3],1,rel_tol,abs_tol))
      b<-min(p[4],p[4]-p[5])
      if(h==1 && b<1) return(.fbst_integrate(function(z) {
        t<-(1-l)*z^(1/b);exp(cc+p[4]*log1p(-l)-log(b)+(p[4]/b-1)*log(z)+(p[3]-1)*log1p(-t)-p[5]*log((1-x)+x*t))
      },0,1,rel_tol,abs_tol))
      .fbst_integrate(function(y)exp(cc+(p[3]-1)*log(y)+(p[4]-1)*log1p(-y)-p[5]*log1p(-x*y)),l,h,rel_tol,abs_tol)
    }
    return(if(lo<.5&&hi==1)direct_piece(lo,.5)+direct_piece(.5,1) else direct_piece(lo,hi))
  }
  trans<-function(y) if(y==1)1 else if(y==0)0 else (1-x)*y/(1-x*y)
  vl<-trans(lo);vh<-trans(hi)
  cst<-(p[1]-1)*log(x)+(p[2]+p[4]-p[5]-1)*log1p(-x)-shift
  logrest<-function(v)cst+(p[5]-p[3]-p[4])*log1p(-x+x*v)
  piece<-function(vl,vh) {
    if(p[3]<1) {
      # v=vh*z^(1/A2) cancels the integrable power singularity exactly.
      return(.fbst_integrate(function(z){v<-vh*z^(1/p[3]);exp(logrest(v)+(p[4]-1)*log1p(-v)+p[3]*log(vh)-log(p[3]))},(vl/vh)^p[3],1,rel_tol,abs_tol))
    }
    if(vh==1 && p[4]<1) {
      return(.fbst_integrate(function(z){w<-(1-vl)*z^(1/p[4]);v<-1-w;exp(logrest(v)+(p[3]-1)*log(v)+p[4]*log1p(-vl)-log(p[4]))},0,1,rel_tol,abs_tol))
    }
    .fbst_integrate(function(v)exp(logrest(v)+(p[3]-1)*log(v)+(p[4]-1)*log1p(-v)),vl,vh,rel_tol,abs_tol)
  }
  tryCatch(if(vl==0&&vh==1&&p[3]<1&&p[4]<1)piece(0,.5)+piece(.5,1) else piece(vl,vh), error=function(e)stop(conditionMessage(e),"; slice x=",format(x,digits=17),", y interval=",paste(format(c(lo,hi),digits=17),collapse=",")))
}
# Integrate the complement directly. Each slice is split at all stationary points.
.fbst_complement_intervals <- function(x,r,target) {
  cst<-.fbst_xlog(r[1],x)+.fbst_xlog(r[2],1-x)
  f<-function(y)cst+.fbst_xlog(r[3],y)+.fbst_xlog(r[4],1-y)-r[5]*log1p(-x*y)-target
  st<-.fbst_roots(x*(r[3]+r[4]-r[5]),-r[3]*(1+x)-r[4]+r[5]*x,r[3])
  pts<-sort(unique(c(0,st[st>0&st<1],1))); roots<-numeric()
  for(j in seq_len(length(pts)-1L)) {
    l<-pts[j];h<-pts[j+1L];fl<-f(l);fh<-f(h)
    if(is.finite(fl)&&fl==0)roots<-c(roots,l)
    if(is.finite(fh)&&fh==0)roots<-c(roots,h)
    if(!is.nan(fl)&&!is.nan(fh)&&sign(fl)*sign(fh)<0) {
      rr<-if(l==0 && r[3]!=0) {
        lf<-f(exp(-744))
        if(sign(lf)*sign(fh)>=0) 0 else exp(uniroot(function(z)f(exp(z)),c(-744,log(h)),tol=1e-12)$root)
      } else uniroot(f,c(l,h),tol=1e-13)$root
      roots<-c(roots,rr)
    }
  }
  b<-sort(unique(c(0,roots,1)))
  ans<-matrix(numeric(),ncol=2)
  for(j in seq_len(length(b)-1L))if(f(mean(b[j:(j+1L)]))<=0)ans<-rbind(ans,b[j:(j+1L)])
  ans
}
evalue <- function(x1,n1,x2,n2,prior,G=64,W=NULL,reference="flat",family="olkin_liu",power=1,
                   method=c("adaptive","grid"),rel_tol=1e-8,abs_tol=1e-10,max_G=1024,strict=TRUE) {
  method<-match.arg(method);reference<-match.arg(reference,c("flat","prior"))
  p<-.fbst_par(x1,n1,x2,n2,prior,family,power)
  if(method=="adaptive" && p[1]/(p[1]+p[2])>p[3]/(p[3]+p[4])) {
    # Exchange parameter coordinates AND their prior shapes; this is a
    # measure-preserving symmetry, not a success/failure complement.
    return(evalue(x2,n2,x1,n1,prior[c(1,3,2)],G=G,reference=reference,family=family,power=power,method=method,rel_tol=rel_tol,abs_tol=abs_tol,max_G=max_G,strict=strict))
  }
  ns<-null_supremum(x1,n1,x2,n2,prior,reference,family,power)
  if(ns$singular)return(structure(1,error_estimate=0,log_evalue=0,diagnostics=list(method="analytic",singular=TRUE)))
  r<-if(reference=="prior")c(power*x1,power*(n1-x1),power*x2,power*(n2-x2),0) else c(p[1:4]-1,p[5])
  if(method=="grid") {
    old<-NA_real_
    for(g in unique(pmin(max_G,G*2^(0:10)))) {
      q<-.fbst_grid(p,g)
      if(exists("fbst_grid_ev",envir=.fbst_env,inherits=FALSE)) val<-.fbst_env$fbst_grid_ev(q$x,q$y,q$w,p,r,ns$log_sup) else {
        L<-outer(q$x,q$y,Vectorize(function(x,y).fbst_logkernel(x,y,c(r[1:4]+1,r[5]))))
        val<-sum(q$w*(L<=ns$log_sup))
      }
      err<-if(is.na(old))Inf else abs(val-old)
      bound<-Inf
      if(val<abs_tol && reference=="flat") {
        norm_bound<-.fbst_normalizer(p,rel_tol=rel_tol,max_G=max_G)
        # A finite partial sum of the positive series gives a lower normalizer;
        # hence exp(sup)/Z_lower bounds the complement (the square has area one).
        logbound<-ns$log_sup-norm_bound$log_lower
        bound<-max(exp(logbound),.Machine$double.xmin)
      }
      isbound<-val<abs_tol && is.finite(bound) && bound<abs_tol
      if(isbound){val<-bound;err<-bound;break}
      if(err<abs_tol&&val>0)break
      old<-val;if(g>=max_G)break
    }
    unresolved<-(!is.finite(err)||err>abs_tol||val==0)
    if(strict&&unresolved)stop(sprintf("E-value grid unresolved: difference %.3g at G=%d",err,g))
    return(structure(val,error_estimate=err,log_evalue=log(val),diagnostics=list(method="grid",G=g,unresolved=unresolved,upper_bound=isbound,log_upper_bound=if(isbound)logbound else NA_real_)))
  }
  norm<-.fbst_normalizer(p,rel_tol=rel_tol/5,max_G=max_G)
  # Flat-reference complement is bounded by exp(log_sup), retaining tiny probabilities.
  shift<-if(reference=="flat")ns$log_sup else norm$logZ
  fun<-function(xx)vapply(xx,function(x) {
    intervals<-.fbst_complement_intervals(x,r,ns$log_sup)
    if(!nrow(intervals))return(0)
    sum(vapply(seq_len(nrow(intervals)),function(j) .fbst_slice_integral(x,intervals[j,1],intervals[j,2],p,shift,rel_tol,abs_tol/10),numeric(1)))
  },numeric(1))
  z<-.fbst_integrate(fun,rel_tol=rel_tol,abs_tol=abs_tol)
  if(z<=0)stop("Complement integral is numerically zero; increase numerical precision")
  logev<-log(z)+shift-norm$logZ
  if(logev<log(.Machine$double.xmin))stop(sprintf("e-value below double range; log(e-value)=%.16g. Use log scale / arbitrary precision",logev))
  ev<-exp(logev)
  if(ev>1+max(1e-7,rel_tol*10))stop("e-value exceeds one: integration inconsistency")
  structure(min(ev,1),error_estimate=ev*(norm$error+rel_tol)+abs_tol*exp(shift-norm$logZ),log_evalue=logev,
            diagnostics=list(method="adaptive",normalizer_G=norm$G,normalizer_error=norm$error,singular=FALSE))
}
enumerate_evalues <- function(n1,n2,prior,G=64,reference="flat",family="olkin_liu",power=1,
                              method="grid",max_G=256,abs_tol=1e-5,rel_tol=1e-8,strict=FALSE,progress=FALSE,checkpoint_file=NULL) {
  stopifnot(n1==as.integer(n1),n2==as.integer(n2),n1>=0,n2>=0)
  started<-proc.time()[3];EV<-matrix(NA_real_,n1+1,n2+1);errors<-EV;unresolved<-matrix(FALSE,n1+1,n2+1);upper_bound<-unresolved
  symmetric<-n1==n2 && identical(unname(prior[2]),unname(prior[3]))
  signature<-list(version=FBST_CORE_VERSION,core_hash=.fbst_source_hash,n1=n1,n2=n2,prior=prior,G=G,reference=reference,family=family,power=power,method=method,max_G=max_G,abs_tol=abs_tol,rel_tol=rel_tol)
  if(!is.null(checkpoint_file)&&file.exists(checkpoint_file)) {
    saved<-readRDS(checkpoint_file)
    if(identical(saved$signature,signature)){EV<-saved$EV;errors<-saved$errors;unresolved<-saved$unresolved;upper_bound<-saved$upper_bound}
  }
  for(i in 0:n1) {
    for(j in 0:n2) {
      if(!is.na(EV[i+1,j+1]))next
      if(symmetric&&j<i){EV[i+1,j+1]<-EV[j+1,i+1];errors[i+1,j+1]<-errors[j+1,i+1];unresolved[i+1,j+1]<-unresolved[j+1,i+1];upper_bound[i+1,j+1]<-upper_bound[j+1,i+1];next}
      z<-evalue(i,n1,j,n2,prior,G=G,reference=reference,family=family,power=power,method=method,
                max_G=max_G,abs_tol=abs_tol,rel_tol=rel_tol,strict=strict)
      EV[i+1,j+1]<-as.numeric(z);errors[i+1,j+1]<-attr(z,"error_estimate");unresolved[i+1,j+1]<-isTRUE(attr(z,"diagnostics")$unresolved);upper_bound[i+1,j+1]<-isTRUE(attr(z,"diagnostics")$upper_bound)
    }
    if(!is.null(checkpoint_file)) {
      dir.create(dirname(checkpoint_file),recursive=TRUE,showWarnings=FALSE)
      tmp<-paste0(checkpoint_file,".tmp");saveRDS(list(signature=signature,EV=EV,errors=errors,unresolved=unresolved,upper_bound=upper_bound),tmp)
      if(!file.rename(tmp,checkpoint_file))stop("Cannot save enumeration checkpoint")
    }
    if(progress&&i%%10==0)message("e-values: ",i,"/",n1)
  }
  attr(EV,"error_estimate")<-errors;attr(EV,"unresolved")<-unresolved;attr(EV,"upper_bound")<-upper_bound
  attr(EV,"diagnostics")<-list(method=method,G=G,max_G=max_G,max_error=max(errors),unresolved_count=sum(unresolved),elapsed=unname(proc.time()[3]-started))
  EV
}
# Positive quadrature over the exact independent representation:
# X~Beta(a1,a0), V~Beta(a2,a0+a1), Y=V/(1-X+X*V).
predictive_A <- function(n1,n2,prior,K=NULL,G=64,max_G=1024,abs_tol=1e-8,family="olkin_liu",strict=TRUE,method=c("series","quadrature"),max_terms=1048576) {
  .fbst_par(0,n1,0,n2,prior,family);a0<-prior[1];a1<-prior[2];a2<-prior[3]
  if(family=="independent") {
    v1<-exp(lchoose(n1,0:n1)+lbeta(a1+0:n1,a0+n1-0:n1)-lbeta(a1,a0))
    v2<-exp(lchoose(n2,0:n2)+lbeta(a2+0:n2,a0+n2-0:n2)-lbeta(a2,a0))
    ans<-outer(v1,v2);attr(ans,"normalization_error")<-abs(sum(ans)-1);attr(ans,"error_estimate")<-0;return(ans)
  }
  method<-match.arg(method)
  if(method=="series" && exists("fbst_predictive_series",envir=.fbst_env,inherits=FALSE)) {
    z<-.fbst_env$fbst_predictive_series(n1,n2,as.numeric(prior),min(1e-10,abs_tol/100),max_terms)
    ans<-z$probability;ne<-abs(sum(ans)-1);err<-sum(z$error)
    if(any(!is.finite(ans))||any(ans<0)||ne>abs_tol||z$unresolved_count>0) {
      if(strict)stop(sprintf("Predictive series unresolved: normalization %.3g, L1 error estimate %.3g, %d cells",ne,err,z$unresolved_count))
    }
    attr(ans,"normalization_error")<-ne;attr(ans,"error_estimate")<-err
    attr(ans,"diagnostics")<-list(method="positive_3F2_series_Richardson6",max_terms=z$max_terms_used,L1_error_estimate=err,unresolved=z$unresolved_count>0)
    return(ans)
  }
  old<-NULL
  for(g in unique(pmin(max_G,G*2^(0:10)))) {
    q1<-.fbst_gauss(a1,a0,g);q2<-.fbst_gauss(a2,a0+a1,g)
    B1<-vapply(q1$nodes,function(x)dbinom(0:n1,n1,x),numeric(n1+1))
    conditional<-vapply(q1$nodes,function(x) {
      y<-q2$nodes/(1-x+x*q2$nodes)
      as.numeric(vapply(0:n2,function(j)sum(q2$weights*dbinom(j,n2,y)),numeric(1)))
    },numeric(n2+1))
    ans<-(B1*rep(q1$weights,each=n1+1))%*%t(conditional)
    normalization_error<-abs(sum(ans)-1)
    err<-if(is.null(old))Inf else sum(abs(ans-old))
    if(err<abs_tol)break
    old<-ans;if(g>=max_G)break
  }
  if(normalization_error>1e-9)stop("Alternative predictive failed normalization: ",normalization_error)
  if(strict&&err>=abs_tol)stop(sprintf("Alternative predictive quadrature unresolved: L1 change %.3g at G=%d",err,g))
  attr(ans,"normalization_error")<-normalization_error;attr(ans,"error_estimate")<-err
  attr(ans,"diagnostics")<-list(method="independent_prior_representation",G=g,L1_difference=err,unresolved=err>=abs_tol)
  ans
}
predictive_H <- function(n1,n2,prior,null="marginal",G=128,rel_tol=1e-10) {
  null<-match.arg(null,c("marginal","arc","lor"));a0<-prior[1];a1<-prior[2];a2<-prior[3]
  if(null=="marginal"&&abs(a1-a2)>1e-12)stop("Main marginal null prior requires lambda1=lambda2")
  S<-outer(0:n1,0:n2,"+");N<-n1+n2;LC<-outer(lchoose(n1,0:n1),lchoose(n2,0:n2),"+")
  if(null=="marginal")ans<-exp(LC+lbeta(a1+S,a0+N-S)-lbeta(a1,a0)) else {
    ca<-a1+a2-if(null=="arc")1 else 0;cb<-a0-if(null=="arc")1 else 0
    if(ca<=0||cb<=0)stop("Improper restriction on difference: requires lambda0>1 and lambda1+lambda2>1")
    logI<-function(a,b) {
      q<-.fbst_gauss(a,b,G);lbeta(a,b)+log(sum(q$weights*exp(-sum(prior)*log1p(q$nodes))))
    }
    vals<-vapply(0:N,function(s)logI(ca+s,cb+N-s),numeric(1))-logI(ca,cb)
    ans<-exp(LC+vals[S+1])
  }
  err<-abs(sum(ans)-1)
  if(err>max(1e-9,rel_tol*10))stop("Null predictive normalization failed: ",err)
  attr(ans,"normalization_error")<-err;ans
}
adaptive_cutoff <- function(EV,pH,pA,tol=.002) {
  if(!identical(dim(EV),dim(pH))||!identical(dim(EV),dim(pA))||anyNA(EV)||any(!is.finite(EV))||any(EV<0|EV>1))stop("Incompatible/invalid evidence and predictive matrices")
  if(abs(sum(pH)-1)>1e-7||abs(sum(pA)-1)>1e-7)stop("Predictives must sum to one before calibration")
  ev<-as.numeric(EV);o<-order(ev);v<-ev[o];ends<-which(c(diff(v)!=0,TRUE));ks<-v[ends]
  alpha<-cumsum(as.numeric(pH)[o])[ends];beta<-sum(pA)-cumsum(as.numeric(pA)[o])[ends]
  beta<-pmax(0,beta)
  if(min(ks)>0){ks<-c(0,ks);alpha<-c(0,alpha);beta<-c(sum(pA),beta)}
  # At k=1 the last region is closed on the right; every other interval is [lo,hi).
  hi<-c(ks[-1],1);risk<-alpha+beta;rmin<-min(risk);j<-which.min(risk)
  intervals<-data.frame(lower=ks,upper=hi,lower_closed=TRUE,upper_closed=c(rep(FALSE,length(ks)-1L),TRUE),alpha=alpha,beta=beta,risk=risk)
  mini<-which(risk==rmin);near<-which(risk<=rmin+tol)
  list(kstar=ks[j],alpha=alpha[j],beta=beta[j],risk=rmin,tie_break="smallest attainable cutoff in [0,1] attaining the numerically smallest risk; numerical near-ties listed separately",ks=ks,
       alpha_k=alpha,beta_k=beta,risk_k=risk,intervals=intervals,minimizers=intervals[mini,,drop=FALSE],numerical_ties=intervals[abs(risk-rmin)<=1e-12,,drop=FALSE],near_optimal=intervals[near,,drop=FALSE],
       k_range=range(c(intervals$lower[near],intervals$upper[near])),alpha_range=range(alpha[near]),no_rejection=list(k=-Inf,alpha=0,beta=sum(pA),risk=sum(pA)))
}
rejection_probability <- function(EV,k,th1,th2) sum(outer(dbinom(0:(nrow(EV)-1),nrow(EV)-1,th1),dbinom(0:(ncol(EV)-1),ncol(EV)-1,th2))*(EV<=k))
ztest_power <- function(n,th1,th2,level) {
  se<-sqrt((th1*(1-th1)+th2*(1-th2))/n);d<-th2-th1
  if(se==0)return(as.numeric(d!=0&&level>0))
  z<-qnorm(1-level/2);pnorm(d/se-z)+pnorm(-d/se-z)
}
# CDF of delta by direct integration, using the same Lebesgue reference throughout.
delta_cdf <- function(d,x1,n1,x2,n2,prior,family="olkin_liu",power=1,rel_tol=1e-7,abs_tol=1e-9,max_G=1024) {
  p<-.fbst_par(x1,n1,x2,n2,prior,family,power);norm<-.fbst_normalizer(p,rel_tol/10,max_G)
  one<-function(v) {
    if(v<=-1)return(0);if(v>=1)return(1)
    f<-function(xx)vapply(xx,function(x) {
      h<-min(1,x+v);if(h<=0)return(0)
      .fbst_slice_integral(x,0,h,p,norm$logZ,rel_tol/5,abs_tol/10)
    },numeric(1))
    split<-sort(unique(c(max(0,-v),min(1,1-v),1)));split<-split[split>=max(0,-v)]
    ans<-sum(vapply(seq_len(length(split)-1L),function(j).fbst_integrate(f,split[j],split[j+1],rel_tol,abs_tol),numeric(1)))
    if(ans< -abs_tol||ans>1+1e-6)stop("Delta CDF violates probability bounds")
    min(1,max(0,ans))
  }
  vapply(d,one,numeric(1))
}
posterior_summary <- function(x1,n1,x2,n2,prior,G=128,family="olkin_liu",power=1,method=c("adaptive","grid"),rel_tol=1e-7,abs_tol=1e-9,max_G=1024) {
  method<-match.arg(method);p<-.fbst_par(x1,n1,x2,n2,prior,family,power)
  q<-.fbst_grid(p,G);d<-outer(q$x,q$y,function(x,y)y-x);mu<-sum(q$w*d)
  q2<-.fbst_grid(p,min(max_G,2*G));mu2<-sum(q2$w*outer(q2$x,q2$y,function(x,y)y-x))
  norm<-.fbst_normalizer(p,rel_tol/10,max_G)
  pp<-p;pp[1]<-pp[1]+1;m1<-exp(.fbst_normalizer(pp,rel_tol/10,max_G)$logZ-norm$logZ)
  pp<-p;pp[3]<-pp[3]+1;m2<-exp(.fbst_normalizer(pp,rel_tol/10,max_G)$logZ-norm$logZ)
  mu2<-m2-m1
  mean_error<-norm$error*3
  if(method=="grid") {
    describe<-function(q) {
      dd<-outer(q$x,q$y,function(x,y)y-x);o<-order(dd);cum<-cumsum(q$w[o])
      c(lo=dd[o][which(cum>=.025)[1]],hi=dd[o][which(cum>=.975)[1]],
        P_increase=sum(q$w*(dd>0))+.5*sum(q$w*(dd==0)),P_decrease=sum(q$w*(dd<0))+.5*sum(q$w*(dd==0)))
    }
    coarse<-describe(q);fine<-describe(q2)
    component_errors<-c(mean=max(abs(mu-mu2),mean_error),abs(fine-coarse))
    return(list(mean=mu2,lo=unname(fine["lo"]),hi=unname(fine["hi"]),P_increase=unname(fine["P_increase"]),P_decrease=unname(fine["P_decrease"]),
      mean_error_estimate=unname(component_errors["mean"]),error_estimate=max(component_errors),component_error_estimates=component_errors,
      diagnostics=list(not_certified=TRUE,coarse_G=G,fine_G=min(max_G,2*G),error_method="successive quadratures; not a mathematical error bound"),method="grid"))
  }
  cdf<-function(d)delta_cdf(d,x1,n1,x2,n2,prior,family,power,rel_tol,abs_tol,max_G)
  qq<-function(v)uniroot(function(d)cdf(d)-v,c(-1,1),tol=max(1e-8,abs_tol))$root
  list(mean=mu2,lo=qq(.025),hi=qq(.975),P_increase=delta_cdf(0,x2,n2,x1,n1,prior[c(1,3,2)],family,power,rel_tol,abs_tol,max_G),P_decrease=cdf(0),error_estimate=max(mean_error,rel_tol),method="adaptive")
}
delta_density <- function(x1,n1,x2,n2,prior,G=96,grid=seq(-1,1,length.out=4001),family="olkin_liu",power=1,rel_tol=1e-7,abs_tol=1e-9,max_G=1024) {
  p<-.fbst_par(x1,n1,x2,n2,prior,family,power);norm<-.fbst_normalizer(p,rel_tol/10,max_G)
  val<-vapply(grid,function(d) {
    lo<-max(0,-d);hi<-min(1,1-d);if(lo>=hi)return(0)
    .fbst_integrate(function(x)exp(.fbst_logkernel(x,x+d,p)-norm$logZ),lo,hi,rel_tol,abs_tol)
  },numeric(1))
  area<-sum(diff(grid)*(head(val,-1)+tail(val,-1))/2)
  if(abs(area-1)>1e-3)stop("Delta density grid failed normalization; refine grid: ",area)
  ans<-data.frame(delta=grid,density=val);attr(ans,"normalization_error")<-abs(area-1);ans
}
contrast <- function(hk,hc) {
  if(!is.data.frame(hk)||!is.data.frame(hc)||!identical(hk$delta,hc$delta))stop("Use delta_density data frames on identical grids")
  dx<-diff(hk$delta)[1];if(max(abs(diff(hk$delta)-dx))>1e-10)stop("Uniform delta grid required")
  wk<-hk$density*dx;wc<-hc$density*dx;wk[c(1,length(wk))]<-wk[c(1,length(wk))]/2;wc[c(1,length(wc))]<-wc[c(1,length(wc))]/2
  nk<-sum(wk);nc<-sum(wc)
  if(max(abs(c(nk,nc)-1))>1e-3)stop("Contrast densities not normalized")
  wk<-wk/nk;wc<-wc/nc;m<-length(wk);N<-2^ceiling(log2(2*m-1))
  mass<-Re(fft(fft(c(wk,rep(0,N-m)))*fft(c(rev(wc),rep(0,N-m))),inverse=TRUE))/N
  mass<-pmax(0,mass[seq_len(2*m-1)]);mass<-mass/sum(mass);x<-(seq_along(mass)-m)*dx;cum<-cumsum(mass)
  qq<-function(v)approx(c(0,cum-.5*mass,1),c(x[1]-dx/2,x,tail(x,1)+dx/2),xout=v,ties="ordered")$y;pg<-sum(mass[x>dx/2])+.5*sum(mass[abs(x)<=dx/2])
  list(mean=sum(x*mass),lo=qq(.025),hi=qq(.975),P_greater=pg,P_positive=pg,P_negative=sum(mass[x< -dx/2])+.5*sum(mass[abs(x)<=dx/2]),error_estimate=max(abs(c(nk,nc)-1),dx),distribution=data.frame(delta=x,mass=mass))
}

posterior_grid <- function(x1,n1,x2,n2,prior,G=128,power=1,family="olkin_liu") {
  q<-.fbst_grid(.fbst_par(x1,n1,x2,n2,prior,family,power),G)
  list(theta1=q$x,theta2=q$y,weights=q$w,log_normalizer=q$logZ)
}
log_posterior_kernel <- function(t1,t2,x1,n1,x2,n2,prior,power=1,family="olkin_liu") .fbst_logkernel(t1,t2,.fbst_par(x1,n1,x2,n2,prior,family,power))
