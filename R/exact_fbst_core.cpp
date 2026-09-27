#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double fbst_grid_ev(NumericVector x, NumericVector y, NumericMatrix w, NumericVector p, NumericVector r, double target) {
  double out=0; int n=x.size(), m=y.size();
  for(int i=0;i<n;i++)for(int j=0;j<m;j++) {
    double xx=x[i],yy=y[j],one=1-xx*yy;
    double L=r[0]*std::log(xx)+r[1]*std::log1p(-xx)+r[2]*std::log(yy)+r[3]*std::log1p(-yy)-r[4]*std::log(one);
    double dx=(i==n-1?1:.5*(x[i]+x[i+1]))-(i==0?0:.5*(x[i-1]+x[i]));
    double dy=(j==m-1?1:.5*(y[j]+y[j+1]))-(j==0?0:.5*(y[j-1]+y[j]));
    double gx=r[0]/xx-r[1]/(1-xx)+r[4]*yy/one;
    double gy=r[2]/yy-r[3]/(1-yy)+r[4]*xx/one;
    double span=std::abs(gx)*dx+std::abs(gy)*dy;
    double frac=span>0?std::min(1.,std::max(0.,.5+(target-L)/span)):(L<=target?1.:0.);
    out+=w(i,j)*frac;
  }
  return out;
}
// [[Rcpp::export]]
Rcpp::List fbst_predictive_series(int n1,int n2,NumericVector prior,double tol=1e-10,int max_terms=131072) {
  double a0=prior[0],a1=prior[1],a2=prior[2],lam=a0+a1+a2;
  NumericMatrix out(n1+1,n2+1),err(n1+1,n2+1);int n_unresolved=0,largest=0;
  for(int i=0;i<=n1;i++)for(int j=0;j<=n2;j++) {
    if(n1==n2&&a1==a2&&j<i){out(i,j)=out(j,i);err(i,j)=err(j,i);continue;}
    double lp=R::lchoose(n1,i)+R::lchoose(n2,j)+R::lgammafn(lam)-R::lgammafn(a0)-R::lgammafn(a1)-R::lgammafn(a2)+R::lbeta(i+a1,n1-i+a0+a2)+R::lbeta(j+a2,n2-j+a0+a1);
    long double term=std::exp((long double)lp),sum=term,prev=0,best=sum,diff=INFINITY;
    double s=a0+n1-i+n2-j;std::vector<long double> old; int next=32,used=0;bool done=false;
    for(int k=1;k<=max_terms;k++) {
      long double z=k-1;
      term*=((lam+z)*(i+a1+z)/(n1+lam+z))*((j+a2+z)/((n2+lam+z)*k));sum+=term;
      if(k==next) {
        std::vector<long double> now;now.push_back(sum);
        int depth=std::min((int)old.size(),6);
        for(int d=1;d<=depth;d++) {
          long double inv=std::exp2(-(long double)(s+d-1));
          now.push_back((now[d-1]-inv*old[d-1])/(1-inv));
        }
        best=now.back();diff=std::abs(best-prev);
        if(old.size()>=6 && diff<=tol*std::abs(best) && best>=0){done=true;used=k;break;}
        if(term==0){best=sum;diff=0;done=true;used=k;break;}
        old=now;prev=best;next*=2;
      }
      used=k;
    }
    if(!done)n_unresolved++;
    largest=std::max(largest,used);out(i,j)=(double)best;err(i,j)=(double)diff;
  }
  return List::create(_["probability"]=out,_["error"]=err,_["unresolved_count"]=n_unresolved,_["max_terms_used"]=largest);
}
// [[Rcpp::export]]
Rcpp::List fbst_normalizer_series(NumericVector p,double tol=1e-11,int max_terms=1048576) {
  double A1=p[0],B1=p[1],A2=p[2],B2=p[3],lam=p[4],s=B1+B2-lam;
  if(s<=0)stop("Non-integrable posterior");
  long double lp=R::lbeta(A1,B1)+R::lbeta(A2,B2);
  long double term=1,sum=term,prev=0,best=sum,diff=INFINITY;
  std::vector<long double> old;int next=32,used=0;bool done=false;
  for(int k=1;k<=max_terms;k++) {
    long double z=k-1;
    term*=((lam+z)*(A1+z)/(A1+B1+z))*((A2+z)/((A2+B2+z)*k));sum+=term;
    if(k==next) {
      std::vector<long double> now;now.push_back(sum);int depth=std::min((int)old.size(),6);
      for(int d=1;d<=depth;d++) {long double inv=std::exp2(-(long double)(s+d-1));now.push_back((now[d-1]-inv*old[d-1])/(1-inv));}
      best=now.back();diff=std::abs(best-prev);
      if(old.size()>=6&&diff<=tol*std::abs(best)&&best>0){done=true;used=k;break;}
      if(term==0){best=sum;diff=0;done=true;used=k;break;}
      old=now;prev=best;next*=2;
    }used=k;
  }
  return List::create(_["logZ"]=(double)(lp+std::log(best)),_["log_lower"]=(double)(lp+std::log(sum))+std::log(1-1e-8),_["error"]=(double)(diff/best),_["G"]=0,_["max_terms_used"]=used,_["unresolved"]=!done);
}
