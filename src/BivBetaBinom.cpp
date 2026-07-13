// BivBetaBinom.cpp
// C++ core for the bivariate beta-binomial model: closed-form posterior
// density, FBST e-values by 2D quadrature, adaptive cutoff k*, and samplers
// for the prior-based and posterior-based error formulations.
// Compile from R with: Rcpp::sourceCpp("src/BivBetaBinom.cpp")

#include <Rcpp.h>
using namespace Rcpp;

// log pFq(U; L; 1) by power series; converges when sum(L) > sum(U).
// [[Rcpp::export]]
double log_genhypergeo_z1(NumericVector U, NumericVector L,
                          int max_iter = 200000, double tol = 1e-16) {
  double log_term = 0.0;
  double sum = 1.0;
  for (int k = 1; k <= max_iter; ++k) {
    double log_ratio = 0.0;
    for (int i = 0; i < U.size(); ++i) log_ratio += std::log(U[i] + k - 1);
    for (int j = 0; j < L.size(); ++j) log_ratio -= std::log(L[j] + k - 1);
    log_ratio -= std::log((double)k);
    log_term += log_ratio;
    double term = std::exp(log_term);
    sum += term;
    if (term < tol * sum) break;
  }
  return std::log(sum);
}

// Precompute every posterior constant that does not depend on (theta1, theta2).
// [[Rcpp::export]]
List bb_constants(int n1, int n2, int x1, int x2,
                  double a0, double a1, double a2) {
  double a = a0 + a1 + a2;
  NumericVector U = NumericVector::create(a, x1 + a1, x2 + a2);
  NumericVector L = NumericVector::create(n1 + a, n2 + a);

  double log_3F2 = log_genhypergeo_z1(U, L);

  double log_C = R::lgammafn(n1 + a) + R::lgammafn(n2 + a)
               - R::lgammafn(x1 + a1) - R::lgammafn(n1 - x1 + a - a1)
               - R::lgammafn(x2 + a2) - R::lgammafn(n2 - x2 + a - a2)
               - log_3F2;

  double e1 = a1 + x1 - 1;              // exponent of theta1
  double f1 = a2 + a0 + (n1 - x1) - 1;  // exponent of (1 - theta1)
  double e2 = a2 + x2 - 1;              // exponent of theta2
  double f2 = a1 + a0 + (n2 - x2) - 1;  // exponent of (1 - theta2)

  return List::create(
    _["log_C"]  = log_C,
    _["e1"] = e1, _["f1"] = f1,
    _["e2"] = e2, _["f2"] = f2,
    _["a"]      = a,
    _["log_3F2"]= log_3F2
  );
}

// Posterior density at a single point.
// [[Rcpp::export]]
double densBB_cpp(double theta1, double theta2, List consts) {
  double log_C = consts["log_C"];
  double e1 = consts["e1"], f1 = consts["f1"];
  double e2 = consts["e2"], f2 = consts["f2"];
  double a  = consts["a"];
  double log_f = log_C
    + e1 * std::log(theta1) + f1 * std::log1p(-theta1)
    + e2 * std::log(theta2) + f2 * std::log1p(-theta2)
    - a  * std::log1p(-theta1 * theta2);
  return std::exp(log_f);
}

// Vectorised posterior density (theta1 and theta2 of equal length).
// [[Rcpp::export]]
NumericVector densBB_vec(NumericVector theta1, NumericVector theta2, List consts) {
  double log_C = consts["log_C"];
  double e1 = consts["e1"], f1 = consts["f1"];
  double e2 = consts["e2"], f2 = consts["f2"];
  double a  = consts["a"];
  int n = theta1.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double t1 = theta1[i], t2 = theta2[i];
    out[i] = std::exp(log_C
      + e1 * std::log(t1) + f1 * std::log1p(-t1)
      + e2 * std::log(t2) + f2 * std::log1p(-t2)
      - a  * std::log1p(-t1 * t2));
  }
  return out;
}

// Posterior density on a grid: returns matrix (length(x) x length(y)).
// [[Rcpp::export]]
NumericMatrix densBB_grid(NumericVector x, NumericVector y, List consts) {
  double log_C = consts["log_C"];
  double e1 = consts["e1"], f1 = consts["f1"];
  double e2 = consts["e2"], f2 = consts["f2"];
  double a  = consts["a"];
  int nx = x.size(), ny = y.size();
  NumericVector lx(nx), ly(ny);
  for (int i = 0; i < nx; ++i)
    lx[i] = e1 * std::log(x[i]) + f1 * std::log1p(-x[i]);
  for (int j = 0; j < ny; ++j)
    ly[j] = e2 * std::log(y[j]) + f2 * std::log1p(-y[j]);

  NumericMatrix Z(nx, ny);
  for (int j = 0; j < ny; ++j) {
    double yj = y[j];
    for (int i = 0; i < nx; ++i) {
      Z(i, j) = std::exp(log_C + lx[i] + ly[j]
                         - a * std::log1p(-x[i] * yj));
    }
  }
  return Z;
}

// Posterior density restricted to H: theta1 == theta2 == theta.
// [[Rcpp::export]]
double densBB_H_cpp(double theta, List consts) {
  double log_C = consts["log_C"];
  double e1 = consts["e1"], f1 = consts["f1"];
  double e2 = consts["e2"], f2 = consts["f2"];
  double a  = consts["a"];
  return std::exp(log_C
    + (e1 + e2) * std::log(theta)
    + (f1 + f2) * std::log1p(-theta)
    - a * std::log1p(-theta * theta));
}

// Vectorised version of densBB_H_cpp.
// [[Rcpp::export]]
NumericVector densBB_H_vec(NumericVector theta, List consts) {
  double log_C = consts["log_C"];
  double e1 = consts["e1"], f1 = consts["f1"];
  double e2 = consts["e2"], f2 = consts["f2"];
  double a  = consts["a"];
  int n = theta.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double t = theta[i];
    out[i] = std::exp(log_C
      + (e1 + e2) * std::log(t)
      + (f1 + f2) * std::log1p(-t)
      - a * std::log1p(-t * t));
  }
  return out;
}

// FBST e-value from posterior draws: share of draws whose density falls
// below the supremum of the posterior on H.
// [[Rcpp::export]]
double ev_FBST(NumericVector theta1_post, NumericVector theta2_post,
               double sup_H, List consts) {
  double log_C = consts["log_C"];
  double e1 = consts["e1"], f1 = consts["f1"];
  double e2 = consts["e2"], f2 = consts["f2"];
  double a  = consts["a"];
  int n = theta1_post.size();
  long count = 0;
  double log_supH = std::log(sup_H);
  for (int i = 0; i < n; ++i) {
    double t1 = theta1_post[i], t2 = theta2_post[i];
    double lf = log_C
      + e1 * std::log(t1) + f1 * std::log1p(-t1)
      + e2 * std::log(t2) + f2 * std::log1p(-t2)
      - a  * std::log1p(-t1 * t2);
    if (lf < log_supH) ++count;
  }
  return (double)count / (double)n;
}

// Supremum of the posterior density on H by a grid sweep over (eps, 1-eps);
// the restricted density is unimodal in typical cases.
// [[Rcpp::export]]
List find_sup_H(List consts, int ngrid = 2001, double eps = 1e-6) {
  double log_C = consts["log_C"];
  double e1 = consts["e1"], f1 = consts["f1"];
  double e2 = consts["e2"], f2 = consts["f2"];
  double a  = consts["a"];

  double h = (1.0 - 2.0 * eps) / (ngrid - 1);
  double best_lf = -INFINITY;
  double best_t  = 0.5;
  for (int i = 0; i < ngrid; ++i) {
    double t = eps + i * h;
    double lf = log_C
      + (e1 + e2) * std::log(t)
      + (f1 + f2) * std::log1p(-t)
      - a * std::log1p(-t * t);
    if (lf > best_lf) { best_lf = lf; best_t = t; }
  }
  return List::create(_["sup_H"] = std::exp(best_lf),
                      _["theta"] = best_t,
                      _["log_sup_H"] = best_lf);
}

// FBST e-value by composite-Simpson 2D quadrature over (eps, 1-eps)^2 on the
// closed-form posterior (no MCMC). Normalising by the total integral absorbs
// the eps truncation and quadrature error, giving P(f(theta|x) <= sup_H | x).
// [[Rcpp::export]]
double ev_quad(List consts, double sup_H, int ngrid = 401, double eps = 1e-6) {
  double log_C = consts["log_C"];
  double e1 = consts["e1"], f1 = consts["f1"];
  double e2 = consts["e2"], f2 = consts["f2"];
  double a  = consts["a"];
  double log_supH = std::log(sup_H);

  if (ngrid % 2 == 0) ++ngrid;  // Simpson needs an odd number of nodes

  double h = (1.0 - 2.0 * eps) / (ngrid - 1);

  std::vector<double> theta(ngrid), lx(ngrid), ly(ngrid), w(ngrid);
  for (int i = 0; i < ngrid; ++i) {
    theta[i] = eps + i * h;
    lx[i] = e1 * std::log(theta[i]) + f1 * std::log1p(-theta[i]);
    ly[i] = e2 * std::log(theta[i]) + f2 * std::log1p(-theta[i]);
    if (i == 0 || i == ngrid - 1) w[i] = 1.0;
    else if (i % 2 == 1)          w[i] = 4.0;
    else                           w[i] = 2.0;
  }
  double scale = (h / 3.0) * (h / 3.0);

  double integral_below = 0.0;
  double integral_total = 0.0;
  for (int i = 0; i < ngrid; ++i) {
    double wi = w[i];
    double ti = theta[i];
    for (int j = 0; j < ngrid; ++j) {
      double lf = log_C + lx[i] + ly[j] - a * std::log1p(-ti * theta[j]);
      double f  = std::exp(lf);
      double wij = wi * w[j];
      integral_total += wij * f;
      if (lf <= log_supH) integral_below += wij * f;
    }
  }
  // `scale` cancels in the ratio; kept for clarity.
  integral_total *= scale;
  integral_below *= scale;
  return integral_below / integral_total;
}

// Convenience wrapper: constants + sup_H + e-value in one call.
// [[Rcpp::export]]
double ev_quad_from_data(int n1, int n2, int x1, int x2,
                         double a0, double a1, double a2,
                         int ngrid_quad = 401, int ngrid_sup = 2001,
                         double eps = 1e-6) {
  List consts = bb_constants(n1, n2, x1, x2, a0, a1, a2);
  List sh = find_sup_H(consts, ngrid_sup, eps);
  double sup_H = sh["sup_H"];
  return ev_quad(consts, sup_H, ngrid_quad, eps);
}

// Tabulated CDF of the prior restricted to H,
// f_H(t) proportional to t^(a1+a2-2) (1-t)^(a0-2) / (1+t)^a.
// Also returns the line-integral normalising constant Z. Note: the integral
// diverges as t -> 1 when a0 <= 1; the eps truncation regularises it.
// [[Rcpp::export]]
List fH_cdf(double a0, double a1, double a2, int ngrid = 10001, double eps = 1e-6) {
  double a = a0 + a1 + a2;
  std::vector<double> tg(ngrid), pdf(ngrid), cdf(ngrid);
  double h = (1.0 - 2.0 * eps) / (ngrid - 1);
  for (int i = 0; i < ngrid; ++i) {
    tg[i] = eps + i * h;
    double lp = (a1 + a2 - 2.0) * std::log(tg[i])
              + (a0 - 2.0)      * std::log1p(-tg[i])
              - a               * std::log1p(tg[i]);
    pdf[i] = std::exp(lp);
  }
  cdf[0] = 0.0;
  for (int i = 1; i < ngrid; ++i)
    cdf[i] = cdf[i-1] + 0.5 * h * (pdf[i] + pdf[i-1]);
  double Iraw = cdf[ngrid-1];
  for (int i = 0; i < ngrid; ++i) cdf[i] /= Iraw;
  double logK = R::lgammafn(a) - R::lgammafn(a0) - R::lgammafn(a1) - R::lgammafn(a2);
  double Z = std::sqrt(2.0) * std::exp(logK) * Iraw;
  return List::create(_["t"] = NumericVector(tg.begin(), tg.end()),
                      _["cdf"] = NumericVector(cdf.begin(), cdf.end()),
                      _["Z"] = Z);
}

// Sample from f_H by inverse-CDF with linear interpolation.
// [[Rcpp::export]]
NumericVector sample_fH(int n, List fH_table) {
  NumericVector tg = fH_table["t"];
  NumericVector cdf = fH_table["cdf"];
  int N = tg.size();
  NumericVector out(n);
  for (int k = 0; k < n; ++k) {
    double u = R::unif_rand();
    int lo = 0, hi = N - 1;
    while (lo < hi) {
      int mid = (lo + hi) / 2;
      if (cdf[mid] < u) lo = mid + 1;
      else              hi = mid;
    }
    if (lo == 0) out[k] = tg[0];
    else {
      double denom = cdf[lo] - cdf[lo-1];
      double frac  = (denom > 0) ? (u - cdf[lo-1]) / denom : 0.5;
      out[k] = tg[lo-1] + frac * (tg[lo] - tg[lo-1]);
    }
  }
  return out;
}

// Sample from the Olkin-Liu bivariate beta prior via its gamma construction:
// V0 ~ Gamma(a0), Vj ~ Gamma(aj) independent; theta_j = Vj / (Vj + V0).
// [[Rcpp::export]]
NumericMatrix sample_prior(int n, double a0, double a1, double a2) {
  NumericMatrix out(n, 2);
  for (int i = 0; i < n; ++i) {
    double v1 = R::rgamma(a1, 1.0);
    double v2 = R::rgamma(a2, 1.0);
    double v0 = R::rgamma(a0, 1.0);
    out(i, 0) = v1 / (v1 + v0);
    out(i, 1) = v2 / (v2 + v0);
  }
  return out;
}

// Prior-based e-value distribution under H: theta ~ f_H, X ~ Bin x Bin(theta),
// e-value by quadrature. Returns the M simulated e-values.
// [[Rcpp::export]]
NumericVector simulate_evs_H(int n1, int n2,
                             double a0, double a1, double a2,
                             int M, int ngrid_quad = 401, int ngrid_sup = 2001,
                             double eps = 1e-6) {
  List fH_table = fH_cdf(a0, a1, a2, 10001, eps);
  NumericVector thetas = sample_fH(M, fH_table);
  NumericVector evs(M);
  for (int i = 0; i < M; ++i) {
    double t = thetas[i];
    int x1 = (int) R::rbinom(n1, t);
    int x2 = (int) R::rbinom(n2, t);
    evs[i] = ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2,
                                ngrid_quad, ngrid_sup, eps);
  }
  return evs;
}

// Prior-based e-value distribution under the full prior (the alternative,
// up to a null set).
// [[Rcpp::export]]
NumericVector simulate_evs_A(int n1, int n2,
                             double a0, double a1, double a2,
                             int M, int ngrid_quad = 401, int ngrid_sup = 2001,
                             double eps = 1e-6) {
  NumericMatrix thetas = sample_prior(M, a0, a1, a2);
  NumericVector evs(M);
  for (int i = 0; i < M; ++i) {
    double t1 = thetas(i, 0), t2 = thetas(i, 1);
    int x1 = (int) R::rbinom(n1, t1);
    int x2 = (int) R::rbinom(n2, t2);
    evs[i] = ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2,
                                ngrid_quad, ngrid_sup, eps);
  }
  return evs;
}

// Averaged error curves on a k grid:
// alpha(k) = mean(ev_H <= k), beta(k) = mean(ev_A > k).
// [[Rcpp::export]]
DataFrame error_curves(NumericVector ev_H, NumericVector ev_A, NumericVector k_grid,
                       double a_w = 1.0, double b_w = 1.0) {
  int K = k_grid.size();
  int MH = ev_H.size(), MA = ev_A.size();
  std::vector<double> sH(ev_H.begin(), ev_H.end());
  std::vector<double> sA(ev_A.begin(), ev_A.end());
  std::sort(sH.begin(), sH.end());
  std::sort(sA.begin(), sA.end());
  NumericVector alpha(K), beta(K), sum(K);
  for (int j = 0; j < K; ++j) {
    double k = k_grid[j];
    int cnt_a = std::upper_bound(sH.begin(), sH.end(), k) - sH.begin();
    int cnt_b = MA - (std::upper_bound(sA.begin(), sA.end(), k) - sA.begin());
    alpha[j] = (double)cnt_a / MH;
    beta[j]  = (double)cnt_b / MA;
    sum[j]   = a_w * alpha[j] + b_w * beta[j];
  }
  return DataFrame::create(_["k"] = k_grid,
                           _["alpha"] = alpha,
                           _["beta"]  = beta,
                           _["sum"]   = sum);
}

// Adaptive cutoff k* = argmin_k a*alpha(k) + b*beta(k). Since alpha and beta
// are step functions, the optimum is attained at one of the simulated
// e-values, so all candidates are enumerated.
// [[Rcpp::export]]
List find_kstar(NumericVector ev_H, NumericVector ev_A,
                double a_w = 1.0, double b_w = 1.0) {
  std::vector<double> sH(ev_H.begin(), ev_H.end());
  std::vector<double> sA(ev_A.begin(), ev_A.end());
  std::sort(sH.begin(), sH.end());
  std::sort(sA.begin(), sA.end());
  int MH = sH.size(), MA = sA.size();

  std::vector<double> ks;
  ks.reserve(MH + MA + 2);
  ks.push_back(0.0);
  ks.push_back(1.0);
  for (auto v : sH) ks.push_back(v);
  for (auto v : sA) ks.push_back(v);
  std::sort(ks.begin(), ks.end());
  ks.erase(std::unique(ks.begin(), ks.end()), ks.end());

  double best_k = 0.0, best_obj = std::numeric_limits<double>::infinity();
  double best_a = 0.0, best_b = 0.0;
  for (size_t j = 0; j < ks.size(); ++j) {
    double k = ks[j];
    int cnt_a = std::upper_bound(sH.begin(), sH.end(), k) - sH.begin();
    int cnt_b = MA - (std::upper_bound(sA.begin(), sA.end(), k) - sA.begin());
    double alpha = (double)cnt_a / MH;
    double beta  = (double)cnt_b / MA;
    double obj   = a_w * alpha + b_w * beta;
    if (obj < best_obj) {
      best_obj = obj; best_k = k; best_a = alpha; best_b = beta;
    }
  }
  return List::create(_["k_star"] = best_k,
                      _["alpha"]  = best_a,
                      _["beta"]   = best_b,
                      _["obj"]    = best_obj);
}

// Tabulated CDF of the posterior restricted to H given the observed counts:
// f_H(t | x_obs) proportional to t^E (1-t)^F / (1+t)^a, with
// E = a1+a2+x1+x2-2 and F = a0+n1+n2-x1-x2-2 (log-sum-exp normalisation).
// [[Rcpp::export]]
List fH_post_cdf(int n1, int n2, int x1, int x2,
                 double a0, double a1, double a2,
                 int ngrid = 10001, double eps = 1e-6) {
  double a = a0 + a1 + a2;
  double E = a1 + a2 + x1 + x2 - 2.0;
  double F = a0 + n1 + n2 - x1 - x2 - 2.0;

  std::vector<double> tg(ngrid), lp(ngrid), pdf(ngrid), cdf(ngrid);
  double h = (1.0 - 2.0 * eps) / (ngrid - 1);
  double max_lp = -INFINITY;
  for (int i = 0; i < ngrid; ++i) {
    tg[i] = eps + i * h;
    lp[i] = E * std::log(tg[i]) + F * std::log1p(-tg[i]) - a * std::log1p(tg[i]);
    if (lp[i] > max_lp) max_lp = lp[i];
  }
  for (int i = 0; i < ngrid; ++i) pdf[i] = std::exp(lp[i] - max_lp);
  cdf[0] = 0.0;
  for (int i = 1; i < ngrid; ++i)
    cdf[i] = cdf[i-1] + 0.5 * h * (pdf[i] + pdf[i-1]);
  double Iraw = cdf[ngrid-1];
  for (int i = 0; i < ngrid; ++i) cdf[i] /= Iraw;

  return List::create(_["t"]   = NumericVector(tg.begin(),  tg.end()),
                      _["cdf"] = NumericVector(cdf.begin(), cdf.end()));
}

// Sample the bivariate posterior by sampling-importance-resampling (SIR)
// with the bivariate beta prior as proposal. Log-scale weights with
// log-sum-exp stabilisation; safe for x = 0 or x = n.
// [[Rcpp::export]]
NumericMatrix sample_posterior(int M, int n1, int n2, int x1, int x2,
                                double a0, double a1, double a2,
                                int N_prop = -1) {
  if (N_prop <= 0) N_prop = std::max(50 * M, 10000);
  NumericMatrix props = sample_prior(N_prop, a0, a1, a2);

  std::vector<double> log_w(N_prop);
  double max_lw = -INFINITY;
  for (int i = 0; i < N_prop; ++i) {
    double t1 = props(i, 0), t2 = props(i, 1);
    double lw = 0.0;
    if (x1 > 0)        lw += x1        * std::log(t1);
    if (n1 - x1 > 0)   lw += (n1 - x1) * std::log1p(-t1);
    if (x2 > 0)        lw += x2        * std::log(t2);
    if (n2 - x2 > 0)   lw += (n2 - x2) * std::log1p(-t2);
    log_w[i] = lw;
    if (lw > max_lw) max_lw = lw;
  }
  std::vector<double> w(N_prop), cw(N_prop);
  double sw = 0.0;
  for (int i = 0; i < N_prop; ++i) { w[i] = std::exp(log_w[i] - max_lw); sw += w[i]; }
  cw[0] = w[0] / sw;
  for (int i = 1; i < N_prop; ++i) cw[i] = cw[i-1] + w[i] / sw;

  NumericMatrix out(M, 2);
  for (int k = 0; k < M; ++k) {
    double u = R::unif_rand();
    int lo = 0, hi = N_prop - 1;
    while (lo < hi) {
      int mid = (lo + hi) / 2;
      if (cw[mid] < u) lo = mid + 1;
      else              hi = mid;
    }
    out(k, 0) = props(lo, 0);
    out(k, 1) = props(lo, 1);
  }
  return out;
}

// Effective sample size of the SIR weights (diagnostic); a low value means
// N_prop is insufficient for that posterior.
// [[Rcpp::export]]
double sir_ess(int n1, int n2, int x1, int x2,
               double a0, double a1, double a2, int N_prop) {
  NumericMatrix props = sample_prior(N_prop, a0, a1, a2);
  std::vector<double> log_w(N_prop);
  double max_lw = -INFINITY;
  for (int i = 0; i < N_prop; ++i) {
    double t1 = props(i, 0), t2 = props(i, 1);
    double lw = 0.0;
    if (x1 > 0)        lw += x1        * std::log(t1);
    if (n1 - x1 > 0)   lw += (n1 - x1) * std::log1p(-t1);
    if (x2 > 0)        lw += x2        * std::log(t2);
    if (n2 - x2 > 0)   lw += (n2 - x2) * std::log1p(-t2);
    log_w[i] = lw;
    if (lw > max_lw) max_lw = lw;
  }
  double s = 0.0, s2 = 0.0;
  for (int i = 0; i < N_prop; ++i) {
    double w = std::exp(log_w[i] - max_lw);
    s  += w;
    s2 += w * w;
  }
  return (s * s) / s2;
}

// Posterior-based e-value distribution under H: theta ~ f_H(. | x_obs),
// X ~ Bin x Bin(theta), e-value by quadrature.
// [[Rcpp::export]]
NumericVector simulate_evs_H_post(int n1, int n2, int x1_obs, int x2_obs,
                                   double a0, double a1, double a2,
                                   int M, int ngrid_quad = 401, int ngrid_sup = 2001,
                                   double eps = 1e-6) {
  List fH_table = fH_post_cdf(n1, n2, x1_obs, x2_obs, a0, a1, a2, 10001, eps);
  NumericVector thetas = sample_fH(M, fH_table);
  NumericVector evs(M);
  for (int i = 0; i < M; ++i) {
    double t = thetas[i];
    int x1 = (int) R::rbinom(n1, t);
    int x2 = (int) R::rbinom(n2, t);
    evs[i] = ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2,
                                ngrid_quad, ngrid_sup, eps);
  }
  return evs;
}

// Posterior-based e-value distribution under the alternative:
// theta ~ full posterior given x_obs, X ~ Bin x Bin(theta).
// [[Rcpp::export]]
NumericVector simulate_evs_A_post(int n1, int n2, int x1_obs, int x2_obs,
                                   double a0, double a1, double a2,
                                   int M, int ngrid_quad = 401, int ngrid_sup = 2001,
                                   double eps = 1e-6, int N_prop = -1) {
  NumericMatrix thetas = sample_posterior(M, n1, n2, x1_obs, x2_obs, a0, a1, a2, N_prop);
  NumericVector evs(M);
  for (int i = 0; i < M; ++i) {
    double t1 = thetas(i, 0), t2 = thetas(i, 1);
    int x1 = (int) R::rbinom(n1, t1);
    int x2 = (int) R::rbinom(n2, t2);
    evs[i] = ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2,
                                ngrid_quad, ngrid_sup, eps);
  }
  return evs;
}
