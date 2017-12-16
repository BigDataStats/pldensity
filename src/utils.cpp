#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 


using namespace Rcpp;
using namespace arma;
using namespace std;

// // This inline function creates one sample of multivariate t-distribution
// // [[Rcpp::export]]
// inline arma::vec rSt(
//     const arma::vec& mu,
//     const arma::mat& Sigma,
//     const double df) {
//   
//   int d = mu.size();
//   vec normal =  chol(Sigma).t() * randn(d);
//   double chi2 = rchisq(1, df)[0];
//   return mu + normal / sqrt(chi2 / df);
// }

double dst(
    const arma::vec& x,
    const arma::vec& mu,
    const arma::mat& Sigma, // Sigma^{-1}
    const double df) {
  int d = mu.n_elem;
  vec xcentered = x - mu;
  mat InvSigmaReg = inv_sympd(Sigma + 1e-24 * eye<mat>(d ,d));
  double ldet;
  double sign;
  log_det(ldet, sign, InvSigmaReg);
  double innerterm, extterm;
  if (df < 35) {
    innerterm =  - 0.5 * (df + d) * log(1.0 + as_scalar(xcentered.t() * InvSigmaReg * xcentered) / df);
    // compute and store the logarithm of determinant
    extterm = log(tgamma(0.5 * (df + d))) - log(tgamma(0.5 * df)) - 0.5 * d * log(df * M_PI) + 0.5 * ldet;
  } else {
    innerterm =  -0.5 * as_scalar(xcentered.t() * InvSigmaReg * xcentered);
    // compute and store the logarithm of determinant
    extterm = - (d / 2.0) * log(2 * M_PI) + 0.5 * ldet;
  }
  if (isnan(extterm)) {
    Rcout << "Error: numeric problems in density estimation" << endl;
    throw "Error: numeric problems in density estimation";
  }
  return exp(extterm + innerterm);
}

// Samples from a multivariate categorical variable
arma::uvec resample(int N, arma::vec prob) {
  vec probsum = cumsum(prob) / sum(prob);
  uvec out(N);
  for (int i = 0; i < N; i++) {
    double u = unif_rand();
    int j = 0;
    while (u > probsum[j]) {j++;}
    out[i] = j;
  }
  return out;
}