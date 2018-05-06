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
  double innerterm =  - 0.5*(df + d) * log(1.0 + as_scalar(xcentered.t() * InvSigmaReg * xcentered) / df);
  double extterm = lgamma(0.5*(df + d)) - lgamma(0.5*df) - 0.5*d * log(df * M_PI) - 0.5*ldet;
  if (isnan(extterm))
    stop("Numeric Problems!");
  return exp(extterm + innerterm);
}

//' @title multivariate t distribution density
//' @descripton Computes the density of a multivariate t
//' @export
// [[Rcpp::export]]
arma::vec dst(
    const arma::mat& X,
    const arma::vec& mu,
    const arma::mat& Sigma, // Sigma^{-1}
    const double df) {
  // Initialise
  vec out(X.n_rows, fill::zeros);
  
  // Sigma Inversion
  int d = mu.n_elem;
  mat InvSigmaReg = inv_sympd(Sigma + 1e-12 * eye<mat>(d ,d));
  double ldet, sign;
  log_det(ldet, sign, InvSigmaReg); 
 
  // Compute for each row
  for (uword i = 0; i < X.n_rows; i++) {
    vec xcentered = X.row(i).t() - mu;
    double innerterm =  - 0.5*(df + d) * log(1.0 + as_scalar(xcentered.t() * InvSigmaReg * xcentered) / df);
    double extterm = lgamma(0.5*(df + d)) - lgamma(0.5*df) - 0.5*d * log(df * M_PI) - 0.5*ldet;
    if (isnan(extterm))
      stop("Numeric Problems!");
    out[i] = exp(extterm + innerterm);
  }
  
  return out;
}

//' @title sample
//' @descripton samples from a categorical variable
//' @export
// [[Rcpp::export]]
arma::uvec sample(int N, arma::vec prob) {
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