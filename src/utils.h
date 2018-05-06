#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 


double dst(
    const arma::vec& x,
    const arma::vec& mu,
    const arma::mat& Sigma, // Sigma^{-1}
    const double df);

arma::vec dst(
    const arma::mat& X,
    const arma::vec& mu,
    const arma::mat& Sigma, // Sigma^{-1}
    const double df);

// Samples from a multivariate categorical variable
arma::uvec sample(
    int N, 
    arma::vec prob);

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

#endif