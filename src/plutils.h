#ifndef PLUTILS_H
#define PLUTILS_H

#define ARMA_64BIT_WORD 1
#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 

double dst(
    const arma::vec& x,
    const arma::vec& mu,
    const arma::mat& Sigma, // Sigma^{-1}
    const double df);

// Samples from a multivariate categorical variable
arma::uvec resample(
    int N, 
    arma::vec prob);


#endif