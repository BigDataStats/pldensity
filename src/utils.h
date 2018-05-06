#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 

double dst(
    const arma::vec& x,
    const arma::vec& mu,
    const arma::mat& Sigma, 
    const double df);

arma::vec dst(
    const arma::mat& X,
    const arma::vec& mu,
    const arma::mat& Sigma, 
    const double df);

arma::uvec resample(
    int N, 
    arma::vec prob);

#endif