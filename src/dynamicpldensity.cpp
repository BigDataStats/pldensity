#define ARMA_64BIT_WORD 1
#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 


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


// Defines a particle, which carries a essensial state
struct DynamicParticle {
  int m; // Total clusters
  vec n; // Observations per cluster
  vec w; // Stick breaking weights
  mat mu; // Mean (dim1) per cluster
  cube S; // SS (dim1, dim2) per cluster
  DynamicParticle (int m_, arma::vec& n_, arma::mat& mu_, arma::cube& S_) 
    : m(m_), n(n_), mu(mu_), S(S_) {};
  DynamicParticle (const DynamicParticle& z) { // copy constructor
    m = z.m;
    n = z.n;
    mu = z.mu;
    S = z.S;
  }; 
  DynamicParticle () {};
  DynamicParticle (const arma::vec& x) {
    m = 1;
    n.resize(1);
    n[0] = 1;
    mu.resize(x.n_elem, 1);
    mu.col(0) = x;
    S.resize(x.n_elem, x.n_elem, 1);
    S.slice(0).fill(0.0);
  }
};
