#define ARMA_64BIT_WORD 1
#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 


using namespace Rcpp;
using namespace arma;
using namespace std;

// 0. Utils ===========================================================================

double dst(
    const arma::vec& x,
    const arma::vec& mu,
    const arma::mat& Sigma, 
    const double df) {
  
  int d = mu.n_elem;
  vec xcentered = x - mu;
  mat InvSigmaReg = inv_sympd(Sigma);
  double ldet;
  double sign;
  log_det(ldet, sign, InvSigmaReg); 
  double innerterm =  - 0.5*(df + d) * log(1.0 + as_scalar(xcentered.t() * InvSigmaReg * xcentered) / df);
  double extterm = lgamma(0.5*(df + d)) - lgamma(0.5*df) - 0.5*d * log(df * M_PI) + 0.5*ldet;
  double out = exp(extterm + innerterm);
  if (isnan(out))
    stop("Numeric Problems!");
  return out;
}


arma::vec dst(
    const arma::mat& X,
    const arma::vec& mu,
    const arma::mat& Sigma, 
    const double df) {
  // Initialise
  vec out(X.n_rows, fill::zeros);
  
  // Sigma Inversion
  int d = mu.n_elem;
  mat InvSigmaReg = inv_sympd(Sigma);
  double ldet, sign;
  log_det(ldet, sign, InvSigmaReg); 
  
  // Compute for each row
  for (uword i = 0; i < X.n_rows; i++) {
    vec xcentered = X.row(i).t() - mu;
    double innerterm =  - 0.5*(df + d) * log(1.0 + as_scalar(xcentered.t() * InvSigmaReg * xcentered) / df);
    double extterm = lgamma(0.5*(df + d)) - lgamma(0.5*df) - 0.5*d * log(df * M_PI) + 0.5*ldet;
    out[i] = exp(extterm + innerterm);
    if (isnan(out[i]))
      stop("Numeric Problems!");
  }
  
  return out;
}


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


// 1. Dynamic Particles ===================================================================

//' @name DynamicParticle
//' @title DynamicParticle class
//' @description This class is used inside C++ only and can't be called from R.
//' It defines a particle, which carries sufficient statistic for the  
//' distribution of a mixture of normals with dynamic stick-breaking weights'.
//' The sufficient statistics are:
//' \itemize{
//'   \item \eqn{M}: total number of clusters in current period
//'   \item \eqn{M_0}: total number of clusters in previous period
//'   \item \eqn{(c_j)_{j=1}^M}: counts of observations in all periods
//'   \item \eqn{(c^t_j)_{j=1}^M}: counts of observations in current period
//'   \item \eqn{(w_{j})_{j=1}^M}: stick-breaking weights of each cluster
//'   \item \eqn{(s_j)_{j=1}^M}: d-dim vector sum of each cluster
//'   \item \eqn{(SS_j)_{j=1}^M}: dxd-matrix of sum-of-squares of each cluster
//' }
//' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
//'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association, 
//'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
struct DynamicParticle {
  
  int M; // number of particles
  int M0; // number of particles in previous cluster
  vec c; // number of counts per cluster in all period
  vec ct; // number of counts per cluster in current period
  vec w; // weights of clusters
  mat s; // dxM matrix of sums per cluster
  cube SS; // dxdxM cube of sum-of-squares outer-product per cluster
  
  // Empty constructor
  DynamicParticle () {};
  
  // Full constructor
  DynamicParticle (
      int M_, 
      int M0_, 
      arma::vec& c_,
      arma::vec& ct_,
      arma::vec& w_, 
      arma::mat& s_, 
      arma::cube& SS_
    ) : M(M_), M0(M0_), c(c_), ct(ct_), w(w_), s(s_), SS(SS_) {};

  // Copy constructor 
  DynamicParticle (const DynamicParticle& z) { 
    M = z.M;
    M0 = z.M0;
    c = arma::vec(z.c);
    ct = arma::vec(z.ct);
    w = arma::vec(z.w);
    s = arma::mat(z.s);
    SS = arma::cube(z.SS);  
  }; 
  
  // Copy constructor with argument for fast copy using pointers
  DynamicParticle (DynamicParticle& z, const bool copy = false) { 
    M = z.M;
    M0 = z.M0;
    int d = z.s.n_rows;
    c = arma::vec(&z.c[0], M, copy);
    ct = arma::vec(&z.ct[0], M, copy);
    w = arma::vec(&z.w[0], M, copy);
    s = arma::mat(&z.s[0], d, M, copy);
    SS = arma::cube(&z.SS[0], d, d, M, copy);  
  }; 
  
  // Create a particle with one cluster from a point
  DynamicParticle (const arma::vec& x, double w_init) {
    M = 0;
    M0 = 0;
    c = {0};
    ct = {0};
    w = {w_init};
    s.resize(x.n_elem, 1);
    s.col(0) = x;
    SS.resize(x.n_elem, x.n_elem, 1);
    SS.slice(0).fill(0.0);
  }
};

//' @name copy_dynamic_particle
//' @title Dynamic Particle Copying
//' @description Interncal C++ function for copying dynamic particles
DynamicParticle copy_dynamic_particle(DynamicParticle& z) { 
  return DynamicParticle(z);
}


//' @name read_dynamic_particle
//' @title From List in R to DynamicParticle in C++
//' @description This class is used inside C++ only and can't be called from R.
//' It reads a list of class DynamicParicle from R a constructs the corresponding object in C++
DynamicParticle read_dynamic_particle(const Rcpp::List& L) {
  if (!L.inherits("DynamicParticle"))
    stop("Not a DynamicParticle list from R!");
  int M = L["M"];
  int M0 = L["M0"];
  vec c = L["c"];
  vec ct = L["ct"];
  vec w = L["w"];
  mat s = L["s"];
  cube SS = L["SS"];
  return DynamicParticle(M, M0, c, ct, w, s, SS);
}

//' @name list_dynamic_particle
//' @title From DynamicParticle in C++ to list in R
//' @description This class is used inside C++ only and can't be called from R.
//' It reads a list of class DynamicParicle from R a constructs the corresponding object in C++
Rcpp::List list_dynamic_particle(const DynamicParticle& z) {
  Rcpp::List out = Rcpp::List::create(
    Named("M") = z.M,
    Named("M0") = z.M0,
    Named("c") = z.c,
    Named("ct") = z.ct,
    Named("w") = z.w,
    Named("s") = z.s,
    Named("SS") = z.SS
  );
  out.attr("class") = "DynamicParticle";
  return out;
}


//' @title Convert a DynamicParticle R list to C++ and back
//' @description This function is used for unit testing that the process of creating DynamicParticles 
//' in C++ from an R-list and listing them back works correctly. Not be used for other purposes.
//' @export
// [[Rcpp::export]]
Rcpp::List iotest_dynamic_particle(Rcpp::List L) {
  DynamicParticle z = read_dynamic_particle(L);
  Rcpp::List out = list_dynamic_particle(z);
  return out;
}


// 2. HyperParameters ===================================================================

//' @name DDPNHyperParam
//' @title DDPNHyperParam class
//' @description This class is used inside C++ only and can't be called from R.
//' It carries the prior distribution for the parameters for Dirichlet Process Mixture
//' \deqn{x_i \mid (\mu, \Sigma) \sim N(\mu, \Sigma)}
//' \deqn{(\mu, \Sigma) \sim \text{DP}(\alpha, \text{NormalInvWishart}(\lambda, \kappa, \Omega, \nu))}
//' The stick-breaking weights follow a \eqn{BAR(\rho)} process. The hyperparameters are
//' \itemize{
//'   \item \eqn{\alpha}: The weights at time \eqn{\tau} have prior 
//'   \item \eqn{\lambda}: mean, \eqn{\mu \sim N(\lambda, \Sigma / \kappa)}
//'   \item \eqn{\kappa}: the shrinking factor
//'   \item \eqn{\Omega}: pseudo sum of squares, \eqn{\Sigma \sim \text{InvWishart}(\nu, \Omega)}
//'   \item \eqn{\nu}: pseudo degrees of freedom
//'   \item \eqn{\rho}: correlation degree from previous stick-breaking weights
//' }
//' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
//'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association, 
//'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
struct DDPNHyperParam { // 

  const arma::vec lambda; // Prior mean of mean
  const double kappa; // Shrink factor
  const double nu; // Deg Freedom/Scale of Inv Wishart prior of Sigma
  const arma::mat Omega; // Shape Param of Inv Wishart prior of Sigma
  const double alpha; // Concentration in Dirichlet Process
  const double rho; // Stick-breaking weights autocorrelation

  // Full spec constructor
  DDPNHyperParam(
    const arma::vec l_,
    const double k_,
    const double n_,
    const arma::mat O_,
    const double a_,
    const double r_)
    : lambda(l_), kappa(k_), nu(n_), Omega(O_), alpha(a_), rho(r_) {};
};


//' @name read_ddpn_hyperparam
//' @title From List in R to DDPNHyperParam in C++
//' @description This class is used inside C++ only and can't be called from R.
//' It reads a list of class DDPNHyperParam from R a constructs the corresponding object in C++
DDPNHyperParam read_ddpn_hyperparam(const Rcpp::List& L) {
  if (!L.inherits("DDPNHyperParam"))
    stop("Not a DDPNHyperParam list from R!");
  double alpha = L["alpha"];
  double kappa = L["kappa"];
  double nu = L["nu"];
  double rho = L["rho"];
  vec lambda = L["lambda"];
  mat Omega = L["Omega"];
  return DDPNHyperParam(lambda, kappa, nu, Omega, alpha, rho);
}

//' @name list_ddpn_hyperparam
//' @title From DDPNHyperParam in C++ to List in R of class DDPNHyperParam
//' @description This class is used inside C++ only and can't be called from R.
//' It reads a list of class DDPNHyperParam from R a constructs the corresponding object in C++
Rcpp::List list_ddpn_hyperparam(const DDPNHyperParam& hp) {
  Rcpp::List out = Rcpp::List::create(
    Named("lambda") = hp.lambda,
    Named("kappa") = hp.kappa,
    Named("nu") = hp.nu,
    Named("Omega") = hp.Omega,
    Named("alpha") = hp.alpha,
    Named("rho") = hp.rho
  );
  out.attr("class") = "DDPNHyperParam";
  return out;
}


//' @title Convert a DDPNHyperParam R list to C++ and back
//' @description This function is used for unit testing that the process of creating DDPNHyperParam 
//' in C++ from an R-list and listing them back works correctly. Not be used for other purposes.
//' @export
// [[Rcpp::export]]
Rcpp::List iotest_ddpn_hyperparam(Rcpp::List L) {
  DDPNHyperParam hp = read_ddpn_hyperparam(L);
  Rcpp::List out = list_ddpn_hyperparam(hp);
  return out;
}

// 2. DDPN ===================================================================

//' @name DDPN
//' @title DDPN class
//' @description This class is used inside C++ only and can't be called from R.
//' It carries the a collection of Dynamic Particles for fitting DP Dynamic Mixture with particle filters,
//' and the information for the prior.
//' \itemize{
//'   \item \eqn{N}: Number of particles
//'   \item \code{hp}: and object of class DDPNHyperParam
//'   \item \code{particle_list}: a vector of size \eqn{N} of objects of class DynamicParticle
//' }
//' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
//'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association, 
//'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
struct DDPN {
  const int N; // Number of particles
  const DDPNHyperParam hp;
  std::vector< DynamicParticle > particle_list;
  DDPN (const int N_,
        const DDPNHyperParam& hp_,
        const std::vector< DynamicParticle > p_)
    : N(N_), hp(hp_), particle_list(p_) {};
};

//' @name read_ddpn
//' @title From List in R to DDPN in C++
//' @description This class is used inside C++ only and can't be called from R.
//' It reads a list of class DDPN from R a constructs the corresponding object in C++
DDPN read_ddpn(const Rcpp::List& L) {
  if (!L.inherits("DDPN"))
    stop("Not a DDPN list from R!");
  const int N = L["nparticles"];
  DDPNHyperParam hp = read_ddpn_hyperparam(L["hyper_param"]);
  Rcpp::List particle_list_ = L["particle_list"];
  std::vector<DynamicParticle> particle_list(N);
  for (int i = 0; i < N; i++) {
    Rcpp::List particle_ = particle_list_[i];
    particle_list[i] = read_dynamic_particle(particle_);
  }
  return DDPN(N, hp, particle_list);
}

// void reset_ddpn(DDPN& mod) {
//   for (int i = 0; i < mod.N; i++)  {
//     mod.particle_list[i].M0 = mod.particle_list[i].M;
//     mod.particle_list[i].c = mod.particle_list[i].M;
//   }
// }


//' @name list_ddpn_hyperparam
//' @title From DDPN in C++ to List in R of class DDPN
//' @description This class is used inside C++ only and can't be called from R.
//' It reads a list of class DDPN from R a constructs the corresponding object in C++
Rcpp::List list_ddpn(DDPN dppn) {
  Rcpp::List particle_list;
  int N = dppn.N;
  for (int i = 0; i < N; i++) {
    DynamicParticle z = dppn.particle_list[i];
    Rcpp::List particle = list_dynamic_particle(z);
    particle_list.push_back(particle);
  }
  Rcpp::List hyper_param = list_ddpn_hyperparam(dppn.hp);
  Rcpp::List out = Rcpp::List::create(
    Named("nparticles") = dppn.N,
    Named("hyper_param") = hyper_param,
    Named("particle_list") = particle_list
  );
  out.attr("class") = "DDPN";
  return out;
}


//' @title Convert a DDPN R list to C++ and back
//' @description This function is used for unit testing that the process of creating DynamicParticles 
//' in C++ from an R-list and listing them back works correctly. Not be used for other purposes
//' @export
// [[Rcpp::export]]
Rcpp::List iotest_ddpn(Rcpp::List L) {
  DDPN model = read_ddpn(L);
  Rcpp::List out = list_ddpn(model);
  return out;
}


// 4. Updates ===========================================================================
 
 //' @name BAR_update
 //' @title Update stick-breaking weights
 //' @description Reweights clusters for the next time period following the
 //' \eqn{\text{BAR}(\rho)} autoregressive process for Beta distributions.
 //' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
 //'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
 //'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
inline void BAR_update(
  DynamicParticle& z,
  const DDPNHyperParam& hp
) {
  // BAR propagation from previous time perod weights
  if (z.M0 > 0) {
    vec u = Rcpp::as<vec>(rbeta(z.M0, hp.alpha, 1 - hp.rho));
    vec v = Rcpp::as<vec>(rbeta(z.M0, hp.rho, 1 - hp.rho));
    z.w.subvec(0, z.M0 - 1) = 1 - u % (1 - v % z.w.subvec(0, z.M0 - 1));
  }
  
  // Stick-breaking weights of current period
  if (z.M0 < z.M) {
    for (int l = z.M0; l < z.M; l++) {
      double csum = (l < z.M - 1) ? sum(z.ct.subvec(l + 1, z.M - 1)) : 0;
      z.w(l) = R::rbeta(1 + z.ct(l), hp.alpha + csum);
    }
  }
}
 
//' @name cluster_predictive
//' @title Computes the cluster predictive distribution of a point
//' @description Computes the predictive distribution of each cluster for a new point \eqn{x}.
//' The updated hyperparameters are
//' @details The notation differs slightly from the paper
//' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
//'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
//'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
inline arma::vec cluster_predictive(
    const arma::vec& x,
    const DynamicParticle& z,
    const DDPNHyperParam& hp
) {
  // Initialise
  int M = z.M;
  double d = x.n_elem;
  vec out(M + 1);

  double total_wt = 0.;
  
  // cluster contributions
  for (int j = 0; j < M; j++) {
    // cluster weight
    double cp_weight = z.w(j);
    if (j >= 1) 
      cp_weight *= prod(1 - z.w.subvec(0, j - 1));

    // cluster posterior predictive
    double cj = z.c[j];
    double kappa_nj = hp.kappa + cj;
    double nu_nj = hp.nu + cj;
    vec meanj = z.s.col(j) / cj;
    vec mu_nj = (hp.kappa*hp.lambda + cj * meanj) / kappa_nj;
    mat Omega_nj = hp.Omega + z.SS.slice(j) - cj * meanj * meanj.t() +
      (hp.kappa * cj / kappa_nj) * (meanj - hp.lambda) * (meanj - hp.lambda).t();
    mat Sigma_nj = (1 + 1 / kappa_nj) * Omega_nj / (nu_nj - d + 1);
    out[j] = cp_weight * dst(x, mu_nj, Sigma_nj, nu_nj);
    total_wt += cp_weight;
  }

  // Prior contribution
  double prior_weight = prod(1 - z.w);
  out[M] = prior_weight * dst(x, hp.lambda, (1 + 1 / hp.kappa) * hp.Omega / (hp.nu - d + 1), hp.nu);

  // Rcout << "clusts: " << total_wt << "  prior: " << prior_weight << "   total: " << total_wt + prior_weight << endl;
  return out;
}

//' @name counts_weight
//' @title Particle counts likelihood
//' @description Computes the portion of the likelihood-based used to resample each particle relative to the counts
inline double counts_weight(
    const arma::vec& xnew,
    const DynamicParticle& z,
    const DDPNHyperParam& hp
) {
  // Allocate space
  double out = 1.0;

  // Likelihood of the counts
  for (int l = 0; l < z.M0; l++) {
    // Rcout << "c(l):  " << z.c(l) << "  n:  " << sum(z.c.subvec(l, z.M - 1)) << "  p:  " << z.w(l) << endl;
    out *= R::dbinom(z.ct(l), sum(z.ct.subvec(l, z.M - 1)), z.w(l), false); 
  }
  
  return out;
}

//' @name update_particle
//' @title Particle state update
//' @description Updates the statistics of a particle after seeing data point \code{xnew}. An allocation
//' of \code{xnew} takes place based on the cluster likelihood.
//' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
//'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
//'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
inline void update_particle(
    DynamicParticle& z,
    const arma::vec& xnew,
    const DDPNHyperParam& hp
) {

  // Choose most likely allocation given observation probabilities
  vec cl_prob = cluster_predictive(xnew, z, hp);
  int k = resample(1, cl_prob)[0];
  // Rcout << "cl weights: " << cl_prob / sum(cl_prob) << "  k: " << k << endl;
  
  // if k == M start a new cluster
  int M = z.M;
  if (k == M) {
    z.w.insert_rows(M, 1); // the value is always updated later
    z.c.insert_rows(M, 1);
    z.ct.insert_rows(M, 1);
    z.s.insert_cols(M, 1);
    z.SS.insert_slices(M, 1);
    z.M += 1;
  } 
  z.c(k) += 1;
  z.ct(k) += 1;
  z.s.col(k) += xnew;
  z.SS.slice(k) += xnew * xnew.t();
}

// 5. Sequential Monte Carlo  =======================================================

//' @title Model Initialisation
//' @description Inialisation function to create a new Dynamic Dirichlet Process Mixture of Normasl (DDPN).
//' @param nparticles Number of particle sused to approximate the posterior distribution. Too few particles will result in particle collapse
//' @param lambda dx1 vector. Prior mean for the mean
//' @param kappa real>0. Shrinking hyperparameter
//' @param nu real>d. Pseudo-degres of freedom
//' @param Omega dxd psd Matrix, Prior sum of squares
//' @param alpha real>0. The DP parameter. Higher values encourage more clusters
//' @param rho 0<real<1. The BAR autocorrelation parameter. Higher values give smoother transitions over time.
//' @details See the reference paper for the full specification of the model.
//' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
//'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
//'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
//' @export
// [[Rcpp::export]]
Rcpp::List ddpn_init(
    const int nparticles,
    const arma::vec& lambda,
    const double kappa,
    const double nu,
    const arma::mat& Omega,
    const double alpha,
    const double rho = 0.0001 // independent fit
) {
  // Hyper Parameters
  DDPNHyperParam hp(lambda, kappa, nu, Omega, alpha, rho);

  // Particle list
  std::vector<DynamicParticle> particle_list(nparticles);
  NumericVector beta = Rcpp::rbeta(nparticles, 1, hp.alpha);
  for (int i = 0; i < nparticles; i++) {
    particle_list[i] = DynamicParticle(lambda, beta[i]);
  }

  // Model
  DDPN mod(nparticles, hp, particle_list);

  return list_ddpn(mod);
}

//' @title Model Training
//' @description Each a DDPN model and new data to learn the distribution.
//' @param model Object of class DDPN
//' @param x nxd matrix where each row is an observation
//' @param epochs int >=1. If it's larger than one, the data will be recycled. Defaults to 1.
//' @param new_period_every int >= -1.
//' @param resample_every int >= -1.
//' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
//'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
//'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
//' @export
// [[Rcpp::export]]
Rcpp::List ddpn_mix(
    const Rcpp::List& model,
    const arma::mat& x,
    const int epochs = 1,
    int new_period_every = 0,
    int resample_every = 1
) {
  // Validate
  if (new_period_every < 0)
    stop("new_period_every must be non-negative. Set to 0 for only one period.");
  if (new_period_every == 0)
    new_period_every = x.n_rows;
  if (resample_every < 1)
    stop("resample_every must be positive.");

  // Model as C++ object
  DDPN mod = read_ddpn(model);
  
  // Initialise Marginal Likelihhood
  vec seq_likelihood(x.n_rows, fill::zeros);
  vec resample_weight(mod.N, fill::ones);
  
  // Sequential Monte Carlo Loop
  for (int b = 0; b < epochs; b++) {
    for (uword t = 0; t < x.n_rows; t++) {
      // Reset M
      if (new_period_every == 1 || ((t % new_period_every) == 0)) 
        for (int i = 0; i < mod.N; i++) {
          mod.particle_list[i].M0 = mod.particle_list[i].M;
          mod.particle_list[i].ct.fill(0);
        }
          
          
      // New observation
      vec xnew = x.row(t).t();

      // Bar update of weights
      for (int i = 0; i < mod.N; i++) 
        BAR_update(mod.particle_list[i], mod.hp);
      
      // Resampling weights
      for (int i = 0; i < mod.N; i++) {
        vec pred_w = cluster_predictive(xnew, mod.particle_list[i], mod.hp);
        double cnts_w = counts_weight(xnew, mod.particle_list[i], mod.hp);
        resample_weight[i] *= cnts_w * sum(pred_w);
        seq_likelihood[t] += sum(pred_w) / mod.N;
      }
        
      if (resample_every == 1 || ((t % resample_every) == 0)) {
        // Resample 
        double sum_weight = sum(resample_weight);
        if (sum_weight > 0)
          resample_weight = resample_weight / sum(resample_weight);
        uvec zeta = resample(mod.N, resample_weight + 1 / pow(mod.N, 2)); // for particle stability
        std::vector<DynamicParticle> temp(mod.N);
        for (int i = 0; i < mod.N; i++) 
          temp[i] =  copy_dynamic_particle(mod.particle_list[i]);
        for (int i = 0; i < mod.N; i++) 
          mod.particle_list[i] =  copy_dynamic_particle(temp[zeta[i]]);
        // Rcout << "Resampling iter " << ((int) t) << endl;
        // Rcout << max(resample_weight + 1/pow(mod.N, 2) / sum(resample_weight + 1/pow(mod.N, 2))) << endl;
        // Rcout << min(resample_weight + 1/pow(mod.N, 2) / sum(resample_weight + 1/pow(mod.N, 2))) << endl;
        resample_weight.fill(1.);  
      }
      
      // Learn and update particles
      for (int i = 0; i < mod.N; i++) 
        update_particle(mod.particle_list[i], xnew, mod.hp);
    }
  }

  return Rcpp::List::create(Named("updated_model") = list_ddpn(mod),
                            Named("sequential_likelihood") = seq_likelihood,
                            Named("marginal_loglikelihood") = sum(log(seq_likelihood)));
}

//' @title Model Training
//' @description Each a DDPN model and new data to learn the distribution.
//' @param model Object of class DDPN
//' @param x nxd matrix where each row is an observation
//' @param epochs int >=1. If it's larger than one, the data will be recycled. Defaults to 1.
//' @param new_period_every int >= 0. When set to zero, all the data corresponds to same time period, and the BAR occurs at the beginning. Otherwise, a BAR update will happen in the frequency specified.
//' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
//'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
//'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
//' @export
// [[Rcpp::export]]
arma::vec ddpn_eval(
    const Rcpp::List& model,
    const arma::mat& x,
    const int nparticles = 50
) {
  // Model as C++ object
  DDPN mod = read_ddpn(model);
  DDPNHyperParam hp = mod.hp;
  
  int N0 = min(mod.N, nparticles);
  double d = x.n_cols;
  
  // Allocate space
  vec out(x.n_rows, fill::zeros);

  // Eval density
  for (int i = 0; i < N0; i++) {
    DynamicParticle z = mod.particle_list[i];
    int M = z.M;
    // double total = 0;
    
    // Cluster contributions: we use the version of the t-student density that only factors Sigma once for all data
    for (int j = 0; j < M; j++) {
      // cluster weight
      
      double cnts_w = z.w(j);
      if (j >= 1) {
        cnts_w *= prod(1 - z.w.subvec(0, j - 1));
      }
      // total += cnts_w;
      
      // cluster posterior paramtersdev
      double cj = z.c[j];
      double kappa_nj = hp.kappa + cj;
      double nu_nj = hp.nu + cj;
      vec meanj = z.s.col(j) / cj;
      vec mu_nj = (hp.kappa*hp.lambda + cj * meanj) / kappa_nj;
      mat Omega_nj = hp.Omega + z.SS.slice(j) - cj * meanj * meanj.t() +
        (hp.kappa * cj / kappa_nj) * (meanj - hp.lambda) * (meanj - hp.lambda).t();
      mat Sigma_nj = (1 + 1 / kappa_nj) * Omega_nj / (nu_nj - d + 1);
      
      // t-distribution for all data points
      vec clusterp_w = dst(x, mu_nj, Sigma_nj, nu_nj); 
      
      // Add to individual
      out += cnts_w * clusterp_w / N0;
    }
    
    // Prior contribution
    double prior_weight = prod(1 - z.w);
    // total += prior_weight;
    out += prior_weight * dst(x, hp.lambda, (1 + 1 / hp.kappa) * hp.Omega / (hp.nu - d + 1), hp.nu) / N0;
    // Rcout << total << endl;
  }

  return out;
}