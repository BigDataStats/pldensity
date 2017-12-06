#include "plutils.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// 1. Dynamic Particles ===================================================================

// Defines a particle, which carries a essensial state
struct DynamicParticle {
  int m; // Total clusters
  vec c; // Observations per cluster
  vec w; // Stick breaking weights
  mat mu; // Mean (dim1) per cluster
  cube S; // SS (dim1, dim2) per cluster
  // Empty constructor
  DynamicParticle () {};
  // Full spec constructor
  DynamicParticle (int m_, arma::vec& c_, arma::vec w_, arma::mat& mu_, arma::cube& S_) 
    : m(m_), c(c_), w(w_), mu(mu_), S(S_) {};
  // copy constructor
  DynamicParticle (const DynamicParticle& z) { 
    m = z.m;
    w = z.w;
    c = z.c;
    mu = z.mu;
    S = z.S;
  }; 
  // Create particle with new center
  DynamicParticle (const arma::vec& x, double x_wt) {
    m = 1;
    c = {1};
    w = {x_wt};
    mu.resize(x.n_elem, 1);
    mu.col(0) = x;
    S.resize(x.n_elem, x.n_elem, 1);
    S.slice(0).fill(0.0);
  }
};

DynamicParticle read_dynamic_particle(const Rcpp::List& L) {
  int m = L["m"];
  vec c = L["c"];
  vec w = L["w"];
  mat mu = L["mu"];
  cube S = L["S"];
  return DynamicParticle(m, c, w, mu, S);
}

Rcpp::List list_dynamic_particle(const DynamicParticle& z) {
  Rcpp::List out = Rcpp::List::create(
    Named("m") = z.m,
    Named("c") = z.c,
    Named("w") = z.w,
    Named("mu") = z.mu,
    Named("S") = z.S
  );
  out.attr("class") = "HyperParam";
  return out;
}

//' @title test dynamic particle io
//' @export
// [[Rcpp::export]]
Rcpp::List iotest_dynamic_particle(Rcpp::List L) {
  DynamicParticle z = read_dynamic_particle(L);
  Rcpp::List out = list_dynamic_particle(z);
  return out;
}


// 2. HyperParameters ===================================================================


// Defines the prior, to avoid passing along all the parameters
struct DDPNHyperParam { // mu ~ N(lambda, S / kappa), S^-1 ~ W(Omega, nu)
  const double alpha; // Concentration in Dirichlet Process
  const arma::vec lambda; // Prior of mean
  const double kappa; // Prior difussion coeff
  const double nu; // Deg Freedom/Scale of Inv Wishart prior of Sigma
  const arma::mat Omega; // Shape Param of Inv Wishart prior of Sigma
  const double rho; // Stick-breaking weights autocorrelation 
  // Full spec constructor
  DDPNHyperParam(
    const double a_,
    const arma::vec l_,
    const double k_, 
    const double n_, 
    const arma::mat O_,
    const double r_) 
    : alpha(a_), lambda(l_), kappa(k_), nu(n_), Omega(O_), rho(r_) {};
};


DDPNHyperParam read_ddpn_hyperparam(const Rcpp::List& L) {
  double alpha = L["alpha"];
  double kappa = L["kappa"];
  double nu = L["nu"];
  double rho = L["rho"];
  vec lambda = L["lambda"];
  mat Omega = L["Omega"];
  return DDPNHyperParam(alpha, lambda, kappa, nu, Omega, rho);
}

Rcpp::List list_ddpn_hyperparam(const DDPNHyperParam& hp) {
  Rcpp::List out = Rcpp::List::create(
    Named("alpha") = hp.alpha,
    Named("lambda") = hp.lambda,
    Named("kappa") = hp.kappa,
    Named("nu") = hp.nu,
    Named("Omega") = hp.Omega,
    Named("rho") = hp.rho
  );
  out.attr("class") = "HyperParam";
  return out;
}


//' @title test hyperparam
//' @export
// [[Rcpp::export]]
Rcpp::List iotest_ddpn_hyperparam(Rcpp::List L) {
  DDPNHyperParam hp = read_ddpn_hyperparam(L);
  Rcpp::List out = list_ddpn_hyperparam(hp);
  return out;
}


struct DDPN {
  const int N; // Number of particles
  const DDPNHyperParam hp;
  std::vector< DynamicParticle > particle_list;
  DDPN (const int N_,
        const DDPNHyperParam& hp_, 
        const std::vector< DynamicParticle > p_)
    : N(N_), hp(hp_), particle_list(p_) {};
};


DDPN read_ddpn(const Rcpp::List& L) {
  const int N = L["N"];
  DDPNHyperParam hp = read_ddpn_hyperparam(L["hyper_param"]);
  Rcpp::List particle_list_ = L["particle_list"];
  std::vector<DynamicParticle> particle_list(N);
  for (int i = 0; i < N; i++) {
    Rcpp::List particle_ = particle_list_[i];
    particle_list[i] = read_dynamic_particle(particle_);
  }
  return DDPN(N, hp, particle_list);
}
  
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
    Named("N") = dppn.N,
    Named("hyper_param") = hyper_param,
    Named("particle_list") = particle_list
  );
  out.attr("class") = "DDPN";
  return out;
}

//' @title test ddpn I/O
//' @export
// [[Rcpp::export]]
Rcpp::List iotest_ddpn(Rcpp::List L) {
  DDPN model = read_ddpn(L);
  Rcpp::List out = list_ddpn(model);
  return out;
}


  
// 4. Code Structure ===========================================================================

inline arma::vec observation_prob(
    const arma::vec& x,
    const DynamicParticle& z,
    const DDPNHyperParam& hp
) {
  // Initialize
  double d = x.n_elem;
  vec out(z.m + 1);
  
  // Cluster contribution
  for (int j = 0; j < z.m; j++) {
    vec aj = (hp.kappa * hp.lambda + z.c[j] * z.mu.col(j)) / (hp.kappa + z.c[j]);
    mat Dj = z.S.slice(j) + hp.kappa * z.c[j] / (hp.kappa +  z.c[j]) * 
      (hp.lambda - z.mu.col(j)) * (hp.lambda - z.mu.col(j)).t();
    double cj = 2 * hp.nu + z.c[j] - d + 1.0;
    mat Bj = 2.0 * (hp.kappa + z.c[j] + 1.0) / (hp.kappa + z.c[j]) / cj * (hp.Omega  + 0.5 * Dj);
    out[j] = z.w[j] * dst(x, aj, Bj, cj);
  }
  
  // Prior contribution
  vec a0 = hp.lambda;
  double c0 = 2.0 * hp.nu - d + 1.0;
  const mat& B0 = 2.0 * (hp.kappa + 1.0) / hp.kappa / c0 * hp.Omega;
  out[z.m] = (1 - sum(z.w)) * dst(x, a0, B0, c0);
  
  return out;
}

// Evaluates the density of a new point given a particle
inline DynamicParticle update_particle(
    const arma::vec& xnew, 
    const DynamicParticle& z,
    const DDPNHyperParam& hp // iteration timestamp
) {
  // Initialize output
  DynamicParticle out(z);
  
  // Sample most likely component according to observation probs
  vec op = observation_prob(xnew, z, hp);
  int k = resample(1, op)[0]; 
  
  // Create new cluster if necessary and update particae
  if (k == z.m) {
    out.c.insert_rows(out.m, 1);
    out.mu.insert_cols(out.m, xnew);
    out.S.insert_slices(out.m, 1);
    out.m += 1;
    out.c[k] += 1;
  } else {
    out.c[k] += 1;
    out.mu.col(k) =  (z.c[k] * z.mu.col(k) + xnew) / out.c[k];
    out.S.slice(k) += (xnew * xnew.t()) + z.c[k] * (z.mu.col(k) * 
      z.mu.col(k).t()) - out.c[k] * (out.mu.col(k) * out.mu.col(k).t());  
  }
  
  return out;
}

