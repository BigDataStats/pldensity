#define ARMA_64BIT_WORD 1
#define _USE_MATH_DEFINES

#include <RcppArmadillo.h>
#include "utils.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. Dynamic Particles ===================================================================

// Defines a particle, which carries a essensial state
struct DynamicParticle {
  int m; // Total clusters
  vec c; // Observations per clusters in current period
  int m0; // Cluster at the end of last period
  vec w; // Stick breaking weights previous and current
  mat mu; // Mean (dim1) per cluster
  cube S; // SS (dim1, dim2) per cluster
  // Empty constructor
  DynamicParticle () {};
  // Full spec constructor
  DynamicParticle (
      int m_, 
      arma::vec& c_, 
      int m0_, 
      arma::vec w_, 
      arma::mat& mu_, 
      arma::cube& S_
    ) 
    : m(m_), c(c_), m0(m0_), w(w_), mu(mu_), S(S_) {};
  // copy constructor
  DynamicParticle (const DynamicParticle& z) { 
    m = z.m;
    w = arma::vec(z.w);
    m0 = z.m0;
    c = arma::vec(z.c);
    mu = arma::mat(z.mu);
    S = arma::cube(z.S);
  }; 
  // Create particle with new center
  DynamicParticle (const arma::vec& x, double x_wt) {
    m = 1;
    c = {1};
    m0 = 1;
    w = {x_wt};
    mu.resize(x.n_elem, 1);
    mu.col(0) = x;
    S.resize(x.n_elem, x.n_elem, 1);
    S.slice(0).fill(0.0);
  }
};

DynamicParticle read_dynamic_particle(const Rcpp::List& L) {
  int m = L["m"];
  int m0 = L["m0"];
  vec c = L["c"];
  vec w = L["w"];
  mat mu = L["mu"];
  cube S = L["S"];
  return DynamicParticle(m, c, m0, w, mu, S);
}

Rcpp::List list_dynamic_particle(const DynamicParticle& z) {
  Rcpp::List out = Rcpp::List::create(
    Named("m") = z.m,
    Named("c") = z.c,
    Named("m0") = z.m0,
    Named("w") = z.w,
    Named("mu") = z.mu,
    Named("S") = z.S
  );
  out.attr("class") = "DynamicParticle";
  return out;
}

// //' @title test dynamic particle io
// //' @export
// // [[Rcpp::export]]
// Rcpp::List iotest_dynamic_particle(Rcpp::List L) {
//   DynamicParticle z = read_dynamic_particle(L);
//   Rcpp::List out = list_dynamic_particle(z);
//   return out;
// }


// 2. HyperParameters ===================================================================


// Defines the prior, to avoid passing along all the parameters
struct DDPNHyperParam { // mu ~ N(lambda, S / kappa), S^-1 ~ W(Omega, nu)
  const double alpha; // Concentration in Dirichlet Process
  const arma::vec lambda; // Prior of mean
  const double kappa; // Prior difussion coeff
  const double nu; // Deg Freedom/Scale of Inv Wishart prior of Sigma
  const arma::mat Omega; // Shape Param of Inv Wishart prior of Sigma
  const double rho; // Stick-breaking weights autocorrelation
  const double thinprob; // Clusters to randomly kill at each stage
  const double discount;
  // Full spec constructor
  DDPNHyperParam(
    const double a_,
    const arma::vec l_,
    const double k_, 
    const double n_, 
    const arma::mat O_,
    const double r_,
    const double th_,
    const double d_) 
    : alpha(a_), lambda(l_), kappa(k_), nu(n_), Omega(O_), 
      rho(r_), thinprob(th_), discount(d_) {};
};


DDPNHyperParam read_ddpn_hyperparam(const Rcpp::List& L) {
  double alpha = L["alpha"];
  double kappa = L["kappa"];
  double nu = L["nu"];
  double rho = L["rho"];
  double thinprob = L["thinprob"];
  double discount = L["discount"];
  vec lambda = L["lambda"];
  mat Omega = L["Omega"];
  return DDPNHyperParam(alpha, lambda, kappa, nu, Omega, rho, thinprob, discount);
}

Rcpp::List list_ddpn_hyperparam(const DDPNHyperParam& hp) {
  Rcpp::List out = Rcpp::List::create(
    Named("alpha") = hp.alpha,
    Named("lambda") = hp.lambda,
    Named("kappa") = hp.kappa,
    Named("nu") = hp.nu,
    Named("Omega") = hp.Omega,
    Named("rho") = hp.rho,
    Named("thinprob") = hp.thinprob,
    Named("discount") = hp.discount
  );
  out.attr("class") = "DHyperParam";
  return out;
}


// //' @title test hyperparam
// //' @export
// // [[Rcpp::export]]
// Rcpp::List iotest_ddpn_hyperparam(Rcpp::List L) {
//   DDPNHyperParam hp = read_ddpn_hyperparam(L);
//   Rcpp::List out = list_ddpn_hyperparam(hp);
//   return out;
// }


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

//' @title test ddpn I/O
//' @export
// [[Rcpp::export]]
Rcpp::List iotest_ddpn(Rcpp::List L) {
  DDPN model = read_ddpn(L);
  Rcpp::List out = list_ddpn(model);
  return out;
}


  
// 4. Update Modules ===========================================================================

inline void reweight_clusters(
  DynamicParticle& z,
  const DDPN& mod
) {
  // BAR propagation from previous time perod weights
  vec u = Rcpp::as<vec>(rbeta(z.m0, mod.hp.alpha, 1 - mod.hp.rho));
  vec v = Rcpp::as<vec>(rbeta(z.m0, mod.hp.rho, 1 - mod.hp.rho));
  z.w.subvec(0, z.m0 - 1) = 1 - u % (1 - v % z.w.subvec(0, z.m0 - 1));
  
  // Stick-breaking of current period
  if (z.m0 < z.m) {
    double cnew_sum = sum(z.w.subvec(z.m0, z.m - 1));
    for (int l = z.m0; l < z.m; l++) {
      z.w(l) = R::rbeta(1 + z.c(l), mod.hp.alpha + cnew_sum);
    } 
  }
}

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
    // Obtain weight
    double cp_weight = z.w(j);
    if (j >= 1) {
      cp_weight *= prod(1 - z.w.subvec(0, j - 1)); 
    }
    
    // Obtain posterior
    vec aj = (hp.kappa * hp.lambda + z.c[j] * z.mu.col(j)) / (hp.kappa + z.c[j]);
    mat Dj = z.S.slice(j) + hp.kappa * z.c[j] / (hp.kappa +  z.c[j]) * 
      (hp.lambda - z.mu.col(j)) * (hp.lambda - z.mu.col(j)).t();
    double cj = 2 * hp.nu + z.c[j] - d + 1.0;
    mat Bj = 2.0 * (hp.kappa + z.c[j] + 1.0) / (hp.kappa + z.c[j]) / cj * (hp.Omega  + 0.5 * Dj);
    out[j] = cp_weight * dst(x, aj, Bj, cj);
  }
  
  // Prior contribution
  double prior_weight = prod(1 - z.w);
  vec a0 = hp.lambda;
  double c0 = 2.0 * hp.nu - d + 1.0;
  const mat& B0 = 2.0 * (hp.kappa + 1.0) / hp.kappa / c0 * hp.Omega;
  out[z.m] = prior_weight * dst(x, a0, B0, c0);
  
  return out;
}

inline double resampling_weight(
    const arma::vec& xnew,
    const DynamicParticle& z,
    const DDPNHyperParam& hp
) {
  // Allcoate space
  double out = 1.0;
  
  // likelihood of p(c_{1,...,m0} | wnew)
  for (int l = 0; l < z.m0; l++) {
    out *= R::dbinom(z.c(l), sum(z.c.subvec(l, z.m0 - 1)), z.w(l), false);
  }
  
  // Observaton probability
  out *= sum(observation_prob(xnew, z, hp));
    
  return out;
}

inline void update_particle(
    DynamicParticle& z,
    const arma::vec& xnew,
    const DDPNHyperParam& hp // iteration timestamp
) {
  // Choose most likely allocation given observation probabilities
  vec op = observation_prob(xnew, z, hp);
  int k = resample(1, op)[0];
  
  // k == z.m means observation starts a new cluster
  if (k == z.m) {
    z.w.insert_rows(z.m, 1); // the value is always updated later
    z.c.insert_rows(z.m, 1);
    z.mu.insert_cols(z.m, xnew);
    z.S.insert_slices(z.m, 1);
    z.m += 1;
    z.c(k) += 1;
  } else {
    z.c(k) += 1;
    vec temp_mu = z.mu.col(k);
    z.mu.col(k) =  ((z.c[k] - 1) * temp_mu + xnew) / z.c[k];
    z.S.slice(k) += (xnew * xnew.t()) + (z.c[k] - 1) * (temp_mu * temp_mu.t()) - z.c[k] * (z.mu.col(k) * z.mu.col(k).t());
  }
}



// This should be run at the beginning of each period to the previous results
inline void thinning(
    DDPN& model
) {
  // Choose most likely allocation given observation probabilities
  for (int i = 0; i < model.N; i++) {
    int m = 0;
    DynamicParticle particle(model.particle_list[i]);
    particle.c = floor(model.hp.discount * particle.c);
    for (int j = 0; j < particle.m; j++) {
      if (particle.c[j] > 0 && particle.c[j] / particle.m > model.hp.thinprob) {
        m++;
      }
    }
    if (m > 0) {
      vec c(particle.c);
      vec w(m);
      mat mu(model.hp.lambda.n_elem, m);
      cube S(model.hp.lambda.n_elem, model.hp.lambda.n_elem, m);
      m = 0;
      for (int j = 0; j < particle.m; j++) {
        if (particle.c[j] > 0 && particle.c[j] / particle.m > model.hp.thinprob) {
          w(m) = particle.w(j);
          mu.col(m) = particle.mu.col(j);
          S.slice(m) = particle.S.slice(j);
          m++;
        }
      }
      model.particle_list[i] = DynamicParticle(m, c, m, w, mu, S);
    } else{
      double beta = R::rbeta(1.0, model.hp.alpha);
      model.particle_list[i] = DynamicParticle(model.hp.lambda, beta);
    }
   
  }
} 
// 
// //' @title thinning test i/o
// //' @export
// // [[Rcpp::export]]
// Rcpp::List thinning(Rcpp::List &model)

// 6. Sequential Monte Carlo  =======================================================

//' @title Model initialisation
//' @export
// [[Rcpp::export]]
Rcpp::List ddpn_init(
    const int nparticles,
    const double alpha,
    const arma::vec& lambda,
    const double kappa,
    const double nu,
    const arma::mat& Omega,
    const double rho = 0.8,
    const double thinprob = 0.001,
    const double discount = 0.99
) {
  // Hyper Parameters
  DDPNHyperParam hp(alpha, lambda, kappa, nu, Omega, rho, thinprob, discount);
  
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

//' @title Dirichlet Process Normal Mixture Kernel Density Estimation
//' @export
// [[Rcpp::export]]
Rcpp::List ddpn_mix(
    const Rcpp::List& model,
    const arma::mat& x,
    const int epochs = 1
) {
  // Model as C++ object
  DDPN mod = read_ddpn(model);
  
  // Drop clusters with small probability
  thinning(mod);
  
  // Sequential Monte Carlo Loop
  for (int b = 0; b < epochs; b++) {
    for (uword t = 0; t < x.n_rows; t++) {
      // New observation
      vec xnew = x.row(t).t();
      
      // Update cluster weights
      for (int i = 0; i < mod.N; i++) {
        reweight_clusters(mod.particle_list[i], mod);
      }
      
      // Obtaining resampling weights from emission probs
      vec weight(mod.N);
      for (int i = 0; i < mod.N; i++) {
        weight[i] = resampling_weight(xnew, mod.particle_list[i], mod.hp);
      }
      
      // Resample 
      uvec zeta = resample(mod.N, weight);
      std::vector<DynamicParticle> temp(mod.particle_list);
      for (int i = 0; i < mod.N; i++) {
        mod.particle_list[i] = DynamicParticle(temp[zeta[i]]);
      }
      
      // Learn and update particles
      for (int i = 0; i < mod.N; i++) {
        update_particle(mod.particle_list[i], xnew, mod.hp);
      }
    }
  }
  return list_ddpn(mod);
}

//' @title Eval Point Density
//' @export
// [[Rcpp::export]]
arma::vec ddpn_eval(
    const Rcpp::List& model,
    const arma::mat& xnew,
    const int nparticles = 50
) {
  // Model as C++ object
  DDPN mod = read_ddpn(model);
  int N0 = min(mod.N, nparticles);
  
  // Allocate space
  vec out(xnew.n_rows, fill::zeros);
  
  // Eval density
  for (int i = 0; i < N0; i++) {
    for (uword t = 0; t < xnew.n_rows; t++) {
      vec xi = xnew.row(t).t();
      out[t] += sum(observation_prob(xi, mod.particle_list[i], mod.hp)) / N0;
    }
  }
  
  return out;
}




// 8. Marginal and Conditional =========================================================

//' @title Take marginals
//' @export
// [[Rcpp::export]]
Rcpp::List ddpn_marginal(
    const Rcpp::List& model,
    const arma::uvec& dims
) {
  // Model to c++ object
  DDPN mod = read_ddpn(model);
  const uvec dims_ = dims - 1;
  
  // Marginal Prior
  DDPNHyperParam mhp(
      mod.hp.alpha, 
      mod.hp.lambda(dims_),
      mod.hp.kappa,
      mod.hp.nu,
      mod.hp.Omega.submat(dims_, dims_),
      mod.hp.rho,
      mod.hp.thinprob,
      mod.hp.discount);
  
  // Marginal Particle List
  std::vector<DynamicParticle> mparticle_list(mod.N);
  for (int i = 0; i < mod.N; i++) {
    DynamicParticle z = mod.particle_list[i];
    mat munew = z.mu.rows(dims_);
    cube S = z.S;
    cube Snew(dims_.n_elem, dims_.n_elem, z.m);
    for (int k = 0; k < z.m; k++) Snew.slice(k) = S.slice(k).submat(dims_, dims_);
    mparticle_list[i] = DynamicParticle(z.m, z.c, z.m0, z.w, munew, Snew);
  }
  
  DDPN marginal(mod.N, mhp, mparticle_list);
  return list_ddpn(marginal);
}

//' @title Eval Point Conditional density
//' @export
// [[Rcpp::export]]
arma::mat ddpn_conditional(
    const Rcpp::List& model,
    const arma::mat& xnew,
    const arma::uvec& eval_dims,
    const arma::uvec& condition_dims,
    const arma::mat& condition_values,
    const int nparticles = 50
) {
  // Model to c++ object
  DDPN mod = read_ddpn(model);
  const int N0 = min(mod.N, nparticles);
  
  // Allocate space
  const uvec edims = eval_dims - 1;
  const uvec cdims = condition_dims - 1;
  const int total_dims = edims.n_elem + cdims.n_elem;
  const uword K = xnew.n_rows, L = condition_values.n_rows;
  mat density(K, L, fill::zeros);
  
  // Marginal model
  DDPN marg = read_ddpn(ddpn_marginal(model, cdims + 1));
  
  //  Eval conditional densities
  for (int i = 0; i < N0; i++) {
    const DynamicParticle &zfull = mod.particle_list[i];
    const DynamicParticle &zmarg = marg.particle_list[i];
    for (uword l = 0; l < L; l++) {
      for (uword k = 0; k < K; k++) {
        vec xfull(total_dims, fill::zeros);
        xfull(edims) = xnew.row(k).t();
        xfull(cdims) = condition_values.row(l).t();
        double djoint = sum(observation_prob(xfull, zfull, mod.hp));
        double dmarg = sum(observation_prob(condition_values.row(l).t(), zmarg, marg.hp));
        density(k, l) += djoint / dmarg / N0;
      }
    }
  }
  
  return density;
}

