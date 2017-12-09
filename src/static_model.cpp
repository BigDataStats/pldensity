#include "pldensity.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. Static Particles =================================================================
  
// Defines a particle, carries an essential state
struct Particle {
  int m; // Total clusters
  vec c; // Observations per cluster
  mat mu; // Mean (dim1) per cluster
  cube S; // SS (dim1, dim2) per cluster
  // Empty constructor
  Particle () {};
  // Full spec constructor
  Particle (int m_, arma::vec& c_, arma::mat& mu_, arma::cube& S_) 
  : m(m_), c(c_), mu(mu_), S(S_) {};
  // copy constructor
  Particle (const Particle& z) { 
    m = z.m;
    c = z.c;
    mu = z.mu;
    S = z.S;
  }; 
  // Create particle with new center
  Particle (const arma::vec& x) {
    m = 1;
    c = {1};
    mu.resize(x.n_elem, 1);
    mu.col(0) = x;
    S.resize(x.n_elem, x.n_elem, 1);
    S.slice(0).fill(0.0);
  }
};

Particle read_particle(const Rcpp::List& L) {
  int m = L["m"];
  vec c = L["c"];
  mat mu = L["mu"];
  cube S = L["S"];
  return Particle(m, c, mu, S);
}

Rcpp::List list_particle(const Particle& z) {
  Rcpp::List out = Rcpp::List::create(
    Named("m") = z.m,
    Named("c") = z.c,
    Named("mu") = z.mu,
    Named("S") = z.S
  );
  out.attr("class") = "Particle";
  return out;
}

// //' @title test dynamic particle io
// //' @export
// // [[Rcpp::export]]
// Rcpp::List iotest_particle(Rcpp::List L) {
//   Particle z = read_particle(L);
//   Rcpp::List out = list_particle(z);
//   return out;
// }


// 2. HyperParameters ===================================================================
  
  // Defines the prior, to avoid passing along all the parameters
struct DPNHyperParam { // mu ~ N(lambda, S / kappa), S^-1 ~ W(Omega, nu)
  const double alpha; // Concentration in Dirichlet Process
  const arma::vec lambda; // Prior of mean
  const double kappa; // Prior difussion coeff
  const double nu; // Deg Freedom/Scale of Inv Wishart prior of Sigma
  const arma::mat Omega; // Shape Param of Inv Wishart prior of Sigma
  // Full spec constructor
  DPNHyperParam(
    const double a_,
    const arma::vec l_,
    const double k_, 
    const double n_, 
    const arma::mat O_)
    : alpha(a_), lambda(l_), kappa(k_), nu(n_), Omega(O_) {};
};


DPNHyperParam read_dpn_hyperparam(const Rcpp::List& L) {
  double alpha = L["alpha"];
  double kappa = L["kappa"];
  double nu = L["nu"];
  vec lambda = L["lambda"];
  mat Omega = L["Omega"];
  return DPNHyperParam(alpha, lambda, kappa, nu, Omega);
}

Rcpp::List list_dpn_hyperparam(const DPNHyperParam& hp) {
  Rcpp::List out = Rcpp::List::create(
    Named("alpha") = hp.alpha,
    Named("lambda") = hp.lambda,
    Named("kappa") = hp.kappa,
    Named("nu") = hp.nu,
    Named("Omega") = hp.Omega
  );
  out.attr("class") = "DPNHyperParam";
  return out;
}

// 
// //' @title test hyperparam
// //' @export
// // [[Rcpp::export]]
// Rcpp::List iotest_dpn_hyperparam(Rcpp::List L) {
//   DPNHyperParam hp = read_dpn_hyperparam(L);
//   Rcpp::List out = list_dpn_hyperparam(hp);
//   return out;
// }


struct DPN {
  const int N; // Number of particles
  const DPNHyperParam hp;
  std::vector< Particle > particle_list;
  // full spec constructor
  DPN (const int N_,
        const DPNHyperParam& hp_, 
        const std::vector<Particle> p_)
  : N(N_), hp(hp_), particle_list(p_) {};
  // copy constructor
  DPN (const DPN& mod)
    : N(mod.N), hp(mod.hp), particle_list(mod.particle_list) {};
};


DPN read_dpn(const Rcpp::List& L) {
  const int N = L["nparticles"];
  DPNHyperParam hp = read_dpn_hyperparam(L["hyper_param"]);
  Rcpp::List particle_list_ = L["particle_list"];
  std::vector<Particle> particle_list(N);
  for (int i = 0; i < N; i++) {
    Rcpp::List particle_ = particle_list_[i];
    particle_list[i] = read_particle(particle_);
  }
  return DPN(N, hp, particle_list);
}

Rcpp::List list_dpn(DPN dpn) {
  Rcpp::List particle_list;
  int N = dpn.N;
  for (int i = 0; i < N; i++) {
    Particle z = dpn.particle_list[i];
    Rcpp::List particle = list_particle(z);
    particle_list.push_back(particle);
  }
  Rcpp::List hyper_param = list_dpn_hyperparam(dpn.hp);
  Rcpp::List out = Rcpp::List::create(
    Named("nparticles") = dpn.N,
    Named("hyper_param") = hyper_param,
    Named("particle_list") = particle_list
  );
  out.attr("class") = "DPN";
  return out;
}

// //' @title test ddpn I/O
// //' @export
// // [[Rcpp::export]]
// Rcpp::List iotest_dpn(Rcpp::List L) {
//   DPN model = read_dpn(L);
//   Rcpp::List out = list_dpn(model);
//   return out;
// }



// 4. Observation Probability ===========================================================
  
inline arma::vec observation_prob(
  const arma::vec& x,
  const Particle& z,
  const DPNHyperParam& hp
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
    out[j] = z.c[j] * dst(x, aj, Bj, cj);
  }
  
  // Prior contribution
  vec a0 = hp.lambda;
  double c0 = 2.0 * hp.nu - d + 1.0;
  const mat& B0 = 2.0 * (hp.kappa + 1.0) / hp.kappa / c0 * hp.Omega;
  out[z.m] = hp.alpha * dst(x, a0, B0, c0);
  out /= hp.alpha + sum(z.c);
  
  return out;
}

// 5. Particle update rule for new data ==================================

// Evaluates the density of a new point given a particle
inline Particle update_particle(
  const arma::vec& xnew, 
  const Particle& z,
  const DPNHyperParam& hp // iteration timestamp
) {
  // Initialize output
  Particle out(z);
  
  // Choose most likely allocation given observation probabilities
  vec op = observation_prob(xnew, z, hp);
  int k = resample(1, op)[0]; 
  
  // k == z.m means observation starts a new cluster
  if (k == z.m) {
    out.c.insert_rows(out.m, 1);
    out.mu.insert_cols(out.m, xnew);
    out.S.insert_slices(out.m, 1);
    out.m += 1;
    out.c[k] += 1;
  } else {
    out.c[k] += 1;
    out.mu.col(k) =  (z.c[k] * z.mu.col(k) + xnew) / out.c[k];
    out.S.slice(k) += (xnew * xnew.t()) + z.c[k] * (z.mu.col(k) * z.mu.col(k).t()) - out.c[k] * (out.mu.col(k) * out.mu.col(k).t());  
  }
  
  return out;
}

// 6. Sequential Monte Carlo  =======================================================

//' @title Model initialisation
//' @export
// [[Rcpp::export]]
Rcpp::List dpn_init(
    const int nparticles,
    const double alpha,
    const arma::vec& lambda,
    const double kappa,
    const double nu,
    const arma::mat& Omega
) {
  // Hyper Parameters
  DPNHyperParam hp(alpha, lambda, kappa, nu, Omega);
  
  // Particle list
  std::vector<Particle> particle_list(nparticles);
  for (int i = 0; i < nparticles; i++) {
    particle_list[i] = Particle(lambda);
  }
  
  // Model
  DPN mod(nparticles, hp, particle_list);
  
  return list_dpn(mod);
}

//' @title Dirichlet Process Normal Mixture Kernel Density Estimation
//' @export
// [[Rcpp::export]]
Rcpp::List dpn_mix(
  const Rcpp::List& model,
  const arma::mat& x,
  const int epochs = 1
) {
  // Model as C++ object
  DPN mod = read_dpn(model);

  // Sequential Monte Carlo Loop
  for (int b = 0; b < epochs; b++) {
    for (uword t = 0; t < x.n_rows; t++) {
      // New observation
      vec xnew = x.row(t).t();
      
      // Resample according to observation probability given particle state
      vec weight(mod.N);
      for (int i = 0; i < mod.N; i++) {
        weight[i] = sum(observation_prob(xnew, mod.particle_list[i], mod.hp));
      }
      uvec resample_id = resample(mod.N, weight);
  
      // Propagate: resample particles and update new states
      std::vector<Particle> temp(mod.particle_list);
      for (int i = 0; i < mod.N; i++) {
        mod.particle_list[i] = update_particle(xnew, temp[resample_id[i]], mod.hp);
      }
    }
  }
  return list_dpn(mod);
}

// 7. Eval the density of the model ============================================

//' @title Eval Point Density
//' @export
// [[Rcpp::export]]
arma::vec dpn_eval(
    const Rcpp::List& model,
    const arma::mat& xnew,
    const int nparticles = 50
) {
  // Model as C++ object
  DPN mod = read_dpn(model);
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
Rcpp::List dpn_marginal(
    const Rcpp::List& model,
    const arma::uvec& dims
) {
  // Model to c++ object
  DPN mod = read_dpn(model);
  const uvec dims_ = dims - 1;
  
  // Marginal Prior
  DPNHyperParam mhp(
      mod.hp.alpha, 
      mod.hp.lambda(dims_),
      mod.hp.kappa,
      mod.hp.nu,
      mod.hp.Omega.submat(dims_, dims_));
    
  // Marginal Particle List
  std::vector<Particle> mparticle_list(mod.N);
  for (int i = 0; i < mod.N; i++) {
    Particle z = mod.particle_list[i];
    mat munew = z.mu.rows(dims_);
    cube S = z.S;
    cube Snew(dims_.n_elem, dims_.n_elem, z.m);
    for (int k = 0; k < z.m; k++) Snew.slice(k) = S.slice(k).submat(dims_, dims_);
    mparticle_list[i] = Particle(z.m, z.c, munew, Snew);
  }
  
  DPN marginal(mod.N, mhp, mparticle_list);
  return list_dpn(marginal);
}

//' @title Eval Point Conditional density
//' @export
// [[Rcpp::export]]
arma::mat dpn_conditional(
    const Rcpp::List& model,
    const arma::mat& xnew,
    const arma::uvec& eval_dims,
    const arma::uvec& condition_dims,
    const arma::mat& condition_values,
    const int nparticles = 50
) {
  // Model to c++ object
  DPN mod = read_dpn(model);
  const int N0 = min(mod.N, nparticles);
  
  // Allocate space
  const uvec edims = eval_dims - 1;
  const uvec cdims = condition_dims - 1;
  const int total_dims = edims.n_elem + cdims.n_elem;
  const uword K = xnew.n_rows, L = condition_values.n_rows;
  mat density(K, L, fill::zeros);
  
  // Marginal model
  DPN marg = read_dpn(dpn_marginal(model, cdims + 1));
  
  //  Eval conditional densities
  for (int i = 0; i < N0; i++) {
    Particle zfull = mod.particle_list[i];
    Particle zmarg = marg.particle_list[i];
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

