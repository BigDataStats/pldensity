// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// iotest_dynamic_particle
Rcpp::List iotest_dynamic_particle(Rcpp::List L);
RcppExport SEXP _pldensity_iotest_dynamic_particle(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(iotest_dynamic_particle(L));
    return rcpp_result_gen;
END_RCPP
}
// iotest_ddpn_hyperparam
Rcpp::List iotest_ddpn_hyperparam(Rcpp::List L);
RcppExport SEXP _pldensity_iotest_ddpn_hyperparam(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(iotest_ddpn_hyperparam(L));
    return rcpp_result_gen;
END_RCPP
}
// iotest_ddpn
Rcpp::List iotest_ddpn(Rcpp::List L);
RcppExport SEXP _pldensity_iotest_ddpn(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(iotest_ddpn(L));
    return rcpp_result_gen;
END_RCPP
}
// ddpn_init
Rcpp::List ddpn_init(const int nparticles, const arma::vec& lambda, const double kappa, const double nu, const arma::mat& Omega, const double alpha, const double rho);
RcppExport SEXP _pldensity_ddpn_init(SEXP nparticlesSEXP, SEXP lambdaSEXP, SEXP kappaSEXP, SEXP nuSEXP, SEXP OmegaSEXP, SEXP alphaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nparticles(nparticlesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(ddpn_init(nparticles, lambda, kappa, nu, Omega, alpha, rho));
    return rcpp_result_gen;
END_RCPP
}
// ddpn_mix
Rcpp::List ddpn_mix(const Rcpp::List& model, const arma::mat& x, const int epochs, int new_period_every);
RcppExport SEXP _pldensity_ddpn_mix(SEXP modelSEXP, SEXP xSEXP, SEXP epochsSEXP, SEXP new_period_everySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< int >::type new_period_every(new_period_everySEXP);
    rcpp_result_gen = Rcpp::wrap(ddpn_mix(model, x, epochs, new_period_every));
    return rcpp_result_gen;
END_RCPP
}
// ddpn_eval
arma::vec ddpn_eval(const Rcpp::List& model, const arma::mat& x, const int nparticles);
RcppExport SEXP _pldensity_ddpn_eval(SEXP modelSEXP, SEXP xSEXP, SEXP nparticlesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type nparticles(nparticlesSEXP);
    rcpp_result_gen = Rcpp::wrap(ddpn_eval(model, x, nparticles));
    return rcpp_result_gen;
END_RCPP
}
// dpn_init
Rcpp::List dpn_init(const int nparticles, const double alpha, const arma::vec& lambda, const double kappa, const double nu, const arma::mat& Omega);
RcppExport SEXP _pldensity_dpn_init(SEXP nparticlesSEXP, SEXP alphaSEXP, SEXP lambdaSEXP, SEXP kappaSEXP, SEXP nuSEXP, SEXP OmegaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nparticles(nparticlesSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    rcpp_result_gen = Rcpp::wrap(dpn_init(nparticles, alpha, lambda, kappa, nu, Omega));
    return rcpp_result_gen;
END_RCPP
}
// dpn_mix
Rcpp::List dpn_mix(const Rcpp::List& model, const arma::mat& x, const int epochs);
RcppExport SEXP _pldensity_dpn_mix(SEXP modelSEXP, SEXP xSEXP, SEXP epochsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type epochs(epochsSEXP);
    rcpp_result_gen = Rcpp::wrap(dpn_mix(model, x, epochs));
    return rcpp_result_gen;
END_RCPP
}
// dpn_eval
arma::vec dpn_eval(const Rcpp::List& model, const arma::mat& xnew, const int nparticles);
RcppExport SEXP _pldensity_dpn_eval(SEXP modelSEXP, SEXP xnewSEXP, SEXP nparticlesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xnew(xnewSEXP);
    Rcpp::traits::input_parameter< const int >::type nparticles(nparticlesSEXP);
    rcpp_result_gen = Rcpp::wrap(dpn_eval(model, xnew, nparticles));
    return rcpp_result_gen;
END_RCPP
}
// dpn_marginal
Rcpp::List dpn_marginal(const Rcpp::List& model, const arma::uvec& dims);
RcppExport SEXP _pldensity_dpn_marginal(SEXP modelSEXP, SEXP dimsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type dims(dimsSEXP);
    rcpp_result_gen = Rcpp::wrap(dpn_marginal(model, dims));
    return rcpp_result_gen;
END_RCPP
}
// dpn_conditional
arma::mat dpn_conditional(const Rcpp::List& model, const arma::mat& xnew, const arma::uvec& eval_dims, const arma::uvec& condition_dims, const arma::mat& condition_values, const int nparticles);
RcppExport SEXP _pldensity_dpn_conditional(SEXP modelSEXP, SEXP xnewSEXP, SEXP eval_dimsSEXP, SEXP condition_dimsSEXP, SEXP condition_valuesSEXP, SEXP nparticlesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xnew(xnewSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type eval_dims(eval_dimsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type condition_dims(condition_dimsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type condition_values(condition_valuesSEXP);
    Rcpp::traits::input_parameter< const int >::type nparticles(nparticlesSEXP);
    rcpp_result_gen = Rcpp::wrap(dpn_conditional(model, xnew, eval_dims, condition_dims, condition_values, nparticles));
    return rcpp_result_gen;
END_RCPP
}
// dst
arma::vec dst(const arma::mat& X, const arma::vec& mu, const arma::mat& Sigma, const double df);
RcppExport SEXP _pldensity_dst(SEXP XSEXP, SEXP muSEXP, SEXP SigmaSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const double >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(dst(X, mu, Sigma, df));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pldensity_iotest_dynamic_particle", (DL_FUNC) &_pldensity_iotest_dynamic_particle, 1},
    {"_pldensity_iotest_ddpn_hyperparam", (DL_FUNC) &_pldensity_iotest_ddpn_hyperparam, 1},
    {"_pldensity_iotest_ddpn", (DL_FUNC) &_pldensity_iotest_ddpn, 1},
    {"_pldensity_ddpn_init", (DL_FUNC) &_pldensity_ddpn_init, 7},
    {"_pldensity_ddpn_mix", (DL_FUNC) &_pldensity_ddpn_mix, 4},
    {"_pldensity_ddpn_eval", (DL_FUNC) &_pldensity_ddpn_eval, 3},
    {"_pldensity_dpn_init", (DL_FUNC) &_pldensity_dpn_init, 6},
    {"_pldensity_dpn_mix", (DL_FUNC) &_pldensity_dpn_mix, 3},
    {"_pldensity_dpn_eval", (DL_FUNC) &_pldensity_dpn_eval, 3},
    {"_pldensity_dpn_marginal", (DL_FUNC) &_pldensity_dpn_marginal, 2},
    {"_pldensity_dpn_conditional", (DL_FUNC) &_pldensity_dpn_conditional, 6},
    {"_pldensity_dst", (DL_FUNC) &_pldensity_dst, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_pldensity(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
