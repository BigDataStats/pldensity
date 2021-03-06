# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @name DynamicParticle
#' @title DynamicParticle class
#' @description This class is used inside C++ only and can't be called from R.
#' It defines a particle, which carries sufficient statistic for the  
#' distribution of a mixture of normals with dynamic stick-breaking weights'.
#' The sufficient statistics are:
#' \itemize{
#'   \item \eqn{M}: total number of clusters in current period
#'   \item \eqn{M_0}: total number of clusters in previous period
#'   \item \eqn{(c_j)_{j=1}^M}: counts of observations in all periods
#'   \item \eqn{(c^t_j)_{j=1}^M}: counts of observations in current period
#'   \item \eqn{(w_{j})_{j=1}^M}: stick-breaking weights of each cluster
#'   \item \eqn{(s_j)_{j=1}^M}: d-dim vector sum of each cluster
#'   \item \eqn{(SS_j)_{j=1}^M}: dxd-matrix of sum-of-squares of each cluster
#' }
#' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
#'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association, 
#'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
NULL

#' @name copy_dynamic_particle
#' @title Dynamic Particle Copying
#' @description Interncal C++ function for copying dynamic particles
NULL

#' @name read_dynamic_particle
#' @title From List in R to DynamicParticle in C++
#' @description This class is used inside C++ only and can't be called from R.
#' It reads a list of class DynamicParicle from R a constructs the corresponding object in C++
NULL

#' @name list_dynamic_particle
#' @title From DynamicParticle in C++ to list in R
#' @description This class is used inside C++ only and can't be called from R.
#' It reads a list of class DynamicParicle from R a constructs the corresponding object in C++
NULL

#' @name DDPNHyperParam
#' @title DDPNHyperParam class
#' @description This class is used inside C++ only and can't be called from R.
#' It carries the prior distribution for the parameters for Dirichlet Process Mixture
#' \deqn{x_i \mid (\mu, \Sigma) \sim N(\mu, \Sigma)}
#' \deqn{(\mu, \Sigma) \sim \text{DP}(\alpha, \text{NormalInvWishart}(\lambda, \kappa, \Omega, \nu))}
#' The stick-breaking weights follow a \eqn{BAR(\rho)} process. The hyperparameters are
#' \itemize{
#'   \item \eqn{\alpha}: The weights at time \eqn{\tau} have prior 
#'   \item \eqn{\lambda}: mean, \eqn{\mu \sim N(\lambda, \Sigma / \kappa)}
#'   \item \eqn{\kappa}: the shrinking factor
#'   \item \eqn{\Omega}: pseudo sum of squares, \eqn{\Sigma \sim \text{InvWishart}(\nu, \Omega)}
#'   \item \eqn{\nu}: pseudo degrees of freedom
#'   \item \eqn{\rho}: correlation degree from previous stick-breaking weights
#' }
#' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
#'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association, 
#'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
NULL

#' @name read_ddpn_hyperparam
#' @title From List in R to DDPNHyperParam in C++
#' @description This class is used inside C++ only and can't be called from R.
#' It reads a list of class DDPNHyperParam from R a constructs the corresponding object in C++
NULL

#' @name list_ddpn_hyperparam
#' @title From DDPNHyperParam in C++ to List in R of class DDPNHyperParam
#' @description This class is used inside C++ only and can't be called from R.
#' It reads a list of class DDPNHyperParam from R a constructs the corresponding object in C++
NULL

#' @name DDPN
#' @title DDPN class
#' @description This class is used inside C++ only and can't be called from R.
#' It carries the a collection of Dynamic Particles for fitting DP Dynamic Mixture with particle filters,
#' and the information for the prior.
#' \itemize{
#'   \item \eqn{N}: Number of particles
#'   \item \code{hp}: and object of class DDPNHyperParam
#'   \item \code{particle_list}: a vector of size \eqn{N} of objects of class DynamicParticle
#' }
#' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
#'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association, 
#'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
NULL

#' @name read_ddpn
#' @title From List in R to DDPN in C++
#' @description This class is used inside C++ only and can't be called from R.
#' It reads a list of class DDPN from R a constructs the corresponding object in C++
NULL

#' @name list_ddpn_hyperparam
#' @title From DDPN in C++ to List in R of class DDPN
#' @description This class is used inside C++ only and can't be called from R.
#' It reads a list of class DDPN from R a constructs the corresponding object in C++
NULL

#' @name cluster_predictive
#' @title Computes the cluster predictive distribution of a point
#' @description Computes the predictive distribution of each cluster for a new point \eqn{x}.
#' The updated hyperparameters are
#' @details The notation differs slightly from the paper
#' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
#'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
#'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
NULL

#' @name counts_weight
#' @title Particle counts likelihood
#' @description Computes the portion of the likelihood-based used to resample each particle relative to the counts
NULL

#' @name update_particle
#' @title Particle state update
#' @description Updates the statistics of a particle after seeing data point \code{xnew}. An allocation
#' of \code{xnew} takes place based on the cluster likelihood.
#' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
#'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
#'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
NULL

#' @title Convert a DynamicParticle R list to C++ and back
#' @description This function is used for unit testing that the process of creating DynamicParticles 
#' in C++ from an R-list and listing them back works correctly. Not be used for other purposes.
#' @export
iotest_dynamic_particle <- function(L) {
    .Call('_pldensity_iotest_dynamic_particle', PACKAGE = 'pldensity', L)
}

#' @title Convert a DDPNHyperParam R list to C++ and back
#' @description This function is used for unit testing that the process of creating DDPNHyperParam 
#' in C++ from an R-list and listing them back works correctly. Not be used for other purposes.
#' @export
iotest_ddpn_hyperparam <- function(L) {
    .Call('_pldensity_iotest_ddpn_hyperparam', PACKAGE = 'pldensity', L)
}

#' @title Convert a DDPN R list to C++ and back
#' @description This function is used for unit testing that the process of creating DynamicParticles 
#' in C++ from an R-list and listing them back works correctly. Not be used for other purposes
#' @export
iotest_ddpn <- function(L) {
    .Call('_pldensity_iotest_ddpn', PACKAGE = 'pldensity', L)
}

#' @title Model Initialisation
#' @description Inialisation function to create a new Dynamic Dirichlet Process Mixture of Normasl (DDPN).
#' @param nparticles Number of particle sused to approximate the posterior distribution. Too few particles will result in particle collapse
#' @param lambda dx1 vector. Prior mean for the mean
#' @param kappa real>0. Shrinking hyperparameter
#' @param nu real>d. Pseudo-degres of freedom
#' @param Omega dxd psd Matrix, Prior sum of squares
#' @param alpha real>0. The DP parameter. Higher values encourage more clusters
#' @param rho 0<real<1. The BAR autocorrelation parameter. Higher values give smoother transitions over time.
#' @details See the reference paper for the full specification of the model.
#' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
#'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
#'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
#' @export
ddpn_init <- function(nparticles, lambda, kappa, nu, Omega, alpha, rho = 0.0001) {
    .Call('_pldensity_ddpn_init', PACKAGE = 'pldensity', nparticles, lambda, kappa, nu, Omega, alpha, rho)
}

#' @title Model Training
#' @description Each a DDPN model and new data to learn the distribution.
#' @param model Object of class DDPN
#' @param x nxd matrix where each row is an observation
#' @param epochs int >=1. If it's larger than one, the data will be recycled. Defaults to 1.
#' @param new_period_every int >= -1.
#' @param resample_every int >= -1.
#' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
#'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
#'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
#' @export
ddpn_mix <- function(model, x, epochs = 1L, new_period_every = 0L, resample_every = 1L) {
    .Call('_pldensity_ddpn_mix', PACKAGE = 'pldensity', model, x, epochs, new_period_every, resample_every)
}

#' @title Model Training
#' @description Each a DDPN model and new data to learn the distribution.
#' @param model Object of class DDPN
#' @param x nxd matrix where each row is an observation
#' @param epochs int >=1. If it's larger than one, the data will be recycled. Defaults to 1.
#' @param new_period_every int >= 0. When set to zero, all the data corresponds to same time period, and the BAR occurs at the beginning. Otherwise, a BAR update will happen in the frequency specified.
#' @references Matthew A. Taddy (2012) Autoregressive Mixture Models for Dynamic Spatial Poisson Processes:
#'  Application to Tracking Intensity of Violent Crime, Journal of the American Statistical Association,
#'  105:492, 1403-1417, DOI: 10.1198/jasa.2010.ap09655
#' @export
ddpn_eval <- function(model, x, nparticles = 50L) {
    .Call('_pldensity_ddpn_eval', PACKAGE = 'pldensity', model, x, nparticles)
}

