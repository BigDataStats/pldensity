library(pldensity)

data("austin")

L <- list(
  alpha = 25, 
  lambda = runif(2), 
  kappa = .001, 
  nu = 2,
  Omega =  0.08 ^ 2 * diag(2),
  rho = 0.8
)
hp <- iotest_ddpn_hyperparam(L)

L <- list(
  m = 2,
  c = c(1,10),
  w = c(.01, .09),
  mu = matrix(0, 3, 2),
  S = array(0, dim = c(3, 3, 2))
)

dp <- iotest_dynamic_particle(L)

L <- list(
  N = 3,
  hyper_param = hp,
  particle_list = list(dp, dp, dp)
)

model <- iotest_ddpn(L)