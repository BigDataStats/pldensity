#' @title Poisson DLM Random Walk
#' @description  Simpel function for filtering with a random walk model
#' @param counts numeric vector with counts
#' @param m0 numeric starting point for state variable mean
#' @param C0 numeric>0 starting point for state variable variance
#' @param nu0 numeric>0 starting degrees of freedom (pseudocounts) for random walk variance
#' @param xi0 numeric>0 starting sum of squares for random walk variance
#' @param delta 0<numeric<=1 correlation parameter
#' @param nparticles number of particles for the filter
#' @export
# [[Rcpp::export]]
firstorder_poisson_dlm <- function(counts, m0, C0, nu0, xi0, delta, nparticles=100) {
  periods <- length(counts)
  
  # reserve space || Keeeping many unnecessary stuff for the moment
  m <- matrix(0, periods + 1, nparticles)
  C <- matrix(0, periods + 1, nparticles)
  nu <- matrix(0, periods + 1, nparticles)
  xi <- matrix(0, periods + 1, nparticles)
  V <- matrix(0, periods + 1, nparticles)
  R <- matrix(0, periods + 1, nparticles)
  Q <- matrix(0, periods + 1, nparticles)
  f <- matrix(0, periods + 1, nparticles)
  logLambda <- matrix(0, periods + 1, nparticles)
  eta <- matrix(0, periods + 1, nparticles)
  
  # initialise
  m[1, ] <- m0
  C[1, ] <- C0
  nu[1, ] <- nu0
  xi[1, ] <- xi0
  V[1, ] <- 1 / rgamma(nparticles, 0.5*nu[1], rate = 0.5*xi[1])
  logLambda[1, ] <- m0
  
  seqlikelihood <- numeric(periods)
  # iterate
  for (t in 1:periods) {
    # 1. propagate and obtain weights
    wts <- numeric(nparticles)
    for (i in 1:nparticles) {
      # forecast moments
      R[t + 1, i] <- C[t, i] / delta
      f[t + 1, i] <- m[t, i]
      Q[t + 1, i] <- R[t + 1, i] + V[t, i]
      
      # propose
      logLambda[t + 1, i] <- rnorm(1, log(counts[t]), 1 / sqrt(counts[t]))
      
      # compute weights
      wts[i] <- dpois(counts[t], exp(logLambda[t + 1, i])) *
        dnorm(logLambda[t + 1, i], f[t + 1, i], sqrt(Q[t + 1, i])) /
        dnorm(logLambda[t + 1, i], log(counts[t]), 1 / sqrt(counts[t]))
    }
    
    # 2. resample particles (brute force way)
    ind <- sample(nparticles, replace = TRUE, prob = wts)
    logLambda <- logLambda[ ,ind]
    m <- m[ ,ind]
    C <- C[ ,ind]
    nu <- nu[ ,ind]
    xi <- xi[ ,ind]
    V <- V[ ,ind]
    R <- R[ ,ind]
    Q <- Q[ ,ind]
    f <- f[ ,ind]
    eta <- eta[, ind]
    
    # 3. Update particles
    for (i in 1:nparticles) {
      eta[t + 1, i] <- rnorm(1, m[t, i], R[t + 1, i])
      nu[t + 1, i] <- nu[t, i] + 1
      xi[t + 1, i] <- xi[t, i] + (logLambda[t + 1, i] - f[t + 1, i])^2
      V[t + 1, i] <- 1 / rgamma(1,  0.5*nu[t + 1, i], rate = 0.5*xi[t + 1, i])
      m[t + 1, i] <- m[t, i] + (R[t + 1, i] / Q[t + 1, i]) * (logLambda[t + 1, i] - f[t + 1, i])
      C[t + 1, i] <- (R[t + 1, i] / Q[t + 1, i]) * V[t, i]
      
      seqlikelihood[t] <-  seqlikelihood[t] + dpois(counts[t], exp(logLambda[t + 1, i])) / nparticles 
    }
  }
  
  list(
    logintensity = logLambda[-1, ],
    seqlikelihood = seqlikelihood,
    logmarginal = sum(log(seqlikelihood)),
    V = V[-1, ]
  )
}