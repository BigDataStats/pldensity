library(pldensity)
library(dplyr)

data(austin)

x <- austin %>% 
  select(start_location_long, start_location_lat) %>% 
  data.matrix()

mod <- ddpn_init(
  nparticles = 250,
  alpha = 50,
  lambda = c(-97.731970, 30.302445),
  kappa = .01,
  nu = 2,
  Omega =  0.01 ^ 2 * diag(2)
) 

mod_thin <- mod %>% 
  ddpn_mix(x)
