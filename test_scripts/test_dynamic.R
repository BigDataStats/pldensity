library(pldensity)
library(dplyr)

data(austin)

x <- austin %>% 
  select(start_location_long, start_location_lat) %>%
  data.matrix()

mod <- ddpn_init(
  nparticles = 500,
  alpha = 8,
  lambda = c(-97.731970, 30.302445),
  kappa = .001,
  nu = 2,
  Omega =  0.01 ^ 2 * diag(2)
) 

mod_trained <- mod %>% 
  ddpn_mix(x)

spatial_plot(mod_trained, x, polygons = TRUE)

# mod_trained %>% 
#   ddpn_mix(x) 
