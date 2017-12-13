library(pldensity)
library(dplyr)

data(austin)

x <- austin %>% 
  select(start_location_long, start_location_lat) %>%
  data.matrix()

mod <- ddpn_init(
  nparticles = 500,
  alpha = 2,
  lambda = c(-97.731970, 30.302445),
  kappa = .001,
  nu = 2,
  Omega =  0.01 ^ 2 * diag(2),
  rho = 0.7,
  thinprob = 0.0001
) 

mod_trained1 <- mod %>% 
  ddpn_mix(x, epochs = 5) 

mod_trained2 <- mod_trained1 %>% 
  ddpn_mix(x)

mod_trained3 <- mod_trained2 %>% 
  ddpn_mix(x)

mod_trained4 <- mod_trained3 %>% 
  ddpn_mix(x)
# mod_trained2 <- mod_trained %>% 
#   ddpn_mix(x)

spatial_plot(mod_trained1, x, polygons = FALSE)
spatial_plot(mod_trained2, x, polygons = FALSE)
spatial_plot(mod_trained3, x, polygons = FALSE)
spatial_plot(mod_trained4, x, polygons = TRUE)

# mod_trained %>% 
#   ddpn_mix(x) 
