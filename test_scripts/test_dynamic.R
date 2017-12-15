library(pldensity)
library(dplyr)

data(austin)

x <- austin %>% 
  select(start_location_long, start_location_lat) %>%
  data.matrix()

mod_init <- ddpn_init(
  nparticles = 500,
  alpha = 25,
  lambda = c(-97.731970, 30.302445),
  kappa = .008,
  nu = 2,
  Omega =  0.008 ^ 2 * diag(2),
  rho = 0.8,
  thinprob = 0
) 

mod_trained1 <- mod_init %>% 
  ddpn_mix(x) 
#
mod_trained2 <- mod_trained1 %>%
  ddpn_mix(x)

mod_trained3 <- mod_trained2 %>%
  ddpn_mix(x)
#
mod_trained4 <- mod_trained3 %>%
  ddpn_mix(x)

modtt <- mod_init
modtt$hyper_param$rho <- 0.5
for (i in 1:16) {
  modtt <- modtt %>% 
    ddpn_mix(x[(10*(i - 1) + 1):(10 * i), , drop = FALSE])
}

spatial_plot(mod_trained1, x, polygons = TRUE)
spatial_plot(mod_trained4, x, polygons = TRUE)
spatial_plot(modtt, x, polygons = TRUE)


# mod_trained %>%
#   ddpn_mix(x)
