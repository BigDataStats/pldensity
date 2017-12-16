library(pldensity)
library(dplyr)
library(ggmap)
library(latex2exp)

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
  thinprob = 0,
  discount = 0.9
) 

mod_trained1 <- mod_init %>% 
  ddpn_mix(x, epochs = 50) 
#
mod_trained2 <- mod_trained1 %>%
  ddpn_mix(x)

mod_trained3 <- mod_trained2 %>%
  ddpn_mix(x)
#
mod_trained4 <- mod_trained3 %>%
  ddpn_mix(x)

spatial_plot(mod_trained1, x, polygons = TRUE)
spatial_plot(mod_trained4, x, polygons = TRUE)


modtt <- ddpn_init(
  nparticles = 500,
  alpha = 25,
  lambda = c(-97.731970, 30.302445),
  kappa = .008,
  nu = 2,
  Omega =  0.008 ^ 2 * diag(2),
  rho = 0.9,
  thinprob = 0.001,
  discount = 0.99
) 

for (i in 1:16) {
  modtt <- modtt %>% 
    ddpn_mix(x[(10*(i - 1) + 1):(10 * i), , drop = FALSE])
}

spatial_plot(modtt, x, polygons = TRUE)

par(mar=c(35,35,100,10))
MyMap <- GetMap(center = c(30.302445, -97.731970), zoom = 10)
PlotOnStaticMap(MyMap, lat = x[ ,2], lon = x[ ,1], size = c(100, 100),
                cex = 1.4, pch = 19, col = "red", FUN = points, add = FALSE)
# title("sdfds")

# lines(1:100, 1:100)
# mod_trained %>%
#   ddpn_mix(x)


# library(rworldmap)
# newmap <- getMap(resolution = "high")
# plot(
#   newmap,
#   xlim = c(min(x[ ,2]), max(x[ ,2])),
#   ylim = c(min(x[ ,1]), max(x[ ,1]))
# )

austinmap <- get_map(location = c(-97.731970, 30.302445), maptype = "toner", zoom = 11)
austinmap_attr <- attr(austinmap, "bb")
resolution <- 128
mesh <- expand.grid(
  lon = seq(austinmap_attr$ll.lon, austinmap_attr$ur.lon, length.out = resolution),
  lat = seq(austinmap_attr$ll.lat, austinmap_attr$ur.lat, length.out = resolution)
)
density <- ddpn_eval(mod_trained1, data.matrix(mesh), nparticles = 20)

austinmap_density <- data.frame(
  lon = mesh[ ,1],
  lat = mesh[ ,2],
  density = density,
  alpha = log(density) / max(log(density))
)
# ggmap(austinmap) + 
#   geom_tile(data = austinmap_density, aes(x = lon, fill = density), alpha = 0.5) +
#   scale_fill_gradientn(colours = heat.colors(50))
#   

g <- ggmap(austinmap) + 
  geom_point(data = data.frame(lon = x[ ,1], lat = x[ ,2]), colour = "red", alpha = 0.9) +
  coord_cartesian() +
  coord_fixed() +
  geom_raster(data = austinmap_density, aes(x = lon, y = lat, fill = density, alpha = alpha)) +
  scale_fill_gradientn(colors = topo.colors(5)) +
  scale_alpha_continuous(guide = FALSE) +
  ggtitle(TeX("$\\alpha = 0.1$"))
g
