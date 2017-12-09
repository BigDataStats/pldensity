
x <- austin %>% 
  select(start_location_long, start_location_lat) %>% 
  data.matrix()

# Define and train the model ===========================

# Note: Best HyperParameters can be obtained with cross-validation
mod <- dpn_init(
  nparticles = 250,
  alpha = 50,
  lambda = c(-97.731970, 30.302445),
  kappa = .01,
  nu = 2,
  Omega =  0.01 ^ 2 * diag(2)
) 

mod_trained <- mod %>% 
  dpn_mix(x, epochs = 2)

density <- mod_trained %>% 
  dpn_eval(x, nparticles = 50)

# Plot some results ========================================

spdata <- sp::SpatialPointsDataFrame(x, data.frame(density))

pal <- colorNumeric("Spectral", domain = density)

spmap <- leaflet(spdata) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>% 
  setView(lat = 30.302445, lng = -97.731970, zoom = 11) %>% 
  addCircleMarkers(
    radius = 4, 
    stroke = FALSE, 
    color = ~pal(density),
    fillOpacity = 1
  ) %>% 
  addLegend(
    "bottomright",
    values = ~density,
    pal = pal) 
spmap

# The entire heatmap ===============================================

# Eval density
resol <- 25
xseq <-  seq(-98.014, -97.58, length.out = resol)
yseq <- seq(30.13, 30.52, length.out = resol)
mesh <- data.matrix(expand.grid(xseq, yseq, length.out = resol))

z <- mod_trained %>% 
  dpn_eval(xnew = mesh[ ,c(1,2)], nparticles = 100)

z <- matrix(z, resol, resol)

plot(x)
contour(xseq, yseq, z, nlevels = 50, add = TRUE, col = "blue")

cl <- contourLines(xseq, yseq, z, nlevels = 500)
for (i in seq_along(cl)) {
  spmap <- spmap %>% addPolygons(
    lng = cl[[i]]$x, 
    lat = cl[[i]]$y, 
    fillColor = "red", 
    fillOpacity = .05,
    stroke = FALSE)  
}
spmap <- spmap %>% 
  addCircleMarkers(
  radius = 4, 
  stroke = FALSE, 
  color = ~pal(density),
  fillOpacity = 1
) 
spmap


# Now let's jointly learn origin and destination and take conditionals ============

x <- austin %>% 
  select(
    start_location_long, 
    start_location_lat,
    end_location_long, 
    end_location_lat
  ) %>% 
  data.matrix()
  

mod2 <- dpn_init(
  nparticles = 1000,
  alpha = 50,
  lambda = c(-97.731970, 30.302445, -97.731970, 30.302445),
  kappa = .01,
  nu = 2,
  Omega =  0.01 ^ 2 * diag(4)
) 

mod2_trained <- mod2 %>% 
  dpn_mix(x, epochs = 1)

airport_coords <- c(-97.666676,  30.204509)

# Eval density
resol <- 25
xseq <-  seq(-98.014, -97.58, length.out = resol)
yseq <- seq(30.13, 30.52, length.out = resol)
mesh <- data.matrix(expand.grid(xseq, yseq, length.out = resol))

z <- mod2_trained %>% 
    dpn_conditional(
      xnew = mesh[ ,c(1,2)], 
      nparticles = 100,
      eval_dims = c(3, 4),
      condition_dims = c(1, 2),
      condition_values = matrix(airtport_coords, nrow = 1)
    )

z <- matrix(z, resol, resol)

spmap_cond <- leaflet(spdata) %>% 
  addProviderTiles(providers$CartoDB.Positron) %>% 
  setView(lat = 30.302445, lng = -97.731970, zoom = 11) %>% 
  addMarkers(lng = ~airport_coords[1], lat = ~airport_coords[2])

cl <- contourLines(xseq, yseq, z, nlevels = 100)
for (i in seq_along(cl)) {
  spmap_cond <- spmap_cond %>% addPolygons(
    lng = cl[[i]]$x, 
    lat = cl[[i]]$y, 
    fillColor = "red", 
    fillOpacity = .03,
    stroke = FALSE)  
}
spmap_cond
