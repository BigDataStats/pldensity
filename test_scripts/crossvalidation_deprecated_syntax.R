# Library ========================================================
library("pldensity") # Our package
library("tidyverse") # Data manipulation
library("progress") # Progress of the cross-validation
library("leaflet") # Geo maps
library("sp") # Spatial data type

# 0. Show how to crossvalidate
t0_train <- lubridate::ymd_hms("2017-04-10 18:00:00")
t1_train <- lubridate::ymd_hms("2017-04-10 18:59:59")
t0_test <- lubridate::ymd_hms("2017-03-27 18:00:00")
t1_test <- lubridate::ymd_hms("2017-03-27 18:59:59")

# 1. Read Data ===================================================
rides <- read_csv("D:/Datasets/Rides_A.csv")

x_train <- rides %>%
  filter(started_on >= t0_train & started_on <= t1_train) %>% 
  select(start_location_lat, start_location_long) %>% 
  data.matrix()

x_test <- rides %>%
  filter(started_on >= t0_test & started_on <= t1_test) %>% 
  select(start_location_lat, start_location_long) %>% 
  data.matrix()

# 2. Parameters Grid =============================================

params <- expand.grid(
  alpha = c(100, 150, 200),
  omega_scale = c(.001, .003, .006, .01),
  kappa = c(.001, .003, .005, .01),
  N = c(500, 750, 200)
)

# 3. Cross-validation ============================================

params$training_time <- 0
params$likelihood <- 0
pb <- progress_bar$new(total = nrow(params)) # progress bar

# Cross-validate
for (i in 1:nrow(params)) {
  # Train on test data
  st <- system.time({
    mod <- dp_normal_mix(
      x_train,
      N = params$N[i],
      alpha = params$alpha[i],
      lambda = c(30.302445, -97.731970),
      kappa = params$kappa[i],
      nu = 2,
      Omega =  (params$omega_scale[i] ^ 2) * diag(2)
    )
  })
  params$training_time[i] <- st["elapsed"]
  # Eval likelihood on test data
  params$likelihood[i] <- mean(dp_normal_deval(mod, x_test, nparticles = params$N[i]))
  pb$tick() # update progress bar
}

# Test-set error
saveRDS(params, "cv_run2.RDS")

# 
# # 4. Visualize best model on training and test set ========================================================
best <- which.min(params$likelihood)
mod <- dp_normal_mix(
  x_train,
  N = params$N[best],
  alpha = 10,
  lambda = c(30.302445, -97.731970),
  kappa = params$kappa[best],
  nu = 2,
  Omega =  (params$omega_scale[best] ^ 2) * diag(2)
)
ll_test <- dp_normal_deval(mod, x_test, nparticles = params$N[best])
ll_train <- dp_normal_deval(mod, x_train, nparticles = params$N[best])
# Run on a mesh grid
resol <- 50
xseq <- seq(30.13, 30.52, length.out = resol)
yseq <-  seq(-98.014, -97.58, length.out = resol)
mesh <- data.matrix(expand.grid(xseq, yseq, length.out = resol))
z <- dp_normal_deval(res2, mesh[ ,c(1,2)], nparticles = params$N[best])

# Training
spdata <- SpatialPointsDataFrame(x_train[ ,2:1], data.frame(density = ll_train))
resol <- 50
xseq <- seq(30.13, 30.52, length.out = resol)
yseq <-  seq(-98.014, -97.58, length.out = resol)
mesh <- data.matrix(expand.grid(xseq, yseq, length.out = resol))
z <- dp_normal_deval(res2, mesh[ ,c(1,2)], nparticles = params$N[best])
z <- matrix(z, resol, resol)
cl <- contourLines(xseq, yseq, z, nlevels = 1000)
pal <- colorNumeric("Spectral", domain = ll_train)
map <- leaflet(spdata) %>%
  addProviderTiles(providers$CartoDB.Positron)%>%
  setView(lat = 30.302445, lng = -97.731970, zoom = 11)
for (i in seq_along(cl)) {
  map <- map %>% addPolygons(
    lng = cl[[i]]$y,
    lat = cl[[i]]$x,
    fillColor = "red",
    fillOpacity = .03,
    stroke = FALSE)
}
map <- map %>%
  addCircleMarkers(
    radius = 3,
    stroke = FALSE,
    color = ~pal(density),
    fillOpacity = 0.8
  )  %>%
  addLegend(
    "bottomright",
    title = "Train density",
    values = ~density,
    pal = pal)
map

# Test
spdata <- SpatialPointsDataFrame(x_test[ ,2:1], data.frame(density = ll_test))
pal <- colorNumeric("Spectral", domain = ll_test)
map <- leaflet(spdata) %>%
  addProviderTiles(providers$CartoDB.Positron)%>%
  setView(lat = 30.302445, lng = -97.731970, zoom = 11)
for (i in seq_along(cl)) {
  map <- map %>% addPolygons(
    lng = cl[[i]]$y,
    lat = cl[[i]]$x,
    fillColor = "red",
    fillOpacity = .03,
    stroke = FALSE)
}
map <- map %>%
  addCircleMarkers(
    radius = 3,
    stroke = FALSE,
    color = ~pal(density),
    fillOpacity = 0.8
  )  %>%
  addLegend(
    "bottomright",
    title = "Test density",
    values = ~density,
    pal = pal)
map
