#' @import leaflet
#' @importFrom magrittr %>%

#' @title Spatial plots 
#' @export
spatial_plot <- function(
    model, 
    newx = NULL, 
    view_center = NULL, 
    zoom = 11,
    polygons = FALSE,
    contour_levels = 250,
    resol = 25, 
    nparticles = 50) 
  {
  # check inputs
  if (is.null(newx) && !polygons) stop("Polygons must be true if now newx is given")
  if (!inherits(newx, c("matrix", "Matrix"))) stop("newx must be a numeric matrix")
  if (is.null(newx) && is.null(view_center)) stop("if newx is null then a view center coordinate must be given")
  
  # Eval density and create spatial data object
  if (!is.null(newx)) {
    z <- dpn_eval(model, newx, nparticles)
    spdata <- sp::SpatialPointsDataFrame(newx, data.frame(density = z))  
  } else {
    spdata <- sp::SpatialPointsDataFrame(matrix(view_center, nrow = 1), data.frame(z = 0))  
  }
  
  # Set view center
  if (is.null(view_center)) {
    view_center <- as.numeric(apply(x, 2, mean))
  }
  
  # Create color pallete for newx
  pal <- colorNumeric("Spectral", domain = spdata$density)  
  
  # Begin map
  map <- leaflet(spdata) %>%
    addProviderTiles(leaflet::providers$CartoDB.Positron)%>%
    setView(lng = view_center[1], lat = view_center[2], zoom = 11)
  
  # Add polygons if necessary
  if (polygons) {
    # Evaluate density on a mesh
    xseq <- seq(min(x[ ,1]), max(x[ ,1]), length.out = resol)
    yseq <-  seq(min(x[ ,2]), max(x[ ,2]), length.out = resol)
    mesh <- data.matrix(expand.grid(xseq, yseq, length.out = resol))[ ,c(1,2)]
    z <- dpn_eval(model, mesh, nparticles)
    z <- matrix(z, resol, resol)
    
    # Create and add contour plots to map
    cl <- contourLines(xseq, yseq, z, nlevels = 250)
    for (i in seq_along(cl)) {
      map <- map %>% addPolygons(
        lng = cl[[i]]$x,
        lat = cl[[i]]$y,
        fillColor = "red",
        fillOpacity = .04,
        stroke = FALSE)
    }
  }
  
  # Add legend if newx was given add markers and circles
  if (!is.null(newx)) {
    map <- map %>% 
      addCircleMarkers(
        radius = 3,
        stroke = FALSE,
        color = ~pal(density),
        fillOpacity = 0.8
      ) %>% 
      addLegend(
        "bottomright",
        title = "Train density",
        values = ~density,
        pal = pal
      ) 
  }
  
  map
}
