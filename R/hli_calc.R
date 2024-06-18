# heat load index function
library(terra)
dem <- terra::rast("data/dem_esv.tif")

get_latitude_raster <- function(dem){
  requireNamespace("terra")
  requireNamespace("sf")
  requireNamespace("dplyr")
  requireNamespace("terra")
  
  dem |>
    as.data.frame(xy=TRUE, cell=T) |>
    sf::st_as_sf(coords = c("x", "y"), crs =sf::st_crs(dem)) |>
    sf::st_transform(crs=4326) %>%
    dplyr::mutate(latitude = sf::st_coordinates(.)[,2]) |>
    sf::st_set_geometry(NULL) |>
    dplyr::left_join(as.data.frame(dem, cell=T, xy = TRUE)) |>
    dplyr::select(x,y, latitude) |>
    terra::rast(type = 'xyz', crs = terra::crs(dem, proj=T))
}

get_hli <- function(dem, out_path){
  #' calculates Heat Load Index (McCune and Keon 2002)
  #' Code adapted from K.C. Rodman
  requireNamespace("terra")
  # getting slope, converting to aspect, setting extreme highs and lows to fall 
  # within the boundaries of the HLI calculation
  slope <- terra::terrain(dem, v="slope")
  slope[slope > 60] <- 60
  slope[slope < 0] <- 0
  slope <- slope * 0.017453293
  aspect <- terra::terrain(dem, v="aspect", unit = "radians")
  latitude_radians <-  get_latitude_raster(dem) * 0.017453293
  cosine_latitude <- cos(latitude_radians)
  sine_latitude <- sin(latitude_radians)
  folded_aspect <- abs(pi - abs(aspect - (5*pi/4)))
  sine_slope <- sin(slope)
  cosine_slope <- cos(slope)
  cosine_fa <- cos(folded_aspect)
  return(
    0.339 + (0.808 * (cosine_latitude * cosine_slope)) - 
      (0.196 *(sine_latitude * sine_slope)) - (0.482 * (cosine_fa * sine_slope))
  )
}

