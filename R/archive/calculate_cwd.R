# calculate CWD at each plot
# we need 30 year normals, spring z-scores, annual z-scores
# site to download polaris data http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/v1.0/theta_s
# need soil available water capacity in the top 150-200 cm of the soil, in mm.
source('R/cwd_function_v3.R')
source('R/a_tc_data_prep.R')
library(sf)
library(terra)
library(tidyverse)
st_read('data/pha_eva_points.gpkg') -> d


ll <- st_read('data/terrain_by_plot.gpkg') |>
  st_transform(4326) |>
  st_coordinates() |>
  as_tibble() |>
  dplyr::rename(longitude = X, latitude = Y)


# polaris calc awc 200 cm ======================================================
depths <- list.files('data/polaris_theta_s/')
filez <- list.files('data/polaris_theta_s/0-5/')



for(ff in filez){
  x0_5 <- terra::rast(paste0('data/polaris_theta_s/0-5/',ff))
  x5_15 <- terra::rast(paste0('data/polaris_theta_s/5-15/',ff))
  x15_30 <- terra::rast(paste0('data/polaris_theta_s/15-30/',ff))
  x30_60 <- terra::rast(paste0('data/polaris_theta_s/30-60/',ff))
  x60_100 <- terra::rast(paste0('data/polaris_theta_s/60-100/',ff))
  x100_200 <- terra::rast(paste0('data/polaris_theta_s/100-200/',ff))
  
  awc_200 <- ((x0_5*50) + (x5_15*100) + (x15_30*150) + (x30_60 * 300) + (x60_100 * 400) + (x100_200 * 1000))/6
  terra::writeRaster(awc_200, filename = paste0('data/polaris_theta_s/awc200_', ff))
}

# extract polaris to points ====================================================
lapply(list.files('data/polaris_theta_s/', pattern = 'tif$', full.names = T), terra::rast) -> rasts
awc200 <- terra::mosaic(rasts[[1]], rasts[[2]], rasts[[3]])
# plot(awc200)

# 800m normals =================================================================
ppt_30y <- terra::rast('data/big/PRISM_ppt_30yr_normal_800mM4_annual_bil/PRISM_ppt_30yr_normal_800mM4_annual_bil.bil')
tmean_30y <- terra::rast('data/big/PRISM_tmean_30yr_normal_800mM5_annual_bil/PRISM_tmean_30yr_normal_800mM5_annual_bil.bil')

tt <- st_read('data/terrain_by_plot.gpkg') |>
  cbind(ll) |>
  st_transform(st_crs(ppt_30y)) %>%
  mutate(ppt_30y = terra::extract(y = ., x=ppt_30y, ID = F)[,1],
         tmean_30y = terra::extract(y = ., x=tmean_30y, ID = F)[,1],
         awc200 = terra::extract(y = ., x=awc200, ID = F)[,1])

# get monthly values for the whole time series for both vars in separate data frames
# monthly 4km prism ppt ============================================================
dirrr <- 'data/big/PRISM_ppt_stable_4kmM3_198101_202406_bil/'
ppt_files <- list.files(dirrr, 
                        pattern = "bil.bil$") |>
  as_tibble() |>
  dplyr::rename(filename = value) |>
  tidyr::separate(filename, sep = "_", remove = F,
                  into = c('product', 'variable', 'stable', 'res', 'yearmonth', 'end')) |>
  dplyr::mutate(year = str_sub(yearmonth, 1,4) |> as.numeric(),
                month = str_sub(yearmonth, 5,6),
                stem = dirrr,
                full_filename = paste0(dirrr, filename))

ppt_brick <- terra::rast(ppt_files$full_filename)

ppt_ext <- tt |>
  dplyr::select(PlotCode) |>
  mutate(terra::extract(ppt_brick, tt, ID=F)) |>
  st_set_geometry(NULL) |>
  pivot_longer(-PlotCode, values_to = 'ppt') |>
  tidyr::separate(name, sep = "_", remove = F,
                  into = c('product', 'variable', 'stable', 'res', 'yearmonth', 'end')) |>
  dplyr::mutate(year = str_sub(yearmonth, 1,4) |> as.numeric(),
                month = str_sub(yearmonth, 5,6) |> as.numeric()) |>
  dplyr::select(PlotCode, year, month, ppt)

tmean_files <- list.files('data/big/PRISM_tmean_stable_4kmM3_198101_202406_bil/',
                          pattern = 'bil.bil$', full.names = T)
tmean_brick <- terra::rast(tmean_files)
tmean_ext <- tt |>
  dplyr::select(PlotCode) |>
  mutate(terra::extract(tmean_brick, tt, ID=F)) |>
  st_set_geometry(NULL) |>
  pivot_longer(-PlotCode, values_to = 'tmean') |>
  tidyr::separate(name, sep = "_", remove = F,
                  into = c('product', 'variable', 'stable', 'res', 'yearmonth', 'end')) |>
  dplyr::mutate(year = str_sub(yearmonth, 1,4) |> as.numeric(),
                month = str_sub(yearmonth, 5,6) |> as.numeric()) |>
  dplyr::select(PlotCode, year, month, tmean)


ttd <- tt |>
  dplyr::select(PlotCode, slope, fa, latitude, awc200) |>
  st_set_geometry(NULL) |>
  right_join(tmean_ext) |>
  left_join(ppt_ext)


cwd <- cwd_function(site = ttd$PlotCode, slope = ttd$slope, latitude = ttd$latitude,
                    foldedaspect = ttd$fa, ppt = ttd$ppt, tmean = ttd$tmean, 
                    year = ttd$year, month = ttd$month, type = "annual",
                    soilawc = ttd$awc200)

cwd_normals_monthly <- 
  cwd |>
  as_tibble() |>
  filter(year > 1989 & year < 2020) |>
  group_by(site, month) |>
  summarise(monthly_normal = mean(cwd)) |>
  ungroup()

cwd_annual_normals <-
  cwd |>
  as_tibble() |>
  filter(year > 1989 & year < 2020) |>
  group_by(site, year) |>
  summarise(annual_cumulative = sum(cwd)) |>
  ungroup() 



cwd |> as_tibble() |>
  mutate(date = as.Date(paste0(year, "-", month, '-01'))) |>
  filter(year>1981) |>
  ggplot(aes(x=date, y=cwd, color = site)) +
  geom_line() 


# look at terraclim ====================
library(ncdf4)
tc_def <- terra::rast('data/TerraClimate_def_2023.nc')

plot(tc_def[[1]])

terra::extract(tc_def, tt) -> tcc

tc_norm <- terra::rast

