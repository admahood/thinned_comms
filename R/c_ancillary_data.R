# get some ancillary data
library(sf)
library(tidyverse)
library(elevatr)
library(terra)
library(topomicro)
files <- list.files("data", pattern = "Export", full.names = TRUE)

plys <- st_read("data/pha_eva_polys.gpkg")
pts <- st_read("data/pha_eva_points.gpkg")

plys |>
  filter(SiteCode == "ESV") |>
  ggplot() +
  geom_sf(aes(fill = TreatmentUnit))

plys |>
  filter(SiteCode == "PHA") |>
  ggplot() +
  geom_sf(aes(fill = TreatmentUnit)) +
  geom_sf(data=pts |>filter(SiteCode == "PHA"))

# annual precip
# vsurf on diversity, fg cover
# beta diversity on the steps
# for both sites, pre-treatment was done august, post-treatment in june
# waiting on heterogeneity
# herbaceous response as far as loading
# cover percentage by class (herbaceous vs shrub)

# response variables used by springer: cohen's d on trt phases ~ cover of natives, non-natives and northern-affinity spp
# GLMMS (binomial w/ ntrials=transect length): total understory species cover, richness and the proportion of northern affinity species
# GLMMS predictors: plot-level (i.e. 30-m spatial resolution) covariates of average climatic water deficit (CWD) 
# from 1991 to 2020, annual CWD in the sampling year (scaled to z-scores relative to the long-term mean and standard deviation at a plot)
# and topographic position index (i.e. TPI; an indicator of topographic exposure of a plot relative to its surroundings),
# as well as the unit-level covariate of treatment (i.e. untreated, 1–5-year post-treatment, 6–10-year post-treatment, or >10-year post-treatment)
# random intercept term of plot to account for repeated sampling within a plot over time, 
# Matérn spatial correlation term (Rousset & Ferdy, 2014) to account for spatial dependence among adjacent plots.

# we have cp_tree, seedling density to add as covariates
# maybe model difference from pre-treatment as gaussian, then have change in seedling/sapling/ba/tpa as predictors
# treatment method as well

plotvisit <- lapply(files, readxl::read_xlsx, "Export_TBL_PlotVisit") |>
  bind_rows() |>
  st_as_sf(coords = c("UTME", "UTMN"), crs =32613) |>
  mutate(project = ifelse(str_sub(PlotCode,1,1)=="E", "ESV", "PHA"))

dem_esv <- elevatr::get_aws_terrain(plotvisit[1,],z = 11, prj = st_crs(plotvisit))
  
dem_pha <- elevatr::get_aws_terrain(plotvisit[200,],z = 11, prj = st_crs(plotvisit))

# terra::writeRaster(dem_esv, 'data/dem_esv.tif')
# terra::writeRaster(dem_pha, 'data/dem_pha.tif')

dem_esv$slope <- terra::terrain(dem_esv[[1]], v = "slope")
dem_esv$aspect <- terra::terrain(dem_esv[[1]], v = "aspect")
dem_esv$twi <- topomicro::get_twi(dem_esv[[1]])$twi
dem_esv$fa <- topomicro::get_folded_aspect(dem_esv$aspect)
dem_esv$fa_x_slope <- sqrt(dem_esv$slope) * dem_esv$fa
dem_esv$hli <- topomicro::get_hli(dem_esv[[1]])
names(dem_esv)[1] <- "elevation"
# plot(dem_esv)

dem_pha$slope <- terra::terrain(dem_pha[[1]], v = "slope")
dem_pha$aspect <- terra::terrain(dem_pha[[1]], v = "aspect")
dem_pha$twi <- topomicro::get_twi(dem_pha[[1]])$twi
dem_pha$fa <- topomicro::get_folded_aspect(dem_pha$aspect)
dem_pha$fa_x_slope <- sqrt(dem_pha$slope) * dem_pha$fa
dem_pha$hli <- topomicro::get_hli(dem_pha[[1]])

names(dem_pha)[1] <- "elevation"
# plot(dem_pha)

esv_terrain <- plotvisit |>
  sf::st_transform(crs = st_crs(dem_esv)) |>
  dplyr::select(Elevation, VisitCode) |>
  dplyr::filter(str_sub(VisitCode, 1,3) == "ESV") %>%
  dplyr::mutate(terra::extract(dem_esv, .),
         Elevation_onsite = Elevation * 0.3048)
pha_terrain <- plotvisit |>
  sf::st_transform(crs = st_crs(dem_pha)) |>
  dplyr::select(Elevation, VisitCode) |>
  dplyr::filter(str_sub(VisitCode, 1,1) == "P") %>%
  dplyr::mutate(terra::extract(dem_pha, .),
                Elevation_onsite = Elevation * 0.3048)

terrain <- bind_rows(esv_terrain, pha_terrain) |>
  tidyr::separate(VisitCode, sep = "\\.", into = c("PlotCode", "phase")) |>
  dplyr::select(-phase, -ID, -Elevation, -Elevation_onsite) |>
  unique()

st_write(terrain, "data/terrain_by_plot.gpkg", delete_dsn = TRUE)

# visualizing to make sure crs's are lined up etc
esv_elv <- plotvisit |>
  dplyr::select(Elevation, project) |>
  filter(project == "ESV") %>%
  mutate(elv = terra::extract(dem_esv, ., ID=F) %>% unlist() %>% as.numeric(),
         Elevation = Elevation * 0.3048)
me<- mblm::mblm(elv ~ Elevation, data = esv_elv)

ggplot(esv_elv, aes(x=Elevation, y=elv)) +
  geom_point() +
  ggtitle("ESV") +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(y = predict(me)), color ="red")


pha_elv <- plotvisit |>
  dplyr::select(Elevation, project) |>
  filter(project == "PHA") %>%
  na.omit() %>%
  mutate(elv = terra::extract(dem_pha, ., ID=F) %>% unlist() %>% as.numeric(),
         Elevation = Elevation * 0.3048)

mp<- mblm::mblm(elv ~ Elevation, data = pha_elv)
ggplot(pha_elv, aes(x=Elevation, y=elv)) +
  geom_point() +
  ggtitle("PHA")+
  geom_abline(slope = 1, intercept = 0)+
  geom_line(aes(y = predict(mp)), color ="red")

glimpse(plotvisit)
plot(dem_esv); plot(plotvisit[0], add=T)
ggplot(plotvisit |> filter(project == "ESV")) +
  geom_sf(aes(color = PlotTreatmentStatus)) 

plot(dem_pha); plot(plotvisit[0], add=T)
ggplot(plotvisit |> filter(project == "PHA")) +
geom_sf(aes(color = PlotTreatmentStatus)) 

 
# usda plant trait data

