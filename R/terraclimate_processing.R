# get terraclimate
source('R/a_tc_data_prep.R')
library(terra)
library(ncdf4)
library(sf)
library(tidyverse)

lll <- st_read('data/terrain_by_plot.gpkg') |>
  st_transform(4326)

z <- list.files('data/big/terraclimate/', full.names = T)[11:43]

result <- list()
for(i in 1:33){
  result[[i]] <- terra::rast(z[i])[[c(3:6)]] |> sum()
  names(result[[i]]) <- str_c("cwd", 1990 + i)
  print(i)
}

brk <- terra::rast(result)

norms <- lll |>
  dplyr::select(PlotCode) |>
  mutate(terra::extract(brk, lll, ID=F)) |>
  st_set_geometry(NULL) |>
  pivot_longer(-PlotCode) |>
  filter(str_sub(name,4,7) |> as.numeric() < 2021) |>
  group_by(PlotCode) |>
  summarise(mean = mean(value),
         sd = sd(value)) |>
  ungroup()

yrs <- plot_visits  |> 
  dplyr::mutate(year = lubridate::year(VisitDate)) |> 
  dplyr::select(PlotCode, year, new_visit_code)

z_scores <- lll |>
  dplyr::select(PlotCode) |>
  mutate(terra::extract(brk, lll, ID=F)) |>
  st_set_geometry(NULL) |>
  pivot_longer(-PlotCode) |>
  mutate(year = str_sub(name,4,7) |> as.numeric()) |>
  right_join(yrs) |>
  left_join(norms) |>
  dplyr::mutate(cwd_z = (value - mean)/sd) |>
  na.omit() |>
  dplyr::select(-name) |>
  dplyr::rename(cwd = value)

print(z_scores, n=199)

write_csv(z_scores,'data/terraclim_cwd_z.csv')
