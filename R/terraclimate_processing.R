# get terraclimate
 source('R/a_tc_data_prep.R')
library(terra)
library(ncdf4)
library(sf)
library(geomtextpath)
library(tidyverse)

lll <- st_read('data/terrain_by_plot.gpkg') |>
  st_transform(4326)

z <- list.files('data/big/terraclimate/', full.names = T)

result <- list()
for(i in 1:33){
  result[[i]] <- terra::rast(z[i])[[c(3:6)]] |> sum()
  names(result[[i]]) <- str_c("cwd", 1990 + i)
  print(i)
}

brk <- terra::rast(result)
monthly_brk <- terra::rast(z)

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

# for climate plot ==========
annual_z_scores <- lll |>
  dplyr::select(PlotCode) |>
  mutate(terra::extract(brk, lll, ID=F)) |>
  st_set_geometry(NULL) |>
  pivot_longer(-PlotCode) |>
  mutate(year = str_sub(name,4,7) |> as.numeric()) |>
  left_join(norms) |>
  dplyr::mutate(cwd_z = (value - mean)/sd) |>
  na.omit() |>
  dplyr::select(-name) |>
  dplyr::rename(cwd = value)

# control vs treatment in plot
# precip bar plots
# facet by site
ggplot(z_scores, aes(x=year, y=cwd_z, group = PlotCode)) +
  geom_line(alpha = 0.5) +
  geom_vline(xintercept = c(2011, 2013, 2017, 2023, 2022), lty=2) +
  geom_vline(xintercept = c(2012), linewidth = 2, color = 'red') +
  ggtitle('Annual CWD Z-Scores')

ggplot(z_scores, aes(x=year, y=cwd, group = PlotCode)) +
  geom_line(alpha =0.5) +
  geom_vline(xintercept = c(2011, 2013, 2017, 2023, 2022), lty=2) +
  geom_vline(xintercept = c(2012), linewidth = 2, color = 'red') +
  ggtitle('Annual CWD (Cumulative)')

ggplot(z_scores |> filter(year>2009), 
       aes(x=as.factor(year), y=cwd)) +
  geom_labelvline(data = yrs, aes(xintercept = as.factor(year)), 
                  lty=2, label = 'sample') +
  geomtextpath::geom_labelvline(xintercept = as.factor(2012), linewidth = 2, color = 'red', label = 'treatment') +
  geom_boxplot(alpha =0.5) +
  ggtitle('Annual CWD (Cumulative)')

# bar plots
ggplot(z_scores |> filter(year>2009) |> mutate(site = str_sub(PlotCode,1,1)), 
       aes(x=as.factor(year), y=cwd_z, fill = site)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = as.factor(c(2011, 2013, 2017, 2023, 2022)), lty=2) +
  geom_vline(xintercept = as.factor(2012), linewidth = 2, color = 'red') +
  geom_boxplot() +
  ggtitle('Annual CWD Z-Scores')

# for glmms ====
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
