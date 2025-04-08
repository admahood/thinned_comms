# annual_z_scores <- read_csv("data/annual_z_scores.csv")

source('R/a_tc_data_prep.R')
library(terra)
library(ncdf4)
library(sf)
library(geomtextpath)
library(tidyverse)

lll <- st_read('data/terrain_by_plot.gpkg') |>
  st_transform(4326) |>
  filter(PlotCode %in% plot_visits_10y$PlotCode)

z <- list.files('data/big/terraclimate', full.names = T)

spring <- list()
# annual <- list()
for(i in 1:33){
  cwd_file <- z[i]
  rst <- terra::rast(cwd_file)[[c(3:6)]] |> sum()
  # ann <- terra::rast(cwd_file) |> sum()
  
  spring[[i]] <- lll |>
    dplyr::select(PlotCode) |>
    mutate(terra::extract(rst, lll, ID=F)) |>
    st_set_geometry(NULL) |>
    dplyr::rename(cwd = sum) |>
    dplyr::mutate(year = str_extract(cwd_file, '\\d{4}'))
  # 
  # annual[[i]] <- lll |>
  #   dplyr::select(PlotCode) |>
  #   mutate(terra::extract(ann, lll, ID=F)) |>
  #   st_set_geometry(NULL) |>
  #   dplyr::rename(cwd = sum) |>
  #   dplyr::mutate(year = str_extract(cwd_file, '\\d{4}'))
  
  print(i)
}

spring_norms <- bind_rows(spring) |>
  filter(year < 2021) |>
  group_by(PlotCode) |>
  summarise(mean = mean(cwd),
            sd = sd(cwd)) |>
  ungroup() 

spring_z <- bind_rows(spring) |>
  left_join(spring_norms) |>
  mutate(cwd_z = (cwd-mean)/sd)

# annual_norms <- bind_rows(annual) |>
#   filter(year < 2021) |>
#   group_by(PlotCode) |>
#   summarise(mean = mean(cwd),
#             sd = sd(cwd)) |>
#   ungroup() 
# 
# annual_z <- bind_rows(annual) |>
#   left_join(annual_norms) |>
#   mutate(cwd_z = (cwd-mean)/sd)


# bar plots ====================================================================
spring_z |> 
  filter(year>2009 & year < 2024) |> 
  mutate(site = ifelse(str_sub(PlotCode,1,1) == "E",
                       "Estes\nValley", "Phantom\nCreek")
         ) |>
  group_by(site, year) |>
  summarise(cwd_z = median(cwd_z)) |>
  ungroup() |>
  ggplot(aes(x=year, y=cwd_z, fill = site)) +
  ylab("Z-Scores") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = as.factor(c(2011, 2013, 2017, 2023, 2022)), lty=2) +
  geom_bar(stat = 'identity', color = 'black', position = 'dodge') +
  geom_vline(xintercept = as.factor(2012), linewidth = 2, color = 'firebrick', linetype =2) +
  ggtitle('D) Spring (March - June) Cumulative Climatic Water Deficit') +
  theme_bw() +
  theme(legend.title = element_blank(),
        plot.title = element_text(face = 'bold', size = 11),
        axis.title.x = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(color = 'black'))
  

ggsave('out/figure_1d_CWD_z.png', width = 7.5, height =3)
