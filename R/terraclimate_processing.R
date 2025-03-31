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

# result <- list()
# for(i in 1:33){
#   result[[i]] <- terra::rast(z[i])[[c(3:6)]] |> sum()
#   names(result[[i]]) <- str_c("cwd", 1990 + i)
#   print(i)
# }

result <- list()
for(i in 1:33){
  result[[i]] <- terra::rast(z[i]) |> sum()
  names(result[[i]]) <- str_c("cwd", 1990 + i)
  print(i)
}

brk <- terra::rast(result)
# monthly_brk <- terra::rast(z)

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
ggplot(annual_z_scores, aes(x=year, y=cwd_z, group = PlotCode)) +
  geom_line(alpha = 0.5) +
  geom_vline(xintercept = c(2011, 2013, 2017, 2023, 2022), lty=2) +
  geom_vline(xintercept = c(2012), linewidth = 2, color = 'red') +
  ggtitle('Annual CWD Z-Scores')

ggplot(annual_z_scores, aes(x=year, y=cwd, group = PlotCode)) +
  geom_line(alpha =0.5) +
  geom_vline(xintercept = c(2011, 2013, 2017, 2023, 2022), lty=2) +
  geom_vline(xintercept = c(2012), linewidth = 2, color = 'red') +
  ggtitle('Annual CWD (Cumulative)')

ggplot(annual_z_scores |> filter(year>2009), 
       aes(x=as.factor(year), y=cwd)) +
  geom_labelvline(data = yrs, aes(xintercept = as.factor(year)), 
                  lty=2, label = 'sample') +
  geomtextpath::geom_labelvline(xintercept = as.factor(2012), linewidth = 2, color = 'red', label = 'treatment') +
  geom_boxplot(alpha =0.5) +
  ggtitle('Annual CWD (Cumulative)')

# bar plots
dodge <- position_dodge(width=.75)
# control vs treatment,
ggplot(annual_z_scores |> 
         filter(year>2009, PlotCode %in% plot_visits_10y$PlotCode) |> 
         mutate(site = str_sub(PlotCode,1,1)) |> 
         mutate(site = ifelse(site == "E", "Estes", "Phantom")) |>
         group_by(year, site) |> summarise(mean_cwd = mean(cwd_z), sd_cwd = sd(cwd_z)) |>
         ungroup(), 
       aes(x=(year), y=mean_cwd, fill = site)) +
  geom_bar(stat = 'identity', position = dodge, color = 'black', width = .75) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(x=year, ymin=mean_cwd-(sd_cwd)*2, ymax = mean_cwd + (sd_cwd)*2, group = site), 
                position = dodge, width = .25) +
  geom_vline(xintercept = (c(2011, 2013, 2017, 2023, 2022)), lty=2) +
  geom_vline(xintercept = (2012), linewidth = 1.5, color = 'firebrick', lty=2) +
  ggtitle('D) Annual CWD Z-Scores') +
  theme_bw() +
  theme(axis.title = element_blank(),
        legend.position = c(0,1),
        legend.title = element_blank(),
        legend.justification  = c(0,1),
        legend.background = element_rect(fill = NA))
ggsave('out/climate_bar_plot.png', width =7, height =2.5, bg = 'white')

# # control vs treatment,
# ggplot(annual_z_scores |> 
#          filter(year>2009, PlotCode %in% plot_visits_10y$PlotCode) |> 
#          left_join(plot_visits_10y |> dplyr::select(PlotCode, PlotTreatmentStatus)) |>
#          group_by(year, PlotTreatmentStatus) |> 
#          mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Treatment", PlotTreatmentStatus)) |>
#          summarise(mean_cwd = mean(cwd_z), sd_cwd = sd(cwd_z)) |>
#          ungroup(), 
#        aes(x=(year), y=mean_cwd, fill = PlotTreatmentStatus)) +
#   geom_bar(stat = 'identity', position = dodge, color = 'black') +
#   geom_hline(yintercept = 0) +
#   # geom_errorbar(aes(x=year, ymin=mean_cwd-(sd_cwd)*2, ymax = mean_cwd + (sd_cwd)*2,
#   #                   group = PlotTreatmentStatus), 
#                 # position = dodge, width = .25) +
#   geom_vline(xintercept = (c(2011, 2013, 2017, 2023, 2022)), lty=2) +
#   geom_vline(xintercept = (2012), linewidth = 1.5, color = 'red', lty=2) +
#   ggtitle('Annual CWD Z-Scores') +
#   theme_bw() +
#   theme(axis.title = element_blank())

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
