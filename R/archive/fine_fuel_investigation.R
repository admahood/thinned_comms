# fuel loading investigation
source("R/a_tc_data_prep.R")
library(emmeans)

ground_cover <- cover_plot %>%
  filter(CodeType %in% non_plant_codes,
         !SpeciesCode %in% c('moss/lichen', 'rock', 'herb veg basal',# 'litter',
                             'coarse fuel in air', 'soil/gravel', 'woody basal')) |>
  dplyr::select(-CodeType) |>
  dplyr::mutate(SpeciesCode = ifelse(SpeciesCode == 'litter', 'litter/duff', SpeciesCode),
                site = str_sub(new_visit_code,1,1)) |>
  tidyr::separate(new_visit_code, into = c('plot', 'phase'), sep = '\\.', remove = F)

ggplot(ground_cover, aes(x=phase, y=cover, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers = F) +
  facet_grid(site~SpeciesCode, scales = 'free_y') 


# litter and duff loading examination ==========================================
sp_groundfuel |>
  ggplot(aes(y=Observation, x=phase_adj, fill = PlotTreatmentStatus)) +
  geom_boxplot() +
  facet_wrap(~Substrate) 


lut_bd<- c('Litter' = 2.37,
       'Duff' = 7.27)

load_df <- sp_groundfuel |>
  group_by(new_visit_code, Substrate) |>
  summarise(n_obs = length(unique(MeterID)),
            value = sum(Observation)/n_obs) |>
  ungroup() |>
  mutate(bd = lut_bd[Substrate],
         loading = value * bd) |>
  group_by(new_visit_code) |>
  summarise(loading = sum(loading)) |>
  ungroup() |>
  left_join(ground_cover |> filter(SpeciesCode == "litter/duff")) |>
  mutate(cover_adj_loading = loading * (cover/100)) 

load4 <- load_df |> filter(phase == 'post04-5')  |> dplyr::rename(load4 = loading)|> dplyr::select(plot, load4)

load_diff <- load_df |>
  dplyr::filter(phase != "post04-5") |>
  left_join(load4) |>
  mutate(diff = loading - load4) 

ggplot(load_diff) +
  geom_histogram(aes(x=diff))

load_df |>
  ggplot(aes(y=cover_adj_loading, x=phase, fill = PlotTreatmentStatus)) +
  geom_boxplot()


load_df|>
  ggplot(aes(x=loading, y=cover_adj_loading)) +
  geom_point() +
  geom_abline()

load_df |>
  ggplot(aes(x=cover, y=cover_adj_loading)) +
  geom_point() +
  geom_smooth()

# fine fuels examination =======================================================
sp_woodyfuel |>
  ggplot(aes(x=phase_adj, y=Calibrated, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers =F)+
  facet_wrap(~Measurement)

woodyfuel <- sp_woodyfuel |>
  dplyr::select(new_visit_code, Measurement, load = Calibrated) |>
  group_by(new_visit_code) |>
  summarise(load= sum(load)) |>
  ungroup() |>
  left_join(ground_cover |> filter(SpeciesCode |> str_sub(1,4) == "1/10")) 

lm(load ~ cover, data = woodyfuel) |> summary()

woodyfuel |>
  ggplot(aes(x=cover, y=load)) +
  geom_point()
