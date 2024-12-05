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
  tidyr::separate(new_visit_code, into = c('plot', 'phase'), sep = '\\.')

ggplot(ground_cover, aes(x=phase, y=cover, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers = F) +
  facet_grid(site~SpeciesCode, scales = 'free_y') 

sp_woodyfuel |>
  ggplot(aes(x=phase_adj, y=Calibrated, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers =F)+
  facet_wrap(~Measurement)

sp_groundfuel |>
  ggplot(aes(y=Observation, x=phase_adj, fill = PlotTreatmentStatus)) +
  geom_boxplot() +
  facet_wrap(~Substrate)
