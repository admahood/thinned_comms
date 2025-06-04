# stuff that doesnt have all 4 years
source('R/a_tc_data_prep.R')
library(ggtext)
dl <- list()

dl$kf <- khr_fuels |>
  mutate(
    r1 = Diameter1/2, r2 = Diameter2/2,
    density_constant = ifelse(RottenSound == "Sound", 0.0125, 0.00935),
    area.multiplier = lut_khr_areas[PlotSize],
    volume_ft3 = (Length*12 * pi / 3 * (r1^2 + r1*r2 + r2^2)) * 0.000578704,
    loading_t_ac = volume_ft3 * density_constant * area.multiplier * 10) |> 
  group_by(new_visit_code, phase_adj, PlotTreatmentStatus, TreatmentUnit) |>
  reframe(loading_t_ac = sum(loading_t_ac, na.rm=T)) |>
  dplyr::select(new_visit_code, PlotTreatmentStatus, value = loading_t_ac)  |>
  dplyr::mutate(value = value*2.52)|>
  dplyr::mutate(name = "Coarse Wood Load (t ha<sup>-1</sup>)")

dl$saps <- cp_sapling |>
  group_by(new_visit_code, phase_adj, PlotTreatmentStatus) |>
  summarise(saplings_per_acre = sum(CalcCol_DensityContribution)) |>
  ungroup() |>
  filter(phase_adj %in% c('post10-11')) |>
  pivot_longer(cols = c(saplings_per_acre)) |>
  dplyr::select(-phase_adj) |>
  dplyr::mutate(name = 'Saplings ha<sup>-1</sup>', 
                value  = value * 2.47)

lut_bd<- c('Litter' = 2.37,
           'Duff' = 7.27)

dl$gf <- sp_groundfuel |>
  group_by(new_visit_code, Substrate,PlotTreatmentStatus) |>
  summarise(n_obs = length(unique(MeterID)),
            value = sum(Observation)/n_obs) |>
  ungroup() |>
  mutate(bd = lut_bd[Substrate],
         loading = value * bd) |>
  group_by(new_visit_code,PlotTreatmentStatus) |>
  summarise(value = sum(loading)) |>
  ungroup() |>
  dplyr::mutate(value = value*2.52) |>
  mutate(name = 'Litter/Duff Load (t ha<sup>-1</sup>)')

dl$wf <- sp_woodyfuel |>
  dplyr::select(new_visit_code, Measurement, load = Calibrated,PlotTreatmentStatus) |>
  group_by(new_visit_code,PlotTreatmentStatus) |>
  summarise(load= sum(load)) |>
  ungroup() |>
  dplyr::rename(value = load) |>
  dplyr::mutate(value = value*2.52,
    name = 'Fine Wood Load (t ha<sup>-1</sup>)')

dl$tcc <- tree_canopy_cover |>
  dplyr::filter(phase_adj == 'post10-11') |>
  group_by(new_visit_code,PlotTreatmentStatus) |>
  summarise(value = sum(PercentCover)) |>
  ungroup() |>
  dplyr::mutate(name = 'Tree Canopy Cover (%)')

dl$sh <- understory_heights |>
  filter(!is.na(Height_inches)) |>
  group_by(new_visit_code,PlotTreatmentStatus) |>
  summarise(value = mean(Height_inches)) |>
  ungroup() |>
  mutate(value = value * 2.54) |>
  mutate(name = 'Shrub Height (cm)')




d <- bind_rows(dl)  |>
  mutate(site = ifelse(str_sub(new_visit_code, 1,1) == "E", "Estes Valley", "Phantom Creek")) 

result <- list()
anovas <- list()
cc <- 1
for(i in unique(d$name)){
  mod <- lm(value ~ PlotTreatmentStatus + site,
                 data = d |> filter(name == i))
  result[[cc]] <- broom::tidy(mod) |> mutate(name = i)
  anovas[[cc]] <- car::Anova(mod) |> broom::tidy()  |> mutate(name = i)
  cc<- cc+1
}


bind_rows(anovas) |>
  filter(term != "Residuals") |>
  dplyr::select(name, term, statistic, p.value) |>
  mutate(term = ifelse(term == "site", "Site", "Treatment")) |>
  mutate_if(is.numeric, signif, 2) |>
  write_csv('out/tenyr_anovas.csv')

stars <- bind_rows(anovas) |>
  filter(term != "Residuals") |>
  dplyr::select(name, term, p.value) |>
  mutate(term = ifelse(term == "site", "Site", "Treatment")) |>
  filter(term == "Treatment") |>
  dplyr::mutate(star = case_when(p.value >= 0.1 ~ "",
                                 p.value < 0.1 & p.value >= 0.05 ~ '.',
                                 p.value < 0.05 & p.value >= 0.01 ~ '*',
                                 p.value < 0.01 & p.value >= 0.001 ~ '**',
                                 p.value < 0.001 ~ "***"))  |>
  dplyr::select(name, star) |>
  mutate(position = c(45, 5000,  45, 30, 70, 15))

p <- bind_rows(dl) |>
  left_join(stars) |>
  mutate(site = ifelse(str_sub(new_visit_code, 1,1) == "E", "Estes Valley", "Phantom Creek")) |>
  ggplot(aes(x=PlotTreatmentStatus, y=value, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers =F) +
  geom_text(aes(label = star, y=position), x=1.5, size=8) +
  facet_wrap(~name, scales ='free') +
  xlab("Site") +
  scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values = c('#FF00C5', '#00A884')) +
  ggtitle("Ten Years Post-Treatment") +
  theme_bw() +
  theme(#legend.position = c(.9,0.1),
        #legend.justification = c(1,0),
        legend.title = element_blank(),
        strip.text = element_markdown(),
        axis.title.y = element_blank()) +
  labs(caption = ". p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001") ;p

ggsave(p, filename = 'out/ten_year_only_variables.png', width = 7, height = 4.5)
