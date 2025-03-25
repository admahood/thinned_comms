# table 2: community metrics

# Table 1: basic metrics
source("R/a_tc_data_prep.R")
library(emmeans)

#todo get rid of the plots without all 4 timepoints

dummy_df <- cover_plot |>
  dplyr::select(-cover, -SpeciesCode) |>
  mutate(cover=0, SpeciesCode = "dummye", NativityL48 = 'exotic') |>
  bind_rows(cover_plot |>
              dplyr::select(new_visit_code) |>
              mutate(cover=0, SpeciesCode = "dummyn", NativityL48 = 'native'))

cp <- cover_plot |>
  filter(!CodeType %in% non_plant_codes) %>%
  left_join(sp_list |> mutate(SpeciesCode = str_to_lower(SpeciesCode))) |>
  tidyr::replace_na(list(NativityL48 = "Unknown")) |>
  bind_rows(dummy_df) |> 
  group_by(new_visit_code, NativityL48,PlotTreatmentStatus) |> 
  summarise(cover = sum(cover)) |>
  ungroup() |>
  group_by(new_visit_code,PlotTreatmentStatus) |>
  mutate(tc = sum(cover),
         rec = (cover/tc)*100) |>
  ungroup() |>
  filter(NativityL48 == "exotic") |>
  dplyr::select(-cover, -NativityL48) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('plot', 'phase')) |>
  mutate( site = str_sub(plot, 1,1))

cp |>
  group_by(phase) |>
  summarise(n=length(unique(plot)))

rp <- cover_plot |>
  filter(!CodeType %in% non_plant_codes) %>%
  left_join(sp_list |> mutate(SpeciesCode = str_to_lower(SpeciesCode))) |>
  tidyr::replace_na(list(NativityL48 = "Unknown")) |>
  group_by(new_visit_code, NativityL48, PlotTreatmentStatus) |>
  summarise(richness = length(unique(SpeciesCode))) |>
  ungroup() |> 
  group_by(new_visit_code,PlotTreatmentStatus) |>
  mutate(tr = sum(richness),
         rr = (richness/tr)*100) |>
  ungroup() |> 
  filter(NativityL48 != "Unknown") |>
  dplyr::select(-richness)|>
  tidyr::pivot_wider(names_from = NativityL48, values_from = rr, 
                     names_glue = "{NativityL48}_rr",
                     values_fill = 0) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('plot', 'phase')) |>
  mutate( site = str_sub(plot, 1,1))
print(rp, n=91)
rp |>
  group_by(phase) |>
  summarise(n=length(unique(plot)))

inv <- rp |>
  mutate(invaded = ifelse(exotic_rr == 0, 0,1),
         site = str_sub(plot, 1,1))

pct_inv <- inv |>
  group_by(phase, PlotTreatmentStatus) |>
  summarise(n_invaded = sum(invaded),
            n_plots = n()) |>
  ungroup() |>
  transmute(phase = phase, Treatment = PlotTreatmentStatus, 
            percent_invaded = round(n_invaded/n_plots * 100, 1))

pfgs <- cover_plot %>%
  filter(!CodeType %in% non_plant_codes) %>%
  filter(!CodeType %in% c("Fuel in Air", "UNK", "Tree")) |>
  left_join(sp_list |> mutate(SpeciesCode = str_to_lower(SpeciesCode))) |>
  tidyr::replace_na(list(NativityL48 = "Unknown")) |>
  group_by(new_visit_code) |>
  mutate(tc = sum(cover)) |>
  ungroup() |>
  group_by(CodeType, new_visit_code, PlotTreatmentStatus, tc) |>
  summarise(cover = sum(cover)) |>
  ungroup() |>
  mutate(rel_cover = (cover/tc)*100) |>
  tidyr::replace_na(list(Lifespan = '')) |>
  mutate(name = paste0(#Lifespan, "_", 
    CodeType)) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('plot', 'phase')) |>
  dplyr::mutate(site = ifelse(str_sub(plot,1,1) == "E", "E", "P"))



# graminoid cover
lmerTest::lmer(cover ~ phase*PlotTreatmentStatus + site + (1|site:plot), 
               data = pfgs |> filter(name == "Graminoid")) -> fit_gram
em_gr <- emmeans(fit_gram, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy(conf.int=T) |> 
  mutate(response = "Graminoid Cover")

# forb cover
lmerTest::lmer(cover ~ phase*PlotTreatmentStatus + site + (1|site:plot), 
               data = pfgs |> filter(name == "Forb")) -> fit_forb
em_fo <- emmeans(fit_forb, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy(conf.int=T) |> 
  mutate(response = "Forb Cover")

# shrub cover
lmerTest::lmer(cover ~ phase*PlotTreatmentStatus + site + (1|site:plot), 
               data = pfgs |> filter(name == "Shrub")) -> fit_sh
em_sh <- emmeans(fit_sh, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy(conf.int=T) |> 
  mutate(response = "Shrub Cover")

# relative ex richness
lmerTest::lmer(exotic_rr ~ phase*PlotTreatmentStatus + site + (1|site:plot), data = rp) -> fit1
emerr <- emmeans(fit1, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy(conf.int=T) |> 
  mutate(response = "Invasion Rate")
# total richness
lme4::lmer(tr ~ phase*PlotTreatmentStatus + site + (1|site:plot), data = rp) -> fit2
emtr <- emmeans(fit2, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy(conf.int=T) |> 
  mutate(response = "Total Richness")
# relative ex cover
lmerTest::lmer(rec ~ phase*PlotTreatmentStatus + site + (1|site:plot), data = cp) -> fit3
emerc <- emmeans(fit3, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy(conf.int=T) |> 
  mutate(response = "Invasion Impact")
# relative total cover
lmerTest::lmer(tc ~ phase*PlotTreatmentStatus + site + (1|site:plot), data = cp) -> fit4
emtc <- emmeans(fit4, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy(conf.int=T) |>
  mutate(response = "Total Cover")

lme4::glmer(invaded ~ phase*PlotTreatmentStatus + site + (1|site:plot), data = inv, family = 'binomial') -> fit5
emin <- emmeans(fit5, ~ PlotTreatmentStatus | phase, predict.type = 'response') |> 
  broom.mixed::tidy(conf.int=T) |>
  dplyr::rename(estimate = prob) |>
  mutate(response = "P(Invasion)")

emmeans(fit5, ~ phase|PlotTreatmentStatus, predict.type = 'response') |>
  pairs() |>
  broom::tidy() |>
  dplyr::select(Treatment = PlotTreatmentStatus, contrast, odds.ratio, std.error, statistic, adj.p.value) |>
  mutate_if(is.numeric, signif, 2) |>
  write.csv('out/trt_contrast_pinvasion.csv')
emmeans(fit5, ~ phase, predict.type = 'response') |> pairs()
emmeans(fit5, ~ PlotTreatmentStatus, predict.type = 'response') |> pairs()

lapply(list(fit1, fit2, fit3, fit4, fit5), performance::r2)


bind_rows(list(emin, emerr, emerc, emtc, emtr, em_gr, em_fo, em_sh)) |>
  dplyr::mutate(conf.low = ifelse(conf.low < 0, 0, conf.low),
    emm = paste0(round(estimate,2), " (", 
                             signif(conf.low,2), ", ",
                             signif(conf.high,2), ")")) |>
  dplyr::select(response, phase, Treatment = PlotTreatmentStatus, emm) |>
  pivot_wider(values_from = emm, names_from = phase) |>
  # left_join(inv) |>
  write_csv('out/table_2_emms.csv')
