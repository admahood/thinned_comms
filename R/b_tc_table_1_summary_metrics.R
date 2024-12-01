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
         rec = cover/tc) |>
  ungroup() |>
  filter(NativityL48 == "exotic") |>
  dplyr::select(-cover, -NativityL48) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('plot', 'phase'))

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
         rr = richness/tr) |>
  ungroup() |> 
  filter(NativityL48 != "Unknown") |>
  dplyr::select(-richness)|>
  tidyr::pivot_wider(names_from = NativityL48, values_from = rr, 
                     names_glue = "{NativityL48}_rr",
                     values_fill = 0) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('plot', 'phase'))
print(rp, n=91)
rp |>
  group_by(phase) |>
  summarise(n=length(unique(plot)))

# relative ex richness
lmerTest::lmer(exotic_rr ~ phase*PlotTreatmentStatus + (1|plot), data = rp) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
pairs(em)

# total richness
lme4::glmer(tr ~ phase*PlotTreatmentStatus + (1|plot), data = rp, family = 'poisson') -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
pairs(em)

# relative ex cover
lmerTest::lmer(rec ~ phase*PlotTreatmentStatus + (1|plot), data = cp) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
pairs(em)

# relative total cover
lmerTest::lmer(tc ~ phase*PlotTreatmentStatus + (1|plot), data = cp) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
pairs(em)

table1 <- cp |>
  left_join(rp) |>
  dplyr::select(-plot) |>
  group_by(phase, PlotTreatmentStatus) |>
  summarise_all(mean) |>
  ungroup() |>
  mutate_if(is.numeric, round, 3) |>
write_csv("out/table_1_summary.csv")

pre_treatment <- cp |>
  left_join(rp)  |>
  filter(phase == '01_Pre') |>
  dplyr::rename(ptc = tc, prec = rec, ptr =tr, perr = exotic_rr) |>
  dplyr::select(-native_rr, -phase)

diffs <- cp |>
  left_join(rp) |>
  left_join(pre_treatment) |>
  mutate(dtc = tc - ptc,
         drec = rec - prec,
         dtr = tr - ptr,
         derr = exotic_rr -perr) |>
  filter(phase != '01_Pre')

lmerTest::lmer(dtc ~ phase*PlotTreatmentStatus + (1|plot), data = diffs) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
pairs(em)

lmerTest::lmer(dtr ~ phase*PlotTreatmentStatus + (1|plot), data = diffs) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
pairs(em)

lmerTest::lmer(derr ~ phase*PlotTreatmentStatus + (1|plot), data = diffs) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
pairs(em)

lmerTest::lmer(drec ~ phase*PlotTreatmentStatus + (1|plot), data = diffs) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
pairs(em)

diffs |>
  dplyr::select(plot, phase, PlotTreatmentStatus, starts_with("d"))
  pivot_longer(cols = c(dtc, drec, derr, dtr)) |>
    dplyr::mutate(phase = case_when(phase == "post0-1" ~ "01",
                                    phase == "post04-5" ~ '04',
                                    phase == "post10-11" ~ "10")
    ) |>
  ggplot() +
  geom_boxplot(aes(x=phase, fill = PlotTreatmentStatus, y=value)) +
  facet_wrap(~name, scales = 'free_y') +
  geom_hline(lintype = 2) +
  ylab("Change from Pre-Treatment") +
  xlab("Years Since Treatment")

