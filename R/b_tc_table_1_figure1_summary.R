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
         rr = (richness/tr)*100) |>
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

inv <- rp |>
  mutate(invaded = ifelse(exotic_rr == 0, 0,1)) |>
  group_by(phase, PlotTreatmentStatus) |>
  summarise(n_invaded = sum(invaded),
            n_plots = n()) |>
  ungroup() |>
  transmute(phase = phase, Treatment = PlotTreatmentStatus, 
            percent_invaded = round(n_invaded/n_plots * 100, 1))

# relative ex richness
lmerTest::lmer(exotic_rr ~ phase*PlotTreatmentStatus + (1|plot), data = rp) -> fit
emerr <- emmeans(fit, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy() |> 
  mutate(response = "Exotic Relative Richness")
# total richness
lme4::lmer(tr ~ phase*PlotTreatmentStatus + (1|plot), data = rp) -> fit
emtr <- emmeans(fit, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy() |> 
  mutate(response = "Total Richness")
# relative ex cover
lmerTest::lmer(rec ~ phase*PlotTreatmentStatus + (1|plot), data = cp) -> fit
emerc <- emmeans(fit, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy() |> 
  mutate(response = "Exotic Relative Cover")
# relative total cover
lmerTest::lmer(tc ~ phase*PlotTreatmentStatus + (1|plot), data = cp) -> fit
emtc <- emmeans(fit, ~ PlotTreatmentStatus | phase) |> 
  broom.mixed::tidy() |>
  mutate(response = "Total Cover")

bind_rows(list(emerr, emerc, emtc, emtr)) |>
  dplyr::mutate(emm = paste0(round(estimate,2), " (", round(std.error, 2), ")")) |>
  dplyr::select(response, phase, Treatment = PlotTreatmentStatus, emm) |>
  pivot_wider(values_from = emm, names_from = response) |>
  left_join(inv) |>
  write_csv('out/table_1_emms.csv')

table1 <- cp |>
  left_join(rp) |>
  dplyr::select(-plot) |>
  group_by(phase, PlotTreatmentStatus) |>
  summarise_all(mean) |>
  ungroup() |>
  mutate_if(is.numeric, signif, 3) |>
  dplyr::rename(treatment = PlotTreatmentStatus, 'Total Cover' = tc, 
                "Relative Exotic Cover" = rec,
                'Total Richness' = tr,
                'Relative Exotic Richness' = exotic_rr) |>
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
  filter(phase != '01_Pre') |>
  mutate(site = ifelse(str_sub(plot,1,1)=="E", "Estes Valley", "Phantom"))

# contrasts ====================================================================
lmerTest::lmer(dtc ~ phase*PlotTreatmentStatus + (1|site/plot), data = diffs) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
m1 <- pairs(em) |> broom::tidy() |> mutate(response = "Total Cover")
 
lmerTest::lmer(dtr ~ phase*PlotTreatmentStatus + (1|site/plot), data = diffs) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
m2 <- pairs(em) |> broom::tidy() |> mutate(response = "Total Richness")

lmerTest::lmer(derr ~ phase*PlotTreatmentStatus + (1|site/plot), data = diffs) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
m3 <- pairs(em) |> broom::tidy() |> mutate(response = "Exotic Relative Richness")

lmerTest::lmer(drec ~ phase*PlotTreatmentStatus + (1|site/plot), data = diffs) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
m4 <- pairs(em) |> broom::tidy() |> mutate(response = "Exotic Relative Cover")

# different than zero =============
responses <- dplyr::select(diffs, starts_with('d')) |> names()
is <- diffs |> dplyr::select(phase, PlotTreatmentStatus) |> unique()
result <- list()
cc <- 1
for(ph in unique(diffs$phase)){
  for(ct in unique(diffs$PlotTreatmentStatus)){
    for(rsp in unique(responses)){
      result[[cc]] <- diffs |>
        filter(PlotTreatmentStatus == ct, phase == ph) |>
        pull(rsp) |>
        wilcox.test() |>
        broom::tidy() |>
        mutate(phase = ph, PlotTreatmentStatus = ct, response = rsp)
      cc <- cc + 1
    }
  }
}
diff0s <- bind_rows(result) |>
  mutate(sig0 = ifelse(p.value < 0.05, "Significant\nChange From\nPre-Treatment", 'ns')) |>
  dplyr::select(phase, response, PlotTreatmentStatus, sig0)

# plotting =====================
contrasts <- bind_rows(m1, m2, m3, m4) |>
  dplyr::select(phase, response, p.value) |>
  mutate(sig = ifelse(p.value < 0.05, "*", ""),
         position = c(50,40,30,20,15,15,14,12,11,.2,.2,.2)) 

diffs |>
  dplyr::select(plot, phase, PlotTreatmentStatus, starts_with("d")) |>
  tidyr::pivot_longer(cols = c(dtc, drec, derr, dtr), 
                      names_to = 'response') |>
  left_join(diff0s) |>
  mutate(response = case_when(response == "dtc" ~ "Total Cover",
                                       response == "drec" ~ "Exotic Relative Cover",
                                       response == "derr" ~ "Exotic Relative Richness",
                                       response == "dtr" ~ "Total Richness")
    )|>
  left_join(contrasts) |>
  dplyr::mutate(phase = case_when(phase == "post0-1" ~ "01",
                                  phase == "post04-5" ~ '04',
                                  phase == "post10-11" ~ "10")) |>
  ggplot() +
  geom_boxplot(aes(x=phase, fill = PlotTreatmentStatus, y=value, color = sig0),
               outliers = F) +
  facet_wrap(~response, scales = 'free_y') +
  geom_hline(linetype = 2, yintercept =0) +
  geom_text(aes(label = sig, x= phase, y=position), size=12) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_manual(values = c("grey",'grey20')) +
  ylab("Change from Pre-Treatment") +
  xlab("Years Since Treatment")+
  labs(caption = "Stars indicate significant treatment differences.") +
  theme_bw() +
  theme(legend.title = element_blank()) 
ggsave('out/figure_x_treatment_diffs.png', width = 7, height = 6, bg = 'white')


# functional group comparisons =================================================
pfgs <- cover_plot %>%
  filter(!CodeType %in% non_plant_codes) %>%
  left_join(sp_list |> mutate(SpeciesCode = str_to_lower(SpeciesCode))) |>
  tidyr::replace_na(list(NativityL48 = "Unknown")) |>
  group_by(CodeType, new_visit_code, Lifespan, PlotTreatmentStatus) |>
  summarise(cover = sum(cover)) |>
  ungroup() |>
  tidyr::replace_na(list(Lifespan = '')) |>
  filter(!CodeType %in% c("Fuel in Air", "UNK")) |>
  mutate(name = paste0(#Lifespan, "_", 
    CodeType)) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('plot', 'phase')) |>
  dplyr::mutate(site = ifelse(str_sub(plot,1,1) == "E", "E", "P"))

pfgs |>
  ggplot() +
  geom_boxplot(aes(x=phase, fill = PlotTreatmentStatus, y = cover)) +
  facet_wrap(~name, scales = 'free_y') +
  theme_bw()

pre_treatment_pfg <- pfgs |>
  filter(phase == '01_Pre') |>
  dplyr::rename(ptc = cover) |>
  dplyr::select(-phase)

diffs_pfg <- pfgs |>
  left_join(pre_treatment_pfg) |>
  filter(phase != '01_Pre') |>
  mutate(dc = cover - ptc)

# contrasts ====================================================================
fgstats <- list()
cc <- 1
for(i in c("Forb", "Graminoid", "Shrub", "Tree")) {
lmerTest::lmer(dc ~ phase*PlotTreatmentStatus + (1|site/plot), 
               data = diffs_pfg |> filter(name == i)) -> fit
em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
fgstats[[cc]] <- pairs(em) |> broom::tidy() |> mutate(response = i)
cc<- cc + 1
}

# t.tests, difference from pre treatment
responses <- unique(diffs_pfg$name)
result <- list()
cc <- 1
for(ph in unique(diffs_pfg$phase)){
  for(ct in unique(diffs_pfg$PlotTreatmentStatus)){
    for(rsp in unique(responses)){
      result[[cc]] <- diffs_pfg |>
        filter(PlotTreatmentStatus == ct, phase == ph, name == rsp) |>
        pull(dc) |>
        wilcox.test() |>
        broom::tidy() |>
        mutate(phase = ph, PlotTreatmentStatus = ct, name = rsp)
      cc <- cc + 1
    }
  }
}
diff0s_pfg <- bind_rows(result) |>
  mutate(sig0 = ifelse(p.value < 0.05, "Significant\nChange From\nPre-Treatment", 'ns')) |>
  dplyr::select(phase, name, PlotTreatmentStatus, sig0)



contrasts_pfg <- bind_rows(fgstats) |>
  mutate(sig = ifelse(p.value < 0.05, "*", "")) |>
  dplyr::select(phase, name = response, sig) |>
  mutate(position = c(5,5,5,10,10,10,4,9,9, 9,9,9))

diffs_pfg |>
  left_join(contrasts_pfg) |>
  left_join(diff0s_pfg) |>
  filter(name != "Tree") |>
  dplyr::mutate(phase = case_when(phase == "post0-1" ~ "01",
                                  phase == "post04-5" ~ '04',
                                  phase == "post10-11" ~ "10")) |>
  ggplot() +
  geom_boxplot(aes(x=phase, fill = PlotTreatmentStatus, y = dc, color = sig0),
               outliers = F) +
  facet_wrap(~name, scales = 'free_y') +
  geom_text(aes(label = sig, x= phase, y=position), size=12) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_manual(values = c("grey", 'grey20')) +
  geom_hline(yintercept = 0, lty =2) +
  ylab("Change From Pre-Treatment") +
  xlab("Years Since Treatment") +
  theme_bw() +
  theme(legend.title = element_blank())
  
ggsave('out/figure_x_fg_treatment_diffs.png', width = 9, height = 4, bg = 'white')

# fwd etc ======================================================================
ground_cover <- cover_plot %>%
  filter(CodeType %in% non_plant_codes,
         !SpeciesCode %in% c('moss/lichen', 'rock', 'herb veg basal', 'litter',
                             'coarse fuel in air', 'soil/gravel', 'woody basal')) |>
  dplyr::select(-CodeType) |>
  dplyr::mutate(SpeciesCode = ifelse(SpeciesCode == 'litter', 'litter/duff', SpeciesCode)) |>
  tidyr::separate(new_visit_code, into = c('plot', 'phase'), sep = '\\.')

ggplot(ground_cover, aes(x=phase, y=cover, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers = F) +
  facet_wrap(~SpeciesCode, scales = 'free_y')

ptgc <- ground_cover |>
  filter(phase == '01_Pre') |>
  dplyr::rename(ptc = cover) |>
  dplyr::select(-phase)

diffs_gc <- ground_cover |>
  left_join(ptgc) |>
  filter(phase != '01_Pre') |>
  mutate(dc = cover - ptc,
         site = str_sub(plot, 1,1)) |>
  dplyr::rename(name = SpeciesCode)

# contrasts ====================================================================
gcstats <- list()
cc <- 1
for(i in unique(diffs_gc$name)) {
  lmerTest::lmer(dc ~ phase*PlotTreatmentStatus + (1|site/plot), 
                 data = diffs_gc |> filter(name == i)) -> fit
  em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
  gcstats[[cc]] <- pairs(em) |> broom::tidy() |> mutate(response = i)
  cc<- cc + 1
}

# t.tests, difference from pre treatment
result <- list()
cc <- 1
for(ph in unique(diffs_gc$phase)){
  for(ct in unique(diffs_gc$PlotTreatmentStatus)){
    for(rsp in unique(diffs_gc$name)){
      result[[cc]] <- diffs_gc |>
        filter(PlotTreatmentStatus == ct, phase == ph, name == rsp) |>
        pull(dc) |>
        wilcox.test() |>
        broom::tidy() |>
        mutate(phase = ph, PlotTreatmentStatus = ct, name = rsp)
      cc <- cc + 1
    }
  }
}
diff0s_gc <- bind_rows(result) |>
  mutate(sig0 = ifelse(p.value < 0.05, "Significant\nChange From\nPre-Treatment", 'ns')) |>
  dplyr::select(phase, name, PlotTreatmentStatus, sig0)



contrasts_gc <- bind_rows(gcstats) |>
  mutate(sig = ifelse(p.value < 0.05, "*", "")) |>
  dplyr::select(phase, name = response, sig) |>
  mutate(position = c(15,12,12,4,4,15,7,7,7))

diffs_gc |>
  left_join(contrasts_gc) |>
  left_join(diff0s_gc) |>
  dplyr::mutate(phase = case_when(phase == "post0-1" ~ "01",
                                  phase == "post04-5" ~ '04',
                                  phase == "post10-11" ~ "10")) |>
  ggplot() +
  geom_boxplot(aes(x=phase, fill = PlotTreatmentStatus, y = dc, color = sig0
                   ), outliers = F) +
  facet_wrap(~name, scales = 'free_y') +
  geom_text(aes(label = sig, x= phase, y=position), size=12) +
  scale_fill_brewer(palette = 'Set1') +
  scale_color_manual(values = c("grey", 'grey20')) +
  geom_hline(yintercept = 0, lty =2) +
  ylab("Change From Pre-Treatment") +
  xlab("Years Since Treatment") +
  theme_bw() +
  theme(legend.title = element_blank())
ggsave('out/figure_x_gc_treatment_diffs.png', width = 9, height = 4, bg = 'white')
