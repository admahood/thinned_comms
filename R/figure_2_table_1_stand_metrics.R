# stand structure changes
library(emmeans)
source("R/a_tc_data_prep.R")
# Seedlings =====================
seeds <- seedling_density |>
  tidyr::separate(new_visit_code, sep = "\\.", remove = T, into = c('plot', "phase")) |>
  dplyr::rename(value = seedlings_per_ha) 


ptsd <- seeds |>
  filter(phase == '01_Pre') |>
  dplyr::rename(ptv = value) |>
  dplyr::select(-phase)

diffs_sd <- seeds |>
  filter(phase != '01_Pre') |>
  left_join(ptsd) |>
  mutate(dv = value - ptv,
         site = str_sub(plot, 1,1),
         name = "seedling_density",
         new_visit_code = str_c(plot, '.', phase))
  

# plot-level stand metrics ====
cp_tree <- cp_tree |>
  dplyr::mutate(TreeStatus = ifelse(TreeStatus == "L", "live", "dead"))

ba_m2ha <- cp_tree |>
  filter(!is.na(BAPrismSize)) |> # filter live trees only?
  group_by(PlotCode, BAPrismSize, PlotTreatmentStatus, phase_adj) |>
  summarise(n=n()) |>
  ungroup() |>
  mutate(ba_ft2peracre = BAPrismSize * n,
         ba_m2pherha = ba_ft2peracre * 0.22956)|>
  group_by(PlotCode, PlotTreatmentStatus, phase_adj) |>
  summarise(ba_m2perha = sum(ba_m2pherha)) |>
  ungroup()

density <- cp_tree |>
  group_by(PlotCode, PlotTreatmentStatus, phase_adj, new_visit_code)|>
  summarise(tpa = sum(CalcCol_DensityContribution, na.rm = T)/0.404686) |>
  ungroup() |>
  mutate(site = str_sub(PlotCode, 1,1)) 

qmd <- cp_tree |>
  mutate(DBH = DBH*2.54) |>
  group_by(PlotCode, PlotTreatmentStatus, phase_adj, new_visit_code, TreatmentUnit)|>
  summarise(n=n(),
            qmd = sqrt(sum(DBH * DBH, na.rm=T)/n)) |>
  ungroup() |>
  dplyr::select(-n,-TreatmentUnit)

# side quest - qmd trend
# qmd_dat <- qmd |>
#   filter(PlotTreatmentStatus == "Treatment") |>
#   mutate(yst = case_when(phase_adj == "post0-1" ~ 1,
#                          phase_adj == "post04-5" ~ 5,
#                          phase_adj == "post10-11" ~ 10,
#                          phase_adj == "01_Pre" ~ -1))
# qmd_mod <- lm(qmd ~ yst, data = qmd_dat)
# qmd_mod1 <- lm(qmd ~ phase_adj, data = qmd_dat)
# 
# 
# summary(qmd_mod1)
# car::Anova(qmd_mod1)

fwd <- fuel_cover_by_type |>
  dplyr::filter(CodeType == 'FWD') |>
  dplyr::select(fwd = cover, new_visit_code)

means <- left_join(density, ba_m2ha) |> 
  left_join(qmd) |>
  left_join(fwd) |>
  pivot_longer(cols = c(tpa, qmd, ba_m2perha, fwd))

# table 1: emmeans =============================================================

allmeans <- means |>
  bind_rows(
    seeds |>
      dplyr::rename(PlotCode = plot, phase_adj = phase) |> 
      dplyr::mutate(name = "seeds",
                    site = str_sub(PlotCode,1,1),
                    new_visit_code = str_c(PlotCode, ".", phase_adj))
  )

ems <- list()
emstat <- list()
cc <- 1
for(i in unique(allmeans$name)) {
  fit <- lmerTest::lmer(value ~ phase_adj*PlotTreatmentStatus + (1|site/PlotCode), 
                 data = allmeans |> filter(name == i))
  emm <- emmeans(fit, ~ PlotTreatmentStatus | phase_adj) 
  ems[[cc]] <- emm |> 
    broom.mixed::tidy() |>
    mutate(response = i)
  emstat[[cc]] <- pairs(emm) |> broom::tidy() |> mutate(response = i)
  cc<- cc + 1
}

em_star <-
  bind_rows(emstat) |>
  dplyr::mutate(star = case_when(p.value >= 0.1 ~ "",
                                 p.value < 0.1 & p.value >= 0.05 ~ '.',
                                 p.value < 0.05 & p.value >= 0.01 ~ '*',
                                 p.value < 0.01 & p.value >= 0.001 ~ '**',
                                 p.value < 0.001 ~ "***")) |>
  dplyr::select(response, phase_adj, star) 

bind_rows(ems) |>
  left_join(em_star) |>
  dplyr::mutate(value = str_c(signif(estimate,3), " (", signif(std.error,3), ") ", star) |> trimws()) |>
  dplyr::select(-statistic, -df, -p.value, -estimate, -std.error, -star) |>
  tidyr::pivot_wider(names_from = phase_adj, values_from = value) |>
  janitor::clean_names() |>
  dplyr::select(response, trt = 1, 3, 4,5,6) |>
  write_csv('out/table1_stand_emms.csv')



# calculating difference from pretreatment for figure ==========================
ptv <- means |>
  dplyr::filter(phase_adj == '01_Pre') |>
  dplyr::select(PlotCode, name, ptv = value)

diffss <- means |>
  dplyr::filter(phase_adj != '01_Pre') |>
  left_join(ptv) |>
  mutate(dv = value - ptv) |>
  dplyr::rename(phase = phase_adj, plot = PlotCode) |>
  bind_rows(diffs_sd)

# contrasts ====================================================================
smstats <- list()
cc <- 1
for(i in unique(diffss$name)) {
  lmerTest::lmer(dv ~ phase*PlotTreatmentStatus + (1|site/plot), 
                 data = diffss |> filter(name == i)) -> fit
  em <- emmeans(fit, ~ PlotTreatmentStatus | phase)
  smstats[[cc]] <- pairs(em) |> broom::tidy() |> mutate(response = i)
  cc<- cc + 1
}

# t.tests, difference from pre treatment
result <- list()
cc <- 1
for(ph in unique(diffss$phase)){
  for(ct in unique(diffss$PlotTreatmentStatus)){
    for(rsp in unique(diffss$name)){
      result[[cc]] <- diffss |>
        filter(PlotTreatmentStatus == ct, phase == ph, name == rsp) |>
        pull(dv) |>
        wilcox.test() |>
        broom::tidy() |>
        mutate(phase = ph, PlotTreatmentStatus = ct, name = rsp)
      cc <- cc + 1
    }
  }
}
diff0ss <- bind_rows(result) |>
  mutate(sig0 = ifelse(p.value < 0.05, "Significant\nChange From\nPre-Treatment", 'ns')) |>
  dplyr::select(phase, name, PlotTreatmentStatus, sig0)



contrastss <- bind_rows(smstats) |>
  mutate(sig = ifelse(p.value < 0.05, "*", "")) |>
  dplyr::select(phase, name = response, sig) |>
  mutate(position = c(50,250,150,4,5,9,1,11,5,13,13,13, 5000, 5000, 5000))

diffss |>
  left_join(diff0ss) |>
  left_join(contrastss)|>
  dplyr::mutate(phase = case_when(phase == "post0-1" ~ "01",
                                  phase == "post04-5" ~ '05',
                                  phase == "post10-11" ~ "10"),
                name = case_when(name == "tpa" ~ "Trees/ha",
                                 name == 'fwd' ~ "Fine Woody Debris (%)",
                                 name == 'ba_m2perha' ~ 'Basal Area (m2/ha)',
                                 name == 'seedling_density' ~ "Seedlings/ha",
                                 name == 'qmd' ~ 'Quadratic Mean Diameter (cm)')) |>
ggplot(aes(x=phase, y=dv, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers = F, aes(color = sig0)) +
  geom_text(aes(label = sig, x= phase, y=position), size=12) +
  facet_wrap(~name, scales = 'free_y') +
  scale_fill_brewer(palette = "Set1") +
  scale_fill_manual(values = c('#FF00C5', '#00A884')) +
  scale_color_manual(values = c("grey", 'grey20')) +
  geom_hline(yintercept = 0, lty=2) +
  ylab("Change From Pre-Treatment") +
  xlab("Years Since Treatment") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(.9, .1),
        legend.justification = c(1,0))

ggsave('out/figure_1_stand_metric_changes.png', width = 9, height =6, bg = 'white')

diffss |> group_by(name, PlotTreatmentStatus) |> summarise(meanv = mean(dv))


