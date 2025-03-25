# brms models
library(tidyverse)
# library(brms)
# library(spaMM)
library(glmmTMB)
# library(ncf)
library(randomForest)
library(splines)
library(ggeffects)
library(ggpubr)
library(car)
library(broom.mixed)
library(broom)
library(cowplot)

forbs <- cover_plot %>%
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
  dplyr::mutate(site = ifelse(str_sub(plot,1,1) == "E", "E", "P")) |>
  dplyr::filter(name == "Forb") |>
  dplyr::select(PlotCode = plot, phase_adj = phase, forb_cover = cover)

plot_level_metrics <- read_csv("data/plot_level_data.csv")  |>
  filter(!PlotCode %in% c('P1-2-T-615-01', 'P1-2-T-615-02 ')) |>
  dplyr::rename(ba_m2perha = ba_m2pherha) |>
  dplyr::mutate(erc = exotic_relative_cover/100,
                err = exotic_relative_richness/100) |>
  dplyr::mutate_if(is.character, as.factor) |>
  left_join(read_csv('data/terraclim_cwd_z.csv')) |>
  mutate(total_cover = native_cover + exotic_cover) |>
  tidyr::replace_na(list(seedlings_per_ha = 0)) |>
  left_join(forbs)
  
# glimpse(plot_level_metrics);summary(plot_level_metrics)

# exotic relative richness =====================================================
# rf_err <- randomForest(err ~ ., data= plot_level_metrics |> 
#                          dplyr::select(-contains('exotic'), -plot, -erc, -phase,-phase_adj,
#                                        -invaded, -TreatmentUnit, -aspect, -year,
#                                        -contains('native'), -new_visit_code, -cwd,
#                                        -PlotCode)) 
# varImpPlot(rf_err)

mod_err <- glmmTMB(I(err+ + 0.000001) ~ 
                     tpa +
                     def_z_trt +
                     # ns(def_norm, 2) +
                     cwd_z +
                     # seedlings_per_ha +
                     # fwd +
                     # quadratic_mean_diameter +
                     # ba_m2perha +
                     # hli +
                     total_cover +
                     (1|plot), 
            data = plot_level_metrics,# map = list(theta = factor(NA)), 
             # start = list(theta = log(10)),
            family = beta_family(), REML =T, na.action = na.fail
              )
# MuMIn::dredge(mod_err)
# summary(mod_err); performance::r2(mod_err); performance::check_collinearity(mod_err); car::Anova(mod_err)


# graminoid cover ========================================================
# rf_gra <- randomForest(graminoid_cover ~ .,ntree=1000, 
#                        data= plot_level_metrics |> 
#                          dplyr::select(-contains('exotic'), -plot,  -phase,-phase_adj, -err, -total_cover,
#                                        -invaded, -TreatmentUnit, -aspect, -year, -treated,
#                                        -contains('native'), -new_visit_code, -cwd,
#                                        -PlotCode)) 
# varImpPlot(rf_gra)

mod_erc <- glmmTMB(I(graminoid_cover/100) ~ 
                     # erc +
                     # total_cover +
                     # def_norm +
                     ns(def_norm,2) +
                     # def_z_trt +
                     cwd_z +
                     # quadratic_mean_diameter +
                     # fwd +
                     hli +
                     tpa +
                     # seedlings_per_ha +
                     ba_m2perha +
                     (1|plot), REML = T, #map = list(theta = factor(NA)), 
                   # start = list(theta = log(10)), #cores = 4, file = 'data/erc_mod',
            data = plot_level_metrics, 
            na.action = na.fail,
            family = beta_family(), 
)

mod_forb <- glmmTMB(I(forb_cover/100) ~ 
                     # erc +
                     total_cover +
                     # def_norm +
                     ns(def_norm,2) +
                     # def_z_trt +
                     cwd_z +
                     # quadratic_mean_diameter +
                     # fwd +
                     hli +
                     tpa +
                     # seedlings_per_ha +
                     # ba_m2perha +
                     (1|plot), REML = T, #map = list(theta = factor(NA)), 
                   # start = list(theta = log(10)), #cores = 4, file = 'data/erc_mod',
                   data = plot_level_metrics, 
                   na.action = na.fail,
                   family = beta_family(), 
)
# MuMIn::dredge(mod_erc)
summary(mod_forb); performance::r2(mod_forb); performance::check_model(mod_forb); car::Anova(mod_forb); AIC(mod_forb)
# plot(ggeffects::ggpredict(mod_erc))
# 
# t0 <- Sys.time()
# bmod_erc <- brm(erc ~ 
#                      # err +
#                      total_cover +
#                      # def_norm +
#                      # s(def_norm) +
#                      def_z_trt +
#                      cwd_z +
#                      # quadratic_mean_diameter +
#                      # fwd +
#                      # hli +
#                      tpa +
#                      # seedlings_per_ha +
#                      # ba_m2perha +
#                      (1|trt_u_adj/plot), 
#                    cores = 4,# save_model = 'data/erc_mod',
#                    data = plot_level_metrics, iter = 3000,
#                    family = 'zero_inflated_beta', 
# );Sys.time()-t0
# 
# summary(bmod_erc)
# conditional_effects(bmod_erc)
# native richness ========================================================
# rf_nr <- randomForest(nspp_native ~ .,ntree=1000, 
#                        data= plot_level_metrics |> 
#                          dplyr::select(-contains('exotic'), -plot,  -phase,-phase_adj, -err, -erc, -total_cover,
#                                        -invaded, -TreatmentUnit, -aspect, -year, -treated, -n,
#                                        -new_visit_code, -cwd, -shannon_native, -native_cover,
#                                        -PlotCode)) 
# varImpPlot(rf_nr)
mod_nr <- glmmTMB(nspp_native ~ 
                  ns(def_norm,2) +
                  cwd_z +
                  hli +
                  # ba_m2perha +
                  tpa +
                  total_cover +
                  # seedlings_per_ha +
                  # quadratic_mean_diameter 
                  (1|plot), #cores = 4, file = 'data/erc_mod',
                data = plot_level_metrics, REML=T,
                family = 'poisson', na.action = na.fail
)
# MuMIn::dredge(mod_nr)

# summary(mod_nr); performance::r2(mod_nr); car::Anova(mod_nr)
# plot(ggeffects::ggpredict(mod_nr))


# native cover ========================================================
# rf_nc <- randomForest(native_cover ~ .,ntree=1000, 
#                       data= plot_level_metrics |> 
#                         dplyr::select(-contains('exotic'), -plot,  -phase,-phase_adj, -err, -erc, -total_cover,
#                                       -invaded, -TreatmentUnit, -aspect, -year, -treated, -n, -nspp_native,
#                                       -new_visit_code, -cwd, -shannon_native, -total_cover,
#                                       -PlotCode)) 
# varImpPlot(rf_nc)

mod_nc <- glmmTMB(I(native_cover/100) ~ 
                    cwd_z +
                    fwd +
                    quadratic_mean_diameter +
                    slope +
                    twi +
                    (1|plot), #cores = 4, file = 'data/erc_mod',
                  data = plot_level_metrics, na.action = na.fail, REML = T,
                  family = beta_family()
)
# MuMIn::dredge(mod_nc)
# summary(mod_nc); performance::r2(mod_nc); car::Anova(mod_nc); performance::check_collinearity(mod_nc); AIC(mod_nc)
# plot(ggeffects::ggpredict(mod_nc))

# invasion =================
# rf_pi <- randomForest(invaded ~ ., data= plot_level_metrics |> 
#                          dplyr::select(-contains('exotic'), -plot, -erc, -err, -phase,-phase_adj,
#                                        -TreatmentUnit, -aspect, -year, -treated,
#                                        -contains('native'), -new_visit_code, -cwd,
#                                        -PlotCode)|>
#                         mutate(invaded = as.factor(invaded))) 
# varImpPlot(rf_pi)

mod_pi <- glmmTMB(invaded ~ 
                    cwd_z +
                    ns(def_norm,2) +
                    def_z_trt +
                    total_cover +
                    # hli +
                    (1|plot), #cores = 4, file = 'data/erc_mod',
                  data = plot_level_metrics, REML = T,
                  family = 'binomial', na.action = na.fail
)
# MuMIn::dredge(mod_pi)
# AIC(mod_pi, mod_pi1)
# summary(mod_pi); performance::r2(mod_pi); car::Anova(mod_pi); performance::check_collinearity(mod_pi)
# plot(ggeffects::ggpredict(mod_pi, show_residuals = T))


# bind_rows(mod_pi |> broom.mixed::tidy() |> mutate(response = "P(Invasion)"),
#           mod_err |> broom.mixed::tidy(exponentiate = F) |> mutate(response = "Exotic Relative Richness"),
#           mod_erc |> broom.mixed::tidy(exponentiate = F) |> mutate(response = "Graminoid Cover"),
#           mod_nr |> broom.mixed::tidy(exponentiate = F) |> mutate(response = "Native Richness"),
#           mod_nc |> broom.mixed::tidy(exponentiate = F) |> mutate(response = "Native Cover")) |> 
#   dplyr::filter(effect != "ran_pars",
#                 !str_sub(term,1,5) %in% c("(Inte", "trt_u", 'phase', 'siteP', 'ns(de', 'hli')) |>
#   dplyr::select(-component, -group, -effect) |>
#   # dplyr::mutate(term = case_when(term == "ba_m2perha" ~ "Basal Area",
#   #                                term == "total_cover" ~ 'Total Cover',
#   #                                term == 'quadratic_mean_diameter' ~ 'QMD',
#   #                                term == 'cwd_z' ~ "CWD Z-Score",
#   #                                term == 'tpa:ba_m2perha' ~ 'tpa:ba_m2perha',
#   #                                term == "tpa" ~ 'Tree Density')) |>
#   # dplyr::mutate(sig = ifelse(p.value < 0.05, "*", "")) |> 
#   ggplot(#aes(color = sig)
#          ) +
#   geom_point(aes(x=estimate, y=term)) +
#   geom_segment(aes(x = estimate-std.error, xend = estimate+std.error,
#                    y = term)) +
#   geom_vline(xintercept = 0, lty = 2, lwd=.5) +
#   facet_wrap(~response, scales = 'free', ncol = 1) +
#   xlab("Estimate") +
#   ylab("Response Variable") +
#   theme_bw() 
# ggsave('out/figure_4_glmm_caterpillar.png', width = 3.5, height = 6, bg = 'white')

# naive models ===========================================================
mod_pic <- glmmTMB(invaded ~ #site + 
                     phase * PlotTreatmentStatus + 
                     (1|trt_u_adj/plot), 
                   data = plot_level_metrics,REML = T,
                  family = 'binomial')#; summary(mod_pic); AIC(mod_pic); car::Anova(mod_pic); performance::r2(mod_pic)

mod_ncc <- glmmTMB(I(native_cover/100) ~ #site + 
                     phase * PlotTreatmentStatus + 
                     (1|plot), 
                   data = plot_level_metrics,REML = T,
                  family = beta_family())
# summary(mod_ncc); performance::r2(mod_ncc)

mod_nrc <- glmmTMB(nspp_native ~#site + 
                     phase * PlotTreatmentStatus + 
                     (1|plot), 
                   data = plot_level_metrics,REML = T,
                  family = 'poisson')
# performance::r2(mod_nrc)
mod_errc <- glmmTMB(I(err + 0.000001) ~ #site + 
                      phase * PlotTreatmentStatus + 
                      (1|trt_u_adj/plot), 
                   data = plot_level_metrics,REML = T,
                   family = beta_family())
# summary(mod_errc); performance::r2(mod_errc)

mod_ercc <- glmmTMB(I(graminoid_cover/100) ~ #site + 
                      phase * PlotTreatmentStatus + 
                      (1|trt_u_adj/plot), 
                   data = plot_level_metrics,REML = T,
                   #map = list(theta = factor(NA)),
                   family = beta_family())

mod_forbc <- glmmTMB(I(forb_cover/100) ~ #site + 
                      phase * PlotTreatmentStatus + 
                      (1|plot), 
                    data = plot_level_metrics,REML = T,
                    family = beta_family())
# summary(mod_ercc);performance::r2(mod_forbc)

# bmod_ercc <- brm(erc ~ site + phase + PlotTreatmentStatus + #treated + 
#                       (1|plot), 
#                     data = plot_level_metrics, 
#                     family = 'zero_inflated_beta')
# summary(bmod_ercc); performance::r2(bmod_ercc)
# r2 table =====================================================================

aic_c <-  lapply(list(mod_pic, mod_errc, mod_ercc, mod_nrc, mod_ncc), MuMIn::AICc) |>
  unlist() |>
  as.data.frame() |>
  dplyr::mutate(response = c("P(Invasion)", "Exotic Relative Richness", 
                             "Graminoid Cover", 'Native Richness',
                             'Native Cover'))|>
  dplyr::select('naive' = 1,'response' = 2)


aic_tab <- lapply(list(mod_pi, mod_err, mod_erc, mod_nr, mod_nc), MuMIn::AICc) |>
  unlist() |>
  as.data.frame() |>
  dplyr::mutate(response = c("P(Invasion)", "Exotic Relative Richness", 
                             "Graminoid Cover", 'Native Richness',
                             'Native Cover'))|>
  dplyr::select('response' = 2,'stand_metrics' = 1) |>
  left_join(aic_c)|>
  mutate(delta_aic = naive - stand_metrics) |>
  mutate_if(is.numeric, round, 2)  |>
  dplyr::select(response, delta_aic)

r2c <-  lapply(list(mod_pic, mod_errc, mod_ercc, mod_nrc, mod_ncc, mod_forbc), 
               function(x)performance::r2(x)) |>
  bind_rows() |>
  
  dplyr::mutate(response = c("P(Invasion)", "Invasion Rate", 
                             "Graminoid Cover", 'Native Richness',
                             'Native Cover', "Forb Cover"))|>
  dplyr::select(3,'naive' = 2)
  
  
r2_tab <- lapply(list(mod_pi, mod_err, mod_erc, mod_nr, mod_nc, mod_forb), 
                 function(x)performance::r2(x)) |>
  bind_rows() |>
  dplyr::mutate(response = c("P(Invasion)", "Invasion Rate", 
                             "Graminoid Cover", 'Native Richness',
                             'Native Cover', "Forb Cover")) |>
  dplyr::select(3, informed = 2) |>
  left_join(r2c) |>
  mutate_if(is.numeric, signif, 2)

write_csv(r2_tab,'out/r2_table.csv')

# coefficients table ===========================================================
resps <-  c("Graminoid Cover", "Invasion Rate", "Native Cover", 
            "Native Richness", "P(Invasion)", "Forb Cover")
naive_models <- list(mod_ercc, mod_errc, mod_ncc, mod_nrc, mod_pic, mod_forbc)
ss_models <- list(mod_erc, mod_err, mod_nc, mod_nr, mod_pi, mod_forb)

lapply(naive_models,performance::model_performance) |> bind_rows() |>
  mutate(response = resps)
lapply(ss_models,performance::model_performance) |> bind_rows() |>
  mutate(response = resps)

ll <- list()
for(i in 1:length(naive_models)){
  ll[[i]] <- car::Anova(naive_models[[i]], test.statistic = "Chisq", type = "II") |>
    broom::tidy() %>%
    mutate(response = resps[i],
           pgroup = 'trt',
           pnum = letters[1:nrow(.)])
}
lll <- list()
for(i in 1:length(ss_models)){
  lll[[i]] <- car::Anova(ss_models[[i]], test.statistic = "Chisq", type = "II") |>
    broom::tidy() %>%
    mutate(response = resps[i],
           pgroup = 'ss',
           pnum = letters[1:nrow(.)])
}
bind_rows(ll, lll)|>
  dplyr::mutate(star = case_when(p.value >= 0.1 ~ "",
                                 p.value < 0.1 & p.value >= 0.05 ~ '.',
                                 p.value < 0.05 & p.value >= 0.01 ~ '*',
                                 p.value < 0.01 & p.value >= 0.001 ~ '**',
                                 p.value < 0.001 ~ "***"))|>
  pivot_wider(names_from = c(pgroup), values_from = c(term, star, statistic, p.value, df)) |>
  dplyr::mutate_if(is.numeric, round, 2) |>
  dplyr::mutate(trtstat = str_c(statistic_trt, star_trt),
                ssstat = str_c(statistic_ss, star_ss),
                term_trt = str_replace_all(term_trt, '\\:', ' x ') |>
                  str_replace_all("PlotTreatmentStatus", "Treatment") |>
                  str_replace_all('phase', "Timestep"),
                term_ss = term_ss |> str_replace_all('tpa', 'Tree Density') |>
                  str_replace_all('def_z_trt', 'CWD(z) Treatment Year') |>
                  str_replace_all('cwd_z', 'CWD(z) Sample Year') |>
                  str_replace_all("total_cover", "Total Cover") |>
                  str_replace_all('ns\\(def_norm, 2\\)', 'CWD normal') |>
                  str_replace_all('hli', 'Heat Load') |>
                  str_replace_all("twi", 'Topographic Wetness') |>
                  str_replace_all('ba_m2perha', "Basal Area") |>
                  str_replace_all('fwd', "Fine Woody Debris Cover") |>
                  str_replace_all('quadratic_mean_d', 'Quadratic Mean D')) |>
  dplyr::select(response,term_trt,trtstat, term_ss,ssstat) |>
  arrange(response) |> #print(n=99)
  write_csv('out/table5_glmms.csv')


# make an effects plot =========================================================
p_df0 <- bind_rows(
  lapply(ggeffects::ggpredict(mod_pi), as.data.frame) |> 
    bind_rows() |> mutate(response = "E. P(Invasion)"),
  lapply(ggeffects::ggpredict(mod_err), as.data.frame) |> 
    bind_rows() |> mutate(response = "F. Invasion Rate"),
  lapply(ggeffects::ggpredict(mod_erc), as.data.frame) |> 
    bind_rows() |> mutate(response = "D. Graminoid Cover"),
  lapply(ggeffects::ggpredict(mod_nr), as.data.frame) |> 
    bind_rows() |> mutate(response = "B. Native Richness"),
  lapply(ggeffects::ggpredict(mod_nc), as.data.frame) |> 
    bind_rows() |> mutate(response = "A. Native Cover"),
  lapply(ggeffects::ggpredict(mod_forb), as.data.frame) |> 
    bind_rows() |> mutate(response = "C. Forb Cover") |>
    filter(conf.high < 20)
)

p_df <- p_df0 |>
  mutate(type = case_when(
                  str_sub(group,1,3) %in% c('def', 'aet', 'cwd') ~ 'Climate',
                  group %in% c('hli', 'twi', 'slope') ~ "Topography",
                  str_sub(group,1,3) %in% c('tot', 'tpa', 'qua', 'fwd', "ba_") ~ "Stand\nStructure"),
         group = case_when(
           group == 'cwd_z' ~ ".CWD Z: Sample Year",
           group == 'def_norm' ~ ".CWD: 30-year Normal",
           group == 'hli' ~ "_Heat Load Index",
           group == 'total_cover' ~ "Total Veg Cover",
           group == 'tpa' ~ "Trees/ha",
           group == 'fwd' ~ "Fine Woody Debris",
           group == 'quadratic_mean_diameter' ~ "QMD",
           group == 'slope' ~ "_Slope",
           group == 'ba_m2perha' ~ 'Basal Area (m2/ha)',
           group == 'twi' ~ "_TWI",
           group == 'def_z_trt' ~ ".CWD Z: Year of Treatment"),
         predicted = ifelse(str_sub(response,1,8) %in% c("Invasion"), predicted * 100, predicted),
         conf.high = ifelse(str_sub(response,1,8) %in% c("Invasion"), conf.high  * 100, conf.high ),
         conf.low = ifelse(str_sub(response,1,8) %in% c("Invasion"), conf.low  * 100, conf.low )
         
  ) |>
  filter(group != "Total Veg Cover" | x < 80)
    

plot_row <- function(p_df, legend = F){
  p <- ggplot(p_df, aes(x=x, y = predicted)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = type)) +
    geom_line(lwd=1) +
    geom_line(aes(y=conf.low), lty = 1, lwd=.15) +
    geom_line(aes(y=conf.high), lty = 1, lwd=.15) +
    facet_wrap(~group, scales = 'free_x', nrow = 1) +
    scale_fill_brewer(palette = "Set2")+
    ylab(p_df$response) +
    theme_bw() +
    theme(axis.title.x = element_blank())
  if(legend){return(ggpubr::get_legend(p + theme(legend.title = element_blank())))}else(return(p + theme(legend.position = 'none')))
}



rpi <- plot_row(filter(p_df,response == "E. P(Invasion)"))
rer <- plot_row(filter(p_df,response == "F. Invasion Rate"))
rgc <- plot_row(filter(p_df,response == "D. Graminoid Cover"))
rnr <- plot_row(filter(p_df,response == "B. Native Richness"))
rnc <- plot_row(filter(p_df,response == "A. Native Cover"))
rf <- plot_row(filter(p_df,response == "C. Forb Cover"))

leg <- plot_row(filter(p_df,response == "B. Native Richness"), legend = TRUE) |> as_ggplot()

full <- cowplot::ggdraw(xlim =c(0,5), ylim =c(0,6)) +  
  draw_plot(rnc, 0, 5, 5, 1) +
  draw_plot(rnr, 0, 4, 5, 1) +  
  draw_plot(rf, 0, 3, 5, 1) +
  draw_plot(rgc, 0, 2, 5, 1) +
  draw_plot(rpi, 0, 1, 4.05, 1) +
  draw_plot(rer, 0, 0, 4.05, 1)+
  draw_plot(leg, 4, 1, 1, 1);full

ggsave(plot = full, filename = 'out/figure_5_effects_plots.png', height = 9, width = 9, bg='white')
