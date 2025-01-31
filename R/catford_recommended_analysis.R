#catford-recommended analysis
library(tidyverse)
library(sf)
library(ggpubr)
library(piecewiseSEM)
# plot-level stand metrics ====

plot_level_metrics <- read_csv("data/plot_level_data.csv") |>
  mutate_if(is.character, as.factor)

# models: explore covariates with glmmTMB, final model with brms
package.list <- c("here", "tidyverse", "glmmTMB", "DHARMa", "MuMIn", "ggplot2", 'brms',
                  "cowplot", "spaMM", "performance", "ncf", "splines", 'ggeffects')
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
sapply(package.list, library, character.only=T)


ggplot(plot_level_metrics, aes(x = tpa, fill = treated)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~site) +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_rect(fill=NA),
        legend.title = element_blank())
ggsave('out/tpa_avgs.png', width=5, height = 5, bg = 'white')

tpa_avgs <- plot_level_metrics |>
  group_by(treated, site) |>
  summarise(tpa = median(tpa),
            tpa_25 = quantile(tpa, 0.25),
            tpa_75 = quantile(tpa, 0.75)) |>
  ungroup()

# Exotic relative richness glmm ==========
err1 <- glmmTMB(exotic_relative_richness ~ tpa * phase_adj  + #PlotTreatmentStatus +
                  #ns(def_norm,2) + def_z_trt + 
                  native_cover + site +
                  (1|trt_u_adj/PlotCode), 
                data = plot_level_metrics, family = tweedie()); summary(err1); performance::r2(err1); diagnose(err1); car::Anova(err1)
p1<- plot(ggpredict(err1, terms = c('tpa', 'phase_adj', 'site'))) +
  ggtitle('Exotic Relative Richness: Marginal R2 = 0.559') 

err2 <- glmmTMB(exotic_relative_richness ~ treated * phase_adj  + 
                  site + native_cover + 
                  (1|trt_u_adj/PlotCode), 
                data = plot_level_metrics, family = tweedie(),
                na.action = "na.fail"); summary(err2); performance::r2(err2); diagnose(err2); car::Anova(err2)
p2 <- plot(ggpredict(err2, terms = c('phase_adj', 'treated', 'site'))) +
  ggtitle("Exotic Relative Richness: Marginal R2 = 0.486");p2
# plot(ggpredict(err2))
AIC(err1, err2)

ggarrange(p1, p2, nrow=2) 
ggsave(filename = "out/err_models.png", height = 8, width = 8, bg='white')


# br1 <- brm(exotic_relative_richness ~ tpa + phase_adj  + 
#               def_z_trt + native_cover +
#               (1|trt_u_adj/PlotCode),file = 'out/err1.rda', 
#             data = plot_level_metrics, family = 'zero_inflated_beta')
# plot(ggpredict(br1, terms = c('tpa', 'phase_adj')))
# 
# summary(br1); plot(conditional_effects(br1), ask = F); performance::r2(br1)
# 
# br2 <- brm(exotic_relative_richness ~ treated + phase_adj  + 
#              native_cover +
#              (1|trt_u_adj/PlotCode), 
#            data = plot_level_metrics, family = 'zero_inflated_beta')
# summary(br2); plot(conditional_effects(br2), ask = F); performance::r2(br2)



# relative cover exotic ==============

erc1 <- glmmTMB(exotic_relative_cover ~ phase_adj + treated + site +
                  #ns(def_norm,2) +#def_z_trt + 
                  (1|trt_u_adj/PlotCode), 
                data = plot_level_metrics, family = tweedie()); summary(erc1); performance::r2(erc1); diagnose(erc1); car::Anova(erc1)
p3 <- plot(ggpredict(erc1, terms = c('phase_adj', 'treated', "site")))

erc2 <- glmmTMB(exotic_relative_cover ~ phase_adj * quadratic_mean_diameter+ site +
                  (1|trt_u_adj/PlotCode), 
                data = plot_level_metrics |> mutate_if(is.character, as.factor), family = tweedie()); summary(erc2); performance::r2(erc2); diagnose(erc2); car::Anova(erc2)
# plot(ggpredict(erc2))
p4 <- plot(ggpredict(erc2, terms = c('quadratic_mean_diameter', 'phase_adj', 'site'))) + scale_y_sqrt()

AIC(erc1, erc2)

ggarrange(p3, p4, nrow=2) 
ggsave(filename = "out/erc_models.png", height = 8, width = 8, bg='white')

# 
glimpse(plot_level_metrics)
m1 <- glmmTMB(exotic_relative_cover ~ phase_adj + quadratic_mean_diameter  + tpa + ba_m2pherha +
                sapling_ba_ft_per_acre  + site +
                (1|trt_u_adj), data = plot_level_metrics, family = tweedie())
summary(m1); car::Anova(m1); performance::check_model(m1); performance::r2(m1)
plot(ggpredict(m1))
plot(ggpredict(m1, terms = c('treated', 'aet_trt')))

bm1 <- brm(exotic_relative_cover ~ quadratic_mean_diameter + 
             sapling_ba_ft_per_acre + phase_adj +
             (1|trt_u_adj), data = plot_level_metrics, family = 'zero_inflated_beta')
summary(bm1)
plot(conditional_effects(bm1), ask=F)

bm2 <- brm(exotic_relative_cover ~ phase_adj*PlotTreatmentStatus + sapling_ba_ft_per_acre + site +
             (1|trt_u_adj), data = plot_level_metrics, family = 'zero_inflated_beta')
summary(bm2); performance::r2(bm2); performance::check_model(bm2)
plot(conditional_effects(bm2), ask=F)
# hood et al reanalysis =============
# brm ============

# mb <- brm(rel_exotic_cover ~ year + ba_m2_ha*burned + (1|block/plot), data = catford, family = 'zero_inflated_beta')
# summary(mb)
# mb1 <- brm(rel_exotic_cover ~ year*trt + (1|block/plot), data = catford, family = 'zero_inflated_beta')
# summary(mb1); plot(conditional_effects(mb1), ask = FALSE)
# mb2 <- brm(rel_exotic_richness ~ year*trt + (1|block/plot), data = catford, family = 'zero_inflated_beta')
# summary(mb2); plot(conditional_effects(mb2), ask = FALSE)
# mb3 <- brm(invaded ~ year*trt + (1|block/plot), data = catford, family = 'bernoulli')
# summary(mb3);plot(conditional_effects(mb3), ask = FALSE)
# save(mb1, mb2, mb3, file = "data/hood_bayesian_models.rda")

load('data/hood_bayesian_models.rda'); read_csv("data/hood_catford.csv")
conditional_effects(mb3)[[3]] |> as.data.frame() |> dplyr::rename(value = 3) |> mutate(name = 'invaded') -> di
conditional_effects(mb2)[[3]] |> as.data.frame() |> dplyr::rename(value = 3) |> mutate(name = 'exotic_relative_richness') -> der
conditional_effects(mb1)[[3]] |> as.data.frame() |> dplyr::rename(value = 3) |> mutate(name = 'exotic_relative_cover') -> dec


bind_rows(di, der, dec) |>
  ggplot(aes(x=estimate__, y = trt, color = trt)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower__, xmax = upper__, y=trt)) +
  facet_wrap(name~year, scales = 'free_x', ncol=2) +
  ggthemes::theme_clean() +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = 'none',
        panel.background = element_rect(color = 'black'),
        axis.title = element_blank()) 
ggsave("out/hood_et_al_plot.png", bg='white', width =7, height =7)

# springer et al 

# explore springer 2023 data
readxl::read_xlsx("../springer_2023/Data/FieldData/LEARN_TransQuadSummaries.xlsx")
read_csv('../springer_2023/Data/Raw/beltTransects.csv')
read_csv('../springer_2023/Data/Raw/lineTransects.csv')
read_csv('../springer_2023/Data/Raw/speciesList.csv')
read_csv('../springer_2023/Data/Raw/treatmentYears.csv')
read_csv('../springer_2023/Data/plantDataWithSpatial.csv') |> glimpse()


springer_d <- read_csv('../springer_2023/Data/plantDataWithSpatial.csv') |>
  mutate(invaded = ifelse(rich_I > 0, 1, 0),
         yst = ifelse(is.na(Year_from_Initial_treatment), 0, Year_from_Initial_treatment),
         yst_c = case_when(yst == 0 ~ '0. pre',
                           yst %in% c(4,5) ~ "1. 4-5 post",
                           yst %in% c(8,9,10) ~ '2. 8-10 post',
                           yst >10 ~ '3. 16-20 post'),
         err = rich_I / rich_Tot,
         erc = freq_I / freq_Tot,
         treated = ifelse(Trt == 2, 'treated', 'control'))

# exotic relative richness of springer data ================
err3 <- glmmTMB(err ~ 
                  liveBA_ha + treated + yst +
                  #ns(avgCWD_1991_2020, 2) +
                  annCWD_z +
                  (1|Site/Block),
                data = springer_d, 
                family = tweedie()); summary(err3); performance::r2(err3); diagnose(err3); car::Anova(err3)
ps1 <- ggarrange(plot(ggpredict(err3, terms = c('treated', 'liveBA_ha'))),
                 plot(ggpredict(err3, terms = c('yst', 'treated'))), nrow=1, ncol=2)

err4 <- glmmTMB(err ~ 
                  liveBA_ha + 
                  yst_c + #treated +
                  # ns(avgCWD_1991_2020, 2) +
                  annCWD_z +
                  (1|Site/Block),
                data = springer_d, 
                family = tweedie()); summary(err4); performance::r2(err4); diagnose(err4); car::Anova(err4)
ps2 <- plot(ggpredict(err4, terms = c('yst_c','liveBA_ha')))

erc3 <- glmmTMB(erc ~ 
                  liveBA_ha + yst_c +
                  ns(avgCWD_1991_2020, 2) +
                  annCWD_z +
                  (1|Site),
                data = springer_d, 
                family = tweedie()); summary(erc3); performance::r2(erc3); diagnose(erc3); car::Anova(erc3)
ps3 <- plot(ggpredict(erc3, terms = c('yst_c','liveBA_ha')))

erc4 <- glmmTMB(erc ~ 
                  liveBA_ha + yst +
                  ns(avgCWD_1991_2020, 2) +
                  annCWD_z +
                  (1|Site),
                data = springer_d, 
                family = tweedie()); summary(erc4); performance::r2(erc4); diagnose(erc4); car::Anova(erc4)
ps4 <- ggarrange(plot(ggpredict(erc4, terms = c('liveBA_ha'))),
                 plot(ggpredict(erc4, terms = c('yst'))), nrow =2, ncol=1)

ggarrange(ps1, ps2, ps4, ps3, nrow=2, ncol=2, labels = 'auto')
ggsave('out/erc_err_springer.png', width=8, height=8, bg='white')

m <- glmmTMB(invaded ~ 
               liveBA_ha + 
               ns(avgCWD_1991_2020, 2) +
               annCWD_z +
               rich_N +
               yst +
               (1|Site),
             data = springer_d, 
             family = binomial(), na.action = 'na.fail') 
MuMIn::dredge(m)
summary(m); performance::check_model(m); car::Anova(m); performance::r2(m); diagnose(m)
ggeffects::ggpredict(m) |> plot()




m <- glmmTMB(invaded ~ 
               liveBA_ha + 
               ns(avgCWD_1991_2020, 2) +
               annCWD_z +
               rich_N +
               yst_c +
               (1|Site),
             data = springer_d, 
             family = binomial(), na.action = 'na.fail') 
MuMIn::dredge(m)
summary(m); performance::check_model(m); car::Anova(m); performance::r2(m); diagnose(m)
ggeffects::ggpredict(m) |> plot()

ggpubr::ggarrange(
  plot(ggeffects::ggpredict(m, terms = c("liveBA_ha")), show_title=F),
  plot(ggeffects::ggpredict(m, terms = c("avgCWD_1991_2020")), show_title=F),
  plot(ggeffects::ggpredict(m, terms = c("annCWD_z")), show_title=F),
  plot(ggeffects::ggpredict(m, terms = c("rich_N")), show_title=F),
  plot(ggeffects::ggpredict(m, terms = c("yst_c")), show_title=F),
  ncol=3, nrow=2) |>
  ggsave(filename = 'out/springer_invaded_cats_glmm_partials.png', width=7, height=5, bg="white")
