# setup
library(tidyverse)
library(brms)
library(spaMM)
library(glmmTMB)
library(randomForest)
library(ncf)
library(splines)
library(ggeffects)
library(ggpubr)

plot_level_metrics <- read_csv("data/plot_level_data.csv")  |>
  filter(!PlotCode %in% c('P1-2-T-615-01', 'P1-2-T-615-02 ')) |>
  dplyr::rename(ba_m2perha = ba_m2pherha) |>
  dplyr::mutate(erc = exotic_relative_cover/100,
                err = exotic_relative_richness/100) |>
  dplyr::mutate_if(is.character, as.factor) |>
  left_join(read_csv('data/terraclim_cwd_z.csv')) |>
  mutate(total_cover = native_cover + exotic_cover) |>
  tidyr::replace_na(list(seedlings_per_ha = 0)) |>
  filter(!phase_adj %in% c('01_Pre', 'post0-1')) |>
  mutate(yst = year - 2012,
         total_nspp = nspp_native + nspp_exotic)

glimpse(plot_level_metrics);summary(plot_level_metrics)

# exotic relative richness =====================================================
# rf
rf_err <- randomForest(err ~ ., data= plot_level_metrics |> 
                         dplyr::select(-contains('exotic'), -plot, -erc, -total_nspp,
                                       -invaded, -TreatmentUnit, -contains('native'))) 
varImpPlot(rf_err)
#glmm
mod_err <- glmmTMB(I(err+0.000001) ~ 
                     # tpa +
                     nspp_native +
                     # elevation +
                     ns(def_norm,2) +
                     cwd_z +
                     def_z_trt +
                     # seedlings_per_ha +
                     # fwd +
                     # trt_u_adj +
                     # quadratic_mean_diameter +
                     # ba_m2perha +
                     # total_cover +
                     (1|plot), 
                   data = plot_level_metrics,
                   family = beta_family(), REML = T, na.action = na.fail
)
# MuMIn::dredge(mod_err)
summary(mod_err); performance::r2(mod_err); performance::check_collinearity(mod_err); car::Anova(mod_err)

# exotic relative cover ========================================================
# rf
rf_erc <- randomForest(erc ~ ., data= plot_level_metrics |> 
                         dplyr::select(-contains('exotic'), -plot, -err, -contains('aet'), -contains('tmin'),
                                       -invaded, -TreatmentUnit, -total_nspp, -trt_u_adj,
                                       -contains('native'), -new_visit_code, 
                                       -PlotCode)) 
varImpPlot(rf_erc)
#glmm
mod_erc <- glmmTMB(I(erc + 0.000001) ~ 
                     # tpa +
                     # total_nspp +
                     # elevation +
                     # twi +
                     def_z_trt +
                     # cwd_z +
                     # def_norm +
                     ns(def_norm,2) +
                     # seedlings_per_ha +
                     # fwd +
                     # trt_u_adj +
                     # quadratic_mean_diameter +
                     # ba_m2perha +
                     # total_cover +
                     (1|plot), 
                   data = plot_level_metrics,
                   family = beta_family(), REML = F, na.action = na.fail
)
# MuMIn::dredge(mod_erc)
summary(mod_erc); performance::r2(mod_erc); performance::check_collinearity(mod_erc); car::Anova(mod_erc)


# p(invasion) ------------------------------------
# rf
rf_pi <- randomForest(invaded ~ ., data= plot_level_metrics |> 
                         dplyr::select(-contains('exotic'), -plot, -err, -erc,
                                       # -invaded, 
                                       -TreatmentUnit,
                                       -contains('native'), -new_visit_code, 
                                       -PlotCode) |> mutate(invaded = as.factor(invaded))) 
varImpPlot(rf_pi)
# glmm
mod_pi <- glmmTMB(invaded ~ 
                    # def_z_trt1 +
                    cwd_z +
                    def_z_trt +
                    # twi +
                    # fwd +
                    ns(def_norm,2) +
                    # def_z_trt +
                    # tpa +
                    # quadratic_mean_diameter +
                    # trt_u_adj +
                    total_cover +
                    # total_nspp +
                    # ba_m2perha +
                    # tpa +
                    (1|plot), #cores = 4, file = 'data/erc_mod',
                  data = plot_level_metrics, REML = T,
                  family = 'binomial', na.action = na.fail
)
# MuMIn::dredge(mod_pi)
# AIC(mod_pi, mod_pi1)
summary(mod_pi); performance::r2(mod_pi); car::Anova(mod_pi); performance::check_collinearity(mod_pi)

# native richness ==============================================================
rf_nr <- randomForest(nspp_native ~ ., data= plot_level_metrics |> 
                         dplyr::select(-contains('exotic'), -plot, -err, -erc, -total_nspp,
                                        -shannon_native, -native_cover, 
                                       -TreatmentUnit, -new_visit_code, 
                                       -PlotCode) |> mutate(invaded = as.factor(invaded))) 
varImpPlot(rf_nr)

mod_nr <- glmmTMB(nspp_native ~ 
                    ns(def_norm,2) +
                    def_z_trt +
                    # cwd_z +
                    hli +
                    # ba_m2perha +
                    tpa +
                    total_cover +
                    # fwd +
                    # twi +
                    # seedlings_per_ha +
                    # quadratic_mean_diameter +
                    # trt_u_adj+
                    (1|plot), #cores = 4, file = 'data/erc_mod',
                  data = plot_level_metrics, REML=T,
                  family = 'poisson', na.action = na.fail
)
#MuMIn::dredge(mod_nr)
summary(mod_nr); performance::r2(mod_nr); car::Anova(mod_nr); performance::check_collinearity(mod_nr)

# native cover =================================================================
rf_nc <- randomForest(native_cover ~ ., data= plot_level_metrics |> 
                        dplyr::select(-contains('exotic'), -plot, -err, -erc, -total_nspp,
                                      -shannon_native, -total_cover, -nspp_native, -contains('phase'),
                                      -TreatmentUnit, -new_visit_code, 
                                      -PlotCode) |> mutate(invaded = as.factor(invaded))) 
varImpPlot(rf_nc)

mod_nc <- glmmTMB(I(native_cover/100) ~ 
                    ns(def_norm,2)+
                    # cwd_z +
                    # elevation +
                    # ba_m2perha +
                    # nspp_native +
                    # seedlings_per_ha +
                    fwd +
                    # quadratic_mean_diameter +
                    # tpa +
                    # +
                    (1|plot), #cores = 4, file = 'data/erc_mod',
                  data = plot_level_metrics, REML = T, na.action = na.fail,
                  # map = list(theta = factor(NA)), 
                  # start = list(theta = log(10)),
                  family = beta_family()
)
# MuMIn::dredge(mod_nc)
summary(mod_nc); performance::r2(mod_nc); car::Anova(mod_nc); performance::check_collinearity(mod_nc)
plot(ggeffects::ggpredict(mod_nc))
