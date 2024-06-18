source("R/a_tc_data_prep.R")
library(sf)
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
  summarise(ba_m2pherha = sum(ba_m2pherha),
            ba_ft2peracre = sum(ba_ft2peracre)) |>
  ungroup()

density <- cp_tree |>
  group_by(PlotCode, PlotTreatmentStatus, phase_adj, new_visit_code)|>
  summarise(tpa = sum(CalcCol_DensityContribution, na.rm = T)) |>
  ungroup() |>
  mutate(site = str_sub(PlotCode, 1,1)) 


# metrics of species diversity, etc ====
colz <- colnames(comm) |>
  as_tibble() |>
  dplyr::rename(SpeciesCode = value) |>
  left_join(sp_list |> dplyr::mutate(SpeciesCode = str_to_lower(SpeciesCode))) %>%
  dplyr::mutate(index = 1:nrow(.))

exotics <- colz |> dplyr::filter(NativityL48 == "exotic") |> pull(SpeciesCode)
natives <- colz |> dplyr::filter(NativityL48 == "native") |> pull(SpeciesCode)

div_indices <- vegan::diversity(comm[natives]) |>
  as_tibble(rownames = "new_visit_code") |>
  dplyr::rename(shannon_native = value) |>
  left_join(
    vegan::diversity(comm[exotics]) |>
      as_tibble(rownames = "new_visit_code") |>
      dplyr::rename(shannon_exotic = value) 
  ) |>
  left_join(
    vegan::specnumber(comm[exotics]) |>
      as_tibble(rownames = "new_visit_code") |>
      dplyr::rename(nspp_exotic = value) 
  )|>
  left_join(
    vegan::specnumber(comm[natives]) |>
      as_tibble(rownames = "new_visit_code") |>
      dplyr::rename(nspp_native = value) 
  )

# exotic vs native cover
species_list <- read_csv("data/species_list_alm_modified.csv") |> 
  unique() |> # one duplicate in species list
  dplyr::select(-notes)
cp <- cover_plot |>
  filter(!CodeType %in% non_plant_codes) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('plot', 'phase'), remove = F) |>
  mutate(site = str_sub(plot, 1,1))
spp <- cp |>
  filter(!CodeType %in% non_plant_codes) |>
  left_join(species_list |> mutate(SpeciesCode = str_to_lower(SpeciesCode))) |>
  left_join(plot_visits)

native_cover <- spp |>
  filter(NativityL48 == "native") |>
  group_by(plot, phase, PlotTreatmentStatus, site, trt_u_adj, new_visit_code) |>
  summarise(native_cover = sum(cover)) |>
  ungroup() |>
  tidyr::replace_na(list('trt_u_adj' = "PHA-1-2"))

exotic_cover <- spp |>
  filter(NativityL48 == "exotic") |>
  group_by(plot, phase, PlotTreatmentStatus, site, trt_u_adj, new_visit_code) |>
  summarise(exotic_cover = sum(cover)) |>
  ungroup() |>
  tidyr::replace_na(list('trt_u_adj' = "PHA-1-2"))

# joining it all together
tt <- st_read("data/terrain_by_plot.gpkg") |>
  st_set_geometry(NULL)

plot_level_metrics <- cp_tree |>
  mutate(site = str_sub(PlotCode, 1,1)) |>
  group_by(PlotCode, PlotTreatmentStatus, phase_adj, new_visit_code, site, TreatmentUnit)|>
  summarise(tpa = sum(CalcCol_DensityContribution, na.rm=T),
            n=n(),
            quadratic_mean_diameter = sqrt(sum(DBH * DBH, na.rm=T)/n)) |>
  left_join(ba_m2ha) |>
  left_join(div_indices)|> 
  left_join(tt) |>
  left_join(saplings_wide) |> 
  left_join(seedling_density) |>
  replace_na(list(sapling_density_per_acre = 0,
                  sapling_ba_ft_per_acre = 0)) |>
  left_join(native_cover) |>
  left_join(exotic_cover) |>
  replace_na(list(exotic_cover = 0)) |>
  mutate(no_exotics = 100 - exotic_cover,
         ec_binary = cbind(exotic_cover, no_exotics),
         log_saplingba = log(sapling_ba_ft_per_acre + 1),
         treated = case_when(PlotTreatmentStatus == "Control" ~ "not_treated",
                             PlotTreatmentStatus == "Treatment" & phase_adj == "01_Pre" ~ "not_treated",
                             PlotTreatmentStatus == "Treatment" & phase_adj != "01_Pre" ~ "treated") |> as.factor(),
         invaded = ifelse(exotic_cover >0, 1, 0)
         ) 

summary(plot_level_metrics)

write_csv(plot_level_metrics, "data/plot_level_data.csv")


# model exploration (Kyle Rodman's code from Springer et al 2023)
package.list <- c("here", "tidyverse", "glmmTMB", "DHARMa", "MuMIn", "ggplot2", 
                  "cowplot", "spaMM", "performance", "ncf", "splines")
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}
## Comparing with and without non-linear CWD term

alt1 <- glmmTMB(nspp_native ~   tpa + phase_adj  + slope + twi + hli  + exotic_cover+ 
                  (1|trt_u_adj/PlotCode), 
                data = plot_level_metrics, family = poisson(),
                na.action = "na.fail"); summary(alt1); performance::r2(alt1); diagnose(alt1); car::Anova(alt1)
alt2 <-  glmmTMB(nspp_native ~   tpa + phase_adj + treated + slope + twi + hli  + exotic_cover+ 
                   (1|trt_u_adj/PlotCode), 
                 data = plot_level_metrics, family = poisson(),
                 na.action = "na.fail"); summary(alt2); performance::r2(alt2); diagnose(alt2)
AICc(alt1, alt2)
dredge(alt1)
simResid <- simulateResiduals(alt1)
plot(simResid)
simResid <- simulateResiduals(alt2)
plot(simResid)


library(ggeffects)

ggpubr::ggarrange(
  plot(ggeffects::ggpredict(alt1, terms = c("tpa")),  show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("twi")),  show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("slope")),  show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("hli")),  show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("exotic_cover")),  show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("phase_adj")),  show_title=F),
ncol=4, nrow=2) |>
  ggsave(filename = 'out/native_richness_glmm_partials.png', width=9, height=5, bg="white")

plot(ggeffects::ggpredict(alt1, terms = c("tpa", "slope" )), show_residuals = TRUE, facets=TRUE, jitter = TRUE)
plot(ggeffects::ggpredict(alt1, terms = c("tpa", "hli")), show_residuals = TRUE, facets=TRUE, jitter = TRUE)
plot(ggeffects::ggpredict(alt1, terms = c("tpa", "phase_adj")), show_residuals = TRUE, facets=TRUE, jitter = TRUE)
plot(ggeffects::ggpredict(alt1, terms = c("tpa", "phase_adj", "slope")), show_residuals = TRUE, facets=TRUE,show_ci=T)

# native richness, bayesian 
library(brms)
brm1 <- brm(nspp_native ~   tpa + phase_adj + PlotTreatmentStatus + slope + twi + sapling_ba_ft_per_acre + hli + 
              (1|trt_u_adj/PlotCode), data = plot_level_metrics, family = 'poisson'); summary(brm1);# performance::r2(brm1); diagnose(brm1)

# probably cover of exotics is more important than species richness
alt2 <- glmmTMB(nspp_exotic ~ quadratic_mean_diameter + phase_adj + site + tpa + fa + fa_x_slope+
                  (1|TreatmentUnit/PlotCode), data = plot_level_metrics, family = poisson(),
                na.action = "na.fail");summary(alt2); performance::r2(alt2); diagnose(alt2)
dredge(alt2)
simResid <- simulateResiduals(alt2)
plot(simResid)

# probably cover of exotics is more important than species richness

ec1 <- glmmTMB(ec_binary ~ phase_adj + treated + tpa + site + slope + native_cover +
                  (1|trt_u_adj/PlotCode), data = plot_level_metrics, family = binomial(),
                na.action = "na.fail");summary(ec1); performance::r2(ec1); diagnose(ec1); car::Anova(ec1)
dredge(ec1) -> ddrd;ddrd

ggpubr::ggarrange(
  plot(ggeffects::ggpredict(alt1, terms = c("tpa")), show_residuals = TRUE, jitter = .08, show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("twi")), show_residuals = TRUE, jitter = .01, show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("slope")), show_residuals = TRUE, jitter = .08, show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("hli")), show_residuals = TRUE, show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("log_saplingba")), show_residuals = TRUE, show_title=F),
  plot(ggeffects::ggpredict(alt1, terms = c("PlotTreatmentStatus"))),
  plot(ggeffects::ggpredict(alt1, terms = c("phase_adj"))),
  ncol=2, nrow=4) |>
  ggsave(filename = 'out/native_richness_glmm_partials.png', width=7, height=10, bg="white")

simResid <- simulateResiduals(alt2)
plot(simResid)

inv1 <- glmmTMB(invaded ~ treated + phase_adj + site + slope + native_cover +
                  shannon_native + 
                 (1|trt_u_adj/PlotCode), data = plot_level_metrics, family = binomial(),
               na.action = "na.fail");summary(inv1); performance::r2(inv1); diagnose(inv1); car::Anova(inv1)

inv2 <- glmmTMB(invaded ~ tpa + ba_m2pherha + quadratic_mean_diameter + site + slope + native_cover +
                  shannon_native + twi + 
                  (1|trt_u_adj/PlotCode), data = plot_level_metrics, family = binomial(),
                na.action = "na.fail");summary(inv2); performance::r2(inv2); diagnose(inv2); car::Anova(inv2)
AIC(inv1, inv2)
simResid <- simulateResiduals(inv1)
plot(simResid)
simResid <- simulateResiduals(inv2)
plot(simResid)

ggpubr::ggarrange(
  plot(ggeffects::ggpredict(inv2, terms = c("tpa")), show_title=F),
  plot(ggeffects::ggpredict(inv2, terms = c("slope")), show_title=F),
  plot(ggeffects::ggpredict(inv2, terms = c("twi")), show_title=F),
  plot(ggeffects::ggpredict(inv2, terms = c("ba_m2pherha")), show_title=F),
  plot(ggeffects::ggpredict(inv2, terms = c("quadratic_mean_diameter")), show_title=F),
  plot(ggeffects::ggpredict(inv2, terms = c("site")), show_title=F),
  plot(ggeffects::ggpredict(inv2, terms = c("native_cover")), show_title=F),
  plot(ggeffects::ggpredict(inv2, terms = c("shannon_native")), show_title=F),
  ncol=4, nrow=2) |>
  ggsave(filename = 'out/invaded_glmm_partials.png', width=10, height=5, bg="white")

ggpubr::ggarrange(
  plot(ggeffects::ggpredict(inv1, terms = c("treated")), show_title=F),
  plot(ggeffects::ggpredict(inv1, terms = c("phase_adj")), show_title=F),
  plot(ggeffects::ggpredict(inv1, terms = c("site")), show_title=F),
  plot(ggeffects::ggpredict(inv1, terms = c("slope")), show_title=F),
  plot(ggeffects::ggpredict(inv1, terms = c("native_cover")), show_title=F),
  plot(ggeffects::ggpredict(inv1, terms = c("shannon_native")), show_title=F),
  ncol=3, nrow=2) |>
  ggsave(filename = 'out/invaded_cats_glmm_partials.png', width=7, height=5, bg="white")

