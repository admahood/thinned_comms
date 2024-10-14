library(tidyverse)
library(lme4)
library(glmmTMB)
library(ggeffects)

# read_csv("data/RDS-2023-0063/Data/fuels.csv")
# read_csv("data/RDS-2023-0063/Data/saplings.csv")
# read_csv("data/RDS-2023-0063/Data/seedlings.csv")
# read_csv("data/RDS-2023-0063/Data/trees.csv")
read_csv("data/RDS-2023-0063/Data/understory.csv") |> filter(ORIGIN == "Exotic", COVER>0) |> 
  group_by(LFCODE, GENUS, YEAR) |>
  summarise(n=n()) |>
  ungroup() |>
  pivot_wider(names_from = YEAR, values_from = n, values_fill = 0) |>
  print(n=99)
read_csv("data/RDS-2023-0063/Data/_variable_descriptions.csv") |> print(n=99)

read_csv("data/RDS-2023-0063/Data/understory.csv") |>
  dplyr::select(1,2,3) |>
  unique() |>
  group_by(BLOCK, TRT) |>
  summarise(n = n())

ba_02 <- read_csv("data/RDS-2023-0063/Data/trees.csv") |>
  janitor::clean_names() |>
  filter(stat02 == 0) |>
  group_by(block, trt, plot) |>
  summarise(ba_m2_ha = ((sum(dbh_01, na.rm=T)^2)/9)* (pi/40000)) |>
  mutate(plot = str_c(block, "_", plot),
         block = as.factor(block),
         year = '2002')#9ha plots

ba_20 <- read_csv("data/RDS-2023-0063/Data/trees.csv") |>
  janitor::clean_names() |>
  filter(stat20 == 0) |>
  group_by(block, trt, plot) |>
  summarise(ba_m2_ha = ((sum(dbh_20, na.rm=T)^2)/9)* (pi/40000)) |>
  ungroup() |>
  mutate(plot = str_c(block, "_", plot),
         block = as.factor(block),
         year = '2022') |>
  bind_rows(ba_02)
  #9ha plots

saplings <- read_csv("data/RDS-2023-0063/Data/saplings.csv") |>
  janitor::clean_names() |>
  group_by(block, plot, trt, year) |>
  summarise(sapling_density_ha = sum(density_ha)) |>
  ungroup() |>
  mutate(plot = str_c(block, "_", plot),
         block = as.factor(block),
         year = as.character(year))

seedlings <- read_csv("data/RDS-2023-0063/Data/seedlings.csv") |>
  janitor::clean_names() |>
  pivot_longer(cols = c('density_ha2002', 'density_ha2022'), names_to = "year", values_to = "density_ha") |>
  mutate(year = str_replace_all(year, 'density_ha', '')) |>
  group_by(block, trt, plot, year) |>
  summarise(seedling_density_ha = sum(density_ha, na.rm = T)) |>
  ungroup() |>
  mutate(plot = str_c(block, "_", plot),
    block = as.factor(block))


exotic_richness <- read_csv("data/RDS-2023-0063/Data/understory.csv") |>
  janitor::clean_names() |>
  group_by(block, trt, year, plot) |>
  filter(origin == 'Exotic', cover > 0 ) |>
  summarise(nspp_exotic = n()) |>
  ungroup() |>
  mutate(trt = str_to_lower(trt),
    year = as.factor(year),
         plot = str_c(block, "_", plot),
         block = as.factor(block)) |>
  left_join(seedlings) |> 
  left_join(saplings) |>
  left_join(ba_20) |>
  tidyr::replace_na(list(seedling_density_ha = 0,
                         sapling_density_ha = 0)) |> print(n=99)

native_richness <- read_csv("data/RDS-2023-0063/Data/understory.csv") |>
  janitor::clean_names() |>
  group_by(block, trt, year, plot) |>
  filter(origin == 'Native', cover > 0 ) |>
  summarise(nspp_native = n()) |>
  ungroup() |>
  mutate(trt = str_to_lower(trt),
         year = as.factor(year),
         plot = str_c(block, "_", plot),
         block = as.factor(block)) |>
  left_join(seedlings) |> 
  left_join(saplings) |>
  left_join(ba_20) |>
  tidyr::replace_na(list(seedling_density_ha = 0,
                         sapling_density_ha = 0)) |> print(n=99)

invaded <- read_csv("data/RDS-2023-0063/Data/understory.csv") |>
  janitor::clean_names() |>
  group_by(block, trt, year, plot) |>
  filter(origin == 'Exotic') |>
  summarise(cover_exotics = sum(cover)) |>
  ungroup() |>
  mutate(invaded = ifelse(cover_exotics >0, 1, 0),
         year = as.factor(year),
         plot = str_c(block, "_", plot),
         block = as.factor(block),
         trt = str_to_lower(trt)) |>
  left_join(seedlings) |> 
  left_join(saplings) |>
  left_join(ba_20) |>
  mutate(burned = ifelse(trt %in% c('bo', 'tb'), 'burned', 'unburned'),
         treated = ifelse(trt == 'co', '0untreated', 'treated')) |>
  tidyr::replace_na(list(seedling_density_ha = 0,
                         sapling_density_ha = 0)) 

native <- read_csv("data/RDS-2023-0063/Data/understory.csv") |>
  janitor::clean_names() |>
  group_by(block, trt, year, plot) |>
  filter(origin == 'Native') |>
  summarise(cover_native = sum(cover)) |>
  ungroup() |>
  mutate(year = as.factor(year),
         plot = str_c(block, "_", plot),
         block = as.factor(block),
         trt = str_to_lower(trt)) |>
  left_join(seedlings) |> 
  left_join(saplings) |>
  left_join(ba_20) |>
  mutate(burned = ifelse(trt %in% c('bo', 'tb'), 'burned', 'unburned'),
         treated = ifelse(trt == 'co', '0untreated', 'treated')) |>
  tidyr::replace_na(list(seedling_density_ha = 0,
                         sapling_density_ha = 0)) |> print(n=99)

catford <-
  left_join(invaded, native) |>
  left_join(native_richness) |>
  left_join(exotic_richness) |>
  tidyr::replace_na(list(nspp_exotic = 0)) |>
  mutate(rel_exotic_richness = nspp_exotic/(nspp_native + nspp_exotic),
         rel_exotic_cover = cover_exotics/(cover_native + cover_exotics),
         trt = case_when(trt == 'co' ~ "Control",
                         trt == 'tb' ~ "TRT: Thin + Burn",
                         trt == "to" ~ 'TRT: Thin',
                         trt == 'bo' ~ 'TRT: Burn'))
# write_csv(catford, 'data/hood_catford.csv')
ggplot(y = trt, x = rel_exotic_richness) +
  geom_boxplot() +
  facet_wrap(~year)
ggplot(y = trt, x = rel_exotic_cover) +
  geom_boxplot() +
  facet_wrap(~year)

# brm ============

mb <- brm(rel_exotic_cover ~ year + ba_m2_ha*burned + (1|block/plot), data = catford, family = 'zero_inflated_beta')
summary(mb)
mb1 <- brm(rel_exotic_cover ~ year*trt + (1|block/plot), data = catford, family = 'zero_inflated_beta')
summary(mb1); plot(conditional_effects(mb1), ask = FALSE)
mb2 <- brm(rel_exotic_richness ~ year*trt + (1|block/plot), data = catford, family = 'zero_inflated_beta')
summary(mb2); plot(conditional_effects(mb2), ask = FALSE)
mb3 <- brm(invaded ~ year*trt + (1|block/plot), data = catford, family = 'bernoulli')
summary(mb3);plot(conditional_effects(mb3), ask = FALSE)

conditional_effects(mb3)[[3]] |> as.data.frame() |> dplyr::rename(value = 3) |> mutate(name = 'invaded') -> di
conditional_effects(mb2)[[3]] |> as.data.frame() |> dplyr::rename(value = 3) |> mutate(name = 'exotic_relative_richness') -> der
conditional_effects(mb1)[[3]] |> as.data.frame() |> dplyr::rename(value = 3) |> mutate(name = 'exotic_relative_cover') -> dec

save(mb1, mb2, mb3, file = "data/hood_bayesian_models.rda")

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

mc <- lme4::glmer(invaded ~ year + ba_m2_ha*burned + (1|block/plot), family = 'binomial', data = invaded)##datawizard::standardise(invaded, exclude = 'invaded'))
summary(mc); performance::check_model(mc); performance::r2(mc)
ggeffects::ggpredict(mc) |> plot() 
ggeffects::ggpredict(mc, terms = 'ba_m2_ha') |> plot() 
ggsave(filename = 'out/hood_ba.png', width = 3.5, height = 3.5, bg='white')

m <- lme4::glmer(invaded ~ year:trt + (1|block/plot), family = 'binomial', data = invaded)##datawizard::standardise(invaded, exclude = 'invaded'))
summary(m); performance::r2(m)
performance::check_model(m)


p1<- ggeffects::ggpredict(m, terms = c('year', 'trt')) |>
  as.data.frame() |>
  mutate(group = case_when(group == 'co' ~ "Control",
                           group == 'tb' ~ "TRT: Thin + Burn",
                           group == "to" ~ 'TRT: Thin',
                           group == 'bo' ~ 'TRT: Burn')) |>
  ggplot(aes(y=group, color = group)) +
  geom_point(aes(x=predicted)) +
  geom_errorbar(aes(y=group, xmin = conf.low, xmax = conf.high), width=.2) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~x) +
  xlab("P(Invasion)") +
  ggthemes::theme_clean() +
  theme(panel.background = element_rect(color = 'black'),
        legend.position = 'none',
        axis.title.y = element_blank()) ;p1

m <- lme4::glmer(nspp ~ year*trt + (1|block/plot), family = 'poisson', data = datawizard::standardise(exotic_richness, exclude = 'nspp'))
summary(m)

p2 <- ggeffects::ggpredict(m, terms = c('year', 'trt')) |>
  as.data.frame() |>
  mutate(group = case_when(group == 'CO' ~ "Control",
                           group == 'TB' ~ "TRT: Thin + Burn",
                           group == "TO" ~ 'TRT: Thin',
                           group == 'BO' ~ 'TRT: Burn')) |>
  ggplot(aes(y=group, color = group)) +
  geom_point(aes(x=predicted)) +
  geom_errorbar(aes(y=group, xmin = conf.low, xmax = conf.high), width=.2) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~x) +
  xlab("Exotic Richness") +
  ggthemes::theme_clean() +
  theme(panel.background = element_rect(color = 'black'),
        legend.position = 'none',
        axis.title.y = element_blank())

ggpubr::ggarrange(p1, p2, ncol = 1, nrow=2) |> ggsave(filename = 'out/hoodetal_preds.png', width =5, height =5, bg="white")

pred_df <- read_csv("data/RDS-2023-0063/Data/understory.csv") |>
  janitor::clean_names() |>
  dplyr::select(block, trt, year, plot) |>
  unique() |>  print(n=99)
  pivot_wider(names_from = year, values_from = pred) |>
  janitor::clean_names()

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
                           yst >10 ~ '3. 16-20 post'))

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