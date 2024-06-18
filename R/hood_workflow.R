# table 1: summary stats
# year, treatment, density, basal area, QMD, SI, Proportion ponderosa pine, 
# Density (seedlings), density (saplings)

source("R/a_tc_data_prep.R")
library(ggtext)

# most of the stuff from table 1

# https://rdrr.io/github/ryanmismith/inventoryfunctions/src/R/ExpansionFactor.R

# EXP.F <- function(DBH, BAF.Area) {
#   
#   TPA <- (BAF.Area/(((DBH/2.54)^2)*0.005454))/0.4046856
#   
#   X <- ifelse(BAF.Area <= 1, 1/BAF.Area, TPA)
#   
#   return(X)
#   
# }
# 
# BAPH <- function(Stand, Plot, BA, EXPF){
#   treebasal <- x <- NULL
#   temp <- tidyr::tibble(Stand, Plot, BA, EXPF)
#   temp <- temp %>%
#     dplyr::mutate(
#       treebasal = BA * EXPF
#     ) %>%
#     dplyr::select(Stand, Plot, treebasal)
#   
#   temp <- temp %>%
#     dplyr::group_by(Stand, Plot) %>%
#     dplyr::mutate(
#       x = sum(treebasal)
#     ) %>%
#     dplyr::select(Stand, Plot, x)
#   
#   temp$x <- round(temp$x, 2)
#   
#   return(temp$x)
# }
cp_tree <- cp_tree |>
  dplyr::mutate(TreeStatus = ifelse(TreeStatus == "L", "live", "dead"))

# ba_ft2ac <- cp_tree |>
#   group_by(PlotCode, TreeStatus, PlotTreatmentStatus, phase_adj) |>
#   summarise(ba = sum(CalcCol_BasalAreaContribution)) |>
#   ungroup()

ba_m2ha <- cp_tree |>
  filter(!is.na(BAPrismSize)) |>
  group_by(PlotCode, BAPrismSize, TreeStatus, PlotTreatmentStatus, phase_adj) |>
  summarise(n=n()) |>
  ungroup() |>
  mutate(ba_ft2peracre = BAPrismSize * n,
         ba_m2pherha = ba_ft2peracre * 0.22956)|>
  group_by(PlotCode, TreeStatus, PlotTreatmentStatus, phase_adj) |>
  summarise(ba_m2pherha = sum(ba_m2pherha),
            ba_ft2peracre = sum(ba_ft2peracre)) |>
  ungroup()

density <- cp_tree |>
  group_by(PlotCode, PlotTreatmentStatus, phase_adj, Species)|>
  summarise(tpa = sum(CalcCol_DensityContribution)) |>
  ungroup() |>
  mutate(site = str_sub(PlotCode, 1,1)) 


density |>
  group_by(PlotTreatmentStatus, phase_adj, Species, site) |>
  summarise(tpa = mean(tpa)) |>
  ungroup() |>
  dplyr::mutate(site = ifelse(site ==  "E", "Estes Valley", "Phantom")) |>
  ggplot(aes(x=phase_adj, y=tpa, fill = Species)) +
  geom_bar(stat = "identity", color = "black") +
  facet_grid(site ~ PlotTreatmentStatus) +
  ylab("Average Trees Per Acre per Plot") +
  theme(axis.title.x = element_blank())

summary_stats <- cp_tree |>
  group_by(new_visit_code, phase_adj, PlotTreatmentStatus) |>
  summarise(tree_density = sum(CalcCol_DensityContribution, na.rm = T),
            basal_area = (BAPrismSize/((DBH/2)^2 * pi)) * 100,
            n=n(),
            quadratic_mean_diameter = sqrt(sum(DBH * DBH, na.rm=T)/n)) |>
  ungroup() |>
  group_by(site = str_sub(new_visit_code,1,1),phase_adj,  PlotTreatmentStatus) |>
  summarise(tree_density = mean(tree_density),
            basal_area = mean(basal_area, na.rm=T),
            QMD = mean(quadratic_mean_diameter)) |>
  print(n=99)

# sapling density: NA's in the plot sizes
# estes valley has all 0.05 plot sizes
sapling_density <- cp_sapling |>
  mutate(FixedPlotSize = ifelse(is.na(FixedPlotSize) & str_sub(new_visit_code, 1, 2) =="ES", 0.005, FixedPlotSize)) |>
  mutate(FixedPlotSize = ifelse(is.na(FixedPlotSize) & str_sub(new_visit_code, 1, 2) =="PH", 0.01, FixedPlotSize)) |>
  mutate(FixedPlotSize = ifelse(is.na(FixedPlotSize) & str_sub(new_visit_code, 1, 2) %in% c("P1", "P2") & phase_adj == "post04-5", 0.0066, FixedPlotSize)) |>
  mutate(FixedPlotSize = ifelse(is.na(FixedPlotSize) & str_sub(new_visit_code, 1, 2) %in% c("P1", "P2") & phase_adj == "post10-11", 0.05, FixedPlotSize)) |>
  group_by(new_visit_code, phase_adj, PlotTreatmentStatus, FixedPlotSize) |>
  summarise(n_saplings = n()) |>
  ungroup() |> 
  mutate(density = n_saplings/FixedPlotSize)|>
  group_by(site = str_sub(new_visit_code,1,1), phase_adj, PlotTreatmentStatus) |>
  summarise(sapling_density = mean(density), sd = sd(density), n=n()) |>
  ungroup()

seedling_density <-
  sp_seedling |>
  filter(Species != "NONE") |>
  mutate(PlotSize_m = as.numeric(PlotSize_m)) |>
  group_by(phase_adj, PlotTreatmentStatus, new_visit_code, PlotSize_m) |>
  summarise(count = sum(TotalCount)) |>
  ungroup() |>
  mutate(density = count/PlotSize_m) |>
  group_by(site = str_sub(new_visit_code,1,1), phase_adj,  PlotTreatmentStatus) |>
  summarise(seedling_density = mean(density)) |>
  ungroup(); seedling_density |> print(n=100)


proportion_pondo <-
  cp_tree |>
  group_by(phase_adj, PlotTreatmentStatus, new_visit_code)|>
  mutate(total_trees = n(),
         PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus)) |>
  ungroup() |>
  group_by(phase_adj, PlotTreatmentStatus, new_visit_code, Species, total_trees) |>
  summarise(n=n()) |>
  ungroup() |>
  mutate(proportion = n/total_trees) |>
  filter(Species == "PIPO") |>
  group_by(site = str_sub(new_visit_code,1,1),phase_adj,  PlotTreatmentStatus) |>
  summarise(proportion_ponderosa = mean(proportion), sd = sd(proportion, na.rm=T), n=n())

ggplot(proportion_pondo, aes(x=phase_adj, y=proportion_ponderosa, fill = PlotTreatmentStatus)) +
  geom_bar(stat = "identity", position = "dodge") 

# figure 1 tree mortality


ba_m2ha |>
  group_by(TreeStatus, PlotTreatmentStatus, phase_adj, site = str_sub(PlotCode, 1,1)) |>
  summarise(ba = mean(ba_m2pherha)) |>
  ungroup() |>
  ggplot(aes(x = phase_adj, y=ba, fill = TreeStatus)) +
  geom_bar(stat = "identity") +
  facet_grid(site~PlotTreatmentStatus) +
  xlab("Phase") +
  ylab("Basal Area m<sup>2</sup> ha<sup>-1</sup>") +
  theme_bw() +
  theme(axis.title.y = element_markdown())
ggsave("out/hood_replication/figure1.png", height = 5, width=7, bg="white" )

ba_ft2ac |>
  group_by(TreeStatus, PlotTreatmentStatus, phase_adj, site = str_sub(PlotCode, 1,1)) |>
  summarise(ba = mean(ba)) |>
  ungroup() |>
  ggplot(aes(x = phase_adj, y=ba, fill = TreeStatus)) +
  geom_bar(stat = "identity") +
  facet_grid(site~PlotTreatmentStatus) +
  xlab("Phase") +
  ylab("Basal Area ft<sup>2</sup> ac<sup>-1</sup>") +
  theme_bw() +
  theme(axis.title.y = element_markdown())
ggsave("out/hood_replication/figure1.png", height = 5, width=7, bg="white" )

ba_m2ha |> 
  pivot_wider(names_from = c(phase_adj), values_from = ba_m2pherha, values_fill = 0) |>
  janitor::clean_names() |>
  mutate(d01 = post0_1 - x01_pre, 
         d05 = post04_5 - x01_pre, 
         d10 = post10_11 - x01_pre) |>
  dplyr::select(-x01_pre, -starts_with("post")) |>
  pivot_longer(-c(plot_code,tree_status, plot_treatment_status)) |>
  mutate(site = str_sub(plot_code, 1,1)) |>
  ggplot(aes(x = name, y = value, fill = tree_status)) +
  geom_hline(yintercept = 0, lty = 3) +
  geom_boxplot() +
  facet_grid(site~plot_treatment_status) +
  ylab("Change from Pre-Treatment m<sup>2</sup> ha<sup>-1</sup>") +
  xlab("Phase") +
  theme_bw() +
  theme(axis.title.y = element_markdown())

ggsave("out/hood_replication/figure1_reimagined.png", height = 7, width=10, bg="white" )

# figure 2 dbh changes

cp_tree |>

cp_tree |>
  filter(!is.na(DBH), TreeStatus == "L") |>
  mutate(site = str_sub(new_visit_code,1,1),
         PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus),
         tpa = EXP.F(BAF.Area = BAPrismSize, DBH=DBH),
         DBH_class = case_when(
           DBH < 7.5 ~ '07.5',
           DBH >= 7.5 & DBH < 12.5 ~ '12.5',
           DBH >= 12.5 & DBH < 17.5 ~ '17.5',
           DBH >= 17.5 & DBH < 22.5 ~ '22.5',
           DBH >= 22.5 & DBH < 27.5 ~ '27.5')) |> 
  ggplot(aes(x=DBH_class, fill = Species, y = tpa)) +
  geom_bar(stat = 'identity', position = "stack")  +
  facet_grid(site + PlotTreatmentStatus~phase_adj) +
  theme_bw() +
  scale_fill_brewer(palette = "Spectral") +
  ylab("Trees ha<sup>-1</sup>") +
  theme(axis.title.y = element_markdown())

ggsave("out/hood_replication/figure2.png", height = 8, width=11, bg="white" )

# figure 3 seedling/sapling: difference from pre-treatment

sp_seedling |>
  filter(Species != "NONE") |>
  mutate(PlotSize_m = as.numeric(PlotSize_m),
         PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus)) |>
  group_by(phase_adj, PlotTreatmentStatus, new_visit_code, PlotSize_m, Species) |>
  summarise(count = sum(TotalCount)) |>
  ungroup() |>
  mutate(density = count/PlotSize_m) |>
  group_by(site = str_sub(new_visit_code,1,1), phase_adj,  PlotTreatmentStatus, Species) |>
  summarise(seedling_density = mean(density)) |>
  ungroup() |>
  ggplot(aes(x=PlotTreatmentStatus, y=seedling_density, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Spectral") +
  facet_grid(site ~phase_adj)  +
  theme_bw()

ggsave("out/hood_replication/figure3a.png", height = 6, width=8, bg="white" )

cp_sapling |>
  filter(!Species %in% c("NONE"), !is.na(Species)) |>
  mutate(FixedPlotSize = ifelse(is.na(FixedPlotSize) & str_sub(new_visit_code, 1, 2) =="ES", 0.005, FixedPlotSize)) |>
  mutate(FixedPlotSize = ifelse(is.na(FixedPlotSize) & str_sub(new_visit_code, 1, 2) =="PH", 0.01, FixedPlotSize)) |>
  mutate(FixedPlotSize = ifelse(is.na(FixedPlotSize) & str_sub(new_visit_code, 1, 2) %in% c("P1", "P2") & phase_adj == "post04-5", 0.0066, FixedPlotSize)) |>
  mutate(FixedPlotSize = ifelse(is.na(FixedPlotSize) & str_sub(new_visit_code, 1, 2) %in% c("P1", "P2") & phase_adj == "post10-11", 0.05, FixedPlotSize)) |>
  group_by(new_visit_code, phase_adj, PlotTreatmentStatus, FixedPlotSize, Species) |>
  summarise(n_saplings = n()) |>
  ungroup() |> 
  mutate(density = n_saplings/FixedPlotSize)|>
  group_by(site = str_sub(new_visit_code,1,1), phase_adj, PlotTreatmentStatus, Species) |>
  summarise(sapling_density = mean(density), sd = sd(density), n=n()) |>
  ungroup() |>
  ggplot(aes(x=PlotTreatmentStatus, y=sapling_density, fill=Species)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Spectral") +
  facet_grid(site ~phase_adj) +
  theme_bw()

ggsave("out/hood_replication/figure3b.png", height = 6, width=8, bg="white" )

  