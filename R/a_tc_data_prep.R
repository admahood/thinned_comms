library(tidyverse)
library(janitor)

# site and plot visits =========================================================
# pha 2011 trt, esv 2012
# not super necessary
site_visits <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                            sheet = "Export_TBL_Site_Visit") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_Site_Visit"))

lut_trt_unit <- c('PHA-1-2' = "PHA-1-2", 'PHA-1-610' = "PHA-1-2", 
                  'PHA-1-615' = "PHA-1-2", 'PHA-1-616'= "PHA-1-2", 
                  'PHA-1-3' = "PHA-1-3",
                  'PHA-1-589' = "PHA-1-3", 'PHA-1-596' = "PHA-1-3", 
                  'PHA-1-601' = "PHA-1-3", 'PHA-1-602' = "PHA-1-3",
                  'PHA-1-603' = "PHA-1-3", 'PHA-1-604' = "PHA-1-3", 
                  'PHA-1-606' = "PHA-1-3", 'PHA-2C1' = 'PHA-2-1', 
                  'PHA-2C2' = 'PHA-2-2', 'PHA-2C3' ='PHA-2-3', 
                  'PHA-2T1' = 'PHA-2-1', 'PHA-2T2' = 'PHA-2-2', 
                  'PHA-2T3' = 'PHA-2-3', 'ESV-13'= 'ESV-13', 
                  'ESV-28'='ESV-28', 'ESV-34' = 'ESV-34')
# this is the main thing
plot_visits <- readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                                           sheet = "Export_TBL_PlotVisit")|>
  bind_rows(readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                 sheet = "Export_TBL_PlotVisit")) |>
  dplyr::mutate(trt_u_adj = lut_trt_unit[TreatmentUnit]) 
 
treatments <- plot_visits |>
  dplyr::select(TreatmentDate1, TreatmentUnit) |>
  unique() |>
  na.omit() |>
  bind_rows(tibble(TreatmentDate1 = rep(2017, 3), TreatmentUnit = c("PHA-2C1", "PHA-2C2", "PHA-2C3"))) |>
  bind_rows(tibble(TreatmentDate1 = rep(2012, 2), TreatmentUnit = c("PHA-1-2", "PHA-1-3")))

lut_trtd <- pull(treatments, TreatmentDate1)
names(lut_trtd) <- pull(treatments, TreatmentUnit)


plot_visits <- plot_visits |>
  mutate(trt_year = lut_trtd[TreatmentUnit],
         visit_year = lubridate::year(VisitDate),
         yst = visit_year - trt_year,
         phase_adj = case_when(yst < 0 ~ "01_Pre",
                               yst >= 0 & yst < 4 ~ "post0-1",
                               yst %in% c(4,5) ~ "post04-5",
                               yst %in% c(10,11) ~ "post10-11"
                               ))  |>
  replace_na(list(phase_adj = "not_measured")) |>
  mutate(new_visit_code = str_c(PlotCode, ".", phase_adj)) 

plot_visits |>
  group_by(trt_year) |>
  summarise(n_plots = length(unique(PlotCode))) |>
  print(n=21)|>
  filter(trt_year != 2017)
# 
plot_visits_10y <- plot_visits|>
  filter(trt_year != 2017);glimpse(plot_visits)

# filter(plot_visits, is.na(yst)) |> glimpse() # 4 plots not sampled

# pull(plot_visits, yst) |> unique()

lut_phase <- pull(plot_visits, VisitCode)
names(lut_phase) <- pull(plot_visits, phase_adj)

new_vcs<- dplyr::select(plot_visits, VisitCode, new_visit_code, phase_adj, trt_year) |> unique()
# table(plot_visits$phase_adj)

# plot_visits |>
#   group_by(phase_adj, trt_year) |>
#   summarise(n_plots = length(unique(PlotCode)))


# khr fuels ====================================================================
khr_fuels <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                 sheet = "Export_TBL_1KhrFuel") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_1KhrFuel")) |>
  left_join(new_vcs) |>
  dplyr::select(-ID, -DataFlag)|>
  filter(trt_year != 2017); glimpse(khr_fuels)
# 
# khr_fuels |>
#   group_by(phase_adj, trt_year, PlotSize) |>
#   summarise(n_plots = length(unique(PlotCode)))
# 
# khr_fuels |>
#   mutate(
#     density_constant = ifelse(RottenSound == "Sound", 0.65, 0.45),
#     mass = (Diameter1/2)*(Diameter2/2)*pi*Length*density_constant) |>
#   group_by(new_visit_code, phase_adj, trt_year,PlotTreatmentStatus,TreatmentUnit) |>
#   reframe(mass = sum(mass, na.rm=T)) |>
#   ggplot(aes(x=TreatmentUnit |> str_sub(1,3), y=mass, fill=PlotTreatmentStatus)) +
#   geom_boxplot() +
#   ylab("Mass (r1 * r2 * pi * length * density)") +
#   xlab("Treatment Units (Phase = 10-11 years post-treatment)")
# ggsave("out/khr_fuels.png", width=15, height=8)

# centerPlotTreeData ============

cp_tree <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                               sheet = "Export_TBL_CenterPlotTreeData") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_CenterPlotTreeData")) |>
  left_join(new_vcs)|>
  filter(trt_year != 2017)|>
  dplyr::select(-DataFlag, -TreeOrder, -PlotSize, -TagID, -DeadTop, -ForkedTrunk, 
                -WitchesBroom, -MistletoeShoots, -FireScar, -LightningScar, -Conks, 
                -HollowBole, -BrokenTop, -NumCavities, -Notes, -BarkBeetle, -MistletoeRating) |>
  mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus))#; glimpse(cp_tree)
# cp_tree |>
#   group_by(phase_adj, trt_year) |>
#   summarise(n_plots = length(unique(PlotCode)))
# dim(cp_tree)



# 
# ggpubr::ggarrange(
# cp_tree |>
#   filter(Species %in% c("PIPO", "PSME")) |>
#   mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus)) |>
#   ggplot(aes(x=phase_adj, y=DBH, fill=PlotTreatmentStatus)) +
#   geom_boxplot() +
#   facet_wrap(trt_year~Species) +
#   ggtitle("CenterPlotTreeData, DBH", "no TagIDs, so we cant track individual trees")
# ,
# cp_tree |>
#   filter(Species %in% c("PIPO", "PSME")) |>
#   mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus)) |>
#   ggplot() +
#   geom_boxplot(aes(x=phase_adj, y=Height, fill=PlotTreatmentStatus)) +
#   facet_wrap(trt_year~Species) +
#   ggtitle("CenterPlotTreeData, Height", "no TagIDs, so we cant track individual trees")
# , common.legend = TRUE, nrow=1)
# ggsave("out/centerPlotTreeData.png", width=15, height=8)

# center plot sapling ==========================================================
cp_sapling <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                             sheet = "Export_TBL_CenterPlotSapling") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_CenterPlotSapling")) |>
  left_join(new_vcs)|>
  filter(trt_year != 2017); glimpse(cp_sapling)

saplings_wide <- cp_sapling |>
  group_by(new_visit_code, trt_year, phase_adj, PlotTreatmentStatus) |>
  summarise(sapling_density_per_acre = mean(CalcCol_DensityContribution, na.rm=T),
            sapling_ba_ft_per_acre = sum(CalcCol_BasalAreaContribution, na.rm=T)) |>
  ungroup() |> 
  dplyr::select(new_visit_code, sapling_density_per_acre, sapling_ba_ft_per_acre)


# cp_sapling |>
#   group_by(phase_adj, trt_year) |>
#   summarise(n_plots = length(unique(PlotCode)))
# dim(cp_sapling)
# 
# # 
# cp_sapling |>
#   group_by(new_visit_code, trt_year, phase_adj, PlotTreatmentStatus) |>
#   summarise(saplings_per_acre = sum(CalcCol_DensityContribution),
#             basal_area_ft_acre = sum(CalcCol_BasalAreaContribution, na.rm=T)) |>
#   ungroup() |>
#   pivot_longer(cols = c(basal_area_ft_acre, saplings_per_acre)) |>
#   ggplot(aes(x=phase_adj, y=value, fill = PlotTreatmentStatus)) +
#   geom_boxplot() +
#   scale_y_sqrt(breaks = c(1, 10, 50, 100, 500, 1000, 5000, 8000)) +
#   facet_grid(name~trt_year, scales = "free_y") +
#   ggtitle("Saplings")
# ggsave("out/centerPlotSaplings.png", width=10, height=6)

# seedlings ====================================================================
sp_seedling <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                sheet = "Export_TBL_SubPlotSeedling") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_SubPlotSeedling")) |> 
  left_join(new_vcs)|>
  filter(trt_year != 2017); glimpse(sp_seedling)

seedling_density <- sp_seedling |>
  mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated",
                                      "Control", PlotTreatmentStatus)) |>
  group_by(new_visit_code, trt_year, phase_adj, PlotTreatmentStatus, PlotSize_ac, MeterID) |>
  summarise(total_count = sum(TotalCount)) |>
  ungroup() |>
  group_by(new_visit_code, trt_year, phase_adj, PlotTreatmentStatus, PlotSize_ac) |> # then mean by transect
  summarise(mean_count = mean(total_count)) |>
  ungroup() |>
  mutate(seedlings_per_acre = mean_count/PlotSize_ac) |>
  dplyr::select(new_visit_code, seedlings_per_acre)
 
# sp_seedling |>
#   group_by(phase_adj, trt_year) |>
#   summarise(n_plots = length(unique(PlotCode)))
# 
# sp_seedling |>
#   mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated",
#                                       "Control", PlotTreatmentStatus)) |>
#   group_by(new_visit_code, trt_year, phase_adj, PlotTreatmentStatus, PlotSize_ac, MeterID) |>
#   summarise(total_count = sum(TotalCount)) |>
#   ungroup() |>
#   group_by(new_visit_code, trt_year, phase_adj, PlotTreatmentStatus, PlotSize_ac) |> # then mean by transect
#   summarise(mean_count = mean(total_count)) |>
#   ungroup() |>
#   mutate(seedlings_per_acre = mean_count/PlotSize_ac)  |>
#   ggplot(aes(y=seedlings_per_acre,x = phase_adj, fill=PlotTreatmentStatus)) +
#   geom_boxplot() +
#   scale_y_sqrt(breaks = c(100, 1000, 5000, 10000, 50000)) +
#   facet_wrap(~trt_year, scales = 'free_y') +
#   ggtitle("Plot seedling counts (all species summed)")
# ggsave("out/SubPlotSeedlings.png", width=10, height=6)

# ground fuel ==================================================================
sp_groundfuel <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                 sheet = "Export_TBL_SubPlotGroundFuel") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_SubPlotGroundFuel")) |>
  left_join(new_vcs)|>
  filter(trt_year != 2017); glimpse(sp_groundfuel)

# sp_groundfuel |>
#   group_by(phase_adj, trt_year) |>
#   summarise(n_plots = length(unique(PlotCode)))
# 
# sp_groundfuel |>
#   group_by(new_visit_code, Substrate,phase_adj, trt_year,PlotTreatmentStatus) |>
#   summarise(depth_nounit = mean(Observation)) |>
#   ungroup() |>
#   ggplot(aes(x=phase_adj, y=depth_nounit, fill = PlotTreatmentStatus)) +
#   geom_boxplot() +
#   facet_grid(Substrate ~ trt_year) +
#   scale_y_continuous(limits = c(0,2.5)) + # one value of 10 that impedes plot interpretibility
#   ggtitle("Ground Fuel")
# ggsave("out/SubPlotGroundFuel.png", width=10, height=6)

# woody fuel ==================================================================
sp_woodyfuel <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                   sheet = "Export_TBL_SubPlotWoodyFuel") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_SubPlotWoodyFuel")) |>
  left_join(new_vcs)|>
  filter(trt_year != 2017); glimpse(sp_woodyfuel)

# sp_woodyfuel |>
#   group_by(phase_adj, trt_year) |>
#   summarise(n_plots = length(unique(PlotCode)))
# 
# sp_woodyfuel |>
#   group_by(new_visit_code, Measurement, phase_adj, trt_year,PlotTreatmentStatus) |>
#   summarise(fuelload_nounit = mean(Calibrated)) |>
#   ungroup() |>
#   ggplot(aes(x=phase_adj, y=fuelload_nounit, fill = PlotTreatmentStatus)) +
#   geom_boxplot() +
#   facet_wrap(~Measurement) +
#   ggtitle("Woody Fuel")
# ggsave("out/SubPlotWoodyFuel.png", width=10, height=6)

# tree canopy cover ==================================================================
tree_canopy_cover <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                  sheet = "Export_TBL_TreeCanopyCover") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_TreeCanopyCover")) |>
  left_join(new_vcs)|>
  filter(trt_year != 2017); glimpse(tree_canopy_cover)

# tree_canopy_cover |>
#   group_by(phase_adj, trt_year) |>
#   summarise(n_plots = length(unique(PlotCode)))

# tree_canopy_cover |>
#   ggplot(aes(x=phase_adj, y=PercentCover, fill = PlotTreatmentStatus)) +
#   geom_boxplot() +
#   ggtitle("Tree Canopy Cover") +
#   facet_wrap(~trt_year)
# ggsave("out/tree_canopy_cover.png", width=10, height=6)

# tree group transect ==================================================================
tree_group_transect <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                       sheet = "Export_TBL_TreeGroupTransect") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_TreeGroupTransect")) |>
  left_join(new_vcs)|>
  filter(trt_year != 2017); glimpse(tree_group_transect)

# tree_group_transect |>
#   group_by(phase_adj, trt_year) |>
#   summarise(n_plots = length(unique(PlotCode)))

# tree_group_transect |>
#   group_by(new_visit_code, TreeGroupClass, phase_adj, trt_year,PlotTreatmentStatus) |>
#   summarise(PercentCover = sum(((EndFt - StartFt)/77)*100)) |>
#   ungroup() |>
#   ggplot(aes(x=phase_adj, y=PercentCover, fill = PlotTreatmentStatus)) +
#   geom_boxplot() +
#   ggtitle("Tree Group Transects") +
#   facet_wrap(~trt_year)
# ggsave("out/tree_group_transects.png", width=10, height=6)

# understory heights ===========================================================
understory_heights <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                         sheet = "Export_TBL_MS_UnderstoryHeights") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_MS_UnderstoryHeights")) |>
  left_join(new_vcs)|>
  filter(trt_year != 2017); glimpse(understory_heights)

understory_heights |>
  group_by(phase_adj, trt_year) |>
  reframe(n_plots = length(unique(PlotCode)))

# understory_heights |>
#   group_by(new_visit_code, PlantType, phase_adj, trt_year,PlotTreatmentStatus) |>
#   summarise(avg_height_in = sum(Height_inches)/4) |>
#   ungroup() |>
#   ggplot(aes(x=phase_adj, y=avg_height_in, fill = PlotTreatmentStatus)) +
#   geom_boxplot() +
#   ggtitle("Understory Heights") +
#   facet_wrap(~trt_year)
# ggsave("out/understory_heights.png", width=10, height=6)

# plant understory =============================================================
esv <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx",
                              sheet = "Export_TBL_MS_Spokes")
pha <- readxl::read_xlsx("data/PHA_Export_20231024.xlsx",
                                  sheet = "Export_TBL_MS_Spokes") 

comm_raw <- bind_rows(esv, pha)|>
  filter(TreatmentUnit %in% plot_visits$TreatmentUnit) |>
  left_join(new_vcs) |>
  mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", 
                                      "Control", PlotTreatmentStatus))|>
  filter(trt_year != 2017)

comm_raw |>
  group_by(phase_adj) |>
  summarise(n=length(unique(PlotCode)))

# glimpse(comm_raw)

# comm_raw |>
#   group_by(phase_adj, trt_year) |>
#   reframe(n_plots = length(unique(PlotCode)))

non_plant_codes <- c("FWD", "Fuel in Air", "Substrate")

meta<- comm_raw %>%
  dplyr::select(VisitCode, PlotTreatmentStatus, TreatmentUnit) %>%
  mutate(PlotTreatmentStatus=ifelse(PlotTreatmentStatus == "NotTreated",
                                    "Control", PlotTreatmentStatus)) |>
  left_join(new_vcs) |>
  tidyr::separate(new_visit_code, into = c("PlotCode", "Phase"), sep = "\\.",remove = F)

cover_plot <- comm_raw %>%
  mutate(SpeciesCode = str_to_lower(SpeciesCode)) %>%
  dplyr::select(-Notes, -DataFlag, -ID, -Cover_class) |>
  tidyr::replace_na(list('Cover_Bottom' = 0, "Cover_Top" = 0)) |>
  unique() %>% # 22113 - 21940 = 173 duplicated rows
  mutate(Cover_Top = ifelse(Transect == "P", .1, Cover_Top), # option 1 convert P plants to trace cover and add to another transect
         Transect = str_replace_all(Transect, "P", "N"),
         cover = Cover_Top + Cover_Bottom) |> 
  dplyr::filter(Cover_Bottom + Cover_Top > 0)  %>%
  dplyr::select(-Cover_Top, -Cover_Bottom) |>
  group_by(trt_year, new_visit_code, Transect, CodeType, SpeciesCode, OriginalSpeciesCode, PlotTreatmentStatus) |>
  summarise(cover = mean(cover)) |>
  ungroup() |>
  # filter(Transect != "P") %>% # option 2: remove trace plants
  group_by(new_visit_code) %>%
  mutate(n_transects = length(unique(Transect))) %>%
  ungroup() %>% #filter(n_transects != 4) |> print(n=999)
  # group_by(new_visit_code, CodeType, SpeciesCode, Transect, n_transects) %>%
  # summarise(cover = (Cover_Top + Cover_Bottom)) %>%
  # ungroup() %>%
  group_by(new_visit_code, CodeType, SpeciesCode, PlotTreatmentStatus) %>%
  summarise(cover = sum(cover)/n_transects) %>%
  ungroup() |>
  unique()#;cover_plot

cover_plot |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('PlotCode', 'phase_adj')) |>
  group_by(phase_adj) |>
  summarise(n=length(unique(PlotCode)))

comm <- cover_plot %>%
  filter(!CodeType %in% non_plant_codes) %>%
  dplyr::select(-CodeType, -PlotTreatmentStatus) %>%
  pivot_wider(names_from = SpeciesCode,
              values_from = cover,# values_fn = function(x)sum(x, na.rm=T),
              values_fill = 0) %>%
   dplyr::mutate_all(list(~ replace(., is.na(.), 0))) %>%
  tibble::column_to_rownames("new_visit_code") %>%
  filter(rowSums(.)>0)
# 
# rowSums(comm)
# dim(comm)
comm <- comm[,colSums(comm)>0];dim(comm) # 20 species have 0 occurrences


plant_cover_by_type <-
  cover_plot %>%
  filter(!CodeType %in% non_plant_codes) |>
  group_by(new_visit_code, CodeType) |>
  summarise(cover = sum(cover, na.rm=T)) |>
  ungroup()

fuel_cover_by_type <- 
  cover_plot %>%
  filter(CodeType %in% non_plant_codes) |>
  group_by(new_visit_code, CodeType, SpeciesCode) |>
  summarise(cover = sum(cover, na.rm=T)) |>
  ungroup()

herbaceous_cover <-
  plant_cover_by_type %>%
  filter(!CodeType %in% non_plant_codes) |>
  mutate(code_modified = ifelse(CodeType %in% c("Forb", "Graminoid", "UNK"), "Herbaceous", CodeType))  |>
  group_by(new_visit_code, code_modified) |>
  summarise(cover = sum(cover, na.rm=T)) |>
  ungroup()

# ggplot(fuel_cover_by_type, aes(x=cover, fill = SpeciesCode)) +
#   geom_histogram() +
#   facet_wrap(~CodeType, scales = "free") +
#   ggtitle("aggregated cover type summaries")
# ggsave("out/fuel_cover.png")
# 
# ggplot(plant_cover_by_type, aes(x=cover, fill = CodeType)) +
#   geom_histogram() +
#   facet_wrap(~CodeType, scales = "free") +
#   ggtitle("aggregated cover type summaries")
# ggsave("out/plant_cover.png")

# investigating repeated rows with different values ============================
# sites_to_investigate <- fuel_cover_by_type |> filter(cover>120) |> pull(new_visit_code)
# 
# comm_raw |> filter(new_visit_code %in% sites_to_investigate, CodeType == "Substrate") |>
#   dplyr::select(-Cover_class, -Notes, -DataFlag, -ID, -TreatmentUnit, -SampleUnit) |>
#   arrange(new_visit_code, Transect, SpeciesCode) |> print(n=999)
# 
# comm_raw |>
#   group_by(new_visit_code, OriginalSpeciesCode, Transect, SpeciesCode) |>
#   summarise(n_dups = n(), 
#             covertops = paste(Cover_Top, collapse = ", "),
#             coverbottoms = paste(Cover_Bottom, collapse = ", ") ) |>
#   filter(n_dups>1, Transect != "P") |>
#   arrange(OriginalSpeciesCode) |>
#   print(n=999) |>
#   write_csv("data/quasiduplicate_investigation.csv")
# 
# unique(cover_by_type$CodeType) 


# 
# write_csv(herbaceous_cover, "data/herb_cover_by_plot.csv")
# write_csv(bind_rows(plant_cover_by_type, fuel_cover_by_type |> 
#                       dplyr::select(-CodeType, CodeType = SpeciesCode)),
#           "data/cover_by_plot.csv")
# 
# ggplot(herbaceous_cover, aes(x=cover, fill = code_modified)) +
#   geom_histogram() +
#   facet_wrap(~code_modified, scales = "free") +
#   ggtitle("aggregated cover type summaries")
# 
# ggsave("out/covertype_aggd.png", width = 7, height=5, bg="white")
# 
# ggplot(cover_by_type, aes(x=cover, fill = CodeType)) +
#   geom_histogram() +
#   facet_wrap(~CodeType, scales = "free") +
#   ggtitle("raw cover type summaries")
# ggsave("out/covertype_raw.png", width = 7, height=5, bg="white")

# big table of data availability ===============================================
bigtab<- bind_rows(
  comm_raw |>
    group_by(phase_adj, trt_year) |>
    reframe(n_plots = length(unique(PlotCode))) |>
    mutate(var = "plant_community"),
  understory_heights |>
    group_by(phase_adj, trt_year) |>
    reframe(n_plots = length(unique(PlotCode))) |>
    mutate(var = "understory_heights"),
  tree_group_transect |>
    group_by(phase_adj, trt_year) |>
    summarise(n_plots = length(unique(PlotCode))) |>
    mutate(var = "tree_group_transects"),
  tree_canopy_cover |>
    group_by(phase_adj, trt_year) |>
    summarise(n_plots = length(unique(PlotCode))) |>
    mutate(var = "tree_canopy_cover"),
  sp_woodyfuel |>
    group_by(phase_adj, trt_year) |>
    summarise(n_plots = length(unique(PlotCode))) |>
    mutate(var = "woody_fuel"),
  sp_groundfuel |>
    group_by(phase_adj, trt_year) |>
    summarise(n_plots = length(unique(PlotCode))) |>
    mutate(var = "ground_fuel"),
  sp_seedling |>
    group_by(phase_adj, trt_year) |>
    summarise(n_plots = length(unique(PlotCode))) |>
    mutate(var = "seedlings"),
  cp_sapling |>
    group_by(phase_adj, trt_year) |>
    summarise(n_plots = length(unique(PlotCode))) |>
    mutate(var = "saplings"),
  cp_tree |>
    group_by(phase_adj, trt_year) |>
    summarise(n_plots = length(unique(PlotCode))) |>
    mutate(var = "center_plot_Tree"),
  khr_fuels |>
    group_by(phase_adj, trt_year) |>
    summarise(n_plots = length(unique(PlotCode)))|>
    mutate(var = "khr_fuels")
) 

bigtab |>
  pivot_wider(names_from = phase_adj, values_from = n_plots, values_fill = 0)
ggplot(bigtab |> mutate(trt_year = paste("Treatment Year:", trt_year)), aes(y=paste(var), x=phase_adj, fill=n_plots)) +
  geom_tile(color="black") +
  scale_fill_viridis_c(direction = -1) +
  facet_wrap(~trt_year)

ggsave("out/data_coverage.png", width=10, height=8)


# creating a species list ======================================================


comm_raw |>
  filter(!CodeType %in% c("Substrate", 
                          "FWD", 
                          # "UNK", 
                          "Fuel in Air")) |>
  pull(SpeciesCode) |>
  unique()

boulder <- readxl::read_xlsx("data/species_attributes_boulder.xlsx") |>
  dplyr::rename(SpeciesCode = FinalCode) |>
  filter(!is.na(FinalName)) |>
  unique()

spp_list <- comm_raw |>
  filter(!CodeType %in% c("Substrate", 
                          "FWD", 
                          "UNK", 
                          "Fuel in Air")) |>
  group_by(SpeciesCode, CodeType) |>
  summarise(prevalence = n()) |>
  ungroup() |> 
  arrange(desc(prevalence)) |>
  left_join(boulder) |> 
  mutate_if(is.character, as.factor) 

summary(spp_list)
write_csv(spp_list,"data/species_list_10.csv")


# creating a species list, with rare species grouped

sp_list <- readr::read_csv("data/species_list_alm_modified.csv") |>
  tidyr::separate(FinalName, into = c("genus", "species"), remove = F, sep = " ") |>
  dplyr::mutate(group = ifelse(prevalence > 10, FinalName, 
                               str_c(NativityL48,"_", Lifespan,"_", CodeType) |> str_remove_all("\\-"))) |>
  dplyr::select(-notes) |>
  unique()
print(sp_list, n=300)

tenyp <- plot_visits_10y |> pull(PlotCode) |> unique()


plotwise_invasion <- cover_plot |>
  left_join(sp_list |> mutate(SpeciesCode = str_to_lower(SpeciesCode))) |>
  filter(!is.na(prevalence),
         NativityL48 == "exotic") |>
  group_by(new_visit_code, NativityL48) |>
  summarise(cover = sum(cover)) |>
  ungroup() |>
  tidyr::separate(new_visit_code, into = c('plot', 'phase'), sep = "\\.") %>%
  pivot_wider(names_from = phase, values_from = cover, values_fill = 0) |>
  mutate(`post10-11` = ifelse(!plot %in% tenyp, NA, `post10-11`)) |>
  left_join(dplyr::select(plot_visits, plot = PlotCode, trt = PlotTreatmentStatus) |> unique())

plotwise_invasion |>
  arrange(`post10-11`) |>
  print(n=99)

plotwise_invasion |>
  mutate_if(is.numeric, function(x) ifelse(x >0, 1,0)) %>%
  pivot_longer(cols = names(.)[c(3:6)]) |>
  group_by(name, trt) |>
  summarise(fraction_invaded = mean(value, na.rm=T)) |>
  ungroup() |>
  ggplot(aes(x=name, y=fraction_invaded, group = trt, color = trt)) +
  geom_line() +
  ggthemes::theme_clean()

plotwise_invasion  %>%
  pivot_longer(cols = names(.)[c(3:6)]) |>
  group_by(name, trt) |>
  summarise(cover_invaded = mean(value, na.rm=T)) |>
  ungroup() |>
  ggplot(aes(x=name, y=cover_invaded, group = trt, color = trt)) +
  geom_line() +
  ggthemes::theme_clean()

ggplot(data.frame(x = c(-1,1,5,10), y = c(49/71, 50/71, 15/71, 4/43))) +
  geom_line(aes(x=x,y=y)) +
  xlab("Years Since Treatment") +
  ylab("Fraction Uninvaded") +
  ggthemes::theme_clean()

# ggsave("out/fraction_uninvaded.png")

# putting together ancillary data for glmm analysis: seedlings, saplings, cp_tree

cp_tree


