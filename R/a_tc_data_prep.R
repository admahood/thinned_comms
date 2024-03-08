library(tidyverse)
library(janitor)

# site and plot visits =========================================================
# pha 2011 trt, esv 2012
site_visits <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                            sheet = "Export_TBL_Site_Visit") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_Site_Visit"))
plot_visits <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                 sheet = "Export_TBL_PlotVisit") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_PlotVisit")) 

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
  print(n=21)

plot_visits_10y <- plot_visits|>
  filter(trt_year != 2017);glimpse(plot_visits)

filter(plot_visits, is.na(yst)) |> glimpse() # 4 plots not sampled

pull(plot_visits, yst) |> unique()

lut_phase <- pull(plot_visits, VisitCode)
names(lut_phase) <- pull(plot_visits, phase_adj)

new_vcs<- dplyr::select(plot_visits, VisitCode, new_visit_code, phase_adj, trt_year) |> unique()
table(plot_visits$phase_adj)

plot_visits |>
  group_by(phase_adj, trt_year) |>
  summarise(n_plots = length(unique(PlotCode)))


# khr fuels ====================================================================
khr_fuels <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                 sheet = "Export_TBL_1KhrFuel") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_1KhrFuel")) |>
  left_join(new_vcs); glimpse(khr_fuels)

khr_fuels |>
  group_by(phase_adj, trt_year, PlotSize) |>
  summarise(n_plots = length(unique(PlotCode)))

# centerPlotTreeData ============

cp_tree <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                               sheet = "Export_TBL_CenterPlotTreeData") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_CenterPlotTreeData")) |>
  left_join(new_vcs); glimpse(cp_tree)

cp_tree |>
  group_by(phase_adj, trt_year) |>
  summarise(n_plots = length(unique(PlotCode)))
dim(cp_tree)

ggpubr::ggarrange(
cp_tree |>
  filter(Species %in% c("PIPO", "PSME")) |>
  mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus)) |>
  ggplot(aes(x=phase_adj, y=DBH, fill=PlotTreatmentStatus)) +
  geom_boxplot() +
  facet_wrap(trt_year~Species) +
  ggtitle("CenterPlotTreeData, DBH", "no TagIDs, so we cant track individual trees")
,
cp_tree |>
  filter(Species %in% c("PIPO", "PSME")) |>
  mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus)) |>
  ggplot() +
  geom_boxplot(aes(x=phase_adj, y=Height, fill=PlotTreatmentStatus)) +
  facet_wrap(trt_year~Species) +
  ggtitle("CenterPlotTreeData, Height", "no TagIDs, so we cant track individual trees")
, common.legend = TRUE, nrow=1)
ggsave("out/centerPlotTreeData.png", width=15, height=8)

# center plot sapling ==========================================================
cp_sapling <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                             sheet = "Export_TBL_CenterPlotSapling") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_CenterPlotSapling")) |>
  left_join(new_vcs); glimpse(cp_sapling)

cp_sapling |>
  group_by(phase_adj, trt_year) |>
  summarise(n_plots = length(unique(PlotCode)))
dim(cp_sapling)


cp_sapling |>
  group_by(new_visit_code, trt_year, phase_adj, PlotTreatmentStatus) |>
  summarise(stems = sum(StemCount),
            mean_height = mean(Height, na.rm=T),
            mean_dbh = mean(DBH, na.rm=T)) |>
  ungroup() |>
  pivot_longer(cols = c(stems, mean_dbh, mean_height)) |>
  ggplot(aes(x=phase_adj, y=value, fill = PlotTreatmentStatus)) +
  geom_boxplot() +
  facet_wrap(trt_year~name, scales = "free_y") +
  ggtitle("Saplings")
ggsave("out/centerPlotSaplings.png", width=10, height=6)

# seedlings ====================================================================
sp_seedling <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                sheet = "Export_TBL_SubPlotSeedling") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_SubPlotSeedling")) |>
  left_join(new_vcs); glimpse(sp_seedling)

sp_seedling |>
  group_by(new_visit_code, trt_year, phase_adj, PlotTreatmentStatus) |>
  summarise(TotalCount = sum(TotalCount)) |>
  ungroup() |>
  mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated",
                                      "Control", PlotTreatmentStatus)) |>
  
  ggplot(aes(x=TotalCount +1, fill=PlotTreatmentStatus)) +
  geom_density(alpha=0.5) +
  facet_grid(trt_year~phase_adj) +
  scale_x_log10() +
  ggtitle("Plot seedling counts (all species summed)")
ggsave("out/SubPlotSeedlings.png", width=10, height=6)

# ground fuel ==================================================================



# plant understory =============================================================
esv <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx",
                              sheet = "Export_TBL_MS_Spokes")
pha <- readxl::read_xlsx("data/PHA_Export_20231024.xlsx",
                                  sheet = "Export_TBL_MS_Spokes")
comm_raw <- bind_rows(esv, pha)|>
  filter(TreatmentUnit %in% plot_visits$TreatmentUnit) |>
  left_join(new_vcs)
glimpse(comm_raw)
non_plant_codes <- c("FWD", "Fuel in Air", "Substrate")

meta<- comm_raw %>%
  dplyr::select(VisitCode, PlotTreatmentStatus, TreatmentUnit) %>%
  mutate(PlotTreatmentStatus=ifelse(PlotTreatmentStatus == "NotTreated",
                                    "Control", PlotTreatmentStatus)) |>
  left_join(new_vcs) |>
  tidyr::separate(new_visit_code, into = c("PlotCode", "Phase"), sep = "\\.",remove = F)

cover_plot <- comm_raw %>%
  mutate(Cover_Top = ifelse(Transect == "P", .1, Cover_Top))%>%
  # filter(Transect != "P") %>%
  group_by(new_visit_code) %>%
  mutate(n_transects = length(unique(Transect))) %>%
  ungroup() %>%
  group_by(new_visit_code, CodeType, SpeciesCode, Transect, n_transects) %>%
  summarise(cover = (Cover_Top + Cover_Bottom)) %>%
  ungroup() %>%
  group_by(new_visit_code, CodeType, SpeciesCode) %>%
  summarise(cover = sum(cover)/n_transects) %>%
  ungroup() %>%
  unique();cover_plot


comm <- cover_plot %>%
  filter(!CodeType %in% non_plant_codes) %>%
  dplyr::select(-CodeType) %>%
  pivot_wider(names_from = SpeciesCode,
              values_from = cover,# values_fn = function(x)sum(x, na.rm=T),
              values_fill = 0) %>%
   dplyr::mutate_all(list(~ replace(., is.na(.), 0))) %>%
  tibble::column_to_rownames("new_visit_code") %>%
  filter(rowSums(.)>0)

rowSums(comm)
dim(comm)
comm <- comm[,colSums(comm)>0];dim(comm) # 20 species have 0 occurrences


cover_by_type <-
  cover_plot |>
  group_by(new_visit_code, CodeType) |>
  summarise(cover = sum(cover, na.rm=T)) |>
  ungroup()

unique(cover_by_type$CodeType)

herbaceous_cover <-
  cover_by_type |>
  mutate(code_modified = ifelse(CodeType %in% c("Forb", "Graminoid"), "Herbaceous", CodeType))  |>
  group_by(new_visit_code, code_modified) |>
  summarise(cover = sum(cover, na.rm=T)) |>
  ungroup()

write_csv(herbaceous_cover, "data/herb_cover_by_plot.csv")
write_csv(cover_by_type, "data/cover_by_plot.csv")
