library(tidyverse)
library(janitor)

# site and plot visits =========================================================
# pha 2011 trt, esv 2012
site_visits <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                            sheet = "Export_TBL_Site_Visit") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_Site_Visit"))

treatments <- plot_visits |>
  dplyr::select(TreatmentDate1, TreatmentUnit) |>
  unique() |>
  na.omit() |>
  bind_rows(tibble(TreatmentDate1 = rep(2017, 3), TreatmentUnit = c("PHA-2C1", "PHA-2C2", "PHA-2C3"))) |>
  bind_rows(tibble(TreatmentDate1 = rep(2012, 2), TreatmentUnit = c("PHA-1-2", "PHA-1-3")))

lut_trtd <- pull(treatments, TreatmentDate1)
names(lut_trtd) <- pull(treatments, TreatmentUnit)


plot_visits <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                         sheet = "Export_TBL_PlotVisit") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_PlotVisit")) |>
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
  summarise(n_plots = length(unique(PlotCode)))

plot_visits <- plot_visits|>
  filter(trt_year != 2017);glimpse(plot_visits)

filter(plot_visits, is.na(yst)) |> glimpse() # 4 plots not sampled

pull(plot_visits, yst) |> unique()

lut_phase <- pull(plot_visits, VisitCode)
names(lut_phase) <- pull(plot_visits, phase_adj)

new_vcs<- dplyr::select(plot_visits, VisitCode, new_visit_code) |> unique()
table(plot_visits$phase_adj)

plot_visits


# khr fuels ====================================================================
khr_fuels <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx", 
                                 sheet = "Export_TBL_1KhrFuel") |>
  bind_rows(readxl::read_xlsx("data/PHA_Export_20231024.xlsx", 
                              sheet = "Export_TBL_1KhrFuel")) |>
  left_join(new_vcs); glimpse(khr_fuels)


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
