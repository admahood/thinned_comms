library(tidyverse)
library(janitor)
esv <- readxl::read_xlsx("data/ESV_Export_20231012.xlsx",
                              sheet = "Export_TBL_MS_Spokes")
pha <- readxl::read_xlsx("data/PHA_Export_20231024.xlsx",
                                  sheet = "Export_TBL_MS_Spokes")
comm_raw <- bind_rows(esv, pha) 
glimpse(comm_raw)
non_plant_codes <- c("FWD", "Fuel in Air", "Substrate")

meta<- comm_raw %>%
  dplyr::select(VisitCode, PlotCode, Phase, PlotTreatmentStatus, TreatmentUnit) %>%
  dplyr::mutate(Phase = factor(Phase, levels = c("Pre", "Post", "Post2", "Post3"),
                               ordered = TRUE)) %>%
  mutate(PlotTreatmentStatus=ifelse(PlotTreatmentStatus == "NotTreated",
                                    "Control", PlotTreatmentStatus))



cover_plot <- comm_raw %>%
  mutate(Cover_Top = ifelse(Transect == "P", .1, Cover_Top))%>%
  # filter(Transect != "P") %>%
  group_by(VisitCode) %>%
  mutate(n_transects = length(unique(Transect))) %>%
  ungroup() %>%
  group_by(VisitCode, CodeType, SpeciesCode, Transect, n_transects) %>%
  summarise(cover = (Cover_Top + Cover_Bottom)) %>%
  ungroup() %>%
  group_by(VisitCode, CodeType, SpeciesCode) %>%
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
  tibble::column_to_rownames("VisitCode") %>%
  filter(rowSums(.)>0)

rowSums(comm)
dim(comm)
comm <- comm[,colSums(comm)>0];dim(comm) # 20 species have 0 occurrences


cover_by_type <-
  cover_plot |>
  group_by(VisitCode, CodeType) |>
  summarise(cover = sum(cover, na.rm=T)) |>
  ungroup()

unique(cover_by_type$CodeType)

herbaceous_cover <-
  cover_by_type |>
  mutate(code_modified = ifelse(CodeType %in% c("Forb", "Graminoid"), "Herbaceous", CodeType))  |>
  group_by(VisitCode, code_modified) |>
  summarise(cover = sum(cover, na.rm=T)) |>
  ungroup()

write_csv(herbaceous_cover, "data/herb_cover_by_plot.csv")
write_csv(cover_by_type, "data/cover_by_plot.csv")
