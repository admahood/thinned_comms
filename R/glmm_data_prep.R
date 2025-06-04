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

graminoid_cover <- spp |>
  filter(CodeType == "Graminoid")|>
  group_by(plot, phase, PlotTreatmentStatus, site, trt_u_adj, new_visit_code) |>
  summarise(graminoid_cover = sum(cover)) |>
  ungroup() |>
  tidyr::replace_na(list('trt_u_adj' = "PHA-1-2"))
  
  
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
clim <- read_csv('data/plot_climate.csv')

tt <- st_read("data/terrain_by_plot.gpkg") |>
  st_set_geometry(NULL) |> 
  left_join(clim)

fwd <- fuel_cover_by_type |>
  dplyr::filter(CodeType == 'FWD') |>
  dplyr::select(fwd = cover, new_visit_code)

plot_level_metrics <- cp_tree |>
  mutate(site = str_sub(PlotCode, 1,1)) |>
  group_by(PlotCode, PlotTreatmentStatus, phase_adj, new_visit_code, site, TreatmentUnit)|>
  summarise(tpa = sum(CalcCol_DensityContribution, na.rm=T),
            n=n(),
            quadratic_mean_diameter = sqrt(sum(DBH * DBH, na.rm=T)/n)) |>
  left_join(ba_m2ha) |>
  left_join(div_indices)|> 
  left_join(fwd) |>
  left_join(tt) |>
  left_join(graminoid_cover) |>
  left_join(saplings_wide) |> 
  left_join(seedling_density) |>
  replace_na(list(sapling_density_per_acre = 0,
                  sapling_ba_ft_per_acre = 0)) |>
  left_join(native_cover) |>
  left_join(exotic_cover) |>
  replace_na(list(exotic_cover = 0)) |>
  mutate(no_exotics = 100 - exotic_cover,
         # ec_binary = cbind(exotic_cover, no_exotics),
         log_saplingba = log(sapling_ba_ft_per_acre + 1),
         treated = case_when(PlotTreatmentStatus == "Control" ~ "not_treated",
                             PlotTreatmentStatus == "Treatment" & phase_adj == "01_Pre" ~ "not_treated",
                             PlotTreatmentStatus == "Treatment" & phase_adj != "01_Pre" ~ "treated") |> as.factor(),
         invaded = ifelse(exotic_cover > 0, 1, 0),
         exotic_relative_cover = 100 * (exotic_cover/(native_cover + exotic_cover)),
         exotic_relative_richness = 100 * (nspp_exotic/(nspp_native + nspp_exotic))
         ) 
# 
# summary(plot_level_metrics)
# glimpse(plot_level_metrics)
write_csv(plot_level_metrics, "data/plot_level_data.csv")


