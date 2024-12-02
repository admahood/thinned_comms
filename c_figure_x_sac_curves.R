# species accumulation curves =============
source("R/a_tc_data_prep.R")
library(vegan)
species_list <- read_csv("data/species_list_alm_modified.csv") |> 
  unique() |> # one duplicate in species list
  dplyr::select(-notes)

exotics <-
  species_list |>
  filter(NativityL48 == "exotic") |>
  mutate(SpeciesCode = str_to_lower(SpeciesCode))

natives <-
  species_list |>
  filter(NativityL48 == "native") |>
  mutate(SpeciesCode = str_to_lower(SpeciesCode))

# defining row and col indexes =================================================
rowz <- rownames(comm) |>
  as_tibble() |>
  dplyr::rename(new_visit_code = value) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c("plot", "phase"), remove = F) %>%
  dplyr::mutate(index = 1:nrow(.)) |>
  left_join(plot_visits)

colz <- colnames(comm) |>
  as_tibble() |>
  dplyr::rename(SpeciesCode = value) |>
  left_join(sp_list |> dplyr::mutate(SpeciesCode = str_to_lower(SpeciesCode))) %>%
  dplyr::mutate(index = 1:nrow(.))

# phases
ph0c <- rowz |> dplyr::filter(phase == "01_Pre", PlotTreatmentStatus == "Control") |> pull(index)
ph1c <- rowz |> dplyr::filter(phase == "post0-1", PlotTreatmentStatus == "Control") |> pull(index)
ph5c <- rowz |> dplyr::filter(phase == "post04-5", PlotTreatmentStatus == "Control") |> pull(index)
ph10c <- rowz |> dplyr::filter(phase == "post10-11", PlotTreatmentStatus == "Control") |> pull(index)

ph0t <- rowz |> dplyr::filter(phase == "01_Pre", PlotTreatmentStatus == "Treatment") |> pull(index)
ph1t <- rowz |> dplyr::filter(phase == "post0-1", PlotTreatmentStatus == "Treatment") |> pull(index)
ph5t <- rowz |> dplyr::filter(phase == "post04-5", PlotTreatmentStatus == "Treatment") |> pull(index)
ph10t <- rowz |> dplyr::filter(phase == "post10-11", PlotTreatmentStatus == "Treatment") |> pull(index)

# exotics 

Exotics <- colz |> dplyr::filter(NativityL48 == "exotic") |> pull(SpeciesCode)
Natives <- colz |> dplyr::filter(NativityL48 == "native") |> pull(SpeciesCode)

# SAC plot (in the least pretty way possible) ==================================

# df 
bind_rows(
  # vegan::specaccum(comm[ph0c,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "total", ph = "0"),
  # vegan::specaccum(comm[ph1c,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "total", ph = "01"),
  # vegan::specaccum(comm[ph5c,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "total", ph = "05"),
  # vegan::specaccum(comm[ph10c,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "total", ph = "10"),
  # vegan::specaccum(comm[ph0t,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "total", ph = "0"),
  # vegan::specaccum(comm[ph1t,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "total", ph = "01"),
  # vegan::specaccum(comm[ph5t,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "total", ph = "05"),
  # vegan::specaccum(comm[ph10t,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "total", ph = "10"),
  vegan::specaccum(comm[ph0c,Natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Control", spp = "Native", ph = "0"),
  vegan::specaccum(comm[ph1c,Natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Control", spp = "Native", ph = "01"),
  vegan::specaccum(comm[ph5c,Natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Control", spp = "Native", ph = "05"),
  vegan::specaccum(comm[ph10c,Natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Control", spp = "Native", ph = "10"),
  vegan::specaccum(comm[ph0t,Natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Treated", spp = "Native", ph = "0"),
  vegan::specaccum(comm[ph1t,Natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Treated", spp = "Native", ph = "01"),
  vegan::specaccum(comm[ph5t,Natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Treated", spp = "Native", ph = "05"),
  vegan::specaccum(comm[ph10t,Natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Treated", spp = "Native", ph = "10"),
  vegan::specaccum(comm[ph0c,Exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Control", spp = "Exotic", ph = "0"),
  vegan::specaccum(comm[ph1c,Exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Control", spp = "Exotic", ph = "01"),
  vegan::specaccum(comm[ph5c,Exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Control", spp = "Exotic", ph = "05"),
  vegan::specaccum(comm[ph10c,Exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Control", spp = "Exotic", ph = "10"),
  vegan::specaccum(comm[ph0t,Exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Treated", spp = "Exotic", ph = "0"),
  vegan::specaccum(comm[ph1t,Exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Treated", spp = "Exotic", ph = "01"),
  vegan::specaccum(comm[ph5t,Exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Treated", spp = "Exotic", ph = "05"),
  vegan::specaccum(comm[ph10t,Exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="Treated", spp = "Exotic", ph = "10"),
) |>
  ggplot(aes(x=sites, color = ph)) +
  geom_line(aes(y=richness)) +
  geom_linerange(aes(ymin=richness-sd, ymax = richness+sd)) +
  scale_color_brewer(palette = "Set1", name = "Years\nSince\nTreatment") +
  facet_grid(spp~trt, scales = "free") +
  ylab("Species Richness") +
  xlab("Number of Plots") +
  theme_bw()
ggsave("out/spp_accumulation.png", width = 7, height=5, bg="white")


# species pool estimations
spp_pools <- 
bind_rows(
  vegan::specpool(comm[ph0c,]) |> as_tibble() |> mutate(trt="Control", spp = "total", ph = "0"),
  vegan::specpool(comm[ph1c,]) |> as_tibble() |> mutate(trt="Control", spp = "total", ph = "01"),
  vegan::specpool(comm[ph5c,]) |> as_tibble() |> mutate(trt="Control", spp = "total", ph = "05"),
  vegan::specpool(comm[ph10c,]) |> as_tibble() |> mutate(trt="Control", spp = "total", ph = "10"),
  vegan::specpool(comm[ph0t,]) |> as_tibble() |> mutate(trt="Treated", spp = "total", ph = "0"),
  vegan::specpool(comm[ph1t,]) |> as_tibble() |> mutate(trt="Treated", spp = "total", ph = "01"),
  vegan::specpool(comm[ph5t,]) |> as_tibble() |> mutate(trt="Treated", spp = "total", ph = "05"),
  vegan::specpool(comm[ph10t,]) |> as_tibble() |> mutate(trt="Treated", spp = "total", ph = "10"),
  vegan::specpool(comm[ph0c,Natives]) |> as_tibble() |> mutate(trt="Control", spp = "Native", ph = "0"),
  vegan::specpool(comm[ph1c,Natives]) |> as_tibble() |> mutate(trt="Control", spp = "Native", ph = "01"),
  vegan::specpool(comm[ph5c,Natives]) |> as_tibble() |> mutate(trt="Control", spp = "Native", ph = "05"),
  vegan::specpool(comm[ph10c,Natives]) |> as_tibble() |> mutate(trt="Control", spp = "Native", ph = "10"),
  vegan::specpool(comm[ph0t,Natives]) |> as_tibble() |> mutate(trt="Treated", spp = "Native", ph = "0"),
  vegan::specpool(comm[ph1t,Natives]) |> as_tibble() |> mutate(trt="Treated", spp = "Native", ph = "01"),
  vegan::specpool(comm[ph5t,Natives]) |> as_tibble() |> mutate(trt="Treated", spp = "Native", ph = "05"),
  vegan::specpool(comm[ph10t,Natives]) |> as_tibble() |> mutate(trt="Treated", spp = "Native", ph = "10"),
  vegan::specpool(comm[ph0c,Exotics]) |> as_tibble() |> mutate(trt="Control", spp = "Exotic", ph = "0"),
  vegan::specpool(comm[ph1c,Exotics]) |> as_tibble() |> mutate(trt="Control", spp = "Exotic", ph = "01"),
  vegan::specpool(comm[ph5c,Exotics]) |> as_tibble() |> mutate(trt="Control", spp = "Exotic", ph = "05"),
  vegan::specpool(comm[ph10c,Exotics]) |> as_tibble() |> mutate(trt="Control", spp = "Exotic", ph = "10"),
  vegan::specpool(comm[ph0t,Exotics]) |> as_tibble() |> mutate(trt="Treated", spp = "Exotic", ph = "0"),
  vegan::specpool(comm[ph1t,Exotics]) |> as_tibble() |> mutate(trt="Treated", spp = "Exotic", ph = "01"),
  vegan::specpool(comm[ph5t,Exotics]) |> as_tibble() |> mutate(trt="Treated", spp = "Exotic", ph = "05"),
  vegan::specpool(comm[ph10t,Exotics]) |> as_tibble() |> mutate(trt="Treated", spp = "Exotic", ph = "10")
) |>
  dplyr::rename(num_plots = n, treatment = trt, phase = ph) |>
  print(n=24)

spp_pools |>
  dplyr::select(boot, boot.se, treatment, spp, phase) |>
  pivot_wider(names_from = spp, values_from = c(boot, boot.se), id_cols = c(phase, treatment)) |>
  dplyr::select(1,2,3,6,4,7,5,8)|>
  # arrange(phase, treatment) |>
  transmute(phase = phase, treatment = treatment,
            total_richness = paste0(round(boot_total), " (", round(boot.se_total), ")"),
            native_richness = paste0(round(boot_Native), " (", round(boot.se_Native), ")"),
            exotic_richness = paste0(round(boot_Exotic), " (", round(boot.se_Exotic,1), ")")) |>
  write_csv("out/ext_spp_pool_estimations.csv")
