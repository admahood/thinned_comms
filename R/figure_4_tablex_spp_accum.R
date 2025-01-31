# basic community analysis

source("R/a_tc_data_prep.R")
library(vegan)
species_list <- read_csv("data/species_list_alm_modified.csv") |> 
  unique() |> # one duplicate in species list
  dplyr::select(-notes)



# which species were there before - ie are introductions after treatment from the regional pool or outside?
exotics <-
  species_list |>
  filter(NativityL48 == "exotic") |>
  mutate(SpeciesCode = str_to_lower(SpeciesCode))

c_vs_t <- plot_visits |>
  dplyr::select(PlotCode, PlotTreatmentStatus) |>
  mutate(PlotTreatmentStatus = ifelse(PlotTreatmentStatus == "NotTreated", "Control", PlotTreatmentStatus)) |>
  unique()

exotic_prevalence <- comm[, names(comm) %in% exotics$SpeciesCode] |>
  tibble::rownames_to_column('new_visit_code') |>
  tidyr::separate(new_visit_code, into = c('PlotCode', 'phase_adj'), sep = "\\.") |>
  pivot_longer(-c(PlotCode, phase_adj)) |>
  mutate(occurrence = ifelse(value >0, 1, 0)) |>
  left_join(c_vs_t) |>
  group_by(phase_adj, PlotTreatmentStatus, name) |>
  summarise(prevalence = sum(occurrence)) |>
  ungroup() |>
  mutate(phase_adj = str_c(phase_adj, "_", PlotTreatmentStatus)) |>
  dplyr::select(-PlotTreatmentStatus) |>
  pivot_wider(names_from = phase_adj, values_from = prevalence) |>
  left_join(exotics |> dplyr::rename('name' = 'SpeciesCode')) |>
  dplyr::select(FinalName, 
                starts_with('01'), starts_with('post'))

write_csv(exotic_prevalence, file = "out/exotic_prevalence.csv")

ggplot(exotic_only, aes(x=phase_adj, y=occurrence, color = name)) +
  geom_boxplot() +
  facet_wrap(~name, scales = "free_y")

# number of plots invaded

invasion_freq <- comm[, names(comm) %in% exotics$SpeciesCode] |>
  tibble::rownames_to_column('new_visit_code') |>
  tidyr::separate(new_visit_code, into = c('PlotCode', 'phase_adj'), sep = "\\.") |>
  pivot_longer(-c(PlotCode, phase_adj)) |>
  mutate(occurrence = ifelse(value >0, 1, 0)) |>
  left_join(c_vs_t) |>
  group_by(phase_adj, PlotTreatmentStatus, PlotCode) |>
  summarise(invaded = ifelse(sum(occurrence)>0, 1, 0)) |>
  ungroup() |>
  group_by(phase_adj, PlotTreatmentStatus) |>
  summarise(n_invaded_plots = sum(invaded),
            n_plots = n()) |>
  ungroup() |>
  mutate(percent_invaded = n_invaded_plots/n_plots *100) ##|>
  # dplyr::select(-n_invaded_plots, -n_plots) |>
  # pivot_wider(names_from = phase_adj, values_from = percent_invaded) |>
  # janitor::clean_names()
write_csv(invasion_freq, "out/invasion_freq.csv")

# do difference from pre-treatment
cp <- cover_plot |>
  filter(!CodeType %in% non_plant_codes) |>
  tidyr::separate(new_visit_code, sep = "\\.", into = c('plot', 'phase'), remove = F) |>
  mutate(site = str_sub(plot, 1,1))

comm |> vegan::diversity("shannon")
comm |> vegan::specnumber()

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

exotics <- colz |> dplyr::filter(NativityL48 == "exotic") |> pull(SpeciesCode)
natives <- colz |> dplyr::filter(NativityL48 == "native") |> pull(SpeciesCode)

# species accumulation curves =========================================

vegan::specaccum(comm[ph5c,]) |> plot(col = "red", main = "Control")
vegan::specaccum(comm[ph0c,]) |> plot(add=T)
vegan::specaccum(comm[ph1c,]) |> plot(col = "gold", add=T)
vegan::specaccum(comm[ph10c,]) |> plot(col = "firebrick", add=T)
vegan::specaccum(comm[ph5t,]) |> plot(col = "red", main = "Treated")
vegan::specaccum(comm[ph0t,]) |> plot(add=T)
vegan::specaccum(comm[ph1t,]) |> plot(col = "gold", add=T)
vegan::specaccum(comm[ph10t,]) |> plot(col = "firebrick", add=T)

vegan::specaccum(comm[ph5c,natives]) |> plot(col = "red", main = "Control, natives")
vegan::specaccum(comm[ph0c,natives]) |> plot(add=T)
vegan::specaccum(comm[ph1c,natives]) |> plot(col = "gold", add=T)
vegan::specaccum(comm[ph10c,natives]) |> plot(col = "firebrick", add=T)
vegan::specaccum(comm[ph5t,natives]) |> plot(col = "red", main = "Treated, natives")
vegan::specaccum(comm[ph0t,natives]) |> plot(add=T)
vegan::specaccum(comm[ph1t,natives]) |> plot(col = "gold", add=T)
vegan::specaccum(comm[ph10t,natives]) |> plot(col = "firebrick", add=T)

vegan::specaccum(comm[ph10c,exotics]) |> plot(col = "firebrick", main = "Control, exotics")
vegan::specaccum(comm[ph5c,exotics]) |> plot(col = "red", add=T)
vegan::specaccum(comm[ph0c,exotics]) |> plot(add=T)
vegan::specaccum(comm[ph1c,exotics]) |> plot(col = "gold", add=T)

vegan::specaccum(comm[ph5t,exotics]) |> plot(col = "red", main = "Treated, exotics")
vegan::specaccum(comm[ph0t,exotics]) |> plot(add=T)
vegan::specaccum(comm[ph1t,exotics]) |> plot(col = "gold", add=T)
vegan::specaccum(comm[ph10t,exotics]) |> plot(col = "firebrick", add=T)

# df 
bind_rows(
vegan::specaccum(comm[ph0c,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "total", ph = "0"),
vegan::specaccum(comm[ph1c,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "total", ph = "01"),
vegan::specaccum(comm[ph5c,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "total", ph = "05"),
vegan::specaccum(comm[ph10c,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "total", ph = "10"),
vegan::specaccum(comm[ph0t,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "total", ph = "0"),
vegan::specaccum(comm[ph1t,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "total", ph = "01"),
vegan::specaccum(comm[ph5t,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "total", ph = "05"),
vegan::specaccum(comm[ph10t,], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "total", ph = "10"),
vegan::specaccum(comm[ph0c,natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "native", ph = "0"),
vegan::specaccum(comm[ph1c,natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "native", ph = "01"),
vegan::specaccum(comm[ph5c,natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "native", ph = "05"),
vegan::specaccum(comm[ph10c,natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "native", ph = "10"),
vegan::specaccum(comm[ph0t,natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "native", ph = "0"),
vegan::specaccum(comm[ph1t,natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "native", ph = "01"),
vegan::specaccum(comm[ph5t,natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "native", ph = "05"),
vegan::specaccum(comm[ph10t,natives], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "native", ph = "10"),
vegan::specaccum(comm[ph0c,exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "exotic", ph = "0"),
vegan::specaccum(comm[ph1c,exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "exotic", ph = "01"),
vegan::specaccum(comm[ph5c,exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "exotic", ph = "05"),
vegan::specaccum(comm[ph10c,exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="control", spp = "exotic", ph = "10"),
vegan::specaccum(comm[ph0t,exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "exotic", ph = "0"),
vegan::specaccum(comm[ph1t,exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "exotic", ph = "01"),
vegan::specaccum(comm[ph5t,exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "exotic", ph = "05"),
vegan::specaccum(comm[ph10t,exotics], permutations = 999)[c(3:5)] |> as_tibble() |> mutate(trt="treated", spp = "exotic", ph = "10"),
) |>
  ggplot(aes(x=sites, color = ph)) +
  geom_line(aes(y=richness)) +
  geom_linerange(aes(ymin=richness-sd, ymax = richness+sd)) +
  scale_color_brewer(palette = "RdYlOr") +
  facet_grid(spp~trt, scales = "free") +
  theme_bw()
ggsave("out/spp_accumulation.png", width = 7, height=5, bg="white")
# species pools ============
cpool <- bind_rows(
  vegan::specpool(comm[ph0c,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "pre"),
  vegan::specpool(comm[ph1c,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01"),
  vegan::specpool(comm[ph5c,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05"),
  vegan::specpool(comm[ph10c,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10")
)

pools <- bind_rows(
  vegan::specpool(comm[ph0t,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "pre"),
  vegan::specpool(comm[ph1t,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01"),
  vegan::specpool(comm[ph5t,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05"),
  vegan::specpool(comm[ph10t,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10")
) |> bind_rows(cpool)

write_csv(pools,"out/specpools.csv")

pools |>
  dplyr::mutate(chao = paste0(round(chao), " (", round(chao.se), ")")) |>
  dplyr::select(PlotTreatmentStatus, phase, chao) |>
  pivot_wider(names_from = PlotTreatmentStatus, values_from = chao) |>
  write_csv("out/specpools_simple.csv")

# div indexes

bind_rows(
  vegan::diversity(comm[ph0t,])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "0pre"),
  vegan::diversity(comm[ph1t,])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01"),
  vegan::diversity(comm[ph5t,])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05"),
  vegan::diversity(comm[ph10t,])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10"),
  vegan::diversity(comm[ph0c,])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "0pre"),
  vegan::diversity(comm[ph1c,])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01"),
  vegan::diversity(comm[ph5c,])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05"),
  vegan::diversity(comm[ph10c,])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10")) |>
  ggplot(aes(x=phase, y=value)) +
  geom_boxplot() +
  facet_wrap(~PlotTreatmentStatus) +
  ggtitle("Shannon-Weiner Alpha Diversity (total)")


bind_rows(
  vegan::diversity(comm[ph0t, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "0pre"),
  vegan::diversity(comm[ph1t, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01"),
  vegan::diversity(comm[ph5t, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05"),
  vegan::diversity(comm[ph10t, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10"),
  vegan::diversity(comm[ph0c, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "0pre"),
  vegan::diversity(comm[ph1c, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01"),
  vegan::diversity(comm[ph5c, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05"),
  vegan::diversity(comm[ph10c, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10")) |>
  ggplot(aes(x=phase, y=value)) +
  geom_boxplot() +
  facet_grid(str_sub(visitcode,1,1)~PlotTreatmentStatus) +
  ggtitle("Shannon-Weiner Alpha Diversity (native)")


bind_rows(
  vegan::diversity(comm[ph0t, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "0pre"),
  vegan::diversity(comm[ph1t, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01"),
  vegan::diversity(comm[ph5t, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05"),
  vegan::diversity(comm[ph10t, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10"),
  vegan::diversity(comm[ph0c, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "0pre"),
  vegan::diversity(comm[ph1c, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01"),
  vegan::diversity(comm[ph5c, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05"),
  vegan::diversity(comm[ph10c, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10")) |>
  ggplot(aes(x=phase, y=value)) +
  geom_boxplot() +
  facet_grid(str_sub(visitcode,1,1)~PlotTreatmentStatus) +
  ggtitle("Shannon-Weiner Alpha Diversity (exotic)")

bind_rows(
  vegan::specnumber(comm[ph0t, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "0pre"),
  vegan::specnumber(comm[ph1t, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01"),
  vegan::specnumber(comm[ph5t, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05"),
  vegan::specnumber(comm[ph10t, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10"),
  vegan::specnumber(comm[ph0c, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "0pre"),
  vegan::specnumber(comm[ph1c, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01"),
  vegan::specnumber(comm[ph5c, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05"),
  vegan::specnumber(comm[ph10c, natives])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10")) |>
  ggplot(aes(x=phase, y=value)) +
  geom_boxplot() +
  facet_wrap(~PlotTreatmentStatus) +
  ggtitle("Richness (native)")

bind_rows(
  vegan::specnumber(comm[ph0t, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "0pre"),
  vegan::specnumber(comm[ph1t, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01"),
  vegan::specnumber(comm[ph5t, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05"),
  vegan::specnumber(comm[ph10t, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10"),
  vegan::specnumber(comm[ph0c, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "0pre"),
  vegan::specnumber(comm[ph1c, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01"),
  vegan::specnumber(comm[ph5c, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05"),
  vegan::specnumber(comm[ph10c, exotics])|> as_tibble(rownames = "visitcode") |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10")) |>
  ggplot(aes(x=phase, y=value)) +
  geom_boxplot() +
  facet_wrap(~PlotTreatmentStatus) +
  ggtitle("Richness (exotic)")

# basic boxplots ===============================================================
cp |>
  ggplot(aes(x=phase, y=cover, fill = CodeType)) +
  geom_boxplot(outliers = F) +
  facet_grid(site ~ PlotTreatmentStatus)

cp |>
  pivot_wider(names_from = phase, values_from = cover, values_fill = 0) |>
  janitor::clean_names() |>
  mutate(d01 = post0_1 - x01_pre, 
         d05 = post04_5 - x01_pre, 
         d10 = post10_11 - x01_pre) |>
  dplyr::select(-x01_pre, -starts_with("post")) |>
  pivot_longer(-c(plot, code_type, species_code, plot_treatment_status, site)) |>
  filter(!code_type %in% non_plant_codes) |>
  left_join(species_list |> 
              mutate(SpeciesCode = str_to_lower(SpeciesCode)) |>
              janitor::clean_names()) |>
  mutate(fg = str_c(nativity_l48, "_", code_type))|>
  filter(str_sub(fg, 1,1) != "U")|>
  ggplot(aes(x=name, y=value, fill = fg)) +
  geom_boxplot(color = "black", position = "dodge", outliers = F) +
  facet_grid(site ~ plot_treatment_status)


spp <- cp |>
  filter(!CodeType %in% non_plant_codes) |>
  left_join(species_list |> mutate(SpeciesCode = str_to_lower(SpeciesCode))) |>
  left_join(plot_visits)


spp |>
  mutate(fg = str_c(NativityL48, "_", CodeType)) |>
  filter(str_sub(fg, 1,1) != "U") |>
  # group_by(fg, phase, site,PlotTreatmentStatus) |>
  # summarise(cover = sum(cover)) |>
  ggplot(aes(x=phase, y=cover, fill = fg)) +
  geom_boxplot(color = "black", position = "dodge", outliers = F) +
  facet_grid(site ~ PlotTreatmentStatus)


# calculate total cover, cover of natives/exotics

# total_cover
total_cover <- spp |>
  group_by(plot, phase, PlotTreatmentStatus, site, trt_u_adj) |>
  summarise(total_cover = sum(cover)) |>
  ungroup() |>
  tidyr::replace_na(list('trt_u_adj' = "PHA-1-2"))

native_cover <- spp |>
  filter(NativityL48 == "native") |>
  group_by(plot, phase, PlotTreatmentStatus, site, trt_u_adj) |>
  summarise(total_cover = sum(cover)) |>
  ungroup() |>
  tidyr::replace_na(list('trt_u_adj' = "PHA-1-2"))

exotic_cover <- spp |>
  filter(NativityL48 == "exotic") |>
  group_by(plot, phase, PlotTreatmentStatus, site, trt_u_adj) |>
  summarise(total_cover = sum(cover)) |>
  ungroup() |>
  tidyr::replace_na(list('trt_u_adj' = "PHA-1-2"))

spp |>
  group_by(plot, phase, PlotTreatmentStatus, site) |>
  summarise(total_cover = sum(cover)) |>
  ungroup() |>
  ggplot(aes(x=phase, y=total_cover, fill=PlotTreatmentStatus)) +
  geom_boxplot() +
  facet_wrap(~site)

spp |>
  filter(NativityL48 != "Unknown") |>
  group_by(plot, phase, PlotTreatmentStatus, site, NativityL48) |>
  summarise(total_cover = sum(cover)) |>
  ungroup() |>
  ggplot(aes(x=phase, y=total_cover, fill = PlotTreatmentStatus)) +
  geom_boxplot(outliers=F) +
  # geom_jitter(aes(shape = site, color = PlotTreatmentStatus), width = .1) +
  facet_wrap(~ NativityL48, scales="free")
ggsave("out/native_exotic_cover.png", height=4, width=7, bg="white")
# cohen's d effect size
library(effsize)
?effsize::cohen.d()

ddf <- data.frame(trt_u = NA, d = NA, cil = NA, ciu=NA, phase = NA, trt = NA, 
                  magnitude = NA)
outlist <- list()
counter = 1
for(i in unique(total_cover$trt_u_adj)){
  tc <- filter(total_cover, trt_u_adj == i)
  for(j in c("post0-1", "post04-5", "post10-11")){
    for(t in c("Control", "Treatment")){
    cohend <- effsize::cohen.d(tc |> filter(phase == j, PlotTreatmentStatus == t) |> 
                              pull(total_cover),
                            tc |> filter(phase == "01_Pre", PlotTreatmentStatus == t) |> 
                              pull(total_cover)) 
    ddf[counter, 1] <- i
    ddf[counter, 2] <- cohend$estimate
    ddf[counter, 3] <- cohend$conf.int[1]
    ddf[counter, 4] <- cohend$conf.int[2]
    ddf[counter, 5] <- j
    ddf[counter, 6] <- t
    ddf[counter, 7] <- cohend$magnitude
    outlist[[counter]] <- cohend
    counter = counter +1
    }
  } 
}

ddfn <- data.frame(trt_u = NA, d = NA, cil = NA, ciu=NA, phase = NA, trt = NA, 
                  magnitude = NA)
outlist <- list()
counter = 1
for(i in unique(native_cover$trt_u_adj)){
  tc <- filter(native_cover, trt_u_adj == i)
  for(j in c("post0-1", "post04-5", "post10-11")){
    for(t in c("Control", "Treatment")){
      cohend <- effsize::cohen.d(tc |> filter(phase == j, PlotTreatmentStatus == t) |> 
                                   pull(total_cover),
                                 tc |> filter(phase == "01_Pre", PlotTreatmentStatus == t) |> 
                                   pull(total_cover)) 
      ddfn[counter, 1] <- i
      ddfn[counter, 2] <- cohend$estimate
      ddfn[counter, 3] <- cohend$conf.int[1]
      ddfn[counter, 4] <- cohend$conf.int[2]
      ddfn[counter, 5] <- j
      ddfn[counter, 6] <- t
      ddfn[counter, 7] <- cohend$magnitude
      outlist[[counter]] <- cohend
      counter = counter +1
    }
  } 
}


ddfe <- data.frame(trt_u = NA, d = NA, cil = NA, ciu=NA, phase = NA, trt = NA, 
                   magnitude = NA)
outlist <- list()
counter = 1
for(i in unique(exotic_cover$trt_u_adj)){
  tc <- filter(exotic_cover, trt_u_adj == i)
  for(j in c("post0-1", "post04-5", "post10-11")){
    for(t in c("Control", "Treatment")){
      cohend <- effsize::cohen.d(tc |> filter(phase == j, PlotTreatmentStatus == t) |> 
                                   pull(total_cover),
                                 tc |> filter(phase == "01_Pre", PlotTreatmentStatus == t) |> 
                                   pull(total_cover)) 
      ddfe[counter, 1] <- i
      ddfe[counter, 2] <- cohend$estimate
      ddfe[counter, 3] <- cohend$conf.int[1]
      ddfe[counter, 4] <- cohend$conf.int[2]
      ddfe[counter, 5] <- j
      ddfe[counter, 6] <- t
      ddfe[counter, 7] <- cohend$magnitude
      outlist[[counter]] <- cohend
      counter = counter +1
    }
  } 
}

dddd<- ddf |>
  mutate(spp = "total") |>
  bind_rows(ddfn |> mutate(spp = "native"))|>
  bind_rows(ddfe |> mutate(spp = "exotic"))

ggplot(dddd |> filter(spp != 'total'), aes(x=phase, y=d)) +
  geom_boxplot() +
  geom_jitter(aes(color = trt_u, size=magnitude), width = .1) +
  geom_hline(yintercept = 0) +
  facet_grid(spp~trt) +
  ylab("Effect Size for Total veg cover (d stat)") +
  theme_bw()

ggsave("out/cohen_d_cover.png", bg="white", width=7, height=5)
