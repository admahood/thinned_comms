# basic community analysis

source("R/a_tc_data_prep.R")
library(vegan)
species_list <- read_csv("data/species_list_alm_modified.csv") |> 
  unique() |> # one duplicate in species list
  dplyr::select(-notes)


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


# df .... yes, this is possible to do with like 5 lines of code instead of whatever this is
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
  dplyr::filter(spp !='total') |>
  dplyr::mutate(ph = ifelse(ph == "0", "Pre", ph)|>
                  str_replace_all("05", "5") |>
                  str_replace_all('01', '1') |>
                  forcats::fct_relevel('Pre', '1', '5', '10'),
                spp = str_to_title(spp),
                trt = str_to_title(trt)) |>
  ggplot(aes(x=sites, color = ph)) +
  geom_line(aes(y=richness)) +
  geom_linerange(aes(ymin=richness-sd, ymax = richness+sd)) +
  scale_color_brewer(palette = "Set1", name = "Years\nSince\nTreatment") +
  facet_grid(spp~trt, scales = "free") +
  ylab("Species Richness") +
  xlab("Number of Plots") +
  theme_bw()
ggsave("out/figure_4_spp_accumulation.png", width = 7, height=5, bg="white")

# species pools ============
cpool <- bind_rows(
  vegan::specpool(comm[ph0c,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "0pre", spp = 'total'),
  vegan::specpool(comm[ph1c,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01", spp = 'total'),
  vegan::specpool(comm[ph5c,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05", spp = 'total'),
  vegan::specpool(comm[ph10c,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10", spp = 'total'),
  vegan::specpool(comm[ph0c,exotics]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "0pre", spp = 'exotic'),
  vegan::specpool(comm[ph1c,exotics]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01", spp = 'exotic'),
  vegan::specpool(comm[ph5c,exotics]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05", spp = 'exotic'),
  vegan::specpool(comm[ph10c,exotics]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10", spp = 'exotic'),
  vegan::specpool(comm[ph0c,natives]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "0pre", spp = 'native'),
  vegan::specpool(comm[ph1c,natives]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_01", spp = 'native'),
  vegan::specpool(comm[ph5c,natives]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_05", spp = 'native'),
  vegan::specpool(comm[ph10c,natives]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Control", phase = "post_10", spp = 'native')
  
)

pools <- bind_rows(
  vegan::specpool(comm[ph0t,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "0pre", spp = 'total'),
  vegan::specpool(comm[ph1t,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01", spp = 'total'),
  vegan::specpool(comm[ph5t,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05", spp = 'total'),
  vegan::specpool(comm[ph10t,]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10", spp = 'total'),
  vegan::specpool(comm[ph0t,exotics]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "0pre", spp = 'exotic'),
  vegan::specpool(comm[ph1t,exotics]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01", spp = 'exotic'),
  vegan::specpool(comm[ph5t,exotics]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05", spp = 'exotic'),
  vegan::specpool(comm[ph10t,exotics]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10", spp = 'exotic'),
  vegan::specpool(comm[ph0t,natives]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "0pre", spp = 'native'),
  vegan::specpool(comm[ph1t,natives]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_01", spp = 'native'),
  vegan::specpool(comm[ph5t,natives]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_05", spp = 'native'),
  vegan::specpool(comm[ph10t,natives]) |> as_tibble() |> dplyr::mutate(PlotTreatmentStatus = "Treatment", phase = "post_10", spp = 'native')
) |> bind_rows(cpool)

write_csv(pools,"out/table_sx_specpools_full.csv")

pools |>
  dplyr::mutate(boot = paste0(round(boot), " (", round(boot.se), ")")) |>
  dplyr::select(phase,PlotTreatmentStatus,  boot, spp) |>
  pivot_wider(names_from = spp, values_from = boot) |>
  dplyr::select(phase,PlotTreatmentStatus, total, native, exotic) |>
  arrange(PlotTreatmentStatus, phase) |>
  write_csv("out/table_3_specpools.csv")

pools |>
  dplyr::mutate(boot = paste0(round(boot), " (", round(boot.se), ")")) |>
  dplyr::select(phase,PlotTreatmentStatus,  boot, spp) |>
  pivot_wider(names_from = PlotTreatmentStatus, values_from = boot) |>
  dplyr::select(phase,spp, Control, Treatment) |>
  arrange(spp, phase) |>
  write_csv("out/table_3_specpools_alt.csv")


pools |>
  dplyr::mutate(boot = paste0(round(boot), " (", round(boot.se), ")")) |>
  dplyr::select(phase, spp,PlotTreatmentStatus,  boot) |>
  pivot_wider(names_from = phase, values_from = boot) |>
  arrange(spp, PlotTreatmentStatus) |>
  write_csv("out/table_3_specpools_final.csv")

# NMDS =========================================================================

nmds <- comm |>
  vegan::decostand(method = 'total') |>
  vegan::metaMDS(trymax = 200)

d_nmds <- nmds$points |>
  as_tibble(rownames = 'meta') |>
  tidyr::separate(meta, sep = "\\.", into = c('plot', 'phase')) |>
  mutate(site = ifelse(str_sub(plot, 1,1)=="E", "Estes Valley", "Phantom Creek")) |>
  left_join(plot_visits_10y |> dplyr::select(PlotTreatmentStatus, PlotCode) |> unique(), 
            by = c('plot'='PlotCode')) |>
  mutate(PlotTreatmentStatus = str_replace_all(PlotTreatmentStatus, "NotTreated", "Treatment"))

d_nmds |>
  ggplot(aes(x=MDS1, y=MDS2, color = phase, shape=phase)) +
  geom_path(aes(group = plot), color = 'grey', linewidth = 1) +
  geom_point(size=3) +
  # stat_ellipse(aes(group = site), level = 0.9) +
  facet_grid(site~PlotTreatmentStatus) +
  theme_bw()

pmnva <- vegan::adonis2(comm |>
                          vegan::decostand(method = 'total') ~ PlotTreatmentStatus+phase + site,
                        data = d_nmds, permutations = 999, 
                        by = 'terms')

# change from pre-treatment

pt_nmds <-
  d_nmds |>
  filter(phase == '01_Pre') |>
  dplyr::select(-phase) |>
  dplyr::rename(pmds1 = MDS1, pmds2 = MDS2)

dd_nmds <-
  d_nmds |>
  filter(phase != '01_Pre') |>
  dplyr::left_join(pt_nmds) |>
  mutate(dMDS1 = MDS1 - pmds1,
         dMDS2 = MDS2 - pmds2)
  
dd_nmds |>
  ggplot(aes(x=dMDS1, y=dMDS2, color = phase, shape=phase)) +
  geom_path(aes(group = plot), color = 'grey', linewidth = 1) +
  geom_point(size=3) +
  # stat_ellipse(aes(group = site), level = 0.9) +
  facet_grid(site~PlotTreatmentStatus) +
  theme_bw() +
  geom_hline(yintercept = 0, lty=2) +
  geom_vline(xintercept = 0, lty=2) +
  stat_ellipse() +
  ggtitle()



