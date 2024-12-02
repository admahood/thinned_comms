# exploratory analysis

source("R/a_tc_data_prep.R")
library(vegan)
library(ggrepel)
library(broom)

# gllvm?

#functions==============================

nms_prep <- function(comm, meta, pa = FALSE) {
  out <- list()
  
  if(pa) nms<- metaMDS(comm %>% decostand(method = "pa")) else nms <- metaMDS(comm)
  out$nms <- nms
  out$stressplot <- stressplot(nms)
  
  out$site_scores <- as.data.frame(vegan::scores(nms)$sites) %>%
    as_tibble(rownames = "new_visit_code")%>%
    left_join(meta) %>%
    dplyr::arrange(Phase) %>%
    dplyr::filter(Phase != "Post2")
  
  ef <- envfit(nms, comm, na.rm = T, permutations = 999)
  sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
  
  out$species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
    tibble::rownames_to_column("species")  %>%
    filter(p <= 0.001)
  return(out)
}



# wrangling
by_grp <- cover_plot %>%
  filter(!CodeType %in% non_plant_codes) %>%
  group_by(CodeType,new_visit_code) %>%
  summarise(cover=sum(cover, na.rm=T)) %>%
  ungroup() %>%
  tidyr::separate(new_visit_code, c("Plot", "phase"), "\\.", remove = FALSE) %>%
  tidyr::separate(Plot, into= c("site", "plot", "trt", "id"), sep = "\\-", remove = FALSE) %>%
  mutate(phase = ifelse(phase == "PreNT", "Pre", phase),
         phase = ifelse(phase == "Post2" & site != "PHA", "Post3", phase),
         phase = ifelse(phase == "Pre", "0Pre", phase),
         site = ifelse(site %in% c("P1", "P2"), "PHA", site),
         trt = ifelse(str_detect(new_visit_code, "C"), "Control", "Treatment")) 

# functional group =============
by_grp %>%
  ggplot(aes(x=phase, y=cover, fill = CodeType)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(trt~site) +
  geom_vline(lty = 2, col = "red", xintercept = 1.5, lwd=1)

# are all the plots monumented?
# fuel model as a response variable, comm as a predictor
# larimer conservation district field dsy, noco fire shed 5m talk, Boulder Valley, longmont field day as well

# nmds ================

nms <- nms_prep(comm, meta)
nms_pa <- nms_prep(comm, meta, pa = TRUE)

# lines connecting plots
# native vs non-native
# sapling/seedling counts with 3m quadrats

p_a <- ggplot(nms$site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed() +
  geom_path(aes(group = PlotCode), color = "grey30") +
  scale_color_brewer(palette = "Set1") +
  geom_point(size=4, stroke=1,
             aes(shape = Phase,
                 color = Phase,
                 group = PlotCode)) +
  scale_shape_manual(values = c(0,1,2,3))+
  # geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
  #              aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  # geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  # ggnewscale::new_scale_color()+
  stat_ellipse(aes(group = paste(TreatmentUnit %>% str_sub(1,3), Phase), 
                   color = Phase), show.legend = F) +
  theme_classic() +
  # scale_color_manual(name = "CRP\nYear",values = c("chocolate4", "turquoise3")) +
  theme(panel.background = element_rect(fill="transparent", color = "black"),
        # legend.position = c(0,0),
        # legend.justification = c(0,0),
        legend.background = element_rect(fill=NA)) +
  facet_grid(PlotTreatmentStatus ~ TreatmentUnit %>% str_sub(1,3))+
  ggtitle("Abundance-Based");p_a

p_o <- ggplot(nms_pa$site_scores, aes(x=NMDS1, y=NMDS2)) +
  coord_fixed()+
  geom_path(aes(group = PlotCode), color = "grey30") +
  scale_color_brewer(palette = "Set1") +
  geom_point(size=4, stroke=1,
             aes(shape = Phase,
                 color = Phase,
                 group = PlotCode)) +
  scale_shape_manual(values = c(0,1,2,3))+
  # geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
  #              aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  # geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  # ggnewscale::new_scale_color()+
  stat_ellipse(aes(group = paste(TreatmentUnit %>% str_sub(1,3), Phase), 
                   color = Phase), show.legend = F) +
  theme_classic() +
  # scale_color_manual(name = "CRP\nYear",values = c("chocolate4", "turquoise3")) +
  theme(panel.background = element_rect(fill="transparent", color = "black"),
        # legend.position = c(0,0),
        # legend.justification = c(0,0),
        legend.background = element_rect(fill=NA)) +
  facet_grid(PlotTreatmentStatus ~ TreatmentUnit %>% str_sub(1,3))+
  ggtitle("Occurrence-Based");p_o

ggsave(filename = "out/nms_occurrence.png", plot = p_o, width =10, height = 4, bg="white")
ggsave(filename = "out/nms_abundance.png", plot = p_a, width =8, height = 4, bg="white")

ggpubr::ggarrange(p_a, p_o, nrow=1, ncol =2, widths = c(1,1)) %>%
  ggsave(filename = "out/nms_panelled.png", width=17, height=8, bg = "white")

# diversity ======================================

div <- data.frame(new_visit_code = rownames(comm)) %>%
  left_join(meta %>% unique())%>%
  mutate(shannon = comm %>% vegan::diversity("shannon"),
         simpson = comm %>% vegan::diversity("shannon"),
         richness = comm %>% vegan::specnumber(),
         site = TreatmentUnit %>% str_sub(1,3)) %>%
  dplyr::filter(Phase != "Post2")

# richmod <- glm(richness ~ PlotTreatmentStatus*Phase + site, data = div, 
#                family = "quasipoisson")
# performance::check_model(richmod)

richmod <- MASS::glm.nb(richness ~ PlotTreatmentStatus*Phase + site, data = div)
summary(richmod)
performance::check_model(richmod)

shanmod<- lm(shannon~ PlotTreatmentStatus*Phase + site, data = div , 
              family = Gamma())
summary(shanmod)
performance::check_model(shanmod)

simpmod<- lm(simpson~ PlotTreatmentStatus*Phase + site, data = div)
summary(simpmod)
performance::check_model(simpmod)
simpmod<- glm(simpson~ PlotTreatmentStatus:Phase, data = div, 
              family = Gamma())
summary(simpmod)

bind_rows(tidy(richmod) %>% mutate(variable = "richness"),
          tidy(shanmod) %>% mutate(variable = "shannon"), 
          tidy(simpmod) %>% mutate(variable = "simpson")) %>%
  dplyr::select(-std.error, -statistic) %>%
  pivot_wider(names_from="variable", values_from = c("estimate", "p.value"))


pdiv <- div %>%
  pivot_longer(cols = c("shannon", "simpson", "richness")) %>%
  arrange(name,Phase) %>%
  ggplot(aes(x=Phase, y = value, fill = PlotTreatmentStatus)) +
  geom_boxplot(position = "dodge") +
  facet_wrap(site~name, scales = "free_y") 

ggsave(plot = pdiv, filename = "out/diversity_indexes.png", bg = "white",
       height = 7, width =10)

vegan::nestedbetasor(comm)
vegan::nestedbetasor(comm[which(div$Phase == "Pre")])
vegan::nestedbetasor(comm[which(div$Phase == "Post")])
vegan::nestedbetasor(comm[which(div$Phase == "Post3")])

# need to make a look up table to separate control vs treatment
nestedchecker(comm)
nestedtemp(comm) %>% plot(kind = "incid")
vegan::nestedtemp(comm[which(div$Phase == "Pre")]) %>% plot(kind = "incid")
vegan::nestedtemp(comm[which(div$Phase == "Post")]) %>% plot(kind = "incid")
vegan::nestedtemp(comm[which(div$Phase == "Post3")]) %>% plot(kind = "incid")
oecosimu(comm, nestedchecker, "quasiswap")
oecosimu(comm[which(div$Phase == "Post3")], nestedchecker, "quasiswap")


