nms_prep <- function(comm, meta, pa = FALSE) {
  out <- list()
  
  if(pa) nms<- metaMDS(comm %>% decostand(method = "pa")) else nms <- metaMDS(comm)
  out$nms <- nms
  out$stressplot <- stressplot(nms)
  
  out$site_scores <- as.data.frame(vegan::scores(nms)$sites) %>%
    as_tibble(rownames = "new_visit_code")%>%
    left_join(meta)
  
  ef <- envfit(nms, comm, na.rm = T, permutations = 999)
  sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
  
  out$species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
    tibble::rownames_to_column("species")  %>%
    filter(p <= 0.001)
  return(out)
}

template <- cover_plot |>
  pivot_wider(names_from = SpeciesCode, values_from = cover, values_fill = 0) |>
  dplyr::select(new_visit_code) |>
  unique()


species_plot_invasion <- cover_plot |>
  filter(!CodeType %in% non_plant_codes) |>
  left_join(sp_list |> mutate(SpeciesCode = str_to_lower(SpeciesCode))) |>
  filter(!is.na(prevalence),
         NativityL48 == "exotic") |>
  left_join(plot_visits) |>
  tibble::rownames_to_column('new_visit_code')

exotics_only <- species_plot_invasion |>
  dplyr::select(new_visit_code, cover, SpeciesCode) |>
  pivot_wider(names_from = SpeciesCode, values_from = cover, values_fill = 0)


eo <- left_join(template, exotics_only) |>
  tidyr::replace_na(list(poco =0,
                         popr =0,   taof =0,   lida =0,   brte =0,   trdu  =0,  
                         veth   =0,canu4 =0,   lase  =0, phpr3  =0, ciar4 =0,  
                         aggi2 =0,  brin2  =0,  civu =0,  juco =0,  livu2 =0))|> 
  tibble::column_to_rownames('new_visit_code')



species_plot_invasion |>
  group_by(trt_u_adj, phase_adj, PlotTreatmentStatus) |>
  summarise(spp = paste(unique(SpeciesCode), collapse = " ")) |>
  ungroup() |>
  print(n=33)

meta <- plot_visits |>
  dplyr::select(new_visit_code, trt_u_adj, PlotTreatmentStatus, phase_adj, TreatmentMethod1, PlotCode) |>
  unique() 
library(vegan)

xx <- nms_prep(eo, meta = meta, pa=T)

ggplot(xx$site_scores |> filter(phase_adj %in% c('post04-5',"post10-11"), trt_u_adj != "PHA-1-2"), aes(x=NMDS1, y=NMDS2)) +
  # coord_fixed() +
  # geom_path(aes(group = PlotCode), color = "grey30") +
  scale_color_brewer(palette = "Set1") +
  geom_point(size=4, stroke=1,
             aes(shape = PlotTreatmentStatus,
                 color = paste(trt_u_adj, phase_adj))) +
  # scale_shape_manual(values = c(0,1,2,3))+
  # geom_segment(data = species,x=0,y=0, color = "grey",arrow = arrow(),
  #              aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  # geom_text_repel(data = species,size=4, aes(label = species), color = "grey40") +
  # ggnewscale::new_scale_color()+
  stat_ellipse(aes(group = paste(trt_u_adj, phase_adj), 
                    color = paste(trt_u_adj, phase_adj)), show.legend = F) +
  theme_classic() +
    # scale_color_manual(name = "CRP\nYear",values = c("chocolate4", "turquoise3")) +
  theme(panel.background = element_rect(fill="transparent", color = "black"),
        # legend.position = c(0,0),
        # legend.justification = c(0,0),
        legend.background = element_rect(fill=NA)) +
  facet_wrap(~trt_u_adj) +
  ggtitle("Abundance-Based")
