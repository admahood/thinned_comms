# Hmsc burned and unburned with b_tectorum as a predictor

# setup ========================================================================
source("R/a_tc_data_prep.R")

library(knitr)
library(ape)
library(MASS)
library(fields)
library(Hmsc)
library(parallel)
library(ggthemes)
library(ggtext)
library(corrplot)
library(ggpubr)
set.seed(1)
theme_set(theme_classic())


# data wrangling ===============================================================
# veg community
Y <- # Hmsc burned and unburned with b_tectorum as a predictor

# setup ========================================================================
source("R/a_data_prep.R")

library(knitr)
library(ape)
library(MASS)
library(fields)
library(Hmsc)
library(parallel)
library(ggthemes)
library(ggtext)
library(corrplot)
library(ggpubr)
set.seed(1)
theme_set(theme_classic())


# data wrangling ===============================================================
# veg community
Y <- comm; dim(Y)
Y[Y>0] <-1;dim(Y)
# Y<- Y[,colSums(Y)>1];dim(Y) # removing species that only occur in one location
colnames(Y)
# colnames(C) <- str_replace_all(colnames(C), " ", "_")
prevalence<- colSums(Y) %>%
  as_tibble(rownames = "Species") %>%
  dplyr::rename(prevalence = value)

coords <-  sites %>%
  arrange(plot) %>%
  st_coordinates %>%
  as.data.frame()  %>%
  dplyr::rename("x-coordinate" = X, "y-coordinate" = Y)

#env data

XData <- sites_w_grazing %>%
  arrange(plot) %>%
  dplyr::select(watershed, elevation_m, folded_aspect,plot,
                grazing_intensity, burned, B_tectorum) %>%
  mutate(watershed = as.factor(watershed),
         burned=as.factor(burned))%>%
  tibble::column_to_rownames("plot") %>%
  as.data.frame %>%
  mutate("y-coordinate" = coords[,1],
         "x-coordinate" = coords[,2])


#species traits
traits <- as.data.frame(veg_traits) %>%
  filter(species %in% colnames(Y))%>%
  dplyr::select(-fg) %>%
  transmute_all((as.factor)) %>%
  tibble::column_to_rownames("species")


tr_form <- ~origin+duration+cots

XFormula <- ~elevation_m+
  folded_aspect+
  grazing_intensity+
  burned+
  B_tectorum

studyDesign <- data.frame(watershed = as.factor(XData$watershed))
studyDesign <- data.frame(sample = as.factor(1:nrow(XData)))

rL = HmscRandomLevel(units = studyDesign$watershed)
rLs = HmscRandomLevel(sData = coords%>% as.matrix()) %>%
  setPriors(nfMin=1, nfMax=1)

# fitting the model ============================================================
mod = Hmsc(Y = Y, XData = XData, XFormula = XFormula,distr="probit",
         TrData = traits,TrFormula = tr_form,
         studyDesign = studyDesign, ranLevels = list(sample = rLs))

nChains = 4
test.run = FALSE
if (test.run){
  #with this option, the vignette evaluates in ca. 1 minute in adam's laptop
  thin = 1
  samples = 100
  transient = ceiling(thin*samples*.5)
}else{
  # with a spatial random effect, evaluates in ---
  # looks like a compute-optimized aws instance is called for, very little ram usage
  thin = 100
  samples = 1000
  transient = ceiling(thin*samples*.5)
}
t0 <- Sys.time()
hmsc_file <- "data/hmsc/hmsc_probit_all_every_sp_w_bromus.Rda"
dir.create("data/hmsc")
if(!file.exists(hmsc_file)){
  m = sampleMcmc(mod, thin = thin, 
                 samples = samples, 
                 transient = transient,
                 adaptNf = rep(ceiling(0.4*samples*thin),1),
                 nChains = nChains, 
                 nParallel = nChains)
  print(Sys.time()-t0)
  save(m, file=hmsc_file)
  system(paste("aws s3 cp",
               hmsc_file,
               file.path("s3://earthlab-amahood/lyb", hmsc_file)))
}else{load(hmsc_file)}

mpost <- convertToCodaObject(m)
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
VP <- computeVariancePartitioning(m)
# mcmc convergence ===================
psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf%>%
  as_tibble() %>% dplyr::rename(psrf_v = `Point est.`)


ess.beta <- effectiveSize(mpost$Beta) %>%
  as_tibble() %>% dplyr::rename(ess_beta = value)

ess.v <- effectiveSize(mpost$V)%>%
  as_tibble() %>% dplyr::rename(ess_v = value)
psrf.beta <- gelman.diag(mpost$Beta, multivariate=FALSE)$psrf%>%
  as_tibble() %>% dplyr::rename(psrf_beta = `Point est.`)

# ess.gamma = effectiveSize(mpost$Gamma)%>%
#   as_tibble() %>% rename(ess_gamma = value)
# psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf%>%
#   as_tibble() %>% rename(psrf_gamma = `Point est.`)
# 
# sppairs = matrix(sample(x = 1:ns^2, size = 100))
# tmp = mpost$Omega[[1]]
# for (chain in 1:length(tmp)){
#   tmp[[chain]] = tmp[[chain]][,sppairs]
# }
# 
# ess.omega = effectiveSize(tmp)%>%
#   as_tibble() %>% rename(ess_omega = value)
# psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf%>%
#   as_tibble() %>% rename(psrf_omega = `Point est.`)

diag_all <- ggarrange(ggplot(ess.beta, aes(x=ess_beta)) + 
                        geom_histogram()+
                        xlab("Effective Sample Size"),
                      ggplot(psrf.beta, aes(x=psrf_beta)) +
                        geom_histogram()+
                        xlab("Gelman Diagnostic"),
          align = "v") +ggtitle("All Plots")+
  ggsave("images/geldman_allplots_wbrom.png", width = 5.5, height=3.5, bg="white")
save(diag_all, file="data/diag_all.Rda")
# fit stuff ===========
MF$TjurR2 %>% mean(na.rm=T)
# explanatory power

ggarrange(
ggplot(as.data.frame(MF),aes(x=(RMSE))) + geom_histogram(),
ggplot(as.data.frame(MF),aes(x=(TjurR2))) + geom_histogram(),
ggplot(as.data.frame(MF),aes(x=(AUC))) + geom_histogram())

# predictive power
# 
# partition = createPartition(m, nfolds = 10) # takes 2.5 hrs with 10 folds
# cv_file <- "data/hmsc/cv_binomial_all_every_sp.Rda"
# # cv_file <- "data/hmsc/cv_binomial_all.Rda"
# 
# if(!file.exists(cv_file)){
#   t0<- Sys.time()
#   print(t0)
#   predYCV = computePredictedValues(m, partition = partition, nParallel = 2)
#   save(predYCV, file = cv_file)
#   print(Sys.time()-t0)
# }else{load(cv_file)}
# # Note that in computePredicted values it is also possible to use the nParallel option
# # Below we construct a plot that compares explanatory power (MF) to predictive power (MFCV)
# # As expected, the explanatory power is higher than the predictive power
# 
# MFCV = evaluateModelFit(hM=m, predY=predYCV)
# plot(MF$AUC, MFCV$AUC)
# abline(0,1)


# plotting variance partitioning ===============================================

mf_df <- data.frame(Species = colnames(m$Y),
                    R2 = MF$TjurR2,
                    AUC = MF$AUC,
                    RMSE = MF$RMSE) %>%
  left_join(prevalence)
mean(mf_df%>% filter(prevalence>7) %>% pull(R2), na.rm=T)
ggplot(mf_df, aes(x=prevalence, y=RMSE)) +
  geom_point()
# VP$R2T$Beta
# VP$R2T$Y
cols<-vector()
cols[1] <- "#FED439"
cols[2] <- RColorBrewer::brewer.pal(12, "Paired")[c(12)]
cols[3:5] <-  RColorBrewer::brewer.pal(3, "Blues")
cols[6]<- "grey"

sbquants <- summary(mpost$Beta)$quantiles %>%
  as_tibble(rownames = "variable") %>% 
  mutate(sign = `2.5%` * `97.5%`) %>%
  filter(sign>0) %>%
  separate(variable,
           into = c("variable", "species"),
           sep = ",") %>%
  mutate(variable = str_sub(variable, 3,nchar(variable)-5),
         species = str_sub(species, 2,nchar(species)-6) %>% trimws) %>%
  filter(variable!= "(Intercept)") %>%
  dplyr::select(variable,species,`2.5%`,`50%`,`97.5%`) %>%
  arrange(variable)


vp_df <- VP$vals%>%
  as_tibble(rownames = "variable") %>%
  pivot_longer(cols=names(.)[2:ncol(.)], 
               names_to = "Species", 
               values_to = "value") %>%
  mutate(origin = lut_all_fg[Species] %>%
           str_sub(1,1),
         fg = lut_all_fg[Species] %>%
           str_sub(2,3))  %>%
  left_join(prevalence) %>%
  na.omit()

vp_summary <- vp_df %>%
  group_by(origin, variable) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  pivot_wider(id_cols = "variable", 
              names_from = "origin", 
              values_from = "value")

vp_order <- vp_df %>%
  filter(variable == "Random: sample") %>%
  arrange(origin,fg, prevalence) %>%
  mutate(Species_f = factor(Species, levels = .$Species)) %>%
  dplyr::select(Species, Species_f, origin, fg) 

# 
# vp_order <- vp_df %>% filter(variable == "Random: sample") %>%
#   filter(origin=="I") %>%
#   left_join(prevalence) %>%
#   arrange(prevalence, origin) %>%
#   mutate(Species_f = factor(Species, levels = .$Species)) %>%
#   dplyr::select(Species, Species_f, origin) %>%
#   rbind(vp_order_n)# %>%
#   #left_join(mf_df)



left_join(vp_df, vp_order) %>% 
  mutate(variable = factor(variable, 
                           levels = c("elevation_m",
                                      "folded_aspect",
                                      "grazing_intensity", 
                                      "burned", 
                                      "B_tectorum",
                                      "Random: sample"),
                           labels = c("Elevation",
                                      "Folded Aspect",
                                      "Grazing Intensity", 
                                      "Fire Occurrence", 
                                      "*B. tectorum* cover",
                                      "Random: Spatial")),
         value = value) %>%
  ggplot(aes(x=value,y=Species_f, fill = variable)) +
  geom_bar(stat="identity")+
  theme_classic() +
  geom_hline(yintercept = table(vp_order$origin)[1]+.5) +
  geom_hline(yintercept = nrow(vp_order)+.5) +
  annotate("text", x = 1.2, y=1, label="Introduced", angle=90, vjust="bottom",
           hjust="left", size=8)+
  annotate("text", x = 1.2, y=nrow(vp_order), label="Native", angle=90, vjust="top",
           hjust="right", size=8)+
  ylab("Species") +
  xlab("Proportion of Variance Explained") +
  scale_fill_manual(values = cols)+
  theme(legend.position = c(1,.315),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.background = element_rect(color="black")) +
  ggtitle("Variance Partitioning, Occurrence Model")+
  ggsave("images/variance_partitioning_allspp_w_bromus.png", height = 10.5, width = 7)


# table(veg_traits$origin)[1]

# species niches ...basically ==================================================

postBeta <- getPostEstimate(m, parName = "Beta")

means <- postBeta$mean %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",VP$groupnames)) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Mean")


supported_all <- postBeta$support %>% 
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",VP$groupnames)) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], 
               names_to = "Species", 
               values_to = "Support") %>%
  filter(Support >0.90|Support<0.1,
         env_var != "intercept") %>%
  left_join(means, by = c("env_var", "Species"))%>%
  mutate(sign = ifelse(Mean>0, "+", "-"))%>%
  mutate(env_var = factor(env_var, 
                          levels = c("elevation_m",
                                     "folded_aspect",
                                     "grazing_intensity", 
                                     "burned", 
                                     "B_tectorum"),
                          labels = c("Elevation (m)",
                                     "Folded Aspect",
                                     "Grazing Intensity", 
                                     "Burned", 
                                     "*B. tectorum* cover"))) %>%
  left_join(vp_order)  %>%
  filter(Species != "Bassia prostrata", # these are there because people planted them...
         Species != "Agropyron cristatum",
         Species != "Unknown perennial grass") %>%
  mutate(fg = lut_all_fg[Species] %>% str_sub(2,3),
         fg = c("W" = "Shrub",
                "AF" = "Annual Forb",
                "AG" = "Annual Grass",
                "PF" = "Perennial Forb",
                "C" = "Shrub",
                "PG" = "Perennial Grass")[fg])

line_df_all <- supported_all %>%
  group_by(Species) %>%
  dplyr::summarise(origin = first(origin)) %>%
  ungroup()

pal<-RColorBrewer::brewer.pal(5, "Accent")
pal <- pal[c(1,2,3,5)]
p_betas_all_plots<-ggplot(supported_all, aes(x=env_var, y = Species_f, 
                                         fill = sign)) +
  
  theme_clean()+
  geom_hline(lwd=11,aes(yintercept=Species, color=fg))+
  geom_hline(lty=3, aes(yintercept=Species), color="black")+
  geom_tile(lwd=.5, color = "black")+#, aes(alpha=Mean)) +
  # scale_color_brewer(name = "Functional Group", palette = "Accent") +
  scale_color_manual(name = "Functional Group", values = pal) +
  # scale_color_colorblind()+
  guides(alpha="none")+
  geom_hline(yintercept = table(line_df_all$origin)[1] + .5, lwd=1) +
  annotate("text", x = 6, y=1, label="Introduced", angle=90, vjust="bottom",
           hjust="left", size=6)+
  annotate("text", x = 6, y=nrow(line_df_all)/1.5, label="Native", angle=90, 
           vjust="bottom",
           hjust="right", size=6)+
  theme(legend.position = "right",#c(1,.315),
        # legend.title = element_blank(),
        # legend.justification = c(1,0),
        legend.background = element_rect(color="black")) +
  theme(axis.text.x = element_text(angle=45, vjust=1,hjust = 1),
        axis.title = element_blank()) +
  ggtitle("All Plot Locations") +
  ggsave("images/betas_binomial_all_spp.png")

save(p_betas_all_plots, supported_all, file = "data/p_beta_all_plots.Rda")


# sanity checks==========
plotBeta(m, post = postBeta, param = "Support",
         supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
plotBeta(m, post = postBeta,param = "Mean",
         supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post=postGamma)
plotGamma(m, post=postGamma, param="Support", supportLevel = 0.95)



# colors for plots =============

pal<-RColorBrewer::brewer.pal(5,"Accent")
# pal[4] <- "#FFC845"
# pal[3] <- "orange"
pal_int <- pal[1:2]
pal_nat <- pal[c(1,3:5)]

# gradient grazing intensity ==========

grazing_gradient = constructGradient(m, focalVariable = "grazing_intensity")

predY_grazing = predict(m, XData=grazing_gradient$XDataNew, 
                        studyDesign=grazing_gradient$studyDesignNew, 
                ranLevels=grazing_gradient$rLNew, expected=TRUE)

n_runs <- nChains*samples

pred_df_grazing <- do.call("rbind", predY_grazing) %>%
  as_tibble() %>%
  mutate(grazing_intensity = rep(grazing_gradient$XDataNew$grazing_intensity,
                                 n_runs)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -grazing_intensity) %>%
  group_by(Species, grazing_intensity) %>%
  dplyr::summarise(mean = median(cover),
                   upr = quantile(cover, 0.975),
                   lwr = quantile(cover, 0.025)) %>%
  ungroup() %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  filter(origin != "IPG" & origin != "IW") %>%
  arrange(desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) %>%
  filter(Species != "unknown_forb") 

pred_df_raw_grazing <- do.call("rbind", predY_grazing) %>%
  as_tibble() %>%
  mutate(grazing_intensity = rep(grazing_gradient$XDataNew$grazing_intensity, 
                                 n_runs),
         run = rep(1:n_runs,each=20)) %>%
  pivot_longer(values_to = "cover", names_to = "Species",
               -c(grazing_intensity,run))%>%
  # mutate(Species = str_replace_all(Species," ", "_")) %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  arrange(origin,desc(prevalence)) %>%
  filter(Species != "unknown_forb")%>%
  mutate(Species_f = factor(Species, levels = unique(.$Species)))

grazing_native <- pred_df_raw_grazing  %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "N") %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin]) %>%
  # mutate(Species_f = lut_sp_nice[Species_f])%>%
  ggplot(aes(x=grazing_intensity, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run, color = origin), key_glyph="rect")+
  geom_line(data = pred_df_grazing  %>%
              # mutate(Species_f = lut_sp_nice[Species_f])%>%
              filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
              filter(str_sub(origin,1,1) == "N") ,
            lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  xlab("Grazing Intensity (AUM/ha)") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  # scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_nat)+
  theme(legend.position = "right",
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.title = element_blank())

grazing_introduced <- pred_df_raw_grazing  %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "I") %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin])%>%
  ggplot(aes(x=grazing_intensity, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run, color = origin), key_glyph="rect")+
  geom_line(data = pred_df_grazing %>%
              filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
              filter(str_sub(origin,1,1) == "I") ,
            lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  xlab("Grazing Intensity (AUM/ha)") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  # scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_int)+
  theme(legend.position = "right",
        legend.justification = c(1,0),
        panel.border = element_rect(fill=NA, size=0.75),
        strip.text = element_markdown(),
        legend.title = element_blank())

# gradient burned ==================
m$XData$burned <- as.factor(m$XData$burned)
burned_gradient = constructGradient(m, focalVariable = "burned")

predY_burned = predict(m, XData=burned_gradient$XDataNew, 
                       studyDesign=burned_gradient$studyDesignNew, 
                ranLevels=burned_gradient$rLNew, expected=TRUE)


n_runs <- nChains*samples


pred_df_burned <- do.call("rbind", predY_burned) %>%
  as_tibble() %>%
  mutate(burned = rep(burned_gradient$XDataNew$burned, n_runs)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -burned) %>%
  group_by(Species, burned) %>%
  dplyr::summarise(mean = median(cover),
                   upr = quantile(cover, 0.975),
                   lwr = quantile(cover, 0.025)) %>%
  ungroup() %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  filter(origin != "IPG" & origin != "IW") %>%
  arrange(origin,desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) %>%
  filter(Species != "unknown forb") 

pred_df_raw_burned <-do.call("rbind", predY_burned) %>%
  as_tibble() %>%
  mutate(burned = rep(burned_gradient$XDataNew$burned, n_runs),
         run = rep(1:n_runs,each=2)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -c(burned,run))%>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  filter(origin != "IPG" & origin != "IW") %>%
  arrange(origin, desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) %>%
  filter(Species != "unknown forb")  %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin])

burned_introduced<-pred_df_raw_burned %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "I")  %>%
  mutate(burned = ifelse(burned == "yes", "B", "U")%>%
           factor(levels = c("U", "B"))) %>%
  ggplot(aes(x=burned, y=cover)) +
  geom_jitter(alpha = 0.03, aes(color = origin))+
  geom_boxplot(fill = "transparent", outlier.shape = NA)+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  xlab("Fire Occurrence") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_int)+
  theme(legend.position = c(1,0),
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())

burned_native<-pred_df_raw_burned %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "N")  %>%
  mutate(burned = ifelse(burned == "yes", "B", "U")%>%
           factor(levels = c("U", "B"))) %>%
  ggplot(aes(x=burned, y=cover)) +
  geom_jitter(alpha = 0.03, aes(color = origin))+
  geom_boxplot(fill = "transparent", outlier.shape = NA)+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  xlab("Fire Occurrence") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_nat)+
  theme(legend.position = c(1,0),
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())
# gradient bromus cover ==========

bromus_gradient = constructGradient(m, focalVariable = "B_tectorum", 
                             non.focalVariables=list(burned = list(1)))

predY_bromus = predict(m, XData=bromus_gradient$XDataNew, 
                       studyDesign=bromus_gradient$studyDesignNew, 
                ranLevels=bromus_gradient$rLNew, expected=TRUE)

n_runs <- nChains*samples

pred_df_bromus <- do.call("rbind", predY_bromus) %>%
  as_tibble() %>%
  mutate(B_tectorum = rep(bromus_gradient$XDataNew$B_tectorum, n_runs)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -B_tectorum) %>%
  group_by(Species, B_tectorum) %>%
  dplyr::summarise(mean = median(cover),
                   upr = quantile(cover, 0.975),
                   lwr = quantile(cover, 0.025)) %>%
  ungroup() %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  filter(origin != "IPG" & origin != "IW") %>%
  arrange(desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) %>%
  filter(Species != "unknown_forb") 

pred_df_raw_bromus <- do.call("rbind", predY_bromus) %>%
  as_tibble() %>%
  mutate(B_tectorum = rep(bromus_gradient$XDataNew$B_tectorum, n_runs),
         run = rep(1:n_runs,each=20)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -c(B_tectorum,run))%>%
  # mutate(Species = str_replace_all(Species," ", "_")) %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  arrange(origin,desc(prevalence)) %>%
  filter(Species != "unknown_forb")%>%
  mutate(Species_f = factor(Species, levels = unique(.$Species)))

bromus_introduced <- pred_df_raw_bromus  %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "I") %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin])%>%
  ggplot(aes(x=B_tectorum, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run, color = origin), key_glyph="rect")+
  geom_line(data = pred_df_bromus %>%
              filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
              filter(str_sub(origin,1,1) == "I"),
            lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  labs(x="*Bromus tectorum* cover (%)",
       y= "Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_int)+
  # scale_color_manual(values = c("#FFC845", "#007DBA"), 
  #                    labels = c("Introduced", "Native"))+
  theme(legend.position = c(1,0),
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        axis.title.x = element_markdown(),
        legend.justification = c(1,0),
        legend.title = element_blank())

bromus_native <- pred_df_raw_bromus  %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "N") %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin])%>%
  ggplot(aes(x=B_tectorum, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run, color = origin), key_glyph="rect")+
  geom_line(data = pred_df_bromus %>%
              filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
              filter(str_sub(origin,1,1) == "N"),
            lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  labs(x="*Bromus tectorum* cover (%)",
       y= "Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_nat)+
  # scale_color_manual(values = c("#FFC845", "#007DBA"), 
  #                    labels = c("Introduced", "Native"))+
  theme(legend.position = c(1,0),
        panel.border = element_rect(fill=NA, size=0.75),
        strip.text = element_markdown(),
        axis.title.x = element_markdown(),
        legend.justification = c(1,0),
        legend.title = element_blank())


# all together =====

nnn<-ggarrange(grazing_native, burned_native,bromus_native,  
               nrow = 3, common.legend = T, labels = c("(a)","(b)","(c)"),
               label.x = 0.01)+
  ggsave("images/probit_preds_nat_3pan.png",
         height = 7, width=12)

iii<-ggarrange(grazing_introduced, burned_introduced,bromus_introduced,  
               nrow = 3, common.legend = T, labels = c("(a)","(b)","(c)"),
               label.x = 0.01)+
  ggsave("images/probit_preds_int_3pan.png", 
         height = 7, width=12)

Y[Y>0] <-1

# Y<- Y[,colSums(Y)>1] # removing species that only occur in one location

# colnames(C) <- str_replace_all(colnames(C), " ", "_")
prevalence<- colSums(Y) %>%
  as_tibble(rownames = "Species") %>%
  dplyr::rename(prevalence = value)

coords <-  sites %>%
  arrange(plot) %>%
  st_coordinates %>%
  as.data.frame()  %>%
  dplyr::rename("x-coordinate" = X, "y-coordinate" = Y)

#env data

XData <- sites_w_grazing %>%
  arrange(plot) %>%
  dplyr::select(watershed, elevation_m, folded_aspect,plot,
                grazing_intensity, burned, B_tectorum) %>%
  mutate(watershed = as.factor(watershed),
         burned=as.factor(burned))%>%
  tibble::column_to_rownames("plot") %>%
  as.data.frame %>%
  mutate("y-coordinate" = coords[,1],
         "x-coordinate" = coords[,2])


#species traits
traits <- as.data.frame(veg_traits) %>%
  filter(species %in% colnames(Y))%>%
  dplyr::select(-fg) %>%
  transmute_all((as.factor)) %>%
  tibble::column_to_rownames("species")


tr_form <- ~origin+duration+cots

XFormula <- ~elevation_m+
  folded_aspect+
  grazing_intensity+
  burned+
  B_tectorum

studyDesign <- data.frame(watershed = as.factor(XData$watershed))
studyDesign <- data.frame(sample = as.factor(1:nrow(XData)))

rL = HmscRandomLevel(units = studyDesign$watershed)
rLs = HmscRandomLevel(sData = coords%>% as.matrix()) %>%
  setPriors(nfMin=1, nfMax=1)

# fitting the model ============================================================
mod = Hmsc(Y = Y, XData = XData, XFormula = XFormula,distr="probit",
           TrData = traits,TrFormula = tr_form,
           studyDesign = studyDesign, ranLevels = list(sample = rLs))

nChains = 4
test.run = FALSE
if (test.run){
  #with this option, the vignette evaluates in ca. 1 minute in adam's laptop
  thin = 1
  samples = 100
  transient = ceiling(thin*samples*.5)
}else{
  # with a spatial random effect, evaluates in ---
  # looks like a compute-optimized aws instance is called for, very little ram usage
  thin = 100
  samples = 1000
  transient = ceiling(thin*samples*.5)
}
t0 <- Sys.time()
hmsc_file <- "data/hmsc/hmsc_probit_all_every_sp_w_bromus.Rda"
dir.create("data/hmsc")
if(!file.exists(hmsc_file)){
  m = sampleMcmc(mod, thin = thin, 
                 samples = samples, 
                 transient = transient,
                 adaptNf = rep(ceiling(0.4*samples*thin),1),
                 nChains = nChains, 
                 nParallel = nChains)
  print(Sys.time()-t0)
  save(m, file=hmsc_file)
  system(paste("aws s3 cp",
               hmsc_file,
               file.path("s3://earthlab-amahood/lyb", hmsc_file)))
}else{load(hmsc_file)}

mpost <- convertToCodaObject(m)
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds)
VP <- computeVariancePartitioning(m)
# mcmc convergence ===================
psrf.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf%>%
  as_tibble() %>% dplyr::rename(psrf_v = `Point est.`)


ess.beta <- effectiveSize(mpost$Beta) %>%
  as_tibble() %>% dplyr::rename(ess_beta = value)

ess.v <- effectiveSize(mpost$V)%>%
  as_tibble() %>% dplyr::rename(ess_v = value)
psrf.beta <- gelman.diag(mpost$Beta, multivariate=FALSE)$psrf%>%
  as_tibble() %>% dplyr::rename(psrf_beta = `Point est.`)

# ess.gamma = effectiveSize(mpost$Gamma)%>%
#   as_tibble() %>% rename(ess_gamma = value)
# psrf.gamma = gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf%>%
#   as_tibble() %>% rename(psrf_gamma = `Point est.`)
# 
# sppairs = matrix(sample(x = 1:ns^2, size = 100))
# tmp = mpost$Omega[[1]]
# for (chain in 1:length(tmp)){
#   tmp[[chain]] = tmp[[chain]][,sppairs]
# }
# 
# ess.omega = effectiveSize(tmp)%>%
#   as_tibble() %>% rename(ess_omega = value)
# psrf.omega = gelman.diag(tmp, multivariate=FALSE)$psrf%>%
#   as_tibble() %>% rename(psrf_omega = `Point est.`)

diag_all <- ggarrange(ggplot(ess.beta, aes(x=ess_beta)) + 
                        geom_histogram()+
                        xlab("Effective Sample Size"),
                      ggplot(psrf.beta, aes(x=psrf_beta)) +
                        geom_histogram()+
                        xlab("Gelman Diagnostic"),
                      align = "v") +ggtitle("All Plots")+
  ggsave("images/geldman_allplots_wbrom.png", width = 5.5, height=3.5, bg="white")
save(diag_all, file="data/diag_all.Rda")
# fit stuff ===========
MF$TjurR2 %>% mean(na.rm=T)
# explanatory power

ggarrange(
  ggplot(as.data.frame(MF),aes(x=(RMSE))) + geom_histogram(),
  ggplot(as.data.frame(MF),aes(x=(TjurR2))) + geom_histogram(),
  ggplot(as.data.frame(MF),aes(x=(AUC))) + geom_histogram())

# predictive power
# 
# partition = createPartition(m, nfolds = 10) # takes 2.5 hrs with 10 folds
# cv_file <- "data/hmsc/cv_binomial_all_every_sp.Rda"
# # cv_file <- "data/hmsc/cv_binomial_all.Rda"
# 
# if(!file.exists(cv_file)){
#   t0<- Sys.time()
#   print(t0)
#   predYCV = computePredictedValues(m, partition = partition, nParallel = 2)
#   save(predYCV, file = cv_file)
#   print(Sys.time()-t0)
# }else{load(cv_file)}
# # Note that in computePredicted values it is also possible to use the nParallel option
# # Below we construct a plot that compares explanatory power (MF) to predictive power (MFCV)
# # As expected, the explanatory power is higher than the predictive power
# 
# MFCV = evaluateModelFit(hM=m, predY=predYCV)
# plot(MF$AUC, MFCV$AUC)
# abline(0,1)


# plotting variance partitioning ===============================================

mf_df <- data.frame(Species = colnames(m$Y),
                    R2 = MF$TjurR2,
                    AUC = MF$AUC,
                    RMSE = MF$RMSE) %>%
  left_join(prevalence)
mean(mf_df%>% filter(prevalence>7) %>% pull(R2), na.rm=T)
ggplot(mf_df, aes(x=prevalence, y=RMSE)) +
  geom_point()
# VP$R2T$Beta
# VP$R2T$Y
cols<-vector()
cols[1] <- "#FED439"
cols[2] <- RColorBrewer::brewer.pal(12, "Paired")[c(12)]
cols[3:5] <-  RColorBrewer::brewer.pal(3, "Blues")
cols[6]<- "grey"

sbquants <- summary(mpost$Beta)$quantiles %>%
  as_tibble(rownames = "variable") %>% 
  mutate(sign = `2.5%` * `97.5%`) %>%
  filter(sign>0) %>%
  separate(variable,
           into = c("variable", "species"),
           sep = ",") %>%
  mutate(variable = str_sub(variable, 3,nchar(variable)-5),
         species = str_sub(species, 2,nchar(species)-6) %>% trimws) %>%
  filter(variable!= "(Intercept)") %>%
  dplyr::select(variable,species,`2.5%`,`50%`,`97.5%`) %>%
  arrange(variable)


vp_df <- VP$vals%>%
  as_tibble(rownames = "variable") %>%
  pivot_longer(cols=names(.)[2:ncol(.)], 
               names_to = "Species", 
               values_to = "value") %>%
  mutate(origin = lut_all_fg[Species] %>%
           str_sub(1,1),
         fg = lut_all_fg[Species] %>%
           str_sub(2,3))  %>%
  left_join(prevalence) %>%
  na.omit()

vp_summary <- vp_df %>%
  group_by(origin, variable) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>%
  pivot_wider(id_cols = "variable", 
              names_from = "origin", 
              values_from = "value")

vp_order <- vp_df %>%
  filter(variable == "Random: sample") %>%
  arrange(origin,fg, prevalence) %>%
  mutate(Species_f = factor(Species, levels = .$Species)) %>%
  dplyr::select(Species, Species_f, origin, fg) 

# 
# vp_order <- vp_df %>% filter(variable == "Random: sample") %>%
#   filter(origin=="I") %>%
#   left_join(prevalence) %>%
#   arrange(prevalence, origin) %>%
#   mutate(Species_f = factor(Species, levels = .$Species)) %>%
#   dplyr::select(Species, Species_f, origin) %>%
#   rbind(vp_order_n)# %>%
#   #left_join(mf_df)



left_join(vp_df, vp_order) %>% 
  mutate(variable = factor(variable, 
                           levels = c("elevation_m",
                                      "folded_aspect",
                                      "grazing_intensity", 
                                      "burned", 
                                      "B_tectorum",
                                      "Random: sample"),
                           labels = c("Elevation",
                                      "Folded Aspect",
                                      "Grazing Intensity", 
                                      "Fire Occurrence", 
                                      "*B. tectorum* cover",
                                      "Random: Spatial")),
         value = value) %>%
  ggplot(aes(x=value,y=Species_f, fill = variable)) +
  geom_bar(stat="identity")+
  theme_classic() +
  geom_hline(yintercept = table(vp_order$origin)[1]+.5) +
  geom_hline(yintercept = nrow(vp_order)+.5) +
  annotate("text", x = 1.2, y=1, label="Introduced", angle=90, vjust="bottom",
           hjust="left", size=8)+
  annotate("text", x = 1.2, y=nrow(vp_order), label="Native", angle=90, vjust="top",
           hjust="right", size=8)+
  ylab("Species") +
  xlab("Proportion of Variance Explained") +
  scale_fill_manual(values = cols)+
  theme(legend.position = c(1,.315),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.justification = c(1,0),
        legend.background = element_rect(color="black")) +
  ggtitle("Variance Partitioning, Occurrence Model")+
  ggsave("images/variance_partitioning_allspp_w_bromus.png", height = 10.5, width = 7)


# table(veg_traits$origin)[1]

# species niches ...basically ==================================================

postBeta <- getPostEstimate(m, parName = "Beta")

means <- postBeta$mean %>%
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",VP$groupnames)) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Mean")


supported_all <- postBeta$support %>% 
  as_tibble() %>%
  rowid_to_column("env_var") %>%
  mutate(env_var = c("intercept",VP$groupnames)) %>%
  pivot_longer(cols=names(.)[2:ncol(.)], 
               names_to = "Species", 
               values_to = "Support") %>%
  filter(Support >0.90|Support<0.1,
         env_var != "intercept") %>%
  left_join(means, by = c("env_var", "Species"))%>%
  mutate(sign = ifelse(Mean>0, "+", "-"))%>%
  mutate(env_var = factor(env_var, 
                          levels = c("elevation_m",
                                     "folded_aspect",
                                     "grazing_intensity", 
                                     "burned", 
                                     "B_tectorum"),
                          labels = c("Elevation (m)",
                                     "Folded Aspect",
                                     "Grazing Intensity", 
                                     "Burned", 
                                     "*B. tectorum* cover"))) %>%
  left_join(vp_order)  %>%
  filter(Species != "Bassia prostrata", # these are there because people planted them...
         Species != "Agropyron cristatum",
         Species != "Unknown perennial grass") %>%
  mutate(fg = lut_all_fg[Species] %>% str_sub(2,3),
         fg = c("W" = "Shrub",
                "AF" = "Annual Forb",
                "AG" = "Annual Grass",
                "PF" = "Perennial Forb",
                "C" = "Shrub",
                "PG" = "Perennial Grass")[fg])

line_df_all <- supported_all %>%
  group_by(Species) %>%
  dplyr::summarise(origin = first(origin)) %>%
  ungroup()

pal<-RColorBrewer::brewer.pal(5, "Accent")
pal <- pal[c(1,2,3,5)]
p_betas_all_plots<-ggplot(supported_all, aes(x=env_var, y = Species_f, 
                                             fill = sign)) +
  
  theme_clean()+
  geom_hline(lwd=11,aes(yintercept=Species, color=fg))+
  geom_hline(lty=3, aes(yintercept=Species), color="black")+
  geom_tile(lwd=.5, color = "black")+#, aes(alpha=Mean)) +
  # scale_color_brewer(name = "Functional Group", palette = "Accent") +
  scale_color_manual(name = "Functional Group", values = pal) +
  # scale_color_colorblind()+
  guides(alpha="none")+
  geom_hline(yintercept = table(line_df_all$origin)[1] + .5, lwd=1) +
  annotate("text", x = 6, y=1, label="Introduced", angle=90, vjust="bottom",
           hjust="left", size=6)+
  annotate("text", x = 6, y=nrow(line_df_all)/1.5, label="Native", angle=90, 
           vjust="bottom",
           hjust="right", size=6)+
  theme(legend.position = "right",#c(1,.315),
        # legend.title = element_blank(),
        # legend.justification = c(1,0),
        legend.background = element_rect(color="black")) +
  theme(axis.text.x = element_text(angle=45, vjust=1,hjust = 1),
        axis.title = element_blank()) +
  ggtitle("All Plot Locations") +
  ggsave("images/betas_binomial_all_spp.png")

save(p_betas_all_plots, supported_all, file = "data/p_beta_all_plots.Rda")


# sanity checks==========
plotBeta(m, post = postBeta, param = "Support",
         supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
plotBeta(m, post = postBeta,param = "Mean",
         supportLevel = 0.95, split=.4, spNamesNumbers = c(T,F))
postGamma = getPostEstimate(m, parName = "Gamma")
plotGamma(m, post=postGamma)
plotGamma(m, post=postGamma, param="Support", supportLevel = 0.95)



# colors for plots =============

pal<-RColorBrewer::brewer.pal(5,"Accent")
# pal[4] <- "#FFC845"
# pal[3] <- "orange"
pal_int <- pal[1:2]
pal_nat <- pal[c(1,3:5)]

# gradient grazing intensity ==========

grazing_gradient = constructGradient(m, focalVariable = "grazing_intensity")

predY_grazing = predict(m, XData=grazing_gradient$XDataNew, 
                        studyDesign=grazing_gradient$studyDesignNew, 
                        ranLevels=grazing_gradient$rLNew, expected=TRUE)

n_runs <- nChains*samples

pred_df_grazing <- do.call("rbind", predY_grazing) %>%
  as_tibble() %>%
  mutate(grazing_intensity = rep(grazing_gradient$XDataNew$grazing_intensity,
                                 n_runs)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -grazing_intensity) %>%
  group_by(Species, grazing_intensity) %>%
  dplyr::summarise(mean = median(cover),
                   upr = quantile(cover, 0.975),
                   lwr = quantile(cover, 0.025)) %>%
  ungroup() %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  filter(origin != "IPG" & origin != "IW") %>%
  arrange(desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) %>%
  filter(Species != "unknown_forb") 

pred_df_raw_grazing <- do.call("rbind", predY_grazing) %>%
  as_tibble() %>%
  mutate(grazing_intensity = rep(grazing_gradient$XDataNew$grazing_intensity, 
                                 n_runs),
         run = rep(1:n_runs,each=20)) %>%
  pivot_longer(values_to = "cover", names_to = "Species",
               -c(grazing_intensity,run))%>%
  # mutate(Species = str_replace_all(Species," ", "_")) %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  arrange(origin,desc(prevalence)) %>%
  filter(Species != "unknown_forb")%>%
  mutate(Species_f = factor(Species, levels = unique(.$Species)))

grazing_native <- pred_df_raw_grazing  %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "N") %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin]) %>%
  # mutate(Species_f = lut_sp_nice[Species_f])%>%
  ggplot(aes(x=grazing_intensity, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run, color = origin), key_glyph="rect")+
  geom_line(data = pred_df_grazing  %>%
              # mutate(Species_f = lut_sp_nice[Species_f])%>%
              filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
              filter(str_sub(origin,1,1) == "N") ,
            lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  xlab("Grazing Intensity (AUM/ha)") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  # scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_nat)+
  theme(legend.position = "right",
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.title = element_blank())

grazing_introduced <- pred_df_raw_grazing  %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "I") %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin])%>%
  ggplot(aes(x=grazing_intensity, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run, color = origin), key_glyph="rect")+
  geom_line(data = pred_df_grazing %>%
              filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
              filter(str_sub(origin,1,1) == "I") ,
            lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  xlab("Grazing Intensity (AUM/ha)") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  # scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_int)+
  theme(legend.position = "right",
        legend.justification = c(1,0),
        panel.border = element_rect(fill=NA, size=0.75),
        strip.text = element_markdown(),
        legend.title = element_blank())

# gradient burned ==================
m$XData$burned <- as.factor(m$XData$burned)
burned_gradient = constructGradient(m, focalVariable = "burned")

predY_burned = predict(m, XData=burned_gradient$XDataNew, 
                       studyDesign=burned_gradient$studyDesignNew, 
                       ranLevels=burned_gradient$rLNew, expected=TRUE)


n_runs <- nChains*samples


pred_df_burned <- do.call("rbind", predY_burned) %>%
  as_tibble() %>%
  mutate(burned = rep(burned_gradient$XDataNew$burned, n_runs)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -burned) %>%
  group_by(Species, burned) %>%
  dplyr::summarise(mean = median(cover),
                   upr = quantile(cover, 0.975),
                   lwr = quantile(cover, 0.025)) %>%
  ungroup() %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  filter(origin != "IPG" & origin != "IW") %>%
  arrange(origin,desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) %>%
  filter(Species != "unknown forb") 

pred_df_raw_burned <-do.call("rbind", predY_burned) %>%
  as_tibble() %>%
  mutate(burned = rep(burned_gradient$XDataNew$burned, n_runs),
         run = rep(1:n_runs,each=2)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -c(burned,run))%>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  filter(origin != "IPG" & origin != "IW") %>%
  arrange(origin, desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) %>%
  filter(Species != "unknown forb")  %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin])

burned_introduced<-pred_df_raw_burned %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "I")  %>%
  mutate(burned = ifelse(burned == "yes", "B", "U")%>%
           factor(levels = c("U", "B"))) %>%
  ggplot(aes(x=burned, y=cover)) +
  geom_jitter(alpha = 0.03, aes(color = origin))+
  geom_boxplot(fill = "transparent", outlier.shape = NA)+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  xlab("Fire Occurrence") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_int)+
  theme(legend.position = c(1,0),
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())

burned_native<-pred_df_raw_burned %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "N")  %>%
  mutate(burned = ifelse(burned == "yes", "B", "U")%>%
           factor(levels = c("U", "B"))) %>%
  ggplot(aes(x=burned, y=cover)) +
  geom_jitter(alpha = 0.03, aes(color = origin))+
  geom_boxplot(fill = "transparent", outlier.shape = NA)+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  xlab("Fire Occurrence") +
  ylab("Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_nat)+
  theme(legend.position = c(1,0),
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        legend.justification = c(1,0),
        legend.background = element_rect(fill="transparent"),
        legend.title = element_blank())
# gradient bromus cover ==========

bromus_gradient = constructGradient(m, focalVariable = "B_tectorum", 
                                    non.focalVariables=list(burned = list(1)))

predY_bromus = predict(m, XData=bromus_gradient$XDataNew, 
                       studyDesign=bromus_gradient$studyDesignNew, 
                       ranLevels=bromus_gradient$rLNew, expected=TRUE)

n_runs <- nChains*samples

pred_df_bromus <- do.call("rbind", predY_bromus) %>%
  as_tibble() %>%
  mutate(B_tectorum = rep(bromus_gradient$XDataNew$B_tectorum, n_runs)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -B_tectorum) %>%
  group_by(Species, B_tectorum) %>%
  dplyr::summarise(mean = median(cover),
                   upr = quantile(cover, 0.975),
                   lwr = quantile(cover, 0.025)) %>%
  ungroup() %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  filter(origin != "IPG" & origin != "IW") %>%
  arrange(desc(prevalence)) %>%
  mutate(Species_f = factor(Species, levels = unique(.$Species))) %>%
  filter(Species != "unknown_forb") 

pred_df_raw_bromus <- do.call("rbind", predY_bromus) %>%
  as_tibble() %>%
  mutate(B_tectorum = rep(bromus_gradient$XDataNew$B_tectorum, n_runs),
         run = rep(1:n_runs,each=20)) %>%
  pivot_longer(values_to = "cover", names_to = "Species", -c(B_tectorum,run))%>%
  # mutate(Species = str_replace_all(Species," ", "_")) %>%
  left_join(prevalence) %>%
  filter(prevalence > 8) %>%
  mutate(origin = lut_all_fg[Species]) %>%
  arrange(origin,desc(prevalence)) %>%
  filter(Species != "unknown_forb")%>%
  mutate(Species_f = factor(Species, levels = unique(.$Species)))

bromus_introduced <- pred_df_raw_bromus  %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "I") %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin])%>%
  ggplot(aes(x=B_tectorum, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run, color = origin), key_glyph="rect")+
  geom_line(data = pred_df_bromus %>%
              filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
              filter(str_sub(origin,1,1) == "I"),
            lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  labs(x="*Bromus tectorum* cover (%)",
       y= "Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_int)+
  # scale_color_manual(values = c("#FFC845", "#007DBA"), 
  #                    labels = c("Introduced", "Native"))+
  theme(legend.position = c(1,0),
        strip.text = element_markdown(),
        panel.border = element_rect(fill=NA, size=0.75),
        axis.title.x = element_markdown(),
        legend.justification = c(1,0),
        legend.title = element_blank())

bromus_native <- pred_df_raw_bromus  %>%
  filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
  filter(str_sub(origin,1,1) == "N") %>%
  mutate(origin = c("IAF" = "Introduced Annual Forb",
                    "IAG" = "Introduced Annual Grass",
                    "NAF" = "Native Annual Forb",
                    "NPF" = "Native Perennial Forb",
                    "NPG" = "Native Perennial Grass",
                    "NW" = "Native Shrub")[origin])%>%
  ggplot(aes(x=B_tectorum, y=cover)) +
  geom_line(alpha = 0.03, aes(group=run, color = origin), key_glyph="rect")+
  geom_line(data = pred_df_bromus %>%
              filter(origin != "IPG" & origin != "IW" & origin != "NC") %>%
              filter(str_sub(origin,1,1) == "N"),
            lwd=1, alpha=0.95, color="black", aes(y=mean))+
  facet_wrap(~Species_f, nrow=1, labeller = labeller(Species_f = lut_sp_nice))+
  labs(x="*Bromus tectorum* cover (%)",
       y= "Probability of Occurrence") +
  guides(color=guide_legend(override.aes = list(alpha=1)))+
  scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = pal_nat)+
  # scale_color_manual(values = c("#FFC845", "#007DBA"), 
  #                    labels = c("Introduced", "Native"))+
  theme(legend.position = c(1,0),
        panel.border = element_rect(fill=NA, size=0.75),
        strip.text = element_markdown(),
        axis.title.x = element_markdown(),
        legend.justification = c(1,0),
        legend.title = element_blank())


# all together =====

nnn<-ggarrange(grazing_native, burned_native,bromus_native,  
               nrow = 3, common.legend = T, labels = c("(a)","(b)","(c)"),
               label.x = 0.01)+
  ggsave("images/probit_preds_nat_3pan.png",
         height = 7, width=12)

iii<-ggarrange(grazing_introduced, burned_introduced,bromus_introduced,  
               nrow = 3, common.legend = T, labels = c("(a)","(b)","(c)"),
               label.x = 0.01)+
  ggsave("images/probit_preds_int_3pan.png", 
         height = 7, width=12)
