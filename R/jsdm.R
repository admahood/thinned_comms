# setup ========================================================================
source("R/a_tc_data_prep.R")

library(Hmsc)
library(parallel)
library(gghmsc)
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
if(!dir.exists("data/hmsc")) dir.create('data/hmsc')
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
 
}else{load(hmsc_file)}
