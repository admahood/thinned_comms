# get some ancillary data
library(sf)
library(tidyverse)
library(elevatr)

files <- list.files("data", pattern = "xlsx$", full.names = TRUE)

# annual precip
# vsurf on diversity, fg cover
# beta diversity on the steps
# for both sites, pre-treatment was done august, post-treatment in june
# waiting on heterogeneity
# herbaceous response as far as loading
# cover percentage by class (herbaceous vs shrub)

plotvisit <- lapply(files, readxl::read_xlsx, "Export_TBL_PlotVisit") |>
  bind_rows() |>
  st_as_sf(coords = c("UTME", "UTMN"), crs =32613) |>
  mutate(project = ifelse(str_sub(PlotCode,1,1)=="E", "ESV", "PHA"))

dem_esv <- elevatr::get_aws_terrain(plotvisit[1,],z = 11, prj = st_crs(plotvisit))
dem_pha <- elevatr::get_aws_terrain(plotvisit[200,],z = 11, prj = st_crs(plotvisit))

# visualizing to make sure crs's are lined up etc
esv_elv <- plotvisit |>
  dplyr::select(Elevation, project) |>
  filter(project == "ESV") %>%
  mutate(elv = terra::extract(dem_esv, ., ID=F) %>% unlist() %>% as.numeric(),
         Elevation = Elevation * 0.3048)
me<- mblm::mblm(elv ~ Elevation, data = esv_elv)

ggplot(esv_elv, aes(x=Elevation, y=elv)) +
  geom_point() +
  ggtitle("ESV") +
  geom_abline(slope = 1, intercept = 0) +
  geom_line(aes(y = predict(me)), color ="red")


pha_elv <- plotvisit |>
  dplyr::select(Elevation, project) |>
  filter(project == "PHA") %>%
  na.omit() %>%
  mutate(elv = terra::extract(dem_pha, ., ID=F) %>% unlist() %>% as.numeric(),
         Elevation = Elevation * 0.3048)

mp<- mblm::mblm(elv ~ Elevation, data = pha_elv)
ggplot(pha_elv, aes(x=Elevation, y=elv)) +
  geom_point() +
  ggtitle("PHA")+
  geom_abline(slope = 1, intercept = 0)+
  geom_line(aes(y = predict(mp)), color ="red")

glimpse(plotvisit)
plot(dem_esv); plot(plotvisit[0], add=T)
ggplot(plotvisit |> filter(project == "ESV")) +
  geom_sf(aes(color = PlotTreatmentStatus)) 

 plot(dem_pha); plot(plotvisit[0], add=T)
 ggplot(plotvisit |> filter(project == "PHA")) +
  geom_sf(aes(color = PlotTreatmentStatus)) 

 
# usda plant trait data
checklist <- read_csv("data/usda_plantlist.csv") %>%
  janitor::clean_names()
read_csv("data/SearchResults.csv", skip = 3) %>%
  janitor::clean_names() %>%
  mutate(growth_habit = "forb") %>%
  write_csv("data/usda/co_forbs.csv")
ptr_search
