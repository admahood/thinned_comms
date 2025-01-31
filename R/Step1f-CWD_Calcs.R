###############################################################################
## Code to calculate climatic water deficit (CWD) following Stephenson (1998) &
  # Lutz et al. (2010). Some code adapted from M. Redmond's script to calculate CWD
  # Adapted to handle both coarse-scale (e.g., 4km/800m) and downscaled 
  # (e.g., 250-m/30-m) climate data

###############################################################################
## Read in packages
package.list <- c("terra", "sf", "here", "tidyverse", "exactextractr")
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###############################################################################
## Reading in the data

# Extent of study area
studySites <- st_read(here("Data", "Spatial", "SiteLocations", "SiteBound.shp"))
studySites <- st_buffer(studySites, dist = -4000)

# Plot locations
plots <- st_read(here("Data", "Spatial", "SiteLocations", "plotBounds_5_11_22_NAD83.shp")) %>%
  select(Tran_Sect)

# Heat load index raster (following McCune and Keon (2002) and Lutz et al. 2010)
hli_raster <- rast(here("Data", "Spatial", "SoilsAndTerrain", "hli_eq3.tif"))

# Raster of latitude values throughout study extent - used to calculate day length
lat_raster <- rast(here("Data", "Spatial", "SoilsAndTerrain", "latitude_deg.tif"))

# Fraction AWC/mean soil depth from Polaris. Both resampled to 30m prior to calculations 
soil_awc <- rast(here("Data", "Spatial", "SoilsAndTerrain", "awc_mm.tif"))

# Get downscaled PRISM data files and filter by file name
ppt_data <- rast(here("Data", "Spatial", "Daymet_downscaled", "PPT.tif"))
tmin_data <- rast(here("Data", "Spatial", "Daymet_downscaled", "TMin.tif"))
tmax_data <- rast(here("Data", "Spatial", "Daymet_downscaled", "TMax.tif"))

## Defining output directory for final files
out_directory <- here("Data", "Spatial", "WaterBalance")

###############################################################################
### Running water balance calcs

## Inputs needed will be: heat load index raster, a raster with latitude values in each cell,
  # soil_awc as mm in top 2 m of soil column, as well as ordered raster stacks of ppt,
  # tmax, and tmin, for each month in the calculation. Lastly, an output directory to
  # store the monthly CWD layer outputs. All files (except climate) should be of 
  # the same resolution and projection, and snapped to the same reference grid.

print("Cropping data to site extent")
## Cropping other layers by extent of climate data
site_ext <- ext(ppt_data)
hli_raster <- crop(hli_raster, site_ext)
lat_raster <- crop(lat_raster, site_ext)
soil_awc <- crop(soil_awc, site_ext)

print("Calculating day length from latitude")
## First, we need to calculate average daylength for every cell in a monthly timestep,
  # relevant code from M. Redmond's CWD script and from daylength function in geosphere package
days_month <- c(31,28,31,30,31,30,31,31,30,31,30,31)
dl_raster_list <- lapply(1:12, function(i){
  ## Modifying default values in this function on each loop. Makes it easier to use with "app"
  worker_fun <- function(x, month = i){
    doy <- c(15,45,74,105,135,166,196,227,258,288,319,349)[month]
    P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.0086 * (doy - 186)))))
    a <- (sin(0.8333 * pi/180) + sin(x * pi/180) * sin(P))/(cos(x * pi/180) * cos(P))
    a <- pmin(pmax(a, -1), 1)
    return(24 - (24/pi) * acos(a))
  }
  return(terra::app(lat_raster, fun = worker_fun))
})

## Looping through each iteration to calculate CWD in several time steps
for(time_step in 1:nlyr(ppt_data)){
  print(paste("Calculating AET/CWD for time step", time_step, "out of", nlyr(ppt_data)))
  ## Calculating tmean for this time step
  tmean_temp <- terra::app(c(tmax_data[[time_step]], tmin_data[[time_step]]), fun = mean)
  
  ## Reprojecting 30-m climate grids to align with 10-m terrain variables using "Average". Doesn't modify values here, just changes resolution to run raster calcs
  tmean_temp <- terra::project(tmean_temp, hli_raster, gdal = T, method = "average")
  ppt_temp <- terra::project(ppt_data[[time_step]], hli_raster, gdal = T, method = "average")
  
  ## Creating estimates of rain fraction for later use. Rain fraction is from Thornthwaite method
    # We assume that the proportion of rain vs. snow received is a function of mean 
    # monthly temperature, with 0 and 6 degrees C being important thresholds
  rain_frac <- tmean_temp; 
  rain_frac[tmean_temp<0] <- 0 
  rain_frac[(tmean_temp>0)&(tmean_temp<6)] <- 
    rain_frac[(tmean_temp>0)&(tmean_temp<6)] * 0.166666666 
  rain_frac[tmean_temp>=6] <- 1
  
  ## Snow and rain are separated into two components from above estimates
  ppt_temp[ppt_temp<0] <- 0 ## Correcting negative values of ppt, sometimes present b/c of downscaling
  snow <- terra::app(c(rain_frac, ppt_temp), fun = function(x){(1-x[[1]])*x[[2]]})
  rain <- terra::app(c(rain_frac, ppt_temp), fun = function(x){x[[1]]*x[[2]]})
  
  ## If we are at the first time step, start with zero initial snowpack and soil moisture
  if(time_step == 1){
    pack <- rain_frac
    values(pack) <- 0
    soilM_tmin1 <- pack
  }
  
  ## Snowmelt and change in new snowpack are follows. Then, remove negative values
  melt <- terra::app(c(rain_frac, snow, pack), fun = function(x){x[[1]]*(x[[2]]+x[[3]])})
  pack <- terra::app(c(rain_frac, ppt_temp, pack), fun = function(x){
    (((1-x[[1]])^2) * x[[2]]) + (1-x[[1]]) * x[[3]]
  })
  pack[pack<0] <- 0
  
  ## Predicting water supply in cell as function of incoming rain and snowmelt
  supply <- terra::app(c(rain, melt), fun = sum)
  
  ## Getting number of days in month for PET calcs
  if(time_step%%12 == 0){dl_index = 12}else{dl_index = time_step%%12}
  ndays <- days_month[dl_index]
  
  ## Calculating saturation vapor pressure (SVP) from temps in each cell - in kPA
    # Formula from Dingman (2002, gave units in pa) and Lutz et al. 2010 (units in kPA)
  ea <- terra::app(tmean_temp, fun = function(x){
    (exp((17.3*x)/(x+237.3)))*0.611
  })
  
  ## Actually calculating PET using day length, number of days in month, heat load index, temp, SVP, and a couple constants
    # worker function for raster calc
  petFunction <- function(x, days = ndays){
    29.8 * days * x[[1]] * x[[2]] * (x[[3]]/(x[[4]] + 273.2))
  }
    # calculating PET 
  pet <- terra::app(c(dl_raster_list[[dl_index]], hli_raster, ea, tmean_temp), fun = petFunction)
    # Setting negative values of PET to zero if they are present
  pet[pet<=0] <- 0 ## Assuming zero PET in negative temps is correct in thornthwaite but not Hamon eq.
  
  ## Next step adjusts soil moisture in given time step based on PET, supply, and AWC
  net <- terra::app(c(supply, pet), fun = function(x){x[[1]]-x[[2]]})
  soilM <- terra::app(c(soilM_tmin1, net), fun = sum)    
  
  ## Removing negative and overly positive values (based on soil AWC). Any soil
    # recharge outside of AWC is assumed to be runoff and unavailable to plants
  soilM[soilM>soil_awc] <- soil_awc[soilM>soil_awc]
  soilM[soilM<0] <- 0
  
  ## Calculating soil water flux
  deltaSoil <- terra::app(c(soilM_tmin1, pet, supply, soil_awc), 
                          fun = function(x){
                            x[[1]] * (1-(exp(-1*(x[[2]]-x[[3]])/(x[[4]]+0.000001))))
                          })
  soil_contrib <- deltaSoil
  soil_contrib[soil_contrib<0] <- 0
  
  ## At end of soil moisture calc, saving soil moisture as soilM from previous month (for next loop)
  soilM_tmin1 <- soilM
  
  ## Calculating AET from incoming supply and change in soil moisture. Basically,
    # the maximum possible AET is the sum of water supply and change in soil moisture.
    # However, aet can not exceed PET, so this is adjusted to equal PET when potential
    # evapotranspiration is less than the combination of supply and deltaSoil.
  aet <- terra::app(c(supply, soil_contrib), fun = sum)
  aet <- terra::app(c(pet, aet), fun = min)
  
  ## Writing out monthly AET and CWD grids to temp. This is done differently for normals and annual data.
    # However, in both cases, CWD is calculated as the difference between PET and AET. We
    # round to two digits to make sure that we don't end up with weird values from subtracting
    # rasters that held a different number of significant digits
  
  ## Getting month # in consecutive order with leading 0s
  number <- str_pad(as.character(time_step), 3, "left", "0")

  ## CWD
  terra::app(c(pet, aet), fun = function(x){
    round(x[[1]] - x[[2]], digits = 2)
  }, filename = here("Data", "Spatial", "WorkingDirectory", 
                    paste("CWD_", "month", number, ".tif", sep = "")), overwrite = T)
  
  ## AET
  terra::app(aet, fun = function(x){
    round(x, digits = 2)
  }, filename = here("Data", "Spatial", "WorkingDirectory", 
                    paste("AET_", "month", number, ".tif", sep = "")), overwrite = T)
}

print("Aggregating monthly data to annual grids, by water year")
## Aggregate monthly values of AET by water year from 1981-2021
  # Sum by year

## Get lists of all CWD and AET files in temp folder
aet_list <- rast(list.files(here("Data", "Spatial", "WorkingDirectory"), pattern = "AET_", full.names = T))
cwd_list <- rast(list.files(here("Data", "Spatial", "WorkingDirectory"), pattern = "CWD_", full.names = T))

###############################################################################
### Taking monthly data and aggregating to annual grids

## Aggregating AET to water year
aetGrids <- lapply(1991:2023, function(year, aets = aet_list){
  # Layer index for end of the water year. First month is Oct 1989. We exclude 
    # the first twelve months because they start with 0 soil moisture and snow. 
    # We then sum over each 12-month water year starting with Oct 1990 to Sep 1991
  eowyIndex <- ((year-1990)*12)+12
  # Calculate yearly totals
  return(terra::app(aets[[c((eowyIndex-11):eowyIndex)]], fun = sum))
})
  # Aggregate list of rasters to stack
aetGrids <- rast(aetGrids)
names(aetGrids) <- paste("AET", 1991:2023, sep = "_")

## Aggregating CWD by water year
cwdGrids <- lapply(1991:2023, function(year, cwds = cwd_list){
  # Layer index for end of the water year. First month is Oct 1989. We exclude 
    # the first twelve months because they start with 0 soil moisture and snow. 
    # We then sum over each 12-month water year starting with Oct 1990 to Sep 1991
  eowyIndex <- ((year-1990)*12)+12
  # Calculate yearly totals
  return(terra::app(cwds[[c((eowyIndex-11):eowyIndex)]], fun = sum))
})
  # Aggregate list of rasters to stack
cwdGrids <- rast(cwdGrids)
names(cwdGrids) <- paste("CWD", 1991:2023, sep = "_")

## If output directory provided, save raster. If not, return to global environment
if(!is.null(out_directory)){
  writeRaster(aetGrids, here(out_directory, "AET_1991-2023.tif"),
              overwrite = T, datatype = 'FLT4S')
  writeRaster(cwdGrids, here(out_directory, "CWD_1991-2023.tif"),
              overwrite = T, datatype = 'FLT4S')
}else{
  ## If not saving outputs, return monthly AET/CWD outputs in two stacks by site
  bothDatasets <- list(aetGrids, cwdGrids)
  names(bothDatasets) <- c("AET_mm", "CWD_mm")
}

###############################################################################
### Extracting annual data, merging with plot list, and exporting as .csv

## Get annual data
aetGrids <- rast(here(out_directory, "AET_1991-2023.tif"))
cwdGrids <- rast(here(out_directory, "CWD_1991-2023.tif"))

## Extracting to plot locations
aetExtract <- exactextractr::exact_extract(aetGrids, plots, "mean")
colnames(aetExtract) <- str_sub(colnames(aetExtract), start = -8, end = -1)
cwdExtract <- exactextractr::exact_extract(cwdGrids, plots, "mean")
colnames(cwdExtract) <- str_sub(colnames(cwdExtract), start = -8, end = -1)

## Merging to one file for export
plots <- cbind(plots, aetExtract, cwdExtract)

## NOTE: Need to add extraction of annual values above...
write_csv(st_drop_geometry(plots), here("Data", "Processed", "AnnualClimate.csv"))

