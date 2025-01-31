#######################################################################################
## This code performs a multi-step procedure that spatially downscales gridded climate
  # data using GIDS (Gradient and Inverse Distance-Squared) of Nalder and Wein (1998)and
  # Flint and Flint (2012). Used to downscale Daymet temp and precip data for Pinyon-
  # Juniper Woodlands north of the San Francisco Peaks, AZ

## Code by Kyle C. Rodman, Ecological Restoration Institute. 
# 12/9/2024

###################################################################################
## Import necessary packages
###################################################################################
package.list <- c("terra", "parallel", "FNN", "tidyverse", "sf", "here")
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###################################################################################
## Reading in data
###################################################################################
## Bounds of study area
studyBound <- st_read(here("Data", "Spatial", "SiteLocations", "SiteBound.shp"))

## DEMs of different resolutions
dem30 <- rast(here("Data", "Spatial", "SoilsAndTerrain", "elev30m_m.tif"))
dem90 <- rast(here("Data", "Spatial", "SoilsAndTerrain", "elev90m_m.tif"))
dem270 <- rast(here("Data", "Spatial", "SoilsAndTerrain", "elev270m_m.tif"))
dem1000 <- rast(here("Data", "Spatial", "SoilsAndTerrain", "elev1km_m.tif"))

## Get Daymet data
pptData <- rast(here("Data", "Spatial", "Daymet_coarse", "daymet_prcp.tif"))
tminData <- rast(here("Data", "Spatial", "Daymet_coarse", "daymet_tmin.tif"))
tmaxData <- rast(here("Data", "Spatial", "Daymet_coarse", "daymet_tmax.tif"))

## NOTE: 372 bands in old dataset vs. 408 now

###################################################################################
## Creating the functions
###################################################################################

## Function to perform GIDS. Called from main function below
gids <- function(coarse_data, fine_data, dist_mat){
  ## Reformatting matrices to allow for model fitting and prediction
  coarse_data <- as.data.frame(coarse_data)
  fine_data <- as.data.frame(fine_data)
  worker <- function(i, fine = fine_data, coarse = coarse_data, dists = dist_mat){
    
    ## First, subset data to get local neighbors
    fine_row <- fine[i,]
    ind <- dists[i,,1]
    dist <- dists[i,,2]
    data <- coarse[c(ind),]
    remove(fine, coarse, dists)
    
    ## Fit lm to data using base-level c function. Quite a bit faster than lm()
    x <- model.matrix(~x + y + anc_layer, data = data)
    y <- data$ds_layer
    fit <- .lm.fit(x = x, y = y) # Fit multiple regression model
    coef <- coefficients(fit)[2:4] # Extract coefficients from regression model
    data <- data[5:50,] # Subset data to exclude cells within nugget distance of 1-
      # cell resolution. Remove “goose egg” effect in interpolation
    x_coord <- fine_row$x; y_coord = fine_row$y; anc = fine_row$anc_layer # renaming
      # variables to make the lines (below) shorter
    
    ## Different portions of the GIDS Formula
    d_sq <- dist[5:50]^2
    sum1 <- sum((data$ds_layer + (x_coord - data$x)*coef[1] + (y_coord -
                data$y)*coef[2] + (anc - data$anc_layer)*coef[3])/(d_sq), na.rm = T)
    sum2 <- sum(1/d_sq, na.rm = T) # Doing this step outside the loop would
      # probably increase speed, but also make things harder for me to follow
    ## And returning output
    return(sum1/sum2)
  }
  ## Running worker function on each point and adding as new column to df
  fine_data$ds_layer <- sapply(seq_len(nrow(fine_data)),worker)
  ## Ouptutting df in matrix form for later
  return(as.matrix(fine_data))
}

## Main function. Wrapper for GIDS that formats data and outputs downscaled grid
multiLevelInterpParrallel <- function(boundary, ds_layer, ancillary_list,
                                      clip_dists = c(20000, 4500, 1350, 450, 150), file_out = NULL){
  
  ## Define number of cores to be used - 6 seems to be about the right balance for minimal export time and reduced processing time
  no_cores <- 6

  ## Creating projected version of boundary .shp that corresponds with output projection
  boundary_proj <- st_transform(boundary, crs(ancillary_list[[1]]))
  
  ## Looping through each of the DEMs
  for(iter in 1:(length(ancillary_list)-1)){
    if(iter == 1){
      
      ## In first iteration of loop, we need to start with the initial ds_layer and
        # ancillary dataset rather than interpolated data
      ds_layer <- terra::project(ds_layer, ancillary_list[[1]]) ## Projecting to
        # align with ancillary data grid
      ds_layer <- terra::crop(ds_layer, extend(ext(boundary_proj), clip_dists[iter]))
      
      ## Cropping to study bounds
      anc_layer <- terra::crop(ancillary_list[[1]], extend(ext(boundary_proj),
                                                    clip_dists[iter])) ## Cropping anc data to study bounds
      
      ## Removing cells with NAs in either layer
      ds_layer[is.na(anc_layer)] <- NA
      anc_layer[is.na(ds_layer)] <- NA
      
      ## Converting to non-referenced matrices
      coarse_pts <- as.data.frame(anc_layer, xy = T)
      coarse_pts <- cbind(coarse_pts, as.data.frame(ds_layer, xy = F)[,1])
      dimnames(coarse_pts)[[2]] <- c("x", "y", "anc_layer", "ds_layer")
      remove(ds_layer, anc_layer)
    
    }else{
      ## More efficient to use previously processed point data in later iterations
      coarse_pts <- out_pts
      remove(out_pts)
    }
    ## Finding nearest neighbors from coarse point layer for each point in finescale ancillary layer
    
    # Cropping raster of finer-scale ancillary data
    new_anc_layer <- crop(ancillary_list[[iter+1]], extend(ext(boundary_proj),
                                                           clip_dists[iter+1]))
    layExt <- ext(new_anc_layer); refer <- crs(new_anc_layer)
    
    # Converting it to a matrix
    new_pts <- as.data.frame(new_anc_layer, xy = T)
    remove(new_anc_layer)
    
    # Defining column names of matrix
    dimnames(new_pts)[[2]] <- c("x", "y", "anc_layer")
    
    # Finding the 50 nearest neighbors and getting distances to each.
      # Based on GIDS typically using 7x7 window
    nns <- get.knnx(coarse_pts[,1:2],new_pts[,1:2], k = 50)
    
    # Reformatting list to 3d array. A little easier to work with later
    nns <- sapply(nns, identity, simplify="array")
    
    ## Splitting data for parrallelization, and initializing cluster
    parts <- split(x = seq_len(nrow(new_pts)), f = 1:no_cores)
    cl <- makeCluster(no_cores)
    clusterExport(cl = cl, varlist = c("coarse_pts", "new_pts", "nns", "parts",
                                       "gids"), envir = environment())
    
    ## Running GIDS, split between number of cores on PC minus 1
    parallelX <- parLapply(cl = cl, X = 1:no_cores,
                           fun = function(t) gids(coarse_data = coarse_pts,
                                                  fine_data = new_pts[parts[[t]],],
                                                  dist_mat = nns[parts[[t]],,]))
    
    ## Terminating cluster
    stopCluster(cl)
    
    ## Merging parallel output to single matrix
    out_pts <- do.call(rbind, parallelX)
    # Write raster to file if we made it through all iterations and if path for outfile provided
    if(iter == (length(ancillary_list)-1)){
      ## Creating raster from matrix
      downscaled <- rast(out_pts[,c(1,2,4)], extent = layExt, crs = refer, type = "xyz")
      if(!is.null(file_out)){
        writeRaster(downscaled, file_out, overwrite = T)
      }
    }
    ## Keeping track of progress
    print(paste("Done with downscale", iter, "out of", (length(ancillary_list)-1)))
  }
  ## Return function output to global environment for later use. Comment out if 
    # just writing output to disk
  return(downscaled)
}

###################################################################################
###################################################################################
### Running code in each study area to perform downscaling
###################################################################################

## Formatting study boundary data. Makes it easier to run function later on
bound <- st_buffer(studyBound, dist = -4000)

## Run the downscaling for each layer in TMin
# Clean up temp workspace
unlink(list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T))

# Run downscaling for each layer
for(x in 373:nlyr(pptData)){
  temp <- multiLevelInterpParrallel(boundary = bound, ds_layer = pptData[[x]],
                                    ancillary_list = list(dem1000, dem270, dem90, dem30),
                                    clip_dists = c(5000, 1350, 450, 150))
  number <- str_pad(as.character(x), 3, "left", "0")
  writeRaster(temp, here("Data", "Spatial", "WorkingDirectory", 
                         paste("PPT_", "month", number, ".tif", sep = "")), overwrite = T)
  print(paste("Done with layer", x, "out of", nlyr(pptData)))
}

# Get names of each produced layer
filesWD <- list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T)

# Get only the ".tif" files
subFiles <- filesWD[str_sub(filesWD, start = -4, end = -1) == ".tif"]

# Read all of them, convert to a multi-band stack, and write output
allGrids <- rast(subFiles)
writeRaster(allGrids, here("Data", "Spatial", "Daymet_downscaled", "PPT.tif"))
print("  DONE WITH TMIN")
unlink(list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T))

## Run the downscaling for each layer in TMin
# Clean up temp workspace
unlink(list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T))

# Run downscaling for each layer
for(x in 373:nlyr(tminData)){
  temp <- multiLevelInterpParrallel(boundary = bound, ds_layer = tminData[[x]],
                                    ancillary_list = list(dem1000, dem270, dem90, dem30),
                                    clip_dists = c(5000, 1350, 450, 150))
  number <- str_pad(as.character(x), 3, "left", "0")
  writeRaster(temp, here("Data", "Spatial", "WorkingDirectory", 
                         paste("TMin_", "month", number, ".tif", sep = "")), overwrite = T)
  print(paste("Done with layer", x, "out of", nlyr(tminData)))
}

# Get names of each produced layer
filesWD <- list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T)

# Get only the ".tif" files
subFiles <- filesWD[str_sub(filesWD, start = -4, end = -1) == ".tif"]

# Read all of them, convert to a multi-band stack, and write output
allGrids <- rast(subFiles)
writeRaster(allGrids, here("Data", "Spatial", "Daymet_downscaled", "TMin.tif"))
print("  DONE WITH TMIN")
unlink(list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T))

## Run the downscaling for each layer in TMax
# Clean up temp workspace
unlink(list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T))

# Run downscaling for each layer
for(x in 373:nlyr(tmaxData)){
  temp <- multiLevelInterpParrallel(boundary = bound, ds_layer = tmaxData[[x]],
                                    ancillary_list = list(dem1000, dem270, dem90, dem30),
                                    clip_dists = c(5000, 1350, 450, 150))
  number <- str_pad(as.character(x), 3, "left", "0")
  writeRaster(temp, here("Data", "Spatial", "WorkingDirectory", 
                         paste("TMax_", "month", number, ".tif", sep = "")), overwrite = T)
  print(paste("Done with layer", x, "out of", nlyr(tmaxData)))
}

# Get names of each produced layer
filesWD <- list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T)

# Get only the ".tif" files
subFiles <- filesWD[str_sub(filesWD, start = -4, end = -1) == ".tif"]

# Read all of them, convert to a multi-band stack, and write output
allGrids <- rast(subFiles)
writeRaster(allGrids, here("Data", "Spatial", "Daymet_downscaled", "TMax.tif"))
print("  DONE WITH TMAX")
unlink(list.files(here("Data", "Spatial", "WorkingDirectory"), full.names = T))
