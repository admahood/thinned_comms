###############################################################################
## Code to calculate heat load index following eq 3 in McCune and Keon (2002)

## Code by Kyle C. Rodman, Ecological Restoration Institute. 
# 10/14/2024

###############################################################################
## Read in packages
package.list <- c("terra", "here")
new.packages <- package.list[!(package.list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## And loading them
for(i in package.list){library(i, character.only = T)}

###############################################################################
## Reading in the data
slope <- rast(here("Data", "Spatial", "SoilsAndTerrain", "slope_deg.tif"))
slope[slope > 60] <- 60 ## HLI Equation 3 only works if slope is <= 60, setting higher values to 60.
aspect <- rast(here("Data", "Spatial", "SoilsAndTerrain", "aspect_deg.tif"))
latitude <- rast(here("Data", "Spatial", "SoilsAndTerrain", "latitude_deg.tif"))
out_path <- here("Data", "Spatial", "SoilsAndTerrain", "hli_eq3.tif")

## Calculating HLI
print("Step 1/11: Transforming latitude to radians")
tmp0 <- terra::app(latitude, fun = function(x) {x * 0.017453293})
remove(latitude)

print("Step 2/11: Calculating cosine of latitude")
cl <- terra::app(tmp0, fun = function(x) {cos(x)})

print("Step 3/11: Calculating sine of latitude")
sl <- terra::app(tmp0, fun = function(x) {sin(x)})
remove(tmp0)

print("Step 4/11: Transforming slope to radians")
tmp1 <- terra::app(slope, fun = function(x) {x * 0.017453293})
remove(slope)

print("Step 5/11: Transforming aspect to radians")
tmp2 <- terra::app(aspect, fun = function(x) {x * 0.017453293})
remove(aspect)

print("Step 6/11: Calculating folded aspect")
tmp3 <- terra::app(tmp2, fun = function(x) {
  abs(pi - abs(x - (5*pi/4)))
}) ## Calculating this based on "southwestness" as this is
remove(tmp2)

print("Step 7/11: Calculating cosine of slope")
tmp4 <- terra::app(tmp1, fun = cos)

print("Step 8/11: Calculating sine of slope")
tmp5 <- terra::app(tmp1, fun = sin)
remove(tmp1)

print("Step 9/11: Calculating cosine of folded aspect")
tmp6 <- terra::app(tmp3, fun = cos)
remove(tmp3)

print("Step 10/11: Final output file")
final_layer <- terra::app(x = c(tmp4, tmp5, tmp6, cl, sl), fun = function(x){
  0.339 + (0.808 * (x[4] * x[1])) - (0.196 *(x[5] * x[2])) - (0.482 * (x[3] * x[2]))
})
remove(tmp4, tmp5, tmp6, cl, sl)

print("Step 11/11: Saving output")
writeRaster(final_layer, out_path, overwrite = T)
