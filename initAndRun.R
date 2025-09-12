source("G:/Home/MyTools/myFunctions.R")
source("G:/Home/MyTests/reclassModel/modules/common/randomInputs.R")
library(SpaDES)
library(SpaDES.core)
library(SpaDES.tools)
library(mapview)
library(terra)
library(data.table)
library(caret)

setBasePath("G:/Home/MyTests/reclassModel")
getPaths() # shows where the 4 relevant paths are

sim <- simInit(times = list(start = 0, end = 0),
               modules = list("WB_VegBasedDrainage"),
               params = list(WB_VegBasedDrainage = list(
                 WB_VegBasedDrainageTimeStep = 1)))

sim <- spades(sim)


terra::writeRaster(sim$WB_VBD_ClayMap, "G:/Home/temp/clay_250m.tif", overwrite = TRUE)
terra::writeRaster(sim$pixelGroupMap, "G:/Home/temp/pgm_250m.tif", overwrite = TRUE)
terra::writeRaster(sim$TWIMap, "G:/Home/temp/TWI_250m.tif", overwrite = TRUE)
