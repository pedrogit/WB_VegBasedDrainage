defineModule(sim, list(
  name = "WB_VegBasedDrainage",
  description = paste("Construct a model for drainage in western boreal canada using plot data and other covariates and reuse it to predict drainage from Biomass_core simulated maps"),
  keywords = c("drainage", "western boreal"),
  authors =  c(
    person("Pierre", "Racine", email= "pierre.racine@sbf.ulaval.ca", role = "aut")
  ),
  childModules = character(0),
  version = list(WB_VegBasedDrainage = "0.0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  # citation = list("citation.bib"),
  # documentation = list("NEWS.md", "README.md", "WB_VegBasedDrainage.Rmd"),
  reqdPkgs = list("data.table", "reproducible", "LandR"),
  parameters = rbind(
    defineParameter("WB_VegBasedDrainageTimeStep", "numeric", 1, NA, NA,
                    "Simulation time at which the drainage map is regenerated.")
  ),
  inputObjects = rbind(
    expectsInput(objectName = "cohortData", objectClass = "data.table", desc = "cohortData from Biomass_core", sourceURL = NA),
    expectsInput(objectName = "pixelGroupMap", objectClass = "SpatRast", desc = "pixelGroupMap from Biomass_core on which to align input and drainage maps", sourceURL = NA),
    expectsInput(objectName = "plotData", objectClass = "data.table", desc = "Plot data to extimate model", sourceURL = NA),
    expectsInput(objectName = "clayMap", objectClass = "SpatRast", desc = "Clay raster from NRCan", sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Clay/Clay_X0_5_cm_100m1980-2000v1.tif"),
    expectsInput(objectName = "sandMap", objectClass = "SpatRast", desc = "Sand raster from NRCan", sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Sand/Sand_X0_5_cm_100m1980-2000v1.tif"),
    expectsInput(objectName = "siltMap", objectClass = "SpatRast", desc = "Silt raster from NRCan", sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Silt/Silt_X0_5_cm_100m1980-2000v1.tif"),
    expectsInput(objectName = "bulkDensityMap", objectClass = "SpatRast", desc = "Bulk Density raster from NRCan", sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Silt/Silt_X0_5_cm_100m1980-2000v1.tif")
  ),
  outputObjects = rbind(
    createsOutput(objectName = "WB_VegBasedDrainageMap", objectClass = "SpatRast", desc = "Drainage map predicted from the model")
  )
))

doEvent.WB_VegBasedDrainage = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, time(sim), "WB_VegBasedDrainage", "ReComputeDrainage")
    },
    ReComputeDrainage = {
      sim <- scheduleEvent(sim, time(sim) + P(sim)$WB_VegBasedDrainageTimeStep, "WB_VegBasedDrainage", "ReComputeDrainage")
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

### template initialization
Init <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

### template for your event1
ReComputeDrainage <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #
  # browser()
  # generate a fake pixelGroupMap if it is not supplied
  if(!suppliedElsewhere("pixelGroupMap", sim)){
    nbGroup <- 200
    pixelGroupRastWidth <- 1000
    pixelGroupRastWidth <- 1000
    message("pixelGrouMap not supplied. Please provide one. Creating random map ", 
            pixelGroupRastWidth, 
            " pixels by ", 
            pixelGroupRastWidth, 
            " pixels with ",
            nbGroup,
            " groups...")
    
    sim$pixelGroupMap <- getRandomPixelGroupMap(origin = c(1541912, 1072021),
                                                width = pixelGroupRastWidth,
                                                crs = "ESRI:102002",
                                                nbPixelGroup = nbGroup)
    # mapView(sim$pixelGroupMap)
  }
  
  # download and cache CANSIS soil data
  mapToProcess <- c("Clay", "Sand", "Silt", "BD")
  baseURL <- "https://sis.agr.gc.ca/cansis/nsdb/psm/"
  nameEnd <- "_X0_5_cm_100m1980-2000v1.tif"
  sapply(mapToProcess, function(mapName){
    varMapName <- paste0("WB_VBD_", mapName, "Map") # e.g. WB_VBD_clayMap
    if (!suppliedElsewhere(varMapName, sim)){
      sim[[varMapName]] <- Cache(
        prepInputs,
        url = paste0(baseURL, mapName, "/", mapName, nameEnd),
        targetFile = paste0(mapName, nameEnd),
        destinationPath = getPaths()$cache,
        fun = terra::rast,
        rasterToMatch = sim$pixelGroupMap,
        maskWithRTM  = TRUE,
        method = "average",
        overwrite = TRUE
      )
    }
  })
  

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
