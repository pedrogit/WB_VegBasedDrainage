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
    # expectsInput(objectName = "cohortData", 
    #              objectClass = "data.table", 
    #              desc = "cohortData from Biomass_core", 
    #              sourceURL = NA),
    # expectsInput(objectName = "pixelGroupMap", 
    #              objectClass = "SpatRast", 
    #              desc = paste("pixelGroupMap from Biomass_core on which to align ",
    #                           "input and drainage maps"), 
    #              sourceURL = NA),
    expectsInput(objectName = "WB_HartJohnstoneForestClassesMap",
                 objectClass = "SpatRast",
                 desc = paste("WB_HartJohnstoneForestClassesMap from the ",
                              "WB_HartJohnstoneForestClasses mudule used as ",
                              "standtype"),
                 sourceURL = NA),
    expectsInput(objectName = "plotPoints", 
                 objectClass = "data.table", 
                 desc = paste("Plot data to fit drainage model. Should contain a ",
                              "latitude, a longitude and at least one column named ",
                              "\"drainage\" with drainage classes. Can also contain a ",
                              "\"standtype\". Standtype will be extracted from the ",
                              "WB_HartJohnstoneMap if it was produced by the corresponding ",
                              "module is present."), 
                 sourceURL = NA),
    expectsInput(objectName = "TWIMap", 
                 objectClass = "SpatRast", 
                 desc = "TWI raster computed from Canada DEM", 
                 sourceURL = NA),
    expectsInput(objectName = "WB_VBD_ClayMap", 
                 objectClass = "SpatRast", 
                 desc = "Clay raster from NRCan", 
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Clay/Clay_X0_5_cm_100m1980-2000v1.tif"),
    expectsInput(objectName = "WB_VBD_SandMap", 
                 objectClass = "SpatRast", 
                 desc = "Sand raster from NRCan", 
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Sand/Sand_X0_5_cm_100m1980-2000v1.tif"),
    expectsInput(objectName = "WB_VBD_SiltMap", 
                 objectClass = "SpatRast", 
                 desc = "Silt raster from NRCan", 
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Silt/Silt_X0_5_cm_100m1980-2000v1.tif"),
    expectsInput(objectName = "WB_VBD_BDMap", 
                 objectClass = "SpatRast", 
                 desc = "Bulk Density raster from NRCan", 
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Silt/Silt_X0_5_cm_100m1980-2000v1.tif")
  ),
  outputObjects = rbind(
    createsOutput(objectName = "WB_VegBasedDrainageMap", 
                  objectClass = "SpatRast", 
                  desc = "Drainage map predicted from the model")
  )
))

doEvent.WB_VegBasedDrainage = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, time(sim), "WB_VegBasedDrainage", "ReComputeDrainageMap")
    },
    ReComputeDrainageMap = {
      sim <- scheduleEvent(sim, time(sim) + P(sim)$WB_VegBasedDrainageTimeStep, "WB_VegBasedDrainage", "ReComputeDrainageMap")
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

### template initialization
Init <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #
  # fit a model based on:
  #    standtype from sim$plotPoints if it has a standtype column or cohortData and pixelGroupMap otherwise.
  #    TWIMap, WB_VBD_ClayMap, WB_VBD_SandMap, WB_VBD_SiltMap, WB_VBD_BDMap (Bulk Density)
  # we assume all maps are projected to the same crs
  covariatesMaps <- c("TWIMap" = "twi", 
                      "WB_VBD_ClayMap" = "clay", 
                      "WB_VBD_SandMap" = "sand",
                      "WB_VBD_SiltMap" = "silt",
                      "WB_VBD_BDMap" = "bulk_den")
  
  sim$plotPoints <- sim$plotPoints[, !(names(sim$plotPoints) %in% unname(covariatesMaps))]
  
  for (i in seq_along(covariatesMaps)) {
    sim$plotPoints <- cbind(sim$plotPoints, extract(sim[[names(covariatesMaps)[i]]], sim$plotPoints)[, -1])  # Remove ID column from extract
    if (! covariatesMaps[i] %in% names(sim$plotPoints)){
      names(sim$plotPoints)[names(sim$plotPoints) == "y"] <- covariatesMaps[i]
    }
  }

  # Keep rows where specified columns are not NA
  covariatesMaps <- c(unname(covariatesMaps), "drainage")
  modelData <- as.data.frame(sim$plotPoints)
  keeps <- complete.cases(modelData[, covariatesMaps])
  modelData <- modelData[keeps, ]
  
  modelData <- modelData[, !(names(modelData) %in% c("X", "plot"))]
  
  
  # split the data frame into training and test data
  inTraining <- createDataPartition(modelData$drainage, p = 0.7, list = FALSE)
  trainSet <- modelData[inTraining, ]
  testSet <- modelData[-inTraining, ]
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = 5,
                             repeats = 5,
                             search = "random")
  
  # fit the model
  modelFit <- train(
    drainage ~ ., # Can consider subsets of covariates here say; e.g. `slope + elevation + age`; `.` considers all the predictors with no interaction terms
    data = trainSet,
    method = "rf", # "rf" = random forest or "xgbTree" = boosted regression trees, # See topepo.github.io/caret for more model tags that can be used here
    trControl = fitControl,
    tuneLength = 10,
    verbose = FALSE
  )
  
  sim$WB_VegBasedDrainageModel <- modelFit
  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

### template for your event1
ReComputeDrainageMap <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  ##############################################################################
  # Generate a fake WB_HartJohnstoneForestClassesMap if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("WB_HartJohnstoneForestClassesMap", sim)){
    rastWidth <- 1000
    message("WB_HartJohnstoneForestClassesMap not supplied. Please couple ",
            "with the WB_HartJohnstoneForestClasses module. Creating random map ",
            rastWidth,
            " pixels by ",
            rastWidth,
            " pixels for the 6 classes (\"deci\", \"mixed\", \"conimix\", ",
            "\"jackpine\", \"larch\" and \"spruce\")...")

    sim$WB_HartJohnstoneForestClassesMap <- Cache(
                           getRandomCategoricalMap,
                           origin = c(-667296, 1758502),
                           width = rastWidth,
                           crs = "ESRI:102002",
                           nbregion = 2000,
                           valuevect = 1:6,
                           seed = 100
    )
    # mapView(sim$WB_HartJohnstoneForestClassesMap)
  }
  
  ##############################################################################
  # Use our own plotPoints data if none is supplied
  ##############################################################################
  if(!suppliedElsewhere("plotPoints", sim) && !suppliedElsewhere("drainageModel", sim)){
    message("plotPoints not supplied. You must provide a CSV table with \"latitude\", ",
            "\"longitude\", \"standtype\" and \"drainage\" following the ",
            "WB_HartJohnstone classification. Loading default plot data points ",
            "as sim$plotPoints...") 
    plotFile <- file.path(getPaths()$modulePath, currentModule(sim), "data/plotData.csv")
    plotDF <- read.csv(plotFile)
    
    # purge standtype of drainage info
    if ("standtype" %in% names(plotDF)){
      plotDF$standtype <- sub("Poorly-drained ", "", plotDF$standtype)
      plotDF$standtype <- sub("Well-drained ", "", plotDF$standtype)
    }
    
    # convert character columns to factor
    plotDF[] <- lapply(plotDF, function(col){
      if (is.character(col)) factor(col) else col
    })
    
  # browser()
    # reassign factors to match WB_HartJohnstone classification used in the WB_WB_HartJohnstoneClasse module
    # labels = c("deci", "mixed", "conimix", "jackpine", "larch", "spruce")
    # levels = c(1L, 2L, 3L, 4L, 5L, 6L)
    new_codes <- c(3L, 1L, 2L, 4L, 5L, 6L)[as.integer(plotDF$standtype)]
    
    # Assign new codes while keeping levels the same
    plotDF$standtype <- factor(new_codes, levels = 1:6, labels = levels(plotDF$standtype))
      
    plotPoints <- vect(plotDF, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")  # WGS84
    sim$plotPoints <- project(plotPoints, crs(sim$WB_HartJohnstoneForestClassesMap))
    # writeVector(sim$plotPoints, "G:/Home/temp/plotPoints.shp", overwrite=TRUE)
  }
 
  ##############################################################################
  # Generate a TWI map from a massive downloaded DEM 
  ##############################################################################
  if(!suppliedElsewhere("TWIMap", sim)){
    library(whitebox)
    
    nbSteps <- 5
    
    message("Computing TWIMap (1/", nbSteps, "): Downloading/cropping/reprojecting/resampling/masking medium resolution MRDEM dem (80GB) to union of studyarea and 100km buffered plot points...") 
    # We will compute a WTI map for the area covered by the plotPoint (buffered 
    # to 100 km) and the studyarea (represented by the WB_HartJohnstoneForestClassesMap)
    
    plotPoints100KmBuffers <- aggregate(buffer(sim$plotPoints, width = 100000))  # 1000 m buffer (if CRS is in meters)
    
    # extract the WB_HartJohnstoneForestClassesMap extent
    WB_HartJohnstoneExtent <- ext(sim$WB_HartJohnstoneForestClassesMap)
    WB_HartJohnstoneExtentPoly <- vect(WB_HartJohnstoneExtent, "polygons")
    crs(WB_HartJohnstoneExtentPoly) <- crs(sim$WB_HartJohnstoneForestClassesMap)
    
    # merge them together
    sim$plotAndPixeGroupArea <- aggregate(rbind(plotPoints100KmBuffers, WB_HartJohnstoneExtentPoly))
    # mapView(sim$plotAndPixeGroupArea)
    # writeVector(sim$plotAndPixeGroupArea, file.path(getPaths()$cache, "plotAndPixeGroupArea.shp"), overwrite = TRUE)
    
    # define the output path
    cachePath <- getPaths()$cachePath
    plotAndPixeGroupAreaDemPath <- file.path(cachePath, "plotAndPixeGroupAreaDem.tif")
    
    # download and process the big thing
    # https://open.canada.ca/data/en/dataset/18752265-bda3-498c-a4ba-9dfe68cb98da

    plotAndPixeGroupAreaDem <- Cache(
      prepInputs,
      url = "https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30/mrdem-30-dtm.tif",
      targetFile = "mrdem-30-dtm.tif",
      destinationPath = cachePath,
      fun = terra::rast,
      cropTo = sim$plotAndPixeGroupArea,
      projectTo = extend(sim$WB_HartJohnstoneForestClassesMap, sim$plotAndPixeGroupArea),
      align_only= TRUE,
      maskTo = sim$plotAndPixeGroupArea,
      method = "bilinear",
      writeTo = plotAndPixeGroupAreaDemPath
    )

    # plotAndPixeGroupAreaDemExtPoly = vect(ext(plotAndPixeGroupAreaDem), "polygons")
    # crs(plotAndPixeGroupAreaDemExtPoly) <- crs(plotAndPixeGroupAreaDem)
    # mapview(plotAndPixeGroupAreaDemExtPoly)+mapview(plotAndPixeGroupArea)
    # writeRaster(plotAndPixeGroupAreaDem, file.path(getPaths()$cache, "plotAndPixeGroupAreaDem.tif"), overwrite = TRUE)
    
    # message("Computing TWIMap (3/", nbSteps, "): Writing plotAndPixeGroupAreaDem as a file in the cache folder for whitebox...")

    # define other paths
    dem_filled_path <- file.path(cachePath, "plotAndPixeGroupAreaDem_filled.tif")
    slope_path <- file.path(cachePath, "plotAndPixeGroupAreaDem_slope.tif")
    flow_acc_path <- file.path(cachePath, "plotAndPixeGroupAreaDem_flowAccum.tif")
    final_twi_path <- file.path(cachePath, "plotAndPixeGroupAreaDem_TWI.tif")
    
    message("Computing TWIMap (2/", nbSteps, "): Filling depressions...")
    Cache(wbt_fill_depressions,
          dem = plotAndPixeGroupAreaDemPath,
          output = dem_filled_path
    )
    
    message("Computing TWIMap (3/", nbSteps, "): Computing slopes...")
    Cache(wbt_slope,
          dem = plotAndPixeGroupAreaDemPath,
          output = slope_path,
          zfactor = 1
    )
    
    message("Computing TWIMap (4/", nbSteps, "): Flow accumulation...")
    Cache(wbt_d8_flow_accumulation,
          input = plotAndPixeGroupAreaDemPath,
          output = flow_acc_path,
          out_type = "specific contributing area"
    )
    message("Computing TWIMap (5/", nbSteps, "): Final step...")
    # Step 4: Wetness Index
    Cache(wbt_wetness_index,
          sca = flow_acc_path,
          slope = slope_path,
          output = final_twi_path
    )
    
    # Crop the 
    sim$TWIMap <- rast(final_twi_path)
  }
  
  ##############################################################################
  # Download, process and cache CANSIS soil data if it is not supplied
  #
  # https://sis.agr.gc.ca/cansis/nsdb/slc/index.html
  # https://sis.agr.gc.ca/cansis/nsdb/psm/index.html
  # https://www.nature.com/articles/s41597-025-05460-4
  # https://open.canada.ca/data/en/dataset/4d39c9f9-a85c-4bf2-b920-138fdd423384
  # https://agriculture.canada.ca/atlas/data_donnees/griddedSoilsCanada/supportdocument_documentdesupport/en/ISO_19131_Soil_Landscape_Grids_of_Canada_100m_%e2%80%93_Data_Product_Specifications.pdf
  # https://agriculture.canada.ca/atlas/data_donnees/griddedSoilsCanada/data_donnees/raster/Silt/
  #
  ##############################################################################
  mapToProcess <- c("Clay", "Sand", "Silt", "BD") # BD is bulk_density
  baseURL <- "https://sis.agr.gc.ca/cansis/nsdb/psm/"
  nameEnd <- "_X0_5_cm_100m1980-2000v1"
  ext <- ".tif"
  sapply(mapToProcess, function(mapName){
    varMapName <- paste0("WB_VBD_", mapName, "Map") # e.g. WB_VBD_clayMap
    # browser()
    if (!suppliedElsewhere(varMapName, sim)){
      message("Downloading/cropping/reprojecting/resampling and masking ", varMapName, " to sim$pixelGroupMap...") 
      fileName <- paste0(mapName, nameEnd, ext)
      sim[[varMapName]] <- Cache(
        prepInputs,
        url = paste0(baseURL, mapName, "/", fileName),
        targetFile = fileName,
        destinationPath = getPaths()$cache,
        fun = terra::rast,
        cropTo = sim$plotAndPixeGroupArea,
        projectTo = extend(sim$WB_HartJohnstoneForestClassesMap, sim$plotAndPixeGroupArea),
        maskTo = sim$plotAndPixeGroupArea,
        writeTo = paste0(mapName, nameEnd, "_processed", ext),
        method = "bilinear"
      )
    }
  })

  ##############################################################################
  # Fit a model of drainage based on :
  #
  #   1) plot data if they are provided (drainage, standtype, latitude, longitude).
  #      If plot data is not provided, use the provided plot data points (cover 
  #      some part of Canada western boreal). 
  #      If plot data is provided but standtype is not, determine standtype from 
  #      sim$WB_HartJohnstoneForestClassesMap.
  #   2) sim$TWIMap computed for the merged area of the plot data and the 
  #      sim$WB_HartJohnstoneForestClassesMap
  #   4) clay, sand, silt and BD (bulk_density). These datasets contains many NA 
  #      areas. Those areas are considered water or rock.
  #
  ##############################################################################

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
