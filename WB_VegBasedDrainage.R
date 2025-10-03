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
    expectsInput(objectName = "MRDEMMap", 
                 objectClass = "SpatRast", 
                 desc = "Medium Resolution Digital Elevation Model (MRDEM) for Canada", 
                 sourceURL = NA),
    expectsInput(objectName = "TWIMap", 
                 objectClass = "SpatRast", 
                 desc = "TWI raster computed from Canada DEM", 
                 sourceURL = NA),
    expectsInput(objectName = "DownslopeDistMap", 
                 objectClass = "SpatRast", 
                 desc = "Downslope distance to water map derived from Canada DEM", 
                 sourceURL = NA),
    expectsInput(objectName = "AspectMap", 
                 objectClass = "SpatRast", 
                 desc = "Aspect raster derived from Canada DEM", 
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
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/Silt/Silt_X0_5_cm_100m1980-2000v1.tif"),
    expectsInput(objectName = "WB_VegBasedDrainageModel", 
                 objectClass = "", 
                 desc = "", 
                 sourceURL = NA)
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
  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

### template for your event1
ReComputeDrainageMap <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  
  # Dummy use of input objects to temporarily stop SpaDES warning
  dummy <- sim$MRDEMMap + sim$TWIMap + sim$AspectMap + sim$DownslopeDistMap +
                                 sim$WB_VBD_ClayMap + sim$WB_VBD_SandMap + sim$WB_VBD_SiltMap + sim$WB_VBD_BDMap
  dummy <- sim$plotPoints
  dummy <- sim$WB_VegBasedDrainageModel
  
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
            " pixels for 6 forest classes (\"deci (1)\", \"mixed (2)\", \"conimix (3)\", ",
            "\"jackpine (4)\", \"larch (5)\" and \"spruce (6)\")...")

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
    plotFile <- file.path(getPaths()$modulePath, currentModule(sim), "data/plotData.csv")
    plotDF <- read.csv(plotFile)
    message("plotPoints not supplied. You must provide a CSV table with \"latitude\", ",
            "\"longitude\", \"standtype\" and \"drainage\" following the ",
            "WB_HartJohnstone classification. Loading default plot data points ",
            "as sim$plotPoints (n=", nrow(plotDF), ")...")
    
    # Purge standtype of drainage info
    if ("standtype" %in% names(plotDF)){
      plotDF$standtype <- sub("Poorly-drained ", "", plotDF$standtype)
      plotDF$standtype <- sub("Well-drained ", "", plotDF$standtype)
    }
    
    # Convert character columns to factor
    plotDF[] <- lapply(plotDF, function(col){
      if (is.character(col)) factor(col) else col
    })
    
    # Reassign factors to match WB_HartJohnstone classification used in the WB_WB_HartJohnstoneClasse module
    # labels = c("deci", "mixed", "conimix", "jackpine", "larch", "spruce", "non forested")
    # levels = c(1L, 2L, 3L, 4L, 5L, 6L, 7L)
    new_codes <- c(3L, 1L, 2L, 4L, 5L, 7L, 6L)[as.integer(plotDF$standtype)]
    
    # Assign new codes while keeping levels the same
    plotDF$standtype <- factor(new_codes, levels = 1:7, labels = terra::levels(plotDF$standtype))
      
    plotPoints <- vect(plotDF, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")  # WGS84
    sim$plotPoints <- project(plotPoints, crs(sim$WB_HartJohnstoneForestClassesMap))
    # writeVector(sim$plotPoints, "G:/Home/temp/plotPoints.shp", overwrite=TRUE)
  }
  
  # Compute the joined area covered by the plot data AND the pixelGroupMap raster
  # (actually WB_HartJohnstoneForestClassesMap here) in order to prepInputs() 
  # covariates to this area to fit the model. Once the model is fitted, we crop 
  # the covariate back to the pixelGroupMap area
  
  # Make a buffer around the plot data
  plotPoints100KmBuffers <- aggregate(buffer(sim$plotPoints, width = 100000))  # 1000 m buffer (if CRS is in meters)
  
  # Extract the WB_HartJohnstoneForestClassesMap extent
  WB_HartJohnstoneExtent <- ext(sim$WB_HartJohnstoneForestClassesMap)
  WB_HartJohnstoneExtentPoly <- vect(WB_HartJohnstoneExtent, crs = crs(sim$WB_HartJohnstoneForestClassesMap))
  crs(WB_HartJohnstoneExtentPoly) <- crs(sim$WB_HartJohnstoneForestClassesMap)
  
  # Merge them together
  plotAndPixelGroupArea <- aggregate(rbind(plotPoints100KmBuffers, WB_HartJohnstoneExtentPoly))
  # mapView(plotAndPixelGroupArea)
  # writeVector(plotAndPixelGroupArea, file.path(getPaths()$cache, "plotAndPixelGroupArea.shp"), overwrite = TRUE)
  
  
  ##############################################################################
  # Download and postProcess the Medium Resolution Digital Elevation Model (MRDEM)
  # for Canada if required
  ##############################################################################
  if(!suppliedElsewhere("MRDEMMap", sim) && ( 
     !suppliedElsewhere("TWIMap", sim) || 
     !suppliedElsewhere("DownslopeDistMap", sim) || 
     !suppliedElsewhere("AspectMap", sim))){
    
    message("Downloading/cropping/reprojecting/resampling/masking medium ", 
            "resolution MRDEM dem (80GB) to union of studyarea and a 100km buffer ", 
            "around  buffered plot points...")
    
    
    # Define the path where to save the dem
    plotAndPixelGroupAreaDemPath <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem.tif")
    
    # Download and process the big thing
    # https://open.canada.ca/data/en/dataset/18752265-bda3-498c-a4ba-9dfe68cb98da
    sim$MRDEMMap <- Cache(
      prepInputs,
      url = "https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30/mrdem-30-dtm.tif",
      targetFile = "mrdem-30-dtm.tif",
      destinationPath = getPaths()$cachePath,
      fun = terra::rast,
      cropTo = plotAndPixelGroupArea,
      projectTo = extend(sim$WB_HartJohnstoneForestClassesMap, plotAndPixelGroupArea),
      align_only = TRUE,
      maskTo = plotAndPixelGroupArea,
      method = "bilinear",
      writeTo = plotAndPixelGroupAreaDemPath,
      purge=7
    )
  }
 
  # Define a wrapper function around whitebox functions to make their arguments
  # and result cacheable
  cacheableWhiteboxFct <- function(fun_name, ...){
    dots <- list(...)
    dots["cacheable_input"] <- NULL
    # output <- dots$output
    # dots["output"] <- NULL
    # do.call(fun_name, output = output, dots)
    do.call(fun_name, dots)
    return (rast(dots$output))
  }
  
  ##############################################################################
  # Generate a TWI map from the MRDEM if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("TWIMap", sim)){
    
    nbSteps <- 4
    
    message("Computing TWIMap (1/", nbSteps, ") from MRDEMMap: Filling depressions...")
    dem_filled_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_filled.tif")
    dep <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = sim$MRDEMMap,
      fun_name = "wbt_fill_depressions",
      dem = plotAndPixelGroupAreaDemPath,
      output = dem_filled_path
    )

    message("Computing TWIMap (2/", nbSteps, ") from MRDEMMap: Computing slopes...")
    slope_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_slope.tif")
    slope <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = dep,
      fun_name = "wbt_slope",
      dem = dem_filled_path,
      output = slope_path,
      zfactor = 1
    )
    
    message("Computing TWIMap (3/", nbSteps, ") from MRDEMMap: Flow accumulation...")
    flow_acc_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_flowAccum.tif")
    flow <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = dep,
      fun_name = "wbt_d8_flow_accumulation",
      input = dem_filled_path,
      output = flow_acc_path,
      out_type = "specific contributing area"
    )
    
    message("Computing TWIMap (4/", nbSteps, ") from MRDEMMap: Final step...")
    final_twi_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_TWI.tif")
    sim$TWIMap <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = flow + slope,
      fun_name = "wbt_wetness_index",
      sca = flow_acc_path,
      slope = slope_path,
      output = final_twi_path
    )
  }
  
  ##############################################################################
  # Generate a downslope distance to water map from the MRDEM if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("DownslopeDistMap", sim)){
 #browser()   
    nbSteps <- 4
    
    message("Computing DownslopeDistMap (1/", nbSteps, ") from MRDEMMap: Breaching depressions...")
    dem_breach_filled_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_breachFilledDep.tif")
    breach_dep <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = sim$MRDEMMap,
      fun_name = "wbt_breach_depressions_least_cost",
      dem = plotAndPixelGroupAreaDemPath,
      dist = 3,
      output = dem_breach_filled_path
    )
    
    message("Computing DownslopeDistMap (2/", nbSteps, ") from MRDEMMap: Flow accumulation from breach filled...")
    bf_flow_acc_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_breachFilledFlowAccum.tif")
    flow <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = breach_dep,
      fun_name = "wbt_d8_flow_accumulation",
      input = dem_breach_filled_path,
      output = bf_flow_acc_path,
      out_type = "cells"
    )
    
    message("Computing DownslopeDistMap (3/", nbSteps, ") from MRDEMMap: Extract streams...")
    streams_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_streams.tif")
    streams <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = flow,
      fun_name = "wbt_extract_streams",
      flow_accum = bf_flow_acc_path,
      output = streams_path,
      threshold = 1000
    )

    message("Computing DownslopeDistMap (4/", nbSteps, ") from MRDEMMap: Final step...")
    downslope_dist_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_downslopeDist.tif")
    sim$DownslopeDistMap <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = breach_dep + streams,
      fun_name = "wbt_downslope_distance_to_stream",
      dem = dem_breach_filled_path,
      streams = streams_path,
      output = downslope_dist_path
    )
  }
  
  ##############################################################################
  # Generate an aspect map from the MRDEM if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("AspectMap", sim)){
    #browser()   
    nbSteps <- 4
    
    message("Computing AspectMap (1/", nbSteps, ") from MRDEMMap: Breaching depressions...")
    aspect_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_aspect.tif")
    sim$AspectMap <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = sim$MRDEMMap,
      fun_name = "wbt_aspect",
      dem = plotAndPixelGroupAreaDemPath,
      output = aspect_path
    )
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
  CANSISMapToProcess <- c("Clay", "Sand", "Silt", "BD") # BD is bulk_density
  equivSoilGridsMaps <- c("clay", "sand", "silt", "bdod") # bdod is bulk_density
  baseURL <- "https://sis.agr.gc.ca/cansis/nsdb/psm/"
  nameEnd <- "_X0_5_cm_100m1980-2000v1"
  ext <- ".tif"
  
  # Define a wrapper around gdalTranslate so it's output becomes cacheable
  cacheableGdalTranslateVRT <- function(mapName, destinationPath){
    baseURL="/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/"
    bbox = c(-13597717,7513741,-10879995,5233135)
    igh = '+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs' # proj string for Homolosine projection
    
    targetPath <- gdal_translate(
      paste0(baseURL, mapName, "/", mapName, "_0-5cm_mean.vrt"),
      file.path(destinationPath, paste0(mapName, '_SoilGrids', "_0-5cm_mean.tif")),
      tr = c(250, 250),
      projwin = bbox,
      projwin_srs = igh
    )
    return(rast(targetPath))
  }
  
  sapply(CANSISMapToProcess, function(mapName){
    varMapName <- paste0("WB_VBD_", mapName, "Map") # e.g. WB_VBD_clayMap
    if (!suppliedElsewhere(varMapName, sim)){
      message("Downloading/cropping/reprojecting/resampling and masking ", varMapName, " to sim$pixelGroupMap...") 
      fileName <- paste0(mapName, nameEnd, ext)
      
      # Assign NULL to dynamically assigned maps so SpaDES stop complaining. 
      if (mapName == "Clay") {
        sim$WB_VBD_ClayMap <- NULL
      } else if (mapName == "Sand") {
        sim$WB_VBD_SandMap <- NULL
      } else if (mapName == "Silt") {
        sim$WB_VBD_SiltMap <- NULL
      } else if (mapName == "BD") {
        sim$WB_VBD_BDMap <- NULL
      }
        
      sim[[varMapName]] <- Cache(
        prepInputs,
        url = paste0(baseURL, mapName, "/", fileName),
        targetFile = fileName,
        destinationPath = getPaths()$cache,
        fun = terra::rast,
        cropTo = plotAndPixelGroupArea,
        projectTo = extend(sim$WB_HartJohnstoneForestClassesMap, plotAndPixelGroupArea),
        maskTo = plotAndPixelGroupArea,
        writeTo = paste0("CANSIS_", mapName, nameEnd, "_postProcessed", ext),
        method = "bilinear"
      )

      ##############################################################################
      # Download, process and cache SoilGrids soil data to patch CANSIS one
      # https://www.isric.org/explore/soilgrids
      # https://files.isric.org/soilgrids/latest/data/
      SGMapName <- equivSoilGridsMaps[[which(CANSISMapToProcess == mapName)]]
      message("Downloading SoilGrids ", SGMapName, "...") 
      rast <- Cache(
        cacheableGdalTranslateVRT, 
        SGMapName, 
        destinationPath = getPaths()$cache
      )
      
      message("Cropping/reprojecting/resampling and masking SoilGrids ", SGMapName, "...") 
      patchRast <- Cache(
        postProcess,
        rast,
        cropTo = plotAndPixelGroupArea,
        projectTo = extend(sim$WB_HartJohnstoneForestClassesMap, plotAndPixelGroupArea),
        maskTo = plotAndPixelGroupArea,
        writeTo = file.path(getPaths()$cache, paste0("SoilGrids_", SGMapName, "_0-5cm_mean_postProcessed.tif")),
        method = "bilinear"
      )
      
      message("Patching CANSIS soil ", mapName, " raster NAs with SoilGrids values...")
      sim[[varMapName]] <- Cache(
        cover,
        sim[[varMapName]],
        patchRast / ifelse(SGMapName == "bdod", 100, 10)
      )
    }
  })

  ##############################################################################
  # Fit a model of drainage based on :
  #
  #   - sim$plotPoints - Must contains the drainage, standtype, latitude, 
  #     longitude attributes.
  #     If plot data is not provided, use the default plot data points (covers 
  #     some part of Canada western boreal) to fit a default model. 
  #     If plot data is provided but standtype is not part of it, determine 
  #     standtype from sim$WB_HartJohnstoneForestClassesMap.
  #
  #   - sim$TWIMap computed for the merged area of the plot data and the 
  #     sim$WB_HartJohnstoneForestClassesMap
  #
  #   - Soild data (sim$WB_VBD_ClayMap, sim$WB_VBD_SandMap, sim$WB_VBD_SiltMap 
  #     and sim$WB_VBD_BDMap).
  #     This dataset contains many NA areas. Those areas are considered water or
  #     rock.
  #
  ##############################################################################
  if (!suppliedElsewhere("drainageModel", sim)){

    nbPLotPoints <- nrow(sim$plotPoints)
    message("drainageModel not supplied. Fitting a model using the provided",
            "plot points (n=", nbPLotPoints, "), soil and TWI maps...") 
    covariatesMaps <- c("TWIMap" = "twi",
                        "DownslopeDistMap" = "downslope_dist",
                        "AspectMap" = "aspect", 
                        "WB_VBD_ClayMap" = "clay", 
                        "WB_VBD_SandMap" = "sand",
                        "WB_VBD_SiltMap" = "silt",
                        "WB_VBD_BDMap" = "bulk_den")
    
    # If standtype is not part of the plot data, get the types from sim$WB_HartJohnstoneForestClassesMap
    if (!"standtype" %in% names(sim$plotPoints)){
      covariatesMaps <- c("WB_HartJohnstoneForestClassesMap" = "standtype", covariatesMaps)
    }
  
    # Extract values from covariate maps
    for (i in seq_along(covariatesMaps)) {
      sim$plotPoints <- cbind(sim$plotPoints, extract(sim[[names(covariatesMaps)[i]]], sim$plotPoints)[, -1])  # Remove ID column from extract
      if (! covariatesMaps[i] %in% names(sim$plotPoints)){
        names(sim$plotPoints)[names(sim$plotPoints) == "y"] <- covariatesMaps[i]
      }
    }
    
    # Keep rows where specified columns are not NA
    modelData <- as.data.frame(sim$plotPoints)
    keeps <- complete.cases(modelData[, c(unname(covariatesMaps), "drainage")])
    modelData <- modelData[keeps, ]
    if (nrow(modelData) < nbPLotPoints){
      message("Removed plot points where covariates could not be extracted. n went from ",
              nbPLotPoints, " to ", nrow(modelData), ". To fix this, make sure ",
              "plot points fall into soil and sim$WB_HartJohnstoneForestClassesMap ",
              "extents and into pixels having values (not NA)...")
    }

    modelData <- modelData[, !(names(modelData) %in% c("X", "plot"))]
    
    
    # Split the data frame into training and test data
    inTraining <- createDataPartition(modelData$drainage, p = 0.7, list = FALSE)
    trainSet <- modelData[inTraining, ]
    testSet <- modelData[-inTraining, ]
    message("plot points (n=", nrow(modelData), ") were split between training (n=", 
            nrow(trainSet), ") and test (n=", nrow(testSet), ")...")
    
    fitControl <- trainControl(method = "repeatedcv",
                               number = 5,
                               repeats = 5,
                               search = "random")
    
    message("Fitting the drainage model...")
    modelFit <- Cache(
      train,
      drainage ~ ., # Can consider subsets of covariates here say; e.g. `slope + elevation + age`; `.` considers all the predictors with no interaction terms
      data = trainSet,
      method = "rf", # "rf" = random forest or "xgbTree" = boosted regression trees, # See topepo.github.io/caret for more model tags that can be used here
      trControl = fitControl,
      tuneLength = 10,
      verbose = FALSE
    )
    
    message("Fitting the drainage model. Done...")
    print(modelFit)
    
    # Crop covariate maps back to groupPixelMap now that the model is fitted and 
    # we don't new to extract covariate values at plot points anymore
    for (i in seq_along(covariatesMaps[names(covariatesMaps) != "WB_HartJohnstoneForestClassesMap"])) {
      message("Cropping sim$", names(covariatesMaps)[i], " from the groupPixelMap + pointPlot area to the groupPixelMap area...")
      sim[[names(covariatesMaps)[i]]] <- Cache(
        postProcessTo,
        sim[[names(covariatesMaps)[i]]],
        cropTo = sim$WB_HartJohnstoneForestClassesMap
      )
    }
    
    sim$WB_VegBasedDrainageModel <- modelFit
  }
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
