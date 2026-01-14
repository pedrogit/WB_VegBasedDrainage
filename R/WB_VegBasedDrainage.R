computeDrainageMap <- function(
    WB_VegBasedDrainageModel,
    HJForestClassesMap,
    TWIMap,
    DownslopeDistMap,
    AspectMap,
    WB_VBD_ClayMap,
    WB_VBD_SandMap,
    WB_VBD_SiltMap,
    WB_VBD_BDMap,
    EcoProvincesMap 
){
  newForestClassesMap <- HJForestClassesMap
  # Merge WB_HartJohnstoneForestClassesMap well drained spruce and poorly drained 
  # spruce together because WB_VegBasedDrainageModel can not be based on any
  # drainage values (only on forest classes values independent of drainage)
  if (any(values(newForestClassesMap) == 7, na.rm = TRUE)){
    message("Some drainage qualified spruce values were found in WB_HartJohnstoneForestClassesMap.")
    message("Merging them together into the spruce category before predicting WB_VegBasedDrainageMap...")
    newForestClassesMap[newForestClassesMap == 7L] <- 6L
    
    levels(newForestClassesMap) <- data.frame(
      value = c(1L, 2L, 3L, 4L, 5L, 6L),
      class = c("deci", "mixed", "conimix", "jackpine", "larch", "spruce")
    )
    names(newForestClassesMap) <- "standtype"
  }
  
  predictors <- c(newForestClassesMap, 
                  TWIMap,
                  DownslopeDistMap,
                  AspectMap,
                  WB_VBD_ClayMap,
                  WB_VBD_SandMap,
                  WB_VBD_SiltMap,
                  WB_VBD_BDMap,
                  EcoProvincesMap
  )
  
  # Recompute the drainage map
  WB_VegBasedDrainageMap <- terra::predict(predictors, WB_VegBasedDrainageModel, na.rm = TRUE)
  
  # Convert to factor and add proper labels
  WB_VegBasedDrainageMap <- terra::as.factor(WB_VegBasedDrainageMap)
  levels(WB_VegBasedDrainageMap) <- data.frame(value = c(1, 2),
                                               class = c("poorly.drained", "well.drained"))
  
  # Assign it a name 
  names(WB_VegBasedDrainageMap) <- "drainage"
  return(WB_VegBasedDrainageMap)
}

##############################################################################
# Function to read and clean the plot data used to fit the WB_VegBasedDrainageModel
##############################################################################
getAndcleanPlotData <- function(plotFile, crs) {
    plotDF <- read.csv(plotFile)

    # Purge standtype of drainage info
    if ("standtype" %in% names(plotDF)){
      plotDF$standtype <- sub("Poorly-drained ", "", plotDF$standtype)
      plotDF$standtype <- sub("Well-drained ", "", plotDF$standtype)
    }

    # Convert character columns to factor and fix values (factor names)
    plotDF[] <- lapply(plotDF, function(col){
    newCol <- col
      if (is.character(newCol)) {
        newCol <- factor(newCol)
        levels(newCol) <- make.names(levels(newCol))
      }
      newCol
    })

    # Reassign factors to match WB_HartJohnstone classification used in the 
    # WB_WB_HartJohnstoneClasse module
    labels = c("deci", "mixed", "conimix", "jackpine", "larch", "spruce", "nonforested")
    # levels = c(1L, 2L, 3L, 4L, 5L, 6L, 7L)
    new_codes <- c(3L, 1L, 2L, 4L, 5L, 7L, 6L)[as.integer(plotDF$standtype)]
    
    # Assign new codes and levels
    # plotDF$standtype <- factor(new_codes, levels = 1:7, labels = terra::levels(plotDF$standtype))
    plotDF$standtype <- factor(new_codes, levels = 1:7, labels = labels)
    
    # Convert the dataframe to a SpatVector object
    plotPoints <- vect(plotDF, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")  # WGS84
    plotPoints <- project(plotPoints, crs)

    return(plotPoints)
  }

##############################################################################
# Define a wrapper function around whitebox functions to make their arguments
# and result cacheable
##############################################################################
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
# Generate a TWI map from a MRDEM map
##############################################################################
generateTWIMap <- function(
  dem,
  dem_path,
  dem_filled_path,
  slope_path,
  flow_acc_path,
  final_twi_path,
  cachePath,
  userTags = NULL
){
  message("Computing TWIMap (1/4) from MRDEMMap: Filling depressions...")
  dep <- Cache(
    cacheableWhiteboxFct,
    cacheable_input = dem,
    fun_name = "wbt_fill_depressions",
    dem = dem_path,
    output = dem_filled_path,
    cachePath = cachePath,
    userTags = c(userTags, "plotAndPixelGroupAreaDem_filled.tif")
  )

  message("------------------------------------------------------------------------------")   
  message("Computing TWIMap (2/4) from MRDEMMap: Computing slopes...")
  slope <- Cache(
    cacheableWhiteboxFct,
    cacheable_input = dep,
    fun_name = "wbt_slope",
    dem = dem_filled_path,
    output = slope_path,
    zfactor = 1,
    cachePath = cachePath,
    userTags = c(userTags, "plotAndPixelGroupAreaDem_slope.tif")
  )
  
  message("------------------------------------------------------------------------------")   
  message("Computing TWIMap (3/4) from MRDEMMap: Flow accumulation...")
  flow <- Cache(
    cacheableWhiteboxFct,
    cacheable_input = dep,
    fun_name = "wbt_d8_flow_accumulation",
    input = dem_filled_path,
    output = flow_acc_path,
    out_type = "specific contributing area",
    cachePath = cachePath,
    userTags = c(userTags, "plotAndPixelGroupAreaDem_flowAccum.tif")
  )
  
  message("------------------------------------------------------------------------------")   
  message("Computing TWIMap (4/4) from MRDEMMap: Final step...")
  # The final result is not cached. It is the responsability of the caller to cache the result.
  TWIMap <- cacheableWhiteboxFct(
    cacheable_input = flow + slope,
    fun_name = "wbt_wetness_index",
    sca = flow_acc_path,
    slope = slope_path,
    output = final_twi_path
  )
  names(TWIMap) <- "twi"
  return (TWIMap)
}

##############################################################################
# Generate a Downslope Distance map from a MRDEM map
##############################################################################
generateDownslopeDistMap <- function(
  dem,
  dem_path,
  dem_breach_filled_path,
  bf_flow_acc_path,
  streams_path,
  downslope_dist_path,
  cachePath,
  userTags = NULL
){
    message("Computing DownslopeDistMap (1/4) from MRDEMMap: Breaching depressions...")
    breach_dep <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = dem,
      fun_name = "wbt_breach_depressions_least_cost",
      dem = dem_path,
      dist = 3,
      output = dem_breach_filled_path,
      cachePath = cachePath,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_breachFilledDep.tif")
    )
    
    message("------------------------------------------------------------------------------")   
    message("Computing DownslopeDistMap (2/4) from MRDEMMap: Flow accumulation from breach filled...")
    flow <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = breach_dep,
      fun_name = "wbt_d8_flow_accumulation",
      input = dem_breach_filled_path,
      output = bf_flow_acc_path,
      out_type = "cells",
      cachePath = cachePath,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_breachFilledFlowAccum.tif")
    )
    
    message("------------------------------------------------------------------------------")   
    message("Computing DownslopeDistMap (3/4) from MRDEMMap: Extract streams...")
    streams <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = flow,
      fun_name = "wbt_extract_streams",
      flow_accum = bf_flow_acc_path,
      output = streams_path,
      threshold = 1000,
      cachePath = cachePath,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_streams.tif")
    )

    message("------------------------------------------------------------------------------")   
    message("Computing DownslopeDistMap (4/4) from MRDEMMap: Final step...")
    # The final result is not cached. It is the responsability of the caller to cache the result.
    downslopeDistMap <- cacheableWhiteboxFct(
      cacheable_input = breach_dep + streams,
      fun_name = "wbt_downslope_distance_to_stream",
      dem = dem_breach_filled_path,
      streams = streams_path,
      output = downslope_dist_path
    )
    names(downslopeDistMap) <- "downslope_dist"
    return (downslopeDistMap)
}

##############################################################################
# Define a wrapper around gdalTranslate for this specific dataset so it's 
# output becomes cacheable
##############################################################################
cacheableGdalTranslateVRT <- function(mapName, destinationPath){
  baseURL="/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/"
  vrt_rast <- rast(paste0(baseURL, mapName, "/", mapName, "_0-5cm_mean.vrt"))
  cropped_rast <- crop(vrt_rast, ext(-13597717, -10879995, 5233135, 7513741))
  
  processedFilePath <- file.path(destinationPath, paste0('SoilGrids_', mapName, "_0-5cm_mean.tif"))
  writeRaster(cropped_rast, processedFilePath, overwrite = TRUE)
  
  return(rast(processedFilePath))
}

##############################################################################
# Download and patch CANSIS soil maps with SoilGrids data
##############################################################################
getAndPatchCANSISSoilMap <- function(
  mapName,
  plotAndPixelGroupArea,
  plotAndPixelGroupAreaRast,
  CANSISMapToProcess,
  equivSoilGridsMaps,
  destinationPath,
  cachePath,
  userTags = NULL
){
  varMapName <- paste0("WB_VBD_", mapName, "Map") # e.g. WB_VBD_clayMap
  message("Downloading/cropping/reprojecting/resampling and masking ", varMapName, "...")

  baseURL = "https://sis.agr.gc.ca/cansis/nsdb/psm"
  nameEnd <- "_X0_5_cm_100m1980-2000v1"
  ext <- ".tif"

  fileName <- paste0(mapName, nameEnd, ext)
  rast <- Cache(
    prepInputs,
    url = paste0(baseURL, "/", mapName, "/", fileName),
    targetFile = fileName,
    destinationPath = destinationPath,
    fun = terra::rast,
    cropTo = plotAndPixelGroupArea,
    projectTo = plotAndPixelGroupAreaRast,
    maskTo = plotAndPixelGroupArea,
    writeTo = paste0("CANSIS_", mapName, nameEnd, "_postProcessed", ext),
    method = "bilinear",
    cachePath = cachePath,
    userTags = c(userTags, paste0("CANSIS_", mapName, nameEnd, "_postProcessed", ext)),
    overwrite = TRUE
  )
  
  # Ensure the raster variable has the right name
  names(rast) <- tolower(mapName)

  #-------------------------------------------------------------------------------
  # Download, process and cache SoilGrids soil data to patch CANSIS one
  # https://www.isric.org/explore/soilgrids
  # https://files.isric.org/soilgrids/latest/data/
  SGMapName <- equivSoilGridsMaps[[which(CANSISMapToProcess == mapName)]]
  message("------------------------------------------------------------------------------")   
  message("Downloading SoilGrids ", SGMapName, "...")
  sgrast <- Cache(
    cacheableGdalTranslateVRT,
    SGMapName,
    destinationPath = destinationPath,
    cachePath = cachePath,
    userTags = c(userTags, paste0('SoilGrids_0-5cm_mean_', SGMapName, ".tif"))
  )

  message("------------------------------------------------------------------------------")   
  message("Cropping/reprojecting/resampling and masking SoilGrids ", SGMapName, "...") 
  patchRast <- Cache(
    postProcess,
    sgrast,
    cropTo = plotAndPixelGroupArea,
    projectTo = plotAndPixelGroupAreaRast,
    maskTo = plotAndPixelGroupArea,
    writeTo = file.path(destinationPath, paste0("SoilGrids_", SGMapName, "_0-5cm_mean_postProcessed.tif")),
    method = "bilinear",
    cachePath = cachePath,
    userTags = c(userTags, paste0('SoilGrids_0-5cm_mean_', SGMapName, "_postProcessed.tif")),
    overwrite = TRUE
  )
  
  message("------------------------------------------------------------------------------")   
  message("Patching CANSIS soil ", mapName, " raster NAs with SoilGrids values...")
  # The final result is not cached. It is the responsability of the caller to cache the result.
  rast <- cover(
    rast,
    patchRast / ifelse(SGMapName == "bdod", 100, 10)
  )

  return(rast)
}

##############################################################################
# Download and patch CANSIS soil maps with SoilGrids data
##############################################################################
fit_WB_VegBasedDrainageModel <- function(
  plotPoints,
  covariateMapList,
  pixelDist = 0
){
  nbPlotPoints <- nrow(plotPoints)
  message("Fitting a model using the plot points provided in the data folder (n=", nbPlotPoints, "), soil (clay, silt, ")
  message("sand and bulked density), aspect, downslope distance to water, ecoprovince and TWI maps...")
  message("------------------------------------------------------------------------------")   
  # Extract values from covariate maps
  for (i in seq_along(covariateMapList)) {
    # Extract the values only if the covariate is not NULL
    if (!is.null(covariateMapList[i])){
      message("Extracting values from ", names(covariateMapList)[i], "...")
      message("  searching at a max distance of ", pixelDist, " pixels for \"with value\" pixels...")
      # For SpatRaster we can:
      # - bind the produced column directly,
      # - add a radius to search for the neared pixel value if the point fall into a NA pixel,
      # - prevent the ID column to be added.
      plotPoints <- terra::extract(
        covariateMapList[[i]], 
        plotPoints, 
        bind = TRUE, 
        method = "bilinear", 
        search_radius = pixelDist * mean(res(covariateMapList[[i]])), 
        ID = FALSE
      )
      
      # Remove the last useless columns
      plotPoints <- plotPoints[, -((ncol(plotPoints) - 1):ncol(plotPoints))]

      # Rename the extracted column
      if (! names(covariateMapList)[i] %in% names(plotPoints)){
        names(plotPoints)[names(plotPoints) == names(covariateMapList)[i]] <- names(covariateMapList)[i]
      }
    }
    else {
      message("WARNING: For some reason,", names(covariateMapList)[i], " does ", 
              "not exist... It will not be taken into account when fittng the model...")
    }
  }
  message("------------------------------------------------------------------------------")   

  # Keep rows where specified columns are not NA
  modelData <- as.data.frame(plotPoints)
  varToKeep <- names(covariateMapList)
  if ("drainage" %in% names(plotPoints)){
    varToKeep <- c(varToKeep, "drainage")
  }
  keeps <- complete.cases(modelData[, varToKeep])
  modelData <- modelData[keeps, ]
  
  if (nrow(modelData) < nbPlotPoints){
    message("------------------------------------------------------------------------------")   
    message("WARNING: Removed ", nbPlotPoints - nrow(modelData), " plot points where ")
    message("covariates could not be extracted. n went from ", nbPlotPoints)
    message(" to ", nrow(modelData), ". To fix this, make sure plot points fall ")
    message("into soil and sim$WB_HartJohnstoneForestClassesMap combined extents ")
    message("and into pixels having values (not NA) or increase the value of the ")
    message("\"searchDistInPixelNb\" parameter...")
    message("Point removed (", paste(as.character(modelData[!keeps, "plot"]), collapse = ", "), 
            ") were saved to a shapefile in the output folder (plotPointsEliminated.shp)")
    message("------------------------------------------------------------------------------")   
    writeVector(plotPoints[!keeps, ], file.path(getPaths()$output, "plotPointsRemoved.shp"), overwrite = TRUE)
  }
  
  # Remove the "X" and the "plot" columns
  modelData <- modelData[, !(names(modelData) %in% c("X", "plot"))]
  
  # Convert the standtype from factor to numeric so the model recognize the standtype raster values
  #modelData$standtype <- as.integer(modelData$standtype)
  
  # Split the plot data into training and test data
  set.seed(1990)
  inTraining <- createDataPartition(modelData$drainage, p = 0.7, list = FALSE)
  trainSet <- modelData[inTraining, ]
  testSet <- modelData[-inTraining, ]
  message("Plot points (n=", nrow(modelData), ") were split between training (n=", 
          nrow(trainSet), ") and test (n=", nrow(testSet), ")...")
  
  # Parametrize the training algorithm
  set.seed(1990)
  fitControl <- trainControl(
    method = "repeatedcv",
    repeats = 5,
    sampling = "down",
    summaryFunction = twoClassSummary,
    classProbs = TRUE
  )

  # Fit the model using all the covariate present in the modeldata frame
  message("Fitting the drainage model...")
  drainageModel <- Cache(
    train,
    drainage ~ .,
    data = trainSet,
    method = "rf", # "rf" = random forest or "xgbTree" = boosted regression trees, 
                    # See topepo.github.io/caret for more model tags that can be used here
    trControl = fitControl,
    tuneLength = 10,
    verbose = FALSE,
    userTags = "drainageModel",
    overwrite = TRUE
  )
  
  message("Fitting the drainage model. Done...")
  print(drainageModel)

  message("------------------------------------------------------------------------------")   
  print(confusionMatrix(predict(drainageModel, testSet), testSet$drainage))

  message("------------------------------------------------------------------------------")   
  print(varImp(drainageModel))
  return(drainageModel)
}

