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
                    "Simulation time at which the drainage map is regenerated."),
    defineParameter("searchDistInPixelNb", "numeric", 1, NA, NA,
                    "Distance, in number of pixel, to search for \"with value\" covariate values when plot points fall into NS values.")
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
  userTags <- c(currentModule(sim), "function:.inputObjects") 
  ##############################################################################
  # Generate a fake WB_HartJohnstoneForestClassesMap if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("WB_HartJohnstoneForestClassesMap", sim)){
    rastWidth <- 1000
    message("##############################################################################")   
    message("WB_HartJohnstoneForestClassesMap not supplied.")   
    message("Please couple with the WB_HartJohnstoneForestClasses module. ")   
    message("Creating random map ", rastWidth, " pixels by ", rastWidth, " pixels for 6 forest ")   
    message("classes (\"deci (1)\", \"mixed (2)\", \"conimix (3)\", \"jackpine (4)\", ")   
    message("\"larch (5)\" and \"spruce (6)\")...")

    sim$WB_HartJohnstoneForestClassesMap <- Cache(
      getRandomCategoricalMap,
      origin = c(-667296, 1758502),
      width = rastWidth,
      crs = "ESRI:102002",
      nbregion = 2000,
      valuevect = 1:6,
      seed = 100,
      userTags = c(userTags, "WB_HartJohnstoneForestClassesMap")
    )
    # mapView(sim$WB_HartJohnstoneForestClassesMap)
  }
  
  ##############################################################################
  # Use our own plotPoints data if none is supplied
  ##############################################################################
  if(!suppliedElsewhere("plotPoints", sim) && !suppliedElsewhere("drainageModel", sim)){
    plotFile <- file.path(getPaths()$modulePath, currentModule(sim), "data/plotData.csv")
    plotDF <- read.csv(plotFile)
    message("##############################################################################")   
    message("plotPoints not supplied.")   
    message("You must provide a CSV table with \"drainage\", \"latitude\", \"longitude\" ")
    message("and \"standtype\" following the WB_HartJohnstone classification.")
    message(" Loading default plot data points as sim$plotPoints (n=", nrow(plotDF), ")...")
    
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
    # labels = c("deci", "mixed", "conimix", "jackpine", "larch", "spruce", "non forested")
    # levels = c(1L, 2L, 3L, 4L, 5L, 6L, 7L)
    new_codes <- c(3L, 1L, 2L, 4L, 5L, 7L, 6L)[as.integer(plotDF$standtype)]
    
    # Assign new codes while keeping levels the same
    plotDF$standtype <- factor(new_codes, levels = 1:7, labels = terra::levels(plotDF$standtype))
    
    # Convert the dataframe to a SpatVector object
    plotPoints <- vect(plotDF, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")  # WGS84
    sim$plotPoints <- project(plotPoints, crs(sim$WB_HartJohnstoneForestClassesMap))
    # writeVector(sim$plotPoints, "G:/Home/temp/plotPoints.shp", overwrite=TRUE)
  }
  
  ##############################################################################
  # Compute the joined area covered by the plot data AND the pixelGroupMap raster
  # (actually WB_HartJohnstoneForestClassesMap here) in order to prepInputs() 
  # covariates to this area to fit the model. Once the model is fitted, we crop 
  # the covariate back to the pixelGroupMap area
  ##############################################################################
  
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
  
  # We create a raster with the merged area so projectTo does not crop it
  # see https://github.com/PredictiveEcology/reproducible/issues/431
  plotAndPixelGroupAreaRast <- extend(sim$WB_HartJohnstoneForestClassesMap, plotAndPixelGroupArea)

  ##############################################################################
  # Download and postProcess the Medium Resolution Digital Elevation Model (MRDEM)
  # for Canada if required
  ##############################################################################
  if(!suppliedElsewhere("MRDEMMap", sim) && ( 
     !suppliedElsewhere("TWIMap", sim) || 
     !suppliedElsewhere("DownslopeDistMap", sim) || 
     !suppliedElsewhere("AspectMap", sim))){
    
    message("##############################################################################")   
    message("MRDEMMap not supplied.")   
    message("Downloading/cropping/reprojecting/resampling/masking medium resolution ")
    message("MRDEM dem (80GB) to union of studyarea and a 100km buffer around buffered plot points...")
    
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
      projectTo = plotAndPixelGroupAreaRast,
      align_only = TRUE,
      maskTo = plotAndPixelGroupArea,
      method = "bilinear",
      writeTo = plotAndPixelGroupAreaDemPath,
      purge=7,
      userTags = c(userTags, "MRDEMMap")
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
    message("##############################################################################")   
    message("TWIMap not supplied.")   
    message("Computing TWIMap (1/4) from MRDEMMap: Filling depressions...")
    dem_filled_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_filled.tif")
    dep <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = sim$MRDEMMap,
      fun_name = "wbt_fill_depressions",
      dem = plotAndPixelGroupAreaDemPath,
      output = dem_filled_path,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_filled.tif")
    )

    message("------------------------------------------------------------------------------")   
    message("Computing TWIMap (2/4) from MRDEMMap: Computing slopes...")
    slope_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_slope.tif")
    slope <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = dep,
      fun_name = "wbt_slope",
      dem = dem_filled_path,
      output = slope_path,
      zfactor = 1,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_slope.tif")
    )
    
    message("------------------------------------------------------------------------------")   
    message("Computing TWIMap (3/4) from MRDEMMap: Flow accumulation...")
    flow_acc_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_flowAccum.tif")
    flow <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = dep,
      fun_name = "wbt_d8_flow_accumulation",
      input = dem_filled_path,
      output = flow_acc_path,
      out_type = "specific contributing area",
      userTags = c(userTags, "plotAndPixelGroupAreaDem_flowAccum.tif")
    )
    
    message("------------------------------------------------------------------------------")   
    message("Computing TWIMap (4/4) from MRDEMMap: Final step...")
    final_twi_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_TWI.tif")
    sim$TWIMap <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = flow + slope,
      fun_name = "wbt_wetness_index",
      sca = flow_acc_path,
      slope = slope_path,
      output = final_twi_path,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_TWI.tif")
    )
  }
  
  ##############################################################################
  # Generate a downslope distance to water map from the MRDEM if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("DownslopeDistMap", sim)){
    message("##############################################################################")   
    message("DownslopeDistMap not supplied.")   
    message("Computing DownslopeDistMap (1/4) from MRDEMMap: Breaching depressions...")
    dem_breach_filled_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_breachFilledDep.tif")
    breach_dep <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = sim$MRDEMMap,
      fun_name = "wbt_breach_depressions_least_cost",
      dem = plotAndPixelGroupAreaDemPath,
      dist = 3,
      output = dem_breach_filled_path,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_breachFilledDep.tif")
    )
    
    message("------------------------------------------------------------------------------")   
    message("Computing DownslopeDistMap (2/4) from MRDEMMap: Flow accumulation from breach filled...")
    bf_flow_acc_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_breachFilledFlowAccum.tif")
    flow <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = breach_dep,
      fun_name = "wbt_d8_flow_accumulation",
      input = dem_breach_filled_path,
      output = bf_flow_acc_path,
      out_type = "cells",
      userTags = c(userTags, "plotAndPixelGroupAreaDem_breachFilledFlowAccum.tif")
    )
    
    message("------------------------------------------------------------------------------")   
    message("Computing DownslopeDistMap (3/4) from MRDEMMap: Extract streams...")
    streams_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_streams.tif")
    streams <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = flow,
      fun_name = "wbt_extract_streams",
      flow_accum = bf_flow_acc_path,
      output = streams_path,
      threshold = 1000,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_streams.tif")
    )

    message("------------------------------------------------------------------------------")   
    message("Computing DownslopeDistMap (4/4) from MRDEMMap: Final step...")
    downslope_dist_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_downslopeDist.tif")
    sim$DownslopeDistMap <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = breach_dep + streams,
      fun_name = "wbt_downslope_distance_to_stream",
      dem = dem_breach_filled_path,
      streams = streams_path,
      output = downslope_dist_path,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_downslopeDist.tif")
    )
  }
  
  ##############################################################################
  # Generate an aspect map from the MRDEM if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("AspectMap", sim)){
    message("##############################################################################")   
    message("AspectMap not supplied.")   
    message("Computing AspectMap from MRDEMMap...")
    aspect_path <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_aspect.tif")
    sim$AspectMap <- Cache(
      cacheableWhiteboxFct,
      cacheable_input = sim$MRDEMMap,
      fun_name = "wbt_aspect",
      dem = plotAndPixelGroupAreaDemPath,
      output = aspect_path,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_aspect.tif")
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
  
  # Define a wrapper around gdalTranslate for this specific dataset so it's 
  # output becomes cacheable
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
      message("##############################################################################")   
      message(varMapName, " not supplied.")   
      message("Downloading/cropping/reprojecting/resampling and masking ", 
              varMapName, " to sim$pixelGroupMap...") 
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
        projectTo = plotAndPixelGroupAreaRast,
        maskTo = plotAndPixelGroupArea,
        writeTo = paste0("CANSIS_", mapName, nameEnd, "_postProcessed", ext),
        method = "bilinear",
        userTags = c(userTags, paste0("CANSIS_", mapName, nameEnd, "_postProcessed", ext))
      )
      
      # Ensure the raster variable has the right name
      names(sim[[varMapName]]) <- mapName

      ##############################################################################
      # Download, process and cache SoilGrids soil data to patch CANSIS one
      # https://www.isric.org/explore/soilgrids
      # https://files.isric.org/soilgrids/latest/data/
      SGMapName <- equivSoilGridsMaps[[which(CANSISMapToProcess == mapName)]]
      message("------------------------------------------------------------------------------")   
      message("Downloading SoilGrids ", SGMapName, "...") 
      rast <- Cache(
        cacheableGdalTranslateVRT, 
        SGMapName, 
        destinationPath = getPaths()$cache,
        userTags = c(userTags, paste0('SoilGrids_0-5cm_mean_', SGMapName, ".tif"))
      )
      
      message("------------------------------------------------------------------------------")   
      message("Cropping/reprojecting/resampling and masking SoilGrids ", SGMapName, "...") 
      patchRast <- Cache(
        postProcess,
        rast,
        cropTo = plotAndPixelGroupArea,
        projectTo = plotAndPixelGroupAreaRast,
        maskTo = plotAndPixelGroupArea,
        writeTo = file.path(getPaths()$cache, paste0("SoilGrids_", SGMapName, "_0-5cm_mean_postProcessed.tif")),
        method = "bilinear",
        userTags = c(userTags, paste0('SoilGrids_0-5cm_mean_', SGMapName, "_postProcessed.tif"))
      )
      
      message("------------------------------------------------------------------------------")   
      message("Patching CANSIS soil ", mapName, " raster NAs with SoilGrids values...")
      sim[[varMapName]] <- Cache(
        cover,
        sim[[varMapName]],
        patchRast / ifelse(SGMapName == "bdod", 100, 10),
        userTags = c(userTags, paste0("CANSIS_", mapName, nameEnd, "_patched", ext))
      )
    }
  })

  ##############################################################################
  # Download, process and cache ecoprovince if it is not supplied
  # https://www.epa.gov/eco-research/ecoregions-north-america
  ##############################################################################
  if(!suppliedElsewhere("EcoProvincesMap", sim)){
    message("##############################################################################")   
    message("EcoProvincesMap not supplied.")   
    message("Downloading and projecting raster NAs with SoilGrids values...")
    sim$EcoProvincesMap <- Cache(
      prepInputs,
      url = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip",
      targetFile = "NA_CEC_Eco_Level3.shp",
      destinationPath = getPaths()$cache,
      projectTo = plotAndPixelGroupAreaRast,
      cropTo = plotAndPixelGroupArea,
      writeTo = file.path(getPaths()$cache, "NA_CEC_Eco_Level3_postProcessed.shp"),
      fun = terra::vect,
      userTags = c(userTags, "NA_CEC_Eco_Level3_postProcessed.shp")
    )
    sim$EcoProvincesMap <- sim$EcoProvincesMap[, c("NA_L3NAME")]
    sim$EcoProvincesMap$NA_L3NAME <- as.factor(sim$EcoProvincesMap$NA_L3NAME)
    
    # Alternative dataset from Canada Open
    # # https://open.canada.ca/data/en/dataset/98fa7335-fbfe-4289-9a0e-d6bf3874b424
    # sim$EcoProvincesMap <- Cache(
    #   prepInputs,
    #   url = "https://agriculture.canada.ca/atlas/data_donnees/nationalEcologicalFramework/data_donnees/geoJSON/ep/nef_ca_ter_ecoprovince_v2_2.geojson",
    #   targetFile = "nef_ca_ter_ecoprovince_v2_2.geojson",
    #   destinationPath = getPaths()$cache,
    #   projectTo = plotAndPixelGroupAreaRast,
    #   cropTo = plotAndPixelGroupArea,
    #   writeTo = file.path(getPaths()$cache, paste0("nef_ca_ter_ecoprovince_v2_2_postProcessed.geojson")),
    #   fun = terra::vect
    # )
    # sim$EcoProvincesMap <- sim$EcoProvincesMap[, c("ECOPROVINCE_NAME_EN")]
    # sim$EcoProvincesMap$ECOPROVINCE_NAME_EN <- as.factor(sim$EcoProvincesMap$ECOPROVINCE_NAME_EN)
  }

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
  if (!suppliedElsewhere("WB_VegBasedDrainageModel", sim)){
    nbPLotPoints <- nrow(sim$plotPoints)
    message("##############################################################################")
    message("WB_VegBasedDrainageModel not supplied. Fitting a model using the ")
    message("provided plot points (n=", nbPLotPoints, "), soil (caly, silt, sand and bulked ")
    message("density), aspect, downslope distance to water, ecoprovince and TWI maps...")
    
    # List the covariates from which tio extract values
    # element's names are the names of sim maps to extract values from (e.g. sim$TWIMap)
    # element values are the names of the column to create in the sim$plotPoints dataframe (e.g clay)
    covariatesMaps <- c("TWIMap" = "twi",
                        "DownslopeDistMap" = "downslope_dist",
                        "AspectMap" = "aspect", 
                        "WB_VBD_ClayMap" = "clay", 
                        "WB_VBD_SandMap" = "sand",
                        "WB_VBD_SiltMap" = "silt",
                        "WB_VBD_BDMap" = "bulk_den",
                        "EcoProvincesMap" = "ecoprov")

    # If standtype is not part of the plot data, get the types from sim$WB_HartJohnstoneForestClassesMap
    if (!"standtype" %in% names(sim$plotPoints)){
      covariatesMaps <- c("WB_HartJohnstoneForestClassesMap" = "standtype", covariatesMaps)
    }
    message("------------------------------------------------------------------------------")   
    
    # Extract values from covariate maps
    for (i in seq_along(covariatesMaps)) {
      # Extract the values only if the covariate exists and it's not NULL
      if (names(covariatesMaps)[i] %in% names(sim) && !is.null(sim[[names(covariatesMaps)[i]]])){
        message("Extracting values from ", names(covariatesMaps)[i], 
                " (", class(sim[[names(covariatesMaps)[i]]]), ")...")
        # Extract properly according to the type of the covariate (SpatRaster or SpatVector)
        if (class(sim[[names(covariatesMaps)[i]]]) == "SpatRaster") {
          message("  searching at a max distance of ", P(sim)$searchDistInPixelNb, " pixels for \"with value\" pixels...")
          # For SpatRaster we can:
          # - bind the produced column directly,
          # - add a radius to search for the neared pixel value if the point fall into a NA pixel,
          # - prevent the ID column to be added.
          sim$plotPoints <- terra::extract(
            sim[[names(covariatesMaps)[i]]], 
            sim$plotPoints, 
            bind = TRUE, 
            method = "bilinear", 
            search_radius = P(sim)$searchDistInPixelNb * mean(res(sim[[names(covariatesMaps)[i]]])), 
            ID=FALSE
          )
          # Remove the last useless columns
          sim$plotPoints <- sim$plotPoints[, -((ncol(sim$plotPoints) - 1):ncol(sim$plotPoints))]
        }
        else {
          # For SpatVector we have to bind manually (as the bind option does not exist for this signature)...
          sim$plotPoints <- cbind(sim$plotPoints, extract(sim[[names(covariatesMaps)[i]]], sim$plotPoints))
          # ...and remove the id column manually as well
          sim$plotPoints <- sim$plotPoints[, !names(sim$plotPoints) %in% "id.y"]
        }

        # Rename the extracted column
        if (! covariatesMaps[i] %in% names(sim$plotPoints)){
          names(sim$plotPoints)[names(sim$plotPoints) == names(sim[[names(covariatesMaps)[i]]])] <- covariatesMaps[i]
        }
      }
      else {
        message("WARNING: For some reason,", names(covariatesMaps)[i], " does ", 
                "not exist... It will not be taken into account when fittng the model...")
      }
    }
    message("------------------------------------------------------------------------------")   
    
    # Keep rows where specified columns are not NA
    modelData <- as.data.frame(sim$plotPoints)
    keeps <- complete.cases(modelData[, c(unname(covariatesMaps), "drainage")])
    modelData <- modelData[keeps, ]
    
    if (nrow(modelData) < nbPLotPoints){
      message("------------------------------------------------------------------------------")   
      message("WARNING: Removed ", nbPLotPoints - nrow(modelData), " plot points where ")
      message("covariates could not be extracted. n went from ", nbPLotPoints)
      message(" to ", nrow(modelData), ". To fix this, make sure plot points fall ")
      message("into soil and sim$WB_HartJohnstoneForestClassesMap combined extents ")
      message("and into pixels having values (not NA) or increase the value of the ")
      message("\"searchDistInPixelNb\" parameter...")
      message("Point removed (", paste(as.character(modelData[!keeps, "plot"]), collapse = ", "), 
              ") were saved to a shapefile in the output folder (plotPointsEliminated.shp)")
      message("------------------------------------------------------------------------------")   
      writeVector(sim$plotPoints[!keeps, ], file.path(getPaths()$output, "plotPointsRemoved.shp"), overwrite = TRUE)
    }
    
    # Remove the "plot" column
    modelData <- modelData[, !(names(modelData) %in% c("X", "plot"))]
    
    # Split the plot data into training and test data
    set.seed(1990)
    inTraining <- createDataPartition(modelData$drainage, p = 0.7, list = FALSE)
    trainSet <- modelData[inTraining, ]
    testSet <- modelData[-inTraining, ]
    message("Plot points (n=", nrow(modelData), ") were split between training (n=", 
            nrow(trainSet), ") and test (n=", nrow(testSet), ")...")
    
    # Parametrize the training algorithm
    fitControl <- trainControl(
      method = "repeatedcv",
      repeats = 5,
      sampling = "down",
      summaryFunction = twoClassSummary,
      classProbs = TRUE
    )

    # Fit the model using all the covariate present in the modeldata frame
    message("Fitting the drainage model...")
    modelFit <- Cache(
      train,
      drainage ~ .,
      data = trainSet,
      method = "rf", # "rf" = random forest or "xgbTree" = boosted regression trees, 
                     # See topepo.github.io/caret for more model tags that can be used here
      trControl = fitControl,
      tuneLength = 10,
      verbose = FALSE,
      userTags = c(userTags, "WB_VegBasedDrainageModel")
    )
    
    message("Fitting the drainage model. Done...")
    print(modelFit)

    message("------------------------------------------------------------------------------")   
    print(confusionMatrix(predict(modelFit, testSet), testSet$drainage))

    message("------------------------------------------------------------------------------")   
    print(varImp(modelFit))
    
    # Crop covariate maps back to groupPixelMap now that the model is fitted and 
    # we don't new to extract covariate values at plot points anymore.
    # From now on the module will, at each iteration step, use the model to predict 
    # drainage from the covariate maps
    for (i in seq_along(covariatesMaps[names(covariatesMaps) != "WB_HartJohnstoneForestClassesMap"])) {
      message("------------------------------------------------------------------------------")   
      message("Cropping sim$", names(covariatesMaps)[i], " from the groupPixelMap + pointPlot area to the groupPixelMap area...")
      sim[[names(covariatesMaps)[i]]] <- Cache(
        postProcessTo,
        sim[[names(covariatesMaps)[i]]],
        cropTo = sim$WB_HartJohnstoneForestClassesMap
      )
    }
    
    message("##############################################################################")
    
    sim$WB_VegBasedDrainageModel <- modelFit
  }
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
