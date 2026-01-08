defineModule(sim, list(
  name = "WB_VegBasedDrainage",
  description = paste("Construct a model for drainage in western boreal canada using plot data and other covariates and reuse it to predict drainage from Biomass_core simulated maps"),
  keywords = c("drainage", "western boreal"),
  authors =  c(
    person("Pierre", "Racine", email= "pierre.racine@sbf.ulaval.ca", role = "cre"),
    person("Andres", "Caseiro Guilhem", email= "andres.caseiro-guilhem.1@ulaval.ca", role = "aut")
  ),
  childModules = character(0),
  version = list(WB_VegBasedDrainage = "0.0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  # citation = list("citation.bib"),
  # documentation = list("NEWS.md", "README.md", "WB_VegBasedDrainage.Rmd"),
  reqdPkgs = list("data.table", "reproducible", "LandR"),
  loadOrder = list(after = c("WB_HartJohnstoneForestClasses")),
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
                              "WB_HartJohnstoneForestClasses module used as ",
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
                 sourceURL = "https://canelevation-dem.s3.ca-central-1.amazonaws.com/mrdem-30/mrdem-30-dtm.tif"),
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
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/psm/BD/BD_X0_5_cm_100m1980-2000v1.tif"),
    expectsInput(objectName = "WB_VegBasedDrainageModel", 
                 objectClass = "", 
                 desc = "", 
                 sourceURL = NA),
    expectsInput(objectName = "EcoProvincesMap", 
                 objectClass = "SpatVector", 
                 desc = "", 
                 sourceURL = "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/NA_CEC_Eco_Level3.zip")
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
      sim <- scheduleEvent(sim, time(sim), "WB_VegBasedDrainage", "reComputeDrainageMap", 2)
    },
    
    reComputeDrainageMap = {
      sim <- reComputeDrainageMap(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$WB_VegBasedDrainageTimeStep, "WB_VegBasedDrainage", "reComputeDrainageMap")
    },
    warning(noEventWarning(sim))
  )
  return(invisible(sim))
}

reComputeDrainageMap <- function(sim) {
  message("Recomputing sim$WB_VegBasedDrainageMap for ", 
          format(ncell(sim$WB_HartJohnstoneForestClassesMap), scientific = FALSE), " pixels..")
  sim$WB_VegBasedDrainageMap <- computeDrainageMap(
    sim$WB_VegBasedDrainageModel,
    sim$WB_HartJohnstoneForestClassesMap,
    sim$TWIMap,
    sim$DownslopeDistMap,
    sim$AspectMap,
    sim$WB_VBD_ClayMap,
    sim$WB_VBD_SandMap,
    sim$WB_VBD_SiltMap,
    sim$WB_VBD_BDMap,
    sim$EcoProvincesMap
  )
  
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
      ncol = rastWidth,
      nrow = rastWidth,
      crs = "ESRI:102002",
      nbregion = 2000,
      valuevect = 1:6,
      seed = 100,
      userTags = c(userTags, "WB_HartJohnstoneForestClassesMap")
    )
    
    # Convert to factor and add proper labels
    sim$WB_HartJohnstoneForestClassesMap <- terra::as.factor(sim$WB_HartJohnstoneForestClassesMap)
    levels(sim$WB_HartJohnstoneForestClassesMap) <- data.frame(
      value = c(1L, 2L, 3L, 4L, 5L, 6L),
      class = c("deci", "mixed", "conimix", "jackpine", "larch", "spruce")
    )
    names(sim$WB_HartJohnstoneForestClassesMap) <- "standtype"
  }

  if (!is.null(sim$WB_HartJohnstoneForestClassesMap)){
    baseRast <- sim$WB_HartJohnstoneForestClassesMap
  }
  else if (!is.null(sim$pixelGroupMap)){
    baseRast <- sim$pixelGroupMap
  }
  else if (!is.null(sim$rasterToMatch)){
    baseRast <- sim$rasterToMatch
  }
  else {
    stop(paste("At least one of WB_HartJohnstoneForestClassesMap, pixelGroupMap or ",
               "rasterToMatch must be defined in sim before WB_VegBasedDrainage can be initialized..."))
  }
  baseExtent <- ext(baseRast)
  baseCRS <- crs(baseRast)
  baseExtentPoly <- vect(baseExtent, crs = baseCRS)

  ##############################################################################
  # Use our own plotPoints data if none is supplied
  ##############################################################################
    if(!suppliedElsewhere("plotPoints", sim) && !suppliedElsewhere("drainageModel", sim)){
    message("##############################################################################")   
    message("plotPoints not supplied.")   
    message("You must provide a CSV table with \"drainage\", \"latitude\", \"longitude\" ")
    message("and \"standtype\" following the WB_HartJohnstone classification in order to ")
    message("fit the WB_VegBasedDrainageModel.")
    plotFilePath <- file.path(getPaths()$modulePath, currentModule(sim), "data/plotData.csv")
    sim$plotPoints <- getAndcleanPlotData(plotFilePath, baseCRS)
    message("Loaded default plot data points as sim$plotPoints (n=", nrow(sim$plotPoints), ")...")
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
  
  # Merge it with the base extent polygon
  plotAndPixelGroupArea <- aggregate(rbind(plotPoints100KmBuffers, baseExtentPoly))
  
  # mapView(plotAndPixelGroupArea)
  # writeVector(plotAndPixelGroupArea, file.path(getPaths()$cache, "plotAndPixelGroupArea.shp"), overwrite = TRUE)
  
  # Create a raster with the merged area so projectTo does not crop it
  # see https://github.com/PredictiveEcology/reproducible/issues/431
  plotAndPixelGroupAreaRast <- terra::extend(baseRast, plotAndPixelGroupArea)
  
  
  ##############################################################################
  # Define a function to download a portion of a VRT raster
  ##############################################################################
  # cropVRTToDisk <- function(vrtURL, cropTo, destinationPath){
  #   vsicurlURL = paste0("/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=", vrtURL)
  #   vrt <- terra::rast(vsicurlURL)
  #   cropTo <- terra::project(cropTo, crs(vrt))
  #   cropped_rast <- terra::crop(vrt, cropTo)
  #   writeRaster(cropped_rast, destinationPath, overwrite = TRUE)
  #   return(rast(destinationPath))
  # }
  
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
    
    # message("------------------------------------------------------------------------------")   
    # message("Downloading a cropped version of MRDEM...")
    # croppedDEMPath <- file.path(getPaths()$cachePath, "mrdem-30-dtm_cropped.tif")
    # sim$MRDEMMap <- Cache(
    #   cropVRTToDisk,
    #   vrtURL = extractURL("MRDEMMap", sim),
    #   cropTo = plotAndPixelGroupArea,
    #   destinationPath = croppedDEMPath
    # )
    # 
    # message("------------------------------------------------------------------------------")   
    # message("Projecting and masking MRDEM...")
    plotAndPixelGroupAreaDemPath <- file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem.tif")
    # sim$MRDEMMap <- Cache(
    #   postProcess,
    #   from = croppedDEMPath,
    #   projectTo = plotAndPixelGroupAreaRast,
    #   align_only = TRUE,
    #   maskTo = plotAndPixelGroupArea,
    #   writeTo = plotAndPixelGroupAreaDemPath,
    #   userTags = c(userTags, "MRDEMMap")
    # )
      
    # Download and process the big thing
      sim$MRDEMMap <- Cache(
        prepInputs,
        url = extractURL("MRDEMMap", sim),
        targetFile = "mrdem-30-dtm.tif",
        destinationPath = getPaths()$cachePath,
        fun = terra::rast,
        cropTo = plotAndPixelGroupArea,
        projectTo = plotAndPixelGroupAreaRast,
        align_only = TRUE,
        maskTo = plotAndPixelGroupArea,
        method = "bilinear",
        writeTo = plotAndPixelGroupAreaDemPath,
        userTags = c(userTags, "MRDEMMap"),
        purge=7,
        verbose = 3,
        overwrite = TRUE
      )
  }
  
  ##############################################################################
  # Generate a TWI map from the MRDEM if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("TWIMap", sim)){
    message("##############################################################################")   
    message("TWIMap not supplied.")   
    sim$TWIMap <- generateTWIMap(
      dem = sim$MRDEMMap,
      dem_path = plotAndPixelGroupAreaDemPath,
      dem_filled_path = file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_filled.tif"),
      slope_path = file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_slope.tif"),
      flow_acc_path = file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_flowAccum.tif"),
      final_twi_path = file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_TWI.tif"),
      cachePath = getPaths()$cache,
      userTags = userTags
    )
  }
  
  ##############################################################################
  # Generate a downslope distance to water map from the MRDEM if it is not supplied
  ##############################################################################
  if(!suppliedElsewhere("DownslopeDistMap", sim)){
    message("##############################################################################")   
    message("DownslopeDistMap not supplied.")
    sim$DownslopeDistMap <- generateDownslopeDistMap(
      dem = sim$MRDEMMap,
      dem_path = plotAndPixelGroupAreaDemPath,
      dem_breach_filled_path = file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_breachFilledDep.tif"),
      bf_flow_acc_path = file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_breachFilledFlowAccum.tif"),
      streams_path = file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_streams.tif"),
      downslope_dist_path = file.path(getPaths()$cachePath, "plotAndPixelGroupAreaDem_downslopeDist.tif"),
      cachePath = getPaths()$cache,
      userTags = userTags
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
      cachePath = getPaths()$cache,
      userTags = c(userTags, "plotAndPixelGroupAreaDem_aspect.tif")
    )
    names(sim$AspectMap) <- "aspect"
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

  sapply(CANSISMapToProcess, function(mapName){
    varMapName <- paste0("WB_VBD_", mapName, "Map") # e.g. WB_VBD_clayMap
    if (!suppliedElsewhere(varMapName, sim)){
      message("##############################################################################")   
      message(varMapName, " not supplied.")   
      
      # Assign NULL to dynamically assigned maps using their explicit names
      # so SpaDES stop complaining. 
      if (mapName == "Clay") {
        sim$WB_VBD_ClayMap <- NULL
      } else if (mapName == "Sand") {
        sim$WB_VBD_SandMap <- NULL
      } else if (mapName == "Silt") {
        sim$WB_VBD_SiltMap <- NULL
      } else if (mapName == "BD") {
        sim$WB_VBD_BDMap <- NULL
      }

      sim[[varMapName]] <- getAndPatchCANSISSoilMap(
        mapName = mapName,
        plotAndPixelGroupArea = plotAndPixelGroupArea,
        plotAndPixelGroupAreaRast = plotAndPixelGroupAreaRast,
        equivSoilGridsMaps = equivSoilGridsMaps,
        CANSISMapToProcess = CANSISMapToProcess,
        destinationPath = getPaths()$cache,
        cachePath = getPaths()$cache,
        userTags = userTags
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
    ecoProv <- Cache(
      prepInputs,
      url = extractURL("EcoProvincesMap", sim),
      targetFile = "NA_CEC_Eco_Level3.shp",
      destinationPath = getPaths()$cache,
      projectTo = plotAndPixelGroupAreaRast,
      cropTo = plotAndPixelGroupArea,
      writeTo = file.path(getPaths()$cache, "NA_CEC_Eco_Level3_postProcessed.shp"),
      fun = terra::vect,
      userTags = c(userTags, "NA_CEC_Eco_Level3_postProcessed.shp"),
      overwrite = TRUE
    )

    ecoProv <- ecoProv[, c("NA_L3NAME")]
    ecoProv$NA_L3NAME <- as.factor(ecoProv$NA_L3NAME)
    names(ecoProv)[names(ecoProv) == "NA_L3NAME"] <- "ecoprov"
    sim$EcoProvincesMap <- rasterize(ecoProv, plotAndPixelGroupAreaRast, field = "ecoprov")
    
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
    nbPlotPoints <- nrow(sim$plotPoints)
    message("##############################################################################")
    message("WB_VegBasedDrainageModel not supplied. Fitting a model using the ")
    message("plot points provided in the data folder (n=", nbPlotPoints, "), soil (clay, silt, ")
    message("sand and bulked density), aspect, downslope distance to water, ecoprovince and TWI maps...")
    
    # List the covariates from which to extract values
    # element's names are the names of sim maps to extract values from (e.g. sim$TWIMap)
    # element values are the names of the column to create in the sim$plotPoints dataframe (e.g clay)
    covariatesMaps <- c("TWIMap" = "twi",
                        "DownslopeDistMap" = "downslope_dist",
                        "AspectMap" = "aspect", 
                        "WB_VBD_ClayMap" = "clay", 
                        "WB_VBD_SandMap" = "sand",
                        "WB_VBD_SiltMap" = "silt",
                        "WB_VBD_BDMap" = "bd",
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
      writeVector(sim$plotPoints[!keeps, ], file.path(getPaths()$output, "plotPointsRemoved.shp"), overwrite = TRUE)
    }
    
    # Remove the "plot" column
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
    sim$WB_VegBasedDrainageModel <- Cache(
      train,
      drainage ~ .,
      data = trainSet,
      method = "rf", # "rf" = random forest or "xgbTree" = boosted regression trees, 
                     # See topepo.github.io/caret for more model tags that can be used here
      trControl = fitControl,
      tuneLength = 10,
      verbose = FALSE,
      userTags = c(userTags, "WB_VegBasedDrainageModel"),
      overwrite = TRUE
    )
    
    message("Fitting the drainage model. Done...")
    print(sim$WB_VegBasedDrainageModel)

    message("------------------------------------------------------------------------------")   
    print(confusionMatrix(predict(sim$WB_VegBasedDrainageModel, testSet), testSet$drainage))

    message("------------------------------------------------------------------------------")   
    print(varImp(sim$WB_VegBasedDrainageModel))
    
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
        cropTo = baseRast
      )
    }
    
    message("##############################################################################")
  }
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
