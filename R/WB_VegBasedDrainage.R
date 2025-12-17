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