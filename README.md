# WB_VegBasedDrainage

### Overview

<!-- Pierre: This two paragraphshould be included, with obvious mods, in each of the four module .mds -->
WB_VegBasedDrainage is a [SpaDES](https://spades.predictiveecology.org/) module 
complementing the [LandR](https://landr-manual.predictiveecology.org/) ecosystem of modules for forest biomass and succession 
simulation. It is part of an ensemble of modules that provide to LandR the statistical prediction of terrestrial lichen biomass from stand type, time-since fire, and terrestrial ecoprovince. The modules are an implementatiom of [Greuel and Degré-Timmons et al (2021)](https://esajournals-onlinelibrary-wiley-com.acces.bibl.ulaval.ca/doi/full/10.1002/ecs2.3481), developed to support lichen biomass modelling for woodland caribou conservation in the Northwest Territories. The geographical zone wherein the model may reasonably be applied should be assessed from Figure 2 of the cited paper. 


The components of the module ensemble are:

- [WB_HartJohnstoneForestClasses](https://github.com/pedrogit/WB_HartJohnstoneForestClasses) - Generates a map classifying LandR forested pixels to 6 (or 7) classes.
- [WB_VegBasedDrainage](https://github.com/pedrogit/WB_VegBasedDrainage) - This module. Generates a drainage class map.
- [WB_NonForestedVegClasses](https://github.com/pedrogit/WB_NonForestedVegClasses) - Generates a map of land cover classes for areas LandR considers to be non-forested.
- [WB_LichenBiomass](https://github.com/pedrogit/WB_LichenBiomass) - Generates a wall-to-wall map of predicted lichen biomass density for forested and non-forested pixels.

The modules are based on extenisve samples of field data collected as described in [Greuel and Degré-Timmons et al (2021)](https://esajournals-onlinelibrary-wiley-com.acces.bibl.ulaval.ca/doi/full/10.1002/ecs2.3481), Casheiro-Guilhem et. al (in prep.) and foundational papers by Harte and Johstone (citations to be added).

At each simulation step, WB_VegBasedDrainage creates a drainage raster map of 
well-drained and poorly-drained pixels for the forested areas covered by LandR 
only. It produces a raster with two classes:

| Class Code | Description |
|-----------|-------------|
| 1 | Poorly-drained |
| 2 | Well-drained |

WB_HartJohnstoneForestClasses (and hence LandR) pixels set to NA are also set to 
NA by WB_VegBasedDrainage.

WB_VegBasedDrainage determines drainage using a machine learning algorithm 
trained on field-observed drainage conditions, against a set of static (TWI, 
downslope distance, aspect, ecoprovince and soil data (clay, sand", silt, 
bulk_density)) and dynamic (forest classes) covariates.

The forest classes are provided by the WB_HartJohnstoneForestClasses module which 
is also based on WB_VegBasedDrainage making an optional cyclic dependency between 
the two modules. Normally a first run of WB_HartJohnstoneForestClasses during 
module initialization will classify the forest to classes 1-6, without taking 
drainage into account. A WB_VegBasedDrainage map will then be computed and used, 
at the next simulation step, by the WB_HartJohnstoneForestClasses module to 
refine the classification of spruce from classes 6 to classes 6 and 7.

1. Initialization:
   - WB_HartJohnstoneForestClasses runs without drainage
   - Forest is classified into classes 1–6

2. Drainage estimation:
   - WB_VegBasedDrainage uses forest classes to compute drainage

3. Simulation steps:
   - WB_HartJohnstoneForestClasses re-runs
   - Spruce pixels are refined into well- or poorly-drained classes (6 or 7)

WB_VegBasedDrainage is dynamic because it depends on forest classes produced by 
WB_HartJohnstoneForestClasses which is itself dynamic. As relative biomass 
changes through time, forest classification is updated, which in turn influences 
the predicted drainage classes.

### Authors and Citation

* Pierre Racine <pierre.racine@sbf.ulaval.ca> [aut, cre]
* Andres Caseiro Guilhem <andres.caseiro-guilhem.1@ulaval.ca> [aut]
* Steven G. Cumming <stevec.boreal@gmail.com> [aut]

Racine, P., Caseiro Guilhem, A., Cumming, S.G. (2026) *WB_VegBasedDrainage: A SpaDES module for drainage classification in boreal forests.* SpaDES Module.


### Module Parameters

| Parameter | Class | Default | Description |
| --- | --- | --- | --- |
| WB_VegBasedDrainageTimeStep | integer | 1 | Simulation time step at which the drainage map is regenerated. |
| searchDistInPixelNb | integer | 1 | Distance, in number of pixel, to search for \"with value\" covariate values when plot points fall into NA pixels. This is generally when plot points are close to water which are set to NA by LandR. |


### Expected Module Inputs

| Input Object | Class | Description |
| --- | --- | --- |
| WB_HartJohnstoneForestClassesMap | SpatRast | Forest classification map produced by the WB_HartJohnstoneForestClasses modules or equivalent. |
| plotPoints | data.table | Plot data used to fit the drainage model. Should contain a latitude, a longitude and at least one column named \"drainage\" with drainage classes. Can also contain a "\"standtype\". Stand type will be extracted from the WB_HartJohnstoneForestClasses map if the corresponding module is present. Users can provide their own plot data in the form of a data.table as long as it follows the same schema. |
| MRDEMMap | SpatRast | Medium Resolution Digital Elevation Model (MRDEM) for Canada or an equivalent DEM. |
| TWIMap | SpatRast | TWI raster map automatically derived from the MRDEM raster or an equivalent TWI raster. |
| DownslopeDistMap | SpatRast | Downslope distance to water map automatically derived from the MRDEM raster or an equivalent downslope distance raster. |
| AspectMap | SpatRast | Aspect raster map automatically derived from the MRDEM raster or an equivalent aspect raster. |
| WB_VBD_ClayMap | SpatRast | Clay raster map automatically derived from the MRDEM raster or an equivalent clay raster. |
| WB_VBD_SandMap | SpatRast | Sand raster map automatically derived from the MRDEM raster or an equivalent sand raster. |
| WB_VBD_SiltMap | SpatRast | Silt raster map automatically derived from the MRDEM raster or an equivalent silt raster. |
| WB_VBD_BDMap | SpatRast | Bulk density raster map automatically derived from the MRDEM raster or an equivalent bulk density raster. |
| WB_VegBasedDrainageModel | TBD | Model used to fit the plot data against all the covariates. |

### Module Outputs

| Output Object | Class | Description |
| --- | --- | --- |
| WB_VegBasedDrainageMap | SpatRast | Raster map classified into 2 drainage classes (1 = poorly-drained, 2 = well-drained). |

### Code

The code is available here: https://github.com/pedrogit/WB_VegBasedDrainage

### Minimal Self Contained Workflow Example

Soon...



