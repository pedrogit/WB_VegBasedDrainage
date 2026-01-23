# WB_VegBaseDrainage

WB_VegBaseDrainage is a SpaDES module associated with the LandR set of 
modules for forest biomass and succession simulation. It is part of a series of 
module modelling boreal forest in Western Canada adding some layers of 
information to the LandR ecosystem. Those modules are:

- [WB_HartJohnstoneForestClasses](https://github.com/pedrogit/WB_HartJohnstoneForestClasses): classify LandR forested pixels to 6 (or 7) classes.
- [WB_VegBaseDrainage](https://github.com/pedrogit/WB_VegBasedDrainage): this module providing two drainage classes.
- [WB_NonForestedVegClasses](https://github.com/pedrogit/WB_NonForestedVegClasses): determine the land cover areas not covered by LandR.
- [WB_LichenBiomass](): generate a map of lichen biomass.

As their names suggest, those modules were developed using field data collected 
in the Western boreal forest of Canada. Their main aim is to model lichen biomass 
for caribou conservation.

At each simulation step, WB_VegBaseDrainage create a drainage raster map of 
well-drained and poorly-drained pixels for the forested areas covered by LandR 
only. It produces a raster with two classes:

- poorly drained (1)
- well drained (2)

WB_HartJohnstoneForestClasses (and hence LandR) pixels set to NA are also set to 
NA by WB_VegBaseDrainage.

WB_VegBaseDrainage use a machine learning algorithm to fit drainage, observed on 
the field, against a set of static covariates: Forest classes, TWI, downslope 
distance, aspect, ecoprovince and soil data (clay, sand", silt, bulk_density).

The forest classes are provided by the WB_HartJohnstoneForestClasses module which 
is also based on WB_VegBasedDrainage making an optional cyclic dependency between 
the two modules. Normally a first run of WB_HartJohnstoneForestClasses during 
module initialization will classify the forest to classes 1-6, without taking 
drainage into account. A WB_VegBasedDrainage maps will then be computed and used, 
at the next simulation step, by the WB_HartJohnstoneForestClasses module to refine 
the classification of spruce from classes 1-6 to classes 1-7.

WB_VegBaseDrainage is dynamic in that WB_HartJohnstoneForestClasses is dynamic: 
At each time step relative biomass changes, modifying how the forested pixels 
are classified and this have an impact on drainage.


### Authors

Pierre Racine <pierre.racine@sbf.ulaval.ca> [aut, cre]
Andres Caseiro Guilhem <andres.caseiro-guilhem.1@ulaval.ca> [aut]

### Module Parameters

| Parameter | Class | Default | Description |
| --- | --- | --- | --- |
| WB_VegBasedDrainageTimeStep | integer | 1 | Simulation time step at which the drainage map is regenerated. |
| searchDistInPixelNb | integer | 1 | Distance, in number of pixel, to search for \"with value\" covariate values when plot points fall into NA pixels. This is generally when plot points are close to water which are set to NA by LandR. |


### Expected Module Inputs

| Input Object | Class | Description |
| --- | --- | --- |
| WB_HartJohnstoneForestClassesMap | SpatRast | Forest classification map produced by the WB_HartJohnstoneForestClasses modules or equivalent. |
| plotPoints | data.table | Plot data used to fit the drainage model. Should contain a latitude, a longitude and at least one column named \"drainage\" with drainage classes. Can also contain a "\"standtype\". Standtype will be extracted from the WB_HartJohnstoneForestClasses map if the corresponding module is present. Users can provide their own plot data in the form of a data.table as long as it follows the same schema. |
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
| WB_VegBasedDrainageMap | SpatRaster | Raster map classified into 2 drainage classes. |

### Code

The code is available here: https://github.com/pedrogit/WB_VegBasedDrainage

### Minimal Self Contained Example



