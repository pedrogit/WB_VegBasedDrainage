library(caret)
library(dplyr)
library(terra)

# plot data
df <- read.csv(".../data.csv") # columns: plot ID, drainage, latitude, latitude (maybe include a vegetation-based indicator)

# covariates (rasters)
twi_raster <- rast("/TWI.tif")
clay <- rast(".../Clay_X0_5_cm_100m1980.tif") # source: Soil Landscape Grids of Canada, 100 m
sand <- rast(".../Sand_X0_5_cm_100m1980-2000v1.tif") # source: Soil Landscape Grids of Canada, 100 m
silt <- rast(".../Silt_X0_5_cm_100m1980-2000v1.tif") # source: Soil Landscape Grids of Canada, 100 m
bulk_density <- rast("/BD_X0_5_cm_100m1980-2000v1.tif") # Soil Canada, 100m 
aspect <- ...
...etc

# extracting covariates data for plot
points <- vect(df, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")  # WGS84

points <- project(points, crs(twi_raster))  # Reproject if necessary
twi_values <- extract(twi_raster, points)
df <- cbind(df, twi_values[, -1])  # Remove ID column from extract
colnames(df)[colnames(df) == "twi_values[, -1]"] <- "TWI_250m_average"

points <- project(points, crs(clay))  # Reproject if necessary
clay_values <- extract(clay, points)
df <- cbind(df, clay_values[, -1])  # Remove ID column from extract
colnames(df)[colnames(df) == "clay_values[, -1]"] <- "clay"

points <- project(points, crs(sand))  # Reproject if necessary
sand_values <- extract(sand, points)
df <- cbind(df, sand_values[, -1])  # Remove ID column from extract
colnames(df)[colnames(df) == "sand_values[, -1]"] <- "sand"

point <- project(points, crs(silt))
silt_values <- extract(silt, points)
df <- cbind(df, silt_values[, -1])  # Remove ID column from extract
colnames(df)[colnames(df) == "silt_values[, -1]"] <- "silt"

point <- project(points, crs(bulk_density))
bulk_density_values <- extract(bulk_density, points)
df <- cbind(df, bulk_density_values[, -1])  # Remove ID column from extract
colnames(df)[colnames(df) == "bulk_density_values[, -1]"] <- "bulk_density"

...etc


df <- na.omit(df)

# Splitting data into train and test set
inTraining <- createDataPartition(df$standtype, p = 0.7, list = FALSE)
trainSet <- df[inTraining, ]
testSet <- df[-inTraining, ]

# Fix class levels BEFORE modeling
levels(trainSet$standtype) <- make.names(levels(trainSet$standtype))
levels(testSet$standtype) <- make.names(levels(testSet$standtype))


# Defining conditions for training such as cross-validation and how to search
# for hyperparameters
fitControl <- trainControl(method = "repeatedcv",
                           number = 5,
                           repeats = 5,
                           search = "random")

# Fitting the model
modelFit <- train(
  standtype ~ ., # Can consider subsets of covariates here say; e.g. `slope + elevation + age`; `.` considers all the predictors with no interaction terms
  data = trainSet,
  method = "rf", # "rf" = random forest or "xgbTree" = boosted regression trees, # See topepo.github.io/caret for more model tags that can be used here
  trControl = fitControl,
  tuneLength = 10,
  verbose = FALSE
)

# Outputs a summary of the trained model
modelFit

# Predict on the test set # but we can also predict other datasets (e.g. Racine's Vegetation map)
predictions <- predict(modelFit, testSet)

# Look at the confusion matrix to evaluate model performance on the test set
confusionMatrix(predictions, testSet$standtype)

# Show variable importance
varImp(modelFit)
