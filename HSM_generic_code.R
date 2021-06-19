#*******************************************************
##  Habitat Suitability Modeling - example code  ##
## == Drew Foster == ##############################
#*******************************************************

#*******************************************************
#*#* Adapted from Guisan et al 2017 (pp. 357-380)*
#* Guisan, A., Thuiller, W., & Zimmermann, N.E. (2017).
#* Habitat Suitability and Distribution Models with Applications in R.
#* Cambridge University Press: Cambridge, UK.

#*******************************************************
#*#* Project: M.S. Thesis Example HSM code
#* Thesis title: A UAS-based approach toward habitat suitability of a rare endemic plant
#* School of Environmental and Forest Sciences, College of Environment, University of Washington
#* 2021
#* In Supplemental Materials

#*******************************************************
#*#* This is generic code, similar to what was run in Andrew Foster's M.S. Thesis
#* No species occurrence data is provided in order to protect the location
#* of the rare plant, Hackelia venusta, used as a focal species in the thesis

#*******************************************************
#*#* This code follows existing biomod2 tutorials and should be easily adaptable
#* for other purposes, and also includes code for importing GBIF records
#* See https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/biomod2/inst/doc/Simple_species_modelling.pdf?root=biomod
#* for another simple worked example


# Load required and other useful packages

library(ade4)
library(sp)
library(sf)
library(raster)
library(dismo)
library(rgbif)
library(ggmap)
library(ggplot2)
library(rgdal)
library(rgeos)
library(biomod2)
library(sdmpredictors)
library(maptools)
library(marmap)
library(rasterVis)
library(dplyr)
library(usdm)
library(RStoolbox)

# set working directory
setwd('D:/WORKING_DIRECTORY') # Change to your appropriate local working directory

####==== Step 1: Load Species Occurrence Data ====####

# # Read in species occurrence data
# Either from Working Drive, GBIF, or other
# Most biomod2 models work with Presence-Only OR Presence-Absence data
# For abundance data, choose different model approach, or convert to PO or PA

# Example for Presence-Absence data, with longitude; latitude; 1/0 for P/A
PA <- read.csv("PAdata.csv", header = TRUE)
PAnames <- c("lon", "lat", "presAbs")
names(PA) <- PAnames

# Check occurrence data, remove NAs if needed
head(PA)

# # Subset data if desired or appropriate for study

# Subset data into 70% train, 30% test
subset30 <- sample(nrow(PA), 30)
PA.train <- PA[-subset30,]
PA.test <- PA[subset30,]

head(PA.train)
head(PA.test)

# # Create response vector from P/A points
# respVec.train <- as.numeric(PA.train$presAbs)
# respVec.test <- as.numeric(PA.test$presAbs)
respVec <- as.numeric(PA$presAbs)

# # Visually Inspect Points

plot(PA$lon, PA$lat,
     col = "green")




#*#* Optionally, Obtain Species Occurrence from GBIF database
#* See https://www.gbif.org/ for more information
#* Note: Consider origin of occurrence records, and associated bias


# # Load species occurrence data from GBIF database
# # First, Identify species of interest: Replace Genus species with species of interest
key <- name_backbone(name = 'Genus species',
                      rank='species')$usageKey

# Run search and download records which contain coordinates, this takes a bit of time
focalSpecies <- occ_search(taxonKey=key,
                   limit=5000,
                   hasCoordinate = TRUE,
                   country=c("US", "ES","FR")) # Can select multiple country codes if you want
# # Note that Data is stored separately for each country

# # Extract Data into data frames, by individual country
# # select columns == "decimalLatitude", "decimalLongitude", "species"
species.data.US <- data.frame(focalSpecies$US$data[,c("decimalLatitude", "decimalLongitude", "species")])
species.data.ES <- data.frame(focalSpecies$ES$data[,c("decimalLatitude", "decimalLongitude", "species")])
species.data.FR <- data.frame(focalSpecies$FR$data[,c("decimalLatitude", "decimalLongitude", "species")])

# Or, combine all occurrences into one data frame
species.data <- rbind(species.data.US,
                   species.data.ES,
                   species.data.FR
                   )

# # Write csv for quicker access next time around
write.csv(species.data, file = "gbifOccurrence.csv")

# # Read csv if records already imported from GBIF
species.data <- read.csv("gbifOccurrence.csv", header = TRUE)

# # GBIF Spp Occurrence Data Q/C

# remove any N/As
species.data <- na.omit(species.data)

# Visualize species occurrences on a simple map
data("wrld_simpl")
plot(wrld_simpl,
     axes = TRUE, col = "light yellow")
box()
points(species.data$decimalLongitude, species.data$decimalLatitude,
       col = "red",
       pch = 16,
       cex = 0.75)

#*#* Inspect and/or remove outliers, address potential issues in occurrence data

####==== Step 2: Prepare Env'tal Predictor Rasters ====####

# Environmental Covariates
# Import single band raster data
# This could be directly imported from bioclim/worldclim depending on study
# Or import your own, can pre-process separately or in R

red <- raster("predictorRasters/REDBAND.tif")
blue <- raster("predictorRasters/BLUEBAND.tif")
green <- raster("predictorRasters/GREENBAND.tif")
nir <- raster("predictorRasters/NIRBAND.tif")
slope <- raster("predictorRasters/SLOPE.tif") # derived from DEM
aspect <- raster("predictorRasters/ASPECT.tif") # derived from DEM
vdvi <- raster("predictorRasters/VDVI.tif") # pre-calculated
ndvi <- raster("predictorRasters/NDVI.tif") # pre-calculated

# Set names for raster layers
names(red) <- 'Red'
names(blue) <- 'Blue'
names(green) <- 'Green'
names(nir) <- 'NIR'
names(slope) <- 'Slope'
names(aspect) <- 'Aspect'
names(vdvi) <- 'VDVI'
names(ndvi) <- 'NDVI'


## Inspect Rasters for QA/QC

# quick visual inspection may be useful
par(mfrow = c(2,2))
plot(red); plot(blue); plot(green); plot(nir)
plot(slope); plot(aspect); plot(ndvi); plot(vdvi)

# Inspect extent and resolution, make sure all are same
bbox(red); ncol(red); nrow(red); res(red)
bbox(blue); ncol(blue); nrow(blue); res(blue)
bbox(green); ncol(green); nrow(green); res(green)
bbox(nir); ncol(nir); nrow(nir); res(nir)
bbox(slope); ncol(slope); nrow(slope); res(slope)
bbox(aspect); ncol(aspect); nrow(aspect); res(aspect)
bbox(vdvi); ncol(vdvi); nrow(vdvi); res(vdvi)
bbox(ndvi); ncol(ndvi); nrow(ndvi); res(ndvi)

# check histograms for outliers
hist(red); hist(blue); hist(green); hist(nir);
hist(slope); hist(slope); hist(ndvi); hist(vdvi)

par(mfrow = c(1,1))

# Check projection
projection(blue); projection(green); projection(red);
projection(nir); projection(ndvi); projection(vdvi);
projection(slope); projection(aspect);
# reset so they all match if necessary

# # RESAMPLE LAYERS AS NEEDED TO MATCH & STACK

# try bilinear interpolation, better for continuous values
# resample all to coarsest resolution raster
blue <- resample(blue, slope, method = "bilinear")
green <- resample(green, slope, method = "bilinear")
red <- resample(red, slope, method = "bilinear")
nir <- resample(nir, slope, method = "bilinear")
ndvi <- resample(ndvi, slope, method = "bilinear")
vdvi <- resample(vdvi, slope, method = "bilinear")
aspect <- resample(aspect, slope, method = "bilinear")

# go back up to inspection to ensure all layers exactly the same

# # check to see if rasters Stack, 
# Do not edit Stack, do nothing to result in a Brick!
# biomod2 will not accept a raster Brick, only a Stack!

envtalStack <- stack(
    blue,
    green,
    red,
    nir,
    ndvi,
    vdvi,
    slope,
    aspect
)

#*#*#* Create new Stack of environmental variables to predict in space or time
#* Ensure same set of variables, 
#* Same resolution & projection as predictor variables
#* Can be different extent

newRed <- raster("projectionRasters/newREDBAND.tif")
newBlue <- raster("projectionRasters/newBLUEBAND.tif")
newGreen <- raster("projectionRasters/newGREENBAND.tif")
newNIR <- raster("projectionRasters/newNIRBAND.tif")
newSlope <- raster("projectionRasters/newSLOPE.tif")
newAspect <- raster("projectionRasters/newASPECT.tif")
newNDVI <- raster("projectionRasters/newNDVI.tif")

#* Above layers could be future climate projections from worldclim,
#* Or different region than predictors

# Set names for raster layers
names(newRed) <- 'Red'
names(newBlue) <- 'Blue'
names(newGreen) <- 'Green'
names(newNIR) <- 'NIR'
names(newSlope) <- 'Slope'
names(newAspect) <- 'Aspect'
names(newVDVI) <- 'VDVI'
names(newNDVI) <- 'NDVI'

# RUN THROUGH SAME QA/QC measures as above if necessary for new rasters

newEnvtalStack <- stack(
    newBlue,
    newGreen,
    newRed,
    newNIR,
    newNDVI,
    newVDVI,
    newSlope,
    newAspect
)



####==== Check for multicollinearity and correlation ====####

# Recognizing some autocorrelation may be acceptable/unavoidable for study

# # check for correlation
envtalStack.corr<-layerStats(envtalStack, 'pearson', na.rm=T)
envtalStack.corr # this checks for Pearson's correlation coefficient
envtalStack.cov <- layerStats(envtalStack, 'cov', na.rm = T)
envtalStack.cov

# Check Variance Inflation Factor
usdm::vifstep(envtalStack)
usdm::vifcor(envtalStack, th = 0.7)

# Remove rasters and restack in above step as necessary

## Alternatively, Run PCA to visualize variable relationships
library(ade4)

# Occurrence data needs to be loaded first
occurPoints <- data.frame(PA[,c("lon","lat")])
points.cell.id <- cellFromXY(subset(envtalStack,1), occurPoints)

envtalStack_df <- na.omit(as.data.frame(envtalStack))
head(envtalStack_df)
pca_Core <- dudi.pca(envtalStack_df, scannf = F, nf = 2)

# Visually check for outliers

plot(pca_Core$li[,1:2])

par(mfrow=c(1,2))
s.class(pca_Core$li[,1:2],
        fac = factor(rownames(envtalStack_df) %in% points.cell.id,
                     levels = c("FALSE", "TRUE"),
                     labels = c("background", "Occurrence")),
        col = c("red", "blue"),
        csta = 0,
        cellipse = 2,
        cpoint = .3,
        pch = 16)

s.corcircle(pca_Core$co, clabel = 1)

# Perfrom Visual check for points within raster extent
plot(ndvi)
points(PA$lon, PA$lat,
     col = "green")

plot(ndviCore)
points(PA$lon, PA$lat,
       col = "blue", 
       pch = 16, cex=2)

# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * 
# After testing for multicollinearity and correlation,
# Evaluate which env'tal variables to bring in to model
# REDO envtalStack and occurrence points as necessary before running model
# * # * # * # * # * # * # * # * # * # * # * # * # * # * # * 


####==== Write Rasters ====####

# # This step may be helpful to write formatted rasters for quick future loading
# # Or convert to ascii format to run Maxent program

# Write rasters to ascii format to run in Maxent
writeRaster(ndvi,
            filename = "maxent/ndvi",
            format = "ascii",
            overwrite = TRUE)
#... continue for all variables or loop through





####==== Step 3: Habitat Suitability Modeling with biomod2 ====####

# # First prepare species data to feed into biomod2

# For Presence-Only data, Generate pseudo-absence points: 
# This is an example using:
# random sampling, repeated 3 times, with 1000 background points
format.biomod <- BIOMOD_FormatingData(
    resp.var = rep(1,nrow(PO)), # could replace PO with species.data from GBIF example above
    expl.var = envtalStack,
    resp.xy = PO[,c("lon","lat")], # or use species.data[,c("decimalLongitude","decimalLatitude")]
    resp.name = "focalSpecies",
    PA.nb.rep = 3,
    PA.nb.absences = 1000,
    PA.strategy = "random"
)

# For Presence-Absence response data
# This example shows using full dataset,
# Optionally, add in independent testing data
format.biomod <- BIOMOD_FormatingData(
    resp.var = respVec,
    expl.var = envtalStack,
    resp.xy = PA[,c("lon","lat")],
    resp.name = "focalSpecies"
    # eval.resp.var = respVec.test, # Independent test data
    # eval.expl.var = envtalStack,
    # eval.resp.xy = PA.test[,c("lon","lat")]
)
# check format
format.biomod

# # Visualise pseudo-absences if generated
plot(format.biomod)

# Set modeling options for individual models
# See ?BIOMOD_ModelingOptions for detailed options and settings;
# this example shows options used in MS Thesis,
biomod.opt <- BIOMOD_ModelingOptions(
    GLM = list(type = "quadratic", interaction.level = 1),
    GBM = list(n.trees = 2500),
    RF = list(ntree = 1000),
    GAM = list(algo = "GAM_mgcv")
)

# Now Run all models! this can take a bit of time...
# See ?BIOMOD_Modeling for detailed arguments
# This example shows 70% subsetting of data for training, and 10 evaluation runs
single.models <- BIOMOD_Modeling(
    data = format.biomod,
    models = c("GLM", "GBM", "RF", "GAM"),
    models.options = biomod.opt,
    NbRunEval = 10,
    DataSplit = 70,
    VarImport = 5,
    models.eval.meth = c("KAPPA","TSS","ROC", "ACCURACY"),
    do.full.models = FALSE,
    modeling.id = "modelsRun01"
)

####==== Single model evaluations ====####

# Examine output, check warnings, view evaluations
single.model.scores <- get_evaluations(single.models)

dim(single.model.scores)
dimnames(single.model.scores)
single.model.scores
single.models.avg <- apply(single.model.scores, 2, mean)
single.models.avg
single.models.df <- as.data.frame(single.model.scores)
single.models.df

# Visual assessment of model performance by type of model
models_scores_graph(single.models, by = "models",
                    metrics = c("ROC", "TSS"), #can change metric on x/y axes
                    xlim = c(0.25,1), # cut off models <0.25
                    ylim = c(0.25,1),
                    main = "Model Performance")

# by evaluation run
models_scores_graph(single.models, by = "cv_run" , metrics = c("ROC","TSS"), 
                    xlim = c(0.25,1), ylim = c(0.25,1))
# by 
models_scores_graph(single.models, by = "data_set" , metrics = c("ROC","TSS"), 
                    xlim = c(0.25,1), ylim = c(0.25,1))

# # See variable importance by algorithm, mean of evaluation runs
# Note: not percent importance, won't add up to 100
(single.models.var.import) <- get_variables_importance(single.models)
apply(single.models.var.import, c(1,2), mean)



# # View variable response plots for probability of occurrence
library(rasterVis)

response_glm <- BIOMOD_LoadModels(single.models, models='GLM')
response_gbm <- BIOMOD_LoadModels(single.models, models='GBM')
response_rf <- BIOMOD_LoadModels(single.models, models='RF')
response_gam <- BIOMOD_LoadModels(single.models, models='GAM')

# Check response plots
glm_eval_strip <- biomod2::response.plot2(
    models  = response_glm,
    Data = get_formal_data(single.models,'expl.var'), 
    show.variables= get_formal_data(single.models,'expl.var.names'),
    do.bivariate = F, # switch to T for bivariate responses, bad idea if lots of variables
    fixed.var.metric = 'mean', # switch to 'median', see if it makes a difference
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(single.models,'resp.var'))

gbm_eval_strip <- biomod2::response.plot2(
    models  = response_gbm,
    Data = get_formal_data(single.models,'expl.var'), 
    show.variables= get_formal_data(single.models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var.metric = 'mean',
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(single.models,'resp.var'))

rf_eval_strip <- biomod2::response.plot2(
    models  = response_rf,
    Data = get_formal_data(single.models,'expl.var'), 
    show.variables= get_formal_data(single.models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var.metric = 'mean',
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(single.models,'resp.var'))

gam_eval_strip <- biomod2::response.plot2(
    models  = response_gam,
    Data = get_formal_data(single.models,'expl.var'), 
    show.variables= get_formal_data(single.models,'expl.var.names'),
    do.bivariate = FALSE,
    fixed.var.metric = 'mean',
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(single.models,'resp.var'))

####==== Alternative Response plot formatting ====####
library(dplyr)
rf.expl.var <- get_formal_data(single.models, 'expl.var')
rf.expl.var.names <- get_formal_data(single.models, 'expl.var.names')

rp.dat <- 
    response.plot2(
        models = response_glm,
        Data = get_formal_data(single.models,'expl.var'), 
        show.variables= get_formal_data(single.models,'expl.var.names'),
        do.bivariate = F, # toggle T: 3D, F: 2D
        fixed.var.metric = 'mean', # mean is default, can choose median, min, max
        legend = F,
        col = c("purple", "green"), # choose colors for bivariate plots
        plot = FALSE,
        display_title = F,
        data_species = get_formal_data(single.models,'resp.var')
    )
rp.dat %>%
    ## transform the pred.name to extract model, cv run and data info
    mutate(
        species = pred.name %>% strsplit('_') %>% sapply(function(x) x[1]),
        pa.dat = pred.name %>% strsplit('_') %>% sapply(function(x) x[2]),
        cv.rep = pred.name %>% strsplit('_') %>% sapply(function(x) x[3]),
        model = pred.name %>% strsplit('_') %>% sapply(function(x) x[4])
    ) %>%
    ggplot(
        aes(
            x = expl.val,
            y = pred.val,
            colour = model,
            group = pred.name
        )
    ) +
    geom_line(size = 1, colour = "black") +
    facet_wrap(~ expl.name, scales = 'free_x') + 
    ylim(0,1) +
    labs(
        x = '',
        y = 'probability of occurence',
        colour = 'model type'
    ) +
    theme_minimal() +
    theme(
        legend.position = 'right'
    )


# get detailed outputs from individual models
bm.form.rf <- lapply(response_rf, function(x) get_formal_model(get(x)))
names(bm.form.rf) <- response_rf
bm.form.rf

bm.form.glm <- lapply(response_glm, function(x) get_formal_model(get(x)))
names(bm.form.glm) <- response_glm
bm.form.glm

#... continue for all single models, or loop

####==== Step 4: Ensemble single models: committee averaging and mean weight ====####

## Determine thresholding for which single models to ensemble
scores_all <- get_evaluations(single.models)
scores_TSS <- as.numeric(scores_all["TSS","Testing.data",,,])
score_thresh <- mean(tail(sort(scores_TSS),10)) # take the mean of the top ten
score_thresh

# This step runs ensemble model(s) based on a chosen threshold metric
ensemble.models <- BIOMOD_EnsembleModeling(
    modeling.output = single.models,
    em.by = 'PA_dataset+repet',
    eval.metric = 'TSS',
    eval.metric.quality.threshold = score_thresh,
    models.eval.meth = c('KAPPA', 'TSS', 'ROC', 'ACCURACY'),
    prob.mean = FALSE,
    prob.cv = TRUE,
    committee.averaging = TRUE,
    prob.mean.weight = TRUE,
    VarImport = 5
)


# Get evaluations for ensemble model(s)
(ensemble.scores <- get_evaluations(ensemble.models))

# Get variable importance for ensemble models
(ensemble.var.import <- get_variables_importance(ensemble.models))




#*#*#* this next step writes big files, bigger for higher resolution, 
#* change directory if needed

####==== Project Ensemble model to new areas or times ====####

# First project single models to new space or time
model.projection <- BIOMOD_Projection(
    modeling.output = single.models,
    new.env = newEnvtalStack, # Provide same set of env'tal variables in new space or time 
    proj.name = "Projection",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
)


# Generate binary presence-absence maps using a threshold method
# This produces binary rasters from ensemble predictions, could save to drive
model.forecast <- BIOMOD_EnsembleForecasting(
    EM.output = ensemble.models,
    projection.output = model.projection,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
)



# Plot Ensemble Models showing Committee Averaging & Weighted Means
plot(model.forecast, 
     str.grep = "EMca|EMwmean")



# #  Optionally, write and save projection rasters to drive

# Binary Ensemble model predictions
ca.binary <- raster("focalSpecies/proj_Current_Prediction/individual_projections/focalSpecies_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
wm.binary <- raster("focalSpecies/proj_Current_Prediction/individual_projections/focalSpecies_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")

# Or, use continuous, non-binary predictions, convert 0-1 probability of occurrence
ca.continuous <- raster("focalSpecies/proj_Current_Prediction/individual_projections/focalSpecies_EMcaByTSS_mergedAlgo_mergedRun_mergedData.img")
wm.continuous <- raster("focalSpecies/proj_Current_Prediction/individual_projections/focalSpecies_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img")
ca.continuous <- ca.continuous/1000
wm.continuous <- wm.continuous/1000

writeRaster(ca.binary, "binaryCAprediction", format = "GTiff")
writeRaster(wm.binary, "binaryWMprediction", format = "GTiff")
writeRaster(ca.continuous, "contCAprediction", format = "GTiff")
writeRaster(wm.continuous, "contWMprediction", format = "GTiff")


