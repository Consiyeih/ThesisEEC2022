###download data from gbif
rm(list=ls())
dev.off()
require(rgbif)
# fill in your gbif.org credentials 
user <- # your gbif.org username 
pwd <- # your gbif.org password
email <-  # your email 
myspecies<- "Eretmochelys imbricata" ##each species name

spp_species <- name_suggest(q=myspecies, rank = "species", limit = 500)
spp_species
species <- spp_species$data
key<- subset(species, species$canonicalName==myspecies)$key

occ_download(
  pred_in("taxonKey", key),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_gte("year", 1945),
  pred_lte("coordinateUncertaintyInMeters", 10000),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

#need to change download key based on above output
d <- occ_download_get('0381401-210914110416597', path="../zipdata/", overwrite = TRUE) %>%
  occ_download_import()

writepath <- paste("../rawdata/",myspecies,".csv",sep="")
write.csv(d,writepath,row.names = FALSE)

######
### Before running code need to setwd to code folder and assign myspecies<-
require(rgbif)
require(robis)
require(dplyr)
require(sf)
require(sp)
require(CoordinateCleaner)

start_time <- Sys.time()
# Creating empty data frame to track the number of records through cleaning and SDM
records_data <- data.frame(process=character(), records=numeric())

##Load and Clean Data from GBIF
spp_species <- name_suggest(q=myspecies, rank = "species", limit = 500)
spp_species #checking species records 

# download GBIF occurrence data for this species; this takes time if there are many data points!
gbif_data <- occ_data(scientificName = myspecies, 
                      hasCoordinate = TRUE, 
                      hasGeospatialIssue = FALSE,limit = 10, 
                      year = '1945,2023')

# take a look at the downloaded data:
gbif_data 
records <- gbif_data$meta$count

# download data from OBIS 
obis_data <- occurrence(myspecies)

# Add records to data frame 
records_data<- rbind(records_data, data.frame(process=c("gbif", "obis"), records=c(records, nrow(obis_data))))

# data downloaded to gbif account 
# fill in your gbif.org credentials 
#user <- "" # your gbif.org username 
#pwd <- "" # your gbif.org password
#email <- "" # your email 
#myspecies<- "Chelonia mydas"

# occ_download(
# pred_in("scientificName", myspecies),
# pred("hasCoordinate", TRUE),
# pred("hasGeospatialIssue", FALSE),
# pred_gte("year", 1945),
# pred_lt("coordinateUncertaintyInMeters", 10000),
# format = "SIMPLE_CSV",
# user=user,pwd=pwd,email=email
# )

filepath <- paste("../rawdata/", myspecies,".csv", sep="")
gbif_data<- read.csv(filepath)

# now let's filter for only human observations, observations and preserved specimens 
gbif_data2<- gbif_data%>%
  dplyr::select("taxonKey", "countryCode", "decimalLatitude","decimalLongitude",
                "basisOfRecord", "occurrenceStatus","individualCount",
                "coordinateUncertaintyInMeters")%>%
  filter(!basisOfRecord%in%c("FOSSIL_SPECIMEN", "LIVING_SPECIMEN"))
gbif_data2$species<- myspecies

names(gbif_data2)
sort(unique(gbif_data2$individualCount))  # notice if some points correspond to zero abundance
sort(unique(gbif_data2$occurrenceStatus))  # check for different indications of "absent", which could be in different languages! and remember that R is case-sensitive

# also eliminate presences with reported coordinate uncertainty (location error, 
#spatial resolution) larger than 5 km (5000 m):
gbif_data3 <- subset(gbif_data2,gbif_data2$coordinateUncertaintyInMeters<10000)
gbif_data3 <-subset(gbif_data3,!is.na(gbif_data3$coordinateUncertaintyInMeters))
nrow(gbif_data3)

# let's do the same for the OBIS data and then combine the datasets 
obis_data2<- obis_data%>%
  filter(!basisOfRecord%in%c("FossilSpecimen", "LivingSpecimen"))%>%
  filter(date_year>=1945)%>%
  filter(coordinateUncertaintyInMeters<10000)%>%
  filter(!is.na(coordinateUncertaintyInMeters))%>%
  dplyr::select("countryCode", "decimalLatitude","decimalLongitude",
                "basisOfRecord", "occurrenceStatus","individualCount",
                "coordinateUncertaintyInMeters")
obis_data2$species<- myspecies
obis_data2$taxonKey<- NA

# Combining 
gbif_data3<- rbind(gbif_data3, obis_data2)

# let's remove duplicates 
dups <- duplicated(gbif_data3[, c('decimalLongitude', 'decimalLatitude')])
sum(dups)
gbif_data4 <- gbif_data3[!dups, ]
nrow(gbif_data4)

rm(list=ls()[! ls() %in% c("gbif_data4","records_data", "myspecies")])

# let's check the data we have is in the IUCN range of the species [NOTE: IUCN shape files need to be downloaded from IUCN separately]
rangepath <- paste("../iucn/",myspecies, "/data_0.shp", sep="")
range <- sf::st_read(rangepath)
spdf <- as(range, "Spatial")
proj4string(spdf)<- CRS("+init=epsg:4326")
spdf@data$species <- myspecies
pc <- spTransform(spdf, CRS("+init=epsg:4326")) 
names(spdf@data)[3] <- "binomial"
gbif_data5<- cc_iucn(
  x=gbif_data4,
  range = spdf,
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "species",
  buffer = 5, 
  value = "clean")

records_data<- rbind(records_data, data.frame(process="combine_cleaned", records=nrow(gbif_data5)))

#Saving Cleaned Data
filename<- paste("../data/",myspecies, "cleaned.csv", sep="")
filename_records <- paste("../results/", myspecies, "records.csv", sep="")
write.csv(gbif_data5,filename,row.names = F)
write.csv(records_data, filename_records,row.names = F)

##### environmental bariable selection


####species distribution code

getwd()
rm(list=ls())
setwd("/cloud/project/sdms/code/")
rm(list=ls()) 
species <- c("Anous minutus", "Anous stolidus", "Chelonia mydas",
             "Dermochelys coriacea", "Eretmochelys imbricata",
             "Galeocerdo cuvier", "Isurus oxyrinchus", "Katsuwonus pelamis",
             "Megaptera novaeangliae", "Mesoplodon europaeus", "Onychoprion fuscatus",
             "Pelagodroma marina", "Phaethon aethereus", "Physeter macrocephalus",
             "Prionace glauca", "Stenella frontalis", "Sula dactylatra", 
             "Sula leucogaster", "Sula sula", "Thunnus albacares", 
             "Thunnus thynnus", "Tursiops truncatus", "Xiphias gladius", "Oceanodroma leucorhoa")
myspecies<-species[2] #change the number 1 and 2 already done 3, 4, 5, 6, 7, 8, 9,
#10, 11, 12, 13, 14, 15, 16, 17, 18, 19 (error), 20 (big error), 21 (big error)

require(raster)
require(rgdal)
require(biomod2)
require(sp)
require(stringr)
require(dplyr)
require(tidyr)
#Surface.Chlorophyll.Min, Surface.Current.Velocity.Min, Present.Surface.Diffuse.attenuation.Min.BOv2_2,
#Surface.Primary.productivity.Max, Surface.Primary.productivity.Min, Surface.Salinity.Max,
#Surface.Temperature.Mean, Surface.Temperature.Range
# Reading Rasters
chloro_min <- raster("../environment/Present.Surface.Chlorophyll.Min.tif")
current_min <- raster("../environment/Present.Surface.Current.Velocity.Min.tif.BOv2_1.tif")
da_min <- raster("../environment/Present.Surface.Diffuse.attenuation.Min.BOv2_2.tif")
pp_max <- raster("../environment/Present.Surface.Primary.productivity.Max.tif")
pp_min <- raster("../environment/Present.Surface.Primary.productivity.Min.tif")
sal_max <- raster("../environment/Present.Surface.Salinity.Max.tif")
temp_mean <- raster("../environment/Present.Surface.Temperature.Mean.tif")
temp_range <- raster("../environment/Present.Surface.Temperature.Range.tif")
enviro<- stack(chloro_min, current_min, da_min,
               pp_max, pp_min, sal_max,
               temp_mean, temp_range)

proj4string(enviro)<- CRS("+init=epsg:4326")

rm(chloro_min, current_min, da_min, 
   pp_max, pp_min, sal_max, temp_mean, temp_range)

# ## Read in gbif species data 
filepath<- paste("../data/",myspecies,"cleaned.csv", sep="")
data<- read.csv(filepath)

## The biomod2 package (the one we are using to do the SDM,
## Thuiller, Georges, Engler, Breiner, 2016) need the
## data to be formatted so we can run the SDM analysis. In this case, we are
## using the spatial points object as response variable, the raster stack as
## explanatory variable and we are creating a folder that is gonna contain 
## our outputs (resp.name). Since we don't have absence points, we are creating
## pseudo-absence points. The number of pseudoabsences is based on this paper:
## (Barbet-Massin, Jiguet, Albert, & Thuiller, 2012)
## Pseudoabsences are meant to avoid the overfitting of the model. 
pseudo<- nrow(data)

myspecies_biomod <- BIOMOD_FormatingData(resp.var = rep(1,nrow(data)),
                                         expl.var = enviro,
                                         resp.xy = data[, c("decimalLongitude", "decimalLatitude")],
                                         resp.name=myspecies,PA.nb.rep=2,
                                         PA.nb.absences=pseudo,
                                         PA.strategy="random")

myspecies_opt <- BIOMOD_ModelingOptions() #default settings

myspecies_models <- BIOMOD_Modeling(
  data = myspecies_biomod, models = c("GLM", "GBM", "RF", "GAM"),
  models.options = myspecies_opt, NbRunEval = 3, DataSplit = 70,
  models.eval.meth = c("ROC", "TSS", "KAPPA"),
  VarImport = 3, do.full.models = F, modeling.id = "ex2", 
  SaveObj = TRUE)

# get model evaluation scores
myspecies_models_scores <- as.data.frame(get_evaluations(myspecies_models))
myspecies_models_scores

path<-paste("../results/",myspecies,"model_scores.csv", sep="")
write.csv(myspecies_models_scores, path)

test_data <- myspecies_models_scores%>%select(starts_with("Testing.data"))
mean_tss<- rowMeans(test_data[2,], na.rm = T)

ensemble_models<- test_data[2,]%>%pivot_longer(names_to = "model", values_to = "tss", 
                                               cols = colnames(test_data))%>%
  filter(tss>mean_tss)

models_for_ensemble_path<-paste("../results/",myspecies,"models_for_ensemble_tss.csv", sep="")
write.csv(ensemble_models, models_for_ensemble_path)

models_scores_graph(myspecies_models, by = "models" ,
                    metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(myspecies_models, by = "cv_run" ,
                    metrics = c("ROC","TSS"), xlim = c(0.5,1), ylim = c(0.5,1))
models_scores_graph(myspecies_models, by = "data_set" ,
                    metrics = c("ROC", "TSS"), xlim = c(0.5, 1), ylim = c(0.5, 1))

# get variable importance 
myspecies_models_var_import <-
  as.data.frame(get_variables_importance(myspecies_models))

mean_var_import<-as.data.frame(apply(myspecies_models_var_import, c(1,2), mean))
var_path<-paste("../results/",myspecies,"model_var_import.csv", sep="")
write.csv(myspecies_models_var_import, var_path)

mean_var_path<-paste("../results/",myspecies,"model_mean_var_import.csv", sep="")
write.csv(mean_var_import, mean_var_path)

# build ensemble model 
myspecies_ensemble_models <- BIOMOD_EnsembleModeling(
  modeling.output = myspecies_models, em.by = "all",
  eval.metric = "TSS", eval.metric.quality.threshold = mean_tss,
  models.eval.meth = c("KAPPA", "TSS", "ROC"), prob.mean = FALSE,
  prob.cv = TRUE, committee.averaging = TRUE,
  prob.mean.weight = TRUE, VarImport = 3)

# get model evaluation scores
myspecies_ensemble_models_scores <- as.data.frame(get_evaluations(myspecies_ensemble_models))
eval_path<-paste("../results/",myspecies,"ensemble_model_scores.csv", sep="")
write.csv(myspecies_ensemble_models_scores, eval_path)

# get variable importance 
myspecies_ensemble_models_var_import <-
  as.data.frame(get_variables_importance(myspecies_ensemble_models))
var_path<-paste("../results/",myspecies,"ensemble_model_var_import.csv", sep="")
write.csv(myspecies_ensemble_models_var_import, var_path)

# model projections for current 
myspecies_models_proj_current <- BIOMOD_Projection(
  modeling.output = myspecies_models, new.env = stack(enviro),
  proj.name = "current", binary.meth = "TSS",
  output.format = ".img", do.stack = FALSE, 
  selected.models = get_needed_models(myspecies_ensemble_models))

myspecies_ensemble_models_proj_current <-
  BIOMOD_EnsembleForecasting(EM.output = myspecies_ensemble_models,
                             projection.output = myspecies_models_proj_current,
                             binary.meth = "TSS", output.format = ".img",
                             do.stack = FALSE)

# produces weighted mean (suggested to be the best: https://doi.org/10.1111/j.1472-4642.2008.00491.x) 
# and committee average projections based on the ensemble models
###plot each species for relative probabilities

bin_path <- paste("../species_proj/", myspecies, "/EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img", sep="")
rel_path <- paste("../species_proj/", myspecies, "/EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img", sep="")
count <- paste("../plastics/ing_count.tif", sep="")
weight <- paste("../plastics/ing_weight.tif",sep="")
data_path <- paste("../species_data/ing.csv",sep="")

#Binary Projections
projections <- raster(bin_path)

#Relative Projections
relative_projections <- raster(rel_path)
relative_projections<- relative_projections/1000 #convert to probabilities 
name<-paste("../species_proj/",myspecies, "/", myspecies, "binary.png", sep="")


png(
  name,
  width     = 11.69,
  #11.69
  height    = 8.27,
  #8.27
  units     = "in",
  res       = 1200,
  pointsize = 2
)

par(oma=c(5,0,5,20))
par(mar=c(0,0,0,15))
plot(relative_projections, legend=FALSE, axes = FALSE, box=FALSE)
#plot(relative_projections,
#    legend= TRUE, cex=0.5, legend.width= 8,legend.shrink= 0.58,
#   legend.only=TRUE, axis.args=list(cex.axis=10),
#  add=TRUE, axes = FALSE, box=FALSE)

dev.off()


######## 4 dose-response ecotoxicology
### adding the plastic avarage and
##polygons and annual ingestion/entanglement rates t the data

###save the following code in a separate r script as dr_data_prep.R
require(raster)
require(sp)
require(sf)
require(exactextractr)
require(rgeos)
require(tidyr)
require(dplyr)
require(MASS)
require(rgdal)

bin_path <- paste("../species_proj/", myspecies, "/EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img", sep="")
rel_path <- paste("../species_proj/", myspecies, "/EMwmeanByTSS_mergedAlgo_mergedRun_mergedData.img", sep="")
count <- paste("../plastics/", ing_ent, "_count.tif", sep="")
weight <- paste("../plastics/", ing_ent, "_weight.tif",sep="")
data_path <- paste("../species_data/", ing_ent, ".csv",sep="")

#Binary Projections
projections <- raster(bin_path)

#Relative Projections
relative_projections <- raster(rel_path)
relative_projections<- relative_projections/1000 #convert to probabilities 

# Plastics data processing 
plastics_count <- raster(count)
plastics_weight <- raster(weight)

#Myspecies Coordinates and migration
data <- read.csv(data_path)
data_mig <- read.csv("../species_data/migration.csv")
data<- subset(data, data$species==myspecies)
data_mig<- subset(data_mig, data_mig$species==myspecies)

#Need to convert to same CRS and then do topological subsetting 
spatial_data <- SpatialPointsDataFrame(c(data[,c('decimalLongitude','decimalLatitude')]), data = data)
proj4string(spatial_data)<- c("+init=epsg:4326")
target_crs <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # This is your string assigned to an object
spatial_data<- spTransform(spatial_data, target_crs) #converting to correct CRS 

#Making Polygons 
projections[projections<1]<- NA # making absences 0 so the polygons are only for presences
projection_polygons <- sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(projections),
                                                   merge = TRUE,na.rm=TRUE)) # requires the sf, sp, raster and stars packages
conversion<- st_as_sf(projection_polygons)
sf::sf_use_s2(FALSE)
buffer<- st_buffer(conversion, dist=data_mig$decimalDist) #distance now in arc decimals 
projection_polygons<- st_make_valid(buffer)
validity<-st_is_valid(projection_polygons) # checks whether polygons are valid 
nrow(validity)-sum(validity) #should be zero 

#Subsetting for the relevant polygons
spatial_data2<-sf::st_as_sf(spatial_data)
projection_polygons2<-sf::st_as_sf(projection_polygons)

##Alternative of assessing the nearest polygon to point 
nearest <- sf::st_nearest_feature(spatial_data2, projection_polygons2, pairwise=FALSE, check_crs = TRUE) 
nearest # finds nearest polygon to every recording
selected_polygons<- projection_polygons2[nearest,] # subset for nearest polygons

#Extracting mean plastic count in each polygon per recording
#reprojecting relative probabilities 
relative_projections2 <- resample(relative_projections,plastics_count,method="bilinear")
plastics_extraction<- exact_extract(plastics_count, selected_polygons)
relatives_extraction<- exact_extract(relative_projections2, selected_polygons)
plastic_values <- lapply(plastics_extraction, function(x) { x[, "value"] })
relative_values <- lapply(relatives_extraction, function(x) { x[, "value"] })
list_mult<-Map('*',relative_values,plastic_values)
data$weighted_mean_count<- unlist(lapply(list_mult, mean, na.rm=TRUE))
data$weighted_median_count<- unlist(lapply(list_mult, median, na.rm=TRUE))
data$weighted_max_count<- unlist(lapply(list_mult, max, na.rm=TRUE))

#Extracting mean plastic weight in each polygon per recording
#reprojecting relative probabilities 
relative_projections2 <- resample(relative_projections,plastics_weight,method="bilinear")
plastics_extraction<- exact_extract(plastics_weight, selected_polygons)
relatives_extraction<- exact_extract(relative_projections2, selected_polygons)
plastic_values <- lapply(plastics_extraction, function(x) { x[, "value"] })
relative_values <- lapply(relatives_extraction, function(x) { x[, "value"] })
list_mult<-Map('*',relative_values,plastic_values)
data$weighted_mean_weight<- unlist(lapply(list_mult, mean, na.rm=TRUE))
data$weighted_median_weight<- unlist(lapply(list_mult, median, na.rm=TRUE))
data$weighted_max_weight<- unlist(lapply(list_mult, max, na.rm=TRUE))

## Adding Columns
data$percent<- data$count/data$total # overall rate 
data$percent_year <- data$percent/data$years # rate per year 
data<-subset(data, data$years>=1)

# Exporting the data 
export_path <- paste("../results/", myspecies, "_", ing_ent, ".csv",sep="")
write.csv(data, export_path, row.names = F)

###### save until here after 

###then run the following make sure you set your wd properly
#entanglement 
ing_ent <- "ent"

species<-c("Chelonia mydas","Dermochelys coriacea","Eretmochelys imbricata",
           "Megaptera novaeangliae","Physeter macrocephalus","Tursiops truncatus")

for(i in species){
  myspecies<- i
  source("dr_data_prep.R")
  print(paste(i, "is done!"))
}

#ingestion 
ing_ent <- "ing"

#ingestion
species<-c("Anous minutus","Anous stolidus","Chelonia mydas","Dermochelys coriacea","Eretmochelys imbricata",
           "Galeocerdo cuvier", "Oceanodroma leucorhoa", "Katsuwonus pelamis", 
           "Mesoplodon europaeus", "Onychoprion fuscatus", "Pelagodroma marina",
           "Phaethon aethereus","Physeter macrocephalus","Prionace glauca",
           "Stenella frontalis","Sula dactylatra","Sula leucogaster","Sula sula",
           "Thunnus albacares","Tursiops truncatus","Xiphias gladius" )

for(i in species){
  myspecies<- i
  source("dr_data_prep.R")
  print(paste(i, "is done!"))
  
}

############


##### now analyze the different lmmm and extract all of them to copy pste in a 
##excel document to have all the values from the models
rm(list=ls())
dev.off()
setwd("../../code/")
files<- list.files("../results/species_df/")
filepath<- paste("../results/species_df/", files[15], sep="") ###change from 1 to 27
##since there is a total of 27 species-interaction files
data<- read.csv(filepath)

# Linear Models with count 
lm_model_count <- lm((percent_year)~0+weighted_mean_count, data=data)
#sqrt
residuals_lm<- simulateResiduals(lm_model_count)
plot(residuals_lm)

summary(lm_model_count)
confint(lm_model_count)

plot(lm_model_count)


# Linear Models with count sqrt
lm_model_count_sqrt<- lm(sqrt(percent_year)~0+weighted_mean_count, data=data)
#sqrt
residuals_lm_sqrt<- simulateResiduals(lm_model_count_sqrt)
plot(residuals_lm_sqrt)

summary(lm_model_count_sqrt)
confint(lm_model_count_sqrt)

plot(lm_model_count_sqrt)

# Linear Mixed Model Method with count 
require(glmmTMB)
require(DHARMa)
lmm_model_count<- glmmTMB((percent_year)~0+weighted_mean_count+(0+weighted_mean_count|ocean)+(0+weighted_mean_count|decade.middle),data=data, family = "gaussian")

residuals_lmm<- simulateResiduals(lmm_model_count)
plot(residuals_lmm)

summary(lmm_model_count)
confint(lmm_model_count)

# Linear Mixed Model Method with count sqrt
require(glmmTMB)
require(DHARMa)
lmm_model_count_sqrt <- glmmTMB(sqrt(percent_year)~0+weighted_mean_count+(0+weighted_mean_count|ocean)+(0+weighted_mean_count|decade.middle),data=data, family = "gaussian")
#sqrt
residuals_lmm_sqrt<- simulateResiduals(lmm_model_count_sqrt)
plot(residuals_lmm_sqrt)

summary(lmm_model_count_sqrt)
confint(lmm_model_count_sqrt)

# Linear Models with weight
lm_model_weight <- lm((percent_year)~0+weighted_mean_weight, data=data)
#sqrt
residuals_lm_w<- simulateResiduals(lm_model_weight)
plot(residuals_lm_w)

summary(lm_model_weight)

confint(lm_model_weight)

plot(lm_model_count)

# Linear Models with weight sqrt
lm_model_weight_sqrt<- lm(sqrt(percent_year)~0+weighted_mean_weight, data=data)
#sqrt
residuals_lm_weight_sqrt<- simulateResiduals(lm_model_weight_sqrt)
plot(residuals_lm_weight_sqrt)

summary(lm_model_weight_sqrt)

confint(lm_model_weight_sqrt)

plot(lm_model_count_sqrt)

# Linear Mixed Model Method with weight
require(glmmTMB)
require(DHARMa)
lmm_model_weight<- glmmTMB((percent_year)~0+weighted_mean_weight+(0+weighted_mean_weight|ocean)+(0+weighted_mean_weight|decade.middle),data=data, family = "gaussian")
#sqrt
residuals_lmm_w<- simulateResiduals(lmm_model_weight)
plot(residuals_lmm_w)

summary(lmm_model_weight)
confint(lmm_model_weight)

# Linear Mixed Model Method with weight sqrt
require(glmmTMB)
require(DHARMa)
lmm_model_weight_sqrt<- glmmTMB(sqrt(percent_year)~0+weighted_mean_weight+(0+weighted_mean_weight|ocean)+(0+weighted_mean_weight|decade.middle),data=data, family = "gaussian")
#sqrt
residuals_lmm_w_sqrt<- simulateResiduals(lmm_model_weight_sqrt)
plot(residuals_lmm_w_sqrt)

summary(lmm_model_weight_sqrt)
confint(lmm_model_weight_sqrt)

## my try to extract values
lmcountestimate<-summary(lm_model_count)$coefficients[1, 1]
lmintervalmin<-confint(lm_model_count)[1,1]
lmintervalmax<-confint(lm_model_count)[1,2]
lmcountestimate_sqrt<-summary(lm_model_count_sqrt)$coefficients[1, 1]
lmintervalmin_sqrt<-confint(lm_model_count_sqrt)[1,1]
lmintervalmax_sqrt<-confint(lm_model_count_sqrt)[1,2]
#LMM how to extract residual???
lmmcountestimate<-summary(lmm_model_count)$coefficients$cond[1, 1]
lmmcountintervalmin<-confint(lmm_model_count)[1,1]
lmmcountintervalmax<-confint(lmm_model_count)[1,2]
lmmcountvarocean<-summary(lmm_model_count)$varcor$cond$ocean[1,1]
lmmcountvar_decade<-summary(lmm_model_count)$varcor$cond$decade[1,1]
#lmm count sqrt
lmmcountestimate_sqrt<-summary(lmm_model_count_sqrt)$coefficients$cond[1, 1]
lmmcountintervalmin_sqrt<-confint(lmm_model_count_sqrt)[1,1]
lmmcountintervalmax_sqrt<-confint(lmm_model_count_sqrt)[1,2]
lmmcountvarocean_sqrt<-summary(lmm_model_count_sqrt)$varcor$cond$ocean[1,1]
lmmcountvar_decade_sqrt<-summary(lmm_model_count_sqrt)$varcor$cond$decade[1,1]

#lm weight 
lm_weightestimate<-summary(lm_model_weight)$coefficients[1, 1]
lm_weightintervalmin<-confint(lm_model_weight)[1,1]
lm_weightintervalmax<-confint(lm_model_weight)[1,2]
lm_weightestimate_sqrt<-summary(lm_model_weight_sqrt)$coefficients[1, 1]
lm_weightintervalmin_sqrt<-confint(lm_model_weight_sqrt)[1,1]
lm_weightintervalmax_sqrt<-confint(lm_model_weight_sqrt)[1,2]

#LMM wight how to extract residual???
lmm_weightestimate<-summary(lmm_model_weight)$coefficients$cond[1, 1]
lmm_weightintervalmin<-confint(lmm_model_weight)[1,1]
lmm_weightintervalmax<-confint(lmm_model_weight)[1,2]
lmm_weightvarocean<-summary(lmm_model_weight)$varcor$cond$ocean[1,1]
lmm_weightvar_decade<-summary(lmm_model_weight)$varcor$cond$decade[1,1]
#lmm wight sqrt
lmm_weightestimate_sqrt<-summary(lmm_model_weight_sqrt)$coefficients$cond[1, 1]
lmm_weightintervalmin_sqrt<-confint(lmm_model_weight_sqrt)[1,1]
lmm_weightintervalmax_sqrt<-confint(lmm_model_weight_sqrt)[1,2]
lmm_weightvarocean_sqrt<-summary(lmm_model_weight_sqrt)$varcor$cond$ocean[1,1]
lmm_weightvar_decade_sqrt<-summary(lmm_model_weight_sqrt)$varcor$cond$decade[1,1]


#lmmcountvar_residuals<-summary(lmm_model_count)$varcor$Residual

hey<-data.frame(lmcountestimate,lmintervalmin, lmintervalmax, "", lmcountestimate_sqrt,
                lmintervalmin_sqrt, lmintervalmax_sqrt, "","",lmmcountestimate,lmmcountintervalmin,
                lmmcountintervalmax,
                lmmcountvarocean, lmmcountvar_decade, "", "", lmmcountestimate_sqrt,
                lmmcountintervalmin_sqrt, lmmcountintervalmax_sqrt, lmmcountvarocean_sqrt,
                lmmcountvar_decade_sqrt, "", "",
                lm_weightestimate,lm_weightintervalmin, lm_weightintervalmax, "",
                lm_weightestimate_sqrt,
                lm_weightintervalmin_sqrt, lm_weightintervalmax_sqrt, "",
                lmm_weightestimate,lmm_weightintervalmin,
                lmm_weightintervalmax,
                lmm_weightvarocean, lmm_weightvar_decade, "", "", lmm_weightestimate_sqrt,
                lmm_weightintervalmin_sqrt, lmm_weightintervalmax_sqrt, lmm_weightvarocean_sqrt,
                lmm_weightvar_decade_sqrt, "", "")


### plot for the model example Dermochelys coriacea ingestion
##first change to lmer
##predicted values
library(lme4)
me<- ggpredict(lm_model_count_sqrt, terms = c("weighted_mean_count"), back.transform = FALSE)


my_y_title <- expression(paste("a) ", "Entanglement (", italic("Dermochelys coriacea"), ")"))
x_axis<-expression(paste("Plastic count density (pieces/", km^2, ")"))
#y_axis<-expression(atop(atop("Rate of affected individuals"),atop("("~sqrt(proportion/year)~")")))
#y_axis<-expression(paste("Rate of affected individuals","(", sqrt(rate/year), ")"))
y_axis<-expression(atop(paste('Rate of affected individuals'), "("~sqrt(proportion/year)~")"))

ggplot(me, aes(x=x, y=predicted))+ 
  geom_line(col="black")+
  geom_point(data, mapping = aes(x=weighted_mean_count, y=sqrt(percent_year)), 
             size=2, col="darkcyan")+
  geom_ribbon(aes(ymin= conf.low, ymax=conf.high), alpha=0.35, fill="grey70", 
              show.legend = FALSE)+ 
  labs(y=y_axis, x=x_axis)+
  #labs(subtitle=text)+
  theme_bw()+
  scale_x_continuous(breaks = seq(from = 0.0, to = 2.5, by = 0.5))+
  scale_y_continuous(breaks = seq(from = 0.0, to = 1, by = 0.10))+
  theme(axis.title = element_text (face = "plain", size = 20, color = "black"),
        axis.text.x = element_text(size= 15),
        axis.text.y = element_text(size= 15),
        title= element_text(size= 25))+
  ggtitle(my_y_title)

### plot for the model example Dermochelys coriacea entanglement
##polt count
##first change to lmer
##predicted values
library(lme4)
me<- ggpredict(lm_model_count_sqrt, terms = c("weighted_mean_count"), back.transform = FALSE)


my_y_title <- expression(paste("a) ", "Entanglement (", italic("Chelonia mydas"), ")"))
x_axis<-expression(paste("Plastic count density (pieces/", km^2, ")"))
#y_axis<-expression(atop(atop("Rate of affected individuals"),atop("("~sqrt(proportion/year)~")")))
#y_axis<-expression(paste("Rate of affected individuals","(", sqrt(rate/year), ")"))
y_axis<-expression(atop(paste('Rate of affected individuals'), "("~sqrt(proportion/year)~")"))

ggplot(me, aes(x=x, y=predicted))+ 
  geom_line(col="black")+
  geom_point(data, mapping = aes(x=weighted_mean_count, y=sqrt(percent_year)), 
             size=2, col="darkcyan")+
  geom_ribbon(aes(ymin= conf.low, ymax=conf.high), alpha=0.35, fill="grey70", 
              show.legend = FALSE)+ 
  labs(y=y_axis, x=x_axis)+
  #labs(subtitle=text)+
  theme_bw()+
  scale_x_continuous(breaks = seq(from = 0.0, to = 2.5, by = 0.5))+
  scale_y_continuous(breaks = seq(from = 0.0, to = 1, by = 0.10))+
  theme(axis.title = element_text (face = "plain", size = 20, color = "black"),
        axis.text.x = element_text(size= 15),
        axis.text.y = element_text(size= 15),
        title= element_text(size= 25))+
  ggtitle(my_y_title)

##polt weight
##first change to lmer
#changge to lmer
lm_model_weight_sqrt<- lmer(sqrt(percent_year)~0+weighted_mean_weight+
                              (0+weighted_mean_weight|ocean),
                            data=data)
##predicted values
me_weight<- ggpredict(lm_model_weight_sqrt, terms = c("weighted_mean_weight"), back.transform = FALSE)
my_y_title <- expression(paste("a) ", "Ingestion (", italic("Chelonia mydas"), ")"))
x_axis<-expression(paste("Plastic count density (kg/", km^2, ")"))
#y_axis<-expression(atop(atop("Rate of affected individuals"),atop("("~sqrt(proportion/year)~")")))
#y_axis<-expression(paste("Rate of affected individuals","(", sqrt(rate/year), ")"))
y_axis<-expression(atop(paste('Rate of affected individuals'), "("~sqrt(proportion/year)~")"))


ggplot(me_weight, aes(x=x, y=predicted))+ 
  geom_line(col="black")+
  geom_point(data, mapping = aes(x=weighted_mean_weight, y=sqrt(percent_year)), 
             size=2, col="darkcyan")+
  geom_ribbon(aes(ymin= conf.low, ymax=conf.high), alpha=0.35, fill="grey70", 
              show.legend = FALSE)+ 
  labs(y=y_axis, x=x_axis)+
  #labs(subtitle=text)+
  theme_bw()+
  scale_x_continuous(breaks = seq(from = 0.0, to = 3, by = 0.5))+
  scale_y_continuous(breaks = seq(from = 0.0, to = 1.2, by = 0.10))+
  theme(axis.title = element_text (face = "plain", size = 20, color = "black"),
        axis.text.x = element_text(size= 15),
        axis.text.y = element_text(size= 15),
        title= element_text(size= 25))+
  ggtitle(my_y_title)

##annual rate avarages




####example for maps of figure EC10 main results
###maps EC110 with percentage

####Chelonia mydas
#1.165515871
rm(list=ls())
dev.off()
require(raster)
require(sp)
require(sf)
require(exactextractr)
require(rgeos)
require(tidyr)
require(dplyr)
require(MASS)
require(rgdal)
require(ggplot2)
require(maptools)
library("RColorBrewer")
#install.packages("sjlabelled")
require(sjlabelled)
count <- paste("../plastics/ing_count.tif", sep="")

data(wrld_simpl)

plastics_count <- raster(count)

####cut to the region to see the plastci raanges
zoom <- extent(-29,-2,-19, -4.5)
st_helena_map_ <- crop(plastics_count, zoom)
# making plot for hotspots example: EC = 0.600578852719709 for count 
chelonia_count <- plastics_count


#par(mfrow=c(1,2), mar=c(1,3,3,1))
st_helena_map <- chelonia_count

##change the 0.849346759 depending on the species values
cuts=c(0, 0.849346759, 1, 2, 3, 4, 5, 6,7.05) #set breaks 

##change depending of the green=safe zones
reds<- brewer.pal(7,"Reds")

#10^0.600578852719709= EC5=3.98
#labelss=c("0", "EC5=3.98", "10", "100", "1000","10000",
#         "100000","1000000", "10000000")
#10^0.600578852719709= EC5=3.98
expression(paste("Plastic count density (kg/", km^2, ")"))

#lb1=expression(paste(10^1, "/0.07"))
#lb2=expression(paste(10^2, "/0.10"))
#lb3=expression(paste(10^3, "/0.66"))
#lb4=expression(paste(10^4, "/1"))
#lb5=expression(paste(10^5, "/1"))
#lb6=expression(paste(10^6, "/1"))
#lb7=expression(paste(10^7, "/1"))
#lbEC5=expression(paste(10^1.2, "/0.10"))

lb1=expression(10^1)
lb2=expression(10^2)
lb3=expression(10^3)
lb4=expression(10^4)
lb5=expression(10^5)
lb6=expression(10^6)
lb7=expression(10^7)
lbEC5=expression(10^0.8) ##modify depending each Sp EC10


labelss=c("0", lbEC5,lb1, lb2,lb3,lb4,lb5,lb6,lb7)


st_helena <- extent(c(-19,-1,-20.5, -3.5))

########
proj4string(st_helena_map)
proj4string(wrld_simpl) # might need to change to same as st_helena map
target_crs <- CRS("+proj=longlat +datum=WGS84 +no_defs")  # This is your string assigned to an object
wrld_simpl2<- spTransform(wrld_simpl, target_crs)
####
pts = data.frame(
  lon = c(-5.7089,
          -14.365076),
  lat = c(-15.9650,
          -7.942932))

conversion_buffers <- st_as_sf(pts,coords=c("lon","lat"))

sf::sf_use_s2(FALSE)
buffer1<- st_buffer(conversion_buffers, dist=3.337) #distance now in arc decimals 200 nautical miles
#also 3.337 degrees km  370.4
projection_polygons1<- st_make_valid(buffer1)
validity<-st_is_valid(projection_polygons1) # checks whether polygons are valid 
nrow(validity)-sum(validity) #should be zero 
projection_polygons3<-sf::st_as_sf(projection_polygons1)

zoom <- extent(-19,-1,-20.5, -3.5)


###### zoomed

st_helena_map_ <- crop(st_helena_map, zoom)


zoom <- extent(-19,-1,-20.5, -3.5)



##Add bathymetry
#install.packages("marmap")
library(marmap)
##-18,-10.5,-11.8, -4
# get bathymetry data
b = getNOAA.bathy(lon1 = -19, lon2 = -1, lat1 = -20.5, lat2 = -3.5, 
                  resolution = 1)
## Querying NOAA database ...
## This may take seconds to minutes, depending on grid size
## Building bathy matrix ...




#install.packages("oce")
#install.packages("ocedata")
library(oce)
library(ocedata)
data("coastlineWorldFine")

# convert bathymetry
bathyLon = as.numeric(rownames(b))
bathyLat = as.numeric(colnames(b))
bathyZ = as.numeric(b)
dim(bathyZ) = dim(b)


###add names of the mounts in the plot (maybe change the color from gray tp other)


# plot coastline (no projection)
#plot(coastlineWorldFine, clon = mlon, clat = mlat, span = span)

###nmes of the island and sea mounts
##https://www.ascension.gov.ac/project/seamounts-project

#hey<-paste("Grattan Seamount")
#wow<- paste("Harris-Stewart Seamount")
#we<- paste("Young Seamount")
#hey1<-paste("Bonaparte Seamount")

hey<-paste("B")
wow<- paste("A")
we<- paste("C")
hey1<-paste("D")
St_Helena<-paste("St. Helena\n Island")

text1 <- data.frame(name=c("Ascension Island", wow, hey, we,
                           St_Helena, hey1, "E"),
                    lon = c(-14.5, -16.7, -13.3, -12.1, -5.7, -7, -6),
                    lat = c(-7.4, -8.5, -9.6, -9.3, -17, -15, -13.5))

#text1 <- data.frame(name=c("Ascension Island",
#                           "St. Helena"),
#                          lon = c(-14.35, -5.7),
#                         lat = c(-7.0, -16.9))


cntr <- c('Saint Helena') 
my_map <- wrld_simpl[wrld_simpl$NAME %in% cntr,]



##########plots
#####Global
png(
  "minPercentageEC10Chelonia_ing_globe.png",
  width     = 11.69,
  #11.69
  height    = 8.27,
  #8.27
  units     = "in",
  res       = 1200,
  pointsize = 2
)

par(oma=c(5,0,5,0))
par(mar=c(0,0,0,0))

## change the colours depending in the number of greens and red for each sp
plot(st_helena_map, breaks=cuts, lab.breaks=labelss, col = c ("green4",reds), 
     legend= FALSE, axes = FALSE, box=FALSE)
#cex=8, axis.args=list(cex.axis=8), legend.width= 4, cex.lab =4, 
#cex.axis = 4,

plot(wrld_simpl2,add=TRUE, col="black", border="black", axes = FALSE, box=FALSE)
#, axes = FALSE, box=FALSE

plot(projection_polygons3, add=TRUE, col=NA, border="white", lwd=3) #can change to whatever
plot(zoom, add=TRUE, col="white", lwd=3)



#raster::scalebar(500, divs=1, type="bar", cex=8, font=2, below = "km", 
#                lonlat = TRUE, label = c(0,500), 
#               xy = c(extent(st_helena_map)[1]+20,
#                     extent(st_helena_map)[3]+20), col= "white")


dev.off()


png(
  "minPercentageEC10Chelonia_ing_zoom.png",
  width     = 11.69,
  #11.69
  height    = 8.27,
  #8.27
  units     = "in",
  res       = 1200,
  pointsize = 2
)

par(oma=c(5,0,5,75))
par(mar=c(0,0,0,0))
#for the letters to change
#par(oma=c(10,10,5,65))
#par(mar=c(0,0,0,0))

## change the colours depending in the number of greens and red for each sp
plot(st_helena_map_, breaks=cuts, lab.breaks=labelss, col = c("green4", reds), 
     legend= FALSE, cex.axis = 8, axes = FALSE, box=FALSE)

#, axes = FALSE, box=FALSE

#Subsetting for the relevant polygons
plot(projection_polygons3, add=TRUE, col=NA, border="black", lwd=3) #can change to whatever

## change the colours depending in the number of greens and red for each sp
plot(st_helena_map, breaks=cuts, lab.breaks=labelss, col = c ("green4", reds), 
     legend= TRUE, cex=0.5, legend.width= 8,legend.shrink= 1,
     legend.only=TRUE, axis.args=list(cex.axis=10),
     add=TRUE, axes = FALSE, box=FALSE)

# plot bathymetry
contour(bathyLon,bathyLat,bathyZ,
        levels = c(-25,-50,-70, -100, -150, -200, -250, -300, -350, -400, -475,-515,-550,-600),
        lwd = c(1, 1,1,1, 1, 1, 1, 1,1,1,1,1,1),
        lty = c(1,1,1, 1, 1, 1, 1, 1,1,1,1,1,1),
        drawlabels = F, add = TRUE, col = 'yellow')
# lwd = c(1, 1,1,1, 2, 2, 3, 1,1,1,1,1,1),
#lty = c(3,1,1, 1, 3, 1, 3, 1,1,1,1,1,1),
##plot island
contour(bathyLon,bathyLat,bathyZ,
        levels= c(0, 25, 50,75, 100, 125,150,175,200,225,250,275,300),
        lwd = c(1, 1, 1, 1, 1, 1,1,1,1,1,1,1,1),
        lty = c(1, 1, 1, 1, 1, 1,1,1,1,1,1,1,1),
        drawlabels = F, add = TRUE, col = 'black')

plot(my_map,add=TRUE, col="black", border="black")


#text1$lon<-as.numeric(unlist(text1$lon))
#text1$latn<-as.numeric(unlist(text1$lat))
#text1 <- st_as_sf(st_h,coords=c("lon","lat"))


text(x = text1$lon, y = text1$lat, 
     text1$name, col=c('black'), cex=12, font=2)

#raster::scalebar(500, divs=1, type="bar", cex=8, font=2, below = "km", 
#        lonlat = TRUE, label = c(0,500), 
#       xy = c(extent(st_helena_map_)[1]+2,
#             extent(st_helena_map_)[3]+2))

raster::scalebar(400, divs=2, type="bar", cex=9, font=2, below = "km", 
                 lonlat = TRUE, label = c(0,200, 400), 
                 xy = c(extent(st_helena_map_)[1]+2,
                        extent(st_helena_map_)[3]+3))



dev.off()

###similar code for each species and type of interation remember entanglement
#is EC5 and ingestion EC10, the vaues of the estmate, CI lower ans max are
##different even if same species.