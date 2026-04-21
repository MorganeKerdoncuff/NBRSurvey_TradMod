##
#### DESCRIPTION ####
##
## Purpose: Cleaning script for secondary data collected for the GrazingFjordScale paper
## Author: Morgane KERDONCUFF
## ORCID: 0000-0003-2223-1857
## github
## Date created: 04/2026
## Last time modified: 04/2026
## Project: TradMod
## Funding: Norges forskningsråd (NFR)
## Institution: University of Bergen, Norway
##
##

#### PACKAGES ####

library(tidyverse) #R language
library(janitor) #Data cleaning
library(readxl) #read xl files
library(lubridate) #standard date
library(terra) #raster files
library(purrr) #file merging

#### DATA ####

# Site description
site_description <- read.csv("data/cleandata/tradmod_nbr_sitedescription.csv", sep=",")
# 1 km^2^ land cover, from Geonorge (https://kartkatalog.geonorge.no/metadata/fkb-ar5/166382b4-82d6-4ea9-a68e-6fd0c87bf788)
landscape_raw <- read_excel(path = "data/rawdata/NBR_RawLandscapeMatrix.xlsx")

#### REGIONAL SCALE ####

#### Climate data ####

#
## Description

# Import and transformation of daily temperature and precipitation data from MET into yearly averages
# Type of data - 1 km x 1km
# Variables
## Daily precipitation
## Daily maximum temperature
## Daily minimum temperature

#
## Prepare coordinates

# Subset from site info file + simplification col names
sitecoord <- subset(siteinfo_full, select = c(SiteID, EPSG.25832_X, EPSG.25832_Y))
names(sitecoord) <- gsub("EPSG.25832_X", "Xsite", names(sitecoord))
names(sitecoord) <- gsub("EPSG.25832_Y", "Ysite", names(sitecoord))

# Create column id which will allow binding with raster analysis results
sitecoord <- sitecoord |> 
  mutate(ID = 1:n()) #|> 
#select(ID, everything())

#
## Precipitation data

# Make function

# Import tiff files according to year
preci_list2010 <- list.files(path = "./data/precipitationraster", pattern= '2010.*\\.tif$', all.files=TRUE, full.names=TRUE) # pattern equal to : contains exactly '2010', then 0 or x characters, then exactly '.tif' and nothing after
preci_list2011 <- list.files(path = "./data/precipitationraster", pattern= '2011.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2012 <- list.files(path = "./data/precipitationraster", pattern= '2012.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2013 <- list.files(path = "./data/precipitationraster", pattern= '2013.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2014 <- list.files(path = "./data/precipitationraster", pattern= '2014.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2015 <- list.files(path = "./data/precipitationraster", pattern= '2015.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2016 <- list.files(path = "./data/precipitationraster", pattern= '2016.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2017 <- list.files(path = "./data/precipitationraster", pattern= '2017.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2018 <- list.files(path = "./data/precipitationraster", pattern= '2018.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2019 <- list.files(path = "./data/precipitationraster", pattern= '2019.*\\.tif$', all.files=TRUE, full.names=TRUE)
preci_list2020 <- list.files(path = "./data/precipitationraster", pattern= '2020.*\\.tif$', all.files=TRUE, full.names=TRUE)

# Import raster files within one object
#preci_daily2010 <- stack(preci_list2010) # base raster, slower
preci_daily2010 <- terra::rast(preci_list2010) #faster, but needs following terra method and not standard raster
preci_daily2011 <- terra::rast(preci_list2011)
preci_daily2012 <- terra::rast(preci_list2012)
preci_daily2013 <- terra::rast(preci_list2013)
preci_daily2014 <- terra::rast(preci_list2014)
preci_daily2015 <- terra::rast(preci_list2015)
preci_daily2016 <- terra::rast(preci_list2016)
preci_daily2017 <- terra::rast(preci_list2017)
preci_daily2018 <- terra::rast(preci_list2018)
preci_daily2019 <- terra::rast(preci_list2019)
preci_daily2020 <- terra::rast(preci_list2020)

# Calculate yearly value
preci_yearly2010 <- terra::app(preci_daily2010, sum, NA.RM = TRUE)
#preci_mean2010 <- calc(preci_daily2010, mean, NA.RM = TRUE)
#plot(preci_yearly2010)
preci_yearly2011 <- terra::app(preci_daily2011, sum, NA.RM = TRUE)
preci_yearly2012 <- terra::app(preci_daily2012, sum, NA.RM = TRUE)
preci_yearly2013 <- terra::app(preci_daily2013, sum, NA.RM = TRUE)
preci_yearly2014 <- terra::app(preci_daily2014, sum, NA.RM = TRUE)
preci_yearly2015 <- terra::app(preci_daily2015, sum, NA.RM = TRUE)
preci_yearly2016 <- terra::app(preci_daily2016, sum, NA.RM = TRUE)
preci_yearly2017 <- terra::app(preci_daily2017, sum, NA.RM = TRUE)
preci_yearly2018 <- terra::app(preci_daily2018, sum, NA.RM = TRUE)
preci_yearly2019 <- terra::app(preci_daily2019, sum, NA.RM = TRUE)
preci_yearly2020 <- terra::app(preci_daily2020, sum, NA.RM = TRUE)

# Check if all sites covered
#precipitation2010 <- as.data.frame(terra::extract(preci_yearly2010, subset(siteinfo_full, select = c(EPSG.25832_X, EPSG.25832_Y))))

# Fill missing cell with average from closest cells
fullpreci_yearly2010 <- terra::focal(preci_yearly2010, w=7, fun = "mean", na.policy = "only")
#plot(fullpreci_yearly2010)
fullpreci_yearly2011 <- terra::focal(preci_yearly2011, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2012 <- terra::focal(preci_yearly2012, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2013 <- terra::focal(preci_yearly2013, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2014 <- terra::focal(preci_yearly2014, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2015 <- terra::focal(preci_yearly2015, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2016 <- terra::focal(preci_yearly2016, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2017 <- terra::focal(preci_yearly2017, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2018 <- terra::focal(preci_yearly2018, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2019 <- terra::focal(preci_yearly2019, w=7, fun = "mean", na.policy = "only")
fullpreci_yearly2020 <- terra::focal(preci_yearly2020, w=7, fun = "mean", na.policy = "only")

# Average annual precipitation from 2010 to 2020
list_allpreci_yearly <- terra::as.list(fullpreci_yearly2010, fullpreci_yearly2011, fullpreci_yearly2012, fullpreci_yearly2013, fullpreci_yearly2014, fullpreci_yearly2015, fullpreci_yearly2016, fullpreci_yearly2017, fullpreci_yearly2018, fullpreci_yearly2019, fullpreci_yearly2020)
allpreci_yearly <- terra::rast(list_allpreci_yearly)
averagepreci <- terra::app(allpreci_yearly, mean, NA.RM = TRUE)

# Extract yearly values for NBR site coordinates
precipitation <- as.data.frame(terra::extract(averagepreci, subset(sitecoord, select = c(Xsite, Ysite))))
names(precipitation) <- gsub("mean", "annualprecipitation", names(precipitation))

#
## Max temperature data

# Import tiff files according to July month
maxtemp_listJuly <- list.files(path = "./data/maxtempraster", pattern= '_07_.*\\.tif$', all.files=TRUE, full.names=TRUE)

# Import raster files within one object
maxtemp_dailyJuly <- terra::rast(maxtemp_listJuly)

# Calculate mean
maxtemp_meanJuly <- terra::app(maxtemp_dailyJuly, mean, NA.RM = TRUE)
#plot(maxtemp_meanJuly)

# Fill missing cell with average from closest cells
fullmaxtemp_meanJuly <- terra::focal(maxtemp_meanJuly, w=7, fun = "mean", na.policy = "only")
plot(fullmaxtemp_meanJuly)

# Extract mean values for NBR site coordinates
maxtempjuly <- as.data.frame(terra::extract(fullmaxtemp_meanJuly, subset(sitecoord, select = c(Xsite, Ysite))))
names(maxtempjuly) <- gsub("focal_mean", "maxtempJuly", names(maxtempjuly))

#
## Min temperature data

# Import tiff files according to January month
mintemp_listJan <- list.files(path = "./data/mintempraster", pattern= '_01_.*\\.tif$', all.files=TRUE, full.names=TRUE)

# Import raster files within one object
mintemp_dailyJan <- terra::rast(mintemp_listJan)

# Calculate mean
mintemp_meanJan <- terra::app(mintemp_dailyJan, mean, NA.RM = TRUE)
#plot(mintemp_meanJan)

# Fill missing cell with average from closest cells
fullmintemp_meanJan <- terra::focal(mintemp_meanJan, w=7, fun = "mean", na.policy = "only")
plot(fullmintemp_meanJan)

# Extract mean values for NBR site coordinates
mintempJan <- as.data.frame(terra::extract(fullmintemp_meanJan, subset(sitecoord, select = c(Xsite, Ysite))))
names(mintempJan) <- gsub("focal_mean", "mintempJan", names(mintempJan))

#
## Average July temperature data

# Import tiff files according to July month
avgtemp_listJuly <- list.files(path = "./data/avgtempraster", pattern= '_07_.*\\.tif$', all.files=TRUE, full.names=TRUE)

# Import raster files within one object
avgtemp_dailyJuly <- terra::rast(avgtemp_listJuly)

# Calculate mean
avgtemp_meanJuly <- terra::app(avgtemp_dailyJuly, mean, NA.RM = TRUE)
#plot(avgtemp_meanJuly)

# Fill missing cell with average from closest cells
fullavgtemp_meanJuly <- terra::focal(avgtemp_meanJuly, w=7, fun = "mean", na.policy = "only")
plot(fullavgtemp_meanJuly)

# Extract mean values for NBR site coordinates
avgtempJuly <- as.data.frame(terra::extract(fullavgtemp_meanJuly, subset(sitecoord, select = c(Xsite, Ysite))))
names(avgtempJuly) <- gsub("focal_mean", "avgtempJuly", names(avgtempJuly))

#
## Average Jan temperature data

# Import tiff files according to Jan month
avgtemp_listJan <- list.files(path = "./data/avgtempraster", pattern= '_01_.*\\.tif$', all.files=TRUE, full.names=TRUE)

# Import raster files within one object
avgtemp_dailyJan <- terra::rast(avgtemp_listJan)

# Calculate mean
avgtemp_meanJan <- terra::app(avgtemp_dailyJan, mean, NA.RM = TRUE)
#plot(avgtemp_meanJan)

# Fill missing cell with average from closest cells
fullavgtemp_meanJan <- terra::focal(avgtemp_meanJan, w=7, fun = "mean", na.policy = "only")
plot(fullavgtemp_meanJan)

# Extract mean values for NBR site coordinates
avgtempJan <- as.data.frame(terra::extract(fullavgtemp_meanJan, subset(sitecoord, select = c(Xsite, Ysite))))
names(avgtempJan) <- gsub("focal_mean", "avgtempJan", names(avgtempJan))

#
## Preparation final dataset

# Merging all climate results
climate_full <- purrr::reduce(list(sitecoord, precipitation, maxtempjuly, mintempJan, avgtempJuly, avgtempJan), dplyr::left_join)

# Rescaling - original data at 0.1 scale - temp conversion to C
climate_full <- climate_full |> 
  mutate(annualprecipitation = annualprecipitation/10,
         maxtempJuly = maxtempJuly/100-273.15,
         mintempJan = mintempJan/100-273.15,
         avgtempJuly = avgtempJuly/10-273.15,
         avgtempJan = avgtempJan/10-273.15)

# Clean data + export in csv
climate_full <- subset(climate_full, select = -c(ID, Xsite, Ysite))
write_csv(climate_full, "data/cleandata/NBR_FullClimate.csv")

# Export summary raster files
writeRaster(averagepreci, "outputs/NBR_maps/NBR_AnnualPreci.tiff", overwrite = TRUE)
writeRaster(fullmaxtemp_meanJuly, "outputs/NBR_maps/NBR_MaxJulyTemp.tiff", overwrite = TRUE)
writeRaster(fullmintemp_meanJan, "outputs/NBR_maps/NBR_MinJanTemp.tiff", overwrite = TRUE)
writeRaster(fullavgtemp_meanJuly, "outputs/NBR_maps/NBR_AvgJulyTemp.tiff", overwrite = TRUE)
writeRaster(fullavgtemp_meanJan, "outputs/NBR_maps/NBR_AvgJanTemp.tiff", overwrite = TRUE)


#### SUMMARY FOR REPOSITORY ####

#
## Site-level data
# Merging of site info full, area 20x20, land use -> NBR_sitescale

# Select variables for repository
repo_site <- purrr::reduce(list(
  subset(siteinfo_full, select = c(SiteID, EcoZone, Year, Habitat)), 
  subset(landuse_full, select = c(SiteID, FieldType, Livestock1, FlockSize1_adults, FlockSize1_young, SelectedFieldArea_ha, FarmInfieldArea_ha, LSU_adult, LSU_young, AvgStockDensity_perha)), 
  subset(area20x20_full, select = -c(HeatLoadIndex, Comments_area20x20))
), dplyr::full_join)




#### LANDSCAPE SCALE ####

# Raw datasets
## landscape_raw

# Desired variables
## siteID - Field identification code for data collection
## fullyCultivatedLandHa - Area of cultivated lands within 1 km radius around the site
## superficiallyCultivatedLandHa - Area of cultivated lands within 1 km radius around the site
## infieldHa - Area of productive grazing fields within 1 km radius around the site
## productiveForestHa - Area of forests within 1 km radius around the site
## nonProductiveForestHa - Area of forests within 1 km radius around the site
## wetlandHa - Area of wetlands within 1 km radius around the site
## outfieldHa - Area of unproductive open lands within 1 km radius around the site
## freshwaterHa - Area of freshwater within 1 km radius around the site
## infrastructureHa - Area of built land and infrastructures within 1 km radius around the site
## seaHa - Area of sea bodies within 1 km radius around the site

# Variable names

## Lower camel with janitor package
landscape_raw <- landscape_raw %>% 
  clean_names("lower_camel")

## Consistent & FAIR
names(landscape_raw) <- gsub("siteId", "siteID", names(landscape_raw))

# Dataset
landscape <- landscape_raw

## Variable types (num/chr)
# str(landscape) #validated
## Duplicate check
# get_dupes(landscape) #validated

## Character variables - desirable categories, NAs, misprints

### Consistent lower camel
landscape <- landscape %>% 
  mutate_if(is.character, tolower)

### Categories & distribution
# table(landscape$siteID) # Unique ID for each site & only grassland sites

## Numeric variables - min/max, distribution, potential outliers

### Min/max
# test <- landscape |>
#   summarise(
#     tibble(
#       across(
#         where(is.numeric),
#         ~min(.x, na.rm = TRUE),
#         .names = "min_{.col}"
#         ),
#       across(
#         where(is.numeric),
#         ~max(.x, na.rm = TRUE),
#         .names = "max_{.col}")
#       )
#     ) |>
#   transpose() #validated

### NA check
# landscape[is.na(landscape$fullyCultivatedLandHa),] #validated
# landscape[is.na(landscape$superficiallyCultivatedLandHa),] #validated
# landscape[is.na(landscape$infieldHa),] #validated
# landscape[is.na(landscape$productiveForestHa),] #validated
# landscape[is.na(landscape$nonProductiveForestHa),] #validated
# landscape[is.na(landscape$wetlandHa),] #validated
# landscape[is.na(landscape$outfieldHa),] #validated
# landscape[is.na(landscape$freshwaterHa),] #validated
# landscape[is.na(landscape$infrastructureHa),] #validated
# landscape[is.na(landscape$seaHa),] #validated

### Variable distribution & outliers
# hist(landscape$fullyCultivatedLandHa) # Poisson distribution, no outlier
# hist(landscape$superficiallyCultivatedLandHa) # Poisson distribution, no outlier
# hist(landscape$infieldHa) # Poisson distribution, no outlier
# hist(landscape$productiveForestHa) # Poisson distribution (a bit uneven), no outlier
# hist(landscape$nonProductiveForestHa) # Poisson distribution (skewed), no outlier
# hist(landscape$wetlandHa) # Poisson distribution (skewed), no outlier
# hist(landscape$outfieldHa) # Poisson distribution, no outlier
# hist(landscape$freshwaterHa) # Poisson distribution (skewed), no outlier
# hist(landscape$infrastructureHa) # Poisson distribution, 2 sites slightly outlying above 60 ha but should not affect analysis
# filter(landscape, infrastructureHa > 60) # oc4 & og2
# hist(landscape$seaHa) # Poisson distribution, no outlier

# Dataset export
write_csv(landscape, "data/cleandata/tradmod_nbr_landscape.csv")

