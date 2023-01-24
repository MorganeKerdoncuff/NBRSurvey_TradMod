# TradMod WP2 - NBR survey cleaning script
#Description of the data
#Date
#Who
#Project
#Funding
#Place

####Package loading####

library(tidyverse) #R language
library(readxl) #read xl files
library(lubridate) #standard date data


#### DATA LOADING ####

siteinfo_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="SiteInfo") #Farm information, field location and habitat type
landuse_raw <- read_excel(path = "data/NBR_FarmerSurvey.xlsx", sheet="Farms Information_R") #Field management data from farmer interview, at site level
#farmer_raw <- read_excel(path = "data/NBR_FarmRepertoire.xlsx") #Farm management data from farmer interview med Margit, at site level
landscape_raw <- read_excel(path = "data/NBR_LandscapeMatrix.xlsx", sheet="MatrixProportion") #Land cover data around the fields from Geonorge, at site level
area20x20_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="20mX20m") #Sampling area description, vegetation cover at the site level
vgcover_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="SoilCover") #Vegetation cover at the quadrat (subplot) level
soilpene_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="SoilPenetration") #Penetration rate in the soil, at subplot level (two collection per subplot)
soilbulk_raw <- read_excel(path = "data/NBR_RawBD2019-2020.xlsx", na="NA") #Soil bulk density, at subplot level (three samples per subplot)
soilmeso_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="Mesofauna") #Height of mesofauna soil core, at subplot level
chem2019 <- read.csv(path = "data/SoilChemistry_NBR2019.txt", sep=";") #2019 Soil chemistry data from Eurofins, at plot level
chem2020 <- read.csv(path = "data/SoilChemistry_NBR2020.txt", sep=";") #2020 soil chemistry data from EUrofins, at plot level
poo_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="Poo", na="NA") #poo data at subplot (quadrat) level
vege_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="PlantRichness") #plant community data, at species level
arthro_raw <- read_excel(path = "data/NBR_RawArthro2019-2020.xlsx", na="NA") #arthropod community data, at family level


#### SITE INFO ####

## Description

## List of variables

# [1] Field identification code for data collection
# [2] Former ecological zonation (adapted from NIBIO) - obsolete, not included in analysis
# [3] Validated ecological zonation (adapted from NIBIO)
# [4] Year of sampling
# [5] Geographical location of the field, hamlet
# [6] Geographical location of the field, municipality
# [7] Geolocation of the field, X-coordinate - CRS EPSG 25832
# [8] Geolocation of the field, y-coordinate - CRS EPSG 25832
# [9] Name of the farmer !! Privacy !!?
# [10] Current grazing livestock, at least for the last 5 years
# [11] Size of the livestock flock - column empty, data collected in another file
# [12] Size of the grazing area - column empty, data collected in another file
# [13] Y/N if the animal was seen on site during the collection
# [14] Poo weight collected on 10% of the sampling area (g) - proxy of grazing intensity on site
# [15] Comments

#
## Summary site info - Check table size, list of variables, variable types (num/chr)

str(siteinfo_raw) # All good

#
## Name & character cleaning Site info

names(siteinfo_raw) <- gsub("%", "percent", names(siteinfo_raw))
names(siteinfo_raw) <- gsub(" ", "", names(siteinfo_raw)) 
names(siteinfo_raw) <-  gsub("\\(", "_", names(siteinfo_raw))
names(siteinfo_raw) <-  gsub("\\)", "", names(siteinfo_raw))
siteinfo_raw <- subset(siteinfo_raw, select = -c(Size_livestock, Surface)) # remove empty variables, which will be included in another dataset

#
## Data cleaning - new R object

siteinfo_full <- siteinfo_raw

#
## Numeric var Site Info - Check min/max, distribution and potential outliers

# Check min/max
test <- siteinfo_full %>% 
  summarise(tibble(across(where(is.numeric), ~min(.x, na.rm = TRUE), .names = "min_{.col}"),across(where(is.numeric), ~max(.x, na.rm = TRUE), .names = "max_{.col}"))) %>% 
  transpose() # All good

# Check distribution of quantitative variable
#hist(siteinfo_full$`Pooestimation_g-10percent`) # Standard distribution of poo data - validated

#
## Char var Site Info - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(siteinfo_full$SiteID) # Unique ID for each site - validated

# Ecological zones
#unique(siteinfo_full$NiBioEcologicalZone) # Three categories for ecological zone - validated

# Geographical locations
#unique(siteinfo_full$Location) # Correct locations and location names

# Muncipalities
#unique(siteinfo_full$Municipality) # missing data + US2 Stordalen -> did you guys went that far ? Yes validated
#siteinfo_full[is.na(siteinfo_full$Municipality),] # NA identified - US2, Stordalen is located in Masfjorden
siteinfo_full <- siteinfo_full |> 
  mutate(Municipality = ifelse(SiteID == "US2", "Masfjorden", Municipality)) # Replace NA by municipality name
#unique(siteinfo_full$Municipality) # NA in Municipality replaced - validated

# Farmers
#unique(siteinfo_full$Farmer) # Missing two full names and one NA
#siteinfo_full[is.na(siteinfo_full$Farmer),] # NA identified - UC1, farmer same as for US1
siteinfo_full <- siteinfo_full |> 
  mutate(Farmer = ifelse(SiteID == "UC1", "Hans Magne Haukeland", Farmer)) |> 
  mutate(Farmer = dplyr::recode(Farmer, "Tormod" = "Tormod Magnesen")) |> 
  mutate(Farmer = dplyr::recode(Farmer, "Dan" = "Ole Mathias Lygren")) # Full farmer names - validated

# Livestock
#unique(siteinfo_full$Type_livestock) # Villsau should be within sheep category
siteinfo_full <- siteinfo_full |> 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "Villsau" = "Sheep")) # Only 3 livestock categories - validated
#table(siteinfo_full$Type_livestock) # correct number of sites per livestock category - validated

# Habitat type
# unique(siteinfo_full$Habitat) # 3 categories + one extra (bog) which will not be considered in the analysis
# table(siteinfo_full$Habitat) # habitat repartition - validated

#unique(siteinfo_full$Animal_on_site)

#
## Geological information # Aborted as raw data were incorrect


#### FIELD MANAGEMENT ####

## Description

## List of variables

# [1] Field identification code for data collection
# [2] Geographical location of the field, hamlet
# [3] Geographical location of the field, postcode
# [4] Geographical location of the field, municipality (both former and current classification)
# [5] Date of the interview with the farmer
# [6] Site productivity - infield high productivity, outfield low productivity
# [7] National identification number of the farm
# [8] Main livestock, currently grazing in the field - ! rotational grazing management
# [9] Other livestock, grazing in other fields - ! rotational grazing management
# [10] Name of cow breed
# [11] Name of goat breed
# [12] Name of sheep breed
# [13] Number of adult animals in the main livestock
# [14] Number of young animals in the main livestock
# [15] Number of adult animals in other livestock
# [16] Number of young animals in other livestock
# [17] Size of the field containing the sampling area (ha) ! rotational grazing management
# [18] FROM GÅRDSKSART - Total grazing area of the farm (ha) ! rotational grazing management
# [19] Grazing density on the field - not collected, empty column
# [20] Period since the farmer has had the current livestock
# [21] Period without grazing management on the field
# [22] Number of months the main livestock grazes in the field in a year
# [23] Number of months other livestock graze in the field in a year ! rotational grazing
# [24] If the animals are kept inside or outside at night
# [25] Farmer impression of grazing pressure on the site
# [26] Former livestock which used to graze in the field (1)
# [27] Former livestock which used to graze in the field (2)
# [28] Former livestock which used to graze in the field (3)
# [29] Former livestock grazing period (1)
# [30] Former livestock grazing period (2)
# [31] Former livestock grazing period (3)
# [32] Type of farm management during the past 10 years (1)
# [33] Type of farm management during the past 10 years (2)
# [34] Type of farm management during the past 10 years (3)
# [35] If applicable, frequency of grass cutting
# [36] If applicable, frequency of mulching
# [37] If applicable, frequency of fertilization with manure
# [38] If applicable, type of manure used
# [39] If applicable, volume of manure used (m3)
# [40] If applicable, in which season the manure is used
# [41] If applicable, frequency of fertilization with artificial fertilizer
# [42] If applicable, weight of artificial fertilizer used (kg)
# [43] If applicable, in which season the artificial fertilizer is used
# [44] If applicable, frequency of fertilization with shell-sand/lime
# [45] If applicable, weight of shell-sand/limer used (kg)
# [46] If applicable, in which season the shell-sand/lime is used
# [47] If applicable, last time the field was sowed
# [48] If applicable, last time the field was plowed
# [49] If applicable, last time the field was drained
# [50] Type of land use prior to grazing
# [51] If applicable, type of production prior to grazing (1)
# [52] If applicable, type of production prior to grazing (2)
# [53] Time period of land use type prior to grazing
# [54] Comments


#### SAMPLING AREA 20X20 ####

## List of variables

# [1] Field identification code for data collection
# [2] Date of data collection
# [3] Names of observers
# [4] Geolocation of corner 1 of the sampling area, Y-coordinate - CRS WPSG84
# [5] Geolocation of corner 1 of the sampling area, X-coordinate - CRS WPSG84
# [6] Geolocation of corner 3 of the sampling area, Y-coordinate - CRS WPSG84
# [7] Geolocation of corner 3 of the sampling area, X-coordinate - CRS WPSG84
# [8] Number of paths in the vegetation created by the animals
# [9] Sum of the lengths of the animal paths (m)
# [10] Elevation at the lowest point of the sampling area (m)
# [11] Elevation at the highest point of the sampling area (m)
# [12] Estimated slope from the highest to the lowest point of the sampling area (degree)
# [13] Aspect of the sampling area - cardinal direction
# [14] Aspect of the sampling area - azimuth degree
# [15] Estimated percentage cover of exposed rock in the sampling area
# [16] Estimated percentage cover of mud in the sampling area
# [17] Estimated percentage cover of trees and shrubs over 1 m in the sampling area
# [18] Estimated percentage cover of shrubs under 1 m in the sampling area
# [19] Estimated percentage cover of herbs in the sampling area
# [20] Estimated percentage cover of monocotyledons (grasses, rushes, sedges) in the sampling area
# [21] Estimated percentage cover of bryophytes (mosses, liverworts) in the sampling area
# [22] Estimated percentage cover of ferns in the sampling area
# [23] Estimated percentage cover of lichens in the sampling area
# [24] Comments

#
## Summary sampling area - Check table size, list of variables, variable types (num/chr)

str(area20x20_raw) # All good, date should be reformatted

#
## Name & character cleaning sampling area

names(area20x20_raw)<-gsub("%", "percent", names(area20x20_raw)) # remove percent signs from names
names(area20x20_raw)<-gsub("&", "_", names(area20x20_raw)) # remove &
names(area20x20_raw)<-gsub("\\(", "", names(area20x20_raw)) # remove (
names(area20x20_raw)<-gsub("\\)", "", names(area20x20_raw)) # remove &

#
## Sampling date standardisation

area20x20_raw$Recording_date <- as.POSIXct(area20x20_raw$Recording_date, format = "%d.%m.%Y")

#
## Data cleaning - New R object

area20x20_full <- area20x20_raw

#
## Char var Site Info - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(area20x20_full$SiteID) # Unique ID for each site - validated

# Aspect
#unique(area20x20_full$Aspect) # One missing value under "/"
#area20x20_full[area20x20_full$Aspect == "/",] # UC1 -> flat bog, will not be used
area20x20_full <- area20x20_full |> 
  mutate(Aspect = dplyr::recode(Aspect, "/" = "NA"))

#
## Numeric var sampling area - Check min/max, distribution and potential outliers

# Check min/max
test <- area20x20_full |> 
  summarise(tibble(across(where(is.numeric), ~min(.x, na.rm = TRUE), .names = "min_{.col}"),across(where(is.numeric), ~max(.x, na.rm = TRUE), .names = "max_{.col}"))) |> 
  transpose() # potential outliers are max herbs 97% and max lichens 20%

# Distribution elevation max
#hist(area20x20_full$Elevation_max) # Sites range from 0 to 900 m elevation, no outlier - validated

# Distribution slope
#hist(area20x20_full$General_slope) # Slopes range from 0 to 35 degrees, no outlier - validated

# Distribution aspect degree
#hist(area20x20_full$AspectDegree) # Aspect covers all spectrum (0 to 360), no outliers - validated

# Distribution rock cover
#hist(area20x20_full$percentRock) # Some sites over 7%, check their location + 3 NAs
#area20x20_full[area20x20_full$percentRock>7,] # 5 sites OV1, UG2, US3, IG3, US5
# UG2, US3 and US5 in subalpine areas with exposed bedrock - validated
# OV1 in coastal heathland habitat, with exposed bedrock - validated
# IG3 in fjord area but at higher elevation, with exposed bedrock - validated
test <- area20x20_full[is.na(area20x20_full$percentRock),] # 3 sites US1, UG1 and OC4 missing all soil cover percentages

# Distribution mud cover
#hist(area20x20_full$percentMud) # One site over 15%
#area20x20_full[area20x20_full$percentMud>10,] # UC1 -> bog site, will not be included into the analysis - validated

# Distribution trees and tall shrubs cover
#hist(area20x20_full$percentTrees_TallShrubs) # No sites over 10% - validated

# Distribution low shrubs
#hist(area20x20_full$percentLowShrubs) # Wide range due to collection in both grassland and heathland habitats. Check that all grassland sites are under 10%
#area20x20_full[area20x20_full$percentLowShrubs>10,] # 11 sites over 10%
# OV1, OV2, OS5, OS7, OS9, IS2 coastal heathlands - validated
# US2, US3, US4, US5, UG2 subalpine heathlands - validated

# Distribution herbs
hist(area20x20_full$percentHerbs) # a few sites over 60%
area20x20_full[area20x20_full$percentHerbs>60,] # OS1, OC1, IG1, IS2, IC1 -> all first year/starting sites, check on vegetation quadrats + site & plot pictures
# OS1 80% - average 20% & no cover over 55% in quadrats, estimation from pictures 35%-40%
# OC1 80% - average 25% & no cover over 50% in quadrats, estimation from pictures 10%-15% 
# IG1 97% - average 20% & no cover over 30% in quadrats, estimation from pictures 15%-20%
# IS2 80% - average 70% & no cover over 90% in quadrats, estimation from pictures 45%-50%
# IC1 70% - average 45% & no cover over 70% in quadrats, estimation from pictures 60%-65%

# Distribution monocotyledons
hist(area20x20_full$percentMonocotyledons) # need to check again sites with weird herb estimations

# Distribution bryophytes
hist(area20x20_full$percentBryophytes) # need to check again sites with weird herb estimations

# Distribution lichens
#hist(area20x20_full$percentLichen) # one site over 10%
#area20x20_full[area20x20_full$percentLichen>15,] # US4 in subalpine area, average of 12% & max 24% in quadrats - validated

#hist(area20x20_full$percentBryophytes)

#Char var 20x20
#unique(area20x20_full$Recording_date)#see all possible values
#table(area20x20_full$SiteID)# no site name appears twice :-)
#unique(area20x20_full$Aspect)# a bit of a mess.

#Heat Load Index
area20x20_full <- area20x20_full |> 
  mutate(HLI = cos(AspectDegree-225)*tan(General_slope))
#hist(area20x20_full$HLI) # 3 outliers: one under 200, two over 100
#area20x20_full[area20x20_full$HLI>100,] #OG4 & IS3 -> both 11 degree slope with SW & SE exposition
#area20x20_full[area20x20_full$HLI<0,] #OS6 -> 11 degree slope with NE exposition

#Last check -> 45 sites
#area20x20_full <- subset(area20x20_full, !SiteID=="UC1")
#unique(area20x20_full$SiteID)
sort(table(area20x20_full$SiteID))
