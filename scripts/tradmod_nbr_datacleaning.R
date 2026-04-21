##
#### DESCRIPTION ####
##
## Purpose: Cleaning script for ecological data collected in 2019 and 2020 in the Nordhordland UNESCO Biosphere Reserve
## Author: Morgane KERDONCUFF
## ORCID: 0000-0003-2223-1857
## github
## Date created: 11/2022
## Last time modified: 03/2026
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

#### RAW DATA ####

# Site location & habitat type
site_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SiteInfo")
# Site management, from farmer interviews
management_raw <- read_excel(path = "data/rawdata/NBR_RawFarmerSurvey.xlsx", sheet="Farms Information_R")
# 1 km^2^ land cover, from Geonorge (https://kartkatalog.geonorge.no/metadata/fkb-ar5/166382b4-82d6-4ea9-a68e-6fd0c87bf788)
landscape_raw <- read_excel(path = "data/rawdata/NBR_RawLandscapeMatrix.xlsx")
# 20 m x 20 m sampling area
sampling_area_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="20mX20m")
# 1 m^2^ ground cover in non-destructive subplots (1 per subplot)
ground_cover_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SoilCover")
# Soil penetration tests in destructive subplots (2 per subplot)
soil_pene_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SoilPenetration")
# Bulk density in destructive subplots (3 per subplot)
soil_bulk_raw <- read_excel(path = "data/rawdata/NBR_RawBD.xlsx", na="NA")
# Soil chemistry 2019 (3 per plot)
soil_chem_2019 <- read.csv("data/rawdata/NBR_RawSoilChemistry2019.txt", sep=";")
# Soil chemistry 2020 (3 per plot)
soil_chem_2020 <- read.csv("data/rawdata/NBR_RawSoilChemistry2020.txt", sep=";")
# Soil chemistry 2020 complementary
soil_chem_2020_DM <- read_excel(path = "data/rawdata/NBR_RawSoilChemistry2020bis.xls")
# 1 m^2^ dung collection in non-destructive subplot
dung_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="Poo", na="NA")
# 1 m^2^ plant community, species level
plant_com_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="PlantRichness")
# Arthropod community, family level (beetles) and order level (others)
arthro_com_raw <- read_excel(path = "data/rawdata/NBR_RawArthro.xlsx", na="NA")
# Arthropod community complementary
arthro_com_sup <- read_excel(path = "data/rawdata/NBR_RawArthroSup.xlsx", na="NA")
# Mesofauna community
meso_com_raw <- read_excel(path = "data/rawdata/NBR_RawMesobio.xlsx", na="NA")
# Soil core height for mesofauna in destructive subplots (1 per subplot)
meso_soil_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="Mesofauna")
# Aboveground biomass in destructive subplots
plant_biomass_raw <- read_excel(path = "data/rawdata/NBR_RawAGB.xlsx", na="NA")

#### SITE DESCRIPTION ####

# Raw datasets
## site_raw
## management_raw

# Desired variables
## siteID - Field identification code for data collection
## samplingYear - Year of data collection (values: 2019; 2020)
## landscapeZone - Landscape zone (values: coastal; fjord; mountain) adapted from Uttakleiv, Lars A., and Trond Simensen. 2011. “Landskapskartlegging i Norge [Landscape Mapping in Norway].” Plan 43 (3-4): 76–83. https://doi.org/10.18261/ISSN1504-3045-2011-03-04-15.
## municipality - Smaller administrative region
## latitude25832 - Latitude in coordinate reference system ETRS89 / UTM zone 32N (EPSG 25832)
## longitude25832 - Longitude in coordinate reference system ETRS89 / UTM zone 32N (EPSG 25832)
## habitatType - Simplified habitat type according to shrub dominance, geographical location and land use (values: coastal heathland; permanent grassland; transition heathland; subalpine heathland)
## fieldType - Land use type according to traditional Norwegian farming system (values: infield; outfield)
## livestockType - Type of livestock grazing on field (values: cattle; goat; sheep)
## numberAnimalsAdult - Number of adult animals grazing on field during the year of sampling
## numberAnimalsYoung - Number of young animals grazing on field during the year of sampling (born in the year for sheep & goats, less than 1 year old for cattle)
## fieldAreaHa - Surface area of the field selected for the survey in hectares
## farmGrazingAreaHa - Total surface area of grazing fields (innmark) of the farm, collected from https://gardskart.nibio.no/search (only valid for fieldType = infield)

# Variable names

## R-friendly with janitor package
site_raw <- site_raw %>% 
  clean_names("lower_camel")
management_raw <- management_raw %>% 
  clean_names("lower_camel")

## Consistent & FAIR
names(site_raw) <- gsub("siteId", "siteID", names(site_raw))
names(site_raw) <- gsub("niBioEcologicalZone", "landscapeZone", names(site_raw))
names(site_raw) <- gsub("year", "samplingYear", names(site_raw))
names(site_raw) <- gsub("epsg25832X", "latitude25832", names(site_raw))
names(site_raw) <- gsub("epsg25832Y", "longitude25832", names(site_raw))
names(site_raw) <- gsub("habitat", "habitatType", names(site_raw))
names(management_raw) <- gsub("siteCode", "siteID", names(management_raw))
names(management_raw) <- gsub("infieldOutfield", "fieldType", names(management_raw))
names(management_raw) <- gsub("typeOfLivestock1", "livestockType", names(management_raw))
names(management_raw) <- gsub("numberOfAnimals1Adults", "numberAnimalsAdult", names(management_raw))
names(management_raw) <- gsub("numberOfAnimals1Young", "numberAnimalsYoung", names(management_raw))
names(management_raw) <- gsub("surveySiteGrazingSurfaceHa", "fieldAreaHa", names(management_raw))
names(management_raw) <- gsub("inmarksbeiteGardskart", "farmGrazingAreaHa", names(management_raw))

# Dataset

## Binding desirable variables
site_description <- left_join(
  subset(site_raw, 
         select = c(siteID, landscapeZone, samplingYear, municipality, latitude25832, longitude25832, habitatType)),
  subset(management_raw,
         select = c(siteID, fieldType, livestockType, numberAnimalsAdult, numberAnimalsYoung, fieldAreaHa, farmGrazingAreaHa))
)

## Variable types (num/chr)
# str(site_description) #validated
## Duplicate check for siteID and geolocation
# get_dupes(site_description, siteID, latitude25832, longitude25832) #validated

## Character variables - desirable categories, NAs, misprints

### Consistent lower camel
site_description <- site_description %>% 
  mutate_if(is.character, tolower)

### Categories & distribution
# table(site_description$landscapeZone) #validated
# unique(site_description$municipality) # NA
site_description[is.na(site_description$municipality),] # NA identified - US2 located in Masfjorden
site_description <- site_description |> 
  mutate(municipality = ifelse(siteID == "us2", "masfjorden", municipality))
# table(site_description$habitatType) # five categories instead of four, bog to be removed
site_description <- filter(site_description, habitatType != "bog")
site_description <- site_description |> 
  mutate(habitatType = dplyr::recode(habitatType, "forest clearing" = "transition heathland"))
# table(site_description$fieldType) #validated
# table(site_description$livestockType) # four categories instead of three
site_description <- site_description %>% 
  mutate(livestockType = dplyr::recode(livestockType, "villsau" = "sheep")) %>% 
  mutate(livestockType = dplyr::recode(livestockType, "cow" = "cattle"))

## Numeric variables - min/max, distribution, potential outliers

### Min/max
# test <- site_description |>  
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
# site_description[is.na(site_description$samplingYear),] #validated
# site_description[is.na(site_description$latitude25832),] #validated
# site_description[is.na(site_description$longitude25832),] #validated
# site_description[is.na(site_description$numberAnimalsAdult),] # Missing values from 2020 farmer interviews, replaced by values from 
site_description <- site_description |> 
  mutate(numberAnimalsAdult = ifelse(siteID == "ic2", 37, numberAnimalsAdult)) |> 
  mutate(numberAnimalsAdult = ifelse(siteID == "us1", 4, numberAnimalsAdult)) |> 
  mutate(numberAnimalsAdult = ifelse(siteID == "us4", 66, numberAnimalsAdult)) |> 
  mutate(numberAnimalsAdult = ifelse(siteID == "is5", 90, numberAnimalsAdult)) |> 
  mutate(numberAnimalsAdult = ifelse(siteID == "os7", 39, numberAnimalsAdult)) |> 
  mutate(numberAnimalsAdult = ifelse(siteID == "ig3", 109, numberAnimalsAdult)) 
  # no information for us6
# site_description[is.na(site_description$numberAnimalsYoung),] # Missing values from 2020 farmer interviews, replaced by values from 
site_description <- site_description |> 
  mutate(numberAnimalsYoung = ifelse(siteID == "ic2", 17, numberAnimalsYoung)) |> 
  mutate(numberAnimalsYoung = ifelse(siteID == "us1", 10, numberAnimalsYoung)) |> 
  mutate(numberAnimalsYoung = ifelse(siteID == "us4", 63, numberAnimalsYoung)) |> 
  mutate(numberAnimalsYoung = ifelse(siteID == "is5", 181, numberAnimalsYoung)) |> 
  mutate(numberAnimalsYoung = ifelse(siteID == "os7", 14, numberAnimalsYoung)) |> 
  mutate(numberAnimalsYoung = ifelse(siteID == "ig3", 135, numberAnimalsYoung)) 
  # no information for US6
# site_description[is.na(site_description$fieldAreaHa),] #validated - No fenced areas in the mountains
# site_description[is.na(site_description$farmGrazingAreaHa),] #validated - Applicable for infields only

### Variable distribution & outliers
# hist(site_description$numberAnimalsAdult) # Poisson distribution, one outlier over 150 animals
# site_description[site_description$numberAnimalsAdult>150,] # is4, no observed inconsistency with farm characteristics
# hist(site_description$numberAnimalsYoung) # Poisson distribution, no outlier
hist(site_description$fieldAreaHa) # One outlier over 1000 ha crushing the distribution
site_description[site_description$fieldAreaHa>1000,] # ug2 outfield site in upland area, value non applicable for analysis
# hist(filter(site_description, siteID != "ug2")$fieldAreaHa) # Poisson distribution
# hist(site_description$farmGrazingAreaHa) # dominance small farms, no visible outliers

# Dataset export
write_csv(site_description, "data/cleandata/tradmod_nbr_sitedescription.csv")

#### LANDSCAPE MATRIX ####

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

#### 20x20 SAMPLING AREA ####

# Raw datasets
## site_raw
## sampling_area_raw

# Desired variables
## siteID - Field identification code for data collection
## samplingDate - Date of data collection
## numberLivestockPath - Number of paths created by the livestock in the sampling area
## lengthLivestockPathM - Total length of livestock path in the sampling area in meters
## elevationMasl - Elevation at the highest point of the sampling area in meters above sea level
## slopeAngleDegree - Slope angle from the highest point to the lowest point of the sampling area in degree
## slopeAspectDegree - Slope aspect from the highest point to the lowest point of the sampling area in degree
## percentRock - Estimated percentage cover of exposed rock in the sampling area
## percentMud - Estimated percentage cover of visible mud in the sampling area
## percentTreeTallShrub - Estimated percentage cover of trees and shrubs over 1 m in the sampling area
## percentLowShrub - Estimated percentage cover of shrubs under 1 m in the sampling area
## percentForb - Estimated percentage cover of forbs in the sampling area
## percentMonocotyledon - Estimated percentage cover of monocotyledons (grasses, rushes, sedges) in the sampling area
## percentBryophyte - Estimated percentage cover of bryophytes (mosses, liverworts) in the sampling area
## percentFern - Estimated percentage cover of ferns in the sampling area
## percentLichen - Estimated percentage cover of lichens in the sampling area

#
## Summary sampling area - Check table size, list of variables, variable types (num/chr)

#str(area20x20_raw) # All good, date should be reformatted

#
## Name & character cleaning sampling area

# R friendly variable names
names(area20x20_raw) <- gsub("%", "percent", names(area20x20_raw)) # remove percent signs from names
names(area20x20_raw) <- gsub("&", "_", names(area20x20_raw)) # remove &
names(area20x20_raw) <- gsub("\\(", "", names(area20x20_raw)) # remove (
names(area20x20_raw) <- gsub("\\)", "", names(area20x20_raw)) # remove )
names(area20x20_raw) <- gsub("Recording_date", "RecordingDate", names(area20x20_raw))
names(area20x20_raw) <- gsub("Number_paths", "NumberAnimalPaths", names(area20x20_raw))
names(area20x20_raw) <- gsub("Total_length_path_m", "AnimalPathsLength_m", names(area20x20_raw))
names(area20x20_raw) <- gsub("Elevation_max", "Elevation", names(area20x20_raw))
names(area20x20_raw) <- gsub("General_slope", "Slope_degree", names(area20x20_raw))
names(area20x20_raw) <- gsub("AspectDegree", "Aspect_degree", names(area20x20_raw))
names(area20x20_raw) <- gsub("percent", "Percent", names(area20x20_raw))
names(area20x20_raw) <- gsub("Trees_Tall", "TreesTall", names(area20x20_raw))
names(area20x20_raw) <- gsub("Fern", "Ferns", names(area20x20_raw))
names(area20x20_raw) <- gsub("Lichen", "Lichens", names(area20x20_raw))
names(area20x20_raw) <- gsub("Comments", "Comments_area20x20", names(area20x20_raw))

#
## Sampling date standardisation

area20x20_raw$RecordingDate <- as.POSIXct(area20x20_raw$RecordingDate, format = "%d.%m.%Y")

#
## Data cleaning - New R object

area20x20_full <- subset(area20x20_raw, select = -c(Team, Latitude1, Longitude1, Latitude2, Longitude2, NumberAnimalPaths, AnimalPathsLength_m, Elevation_min, Aspect))

#
## Char var Site Info - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(area20x20_full$SiteID) # Unique ID for each site - validated

#
## Numeric var sampling area - Check min/max, distribution and potential outliers

# Check min/max
test <- area20x20_full |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # potential outliers are max herbs 97% and max lichens 20%

# Distribution elevation max
#hist(area20x20_full$Elevation) # Sites range from 0 to 900 m elevation, no outlier - validated

# Distribution slope
#hist(area20x20_full$Slope_degree) # Slopes range from 0 to 35 degrees, no outlier - validated

# Distribution aspect degree
#hist(area20x20_full$Aspect_degree) # Aspect covers all spectrum (0 to 360), no outliers - validated

# Distribution distance to sea
#area20x20_full[is.na(area20x20_full$DistanceToSea_m),] # 9 NA, corresponding to upland sites -> validated
#hist(area20x20_full$DistanceToSea_m) # One outlier, above 10 km distance
#filter(area20x20_full, DistanceToSea_m>10000) # IG3 in Modalen

# Distribution rock cover
#hist(area20x20_full$PercentRock) # Some sites over 7%, check their location + 3 NAs
#area20x20_full[area20x20_full$PercentRock>7,] # 5 sites OV1, UG2, US3, IG3, US5
# UG2, US3 and US5 in subalpine areas with exposed bedrock - validated
# OV1 in coastal heathland habitat, with exposed bedrock - validated
# IG3 in fjord area but at higher elevation, with exposed bedrock - validated
test <- area20x20_full[is.na(area20x20_full$PercentRock),] # 3 sites US1, UG1 and OC4 missing all soil cover percentages

# Distribution mud cover
#hist(area20x20_full$PercentMud) # One site over 15%
#area20x20_full[area20x20_full$PercentMud>10,] # UC1 -> bog site, will not be included into the analysis - validated

# Distribution trees and tall shrubs cover
#hist(area20x20_full$PercentTreesTallShrubs) # No sites over 10% - validated

# Distribution low shrubs
#hist(area20x20_full$PercentLowShrubs) # Wide range due to collection in both grassland and heathland habitats. Check that all grassland sites are under 10%
#area20x20_full[area20x20_full$PercentLowShrubs>10,] # 11 sites over 10%
# OV1, OV2, OS5, OS7, OS9, IS2 coastal heathlands - validated
# US2, US3, US4, US5, UG2 subalpine heathlands - validated

# Distribution herbs
hist(area20x20_full$PercentHerbs) # a few sites over 60%
area20x20_full[area20x20_full$PercentHerbs>60,] # OS1, OC1, IG1, IS2, IC1 -> all first year/starting sites, check on vegetation quadrats + site & plot pictures
# OS1 80% - average 20% & no cover over 55% in quadrats, estimation from pictures 35%-40%
# OC1 80% - average 25% & no cover over 50% in quadrats, estimation from pictures 10%-15% 
# IG1 97% - average 20% & no cover over 30% in quadrats, estimation from pictures 15%-20%
# IS2 80% - average 70% & no cover over 90% in quadrats, estimation from pictures 45%-50%
# IC1 70% - average 45% & no cover over 70% in quadrats, estimation from pictures 60%-65%

# Distribution monocotyledons
hist(area20x20_full$PercentMonocotyledons) # need to check again sites with weird herb estimations

# Distribution bryophytes
hist(area20x20_full$PercentBryophytes) # need to check again sites with weird herb estimations

# Distribution lichens
#hist(area20x20_full$PercentLichens) # one site over 10%
#area20x20_full[area20x20_full$PercentLichens>15,] # US4 in subalpine area, average of 12% & max 24% in quadrats - validated

#
## Calculation new variables

# Heat Load Index
area20x20_full <- area20x20_full |> 
  mutate(HeatLoadIndex = cos(Aspect_degree-225)*tan(Slope_degree))
#hist(area20x20_full$HLI) # 3 outliers: one under 200, two over 100
#area20x20_full[area20x20_full$HLI>100,] #OG4 & IS3 -> both 11 degree slope with SW & SE exposition
#area20x20_full[area20x20_full$HLI<0,] #OS6 -> 11 degree slope with NE exposition


## Export clean data in new excel file

write_csv(area20x20_full, "data/cleandata/NBR_FullArea20x20.csv")


#### Ground cover quadrats ####

## Description

## List of variables

# [1] Field identification code for data collection
# [2] Date of data collection
# [3] Sample identification code
# [4] Percent cover of bare soil in the quadrat
# [5] Percent cover of exposed rock in the quadrat
# [6] Percent cover of litter in the quadrat
# [7] Percent cover of dead wood in the quadrat
# [8] Percent cover of bryophytes in the quadrat
# [9] Percent cover of lichens in the quadrat
# [10] Percent cover of vascular plants in the quadrat
# [11] Percent cover of blossom in the quadrat
# [12] Species in blossom (1)
# [13] Species in blossom (2)
# [14] Species in blossom (3)
# [15] Species in blossom (4)
# [16] Species in blossom (5)
# [17] Percent cover of dung in the quadrat
# [18] Vegetation mean height in the quadrat, in cm
# [19] Vegetation max height in the quadrat, in cm
# [20] Plant species richness in the quadrat
# [21] Comments

#
## Summary sampling area - Check table size, list of variables, variable types (num/chr)

#str(groundcover_raw) # Date should be reformatted, plotID renamed as sampleID and plotID created

#
## Name & character cleaning sampling area

# R friendly variable names
names(groundcover_raw) <-gsub ("Site", "SiteID", names(groundcover_raw)) # rename in SiteID so it matches with other files
names(groundcover_raw)<- gsub ("PlotID", "SampleID", names(groundcover_raw)) # Rename plotID as sampleID
names(groundcover_raw)<- gsub ("Date", "Recording_date", names(groundcover_raw)) # Rename so it matches with other files
names(groundcover_raw)<- gsub ("\\(", "", names(groundcover_raw)) # remove (
names(groundcover_raw)<- gsub ("\\)", "", names(groundcover_raw)) # remove )
names(groundcover_raw) <- gsub("Comments", "Comments_soilcover", names(groundcover_raw))
groundcover_raw$PlotID <- substr(groundcover_raw$SampleID, 1,6) #create PlotID column

#
## Sampling date standardisation

groundcover_raw$Recording_date <- as.POSIXct(groundcover_raw$Recording_date, format = "%d.%m.%Y")

#
## Data cleaning - New R object

groundcover_full <- groundcover_raw

#
## Char var - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(groundcover_full$SiteID) # 15 samples per site - validated

# Plot ID
#table(groundcover_full$PlotID) # 5 samples per plot - validated

# Sample ID
#groundcover_full[duplicated(groundcover_full$SampleID),] # Unique ID for sample - validated

# Blossom sp1
#table(groundcover_full$Blossom_sp1) # two latin names for Cirsium palustre + "Alchemilla millefolium" either Achillea millefolium or Alchemilla vulgaris
groundcover_full <- groundcover_full |> 
  mutate(Blossom_sp1 = dplyr::recode(Blossom_sp1, "Cirsium palustris" = "Cirsium palustre"))
#filter(groundcover_full, Blossom_sp1 == "Alchemilla millefolium") # IS1-P3-N3 -> check on field sheet -> species is Achillea millefolium
groundcover_full <- groundcover_full |> 
  mutate(Blossom_sp1 = dplyr::recode(Blossom_sp1, "Alchemilla millefolium" = "Achillea millefolium"))

# Blossom sp2
#table(groundcover_full$Blossom_sp2) # Bad ID Leontodon saxatile, should be Leontodon autumnalis
groundcover_full <- groundcover_full |> 
  mutate(Blossom_sp2 = dplyr::recode(Blossom_sp2, "Leontodon saxatile" = "Leontodon autumnalis"))

# Blossom sp3
#table(groundcover_full$Blossom_sp3) # All good

# Blossom sp4
#table(groundcover_full$Blossom_sp4) # All good

# Blossom sp5
#table(groundcover_full$Blossom_sp5) # All good

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- groundcover_full |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # no visible outlier (maybe lichen at 80%?), no percentage above 100%

# Bare soil - NA + Distribution
#groundcover_full[is.na(groundcover_full$Bare_soil),] # No NA
#hist(groundcover_full$Bare_soil) # Samples range from 0 to 40% -> check the maximum
#filter(groundcover_full, Bare_soil>30) # Two samples on OV1, recently burnt coastal heathland -> validated

# Rocks - NA + Distribution
#groundcover_full[is.na(groundcover_full$Rocks),] # No NA
#hist(groundcover_full$Rocks) # Samples range from 0 to 8% -> validated

# Litter - NA + Distribution
#groundcover_full[is.na(groundcover_full$Litter),] # No NA
#hist(groundcover_full$Litter) # Samples range from 0 to 60% -> check above 30%
#filter(groundcover_full, Litter>30) # Samples either from OV1, recently burnt coastal heathland, or mountain sites (US2, UG1) -> validated

# Dead wood - NA + Distribution
#groundcover_full[is.na(groundcover_full$Dead_wood),] # No NA
#hist(groundcover_full$Dead_wood) # Samples range from 0 to 10% -> validated

# Bryophytes - NA + Distribution
#groundcover_full[is.na(groundcover_full$Bryophytes),] # No NA
#hist(groundcover_full$Bryophytes) # Samples range from 0 to 100% -> validated

# Lichen - NA + Distribution
#groundcover_full[is.na(groundcover_full$Lichen),] # No NA
#hist(groundcover_full$Lichen) # Samples range from 0 to 80%% -> check outliers above 40%
#filter(groundcover_full, Lichen>40) # 1 sample from mountain site (US1-P1-N5) -> validated from the field sheet

# Vascular - NA + Distribution
#groundcover_full[is.na(groundcover_full$Vascular),] # No NA
#hist(groundcover_full$Vascular) # Samples range from 0 to 100% -> validated

# Blossom cover - NA + Distribution
#groundcover_full[is.na(groundcover_full$Blossom_cover),] # No NA
#hist(groundcover_full$Blossom_cover) # Samples range from 0 to 15% -> validated

# Dung cover - NA + Distribution
#groundcover_full[is.na(groundcover_full$Dung),] # No NA
#hist(groundcover_full$Dung) # Samples range from 0 to 15% -> check maximum
#filter(groundcover_full, Dung>10) # Samples from cow site (IC1) -> validated

# Vegetation maximum height - NA + Distribution
#groundcover_full[is.na(groundcover_full$VG_max_height_cm),] # No NA
#hist(groundcover_full$VG_max_height_cm) # Samples range from 0 to 160 cm with normal distribution -> validated

# Vegetation mean height - NA + Distribution
#groundcover_full[is.na(groundcover_full$VG_mean_height_cm),] # No NA
#hist(groundcover_full$VG_mean_height_cm) # Samples range from 0 to 90 cm with Poisson distribution -> validated

# Vegetation species richness - NA + Distribution
#groundcover_full[is.na(groundcover_full$Plant_species_richness),] # No NA
#hist(groundcover_full$Plant_species_richness) # Samples range from 0 to 35 species with Normal distribution -> validated


## Export clean data in new excel file

write_csv(groundcover_full, "data/cleandata/NBR_FullGroundCover.csv")



#### Soil penetration quadrat ####

## Description

## List of variables

# [1] Field identification code for data collection
# [2] Date of data collection
# [3] Sample identification code
# [4] Stick height remaining above the ground on the left side of the quadrat (cm)
# [5] Stick height which penetrated the ground on the left side of the quadrat (cm) -> calculated in Excel
# [6] If the stick hit a rock on the left side of the quadrat Y/N
# [7] Stick height remaining above the ground on the right side of the quadrat (cm)
# [8] Stick height which penetrated the ground on the right side of the quadrat (cm) -> calculated in Excel
# [9] If the stick hit a rock on the left side of the quadrat Y/N
# [10] Average soil penetration for quadrat (mean left and right)
# [11] Stick total length (cm) - the stick would wear out along with repetitive use, so its length would decrease over time
# [12] Comments

#
## Summary - Check table size, list of variables, variable types (num/chr)

#str(soilpene_raw) # Date should be reformatted, plotID renamed as sampleID and plotID created

#
## Name & character cleaning

# R friendly variable names
names(soilpene_raw) <- gsub("Site", "SiteID", names(soilpene_raw)) # rename in SiteID so it matches with other files
names(soilpene_raw) <- gsub("PlotID", "SampleID", names(soilpene_raw)) # Rename plotID as sampleID
names(soilpene_raw) <- gsub("Date", "Recording_date", names(soilpene_raw)) # Rename so it matches with other files
names(soilpene_raw) <- gsub("\\(", "", names(soilpene_raw)) # remove (
names(soilpene_raw) <- gsub("\\)", "", names(soilpene_raw)) # remove )
names(soilpene_raw) <- gsub("Comments", "Comments_soilpene", names(soilpene_raw))
soilpene_raw$PlotID <- substr(soilpene_raw$SampleID, 1,6) #create PlotID column

#
## Sampling date standardisation

soilpene_raw$Recording_date <- as.POSIXct(soilpene_raw$Recording_date, format = "%d.%m.%Y")

#
## Data cleaning - New R object

soilpene_full <- soilpene_raw

#
## Char var - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(soilpene_full$SiteID) # 12 samples per site - validated

# Plot ID
#table(soilpene_full$PlotID) # 4 samples per plot - validated

# Sample ID
#soilpene_full[duplicated(soilpene_full$SampleID),] # Unique ID for sample - validated

# If rock hit on the left or right
#unique(soilpene_full$Left_rock_hit) # only Y & N -> validated
#unique(soilpene_full$Right_rock_hit) # only Y & N -> validated

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- soilpene_full |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # no visible outlier, no penetration or standing part above stick full length

# Left standing height
#soilpene_full[is.na(soilpene_full$Left_standing_part_cm),] # No NA
#hist(soilpene_full$Left_standing_part_cm) # Heights range from 28 to 42 cm in a normal distribution -> validated

# Left penetration height
#soilpene_full[is.na(soilpene_full$Left_PT_cm),] # No NA
#hist(soilpene_full$Left_PT_cm) # Heights range from 0 to 14 cm in a normal distribution -> validated

# Right standing height
#soilpene_full[is.na(soilpene_full$Right_standing_part_cm),] # No NA
#hist(soilpene_full$Right_standing_part_cm) # Heights range from 25 to 45 cm in a normal distribution -> validated

# Right penetration height
#soilpene_full[is.na(soilpene_full$Right_PT_cm),] # No NA
#hist(soilpene_full$Right_PT_cm) # Heights range from 0 to 16 cm in a normal distribution -> validated

# Total stick length
#soilpene_full[is.na(soilpene_full$Stick_height),] # No NA
#hist(soilpene_full$Stick_height) # Heights range from 42.7 to 43.4 cm -> validated


## Export clean data in new excel file

write_csv(soilpene_full, "data/cleandata/NBR_FullSoilPene.csv")



#### Soil bulk density - quadrats ####

## Description

## List of variables

# [1] Sample identification code
# [2] Field identification code for data collection
# [3] Plot identification code
# [4] Height of the core fully completed with soil (cm)
# [5] Volume of the soil core before correction for holes or slopes (cm3) -> not to be used in the analysis
# [6] Volume of the soil core manually corrected for holes and slopes if applicable (cm3) -> not to be used in the analysis
# [7] Best estimation of the volume of the soil core, with correction for holes or slopes if needed (cm3)
# [8] Core weight (including soil + PVC core + cheesecloth) on fresh soil, before water saturation (g)
# [9] Core weight (including soil + PVC core + cheesecloth) after water saturation (g)
# [10] Core weight (including soil + PVC core + cheesecloth) after 24h of drying (g)
# [11] Core weight (including soil + PVC core + cheesecloth) after 48h of drying (g)
# [12] Core weight (including soil + PVC core + cheesecloth) after drying at 105C in oven (g)
# [13] Weight of the cheesecloth (g)
# [14] Weight of the PVC core with the cheesecloth (g)
# [15] Percentage of water loss over 24h -> calculated from W0 and W24
# [16] Percentage of water loss over 48h -> calculated from W0 and W48
# [17] Bulk density calculated from the best volume estimation [7]
# [18] Weight of percentage of soil moisture (g)
# [19] Volume of percentage of soil moisture (cm3) -> calculated from the BD
# [20] Percentage of soil porosity -> calculated from the BD
# [21] Percentage of WFPS
# [22] Comments during the lab processing of the soil
# [23] Other comments
# [24] If the samples are concerned by scale calibration issue
# [25] Who performed the task

#
## Summary - Check table size, list of variables, variable types (num/chr)

#str(soilbulk_raw) # missing sample ID, plotID to be reformated

#
## Name & character cleaning

# R friendly variable names
names(soilbulk_raw) <- gsub("\\(", "", names(soilbulk_raw)) # remove (
names(soilbulk_raw) <- gsub("\\)", "", names(soilbulk_raw)) # remove )
names(soilbulk_raw) <- gsub(" ", "", names(soilbulk_raw)) # remove spaces
names(soilbulk_raw) <- gsub("Site", "SiteID", names(soilbulk_raw)) # rename in SiteID so it matches with other files
names(soilbulk_raw) <- gsub("cm", "_cm", names(soilbulk_raw))
names(soilbulk_raw) <- gsub("%", "percent_", names(soilbulk_raw))
names(soilbulk_raw) <- gsub("Comments_processing", "CommentsProcessing_soilbulk", names(soilbulk_raw))
names(soilbulk_raw) <- gsub("Other_comments", "OtherComments_soilbulk", names(soilbulk_raw))

# New ID variables
soilbulk_raw$PlotID <- substr(soilbulk_raw$BDcoreID, 1,6) # recreate PlotID column
soilbulk_raw$SampleID <- substr(soilbulk_raw$BDcoreID, 1,9) # recreate SampleID column

#
## Data cleaning - New R object

soilbulk_full <- soilbulk_raw

#
## Char var - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(soilbulk_full$SiteID) # 36 samples per site - validated

# Plot ID
#table(soilbulk_full$PlotID) # 12 samples per plot - validated

# Bulk density core ID
#soilbulk_full[duplicated(soilbulk_full$BDcoreID),] # Unique ID for core - validated

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- soilbulk_full |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # min core volume quite low, negative values for water loss 24h and 48h, negative BD

# Best estimation soil core volume
#soilbulk_full[is.na(soilbulk_full$CoreVol),] # two NA in IS1 & OG1 + 24 NAs in IC3 -> check datasheet -> lab incident, the 2 cores were discarded - Samples missing for IC3, never found
#hist(soilbulk_full$CoreVol) # Volumes range from 15 to 60 cm3, most above 45-50 cm3 -> low volume = bad estimation of BD, too low volumes should be discarded

# W0 - Weight of fresh soil before saturation
#soilbulk_full[is.na(soilbulk_full$W0g),] # same NAs -> validated
#hist(soilbulk_full$W0g) # Weights range from 50 to 180 g in a normal distribution -> very low weights likely to be linked to low volumes

# WSAT - Soil weight after water saturation
#soilbulk_full[is.na(soilbulk_full$WSAT),] # same NAs -> validated
#hist(soilbulk_full$WSAT) # Weights range from 60 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes

# W24H - Soil weight after 24h of drying
#soilbulk_full[is.na(soilbulk_full$W24H),] # same NAs -> validated
#hist(soilbulk_full$W24H) # Weights range from 50 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes, not so much difference compared to WSAT

# W48H - Soil weight after 48h of drying
#soilbulk_full[is.na(soilbulk_full$W48H),] # same NAs -> validated
#hist(soilbulk_full$W48H) # Weights range from 50 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes

# WDRY - Soil weight after over at 105C
#soilbulk_full[is.na(soilbulk_full$WDRY),] # same NAs -> validated
#hist(soilbulk_full$WDRY) # Weights range from 20 to 140 g in a normal distribution -> very low weights likely to be linked to low volumes

# Percent water loss in 24h
#soilbulk_full[is.na(soilbulk_full$percent_Waterloss24h),] # same NAs -> validated
#hist(soilbulk_full$percent_Waterloss24h) # % range from -40% to 40%, main between 0 and 10% -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Percent water loss in 48h
#soilbulk_full[is.na(soilbulk_full$percent_Waterloss48h),] # same NAs -> validated
#hist(soilbulk_full$percent_Waterloss48h) # % range from -70% to 70%, main between 0 and 20% -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Bulk density
#soilbulk_full[is.na(soilbulk_full$BD),] # same NAs -> validated
#hist(soilbulk_full$BD) # % range from -0.2% to 3, normal distribution -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Soil moisture in percentage weight = gravimetric water content
#soilbulk_full[is.na(soilbulk_full$Weightpercent_Soilmoisture),] # same NAs -> validated
#hist(soilbulk_full$Weightpercent_Soilmoisture) # % range from -50% to 100%, normal distribution -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Soil moisture in percentage volume
#soilbulk_full[is.na(soilbulk_full$Volpercent_Soilmoisture),] # same NAs -> validated
#hist(soilbulk_full$Volpercent_Soilmoisture) # % range from -50% to 100%, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

# Soil porosity
#soilbulk_full[is.na(soilbulk_full$percent_Soilporosity),] # same NAs -> validated
#hist(soilbulk_full$percent_Soilporosity) # % range from -70% to 80%, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

# WFPS
#soilbulk_full[is.na(soilbulk_full$percent_WFPS),] # same two NAs -> validated
#hist(soilbulk_full$percent_WFPS) # % range from -50% to one outlier over 10000, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

#
## New variables - gravimetric & volumetric water content from standardised W+48h dried soil

# New variables
soilbulk_full <- soilbulk_full |> 
  mutate(GWC_48 = (W48H - WDRY)/W48H*100) |> 
  mutate(VWC_48 = GWC_48*BD)

# Distribution
hist(soilbulk_full$GWC_48) # Normal distribution, from 20% to 80% -> some very high values
hist(soilbulk_full$VWC_48) # Normal distribution, from 5% to 60%

#
## Data filtering

# Min soil core volume
#filter(soilbulk_full, CoreVol<40 & !is.na(CoreVol)) # 40 or 2% samples unfit
#filter(soilbulk_full, CoreVol<45 & !is.na(CoreVol)) # 111 or 7% samples unfit
#filter(soilbulk_full, CoreVol<50 & !is.na(CoreVol)) # 344 or 21% samples unfit -> cores should be minimum vol of 50 cm3

# Negative water loss values
#filter(soilbulk_full, percent_Waterloss24h<0 & !is.na(percent_Waterloss24h)) #132 samples with water loss 24h negative
#filter(soilbulk_full, percent_Waterloss48h<0 & !is.na(percent_Waterloss48h)) #82 samples with water loss 48h negative

# Selection data with min 50 cm3 soil volume and positive water loss
soilbulk_full <- subset(soilbulk_full, CoreVol>50)
soilbulk_full <- subset(soilbulk_full, percent_Waterloss24h>0)
soilbulk_full <- subset(soilbulk_full, percent_Waterloss48h>0)

# Check new variable distribution - water loss 24h
#hist(soilbulk_full$percent_Waterloss24h) # still some extreme values over 20%
#filter(soilbulk_full, percent_Waterloss24h>20) # 6 cores with more than 20% over 24h
# 2 cores from UC1, which is excluded from the analysis -> should be removed
# 1 cores from OC3, concerned with scale issue (lots of negative values which are already removed). Water loss between 0-24 and 24-48 not coherent -> should be removed
# 2 cores from OC2, concerned with scale issue. Water loss between 0-24 and 24-48 not coherent with other samples from same plot (W48h>W24h for P1-D1_2) -> should be removed
# 1 cores from OC5, concerned with scale issue. Water loss between 0-24 and 24-48 not coherent with other samples from same site -> should be removed
soilbulk_full <- subset(soilbulk_full, percent_Waterloss24h<20)

# Check new variable distribution - water loss 48h
#hist(soilbulk_full$percent_Waterloss48h) # still some extreme values over 25%
#filter(soilbulk_full, percent_Waterloss48h>25) # 2 cores with more than 25% over 48h
# OG6-P1-D3_1, not concerned by the scale issue and with values from other cores coherent -> to be kept
# OC2-P2-D1_3, concerned with scale issue - value not coherent with water loss 24h and with other cores -> to be removed
soilbulk_full <- subset(soilbulk_full, BDcoreID != "OC2-P2-D1_3")

# Check new variable distribution - bulk density
#hist(soilbulk_full$BD) # no negative values anymore, quite nice normal distribution -> validated

# Check new variable distribution - soil moisture in percent weight
#hist(soilbulk_full$Weightpercent_Soilmoisture) # still some negative and extreme values (100%)
#filter(soilbulk_full, Weightpercent_Soilmoisture<20) # 5 cores with less than 20% soil moisture
# 3 cores from OC2, concerned with scale issue. OC2-P2-D2_2 negative value, OC2-P1-D1_3 very low not coherent with other samples from the plot -> to be removed - OC2-P3-D3_3 just under 20, not extreme compared with the other samples -> to be kept
# OG4-P3-D3_1, concerned with scale issue -> values are coherent within the plot and relatively close to what is find in other plots (10%-30%) -> to be kept
soilbulk_full <- subset(soilbulk_full, BDcoreID != "OC2-P2-D2_2")
soilbulk_full <- subset(soilbulk_full, BDcoreID != "OC2-P1-D1_3")
#filter(soilbulk_full, Weightpercent_Soilmoisture>90) # IS3-P3-D4_1, with P3 concerned with scale issue. Incoherent with other samples from same plot -> to be removed
soilbulk_full <- subset(soilbulk_full, Weightpercent_Soilmoisture<90)

# Check new variable distribution - soil moisture in percent volume
#hist(soilbulk_full$Volpercent_Soilmoisture) # no extreme. nice normal distribution

# Check new variable distribution - standardised gravimetric water content
#hist(soilbulk_full$GWC_48) # no extreme. nice normal distribution

# Check new variable distribution - standardised volumetric water content
#hist(soilbulk_full$VWC_48) # no extreme. nice normal distribution

# Check new variable distribution - WFPS
#hist(soilbulk_full$percent_WFPS) # no extremes, nice normal distribution

# Check new number of replicates per site
#sort(table(soilbulk_full$SiteID)) 
# 9 sites with less than 20 replicates and lowest IC3 with 9 replicates (due to missing values) -> validated


## Export clean data in new excel file

write_csv(soilbulk_full, "data/cleandata/NBR_FullSoilBulk.csv")



#### Soil chemistry ####

## Description
# 3 datasets: one from 2019, one from 2020, and one complementary for a missing variable in 2020

## List of variables

# [1] Eurofins protocole -> not to be used in the analysis
# [2] Eurofins protocole -> not to be used in the analysis
# [3] Eurofins protocole -> not to be used in the analysis
# [4] Eurofins protocole -> not to be used in the analysis
# [5] Eurofins protocole -> not to be used in the analysis
# [6] Eurofins protocole -> not to be used in the analysis
# [7] Eurofins protocole -> not to be used in the analysis
# [8] Eurofins protocole -> not to be used in the analysis
# [9] Eurofins protocole -> not to be used in the analysis
# [10] Eurofins protocole -> not to be used in the analysis
# [11] Eurofins protocole -> not to be used in the analysis
# [12] Eurofins protocole -> not to be used in the analysis
# [13] Plot identification code
# [14] Eurofins protocole -> not to be used in the analysis
# [15] Eurofins protocole -> not to be used in the analysis
# [16] Soil type code (e.g. clay, sand, humus)
# [17] Clay class code
# [18] Loss of Ignition (LOI)
# [19] Soil density (kg/L)
# [20] Humus quantity in percent dry matter
# [21] Humus class code
# [22] pH
# [23] Rate of phosphorus (mg/100g)
# [24] Rate of potassium (mg/100g)
# [25] Rate of magnesium (mg/100g)
# [26] Rate of calcium (mg/100g)
# [27] Rate of potassium nitrate (mg/100g)
# [28] Rate of copper (mg/100g)
# [29] Rate of boron (mg/100g)
# [30] Rate of sodium (mg/100g)
# [31] Rate of sulfur (mg/100g)
# [32] Rate of iron (mg/100g)
# [33] Rate of manganese (mg/100g)
# [34] Rate of zinc (mg/100g)
# [35] Rate of molybdenum (mg/100g)
# [36] Rate of selenium (mg/100g)

#
## Summary - Check table size, list of variables, variable types (num/chr)

#str(chem2019) # two missing variables (dry matter and total N) available in the PDF version of the document
#str(chem2020) # two missing variables (dry matter and total N) respectively available in the complementary dataset and in the PDF version of the document
#str(chem2020_DM) # dry matter as character

#
## Common ID for merging

names(chem2019) <- gsub("Provenummer", "PlotID", names(chem2019)) # common ID
names(chem2020) <- gsub("Provenummer", "PlotID", names(chem2020)) # common ID
names(chem2020_DM) <- gsub("Merking", "PlotID", names(chem2020_DM)) # common ID

#
## Filling missing variables

# Soil chemistry 2019 - dry matter and total N
extra2019 <- data.frame(
  PlotID = c("OC11", "OC12", "OC13", "OS11", "OS12", "OS13", "OG11", "OG12", "OG13", "OS21", "OS22", "OS23", "OV21", "OV22", "OV23", "IC11", "IC12", "IC13", "IG11", "IG12", "IG13", "IG21", "IG22", "IG23", "IS11", "IS12", "IS13", "IS21", "IS22", "IS23", "IC21", "IC22", "IC23", "IV11", "IV12", "IV13", "UG11", "UG12", "UG13", "UG21", "UG22", "UG23", "US11", "US12", "US13", "US21", "US22", "US23", "US31", "US32", "US33", "US41", "US42", "US43"),
  DryMatter_percent = c(68.9, 66.7, 70.8, 96.8, 97.4, 97, 94.2, 94.9, 95.7, 95.6, 95.6, 95.5, 88.7, 93, 93.5, 96.2, 95.4, 95.8, 94.8, 95.4, 95.1, 93.6, 92.9, 94.4, 94.8, 94.4, 94.9, 94.7, 93.5, 93.7, 95.1, 95, 93.3, 96.7, 96.4, 97.1, 83.8, 89.1, 90.6, 94.8, 94.3, 90.6, 95.2, 94.5, 89, 94.6, 90.7, 94.3, 86.7, 93.5, 91.3, 92.9, 94.4, 93),
  TotalN_percentDM = c(1.98, 2.08, 1.46, 0.28, 0.29, 0.34, 0.88, 0.65, 0.59, 0.54, 0.63, 0.62, 1.14, 0.61, 0.72, 0.5, 0.62, 0.5, 0.4, 0.45, 0.51, 0.76, 0.83, 0.74, 0.64, 0.62, 0.67, 0.59, 0.64, 0.92, 0.64, 0.7, 0.92, 0.52, 0.55, 0.71, 1.78, 1.53, 1.41, 0.71, 0.78, 1.26, 0.66, 0.63, 1.57, 0.72, 1.2, 0.94, 1.6, 0.62, 0.86, 1.08, 0.84, 0.67)
)
chem2019 <- full_join(chem2019, extra2019)

# Soil chemistry 2020 - total N
extra2020 <- data.frame(
  PlotID = c("OC21", "OC22", "OC23", "OC31", "OC32", "OC33", "OC41", "OC42", "OC43", "OC51", "OC52", "OC53", "OG21", "OG22", "OG23", "OG31", "OG32", "OG33", "OG51", "OG52", "OG53", "OG61", "OG62", "OG63", "IG31", "IG32", "IG33", "IS31", "IS32", "IS33", "IS41", "IS42", "IS43", "IS51", "IS52", "IS53", "OS31", "OS32", "OS33", "OS41", "OS42", "OS43", "OS51", "OS52", "OS53", "OS61", "OS62", "OS63", "OS71", "OS72", "OS73", "OS81", "OS82", "OS83", "OS91", "OS92", "OS93", "IC31", "IC32", "IC33", "IC41", "IC42", "IC43", "IC51", "IC52", "IC53", "OG41", "OG42", "OG43", "US51", "US52", "US53", "US61", "US62", "US63"),
  TotalN_percentDM = c(0.25, 0.32, 0.33, 0.88, 1.13, 0.84, 0.54, 0.63, 0.63, 0.63, 0.76, 1.26, 0.69, 0.34, 0.47, 1.97, 1.4, 2.04, 1.57, 0.87, 0.42, 0.53, 0.5, 0.4, 0.27, 0.73, 0.5, 0.68, 0.95, 0.69, 0.41, 0.93, 0.42, 1.04, 1.21, 1.66, 0.69, 0.53, 0.5, 0.42, 0.3, 0.35, 0.61, 0.48, 0.42, 0.39, 0.4, 0.41, 1.21, 1.5, 1.82, 1.77, 1.96, 1.51, 0.76, 0.84, 0.63, 0.61, 0.77, 0.68, 0.31, 0.45, 0.39, 0.55, 1.33, 1.27, 0.56, 0.45, 0.27, 0.79, 0.85, 0.97, 0.65, 0.56, 1.22)
)
chem2020 <- full_join(chem2020, extra2020)

# Soil chemistry 2020 - dry matter
names(chem2020_DM) <- gsub("Torrstoff.g.100g", "DryMatter_percent", names(chem2020_DM))
chem2020_DM$DryMatter_percent <- as.numeric(chem2020_DM$DryMatter_percent)
chem2020 <- left_join(chem2020, chem2020_DM)

# Merging 2019 and 2020 tables
soilchem_raw <- full_join(chem2019, chem2020)

#
## Name & character cleaning

# R friendly variable names
names(soilchem_raw) <- gsub("Jordart", "SoilType", names(soilchem_raw))
names(soilchem_raw) <- gsub("Leirklasse", "ClayCategory", names(soilchem_raw))
names(soilchem_raw) <- gsub("Glodetap", "LOI", names(soilchem_raw))
names(soilchem_raw) <- gsub("Volumvekt", "SoilDensity_kg.L", names(soilchem_raw))
names(soilchem_raw) <- gsub("Mold", "Humus_percentDM", names(soilchem_raw))
names(soilchem_raw) <- gsub("Humus_percentDMklasse", "HumusCategory", names(soilchem_raw))
names(soilchem_raw) <- gsub("P.Al", "P.Al_mg.100g", names(soilchem_raw))
names(soilchem_raw) <- gsub("K.Al", "K.Al_mg.100g", names(soilchem_raw))
names(soilchem_raw) <- gsub("Mg.Al", "Mg.Al_mg.100g", names(soilchem_raw))
names(soilchem_raw) <- gsub("Ca.Al", "Ca.Al_mg.100g", names(soilchem_raw))
names(soilchem_raw) <- gsub("KHNO3", "KHNO3_mg.100g", names(soilchem_raw))
names(soilchem_raw) <- gsub("Na.Al", "Na.Al_mg.100g", names(soilchem_raw))

# Removal dummy and/or empty variables from Eurofins protocol
soilchem_raw <- soilchem_raw |>   
  discard(~all(is.na(.) | . =="")) # remove all empty columns
soilchem_raw <- subset(soilchem_raw, select = -c(Arstall, Journalnr, Navn, Postnr, Poststed, Registreringsdato, batchCode, contactName, Adresse, samplePartnerCode)) # remove useless columns

# New ID variables
soilchem_raw$PlotID <- paste(substr(soilchem_raw$PlotID, 1, 3), substr(soilchem_raw$PlotID, 4, 4), sep= "-P") #Create PlotID same format as other sheets
soilchem_raw$SiteID <- substr(soilchem_raw$PlotID, 1, 3) #Create SiteID

#
## Data cleaning - New R object

soilchem_full <- soilchem_raw

#
## Char var - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(soilchem_full$SiteID) # 3 samples per site - validated - but missing OV1 site

# Adding OV1-P1 replicate
soilchem_full <- soilchem_full |> 
  add_row(PlotID = 'OV1-P1',
          SoilType = 13,
          ClayCategory = 1,
          LOI = 25.8,
          SoilDensity_kg.L = 0.59,
          Humus_percentDM = 25.8,
          HumusCategory = 5,
          pH = 5.4,
          P.Al_mg.100g = 3,
          K.Al_mg.100g = 14,
          Mg.Al_mg.100g = 20,
          Ca.Al_mg.100g = 31,
          Na.Al_mg.100g = 6,
          DryMatter_percent = 0.5,
          TotalN_percentDM = 0.44,
          SiteID = 'OV1'
          )

# Adding OV1-P2 replicate
soilchem_full <- soilchem_full |> 
  add_row(PlotID = 'OV1-P2',
          SoilType = 13,
          ClayCategory = 1,
          LOI = 26.0,
          SoilDensity_kg.L = 0.74,
          Humus_percentDM = 26.0,
          HumusCategory = 5,
          pH = 4.9,
          P.Al_mg.100g = 3,
          K.Al_mg.100g = 13,
          Mg.Al_mg.100g = 15,
          Ca.Al_mg.100g = 21,
          Na.Al_mg.100g = 5,
          DryMatter_percent = 91.4,
          TotalN_percentDM = 0.51,
          SiteID = 'OV1'
  )

# Adding OV1-P3 replicate
soilchem_full <- soilchem_full |> 
  add_row(PlotID = 'OV1-P3',
          SoilType = 14,
          ClayCategory = 1,
          LOI = 43.2,
          SoilDensity_kg.L = 0.38,
          Humus_percentDM = 43.2,
          HumusCategory = 6,
          pH = 5.0,
          P.Al_mg.100g = 3,
          K.Al_mg.100g = 12,
          Mg.Al_mg.100g = 14,
          Ca.Al_mg.100g = 16,
          Na.Al_mg.100g = 5,
          DryMatter_percent = 85.4,
          TotalN_percentDM = 0.72,
          SiteID = 'OV1'
  )

# Plot ID
#soilchem_full[duplicated(soilchem_full$PlotID),] # Unique plot ID - validated

# Soil type - make categories explicit
#unique(soilchem_full$SoilType) # 6 categories
soilchem_full <- soilchem_full |>  
  mutate(SoilType = ifelse(
    SoilType == 2, "Medium_sand", ifelse(
      SoilType == 3, "Fine_sand", ifelse(
        SoilType == 5, "Silty_medium_sand", ifelse(
          SoilType == 6, "Silty_fine_sand", ifelse(
            SoilType == 13, "Mineral_mixed_humus_soil", "Organic_soil"
            ))))))

# Clay category - make categories explicit
#unique(soilchem_full$ClayCategory) # 2 categories
soilchem_full <- soilchem_full |>  
  mutate(ClayCategory = ifelse(
    ClayCategory == 1, "0-5%", "5-10%"
  ))

# Humus category - make categories explicit
#unique(soilchem_full$HumusCategory) # 4 categories
soilchem_full <- soilchem_full %>% 
  mutate(HumusCategory = ifelse(
    HumusCategory == 5, "Mineral_mixed_humus", ifelse(
      HumusCategory == 6, "Organic_soil", "Moderately_humus-rich"
    )))

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- soilchem_full |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # max calcium over 1000 mg/100g ?

# LOI
#soilchem_full[is.na(soilchem_full$LOI),] # no NA
#hist(soilchem_full$LOI) # range from 0 to 90 -> very wide, but include both heathland and grassland. No visible outlier. Distribution a bit hectic

# Soil density
#soilchem_full[is.na(soilchem_full$SoilDensity_kg.L),] # no NA
#hist(soilchem_full$SoilDensity_kg.L) # range from 0 to 1.4 -> quite wide, but include both heathland and grassland. No visible outlier. Distribution a bit hectic

# Percent of humus in dry matter
#soilchem_full[is.na(soilchem_full$Humus_percentDM),] # no NA
#hist(soilchem_full$Humus_percentDM) # range from 0 to 90 -> matching with LOI

# pH
#soilchem_full[is.na(soilchem_full$pH),] # no NA
#hist(soilchem_full$pH) # range from 4 to 7, Normal distribution -> one outlier over 6.5
#filter(soilchem_full, pH>6.5) # OC2-P1 with the calcium outlier -> should be removed

# Phosphorus
#soilchem_full[is.na(soilchem_full$P.Al_mg.100g),] # no NA
#hist(soilchem_full$P.Al_mg.100g) # range from 0 to 40, Poisson distribution -> check high values
#filter(soilchem_full, P.Al_mg.100g>20) # 10 plots among 5 sites over 20 mg/100g
# 2 sites with all values over 20 (IS4, OC2)
# other sites (IC2, OC5, OG4), values not to far from other plots

# Potassium
#soilchem_full[is.na(soilchem_full$K.Al_mg.100g),] # no NA
#hist(soilchem_full$K.Al_mg.100g) # range from 0 to 30, Normal distribution -> check high values
#filter(soilchem_full, K.Al_mg.100g>20) # 4 plots among 2 sites (IC2, OC1) over 20 mg/100g -> coherent with other values

# Magnesium
#soilchem_full[is.na(soilchem_full$Mg.Al_mg.100g),] # no NA
#hist(soilchem_full$Mg.Al_mg.100g) # range from 0 to 35, Normal distribution -> check high values
#filter(soilchem_full, Mg.Al_mg.100g>20) # 4 plots among 2 sites (IC2, OC1) over 20 mg/100g, same as for Potassium

# Calcium
#soilchem_full[is.na(soilchem_full$Ca.Al_mg.100g),] # no NA
#hist(soilchem_full$Ca.Al_mg.100g) # range from 0 to 1000, one clear outlier
#filter(soilchem_full, Ca.Al_mg.100g>1000) # OC2-P1, not coherent with other samples -> to be removed
#filter(soilchem_full, Ca.Al_mg.100g>200) # 2 plots from same site (IC5) over 200 mg/100g

# Sodium
#soilchem_full[is.na(soilchem_full$Na.Al_mg.100g),] # no NA
#hist(soilchem_full$Na.Al_mg.100g) # range from 0 to 21, one clear outlier over 20
#filter(soilchem_full, Na.Al_mg.100g>20) # OC2-P1, same as Calcium -> to be removed
#filter(soilchem_full, Na.Al_mg.100g>12) # IS4-P3 & IS5-P2 -> coherent with rest of the samples

# Percent Dry Matter
#soilchem_full[is.na(soilchem_full$DryMatter_percent),] # no NA
#hist(soilchem_full$DryMatter_percent) # range from 10 to 100, distribution a bit hectic

# Total N in percent dry matter
#soilchem_full[is.na(soilchem_full$TotalN_percentDM),] # no NA
#hist(soilchem_full$TotalN_percentDM) # range from 0.2 to 2.2, distribution a bit hectic

#
## Data filtering/removal

# OC2-P1 outlier in several parameter -> farmer fertilizes in spring and summmer, maybe samples taken on a chunk
soilchem_full <- subset(soilchem_full, !PlotID == "OC2-P1")


## Export clean data in new excel file

write_csv(soilchem_full, "data/cleandata/NBR_FullSoilChem.csv")



#### Plant species community ####

## Description

# Plant community table, with species as rows and quadrats as columns.

#
## Check species names

#table(vege_raw$Species) # several name repetitions + some unknown species

# Liste species names
unique(vege_raw$Species)

# Correction latin name
vege_raw <- vege_raw |>
  mutate(Species = dplyr::recode(Species, "Antittrichia curtipendula" = "Antitrichia curtipendula")) |> 
  mutate(Species = dplyr::recode(Species, "Antrichum undulatum" = "Atrichum undulatum")) |> 
  mutate(Species = dplyr::recode(Species, "Braen sp US5-P2" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Bryum sp" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Carex ovalis" = "Carex leporina")) |>
  mutate(Species = dplyr::recode(Species, "Carex sp2 US2-P2" = "Carex sp")) |> 
  mutate(Species = dplyr::recode(Species, "Chamaepericlymenum suecica" = "Chamaepericlymenum suecicum")) |> 
  mutate(Species = dplyr::recode(Species, "Chiloscytinus pallenscens" = "Chiloscyphus pallescens")) |> 
  mutate(Species = dplyr::recode(Species, "Cuppy moss UC1" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Gnaphlium sylvaticum" = "Gnaphalium sylvaticum")) |> 
  mutate(Species = dplyr::recode(Species, "Golden curly fern UG2-P1" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Hieracium sp IG3-P1" = "Hieracium sp")) |> 
  mutate(Species = dplyr::recode(Species, "Kystkransmose" = "Rhytidiadelphus loreus")) |> 
  mutate(Species = dplyr::recode(Species, "Liverwort UG2-P2" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Luminous liverwort US1-P3" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Moss ND1 UG1-P1" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Mystery Carex OV1-P1" = "Carex sp")) |> 
  mutate(Species = dplyr::recode(Species, "Mystery grass OC5-P1" = "Poaceae sp")) |> 
  mutate(Species = dplyr::recode(Species, "Orchid sp" = "Orchidaceae sp")) |> 
  mutate(Species = dplyr::recode(Species, "Plagomnium medium" = "Plagiomnium medium")) |> 
  mutate(Species = dplyr::recode(Species, "Polytrichiastrum alpinum" = "Polytrichastrum alpinum")) |> 
  mutate(Species = dplyr::recode(Species, "Polytrichum juniperirinum" = "Polytrichum juniperinum")) |> 
  mutate(Species = dplyr::recode(Species, "Rhacomitrium lanuginosum" = "Racomitrium lanuginosum")) |> 
  mutate(Species = dplyr::recode(Species, "Rhacomitrium canescens" = "Racomitrium canescens")) |> 
  mutate(Species = dplyr::recode(Species, "Rhacomitrium sudeticum" = "Racomitrium sudeticum")) |> 
  mutate(Species = dplyr::recode(Species, "Rhacomitrium affine" = "Racomitrium affine")) |> 
  mutate(Species = dplyr::recode(Species, "Saniona uncinata" = "Sanionia uncinata")) |> 
  mutate(Species = dplyr::recode(Species, "Skinny pleurocarp IS5-P1" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Sphagunm compactum" = "Sphagnum compactum")) |>
  mutate(Species = dplyr::recode(Species, "Sphagnum sp UG2-P2" = "Bryophyte sp")) |> 
  mutate(Species = dplyr::recode(Species, "Straminergon straminergon" = "Straminergon stramineum")) |> 
  mutate(Species = dplyr::recode(Species, "Unknown shrub US2-P1" = "Woody sp"))

#
## Data cleaning and manipulation

# New R object
vege_full <- vege_raw

# Remove full NA rows
vege_full <- vege_full[rowSums(is.na(vege_full)) != ncol(vege_full), ] #remove full NA rows

# Aggregate repetitions
vege_full <- vege_full |>  
  group_by(Species) |> 
  summarise_if(is.numeric, sum, na.rm = TRUE) |>  
  distinct()
#sort(table(vege_full$Species)) #no doubletons remaining
# from 676 to 674 variables -> 2 quadrats removed?

# Remove non-identified species
vege_full <- subset(vege_full, Species != "Bryophyte sp" &
           Species != "Woody sp" &
           Species != "Hieracium sp" &
           Species != "Carex sp" &
           Species != "Poaceae sp" &
           Species != "Orchidaceae sp")

#
## Table formatting for vegan (long table)
vege_full <- vege_full |> 
  pivot_longer(cols = c(-Species), names_to = "SampleID", values_to = "Abundance") |> 
  mutate(PlotID = substr(SampleID, 1, 6)) |>
  mutate(SiteID = substr(SampleID, 1, 3))


## Export clean data in new excel file

write_csv(vege_full, "data/cleandata/NBR_FullPlantComm.csv")



#### Beetle families community ####

## Description

# Arthropod community table, with families/orders as columns and pitfall traps as rows

#
## Summary - Check table size, list of variables, variable types (num/chr)

#str(arthro_main) # Sampling period a bit messy
#str(arthro_sup) # Need to check Latin names for families

#
## Character cleaning, Common ID and correct Latin names for merging

# Character cleaning
names(arthro_sup) <- gsub(" ", "", names(arthro_sup))
names(arthro_main) <- gsub("-", "", names(arthro_main))

# Common ID
names(arthro_main) <- gsub("Site", "SiteID", names(arthro_main))
names(arthro_main) <- gsub("Pitfall_ID", "SampleID", names(arthro_main))
names(arthro_sup) <- gsub("Site", "SiteID", names(arthro_sup))
names(arthro_sup) <- gsub("PitfallID", "SampleID", names(arthro_sup))
names(arthro_sup) <- gsub("Samplingperiod", "Sampling_period", names(arthro_sup))

# Latin names
names(arthro_main) <- gsub("Scarabidae", "Scarabaeidae", names(arthro_main)) #rename families
names(arthro_main) <- gsub("Ptilidae", "Ptiliidae", names(arthro_main)) #rename families
names(arthro_sup) <- gsub("Curculinoideae", "Curculionidae", names(arthro_sup)) #rename families
names(arthro_sup) <- gsub("Chysomelidae", "Chrysomelidae", names(arthro_sup)) #rename families

# Binding main and complementary tables
arthro_raw <- full_join(arthro_main, arthro_sup)

#
## Summary sampling area - Check table size, list of variables, variable types (num/chr)

#str(arthro_raw) # all good except from period of sampling which has messy format - also missing plotID
arthro_raw$PlotID <- substr(arthro_raw$SampleID, 1, 6) #arrange plot ID to match other files

#
## Data cleaning - New R object + removal variables redundant with other dataset

arthro_full <- subset(arthro_raw, select = -c(Year, Habitat, Livestock, Storage))

#
## Char var - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
length(table(arthro_full$SiteID)) #39 sites, missing 6
table(arthro_full$SiteID) # 12 samples per site - sites missing are OV1, US3, US4, UG1, UG2, UC1

# Plot ID
#table(arthro_full$PlotID) # 4 samples per plot - validated

# Bulk density core ID
#arthro_full[duplicated(arthro_full$SampleID),] # Unique ID for core - validated

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- arthro_full |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # no obvious outliers, some loners

# Tot Beetle
#arthro_full[is.na(arthro_full$Beetle),] # 9 NAs -> pitfall traps crushed in OC5 (5), OC4 (1), IC5 (1) or water overloaded in US2 (2)
#hist(arthro_full$Beetle) # Poisson distribution, skewed.

# Distribution of beetle families
#hist(arthro_full$Staphylinidae) # Poisson, skewed from 0 to over 1000
#hist(arthro_full$Carabidae) # Poisson, skewed from 0 to 15
#hist(arthro_full$Hydrophilidae) # Poisson, skewed, from 0 to 130
#hist(arthro_full$Scarabaeidae) # Poisson, skewed, from 0 to 35
#hist(arthro_full$Ptiliidae) # Poisson, skewed from 0 to 400
#hist(arthro_full$Curculionidae) # Low abundance (from 1 to 4), but 50 samples have at least one
#hist(arthro_full$Elateridae) # Low abundance (from 1 to 4), but 50 samples have at least one
#hist(arthro_full$Leiodidae) # Only two loners -> not to be taken in account in the analysis
#hist(arthro_full$Rhizophagidae) # Only two loners -> not to be taken in account in the analysis
#hist(arthro_full$Silphidae) # Poisson, skewed from 0 to 35
#hist(arthro_full$Histeridae) # Very low abundance (10 pitfall traps with 1-2) -> not to be taken in account
#hist(arthro_full$Geotrupidae) # Low abundance (from 1 to 3), but 20 samples have at least one
#hist(arthro_full$Chrysomelidae) # Only four loners -> not to be taken in account
#hist(arthro_full$Dascillidae) # Only one loner -> not to be taken in account
#hist(arthro_full$Erotylidae) # Only one loner -> not to be taken in account

#
## Table formatting for vegan

# Replace NAs by zeros
arthro_full <-  mutate_if(is.numeric, ~replace(., is.na(.), 0))

# Long table beetle only
beetle_full <- subset(arthro_full, select = c(SiteID, SampleID, Staphylinidae, Carabidae, Hydrophilidae, Scarabaeidae, Ptiliidae, Curculionidae, Elateridae, Rhizophagidae, Leiodidae, Silphidae, Histeridae, Geotrupidae, Chrysomelidae, Dascillidae, Other, Erotylidae))
beetle_full <- beetle_full |> 
  pivot_longer(
    cols = c(Staphylinidae, Carabidae, Hydrophilidae, Scarabaeidae, Ptiliidae, Curculionidae, Elateridae, Rhizophagidae, Leiodidae, Silphidae, Histeridae, Geotrupidae, Chrysomelidae, Dascillidae, Other, Erotylidae),
    names_to = "BeetleFamilies", 
    values_to = "BeetleFam_abundance") |> 
  mutate(PlotID = substr(SampleID, 1, 6)) |>
  mutate(SiteID = substr(SampleID, 1, 3)) |> 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) # replace NAs by zeros


# Long table all arthropods
arthro_full <- arthro_full |> 
  pivot_longer(
    cols = c(Staphylinidae, Carabidae, Hydrophilidae, Scarabaeidae, Ptiliidae, Curculionidae, Elateridae, Rhizophagidae, Leiodidae, Silphidae, Histeridae, Geotrupidae, Chrysomelidae, Dascillidae, Other, Erotylidae),
    names_to = "BeetleFamilies", 
    values_to = "BeetleFam_abundance") |> 
  pivot_longer(
    cols = c(Beetle, Spider, Diptera, Hemiptera, Opilion, Worm, Slug, Snail, Cloporte, Millipoda, Orthoptera, Hymenoptera),
    names_to = "Orders", 
    values_to = "Order_abundance") |>
  mutate(PlotID = substr(SampleID, 1, 6)) |>
  mutate(SiteID = substr(SampleID, 1, 3)) |> 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))

# Export clean data in new excel file
write_csv(beetle_full, "data/cleandata/NBR_FullBeetleComm.csv")
write_csv(arthro_full, "data/cleandata/NBR_FullArtComm.csv")



#### Mesofauna abundance data ####

## Description

## List of variables

# [1] Name of the counter
# [2] Sample ID
# [3] Acari abundance
# [4] Collembola abundance
# [5] Site ID
# [6] Other observer ?

#
## Summary - Check table size, list of variables, variable types (num/chr)

#str(mesobio_raw) # need to rename comments, siteID and plotID should also be added
#str(soilmeso_raw) # need to rename comments, siteID and plotID should also be added

#
## Character cleaning, Common ID and correct Latin names for merging

# Common ID
names(mesobio_raw) <- gsub("Sample_name_on_pot", "SampleID", names(mesobio_raw))
names(mesobio_raw) <- gsub("...5", "Comments", names(mesobio_raw))
names(soilmeso_raw) <- gsub("Site", "SiteID", names(soilmeso_raw))
names(soilmeso_raw) <- gsub("PlotID", "SampleID", names(soilmeso_raw))
names(soilmeso_raw) <- gsub("Core_depth", "CoreDepth", names(soilmeso_raw))
names(soilmeso_raw) <- gsub("\\)", "", names(soilmeso_raw))
names(soilmeso_raw) <- gsub("\\(", "", names(soilmeso_raw))

# Check ID coding
table(mesobio_raw$SampleID) # ØY experiment samples merged in, need to extract ØY-GK as OV1
mesobio_raw <- mesobio_raw |> 
  mutate(SampleID = dplyr::recode(SampleID, "ØY-R1-T2-D1" = "OV1-P1-D1")) |> 
  mutate(SampleID = dplyr::recode(SampleID, "ØY-R1-T3-D1" = "OV1-P2-D1")) |> 
  mutate(SampleID = dplyr::recode(SampleID, "ØY-GK-T4-D1" = "OV1-P3-D1"))
#table(soilmeso_raw$SiteID) # all good, 12 cores per site

# Add PlotID and SiteID
mesobio_raw <- mesobio_raw |> 
  mutate(PlotID = substr(SampleID, 1, 6)) |>
  mutate(SiteID = substr(SampleID, 1, 3))
soilmeso_raw <- soilmeso_raw |> 
  mutate(PlotID = substr(SampleID, 1, 6))

#
## Char var Site Info - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(mesobio_raw$SiteID) # at least 3 per sheep (S or V) sites -> validated

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- mesobio_raw |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # All good

# Check distribution of quantitative variable
#hist(mesobio_raw$Acari) # Poisson distribution
#hist(mesobio_raw$Collembola) # Poisson distribution
#hist(soilmeso_raw$CoreDepth_cm) # Most around 14 cm -> validated

#
## New variable - abundance per soil area with correction soil volume

# Extraction survey data only
mesobio_full <- mesobio_raw |> 
  filter(SiteID != "ØY-")

# Merging datasets according to sorted fauna
mesobio_full <- left_join(mesobio_full, soilmeso_raw)
mesobio_full <- subset(mesobio_full, select = c(SampleID, Acari, Collembola, PlotID, SiteID, CoreDepth_cm))

# New variable with correction for soil volume
mesobio_full <- mesobio_full |> 
  # corrected abundance = (measured_abundance*standard_coreheight)/measured_coreheight
  mutate(Acari.m2 = ((Acari*14)/CoreDepth_cm)/(3.14*(0.105/2)^2)) |> 
  mutate(Collembola.m2 = ((Collembola*14)/CoreDepth_cm)/(3.14*(0.105/2)^2))

#
## Export clean data in new excel file

write_csv(mesobio_full, "data/cleandata/NBR_FullMesobio.csv")



#### ABOVEGROUND BIOMASS ####

## Description

## List of variables

# [1] Name of the observer who sorted and weighted the biomass
# [2] Sample ID
# [3] Plant functional type
# [4] Dry weight of the biomass with the bag
# [5] Dry weight of the biomass without the bag
# [6] Type of bag
# [7] Comments

#
## Summary  - Check table size, list of variables, variable types (num/chr)

#str(biomass_raw) # all good

#
## Character cleaning, Common ID and correct Latin names for merging

# Character cleaning & change variable name
names(biomass_raw) <- gsub(" ", "", names(biomass_raw))
names(biomass_raw) <-  gsub("\\(", "_", names(biomass_raw))
names(biomass_raw) <-  gsub("\\)", "", names(biomass_raw))
names(biomass_raw) <-  gsub("\\+", "And", names(biomass_raw))
names(biomass_raw) <- gsub("Who", "ProcessedBy", names(biomass_raw))
names(biomass_raw) <- gsub("Dryweight_g", "DWbiomass_g", names(biomass_raw))


# Check ID coding
#table(biomass_raw$SampleID) # all good

# Add PlotID and SiteID
biomass_raw <- biomass_raw |> 
  mutate(PlotID = substr(SampleID, 1, 6)) |>
  mutate(SiteID = substr(SampleID, 1, 3))

#
## Data cleaning - New R object

biomass_full <- subset(biomass_raw, select = -c(BiomassAndbag_g, BagType, ProcessedBy, Comments))

#
## Char var Site Info - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
table(biomass_full$SiteID) # all sheep sites with between 12 and 77 bags -> need to check functional types

# Functional type
#unique(biomass_full$FunctionalType) # Need clear categories
biomass_full <- biomass_full |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Grasses" = "monocotyledons")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Graminoids" = "monocotyledons")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "graminoids" = "monocotyledons")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Graminioids" = "monocotyledons")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Graminiods" = "monocotyledons")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Herbs" = "forbs")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Forbs" = "forbs")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Ferns" = "ferns")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "fern" = "ferns")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Fern" = "ferns")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Woody" = "woody")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Mosses" = "bryophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "moss" = "bryophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Moss" = "bryophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Bryo" = "bryophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "mosses" = "bryophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "leaf litter" = "litter")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Litter" = "litter")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Lichens" = "lichens")) |>
  mutate(FunctionalType = dplyr::recode(FunctionalType, "lichen" = "lichens")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Club mosses" = "lycophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "lycopodium" = "lycophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Lycopodium" = "lycophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Huperzia selago" = "lycophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Diphasiastrum alpinum" = "lycophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "Selaginella selaginoides" = "lycophytes")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "lycophytes" = "cryptogams")) |> 
  mutate(FunctionalType = dplyr::recode(FunctionalType, "bryophytes" = "cryptogams"))
  
  # 9 categories

# Check categories
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "IS1")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "IS2")) # 3 D1s, all functional groups but woody 0 for IS2-P1-D1 -> no samples found for woody
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "IS3")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "IS4")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "IS5")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "IV1")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS1")) # 3 D1s, all functional groups -> validated
xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS2")) # all site, but some samples including D1s not finished sorting -> need to standardise or exclude samples with mix bryo/litter
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS3")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS4")) # 3 D1s, all functional groups -> validated
xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS5")) # one non ID woody sample which should be removed - even if not all D1s, all samples with same functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS6")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS7")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS8")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OS9")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OV1")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "OV2")) # all Ds, all functional groups
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "US1")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "US2")) # 3 D1s, all functional groups -> validated 
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "US3")) # all Ds, all functional groups
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "US4")) # all Ds, all functional groups
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "US5")) # 3 D1s, all functional groups -> validated
#xtabs(DWbiomass_g ~ SampleID + FunctionalType, data = filter(biomass_full, SiteID == "US6")) # 3 D1s, all functional groups -> validated

#
## Remove mislabelled & non fully sorted samples

biomass_full <- filter(biomass_full, SampleID != "OS5-?-?" &
                         SampleID != "OS2-P2-D1" &
                         SampleID != "OS2-P3-D2")

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- biomass_full |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # some negative values -> should be 

# Check distribution of quantitative variable
hist(biomass_full$DWbiomass_g) # quite a few negative values
#filter(biomass_full, DWbiomass_g < 1) # 5 negative values

# Check for duplicates in weights
dupli <- biomass_full |>
  group_by(SampleID, PlotID, SiteID, FunctionalType) |> 
  summarise(n = n(), .groups = "drop") |> 
  filter(n > 1L)
unique(dupli$SampleID) # 8 duplicates

# OS2-P1-D1 monocotyledons, OS2-P2-D3 forbs & monocotyledons
# filter(biomass_full, SampleID == "OS2-P1-D1")
# filter(biomass_full, SampleID == "OS2-P2-D3")
# filter(biomass_full, SiteID == "OS2" & FunctionalType == "monocotyledons")
# filter(biomass_full, SiteID == "OS2" & FunctionalType == "forbs")
# Check on comments -> sample sorted in two half separatly -> weights should be summed

# US4 lycophytes
# filter(biomass_full, SampleID == "US4-P1-D1")
# filter(biomass_full, SampleID == "US4-P1-D2")
# filter(biomass_full, SampleID == "US4-P1-D3")
# filter(biomass_full, SampleID == "US4-P1-D4")
# filter(biomass_full, SampleID == "US4-P2-D4")
# filter(biomass_full, SampleID == "US4-P3-D1")
# filter(biomass_full, SiteID == "US4" & FunctionalType == "lycophytes") # duplicates due to renaming -> need to sum up all rows

# Summarise duplicates by sum
biomass_full <- biomass_full |> 
  group_by(SampleID, PlotID, SiteID, FunctionalType) |> 
  summarise(DWbiomass_g = sum(DWbiomass_g, na.rm = TRUE)) |> 
  ungroup()

# New biomass value per m^2^
biomass_full <- biomass_full |> 
  mutate(Biomass.m2 = DWbiomass_g*4)
biomass_full <- subset(biomass_full, select = -c(DWbiomass_g))

#
## Prepare data for vegan

# Select 3 replicates per site
biomass_full <- biomass_full |> 
  group_by(SiteID, PlotID) |>
  pivot_wider(names_from = FunctionalType, values_from = Biomass.m2) |> 
  # Select randomly one row which match unique combination of site & plot IDs
  slice(1) |>
  ungroup() |>
  pivot_longer(cols = c(-SiteID, -PlotID, -SampleID), names_to = "FunctionalType", values_to = "Biomass.m2")

# Export clean data in new excel file
write_csv(biomass_full, "data/cleandata/NBR_FullBiomass.csv")


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


