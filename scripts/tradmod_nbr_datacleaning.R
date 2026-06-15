##
#### DESCRIPTION ####
##
## Purpose: Cleaning script for data paper "Grazing-dependent communities in the fjord region - data from an observational study in the Nordhordland UNESCO Biosphere Reserve"
## Author: Morgane KERDONCUFF
## ORCID: 0000-0003-2223-1857
## github
## Date created: 11/2022
## Last time modified: 04/2026
## Project: TradMod
## Funding: Norges forskningsråd (NFR)
## Institution: University of Bergen, Norway
##
##

#### PACKAGES ####

library(tidyverse) # R language
library(janitor) # Data cleaning
library(readxl) # Read xl files
library(lubridate) # Standard date
library(purrr) # Merge tables

#### RAW DATA ####

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
## General information about grazing fields selected as study sites
site_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SiteInfo")
## Management data on each study site collected during farmer interviews
management_raw <- read_excel(path = "data/rawdata/NBR_RawFarmerSurvey.xlsx", sheet="Farms Information_R")

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
# colnames(site_description)[apply(site_description, 2, anyNA)] # numberAnimalAdult; numberAnimalYoung; fieldAreaHA; farmGrazingAreaHa
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

#### 20x20 SAMPLING AREA ####

# Raw datasets
## Description of the representative 20 m x 20 m sampling areas
sampling_area_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="20mX20m")

# Desired variables
## siteID - Field identification code for data collection
## eventDate - Date of data collection
## numberLivestockPaths - Number of paths created by the livestock in the sampling area
## lengthLivestockPathM - Total length of livestock path in the sampling area in meters
## elevationMasl - Elevation at the highest point of the sampling area in meters above sea level
## slopeAngleDegree - Slope angle from the highest point to the lowest point of the sampling area in degree
## slopeAspectDegree - Slope aspect from the highest point to the lowest point of the sampling area in degree
## percentRock - Estimated percentage cover of exposed rock in the sampling area
## percentMud - Estimated percentage cover of visible mud in the sampling area
## percentTreeTallShrubs - Estimated percentage cover of trees and shrubs over 1 m in the sampling area
## percentLowShrubs - Estimated percentage cover of shrubs under 1 m in the sampling area
## percentForbs - Estimated percentage cover of forbs in the sampling area
## percentMonocotyledons - Estimated percentage cover of monocotyledons (grasses, rushes, sedges) in the sampling area
## percentBryophytes - Estimated percentage cover of bryophytes (mosses, liverworts) in the sampling area
## percentFerns - Estimated percentage cover of ferns in the sampling area
## percentLichens - Estimated percentage cover of lichens in the sampling area

# Variable names

## R-friendly with janitor package
sampling_area_raw <- sampling_area_raw %>% 
  clean_names("lower_camel")

## Consistent & FAIR
names(sampling_area_raw) <- gsub("siteId", "siteID", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("recordingDate", "eventDate", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("numberPaths", "numberLivestockPaths", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("totalLengthPathM", "lengthLivestockPathM", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("elevationMax", "elevationMasl", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("generalSlope", "slopeAngleDegree", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("aspectDegree", "slopeAspectDegree", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("percentHerbs", "percentForbs", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("percentFern", "percentFerns", names(sampling_area_raw))
names(sampling_area_raw) <- gsub("percentLichen", "percentLichens", names(sampling_area_raw))

# Dataset

## Selection desirable variables
sampling_area <- subset(
  sampling_area_raw, select = -c(
    # anonymous
    team,
    # redundant with site_description
    latitude1,
    longitude1,
    latitude2,
    longitude2,
    # only one elevation variable necessary
    elevationMin,
    # only aspect degree necessary
    aspect,
    comments
  )
)

## Variable types (num/chr)
# str(sampling_area) #validated
## Standard date according to ISO 8601-1:2019
sampling_area$eventDate <- as.POSIXct(sampling_area$eventDate, format = "%d.%m.%Y")
## Duplicate check
# get_dupes(sampling_area) #validated

## Character variables - desirable categories, NAs, misprints

### Consistent lower camel
sampling_area <- sampling_area %>% 
  mutate_if(is.character, tolower)

### Categories & distribution
# table(sampling_area$siteID) #validated - uc1 site (bog) to be removed
sampling_area <- filter(sampling_area, siteID != "uc1")

## Numeric variables - min/max, distribution, potential outliers

## Min/max
# test <- sampling_area |>
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
#   transpose() #validated - need further check for max herbs 97% & max lichens 20%

### NA check
# colnames(sampling_area)[apply(sampling_area, 2, anyNA)] # col: numberLivestockPaths, lengthLivestockPaths & all percent cover
sampling_area[!complete.cases(sampling_area),] # rows: us1, ug1, oc4

### Variable distribution & outliers
table(sampling_area$numberLivestockPaths) # dominance 0, variable to be taken out
hist(sampling_area$lengthLivestockPathM) # dominance 0, variable to be taken out
# hist(sampling_area$elevationMasl) # Poisson distribution, no outlier
# hist(sampling_area$slopeAngleDegree) # Normal distribution, no outlier
hist(sampling_area$slopeAspectDegree) # Uneven distribution, no outlier
# hist(sampling_area$percentRock) # Poisson distribution, further check for sites over 7%
# sampling_area[sampling_area$percentRock>7,] # subalpine heathlands (ug2, us3, us5), coastal heatland (ov1) & fjord at higher elevation (ig3)
# hist(sampling_area$percentMud) # very skewed Poisson distribution, no outlier
# hist(sampling_area$percentTreesTallShrubs) # Skewed Poisson distribution, no outlier
# hist(sampling_area$percentLowShrubs) # Poisson distribution, further check all grassland sites should be under 10%
# sampling_area[sampling_area$percentLowShrubs>10,] # 11 heathland sites
hist(sampling_area$percentForbs) # Poisson distribution, further check for sites >50%
sampling_area[sampling_area$percentForbs>50,] # os1, oc1, ig1, ig2, is2, iv1, ic1, og2, ic4 -> all first year/starting sites, check on vegetation quadrats + site & plot pictures
# OS1 80% - average 20% & no cover over 55% in quadrats, estimation from pictures 35%-40%
# OC1 80% - average 25% & no cover over 50% in quadrats, estimation from pictures 10%-15% 
# IG1 97% - average 20% & no cover over 30% in quadrats, estimation from pictures 15%-20%
# IS2 80% - average 70% & no cover over 90% in quadrats, estimation from pictures 45%-50%
# IC1 70% - average 45% & no cover over 70% in quadrats, estimation from pictures 60%-65%
hist(sampling_area$percentMonocotyledons) # further check for sites with odd forb distribution
hist(sampling_area$percentBryophytes) # uneven distribution, further check for sites with odd forb distribution
# hist(sampling_area$percentLichens) # highly skewed Poisson distribution, one site over 10%
# sampling_area[sampling_area$percentLichens>10,] # US4 in subalpine area, average of 12% & max 24% in quadrats - validated

# Add/remove variables

## Removal numberLivestockPaths and lengthLivestockPathM
sampling_area <- subset(sampling_area, select = -c(numberLivestockPaths, lengthLivestockPathM))

## Heat Load Index
# sampling_area <- sampling_area |> 
#   mutate(heatLoadIndex = cos(slopeAspectDegree-225)*tan(slopeAngleDegree))
# hist(sampling_area$heatLoadIndex) # 3 outliers: one under 200, two over 100
# sampling_area[sampling_area$heatLoadIndex>100,] #OG4 & IS3 -> both 11 degree slope with SW & SE exposition
# sampling_area[sampling_area$heatLoadIndex<0,] #OS6 -> 11 degree slope with NE exposition

# Dataset export
write_csv(sampling_area, "data/cleandata/tradmod_nbr_samplingarea.csv")

#### NON-DESTRUCTIVE SUBPLOTS - GROUND COVER ####

# Raw datasets
## 1 m^2^ ground cover in non-destructive subplots (15 per site)
ground_cover_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SoilCover")

# Desired variables
## siteID - Field identification code for data collection
## plotID - Plot identification code for data collection
## subplotID - Subplot identification code for data collection
## percentBareGround - Estimated percentage cover of bare ground in the subplot
## percentRock - Estimated percentage cover of exposed rock in the subplot
## percentLitter - Estimated percentage cover of litter in the subplot
## percentDeadWood - Estimated percentage cover of dead wood in the subplot
## percentBryophytes - Estimated percentage cover of bryophytes in the subplot
## percentVascular - Estimated percentage cover of vascular plants in the subplot
## percentLichens - Estimated percentage cover of lichens in the subplot
## percentBlossom - Estimated percentage cover of flower blossom in the subplot
## blossomSp1-blossomSp5 - Plant species with blossoms in the subplot when applicable
## percentDung - Estimated percentage cover of dung pellets or pads in the subplot
## avgVegetationHeightCm - Average vegetation height in the subplot in cm
## maxVegetationHeightCm - Maximum vegetation height in the subplot in cm

# Variable names

## R-friendly with janitor package
ground_cover_raw <- ground_cover_raw %>% 
  clean_names("lower_camel")

## Consistent & FAIR
names(ground_cover_raw) <- gsub("site", "siteID", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("plotId", "subplotID", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("bareSoil", "percentBareGround", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("rocks", "percentRock", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("litter", "percentLitter", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("deadWood", "percentDeadWood", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("bryophytes", "percentBryophytes", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("vascular", "percentVascular", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("lichen", "percentLichens", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("blossomCover", "percentBlossom", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("dung", "percentDung", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("vgMeanHeightCm", "avgVegetationHeightCm", names(ground_cover_raw))
names(ground_cover_raw) <- gsub("vgMaxHeightCm", "maxVegetationHeightCm", names(ground_cover_raw))

## New plotID variable
ground_cover_raw <- ground_cover_raw %>% 
  mutate(plotID = substr(subplotID, 1, 6))

# Dataset

## Selection desirable variables
ground_cover <- subset(
  ground_cover_raw, select = -c(
    # redundant with sampling_area
    date,
    # species richness to be calculated from community data
    plantSpeciesRichness,
    comment
  )
)

## Variable types (num/chr)
# str(ground_cover) #validated
## Duplicate check
# get_dupes(ground_cover) #validated

## Character variables - desirable categories, NAs, misprints

### Consistent lower typo
ground_cover <- ground_cover %>%
  mutate_if(is.character, tolower)

### Categories & distribution
# table(ground_cover$siteID) #validated - uc1 site (bog) to be removed
ground_cover <- filter(ground_cover, siteID != "uc1")
# table(ground_cover$plotID) #validated
# table(ground_cover$blossomSp1) # latin names check
ground_cover <- ground_cover |> 
  # Cirsium palustre misprint
  mutate(blossomSp1 = dplyr::recode(blossomSp1, "cirsium palustris" = "cirsium palustre")) %>% 
  # Only Achillea millefolium was recorded in is1-p3-n3
  mutate(blossomSp1 = dplyr::recode(blossomSp1, "alchemilla millefolium" = "achillea millefolium"))
table(ground_cover$blossomSp2) # latin names check
ground_cover <- ground_cover |> 
  # Leontodon saxatile wrong id
  mutate(blossomSp2 = dplyr::recode(blossomSp2, "leontodon saxatile" = "leontodon autumnalis"))
# table(ground_cover$blossomSp3) #validated
# table(ground_cover$blossomSp4) #validated
# table(ground_cover$blossomSp5) #validated

## Numeric variables - min/max, distribution, potential outliers

## Min/max
# test <- ground_cover |>
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
#   transpose() #validated - need further check for max lichens 80%

### NA check
# colnames(ground_cover)[apply(ground_cover, 2, anyNA)] #validated - only blossom species ID

### Variable distribution & outliers
# hist(ground_cover$percentBareGround) # skewed Poisson distribution, check subplots > 20%
# filter(ground_cover, percentBareGround>20) #validated - subplots from recently burnt heathland ov1
# hist(ground_cover$percentRock) # skewed Poisson distribution
# hist(ground_cover$percentLitter) # skewed Poisson distribution, check subplots > 30%
# filter(ground_cover, percentLitter>30) #validated - subplots from burnt (ov1) & mountain heathlands (us2, ug1)
# hist(ground_cover$percentDeadWood) # skewed Poisson distribution, check subplots > 2%
# filter(ground_cover, percentDeadWood>2) #validated - subplots from is5, in the middle of a wood clearing
# hist(ground_cover$percentBryophytes) # Poisson distribution, no outliers
# hist(ground_cover$percentLichens) # Poisson distribution, one outlier above 40%
# filter(ground_cover, percentLichens>40) #validated - mountain site (us1-p1-n5) with high lichen cover
# hist(ground_cover$percentVascular) # Exponential distribution, no outlier
# hist(ground_cover$percentBlossom) # Skewed Poisson distribution, no outlier
# hist(ground_cover$percentDung) # Skewed Poisson distribution, check subplots > 10%
# filter(ground_cover, percentDung>10) #validated - cattle site (ic1)
# hist(ground_cover$avgVegetationHeightCm) # Poisson distribution, no outliers
# hist(ground_cover$maxVegetationHeightCm) # Normal distribution, no outliers

# Dataset export
write_csv(ground_cover, "data/cleandata/tradmod_nbr_groundcover.csv")

#### DESTRUCTIVE SUBPLOTS - SOIL PENETRATION TESTS ####

# Raw datasets
## Soil penetration tests in destructive subplots (2 per subplot, 24 per site)
soil_pene_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SoilPenetration")

# Desired variables
## siteID - Field identification code for data collection
## plotID - Plot identification code for data collection
## subplotID - Subplot identification code for data collection
## recordID - Identification code for individual penetration tests
## soilPenetrationDepthCm - Soil penetration depth of a sharpened metal rod of 43.4 cm length (diam. 2.2 cm; weight 1.3 kg) dropped from 1 m above ground in cm
## bedrockHit - If the penetration stick hit the bedrock during the test (yes/no)

# Variable names & structure

## R-friendly with janitor package
soil_pene_raw <- soil_pene_raw %>% 
  clean_names("lower_camel")

## Consistent & FAIR
names(soil_pene_raw) <- gsub("site", "siteID", names(soil_pene_raw))
names(soil_pene_raw) <- gsub("plotId", "subplotID", names(soil_pene_raw))

## New plot & record ID variables
soil_pene_raw <- soil_pene_raw %>%
  mutate(plotID = substr(subplotID, 1, 6)) %>%
  # replace left/right test by a record ID
  mutate(recordID = ifelse(
    leftRight == "left",
    paste(subplotID, "r1", sep = "-"),
    paste(subplotID, "r2", sep = "-")
  ))

# Dataset

## Removal redundant or unnecessary variables
soil_pene <- subset(
  soil_pene_raw, select = -c(
    # redundant with sampling_area
    date,
    # replaced by recordID
    leftRight,
    comments
  )
)

## Variable types (num/chr)
# str(soil_pene) #validated
## Duplicate check
# get_dupes(soil_pene) #validated

## Character variables - desirable categories, NAs, misprints

### Consistent lower typo
soil_pene <- soil_pene %>%
  mutate_if(is.character, tolower)

### Categories & distribution
# table(soil_pene$siteID) #validated - uc1 site (bog) to be removed
soil_pene <- filter(soil_pene, siteID != "uc1")
# table(soil_pene$plotID) #validated
# table(soil_pene$bedrockHit) # 39 failed tests over 1017 due to bedrock hit
table(filter(soil_pene, bedrockHit == "y")$siteID) # 8 sites with up to 12 failures

## Numeric variables - min/max, distribution, potential outliers

### Min/max
# test <- soil_pene |>
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
#   transpose() #validated - no visible height above maximum stick length

### NA check
# colnames(soil_pene)[apply(soil_pene, 2, anyNA)] #validated

### Variable distribution & outliers
# hist(soil_pene$visibleHeightCm) # Normal distribution, no outliers

# Add/remove variables

## New variable soilPenetrationDepth
soil_pene <- soil_pene %>% 
  mutate(soilPeneDepthCm = stickHeight - visibleHeightCm)

## Removal stickHeight and visibleHeightCm
soil_pene <- subset(soil_pene, select = -c(stickHeight, visibleHeightCm))

# Dataset export
write_csv(soil_pene, "data/cleandata/tradmod_nbr_soilpene.csv")

#### DESTRUCTIVE SUBPLOTS - BULK DENSITY & GRAVIMETRIC WATER CONTENT ####

# Raw datasets
# Bulk density in destructive subplots (3 per subplot, 36 per site)
soil_bulk_raw <- read_excel(path = "data/rawdata/NBR_RawBD.xlsx", na="NA")

# Desired variables
## siteID - Field identification code for data collection
## plotID - Plot identification code for data collection
## subplotID - Subplot identification code for data collection
## recordID - Identification code for individual bulk density cores
## coreVolCm3 - Core volume corrected for gaps in cm3
## weight0hG - Initial core weight after collection in g
## weightSatG - Water saturated weight of soil cores in g
## weight24hG - Core weight after 24h air drying in g
## weight48hG - Core weight after 48h air drying in g
## weightDryG - Core weight after 48h oven drying at 105C in g

# Variable names & structure

## R-friendly with janitor package
soil_bulk_raw <- soil_bulk_raw %>% 
  clean_names("lower_camel")

## Consistent & FAIR
names(soil_bulk_raw) <- gsub("site", "siteID", names(soil_bulk_raw))
names(soil_bulk_raw) <- gsub("plotId", "plotID", names(soil_bulk_raw))
names(soil_bulk_raw) <- gsub("bdCoreId", "recordID", names(soil_bulk_raw))
names(soil_bulk_raw) <- gsub("coreVol", "coreVolCm3", names(soil_bulk_raw))
names(soil_bulk_raw) <- gsub("w0", "weight0h", names(soil_bulk_raw))
names(soil_bulk_raw) <- gsub("wsat", "weightSatG", names(soil_bulk_raw))
names(soil_bulk_raw) <- gsub("w24H", "weight24hG", names(soil_bulk_raw))
names(soil_bulk_raw) <- gsub("w48H", "weight48hG", names(soil_bulk_raw))
names(soil_bulk_raw) <- gsub("wdry", "weightDryG", names(soil_bulk_raw))

## Consistent plot, subplot & record ID variables
soil_bulk_raw <- soil_bulk_raw %>%
  mutate(plotID = substr(recordID, 1, 6)) %>% 
  mutate(subplotID = substr(recordID, 1, 9)) %>% 
  group_by(subplotID) %>% 
  mutate(recordID = paste(subplotID, row_number(), sep = "-r")) %>% 
  ungroup()

# Dataset

## Removal redundant or unnecessary variables
soil_bulk <- subset(
  soil_bulk_raw, select = c(
    siteID,
    plotID,
    subplotID,
    recordID,
    # corrected volume of soil cores in cm3
    coreVolCm3,
    # soil core weight measurements
    weight0hG,
    weightSatG,
    weight24hG,
    weight48hG,
    weightDryG
  )
)

## Variable types (num/chr)
# str(soil_bulk) #validated
## Duplicate check
# get_dupes(soil_bulk) #validated

## Character variables - desirable categories, NAs, misprints

### Consistent lower typo
soil_bulk <- soil_bulk %>%
  mutate_if(is.character, tolower)

### Categories & distribution
# table(soil_bulk$siteID) #validated - uc1 site (bog) to be removed
soil_bulk <- filter(soil_bulk, siteID != "uc1")
# table(soil_bulk$plotID) #validated

## Numeric variables - min/max, distribution, potential outliers

# Min/max
# test <- soil_bulk |>
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
#   transpose() #quality check needed - min core volume low

 # min core volume quite low, negative values for water loss 24h and 48h, negative BD

### NA check
# colnames(soil_bulk)[apply(soil_bulk, 2, anyNA)] # all variable, check row identification
# soil_bulk[!complete.cases(soil_bulk),] # two missing records (is1-p3-d4-r2 & og1-p3-d2-r1) -> discarded due to lab incident
soil_bulk <- filter(soil_bulk, recordID != "is1-p3-d4-r2" & recordID != "og1-p3-d2-r1")

### Quality check

#### Soil core volume - low core volume affect weight measurements
# hist(soil_bulk$coreVolCm3) # Visible threshold around 50 cm3
# filter(soil_bulk, coreVolCm3<40) # 37 or 2% samples unfit
# filter(soil_bulk, coreVolCm3<45) # 105 or 7% samples unfit
# filter(soil_bulk, coreVolCm3<50) # 330 or 21% samples unfit -> cores should be minimum vol of 50 cm3
soil_bulk <- filter(soil_bulk, coreVolCm3 >= 50)

### Soil core weight - weight loss should be consistent with drying processes
qualitycheck <- filter(soil_bulk,
                        weightSatG - weight0hG < 0 |
                        weight24hG - weightSatG > 0 |
                        weight48hG - weight24hG > 0 |
                        weightDryG - weight0hG > 0)

### Variable distribution & outliers

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



#### DESTRUCTIVE SUBPLOTS - SOIL CHEMISTRY ####

# Raw datasets
# Soil chemistry 2019 (3 per plot)
soil_chem_2019 <- read.csv("data/rawdata/NBR_RawSoilChemistry2019.txt", sep=";")
# Soil chemistry 2020 (3 per plot)
soil_chem_2020 <- read.csv("data/rawdata/NBR_RawSoilChemistry2020.txt", sep=";")
# Soil chemistry 2020 complementary dry matter content
soil_chem_2020_DM <- read_excel(path = "data/rawdata/NBR_RawSoilChemistry2020bis.xls")

# Desired variables
## siteID - Field identification code for data collection
## plotID - Plot identification code for data collection
## soilType - Soil type according to sand, clay & humus content
## clayContentClass - Content of clay in soil in percent classes
## lossOnIgnitionPercentDM - Labile carbon content in soil in percent dry matter content
## soilDensityKgL - Soil density in kg per liter
## dryMatterPercent - Dry matter content in the soil in percent weight
## humusPercentDM - Content of humus in soil in percent dry matter content
## humusContentClass - Content of humus in soil in classes
## pH - Soil pH
## availablePMg100g - Concentration of phosphorus in mg in 100 g of soil
## availableKMg100g - Concentration of potassium in mg in 100 g of soil
## availableMgMg100g - Concentration of magnesium in mg in 100 g of soil
## availableCaMg100g - Concentration of calcium in mg in 100 g of soil
## availableNaMg100g - Concentration of sodium in mg in 100 g of soil
## totalNPercentDM - Nitrogen content in percent dry matter content

# Variable names & structure

## R-friendly with janitor package
soil_chem_2019 <- soil_chem_2019 %>% 
  clean_names("lower_camel")
soil_chem_2020 <- soil_chem_2020 %>% 
  clean_names("lower_camel")
soil_chem_2020_DM <- soil_chem_2020_DM %>% 
  clean_names("lower_camel")

## Consistent plotID
names(soil_chem_2019) <- gsub("provenummer", "plotID", names(soil_chem_2019))
names(soil_chem_2020) <- gsub("provenummer", "plotID", names(soil_chem_2020))
names(soil_chem_2020_DM) <- gsub("merking", "plotID", names(soil_chem_2020_DM))

# Dataset

## Variable list & types
# str(soil_chem_2019) # dryMatterPercent & totalNPercentDM missing -> in the PDF version of the document
# str(soil_chem_2020) # dryMatterPercent & totalNPercentDM missing -> dryMatterPercent in soil_chem_2020_DM & totalNPercentDM in the PDF version of the document
# str(soil_chem_2020_DM) # dryMatterPercent as character

## Missing variables

### Soil chemistry 2019 complement from PDF document
complement_2019 <- data.frame(
  plotID = c("OC11", "OC12", "OC13", "OS11", "OS12", "OS13", "OG11", "OG12", "OG13", "OS21", "OS22", "OS23", "OV21", "OV22", "OV23", "IC11", "IC12", "IC13", "IG11", "IG12", "IG13", "IG21", "IG22", "IG23", "IS11", "IS12", "IS13", "IS21", "IS22", "IS23", "IC21", "IC22", "IC23", "IV11", "IV12", "IV13", "UG11", "UG12", "UG13", "UG21", "UG22", "UG23", "US11", "US12", "US13", "US21", "US22", "US23", "US31", "US32", "US33", "US41", "US42", "US43"),
  dryMatterPercent = c(68.9, 66.7, 70.8, 96.8, 97.4, 97, 94.2, 94.9, 95.7, 95.6, 95.6, 95.5, 88.7, 93, 93.5, 96.2, 95.4, 95.8, 94.8, 95.4, 95.1, 93.6, 92.9, 94.4, 94.8, 94.4, 94.9, 94.7, 93.5, 93.7, 95.1, 95, 93.3, 96.7, 96.4, 97.1, 83.8, 89.1, 90.6, 94.8, 94.3, 90.6, 95.2, 94.5, 89, 94.6, 90.7, 94.3, 86.7, 93.5, 91.3, 92.9, 94.4, 93),
  totalNPercentDM = c(1.98, 2.08, 1.46, 0.28, 0.29, 0.34, 0.88, 0.65, 0.59, 0.54, 0.63, 0.62, 1.14, 0.61, 0.72, 0.5, 0.62, 0.5, 0.4, 0.45, 0.51, 0.76, 0.83, 0.74, 0.64, 0.62, 0.67, 0.59, 0.64, 0.92, 0.64, 0.7, 0.92, 0.52, 0.55, 0.71, 1.78, 1.53, 1.41, 0.71, 0.78, 1.26, 0.66, 0.63, 1.57, 0.72, 1.2, 0.94, 1.6, 0.62, 0.86, 1.08, 0.84, 0.67)
)
soil_chem_2019 <- full_join(soil_chem_2019, complement_2019)

### Soil chemistry 2020 complement from PDF document & soil_chem_2020_DM
#### totalNPercentDM
complement_2020 <- data.frame(
  plotID = c("OC21", "OC22", "OC23", "OC31", "OC32", "OC33", "OC41", "OC42", "OC43", "OC51", "OC52", "OC53", "OG21", "OG22", "OG23", "OG31", "OG32", "OG33", "OG51", "OG52", "OG53", "OG61", "OG62", "OG63", "IG31", "IG32", "IG33", "IS31", "IS32", "IS33", "IS41", "IS42", "IS43", "IS51", "IS52", "IS53", "OS31", "OS32", "OS33", "OS41", "OS42", "OS43", "OS51", "OS52", "OS53", "OS61", "OS62", "OS63", "OS71", "OS72", "OS73", "OS81", "OS82", "OS83", "OS91", "OS92", "OS93", "IC31", "IC32", "IC33", "IC41", "IC42", "IC43", "IC51", "IC52", "IC53", "OG41", "OG42", "OG43", "US51", "US52", "US53", "US61", "US62", "US63"),
  totalNPercentDM = c(0.25, 0.32, 0.33, 0.88, 1.13, 0.84, 0.54, 0.63, 0.63, 0.63, 0.76, 1.26, 0.69, 0.34, 0.47, 1.97, 1.4, 2.04, 1.57, 0.87, 0.42, 0.53, 0.5, 0.4, 0.27, 0.73, 0.5, 0.68, 0.95, 0.69, 0.41, 0.93, 0.42, 1.04, 1.21, 1.66, 0.69, 0.53, 0.5, 0.42, 0.3, 0.35, 0.61, 0.48, 0.42, 0.39, 0.4, 0.41, 1.21, 1.5, 1.82, 1.77, 1.96, 1.51, 0.76, 0.84, 0.63, 0.61, 0.77, 0.68, 0.31, 0.45, 0.39, 0.55, 1.33, 1.27, 0.56, 0.45, 0.27, 0.79, 0.85, 0.97, 0.65, 0.56, 1.22)
)
soil_chem_2020 <- full_join(soil_chem_2020, complement_2020)
#### dryMatterPercent
names(soil_chem_2020_DM) <- gsub("torrstoffG100G", "dryMatterPercent", names(soil_chem_2020_DM)) # Consistent variable name
soil_chem_2020_DM$dryMatterPercent <- as.numeric(soil_chem_2020_DM$dryMatterPercent) # Make variable as numeric
soil_chem_2020 <- left_join(soil_chem_2020, subset(soil_chem_2020_DM, select = c(plotID, dryMatterPercent)))

## Binding 2019 and 2020 datasets
soil_chem <- full_join(soil_chem_2019, soil_chem_2020)

## Remove empty columns
soil_chem <- soil_chem %>% 
  remove_empty()

## Consistent ID variables
soil_chem$plotID <- paste(substr(soil_chem$plotID, 1, 3), substr(soil_chem$plotID, 4, 4), sep= "-p")
soil_chem$siteID <- substr(soil_chem$plotID, 1, 3)

## Consistent & FAIR variable names
names(soil_chem) <- gsub("jordart", "soilType", names(soil_chem))
names(soil_chem) <- gsub("leirklasse", "clayContentClass", names(soil_chem))
names(soil_chem) <- gsub("glodetap", "lossOnIgnitionPercentDM", names(soil_chem))
names(soil_chem) <- gsub("volumvekt", "soilDensityKgL", names(soil_chem))
names(soil_chem) <- gsub("mold", "humusPercentDM", names(soil_chem))
names(soil_chem) <- gsub("humusPercentDMklasse", "humusContentClass", names(soil_chem))
names(soil_chem) <- gsub("pAl", "availablePMg100g", names(soil_chem))
names(soil_chem) <- gsub("kAl", "availableKMg100g", names(soil_chem))
names(soil_chem) <- gsub("mgAl", "availableMgMg100g", names(soil_chem))
names(soil_chem) <- gsub("caAl", "availableCaMg100g", names(soil_chem))
names(soil_chem) <- gsub("naAl", "availableNaMg100g", names(soil_chem))

## Selection desirable variables
soil_chem <- subset(
  soil_chem, select = -c(
    arstall,
    # anonymous
    journalnr,
    navn,
    adresse,
    postnr,
    poststed,
    registreringsdato
  )
)

## Duplicate check
# get_dupes(soil_chem) #validated

## Site distribution check
# table(soil_chem$siteID) # ov1 missing - in another PDF document

### Adding ov1 plots from PDF document
soil_chem <- soil_chem |> 
  add_row(plotID = "ov1-p1",
          soilType = 13,
          clayContentClass = 1,
          lossOnIgnitionPercentDM = 25.8,
          soilDensityKgL = 0.59,
          humusPercentDM = 25.8,
          humusContentClass = 5,
          pH = 5.4,
          availablePMg100g = 3,
          availableKMg100g = 14,
          availableMgMg100g = 20,
          availableCaMg100g = 31,
          availableNaMg100g = 6,
          dryMatterPercent = 0.5,
          totalNPercentDM = 0.44,
          siteID = "ov1"
  ) %>% 
  add_row(plotID = "ov1-p2",
          soilType = 13,
          clayContentClass = 1,
          lossOnIgnitionPercentDM = 26.0,
          soilDensityKgL = 0.74,
          humusPercentDM = 26.0,
          humusContentClass = 5,
          pH = 4.9,
          availablePMg100g = 3,
          availableKMg100g = 13,
          availableMgMg100g = 15,
          availableCaMg100g = 21,
          availableNaMg100g = 5,
          dryMatterPercent = 91.4,
          totalNPercentDM = 0.51,
          siteID = "ov1"
  ) %>% 
  add_row(plotID = "ov1-p3",
          soilType = 14,
          clayContentClass = 1,
          lossOnIgnitionPercentDM = 43.2,
          soilDensityKgL = 0.38,
          humusPercentDM = 43.2,
          humusContentClass = 6,
          pH = 5.0,
          availablePMg100g = 3,
          availableKMg100g = 12,
          availableMgMg100g = 14,
          availableCaMg100g = 16,
          availableNaMg100g = 5,
          dryMatterPercent = 85.4,
          totalNPercentDM = 0.72,
          siteID = "ov1"
  )
  
## Character variables - desirable categories, NAs, misprints

### Consistent lower camel
soil_chem <- soil_chem %>% 
  mutate_if(is.character, tolower)

### Explicit category names for class variables from Eurofins classification
soil_chem <- soil_chem %>% 
  mutate(
    soilType = ifelse(
      soilType == 2, "medium-sand", ifelse(
        soilType == 3, "fine-sand", ifelse(
          soilType == 5, "silty-medium-sand", ifelse(
            soilType == 6, "silty-fine-sand", ifelse(
              soilType == 13, "mineral-mixed-humus-soil", "organic-soil"
            )
          )
        )
      )
    )
  ) %>% 
  mutate(
    clayContentClass = ifelse(clayContentClass == 1, "0-5%", "5-10%")
  ) %>% 
  mutate(
    humusContentClass = ifelse(
      humusContentClass == 5, "mineral-mixed-humus", ifelse(
        humusContentClass == 6, "organic-soil", "moderately-humus-rich"
      )
    )
  )

### Categories & distribution
# table(soil_chem$soilType) #validated
# table(soil_chem$clayContentClass) #validated
# table(soil_chem$humusContentClass) #validated

## Numeric variables - min/max, distribution, potential outliers

### Min/max
# test <- soil_chem |>
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
#   transpose() # Check max Ca over 1000 mg/100g

### NA check
# colnames(soil_chem)[apply(soil_chem, 2, anyNA)] #validated

### Variable distribution & outliers
hist(soil_chem$lossOnIgnitionPercentDM) # 

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


