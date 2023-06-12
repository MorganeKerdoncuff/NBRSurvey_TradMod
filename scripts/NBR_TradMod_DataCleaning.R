# TradMod WP2 - NBR survey cleaning script
#Description of the data
#Date
#Who
#Project
#Funding
#Place

#### PACKAGE LOADING ####

library(tidyverse) #R language
library(readxl) #read xl files
library(lubridate) #standard date data
#library(raster) # work with raster files -> soon outdated, replaced by terra
#library(rgdal) # work with raster files -> soon outdated
library(terra) # work with raster files
library(purrr) # merging several files at the same time

#### DATA LOADING ####

siteinfo_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SiteInfo") # Farm information, field location and habitat type
landuse_raw <- read_excel(path = "data/rawdata/NBR_RawFarmerSurvey.xlsx", sheet="Farms Information_R") # Field management data from farmer interview, at site level
#farmer_raw <- read_excel(path = "data/rawdata/NBR_RawFarmRepertoire.xlsx") # Farm management data from farmer interview med Margit, at site level
landscape_raw <- read_excel(path = "data/rawdata/NBR_RawLandscapeMatrix.xlsx", sheet="MatrixProportion") #Land cover data around the fields from Geonorge, at site level
area20x20_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="20mX20m") # Sampling area description, vegetation cover at the site level
groundcover_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SoilCover") # Vegetation cover at the quadrat (subplot) level
soilpene_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="SoilPenetration") # Penetration rate in the soil, at subplot level (two collection per subplot)
soilbulk_raw <- read_excel(path = "data/rawdata/NBR_RawBD.xlsx", na="NA") # Soil bulk density, at subplot level (three samples per subplot)
soilmeso_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="Mesofauna") # Height of mesofauna soil core, at subplot level
chem2019 <- read.csv("data/rawdata/NBR_RawSoilChemistry2019.txt", sep=";") # 2019 Soil chemistry data from Eurofins, at plot level
chem2020 <- read.csv("data/rawdata/NBR_RawSoilChemistry2020.txt", sep=";") # 2020 soil chemistry data from Eurofins, at plot level
chem2020_DM <- read_excel(path = "data/rawdata/NBR_RawSoilChemistry2020bis.xls") # 2020 complementary soil chemistry data from Eurofins, at plot level
poo_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="Poo", na="NA") # poo data at subplot (quadrat) level
vege_raw <- read_excel(path = "data/rawdata/NBR_RawAll.xlsx", sheet="PlantRichness") # plant community data, at species level
arthro_main <- read_excel(path = "data/rawdata/NBR_RawArthro.xlsx", na="NA") # arthropod community data, at family level for beetles and order level for other arthropods
arthro_sup <- read_excel(path = "data/rawdata/NBR_RawArthroSup.xlsx", na="NA") # complementary arthropod community data, at family level for beetles and order level for other arthropods
mesobio_raw <- read_excel(path = "data/rawdata/NBR_RawMesobio.xlsx", na="NA")


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
# [9] Current grazing livestock, at least for the last 5 years
# [10] Size of the livestock flock - column empty, data collected in another file
# [11] Size of the grazing area - column empty, data collected in another file
# [12] Y/N if the animal was seen on site during the collection
# [13] Poo weight collected on 10% of the sampling area (g) - proxy of grazing intensity on site
# [14] Comments

#
## Summary site info - Check table size, list of variables, variable types (num/chr)

#str(siteinfo_raw) # All good

#
## Name & character cleaning

# R friendly variable names
names(siteinfo_raw) <- gsub("%", "percent", names(siteinfo_raw))
names(siteinfo_raw) <- gsub(" ", "", names(siteinfo_raw)) 
names(siteinfo_raw) <-  gsub("\\(", "_", names(siteinfo_raw))
names(siteinfo_raw) <-  gsub("\\)", "", names(siteinfo_raw))
names(siteinfo_raw) <- gsub("Comments", "Comments_siteinfo", names(siteinfo_raw))

# Removal empty columns
#siteinfo_raw <- subset(siteinfo_raw, select = -c(Size_livestock, Surface)) # remove empty variables, which will be included in another dataset

#
## Data cleaning - new R object with removal empty columns or variables redundant with other datasets

siteinfo_full <- subset(siteinfo_raw, select = -c(Size_livestock, Surface))

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- siteinfo_full |>  
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
#hist(siteinfo_full$`Pooestimation_g-10percent`) # Standard distribution of poo data - validated

#
## Char var Site Info - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(siteinfo_full$SiteID) # Unique ID for each site - validated

# Ecological zones
#unique(siteinfo_full$NiBioEcologicalZone) # Three categories for ecological zone - validated
# Homogeneous character writing - removing capital letters for common nouns
siteinfo_full <- siteinfo_full |> 
  mutate(NiBioEcologicalZone = dplyr::recode(NiBioEcologicalZone, "Coastal" = "coastal")) |> 
  mutate(NiBioEcologicalZone = dplyr::recode(NiBioEcologicalZone, "Fjord" = "fjord")) |> 
  mutate(NiBioEcologicalZone = dplyr::recode(NiBioEcologicalZone, "Mountain" = "mountain"))

# Geographical locations
#unique(siteinfo_full$Location) # Correct locations and location names

# Muncipalities
#unique(siteinfo_full$Municipality) # missing data + US2 Stordalen -> did you guys went that far ? Yes validated
#siteinfo_full[is.na(siteinfo_full$Municipality),] # NA identified - US2, Stordalen is located in Masfjorden
siteinfo_full <- siteinfo_full |> 
  mutate(Municipality = ifelse(SiteID == "US2", "Masfjorden", Municipality)) # Replace NA by municipality name
#unique(siteinfo_full$Municipality) # NA in Municipality replaced - validated

# Livestock
#unique(siteinfo_full$Type_livestock) # Villsau should be within sheep category
siteinfo_full <- siteinfo_full |> 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "Villsau" = "sheep")) # Only 3 livestock categories - validated
# Homogeneous character writing - removing capital letters for common nouns
siteinfo_full <- siteinfo_full |> 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "Sheep" = "sheep")) |> 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "Cows" = "cow")) |> 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "Goats" = "goat"))
#table(siteinfo_full$Type_livestock) # correct number of sites per livestock category - validated

# Habitat type
#unique(siteinfo_full$Habitat) # 3 categories + one extra (bog) which will not be considered in the analysis
# Homogeneous character writing - removing capital letters for common nouns
siteinfo_full <- siteinfo_full |> 
  mutate(Habitat = dplyr::recode(Habitat, "Coastal heathland" = "coastal heathland")) |> 
  mutate(Habitat = dplyr::recode(Habitat, "Permanent grassland" = "permanent grassland")) |> 
  mutate(Habitat = dplyr::recode(Habitat, "Subalpine heathland" = "subalpine heathland"))
#table(siteinfo_full$Habitat) # habitat repartition - validated

#unique(siteinfo_full$Animal_on_site)


## Export clean data in new excel file

write_csv(siteinfo_full, "data/cleandata/NBR_FullSiteInfo.csv")


#### FIELD MANAGEMENT ####

## Description

## List of variables

# [1] Field identification code for data collection
# [2] Geographical location of the field, hamlet
# [3] Geographical location of the field, postcode
# [4] Geographical location of the field, municipality (both former and current classification)
# [5] Date of the interview with the farmer
# [6] Site productivity - infield high productivity, outfield low productivity
# [7] Main livestock, currently grazing in the field - ! rotational grazing management
# [8] Other livestock, grazing in other fields - ! rotational grazing management
# [9] Name of cow breed
# [10] Name of goat breed
# [11] Name of sheep breed
# [12] Number of adult animals in the main livestock
# [13] Number of young animals in the main livestock
# [14] Number of adult animals in other livestock
# [15] Number of young animals in other livestock
# [16] Size of the field containing the sampling area (ha) ! rotational grazing management
# [17] FROM GÅRDSKSART - Total grazing area of the farm (ha) ! rotational grazing management
# [18] Grazing density on the field - not collected, empty column
# [19] Period since the farmer has had the current livestock
# [20] Period without grazing management on the field
# [21] Number of months the main livestock grazes in the field in a year
# [22] Number of months other livestock graze in the field in a year ! rotational grazing
# [23] If the animals are kept inside or outside at night
# [24] Farmer impression of grazing pressure on the site
# [25] Former livestock which used to graze in the field (1)
# [26] Former livestock which used to graze in the field (2)
# [27] Former livestock which used to graze in the field (3)
# [28] Former livestock grazing period (1)
# [29] Former livestock grazing period (2)
# [30] Former livestock grazing period (3)
# [31] Type of farm management during the past 10 years (1)
# [32] Type of farm management during the past 10 years (2)
# [33] Type of farm management during the past 10 years (3)
# [34] If applicable, frequency of grass cutting
# [35] If applicable, frequency of mulching
# [36] If applicable, frequency of fertilization with manure
# [37] If applicable, type of manure used
# [38] If applicable, volume of manure used (m3)
# [39] If applicable, in which season the manure is used
# [40] If applicable, frequency of fertilization with artificial fertilizer
# [41] If applicable, weight of artificial fertilizer used (kg)
# [42] If applicable, in which season the artificial fertilizer is used
# [43] If applicable, frequency of fertilization with shell-sand/lime
# [44] If applicable, weight of shell-sand/limer used (kg)
# [45] If applicable, in which season the shell-sand/lime is used
# [46] If applicable, last time the field was sowed
# [47] If applicable, last time the field was plowed
# [48] If applicable, last time the field was drained
# [49] Type of land use prior to grazing
# [50] If applicable, type of production prior to grazing (1)
# [51] If applicable, type of production prior to grazing (2)
# [52] Time period of land use type prior to grazing
# [53] Comments

#
## Summary land use - Check table size, list of variables, variable types (num/chr)

#str(landuse_raw) # All good

#
## Name & character cleaning land use

# R friendly variable names
names(landuse_raw)<-  gsub(" ", "", names(landuse_raw))
names(landuse_raw)<-  gsub("\\)", "", names(landuse_raw))
names(landuse_raw)<-  gsub("\\(", "_", names(landuse_raw))
names(landuse_raw)<-  gsub("/", "_", names(landuse_raw))
names(landuse_raw)<-  gsub(";", "_", names(landuse_raw))
names(landuse_raw)<-  gsub("-", "_", names(landuse_raw))
names(landuse_raw)<-  gsub("\\?", "", names(landuse_raw))
names(landuse_raw)<-  gsub(":", "", names(landuse_raw))
names(landuse_raw) <- gsub("Sitecode", "SiteID", names(landuse_raw)) #rename siteID so it matches with other sheets
names(landuse_raw) <- gsub("Typeoflivestock", "Livestock", names(landuse_raw))
names(landuse_raw) <- gsub("_ifapplicable", "", names(landuse_raw))
names(landuse_raw) <- gsub("_villsaubreed", "breed", names(landuse_raw))
names(landuse_raw) <- gsub("Numberofanimals", "FlockSize", names(landuse_raw))
names(landuse_raw) <- gsub("Surveysitegrazingsurface_ha", "GrazingSurface_ha", names(landuse_raw))
names(landuse_raw) <- gsub("Inmarksbeite_Gardskart", "TotalInfieldSurface", names(landuse_raw))
names(landuse_raw) <- gsub("Currentlivestock_fromdate_year", "LivestockFrom_year", names(landuse_raw))
names(landuse_raw) <- gsub("Gap_swithoutgrazing_timeperiod", "NoGrazing_period", names(landuse_raw))
names(landuse_raw) <- gsub("Yearlygrazingtimesurvey_year1_currentlivestock_monthsperyear", "YearlyGrazing1_month", names(landuse_raw))
names(landuse_raw) <- gsub("Yearlygrazingtimesurvey_year2_currentlivestock_monthsperyear", "YearlyGrazing2_month", names(landuse_raw))
names(landuse_raw) <- gsub("Impressionofgrazingpressureonsurveysite", "FarmerImpression_GrazingPressure", names(landuse_raw))
names(landuse_raw) <- gsub("_type", "", names(landuse_raw))
names(landuse_raw) <- gsub("_Ifseveral", "", names(landuse_raw))
names(landuse_raw) <- gsub("Type_soffarmmanagement_soiltreatmentusedthelast10years", "FieldManagement1", names(landuse_raw))
names(landuse_raw) <- gsub("1...32", "2", names(landuse_raw))
names(landuse_raw) <- gsub("1...33", "3", names(landuse_raw))
names(landuse_raw) <- gsub("frequencyoffertilizing", "_freq", names(landuse_raw))
names(landuse_raw) <- gsub("frequency", "_freq", names(landuse_raw))
names(landuse_raw) <- gsub("type", "_type", names(landuse_raw))
names(landuse_raw) <- gsub("amount_vol_m3", "_volm3", names(landuse_raw))
names(landuse_raw) <- gsub("Artificialfertilizer", "ArtificialFert", names(landuse_raw))
names(landuse_raw) <- gsub("Artificialfertiliser", "ArtificialFert", names(landuse_raw))
names(landuse_raw) <- gsub("Art.fertilizer", "ArtificialFert", names(landuse_raw))
names(landuse_raw) <- gsub("amount_mass_kg", "_masskg", names(landuse_raw))
names(landuse_raw) <- gsub("season", "_season", names(landuse_raw))
names(landuse_raw) <- gsub("Shell_sand_lime", "ShellSandLime", names(landuse_raw))
names(landuse_raw) <- gsub("...53", "", names(landuse_raw))
names(landuse_raw) <- gsub("...54", "2", names(landuse_raw))
names(landuse_raw) <- gsub("Comments", "Comments_landuse", names(landuse_raw))

# Removal empty columns
#landuse_raw <- subset(landuse_raw, select = -c(Grazingdensity_perha)) # will be calculated later in the script

#
## Data cleaning - new R object

landuse_full <- landuse_raw

#
## Char var land use - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(landuse_full$SiteID) # Unique ID for each site - validated

# Livestock1
#unique(landuse_full$Livestock1) # villsau should be in sheep category
landuse_full <- landuse_full |> 
  mutate(Livestock1 = dplyr::recode(Livestock1, "villsau" = "sheep"))
#table(landuse_full$Livestock1) # correct number of sites per livestock category - validated

# Livestock2
#unique(landuse_full$Livestock2) # villsau should be in sheep category
landuse_full <- landuse_full |> 
  mutate(Livestock2 = dplyr::recode(Livestock2, "villsau" = "sheep"))
#table(landuse_full$Livestock2) # Only 5 farmers with another livestock -> not to be included in the analysis

# Period without grazing
#unique(landuse_full$NoGrazing_period) # NAs ?
#landuse_full[is.na(landuse_full$NoGrazing_period),] # 8 NAs, including the 6 farmers who did not respond to the survey
#table(landuse_full$NoGrazing_period) # Only 7 sites with potential grazing interruption, in a quite hectic way (some 2-3 years, other 10 years), 30 sites with no interruption -> not to be included into the analysis, but interesting for the material and methods
#landuse_full[landuse_full$NoGrazing_period != "no",] # Interruption are on OS1, IG1, IS2, OS5, OS7, OC2 and OS3

# If the animals are kept inside or outside at night
#unique(landuse_full$Inside_outsideatnight) # some NAs
#landuse_full[is.na(landuse_full$Inside_outsideatnight),] #18 sites with missing values -> not to be included in the analysis

# Farmer impression of grazing pressure on the site
#unique(landuse_full$Impressionofgrazingpressureonsurveysite) # NAs
#landuse_full[is.na(landuse_full$Impressionofgrazingpressureonsurveysite),] # 31 missing values -> not to be included in the analysis

# Former livestock (1), before the current livestock
#unique(landuse_full$Formerlivestock1) # NAs + villsau should be sheep category
landuse_full <- landuse_full |> 
  mutate(Formerlivestock1 = dplyr::recode(Formerlivestock1, "villsau" = "sheep"))
#table(landuse_full$Formerlivestock1) # 11 former cows, 13 former sheep, 2 former horse and 1 "other"
#landuse_full[is.na(landuse_full$Formerlivestock1),] # 18 missing values, including the 6 farmers which did not reply the survey -> 12 with no former livestock known -> could be used in the analysis?

# Former livestock (2), before livestock (1)
#unique(landuse_full$Formerlivestock2) # NAs + villsau should be sheep category
landuse_full <- landuse_full |> 
  mutate(Formerlivestock2 = dplyr::recode(Formerlivestock2, "villsau" = "sheep"))
#table(landuse_full$Formerlivestock2) # 5 former cows, 2 former sheep, 2 former horse
#landuse_full[is.na(landuse_full$Formerlivestock2),] # 36 missing values -> not to be included in the analysis

# Former livestock (3), before livestock (2)
#unique(landuse_full$Formerlivestock3) # NAs
#table(landuse_full$Formerlivestock3) # 1 former cows, 1 other
#landuse_full[is.na(landuse_full$Formerlivestock3),] # 43 missing values -> not to be included in the analysis

# Former livestock (1) grazing period
#unique(landuse_full$Formerlivestock1_timeperiod) # quite hectic, need to be transformed into how many years the former livestock has grazed -> not to be included in the analysis

# Management or soil treatment used in the past 10 years
unique(landuse_full$FieldManagement1) # NAs
unique(landuse_full$FieldManagement2)
unique(landuse_full$FieldManagement3)
#landuse_full[is.na(landuse_full$FieldManagement1),] # 6 missing values -> farmers who did not reply to the survey
#landuse_full[is.na(landuse_full$FieldManagement2),] # 36 missing values
#landuse_full[is.na(landuse_full$FieldManagement3),] # 41 missing values
#table(landuse_full$FieldManagement1) # 9 sites with no treatments at all, 21 with at least fertilization
#table(landuse_full$FieldManagement2) # 8 sites with at least 2 treatments
#table(landuse_full$FieldManagement3) # 3 sites with 3 treatments

# Sites with no treatment
#subset(landuse_full, FieldManagement1 == "none") # 4 mountain sites (US2, US3, US5, US6), 2 coastal heathlands (IS2, OS5), 3 grasslands (IG1, IG2, IS3)

# Sites with burning
#subset(landuse_full, FieldManagement1 == "burning") # 3 coastal heathlands only treatment (OV1, OV2, OS7) -> validated
#subset(landuse_full, FieldManagement3 == "burning") # 1 grassland (IS1), with also mulching and fertilisation -> validated

# Sites with tree cutting
#subset(landuse_full, FieldManagement1 == "sitka spruce removal") # 1 coastal heathland only treatment (OS9) -> validated
#subset(landuse_full, FieldManagement1 == "tree cutting") # 1 mountain site only treatment (UG1), 1 mountain site with also grass cutting (UG2) -> validated

# Sites with herbicide
#subset(landuse_full, FieldManagement1 == "herbicides") # 1 grassland (OS5) with also tree cutting as second treatment
#subset(landuse_full, FieldManagement3 == "herbicides") # 1 grassland (OC2)

# Sites with mulching
#subset(landuse_full, FieldManagement1 == "mulching") # 2 grassland (OS3, OS1) only treatment -> validated
#subset(landuse_full, FieldManagement2 == "mulching") # 6 grasslands (OC1, IS1, IV1, IC5, OS4, OC2) with mulching + fertilization. OC2 also has herbicides.

# Frequency of mulching
#table(landuse_full$Mulching_freq) # 8 sites with mulching
#subset(landuse_full, Mulching_freq == "every year") # 3 sites with annual mulching (IS1, OS4, OS3) -> validated
#subset(landuse_full, Mulching_freq == "sometimes") # 5 sites with occasional mulching (OC1, OS1, IV1, IC5, OC2) -> validated

# Types of fertilization
#unique(landuse_full$Fertilizing_type1)
#unique(landuse_full$Fertilizing_type2)
#unique(landuse_full$Fertilizing_type3) # 3 types of fertilization: artificial fertilizer, organic (manure) and inorganic (shell-sand/lime)
#table(landuse_full$Fertilizing_type3) # 5 sites with 3 types of fertilization
#subset(landuse_full, Fertilizing_type3 != is.na(Fertilizing_type3)) # 5 grassland sites (IC1, OC3, IC4, IC5, OC2)
#table(landuse_full$Fertilizing_type2) # 11-5 = 6 sites with 2 types of fertilization
#filter(landuse_full, !is.na(Fertilizing_type2) & is.na(Fertilizing_type3)) # 6 grasslands (OC1, IS1, OC4, OG6, OC5, OS8)
#table(landuse_full$Fertilizing_type1) # 21-11 = 10 sites with only 1 type of fertilization
#filter(landuse_full, !is.na(Fertilizing_type1) & is.na(Fertilizing_type2)) # 10 grasslands (OG1, OS2, IV1, OG3, IC3, OG2, OG4, OS4, OS6, IS4)

# Frequency of manure
#unique(landuse_full$Manure_freq) # possibilities are "yearly" or "occasionally"
#subset(landuse_full, Manure_freq == "every year") # 9 grasslands (OC1, IC1, IS1, OC3, OC4, OG6, IC4, IC5, OC5) with annual manure
#subset(landuse_full, Manure_freq == "sometimes") # 2 grasslands (OC2, OS8) with occasional manure

# Frequency of artificial fertilization
#unique(landuse_full$ArtificialFert_freq) # possibilities are "yearly" or "occasionally"
#subset(landuse_full, ArtificialFert_freq == "every year") # 15 grasslands (OC1, OS2, IC1, OG3, IC3, OC3, OC4, OG2, OG4, IC4, IC5, OC5, OS6, OC2, IS4) with annual fertilization with artificial fertiliser
#subset(landuse_full, ArtificialFert_freq == "sometimes") # 4 grasslands (OC2, OG6, IV1, OS8) with occasional fertilization with artificial fertiliser

# Frequency of inorganic fertilization with shell-sand/lime
#unique(landuse_full$ShellSandLime_freq) # possibilities are "occasionally" only
#subset(landuse_full, ShellSandLime_freq == "sometimes") # 7 grassland (IC1, OG3, IC3, OC3, IC4, IC5, OC2)

#
## New variables for field management/treatments for cleaner and standardised output

# Herbicide
landuse_full <- landuse_full |> 
  mutate(FM_herbicide = ifelse(
      FieldManagement1 == "herbicides" |
        FieldManagement3 == "herbicides",
      "yes",
      "no"))

# Tree cutting
landuse_full <- landuse_full |> 
  mutate(FM_treecutting = ifelse(
      FieldManagement1 == "tree cutting" |
        FieldManagement1 == "sitka spruce removal" |
        FieldManagement2 == "tree cutting" |
        FieldManagement3 == "tree cutting",
      "yes",
      "no"))

# Burning
landuse_full <- landuse_full |> 
  mutate(FM_burning = ifelse(
      FieldManagement1 == "burning" |
        FieldManagement3 == "burning",
      "yes",
      "no"))

# Mowing
landuse_full <- landuse_full |> 
  mutate(FM_mowing = ifelse(FieldManagement2 == "grass cutting", "yes", "no"))

# Mulching
landuse_full <- landuse_full |> 
  mutate(FM_mulching = ifelse(
      Mulching_freq == "every year", "yearly", ifelse(
        Mulching_freq == "sometimes", "occasionally", "no"
      )))

# Manure
landuse_full <- landuse_full |> 
  mutate(FM_manure = ifelse(
    Manure_freq == "every year", "yearly", ifelse(
      Manure_freq == "sometimes", "occasionally", "no"
  )))

# Artificial fertilizer
landuse_full <- landuse_full |> 
  mutate(FM_artificialfert = ifelse(
    ArtificialFert_freq == "every year", "yearly", ifelse(
      ArtificialFert_freq == "sometimes", "occasionally", "no"
    )))

# Shell-sand-lime
landuse_full <- landuse_full |> 
  mutate(FM_shellsandlime = ifelse(
    ShellSandLime_freq == "every year", "yearly", ifelse(
      ShellSandLime_freq == "sometimes", "occasionally", "no"
    )))

# Replace NAs by "no" in summary management columns only
landuse_full <- landuse_full |> 
  mutate_at(vars(starts_with("FM_")) & filter(
                   SiteID != "IC2" |
                SiteID != "US1" |
                SiteID != "US2" |
                SiteID != "US4" |
                SiteID != "IS5" |
                SiteID != "IG3"),
            ~replace(is.na(.), "no"))

landuse_full <- landuse_full |> 
  mutate_at(vars(starts_with("FM_")), ~replace(., is.na(.), "no"))

#
## Numeric var land use - Check min/max, distribution and potential outliers

# Check min/max
test <- landuse_full |>  
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

# Flock size adults - missing values
#landuse_full[is.na(landuse_full$FlockSize1_adults),] # some farmers would not or could not respond to our interview. The missing value were therefore taken from another series of interviews made in 2018
landuse_full <- landuse_full |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "IC2", 37, FlockSize1_adults)) |> 
  mutate(Cowbreed = ifelse(SiteID == "IC2", "norsk rødt fe", Cowbreed)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "UC1", 13, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "US1", 4, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "US4", 66, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "IS5", 90, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "OS7", 39, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "IG3", 109, FlockSize1_adults)) # no information for US6

# Flock size adults - distribution
#hist(landuse_full$FlockSize1_adults) # wide range but most farms under 50 animals - one farm over 150 animals
#landuse_full[landuse_full$FlockSize1_adults>150,] # IS4 over 150 animals + US6 as NA

# Flock size young - missing values
landuse_full <- landuse_full |> 
  mutate(FlockSize1_young = ifelse(SiteID == "IC2", 17, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "UC1", 13, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "US1", 10, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "US4", 63, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "IS5", 181, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "OS7", 14, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "IG3", 135, FlockSize1_young)) # no information for US6

# Flock size young - distribution
#hist(landuse_full$FlockSize1_young)
#landuse_full[landuse_full$FlockSize1_young>150,] # US3 & IS5 + 4 farms with no young animals + US6 as NA

# Grazing surface
#landuse_full[is.na(landuse_full$GrazingSurface_ha),] # No missing value - validated
#hist(landuse_full$GrazingSurface_ha) # big range, one outlier over 1000 ha
#landuse_full[landuse_full$GrazingSurface_ha>1000,] # UG2 in the mountain -> extended area which serves for several farmers

# Total infield surface - missing values & distribution
#landuse_full[is.na(landuse_full$TotalInfieldSurface),] # All missing values should be outfields/heathland sites - validated
hist(landuse_full$TotalInfieldSurface) # big range but no visible outlier

# Date since the current livestock has been grazing
#landuse_full[is.na(landuse_full$LivestockFrom_year),] # 7 missing values from farmers who did not reply the survey
#hist(landuse_full$LivestockFrom_year) # ranges from 1920 to 2020, but ! not equivalent to no grazing, neither to grazing with another animal -> too hectic, not to be included in the analysis

# How many months main livestock graze on site during the year
#landuse_full[is.na(landuse_full$YearlyGrazing1_month),] # 6 missing values from farmers who did not reply the survey
#hist(landuse_full$YearlyGrazing1_month) # coherent values (between 2 and 12), no visible outliers -> validated

#
## New variable - grazing intensity

# Livestock unit for adults - depends on type of livestock and breed (for cows)
#unique(landuse_full$Cowbreed) # Milking cows (norsk rødt fe) are 1 SSU - Beef cattle (Limousin, Aberdeen angus, Highland) are 0.8 LSU - Sheep and goats are 0.1 LSu
landuse_full <- landuse_full |> 
  mutate(LSU_adultind = ifelse(Livestock1 == "sheep" | Livestock1 == "goat", 0.1,
                               ifelse(Cowbreed == "norsk rødt fe", 1, 0.8)))
#hist(landuse_full$LSU_adultind) # validated

# Livestock unit for young - half for goats and sheep (0.05), young cows on field (0.7)
landuse_full <- landuse_full |> 
  mutate(LSU_youngind = ifelse(Livestock1 == "sheep" | Livestock1 == "goat", 0.05, 0.7))
#hist(landuse_full$LSU_youngind) # validated

# Grazing intensity over the year - two goat flocks in Osterøy were known to be moved to summer farm (outfield)
landuse_full <- landuse_full |>  
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
  mutate(Grazingdensity_perha = ifelse(
    Livestock1 == "goat" & Municipality_old == "Osterøy", 
    ((FlockSize1_adults*LSU_adultind+FlockSize1_young*LSU_youngind)/TotalInfieldSurface)*(YearlyGrazing1_month/12), 
    (FlockSize1_adults*LSU_adultind+FlockSize1_young*LSU_youngind)/TotalInfieldSurface)
    )
#hist(landuse_full$Grazingdensity_perha) # Range from 0 to 0.6, no visible outlier
#landuse_full[is.na(landuse_full$Grazingdensity_perha),] # US6 as NA


## Export clean data in new excel file

write_csv(landuse_full, "data/cleandata/NBR_FullLanduse.csv")



#### LANDSCAPE MATRIX ####

## Description

## List of variables

# [1] Field identification code for data collection
# [2] Area cover of fully cultivated land within 1 km buffer zone around the field
# [3] Area cover of partially cultivated land within 1 km buffer zone around the field
# [4] Area cover of grazed infields within 1 km buffer zone around the field
# [5] Area cover of productive forest within 1 km buffer zone around the field
# [6] Area cover of unproductive forest within 1 km buffer zone around the field
# [7] Area cover of wetland within 1 km buffer zone around the field
# [8] Area cover of outfield within 1 km buffer zone around the field
# [9] Area cover of freshwater within 1 km buffer zone around the field
# [10] Area cover of infrastructure within 1 km buffer zone around the field
# [11] Area cover of sea within 1 km buffer zone around the field

#
## Summary landscape - Check table size, list of variables, variable types (num/chr)

#str(landscape_raw) # All good

#
## Data cleaning - new R object

landscape_full <- landscape_raw

#
## Char var land use - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(landscape_full$SiteID) # Unique ID for each site & only grassland sites

#
## Numeric var - Check min/max, distribution and potential outliers

# Check min/max
test <- landscape_full |>  
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
  transpose() # no visible outlier, no data over 100%

# Fully cultivated land
#landscape_full[is.na(landscape_full$FullyCultivatedLand_percent),] # No NA
#hist(landscape_full$FullyCultivatedLand_percent) # Range between 0 and 20%, no outlier
#filter(landscape_full, FullyCultivatedLand_percent == 0) # No field without cultivated land around

# Partially cultivated land
#landscape_full[is.na(landscape_full$SuperficiallyCultivatedLand_percent),] # No NA
#hist(landscape_full$SuperficiallyCultivatedLand_percent) # Range between 0 and 20%, no outlier
#filter(landscape_full, SuperficiallyCultivatedLand_percent == 0) # No field without Partially land around

# Grazed infield
#landscape_full[is.na(landscape_full$Infield_percent),] # No NA
#hist(landscape_full$Infield_percent) # Range between 0 and 35%, no outlier
#filter(landscape_full, Infield_percent == 0) # No field without infield around

# Productive forest
#landscape_full[is.na(landscape_full$ProductiveForest_percent),] # No NA
#hist(landscape_full$ProductiveForest_percent) # Range between 0 and 70%, no outlier
#filter(landscape_full, ProductiveForest_percent == 0) # No field without productive forest around

# Unproductive forest
#landscape_full[is.na(landscape_full$NonProductiveForest_percent),] # No NA
#hist(landscape_full$NonProductiveForest_percent) # Range between 0 and 25%, no outlier
#filter(landscape_full, ProductiveForest_percent == 0) # No field without unproductive forest around

# Wetlands
#landscape_full[is.na(landscape_full$Wetland_percent),] # No NA
#hist(landscape_full$Wetland_percent) # Range between 0 and 20%, no outlier
#filter(landscape_full, Wetland_percent == 0) # 1 site without wetland around (IG3)

# Outfield
#landscape_full[is.na(landscape_full$Outfield_percent),] # No NA
#hist(landscape_full$Outfield_percent) # Range between 0 and 45%, no outlier
#filter(landscape_full, Outfield_percent == 0) # No field without outfield around

# Freshwater
#landscape_full[is.na(landscape_full$Freshwater_percent),] # No NA
#hist(landscape_full$Freshwater_percent) # Range between 0 and 16%, no outlier
#filter(landscape_full, Freshwater_percent == 0) # No field without freshwater around

# Infrastructure
#landscape_full[is.na(landscape_full$Infrastructure_percent),] # No NA
#hist(landscape_full$Infrastructure_percent) # Range between 0 and 25%, no outlier
#filter(landscape_full, Infrastructure_percent == 0) # No field without infrastructure around

# Sea
#landscape_full[is.na(landscape_full$Sea_percent),] # No NA
#hist(landscape_full$Sea_percent) # Range between 0 and 70%, no outlier
#filter(landscape_full, Sea_percent == 0) # 11 fields without sea (ocean or fjord) around


## Export clean data in new excel file

write_csv(landscape_full, "data/cleandata/NBR_FullLandscape.csv")



#### SAMPLING AREA 20X20 ####

## Description

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
# [15] Distance to sea, either ocean or fjord (m)
# [16] Estimated percentage cover of exposed rock in the sampling area
# [17] Estimated percentage cover of mud in the sampling area
# [18] Estimated percentage cover of trees and shrubs over 1 m in the sampling area
# [19] Estimated percentage cover of shrubs under 1 m in the sampling area
# [20] Estimated percentage cover of herbs in the sampling area
# [21] Estimated percentage cover of monocotyledons (grasses, rushes, sedges) in the sampling area
# [22] Estimated percentage cover of bryophytes (mosses, liverworts) in the sampling area
# [23] Estimated percentage cover of ferns in the sampling area
# [24] Estimated percentage cover of lichens in the sampling area
# [25] Comments

#
## Summary sampling area - Check table size, list of variables, variable types (num/chr)

#str(area20x20_raw) # All good, date should be reformatted

#
## Name & character cleaning sampling area

# R friendly variable names
names(area20x20_raw)<-gsub("%", "percent", names(area20x20_raw)) # remove percent signs from names
names(area20x20_raw)<-gsub("&", "_", names(area20x20_raw)) # remove &
names(area20x20_raw)<-gsub("\\(", "", names(area20x20_raw)) # remove (
names(area20x20_raw)<-gsub("\\)", "", names(area20x20_raw)) # remove )
names(area20x20_raw) <- gsub("Comments", "Comments_area20x20", names(area20x20_raw))

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
#hist(area20x20_full$Elevation_max) # Sites range from 0 to 900 m elevation, no outlier - validated

# Distribution slope
#hist(area20x20_full$General_slope) # Slopes range from 0 to 35 degrees, no outlier - validated

# Distribution aspect degree
#hist(area20x20_full$AspectDegree) # Aspect covers all spectrum (0 to 360), no outliers - validated

# Distribution distance to sea
#area20x20_full[is.na(area20x20_full$DistanceToSea_m),] # 9 NA, corresponding to upland sites -> validated
#hist(area20x20_full$DistanceToSea_m) # One outlier, above 10 km distance
#filter(area20x20_full, DistanceToSea_m>10000) # IG3 in Modalen

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

#
## Calculation new variables

# Heat Load Index
area20x20_full <- area20x20_full |> 
  mutate(HLI = cos(AspectDegree-225)*tan(General_slope))
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
# [17] Bulk density calculated from the non-corrected volumes -> not to be used in the analysis
# [18] Bulk density calculated for the corrected volumes only -> not to be used in the analysis
# [19] Bulk density calculated from the best volume estimation
# [20] Weight of percentage of soil moisture (g)
# [21] Volume of percentage of soil moisture (cm3) -> calculated from the BD
# [22] Percentage of soil porosity -> calculated from the BD
# [23] Percentage of WFPS
# [24] Comments during the lab processing of the soil
# [25] Other comments
# [26] If the samples are concerned by scale calibration issue
# [27] Who performed the task

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
#soilbulk_full[is.na(soilbulk_full$CoreVol),] # two NA in IS1 & OG1 -> check datasheet -> lab incident, the 2 cores were discarded
#hist(soilbulk_full$CoreVol) # Volumes range from 120 to 450 cm3, most above 360 cm3 -> low volume = bad estimation of BD, too low volumes should be discarded

# W0 - Weight of fresh soil before saturation
#soilbulk_full[is.na(soilbulk_full$W0g),] # same two NAs -> validated
#hist(soilbulk_full$W0g) # Weights range from 50 to 180 g in a normal distribution -> very low weights likely to be linked to low volumes

# WSAT - Soil weight after water saturation
#soilbulk_full[is.na(soilbulk_full$WSAT),] # same two NAs -> validated
#hist(soilbulk_full$WSAT) # Weights range from 60 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes

# W24H - Soil weight after 24h of drying
#soilbulk_full[is.na(soilbulk_full$W24H),] # same two NAs -> validated
#hist(soilbulk_full$W24H) # Weights range from 50 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes, not so much difference compared to WSAT

# W48H - Soil weight after 48h of drying
#soilbulk_full[is.na(soilbulk_full$W48H),] # same two NAs -> validated
#hist(soilbulk_full$W48H) # Weights range from 50 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes

# WDRY - Soil weight after over at 105C
#soilbulk_full[is.na(soilbulk_full$WDRY),] # same two NAs -> validated
#hist(soilbulk_full$WDRY) # Weights range from 20 to 140 g in a normal distribution -> very low weights likely to be linked to low volumes

# Percent water loss in 24h
#soilbulk_full[is.na(soilbulk_full$percent_Waterloss24h),] # same two NAs -> validated
#hist(soilbulk_full$percent_Waterloss24h) # % range from -20% to 40%, main between 0 and 5% -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Percent water loss in 48h
#soilbulk_full[is.na(soilbulk_full$percent_Waterloss48h),] # same two NAs -> validated
#hist(soilbulk_full$percent_Waterloss48h) # % range from -70% to 60%, main between 0 and 5% -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Bulk density
#soilbulk_full[is.na(soilbulk_full$BD),] # same two NAs -> validated
#hist(soilbulk_full$BD) # % range from -0.1% to 0.4%, normal distribution -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Soil moisture in percentage weight
#soilbulk_full[is.na(soilbulk_full$Weightpercent_Soilmoisture),] # same two NAs -> validated
#hist(soilbulk_full$Weightpercent_Soilmoisture) # % range from -50% to 100%, normal distribution -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Soil moisture in percentage volume
#soilbulk_full[is.na(soilbulk_full$Volpercent_Soilmoisture),] # same two NAs -> validated
#hist(soilbulk_full$Volpercent_Soilmoisture) # % range from -50% to 100%, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

# Soil porosity
#soilbulk_full[is.na(soilbulk_full$percent_Soilporosity),] # same two NAs -> validated
#hist(soilbulk_full$percent_Soilporosity) # % range from -50% to 100%, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

# WFPS
#soilbulk_full[is.na(soilbulk_full$percent_WFPS),] # same two NAs -> validated
#hist(soilbulk_full$percent_WFPS) # % range from -50% to 100%, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

#
## Data filtering

# Min soil core volume
#filter(soilbulk_full, CoreVol<350 & !is.na(CoreVol)) # 136 or 9% samples unfit
#filter(soilbulk_full, CoreVol<380 & !is.na(CoreVol)) # 366 or 25% samples unfit
#filter(soilbulk_full, CoreVol<400 & !is.na(CoreVol)) # 690 or 46% samples unfit

# Negative water loss values
#filter(soilbulk_full, percent_Waterloss24h<0 & !is.na(percent_Waterloss24h)) #138 samples with water loss 24h negative
#filter(soilbulk_full, percent_Waterloss48h<0 & !is.na(percent_Waterloss48h)) #84 samples with water loss 48h negative

# Selection data with min 350 cm3 soil volume and positive water loss
soilbulk_full <- subset(soilbulk_full, CoreVol>350)
soilbulk_full <- subset(soilbulk_full, percent_Waterloss24h>0)
soilbulk_full <- subset(soilbulk_full, percent_Waterloss48h>0)

# Check new variable distribution - water loss 24h
#hist(soilbulk_full$percent_Waterloss24h) # still some extreme values over 20%
#filter(soilbulk_full, percent_Waterloss24h>20) # 8 cores with more than 20% over 24h
# 2 cores from UC1, which is excluded from the analysis -> should be removed
# 2 cores from OC3, concerned with scale issue (lots of negative values which are already removed). Water loss between 0-24 and 24-48 not coherent -> should be removed
# 2 cores from OC2, concerned with scale issue. Water loss between 0-24 and 24-48 not coherent with other samples from same plot (W48h>W24h for P1-D1_2) -> should be removed
# 1 cores from OC5, concerned with scale issue. Water loss between 0-24 and 24-48 not coherent with other samples from same site -> should be removed
# 1 core from OG6, not concerned by scale issue. High value coherent with other samples from same plot -> to be kept
soilbulk_full <- subset(soilbulk_full, percent_Waterloss24h<20 | BDcoreID == "OG6-P2-D4_2")

# Check new variable distribution - water loss 48h
#hist(soilbulk_full$percent_Waterloss48h) # still some extreme values over 30%
#filter(soilbulk_full, percent_Waterloss48h>30) # 2 cores with more than 30% over 48h
# 1 core from OG6, same one than with more than 20% reduction in 24h -> to be kept
# 1 core from IS3, with P3 concerned with scale issue. Water loss ~60%, completely incoherent with other samples from same plot -> to be removed
soilbulk_full <- subset(soilbulk_full, percent_Waterloss48h<40)

# Check new variable distribution - bulk density
#hist(soilbulk_full$BD) # no negative values anymore, quite nice normal distribution -> validated

# Check new variable distribution - soil moisture in percent weight
#hist(soilbulk_full$Weightpercent_Soilmoisture) # still some negative and extreme values (100%)
#filter(soilbulk_full, Weightpercent_Soilmoisture<20) # 8 cores with less than 20% soil moisture
# 4 cores from OC2, concerned with scale issue. 2 negative values, 2 other very low not coherent with other samples from the plot -> to be removed
# 1 core from OC4, concerned with scale issue -> weight loss after 24h almost 0, value not coherent with other samples from the plot -> to be removed
# 3 cores from OG4, concerned with scale issue -> values are coherent within the plot and relatively close to what is find in other plots (30%-40%) -> to be kept
soilbulk_full <- subset(soilbulk_full, Weightpercent_Soilmoisture>20 | SiteID == "OG4")
#filter(soilbulk_full, Weightpercent_Soilmoisture>90) # 1 core from IS3, with P3 concerned with scale issue. Incoherent with other samples from same plot -> to be removed
soilbulk_full <- subset(soilbulk_full, Weightpercent_Soilmoisture<90) # Distribution a bit hectic

# Check new variable distribution - soil moisture in percent volume
#hist(soilbulk_full$Volpercent_Soilmoisture) # still some very low values, but no extremes, nice normal distribution

# Check new variable distribution - WFPS
#hist(soilbulk_full$percent_WFPS) # no extremes, nice normal distribution

# Check new number of replicates per site
#sort(table(soilbulk_full$SiteID)) # 3 sites with less than 20 replicates (OS5, OC3, OC2) -> validated


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
chem2020_DM$`Tørrstoff/g/100g` <- as.numeric(chem2020_DM$`Tørrstoff/g/100g`)

#
## Common ID for merging

names(chem2019) <- gsub("Prøvenummer", "PlotID", names(chem2019)) # common ID
names(chem2020) <- gsub("Prøvenummer", "PlotID", names(chem2020)) # common ID
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
names(chem2020_DM) <- gsub("Tørrstoff/g/100g", "DryMatter_percent", names(chem2020_DM))
chem2020 <- left_join(chem2020, chem2020_DM)

# Merging 2019 and 2020 tables
soilchem_raw <- full_join(chem2019, chem2020)

#
## Name & character cleaning

# R friendly variable names
names(soilchem_raw) <- gsub("Jordart", "SoilType", names(soilchem_raw))
names(soilchem_raw) <- gsub("Leirklasse", "ClayCategory", names(soilchem_raw))
names(soilchem_raw) <- gsub("Glødetap", "LOI", names(soilchem_raw))
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
soilchem_raw <- subset(soilchem_raw, select = -c(Årstall, Journalnr, Navn, Postnr, Poststed, Registreringsdato, batchCode, contactName, Adresse, samplePartnerCode)) # remove useless columns

# New ID variables
soilchem_raw$PlotID <- paste(substr(soilchem_raw$PlotID, 1, 3), substr(soilchem_raw$PlotID, 4, 4), sep= "-P") #Create PlotID same format as other sheets
soilchem_raw$SiteID <- substr(soilchem_raw$PlotID, 1, 3) #Create SiteID

#
## Data cleaning - New R object

soilchem_full <- soilchem_raw

#
## Char var - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(soilchem_full$SiteID) # 3 samples per site - validated

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

# Correction latin name
vege_raw <- vege_raw |>  
  mutate(Species = dplyr::recode(Species, "Carex ovalis" = "Carex leporina")) |> 
  mutate(Species = dplyr::recode(Species, "Polytrichiastrum" = "Polytrichastrum")) |> 
  mutate(Species = dplyr::recode(Species, "Polytrichum juniperirinum" = "Polytrichum juniperinum")) |> 
  mutate(Species = dplyr::recode(Species, "Rhacomitrium" = "Racomitrium")) |> 
  mutate(Species = dplyr::recode(Species, "Kystkransmose" = "Rhytidiadelphus loreus")) |> 
  mutate(Species = dplyr::recode(Species, "Sphagunm" = "Sphagnum"))

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

# from 676 to 673 variables -> 3 quadrats removed?

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
## Table formatting for vegan beetle only (long table)
beetle_full <- subset(arthro_full, select = c(SiteID, SampleID, Staphylinidae, Carabidae, Hydrophilidae, Scarabaeidae, Ptiliidae, Curculionidae, Elateridae, Rhizophagidae, Leiodidae, Silphidae, Histeridae, Geotrupidae, Chrysomelidae, Dascillidae, Other, Erotylidae))
beetle_full <- beetle_full |> 
  pivot_longer(
    cols = c(Staphylinidae, Carabidae, Hydrophilidae, Scarabaeidae, Ptiliidae, Curculionidae, Elateridae, Rhizophagidae, Leiodidae, Silphidae, Histeridae, Geotrupidae, Chrysomelidae, Dascillidae, Other, Erotylidae),
    names_to = "BeetleFamilies", 
    values_to = "BeetleFam_abundance") |> 
  mutate(PlotID = substr(SampleID, 1, 6)) |>
  mutate(SiteID = substr(SampleID, 1, 3))

#
## Table formatting for vegan (long table)
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
  mutate(SiteID = substr(SampleID, 1, 3))

#
## Export clean data in new excel file

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

#str(mesobio_raw) # RAS

#
## Character cleaning, Common ID and correct Latin names for merging

# Common ID
names(mesobio_raw) <- gsub("Sample_name_on_pot", "SampleID", names(mesobio_raw))

#
## Char var Site Info - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
table(mesobio_raw$SiteID) # Missing IG3, OG2, OG3, OG5; IC3; tous les OC

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
hist(mesobio_raw$Acari) # Poisson distribution
hist(mesobio_raw$Collembola) # Poisson distribution

#
## Extract survey data only

mesobio_full <- mesobio_raw |> 
  filter(SiteID != "OY")

#
## Export clean data in new excel file

write_csv(mesobio_full, "data/cleandata/NBR_FullMesobio.csv")


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
names(precipitation) <- gsub("mean", "AnnualPrecipitation", names(precipitation))

#
## Max temperature data

# Import tiff files according to July month
maxtemp_listJuly <- list.files(path = "./data/maxtempraster", pattern= '_07_.*\\.tif$', all.files=TRUE, full.names=TRUE)

# Import raster files within one object
maxtemp_dailyJuly <- terra::rast(maxtemp_listJuly)

# Calculate mean
maxtemp_meanJuly <- terra::app(maxtemp_dailyJuly, mean, NA.RM = TRUE)
plot(maxtemp_meanJuly)

# Fill missing cell with average from closest cells
fullmaxtemp_meanJuly <- terra::focal(maxtemp_meanJuly, w=7, fun = "mean", na.policy = "only")
plot(fullmaxtemp_meanJuly)

# Extract mean values for NBR site coordinates
maxtempjuly <- as.data.frame(terra::extract(fullmaxtemp_meanJuly, subset(sitecoord, select = c(Xsite, Ysite))))
names(maxtempjuly) <- gsub("focal_mean", "MaxTempJuly", names(maxtempjuly))

#
## Min temperature data

# Import tiff files according to January month
mintemp_listJan <- list.files(path = "./data/mintempraster", pattern= '_01_.*\\.tif$', all.files=TRUE, full.names=TRUE)

# Import raster files within one object
mintemp_dailyJan <- terra::rast(mintemp_listJan)

# Calculate mean
mintemp_meanJan <- terra::app(mintemp_dailyJan, mean, NA.RM = TRUE)
plot(mintemp_meanJan)

# Fill missing cell with average from closest cells
fullmintemp_meanJan <- terra::focal(mintemp_meanJan, w=7, fun = "mean", na.policy = "only")
plot(fullmintemp_meanJan)

# Extract mean values for NBR site coordinates
mintempJan <- as.data.frame(terra::extract(fullmintemp_meanJan, subset(sitecoord, select = c(Xsite, Ysite))))
names(mintempJan) <- gsub("focal_mean", "MinTempJan", names(mintempJan))

#
## Preparation final dataset

# Merging all climate results
climate_full <- purrr::reduce(list(sitecoord, precipitation, maxtempjuly, mintempJan), dplyr::left_join)

# Rescaling - original data at 0.1 scale - temp conversion to C
climate_full <- climate_full |> 
  mutate(annualprecipitation = AnnualPrecipitation/10,
         MaxTempJuly = MaxTempJuly/100-273.15,
         MinTempJan = MinTempJan/100-273.15)

# Clean data + export in csv
climate_full <- subset(climate_full, select = -c(ID, Xsite, Ysite))
write_csv(climate_full, "data/cleandata/NBR_FullClimate.csv")

