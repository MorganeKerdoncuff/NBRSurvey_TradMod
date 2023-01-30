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

# R friendly variable names
names(siteinfo_raw) <- gsub("%", "percent", names(siteinfo_raw))
names(siteinfo_raw) <- gsub(" ", "", names(siteinfo_raw)) 
names(siteinfo_raw) <-  gsub("\\(", "_", names(siteinfo_raw))
names(siteinfo_raw) <-  gsub("\\)", "", names(siteinfo_raw))

# Removal empty columns
#siteinfo_raw <- subset(siteinfo_raw, select = -c(Size_livestock, Surface)) # remove empty variables, which will be included in another dataset

#
## Data cleaning - new R object

siteinfo_full <- siteinfo_raw

#
## Numeric var Site Info - Check min/max, distribution and potential outliers

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
unique(siteinfo_full$NiBioEcologicalZone) # Three categories for ecological zone - validated
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

#
## Geological information # Aborted as raw data were incorrect


#### FIELD MANAGEMENT ####

## Description

## List of variables

# [1] Field identification code for data collection
# [2] Geographical location of the field, hamlet
# [3] Geographical location of the field, postcode
# [4] Geographical location of the field, municipality (both former and current classification)
# [5] Name of the farmer
# [6] Date of the interview with the farmer
# [7] Site productivity - infield high productivity, outfield low productivity
# [8] National identification number of the farm
# [9] Main livestock, currently grazing in the field - ! rotational grazing management
# [10] Other livestock, grazing in other fields - ! rotational grazing management
# [11] Name of cow breed
# [12] Name of goat breed
# [13] Name of sheep breed
# [14] Number of adult animals in the main livestock
# [15] Number of young animals in the main livestock
# [16] Number of adult animals in other livestock
# [17] Number of young animals in other livestock
# [18] Size of the field containing the sampling area (ha) ! rotational grazing management
# [19] FROM GÅRDSKSART - Total grazing area of the farm (ha) ! rotational grazing management
# [20] Grazing density on the field - not collected, empty column
# [21] Period since the farmer has had the current livestock
# [22] Period without grazing management on the field
# [23] Number of months the main livestock grazes in the field in a year
# [24] Number of months other livestock graze in the field in a year ! rotational grazing
# [25] If the animals are kept inside or outside at night
# [26] Farmer impression of grazing pressure on the site
# [27] Former livestock which used to graze in the field (1)
# [28] Former livestock which used to graze in the field (2)
# [29] Former livestock which used to graze in the field (3)
# [30] Former livestock grazing period (1)
# [31] Former livestock grazing period (2)
# [32] Former livestock grazing period (3)
# [33] Type of farm management during the past 10 years (1)
# [34] Type of farm management during the past 10 years (2)
# [35] Type of farm management during the past 10 years (3)
# [36] If applicable, frequency of grass cutting
# [37] If applicable, frequency of mulching
# [38] If applicable, frequency of fertilization with manure
# [39] If applicable, type of manure used
# [40] If applicable, volume of manure used (m3)
# [41] If applicable, in which season the manure is used
# [42] If applicable, frequency of fertilization with artificial fertilizer
# [43] If applicable, weight of artificial fertilizer used (kg)
# [44] If applicable, in which season the artificial fertilizer is used
# [45] If applicable, frequency of fertilization with shell-sand/lime
# [46] If applicable, weight of shell-sand/limer used (kg)
# [47] If applicable, in which season the shell-sand/lime is used
# [48] If applicable, last time the field was sowed
# [49] If applicable, last time the field was plowed
# [50] If applicable, last time the field was drained
# [51] Type of land use prior to grazing
# [52] If applicable, type of production prior to grazing (1)
# [53] If applicable, type of production prior to grazing (2)
# [54] Time period of land use type prior to grazing
# [55] Comments

#
## Summary land use - Check table size, list of variables, variable types (num/chr)

str(landuse_raw) # All good

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
names(landuse_raw) <- gsub("Farmer_interviewee", "Farmer", names(landuse_raw))
names(landuse_raw) <- gsub("Farm_name_nr.orbeitelag", "FarmID", names(landuse_raw))
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
names(landuse_raw) <- gsub("1...34", "2", names(landuse_raw))
names(landuse_raw) <- gsub("1...35", "3", names(landuse_raw))
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
#unique(landuse_full$FieldManagement2)
#unique(landuse_full$FieldManagement3)
#landuse_full[is.na(landuse_full$FieldManagement1),] # 6 missing values -> farmers who did not reply to the survey
#landuse_full[is.na(landuse_full$FieldManagement2),] # 36 missing values
#landuse_full[is.na(landuse_full$FieldManagement3),] # 41 missing values
table(landuse_full$FieldManagement1) # 9 sites with no treatments at all, 21 with at least fertilization
table(landuse_full$FieldManagement2) # 8 sites with at least 2 treatments
table(landuse_full$FieldManagement3) # 3 sites with 3 treatments

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
filter(landuse_full, Fertilizing_type2 != is.na(Fertilizing_type2) & Fertilizing_type3 == is.na(Fertilizing_type3)) # 5 grasslands (OC1, IS1, OC4, OG6, OC5) and 1 coastal heathland (OS8)

table(landuse_full$Fertilizing_type1) # 21-11 = 10 sites with only 1 type of fertilization
subset(landuse_full, Fertilizing_type1 != is.na(Fertilizing_type1)) # 

# Frequency of manure
#unique(landuse_full$Manure_freq) # possibilities are "yearly" or "occasionally"
#subset(landuse_full, Manure_freq == "every year") # 9 grasslands (OC1, IC1, IS1, OC3, OC4, OG6, IC4, IC5, OC5) with annual manure
#subset(landuse_full, Manure_freq == "sometimes") # 1 grassland (OC2) and 1 heathland (OS8) with occasional manure

# Frequency of artificial fertilization
#unique(landuse_full$ArtificialFert_freq) # possibilities are "yearly" or "occasionally"
#subset(landuse_full, ArtificialFert_freq == "every year") # 15 grasslands (OC1, OS2, IC1, OG3, IC3, OC3, OC4, OG2, OG4, IC4, IC5, OC5, OS6, OC2, IS4) with annual fertilization with artificial fertiliser
#subset(landuse_full, ArtificialFert_freq == "sometimes") # 3 grasslands (OC2, OG6, IV1) and 1 heathland (OS8) with occasional fertilization with artificial fertiliser

# Frequency of inorganic fertilization with shell-sand/lime
#unique(landuse_full$ShellSandLime_freq) # possibilities are "occasionally" only
#subset(landuse_full, ShellSandLime_freq == "sometimes") # 7 grassland (IC1, OG3, IC3, OC3, IC4, IC5, OC2)

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
hist(landuse_full$FlockSize1_adults) # wide range but most farms under 50 animals - one farm over 150 animals
landuse_full[landuse_full$FlockSize1_adults>150,] # IS4 over 150 animals + US6 as NA

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
hist(landuse_full$FlockSize1_young)
landuse_full[landuse_full$FlockSize1_young>150,] # US3 & IS5 + 4 farms with no young animals + US6 as NA

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
landuse_full[is.na(landuse_full$YearlyGrazing1_month),] # 6 missing values from farmers who did not reply the survey
hist(landuse_full$YearlyGrazing1_month) # coherent values (between 2 and 12), no visible outliers -> validated

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
hist(landuse_full$Grazingdensity_perha)
landuse_full[is.na(landuse_full$Grazingdensity_perha),]


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

# R friendly variable names
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

#
## Calculation new variables

# Heat Load Index
area20x20_full <- area20x20_full |> 
  mutate(HLI = cos(AspectDegree-225)*tan(General_slope))
#hist(area20x20_full$HLI) # 3 outliers: one under 200, two over 100
#area20x20_full[area20x20_full$HLI>100,] #OG4 & IS3 -> both 11 degree slope with SW & SE exposition
#area20x20_full[area20x20_full$HLI<0,] #OS6 -> 11 degree slope with NE exposition