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
landuse_raw <- read_excel(path = "data/Farmer_Survey.xlsx", sheet="Farms Information_R") #Field management data from farmer interview, at site level
farmer_raw <- read_excel(path = "data/NBR_FarmRepertoire.xlsx") #Farm management data from farmer interview med Margit, at site level
landscape_raw <- read_excel(path = "data/NBR_LandscapeMatrix.xlsx", sheet="MatrixProportion") #Land cover data around the fields from Geonorge, at site level
site20x20_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="20mX20m") #Sampling area description, vegetation cover at the site level
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
## Numeric var Site Info - Check min/max, distribution and potential outliers

# Check min/max
# test <- siteinfo_raw %>% 
  # summarise(tibble(across(where(is.numeric), ~min(.x, na.rm = TRUE), .names = "min_{.col}"),across(where(is.numeric), ~max(.x, na.rm = TRUE), .names = "max_{.col}"))) %>% 
  # transpose() # All good

# Check distribution of quantitative variable
# hist(siteinfo_raw$`Pooestimation_g-10percent`) # Standard distribution of poo data - validated

#
## Char var Site Info - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# table(siteinfo_raw$SiteID) # Unique ID for each site - validated
# unique(siteinfo_raw$NiBioEcologicalZone) # Three categories for ecological zone - validated
# unique(siteinfo_raw$Location) # Correct locations and location names
# unique(siteinfo_raw$Municipality) # missing data + US2 Stordalen -> did you guys went that far ? Yes validated
# siteinfo_raw[is.na(siteinfo_raw$Municipality),] # NA identified - US2, Stordalen is located in Masfjorden
siteinfo_raw <- siteinfo_raw |> 
  mutate(Municipality=ifelse(SiteID == "US2", "Masfjorden", Municipality)) # Replace NA by municipality name
# unique(siteinfo_raw$Municipality) # NA in Municipality replaced - validated
unique(siteinfo_raw$Farmer) # Missing two full names and one NA
#siteinfo_raw[is.na(siteinfo_raw$Farmer)]
#unique(siteinfo_raw$Type_livestock)
siteinfo_raw <- siteinfo_raw %>% 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "Villsau" = "Sheep")) # Villsau as sheep category
table(siteinfo_raw$Type_livestock) #10 cows, 11 goats, 23 sheep -> ask Amy again about UG1 & UG2
#unique(siteinfo_raw$Habitat) #needs to standardise categories, no more than 5 - burnt/non-burnt? alpine/subalpine? 
#unique(siteinfo_raw$Animal_on_site)

#Binding geology
Geol_raw <- read_xlsx("..\\06_SurveyDataSets\\Tradmod site geology.xlsx")
#length(Geol_raw$SiteID) #45 sites
#sort(table(Geol_raw$SiteID)) # 1 value per site
#unique(Geol_raw$Geology) # 2 types of bedrock, no doubletons
Site_Info_raw <- left_join(Site_Info_raw, Geol_raw)