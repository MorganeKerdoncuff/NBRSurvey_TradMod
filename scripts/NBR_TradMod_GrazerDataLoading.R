# TradMod WP2 - NBR survey grazer loading script
#Description of the data
#Date
#Who
#Project
#Funding
#Place

#### PACKAGE LOADING ####

library(tidyverse) # R language
library(purrr) # Data manipulation: function "reduce" to bind several tables at the same time



#### DATA LOADING ####

siteinfo_full <- read.csv("data/cleandata/NBR_FullSiteInfo.csv", sep=",") # Clean site info data
landscape_full <- read.csv("data/cleandata/NBR_FullLandscape.csv", sep=",") # Clean landscape matrix data
landuse_full <- read.csv("data/cleandata/NBR_FullLanduse.csv", sep=",") # Clean field management data
area20x20_full <- read.csv("data/cleandata/NBR_FullArea20x20.csv", sep=",") # Clean sampling area 20mx20m data
groundcover_full <- read.csv("data/cleandata/NBR_FullGroundCover.csv", sep=",") # Clean aboveground cover data
soilbulk_full <- read.csv("data/cleandata/NBR_FullSoilBulk.csv", sep=",") # Clean soil bulk density data
soilchem_full <- read.csv("data/cleandata/NBR_FullSoilChem.csv", sep=",") # Clean soil chemistry data
soilpene_full <- read.csv("data/cleandata/NBR_FullSoilPene.csv", sep=",") # Clean soil bulk density data
vege_full <- read.csv("data/cleandata/NBR_FullPlantComm.csv", sep=",") # Clean plant community data
arthro_full <- read.csv("data/cleandata/NBR_FullArtComm.csv", sep=",") # Clean arthropod community data


#### DATA PREPARATION ####

#
## Filter grassland sites

# Site selection
siteinfo_grass <- siteinfo_full |>  
  filter(Habitat == "permanent grassland") |> 
  dplyr::select(SiteID, Type_livestock)

# Extraction in other datasets
landscape_grass <- filter(landscape_full, SiteID %in% siteinfo_grass$SiteID)
landuse_grass <- filter(landuse_full, SiteID %in% siteinfo_grass$SiteID)
area20x20_grass <- filter(area20x20_full, SiteID %in% siteinfo_grass$SiteID)
groundcover_grass <- filter(groundcover_full, SiteID %in% siteinfo_grass$SiteID)
soilbulk_grass <- filter(soilbulk_full, SiteID %in% siteinfo_grass$SiteID)
soilchem_grass <- filter(soilchem_full, SiteID %in% siteinfo_grass$SiteID)
soilpene_grass <- filter(soilpene_full, SiteID %in% siteinfo_grass$SiteID)
vege_grass <- filter(vege_full, SiteID %in% siteinfo_grass$SiteID)
arthro_grass <- filter(arthro_full, SiteID %in% siteinfo_grass$SiteID)

#
## Summarise data at site level -> should be 30 observations for non community data

# Site info - validated

# Landscape matrix - validated

# Field management - validated

# Area 20 x 20 - validated

# Ground cover - current at sample level -> summary by average
groundcover_grass <- groundcover_grass |>  
  group_by(SiteID) |>  
  summarise(MeanBareSoil = mean(Bare_soil, na.rm=TRUE),
            MeanRocks = mean(Rocks, na.rm=TRUE),
            MeanLitter = mean(Litter, na.rm=TRUE),
            MeanDeadWood = mean(Dead_wood, na.rm=TRUE),
            MeanBryo = mean(Bryophytes, na.rm=TRUE),
            MeanLichen = mean(Lichen, na.rm=TRUE),
            MeanVasc = mean(Vascular, na.rm=TRUE),
            MeanBlossom = mean(Blossom_cover, na.rm=TRUE),
            MeanHeight = mean(VG_mean_height_cm, na.rm=TRUE),
            MaxHeight = mean(VG_max_height_cm, na.rm=TRUE),
            MeanRichness = mean(Plant_species_richness, na.rm=TRUE))

# Bulk density - current at sample level -> summary by average
soilbulk_grass <- soilbulk_grass |> 
  group_by(SiteID) |> 
  summarise(MeanBD = mean(BD),
            MeanMoisture = mean(Weightpercent_Soilmoisture))

# Soil chemistry - current at plot level -> summary by average
soilchem_grass <- soilchem_grass |> 
  group_by(SiteID) |> 
  summarise(MeanLOI = mean(LOI),
            MeanHumus = mean(Humus_percentDM),
            MeanpH = mean(pH),
            MeanPhosphorus = mean(P.Al_mg.100g),
            MeanPotassium = mean(K.Al_mg.100g),
            MeanNitrogen = mean(TotalN_percentDM),
            MeanSodium = mean(Na.Al_mg.100g))

# Soil penetration - current at sample level -> summary by average
soilpene_grass <- soilpene_grass |> 
  group_by(SiteID) |> 
  summarise(MeanPT = mean(AveragePT_cm))

# Plant community - current at quadrat level -> summary by average
vege_grass <- vege_grass |> 
  group_by(SiteID, Species) |> 
  summarise(PlantSp_Cover = mean(Abundance))

# Beetle community - current at pitfall level -> summary by sum only main families
arthro_grass <- arthro_grass |> 
  group_by(SiteID, BeetleFamilies) |> 
  summarise(BeetleFam_Abundance = sum(BeetleFam_Abundance, na.rm = TRUE))

#
## Transformation community data

# Data distribution
hist(vege_grass$PlantSp_Cover) # Poisson, highly skewed -> should ren\move rare species
hist(arthro_grass$BeetleFam_Abundance) # Poisson, highly skewed