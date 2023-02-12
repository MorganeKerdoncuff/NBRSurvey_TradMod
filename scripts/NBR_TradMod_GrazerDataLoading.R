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
  summarise(PlantSp_cover = mean(Abundance))

# Beetle community - current at pitfall level -> summary by sum only main families
arthro_grass <- arthro_grass |> 
  group_by(SiteID, BeetleFamilies) |> 
  summarise(BeetleFam_abundance = sum(BeetleFam_abundance, na.rm = TRUE))

#
## Transformation community data

# Data distribution
hist(vege_grass$PlantSp_cover) # Poisson, highly skewed -> should remove rare species
hist(arthro_grass$BeetleFam_abundance) # Poisson, highly skewed

# Removal rare plant species with total mean cover across sites under 0.5% -> from 275 to 88 species
vege_grass <- vege_grass |> 
  pivot_wider(names_from = Species, values_from = PlantSp_cover) |>  
  select_if(negate(function(col) is.numeric(col) && sum(col) < 0.5)) |> 
  pivot_longer(cols = c(-SiteID), names_to = "PlantSp", values_to = "PlantSp_cover")

# Removal rare beetle families (determined in data cleaning script)
arthro_grass <- arthro_grass |> 
  pivot_wider(names_from = BeetleFamilies, values_from = BeetleFam_abundance)
arthro_grass <- subset(arthro_grass, select = c(SiteID, Carabidae, Staphylinidae, Hydrophilidae, Ptiliidae, Scarabaeidae, Curculionidae, Elateridae))
arthro_grass <- arthro_grass |>  
  pivot_longer(cols = c(-SiteID), names_to = "BeetleFam", values_to = "BeetleFam_abundance")

# Log transformation
vege_grass <- vege_grass |> 
  mutate(PlantSp_logcover = log1p(PlantSp_cover))
arthro_grass <- arthro_grass |> 
  mutate(BeetleFam_logabundance = log1p(BeetleFam_abundance))

# Contingency tables
contin_vege <- xtabs(formula = PlantSp_logcover ~ SiteID + PlantSp, data = vege_grass)
contin_beetle <- xtabs(formula = BeetleFam_logabundance ~ SiteID + BeetleFam, data = arthro_grass)

#
## Creation fjord system matrix

# Desired variables
## Elevation (num), from area 20x20
## Slope angle (num), from area 20x20
## Aspect angle (num), from area 20x20
## Distance to sea (num), from area 20x20

# Selection variables
fjordsys <- subset(area20x20_grass, select = c(SiteID, Elevation_max, General_slope, AspectDegree, DistanceToSea_m))

# Contingency table
fjordsys <- fjordsys |> 
  pivot_longer(
    cols = c(Elevation_max, General_slope, AspectDegree, DistanceToSea_m),
    names_to = "Factors",
    values_to = "Values")
contin_fjordsys <- xtabs(formula = Values ~ SiteID + Factors, data = fjordsys)

# Data scaling
contin_fjordsys <- scale(contin_fjordsys)

#
## Creation landscape matrix

# Desired variables
## Sum forest area (num), from landscape
## Sum cultivated area (num), from landscape
## Sum infield area (num), from landscape
## Sum outfield area (num), from landscape
## Sum wetland area (num), from landscape

# New variables - sum of cultivated land & forest
landscape_grass <- landscape_grass |> 
  mutate(TotCultivatedLand_percent = FullyCultivatedLand_percent + SuperficiallyCultivatedLand_percent) |> 
  mutate(TotForest_percent = ProductiveForest_percent + NonProductiveForest_percent)

# Selection variables
landscape <- subset(landscape_grass, select = c(SiteID, TotCultivatedLand_percent, TotForest_percent, Infield_percent, Outfield_percent, Wetland_percent))

# Contingency table
landscape <- landscape |> 
  pivot_longer(
    cols = c(TotCultivatedLand_percent, TotForest_percent, Infield_percent, Outfield_percent, Wetland_percent),
    names_to = "Factors",
    values_to = "Values")
contin_landscape <- xtabs(formula = Values ~ SiteID + Factors, data = landscape)

# Data scaling
contin_landscape <- scale(contin_landscape)

#
## Creation grazing management matrix

# Desired variables
## Type of livestock (char), from landuse
## Number of adult animals (num), from landuse
## Grazing surface of the collected field (num), from landuse
## Total infield surface (num), from landuse
## Standardised grazing density, from landuse

# Dummy numeric for original character variables
landuse_grass <- landuse_grass |> 
  mutate(Sheep = ifelse(Livestock1 == "sheep",1,0)) |> 
  mutate(Cow = ifelse(Livestock1 == "cow", 1,0))

# Selection variables
grazing <- subset(landuse_grass, select = c(SiteID, Sheep, Cow, FlockSize1_adults, GrazingSurface_ha, TotalInfieldSurface, Grazingdensity_perha))

# Contingency table
grazing <- grazing |> 
  pivot_longer(
    cols = c(Sheep, Cow, FlockSize1_adults, GrazingSurface_ha, TotalInfieldSurface, Grazingdensity_perha),
    names_to = "Factors",
    values_to = "Values")
contin_grazing <- xtabs(formula = Values ~ SiteID + Factors, data = grazing)

# Data scaling
contin_grazing <- scale(contin_grazing)

#
## Creation plant local environment matrix

# Desired variables
## Bulk density (num), from soilbulk
## Moisture content (num), from soilbulk
## Penetration rate (num), from soilpene
## LOI (num), from soilchem
## Nitrogen content (num), from soilchem
## Phosphorous content (num), from soilchem
## pH (num), from soilchem
## Humus content (num), from soilchem

# Selection variables
plantenvi <- purrr::reduce(list(soilbulk_grass, soilpene_grass, soilchem_grass), dplyr::left_join)
plantenvi <- subset(plantenvi, select = -c(MeanPotassium, MeanSodium))

# Contingency table
plantenvi <- plantenvi |> 
  pivot_longer(
    cols = c(MeanBD, MeanMoisture, MeanPT, MeanLOI, MeanHumus, MeanpH, MeanPhosphorus, MeanNitrogen),
    names_to = "Factors",
    values_to = "Values")
contin_plantenvi <- xtabs(formula = Values ~ SiteID + Factors, data = plantenvi)

# Data scaling
contin_plantenvi <- scale(contin_plantenvi)

#
## Creation beetle local environment matrix

# Desired variables
## Exposed ground (bare soil + rock) cover (num), from groundcover
## Litter cover (num), from groundcover
## Mean vegetation height (num), from groundcover
## Mean vegetation richness (num), from groundcover
## Bulk density (num), from soilbulk
## Moisture content (num), from soilbulk
## Humus content (num), from soilchem

# Selection variables
beetlenvi <- purrr::reduce(list(groundcover_grass, soilbulk_grass, soilchem_grass), dplyr::left_join)
beetlenvi <- beetlenvi |> 
  mutate(MeanExposedGround = MeanBareSoil + MeanRocks)
beetlenvi <- subset(beetlenvi, select = c(SiteID, MeanExposedGround, MeanLitter, MeanHeight, MeanRichness, MeanBD, MeanMoisture, MeanHumus))

# Contingency beetlenvi
beetlenvi <- beetlenvi |> 
  pivot_longer(
    cols = c(MeanExposedGround, MeanLitter, MeanHeight, MeanRichness, MeanBD, MeanMoisture, MeanHumus),
    names_to = "Factors",
    values_to = "Values")
contin_beetlenvi <- xtabs(formula = Values ~ SiteID + Factors, data = beetlenvi)

# Data scaling
contin_beetlenvi <- scale(contin_beetlenvi)

#
## Check variable colinearity
