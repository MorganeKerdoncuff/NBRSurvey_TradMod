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
library(vegan) # Community ecology analysis
library(ggplot2) # Visual representation
library(GGally) # Extension ggplot
library(car) # Visualisation lm residuals plots
library(CCA) # Canonical correlation
library(CCP) # Statistical test canonical correlation

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
## Transformation plant community data

# Data distribution
#hist(vege_grass$PlantSp_cover) # Poisson, highly skewed -> should remove rare species

# Species frequencies all sites included
Vege_freq <- vege_grass |> 
  group_by(Species) |> 
  summarise_if(is.numeric, sum, na.rm = TRUE) |> 
  dplyr::arrange(desc(PlantSp_cover)) # Main grass species are Agrostis capillaris, Festuca rubra, Holcus lanatus, Poa pratensis, Deschampsia cespitosa, Anthoxantum odoratum, Lolium perenne, Deschampsia flexuosa and Poa trivialis.

# Removal rare plant species with total mean cover across sites under 15% -> from 274 to 28 species -> ABORTED
#vege_grass <- vege_grass |> 
  #pivot_wider(names_from = Species, values_from = PlantSp_cover) |>  
  #select_if(negate(function(col) is.numeric(col) && sum(col) < 30)) |> 
  #pivot_longer(cols = c(-SiteID), names_to = "PlantSp", values_to = "PlantSp_cover")

# Selection main grass species
vege_grass <- subset(vege_grass,
                       Species == "Agrostis capillaris" |
                       Species == "Festuca rubra" | 
                       Species == "Holcus lanatus" | 
                       Species == "Poa pratensis" | 
                       Species == "Deschampsia cespitosa" | 
                       Species == "Anthoxantum odoratum" | 
                       Species == "Lolium perenne" | 
                       Species == "Deschampsia flexuosa" |
                       Species == "Poa trivialis")
  
# Log transformation
vege_grass <- vege_grass |> 
  mutate(PlantSp_logcover = log1p(PlantSp_cover))

# Contingency tables & wide tables
contin_vege <- xtabs(formula = PlantSp_logcover ~ SiteID + Species, data = vege_grass)
vege <- subset(vege_grass, select = -c(PlantSp_cover))
vege <- vege |> 
  pivot_wider(names_from = Species, values_from = PlantSp_logcover)
vege <- as.data.frame(vege)

#
## Transformation beetle assemblage data

#hist(arthro_grass$BeetleFam_abundance) # Poisson, highly skewed

# Removal rare beetle families (determined in data cleaning script)
arthro_grass <- arthro_grass |> 
  pivot_wider(names_from = BeetleFamilies, values_from = BeetleFam_abundance)
arthro_grass <- subset(arthro_grass, select = c(SiteID, Carabidae, Staphylinidae, Hydrophilidae, Ptiliidae, Scarabaeidae, Curculionidae, Elateridae))
arthro_grass <- arthro_grass |>  
  pivot_longer(cols = c(-SiteID), names_to = "BeetleFam", values_to = "BeetleFam_abundance")

# Log transformation
arthro_grass <- arthro_grass |> 
  mutate(BeetleFam_logabundance = log1p(BeetleFam_abundance))

# Contingency tables & wide tables
contin_beetle <- xtabs(formula = BeetleFam_logabundance ~ SiteID + BeetleFam, data = arthro_grass)
beetle <- pivot_wider(arthro_grass, names_from = BeetleFam, values_from = BeetleFam_logabundance)
beetle <- as.data.frame(beetle)

#
## Creation fjord system matrix

# Desired variables
## Elevation (num), from area 20x20
## Slope angle (num), from area 20x20
## Aspect angle (num), from area 20x20
## Distance to sea (num), from area 20x20

# Selection variables & scaling
fjordsys <- subset(area20x20_grass, select = c(SiteID, Elevation_max, General_slope, AspectDegree, DistanceToSea_m))
fjordsys <- fjordsys |> 
  mutate(across(where(is.numeric), scale))

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

# Selection variables & scaling
landscape <- subset(landscape_grass, select = c(SiteID, TotCultivatedLand_percent, TotForest_percent, Infield_percent, Outfield_percent, Wetland_percent))
landscape <- landscape |> 
  mutate(across(where(is.numeric), scale))

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

# Selection variables & scaling
grazing <- subset(landuse_grass, select = c(SiteID, Sheep, Cow, FlockSize1_adults, GrazingSurface_ha, TotalInfieldSurface, Grazingdensity_perha))
grazing <- grazing |> 
  mutate(across(where(is.numeric), scale))

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

# Selection variables & scaling
locenvi_vege <- purrr::reduce(list(soilbulk_grass, soilpene_grass, soilchem_grass), dplyr::left_join)
locenvi_vege <- subset(locenvi_vege, select = -c(MeanPotassium, MeanSodium))
locenvi_vege <- locenvi_vege |> 
  mutate(across(where(is.numeric), scale))

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

# Selection variables & scaling
locenvi_beetle <- purrr::reduce(list(groundcover_grass, soilbulk_grass, soilchem_grass), dplyr::left_join)
locenvi_beetle <- locenvi_beetle |> 
  mutate(MeanExposedGround = MeanBareSoil + MeanRocks)
locenvi_beetle <- subset(locenvi_beetle, select = c(SiteID, MeanExposedGround, MeanLitter, MeanBryo, MeanHeight, MeanRichness, MeanBD, MeanMoisture, MeanHumus))
locenvi_beetle <- locenvi_beetle |> 
  mutate(across(where(is.numeric), scale))


#### VERIFICATION ASSUMPTIONS ####

#
## Variables colinearity

# Function for paired panels
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

# Fjord system variables
#paircor_fjordsys <- ggpairs(select_if(fjordsys, is.numeric))
#paircor_fjordsys # no colinearity between variables, one outlier on distance to sea
pairs(select_if(fjordsys, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # no colinearity between variables -> validated

# Landscape variables
#paircor_landscape <- ggpairs(select_if(landscape, is.numeric))
#paircor_landscape # one correlation a bit strong between total cultivated land and infield percent
pairs(select_if(landscape, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # one correlation a bit strong (0.74) between total cultivated land and infield percent

# Grazing variables
pairs(select_if(grazing, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # no correlation stronger than 0.64 -> validated

# Local plant environmental variables
pairs(select_if(locenvi_vege, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # 5 variables highly correlated (>0.8): BD, moisture, LOI, humus, Nitrogen

# Local beetle environmental variables
pairs(select_if(locenvi_beetle, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # 3 variables highly correlated (>0.8): BD, moisture & humus

# Removal strongly correlated variables - BD kept as with higher number of replication for each site
locenvi_vege <- subset(locenvi_vege, select = -c(MeanMoisture, MeanLOI, MeanHumus, MeanNitrogen))
locenvi_beetle <- subset(locenvi_beetle, select = -c(MeanMoisture, MeanHumus))

#
## Linear relationships of residuals

# Ordination community data
DCA_vege <- decorana(contin_vege)
DCA_vege # Axis length of 3.5 -> keep DCA
DCA_beetle <- decorana(contin_beetle)
DCA_beetle # Axis length under 1 -> run PCA
PCA_beetle <- prcomp(contin_beetle)

# Extraction plant species scores from DCA
scores_vege <- as.data.frame(scores(DCA_vege, choices = c(1,2), display = "sites"))
names(scores_vege) <- gsub("DCA1", "VegDCA1", names(scores_vege))
names(scores_vege) <- gsub("DCA2", "VegDCA2", names(scores_vege))
scores_vege <- scores_vege |>  
  mutate(SiteID = row.names(scores_vege)) # PlotID as column for future binding

# Extraction beetle family scores from PCA
scores_beetle <- as.data.frame(scores(PCA_beetle, choices = c(1,2), display = "sites"))
names(scores_beetle) <- gsub("PC1", "BeetlePCA1", names(scores_beetle))
names(scores_beetle) <- gsub("PC2", "BeetlePCA2", names(scores_beetle))
scores_beetle <- scores_beetle |>  
  mutate(SiteID = row.names(scores_beetle)) # PlotID as column for future binding

# Table binding
allvar_vege <- purrr::reduce(list(scores_vege, fjordsys, landscape, grazing, locenvi_vege), dplyr::left_join)
allvar_beetle <- purrr::reduce(list(scores_beetle, fjordsys, landscape, grazing, locenvi_beetle), dplyr::left_join)

# Residuals LM for plants - 28 species -> ABORTED
#residualPlots(lm(VegDCA1~Elevation_max, data = allvar_vege)) # Tukey=1.26 ; p-value=0.21 -> validated
#residualPlots(lm(VegDCA1~General_slope, data = allvar_vege)) # Tukey=0.35 ; p-value=0.72 -> validated
#residualPlots(lm(VegDCA1~AspectDegree, data = allvar_vege)) # Tukey=1.48 ; p-value=0.14 -> validated
#residualPlots(lm(VegDCA1~DistanceToSea_m, data = allvar_vege)) # Tukey=-1.24 ; p=0.21 -> validated
#residualPlots(lm(VegDCA1~TotCultivatedLand_percent, data = allvar_vege)) # Tukey=-0.027 ; p=0.98 -> validated
#residualPlots(lm(VegDCA1~TotForest_percent, data = allvar_vege)) # Tukey=-0.55 ; p=0.58 -> validated
#residualPlots(lm(VegDCA1~Infield_percent, data = allvar_vege)) # Tukey=-1.73 ; p=0.083 -> a bit tight
#residualPlots(lm(VegDCA1~Outfield_percent, data = allvar_vege)) # Tukey=0.83 ; p=0.41 -> validated
#residualPlots(lm(VegDCA1~Wetland_percent, data = allvar_vege)) # Tukey=0.79 ; p=0.43 -> validated
#residualPlots(lm(VegDCA1~Sheep, data = allvar_vege)) # binary, cannot be checked
#residualPlots(lm(VegDCA1~Cow, data = allvar_vege)) # binary, cannot be checked
#residualPlots(lm(VegDCA1~FlockSize1_adults, data = allvar_vege)) # Tukey=2.27 ; p=0.023 -> rejected (one outlier? -> no)
#residualPlots(lm(VegDCA1~GrazingSurface_ha, data = allvar_vege)) # Tukey=0.038 ; p=0.97 -> validated
#residualPlots(lm(VegDCA1~TotalInfieldSurface, data = allvar_vege)) # Tukey=1.65 ; p=0.096 -> a bit tight
#residualPlots(lm(VegDCA1~Grazingdensity_perha, data = allvar_vege)) # Tukey=0.82 ; p=0.41 -> validated
#residualPlots(lm(VegDCA1~MeanBD, data = allvar_vege)) # Tukey=0.95 ; p=0.34 -> validated
#residualPlots(lm(VegDCA1~MeanPT, data = allvar_vege)) # Tukey=-1.24 ; p=0.22 -> validated
#residualPlots(lm(VegDCA1~MeanpH, data = allvar_vege)) # Tukey=-0.12 ; p=0.9 -> validated
#residualPlots(lm(VegDCA1~MeanPhosphorus, data = allvar_vege)) # Tukey=-1.42 ; p=0.16 -> validated

# Residuals LM for grasses
#residualPlots(lm(VegDCA1~Elevation_max, data = allvar_vege)) # Tukey=0.98 ; p-value=0.33 -> validated
#residualPlots(lm(VegDCA1~General_slope, data = allvar_vege)) # Tukey=-0.79 ; p-value=0.43 -> validated
#residualPlots(lm(VegDCA1~AspectDegree, data = allvar_vege)) # Tukey=0.5 ; p-value=0.62 -> validated
#residualPlots(lm(VegDCA1~DistanceToSea_m, data = allvar_vege)) # Tukey=0.34 ; p=0.74 -> validated
#residualPlots(lm(VegDCA1~TotCultivatedLand_percent, data = allvar_vege)) # Tukey=0.3 ; p=0.76 -> validated
#residualPlots(lm(VegDCA1~TotForest_percent, data = allvar_vege)) # Tukey=-0.55 ; p=0.58 -> validated
#residualPlots(lm(VegDCA1~Infield_percent, data = allvar_vege)) # Tukey=-0.45 ; p=0.65 -> validated
#residualPlots(lm(VegDCA1~Outfield_percent, data = allvar_vege)) # Tukey=0.71 ; p=0.48 -> validated
residualPlots(lm(VegDCA1~Wetland_percent, data = allvar_vege)) # Tukey=2.05 ; p=0.04 -> rejected
#residualPlots(lm(VegDCA1~Sheep, data = allvar_vege)) # binary, cannot be checked
#residualPlots(lm(VegDCA1~Cow, data = allvar_vege)) # binary, cannot be checked
residualPlots(lm(VegDCA1~FlockSize1_adults, data = allvar_vege)) # Tukey=2.66 ; p=0.0078 -> rejected (one outlier? -> no)
#residualPlots(lm(VegDCA1~GrazingSurface_ha, data = allvar_vege)) # Tukey=0.073 ; p=0.94 -> validated
#residualPlots(lm(VegDCA1~TotalInfieldSurface, data = allvar_vege)) # Tukey=1.58 ; p=0.11 -> validated
#residualPlots(lm(VegDCA1~Grazingdensity_perha, data = allvar_vege)) # Tukey=0.66 ; p=0.51 -> validated
#residualPlots(lm(VegDCA1~MeanBD, data = allvar_vege)) # Tukey=0.34 ; p=0.73 -> validated
residualPlots(lm(VegDCA1~MeanPT, data = allvar_vege)) # Tukey=-1.7 ; p=0.088 -> a bit tight
#residualPlots(lm(VegDCA1~MeanpH, data = allvar_vege)) # Tukey=0.35 ; p=0.73 -> validated
#residualPlots(lm(VegDCA1~MeanPhosphorus, data = allvar_vege)) # Tukey=0.48 ; p=0.63 -> validated

# Residuals LM for beetles
#residualPlots(lm(BeetlePCA1~Elevation_max, data = allvar_beetle)) # Tukey=1.14 ; p=0.25 -> validated
#residualPlots(lm(BeetlePCA1~General_slope, data = allvar_beetle)) # Tukey=-0.57 ; p=0.57 -> validated
#residualPlots(lm(BeetlePCA1~AspectDegree, data = allvar_beetle)) # Tukey=-0.61 ; p=0.54 -> validated
#residualPlots(lm(BeetlePCA1~DistanceToSea_m, data = allvar_beetle)) # Tukey=0.66 ; p=0.51 -> validated
#residualPlots(lm(BeetlePCA1~TotCultivatedLand_percent, data = allvar_beetle)) # Tukey=0.35 ; p=0.73 -> validated
#residualPlots(lm(BeetlePCA1~TotForest_percent, data = allvar_beetle)) # Tukey=0.15 ; p=0.88 -> validated
#residualPlots(lm(BeetlePCA1~Infield_percent, data = allvar_beetle)) # Tukey=0.26 ; p=0.8 -> validated
residualPlots(lm(BeetlePCA1~Outfield_percent, data = allvar_beetle)) # Tukey=-1.94 ; p=0.052 -> a bit tight
#residualPlots(lm(BeetlePCA1~Wetland_percent, data = allvar_beetle)) # Tukey=0.71 ; p=0.48 -> validated
#residualPlots(lm(BeetlePCA1~Sheep, data = allvar_beetle)) # binary, cannot be checked
#residualPlots(lm(BeetlePCA1~Cow, data = allvar_beetle)) # binary, cannot be checked
#residualPlots(lm(BeetlePCA1~FlockSize1_adults, data = allvar_beetle)) # Tukey=0.82 ; p=0.41 -> validated
#residualPlots(lm(BeetlePCA1~GrazingSurface_ha, data = allvar_beetle)) # Tukey=0.87 ; p=0.39 -> validated
#residualPlots(lm(BeetlePCA1~TotalInfieldSurface, data = allvar_beetle)) # Tukey=-0.88 ; p=0.38 -> validated
#residualPlots(lm(BeetlePCA1~Grazingdensity_perha, data = allvar_beetle)) # Tukey=0.31 ; p=0.76 -> validated
#residualPlots(lm(BeetlePCA1~MeanBD, data = allvar_beetle)) # Tukey=-0.6 ; p=0.55 -> validated
residualPlots(lm(BeetlePCA1~MeanExposedGround, data = allvar_beetle)) # Tukey=-2.28 ; p=0.023 -> rejected (one outlier? -> not really after verification)
residualPlots(lm(BeetlePCA1~MeanLitter, data = allvar_beetle)) # Tukey=-2.06 ; p=0.039 -> rejected (one outlier? -> not really after verification)
#residualPlots(lm(BeetlePCA1~MeanBryo, data = allvar_beetle)) # Tukey=-0.66 ; p=0.51 -> validated
#residualPlots(lm(BeetlePCA1~MeanHeight, data = allvar_beetle)) # Tukey=0.21 ; p=0.83 -> validated
#residualPlots(lm(BeetlePCA1~MeanRichness, data = allvar_beetle)) # Tukey=0.34 ; p=0.73 -> validated

# Removal variables with non-linear relationships with residuals
#fjordsys validated for both plants and beetles
landscape_vege <- subset(landscape, select = -c(Wetland_percent)) # linear relationship rejected
landscape_beetle <- subset(landscape, select = -c(Outfield_percent)) # very tight on linear relationship assumption
grazing_vege <- subset(grazing, select = -c(FlockSize1_adults)) # linear relationship rejected
grazing_beetle <- grazing
#locenvi_vege validated
locenvi_beetle <- subset(locenvi_beetle, select = -c(MeanExposedGround, MeanLitter)) # linear relationship rejected


#### Canonical correlation preparation ####

#
## Data preparation - Contingency tables

# Fjord system
fjordsys_long <- fjordsys |> 
  pivot_longer(
    cols = c(Elevation_max, General_slope, AspectDegree, DistanceToSea_m),
    names_to = "Factors",
    values_to = "Values")
contin_fjordsys <- xtabs(formula = Values ~ SiteID + Factors, data = fjordsys_long)

# Landscape matrix - vege
landscape_long_vege <- landscape_vege |> 
  pivot_longer(
    cols = c(TotCultivatedLand_percent, TotForest_percent, Outfield_percent, Infield_percent),
    names_to = "Factors",
    values_to = "Values")
contin_landscape_vege <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long_vege)

# Landscape matrix - beetle
landscape_long_beetle <- landscape_beetle |> 
  pivot_longer(
    cols = c(TotCultivatedLand_percent, TotForest_percent, Infield_percent, Wetland_percent),
    names_to = "Factors",
    values_to = "Values")
contin_landscape_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long_beetle)

# Grazing - vege
grazing_long_vege <- grazing_vege |> 
  pivot_longer(
    cols = c(Sheep, Cow, GrazingSurface_ha, TotalInfieldSurface, Grazingdensity_perha),
    names_to = "Factors",
    values_to = "Values")
contin_grazing_vege <- xtabs(formula = Values ~ SiteID + Factors, data = grazing_long_vege)

# Grazing - beetle
grazing_long_beetle <- grazing_beetle |> 
  pivot_longer(
    cols = c(Sheep, Cow, FlockSize1_adults, GrazingSurface_ha, TotalInfieldSurface, Grazingdensity_perha),
    names_to = "Factors",
    values_to = "Values")
contin_grazing_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = grazing_long_beetle)

# Local environment - vege
locenvi_long_vege <- locenvi_vege |> 
  pivot_longer(
    cols = c(MeanBD, MeanPT, MeanpH, MeanPhosphorus),
    names_to = "Factors",
    values_to = "Values")
contin_locenvi_vege <- xtabs(formula = Values ~ SiteID + Factors, data = locenvi_long_vege)

# Local environment - beetle
locenvi_long_beetle <- locenvi_beetle |> 
  pivot_longer(
    cols = c(MeanHeight, MeanRichness, MeanBD, MeanBryo),
    names_to = "Factors",
    values_to = "Values")
contin_locenvi_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = locenvi_long_beetle)

#
## Graphical representation of correlations plant community

# Fjord effect
correl_fjordxlandscape_vege <- matcor(contin_fjordsys, contin_landscape_vege)
img.matcor(correl_fjordxlandscape_vege, type = 2)
correl_fjordxgrazing_vege <- matcor(contin_fjordsys, contin_grazing_vege)
img.matcor(correl_fjordxgrazing_vege, type = 2)
correl_fjordxlocenvi_vege <- matcor(contin_fjordsys, contin_locenvi_vege)
img.matcor(correl_fjordxlocenvi_vege, type = 2)
correl_fjordxvege <- matcor(contin_fjordsys, contin_vege)
img.matcor(correl_fjordxvege, type = 2)

# Landscape effect
correl_landscapexgrazing_vege <- matcor(contin_landscape_vege, contin_grazing_vege)
img.matcor(correl_landscapexgrazing_vege, type = 2)
correl_landscapexlocalenvi_vege <- matcor(contin_landscape_vege, contin_locenvi_vege)
img.matcor(correl_landscapexlocalenvi_vege, type = 2)
correl_landscapexvege <- matcor(contin_landscape_vege, contin_vege)
img.matcor(correl_landscapexvege, type = 2)

# Grazing effect
correl_grazingxlocenvi_vege <- matcor(contin_grazing_vege, contin_locenvi_vege)
img.matcor(correl_grazingxlocenvi_vege, type = 2)
correl_grazingxvege <- matcor(contin_grazing_vege, contin_vege)
img.matcor(correl_grazingxlocenvi_vege, type = 2)

# Local environment effect
correl_locenvixvege <- matcor(contin_locenvi_vege, contin_vege)
img.matcor(correl_grazingxlocenvi_vege, type = 2)

#
## Graphical representation of correlations beetle community

# Fjord effect
correl_fjordxlandscape_beetle <- matcor(contin_fjordsys, contin_landscape_beetle)
img.matcor(correl_fjordxlandscape_beetle, type = 2)
correl_fjordxgrazing_beetle <- matcor(contin_fjordsys, contin_grazing_beetle)
img.matcor(correl_fjordxgrazing_beetle, type = 2)
correl_fjordxlocenvi_beetle <- matcor(contin_fjordsys, contin_locenvi_beetle)
img.matcor(correl_fjordxlocenvi_beetle, type = 2)
correl_fjordxbeetle <- matcor(contin_fjordsys, contin_beetle)
img.matcor(correl_fjordxbeetle, type = 2)

# Landscape effect
correl_landscapexgrazing_beetle <- matcor(contin_landscape_beetle, contin_grazing_beetle)
img.matcor(correl_landscapexgrazing_beetle, type = 2)
correl_landscapexlocalenvi_beetle <- matcor(contin_landscape_beetle, contin_locenvi_beetle)
img.matcor(correl_landscapexlocalenvi_beetle, type = 2)
correl_landscapexbeetle <- matcor(contin_landscape_beetle, contin_beetle)
img.matcor(correl_landscapexbeetle, type = 2)

# Grazing effect
correl_grazingxlocenvi_beetle <- matcor(contin_grazing_beetle, contin_locenvi_beetle)
img.matcor(correl_grazingxlocenvi_beetle, type = 2)
correl_grazingxbeetle <- matcor(contin_grazing_vege, contin_beetle)
img.matcor(correl_grazingxbeetle, type = 2)

# Local environment effect
correl_locenvixbeetle <- matcor(contin_locenvi_beetle, contin_beetle)
img.matcor(correl_locenvixbeetle, type = 2)


#### Canonical correlation analysis plant community ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape_vege <- length(select_if(landscape_vege, is.numeric))
nvar_grazing_vege <- length(select_if(grazing_vege, is.numeric))
nvar_locenvi_vege <- length(select_if(locenvi_vege, is.numeric))
nvar_vege <- length(select_if(vege, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor_fjordxlandscape_vege <- cc(contin_fjordsys, contin_landscape_vege)
rho_fjordxlandscape_vege <- cancor_fjordxlandscape_vege$cor
rho_fjordxlandscape_vege # 1st axis correlation 0.78
p.asym(rho_fjordxlandscape_vege, nobs, nvar_fjordsys, nvar_landscape_vege, tstat = "Hotelling") # 1st dimension significant pval = 1.46.10-4

# Fjord x grazing
cancor_fjordxgrazing_vege <- cc(contin_fjordsys, contin_grazing_vege)
rho_fjordxgrazing_vege <- cancor_fjordxgrazing_vege$cor
rho_fjordxgrazing_vege # 1st axis correlation 0.53
p.asym(rho_fjordxgrazing_vege, nobs, nvar_fjordsys, nvar_grazing_vege, tstat = "Hotelling") # NS pval = 0.95

# Fjord x local environment
cancor_fjordxlocenvi_vege <- cc(contin_fjordsys, contin_locenvi_vege)
rho_fjordxlocenvi_vege <- cancor_fjordxlocenvi_vege$cor
rho_fjordxlocenvi_vege # 1st axis correlation 0.52
p.asym(rho_fjordxlocenvi_vege, nobs, nvar_fjordsys, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.8

# Fjord x plant community
cancor_fjordxvege <- cc(contin_fjordsys, contin_vege)
rho_fjordxvege <- cancor_fjordxvege$cor
rho_fjordxvege # if 15 species, 1st axis correlation 0.93 ; if 8 main grasses, 1st axis correlation 0.76
p.asym(rho_fjordxvege, nobs, nvar_fjordsys, nvar_vege, tstat = "Hotelling") # if 15 species, NS with pval=0.21 ; if 8 main grasses, NS with pval=0.35

#
## Landscape effect

# Landscape x grazing
cancor_landscapexgrazing_vege <- cc(contin_landscape_vege, contin_grazing_vege)
rho_landscapexgrazing_vege <- cancor_landscapexgrazing_vege$cor
rho_landscapexgrazing_vege # 1st axis correlation 0.55
p.asym(rho_landscapexgrazing_vege, nobs, nvar_landscape_vege, nvar_grazing_vege, tstat = "Hotelling") # NS pval = 0.76

# Landscape x local environment
cancor_landscapexlocalenvi_vege <- cc(contin_landscape_vege, contin_locenvi_vege)
rho_landscapexlocalenvi_vege <- cancor_landscapexlocalenvi_vege$cor
rho_landscapexlocalenvi_vege # 1st axis correlation 0.43
p.asym(rho_landscapexlocalenvi_vege, nobs, nvar_landscape_vege, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.84

# Landscape x plant community
cancor_landscapexvege <- cc(contin_landscape_vege, contin_vege)
rho_landscapexvege <- cancor_landscapexvege$cor
rho_landscapexvege # if 15 species, 1st axis correlation 0.93 ; if 8 main grasses, 1st axis correlation 0.72
p.asym(rho_landscapexvege, nobs, nvar_landscape_vege, nvar_vege, tstat = "Hotelling") # if 15 species, NS with pval=0.14 ; if 8 main grasses, NS with pval=0.72

#
## Grazing effect

# Grazing x local environment
cancor_grazingxlocenvi_vege <- cc(contin_grazing_vege, contin_locenvi_vege)
rho_grazingxlocenvi_vege <- cancor_grazingxlocenvi_vege$cor
rho_grazingxlocenvi_vege # 1st axis correlation 0.58
p.asym(rho_grazingxlocenvi_vege, nobs, nvar_grazing_vege, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.66

cancor_grazingxvege <- cc(contin_grazing_vege, contin_vege)
rho_grazingxvege <- cancor_grazingxvege$cor
rho_grazingxvege # if 15 species, 1st axis correlation 0.94 ; if 8 main grasses, 1st axis correlation 0.89
p.asym(rho_grazingxvege, nobs, nvar_grazing_vege, nvar_vege, tstat = "Hotelling") # if 15 species, 1st axis significant with pval=0.033 ; if 8 main grasses, 1st axis significant with pval=9.87.10-5

#
## Local environment effect

cancor_locenvixvege <- cc(contin_locenvi_vege, contin_vege)
rho_locenvixvege <- cancor_locenvixvege$cor
rho_locenvixvege # if 15 species, 1st axis correlation 0.93 ; if 8 main grasses, 1st axis correlation 0.83
p.asym(rho_locenvixvege, nobs, nvar_locenvi_vege, nvar_vege, tstat = "Hotelling") # if 15 species, 1st axis significant with pval=0.018 ; if 8 main grasses, 1st axis significant with pval=6.52.10-3


#### Canonical correlation analysis beetle assemblage ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape_beetle <- length(select_if(landscape_beetle, is.numeric))
nvar_grazing_beetle <- length(select_if(grazing_beetle, is.numeric))
nvar_locenvi_beetle <- length(select_if(locenvi_beetle, is.numeric))
nvar_beetle <- length(select_if(beetle, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor_fjordxlandscape_beetle <- cc(contin_fjordsys, contin_landscape_beetle)
rho_fjordxlandscape_beetle <- cancor_fjordxlandscape_beetle$cor
rho_fjordxlandscape_beetle # 1st axis correlation 0.73, 2nd axis correlation 0.67
p.asym(rho_fjordxlandscape_beetle, nobs, nvar_fjordsys, nvar_landscape_beetle, tstat = "Hotelling") # 1st axis significant pval 0.002 ; 2nd axis significant pval 0.028

# Fjord x grazing
cancor_fjordxgrazing_beetle <- cc(contin_fjordsys, contin_grazing_beetle)
rho_fjordxgrazing_beetle <- cancor_fjordxgrazing_beetle$cor
rho_fjordxgrazing_beetle # 1st axis correlation 0.54
p.asym(rho_fjordxgrazing_beetle, nobs, nvar_fjordsys, nvar_grazing_beetle, tstat = "Hotelling") # NS pval=0.94

# Fjord x local environment
cancor_fjordxlocenvi_beetle <- cc(contin_fjordsys, contin_locenvi_beetle)
rho_fjordxlocenvi_beetle <- cancor_fjordxlocenvi_beetle$cor
rho_fjordxlocenvi_beetle # 1st axis correlation 0.5
p.asym(rho_fjordxlocenvi_beetle, nobs, nvar_fjordsys, nvar_locenvi_beetle, tstat = "Hotelling") # NS pval=0.65

# Fjord x beetle assemblage
cancor_fjordxbeetle <- cc(contin_fjordsys, contin_beetle)
rho_fjordxbeetle <- cancor_fjordxbeetle$cor
rho_fjordxbeetle # 1st axis correlation 0.68
p.asym(rho_fjordxbeetle, nobs, nvar_fjordsys, nvar_beetle, tstat = "Hotelling") # NS pval=0.92

#
## Landscape effect

# Landscape x grazing
cancor_landscapexgrazing_beetle <- cc(contin_landscape_beetle, contin_grazing_beetle)
rho_landscapexgrazing_beetle <- cancor_landscapexgrazing_beetle$cor
rho_landscapexgrazing_beetle # 1st axis correlation 0.7
p.asym(rho_landscapexgrazing_beetle, nobs, nvar_landscape_beetle, nvar_grazing_beetle, tstat = "Hotelling") # NS pval=0.34

# Landscape x local environment
cancor_landscapexlocalenvi_beetle <- cc(contin_landscape_beetle, contin_locenvi_beetle)
rho_landscapexlocalenvi_beetle <- cancor_landscapexlocalenvi_beetle$cor
rho_landscapexlocalenvi_beetle # 1st axis correlation 0.48
p.asym(rho_landscapexlocalenvi_beetle, nobs, nvar_landscape_beetle, nvar_locenvi_beetle, tstat = "Hotelling") # NS pval=0.74

# Landscape x beetle assemblage
cancor_landscapexbeetle <- cc(contin_landscape_beetle, contin_beetle)
rho_landscapexbeetle <- cancor_landscapexbeetle$cor
rho_landscapexbeetle # 1st axis correlation 0.64
p.asym(rho_landscapexbeetle, nobs, nvar_landscape_beetle, nvar_beetle, tstat = "Hotelling") # NS pval=0.84

#
## Grazing effect

# Grazing x local environment
cancor_grazingxlocenvi_beetle <- cc(contin_grazing_beetle, contin_locenvi_beetle)
rho_grazingxlocenvi_beetle <- cancor_grazingxlocenvi_beetle$cor
rho_grazingxlocenvi_beetle # 1st axis correlation 0.77
p.asym(rho_grazingxlocenvi_beetle, nobs, nvar_grazing_beetle, nvar_locenvi_beetle, tstat = "Hotelling") # 1st axis significant with pval=0.048

# Grazing x beetle community
cancor_grazingxbeetle <- cc(contin_grazing_beetle, contin_beetle)
rho_grazingxbeetle <- cancor_grazingxbeetle$cor
rho_grazingxbeetle # 1st axis correlation 0.75
p.asym(rho_grazingxbeetle, nobs, nvar_grazing_beetle, nvar_beetle, tstat = "Hotelling") # NS pval=0.68

#
## Local environment effect

cancor_locenvixbeetle <- cc(contin_locenvi_beetle, contin_beetle)
rho_locenvixbeetle <- cancor_locenvixbeetle$cor
rho_locenvixbeetle # 1st axis correlation 0.76
p.asym(rho_locenvixbeetle, nobs, nvar_locenvi_beetle, nvar_beetle, tstat = "Hotelling") # NS pval=0.12