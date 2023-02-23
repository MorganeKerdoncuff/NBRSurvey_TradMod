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
#library(candisc) # visual representation canonical correlation


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
grass <- subset(vege_grass,
                Species == "Agrostis capillaris" |
                Species == "Festuca rubra" | 
                Species == "Holcus lanatus" | 
                Species == "Poa pratensis" | 
                Species == "Deschampsia cespitosa" | 
                Species == "Anthoxantum odoratum" | 
                Species == "Lolium perenne" | 
                Species == "Deschampsia flexuosa" |
                Species == "Poa trivialis")

# Selection main forb species
forb <- subset(vege_grass,
                Species == "Trifolium repens" |
                  Species == "Rumex acetosa" |
                  Species == "Galium saxatile" |
                  Species == "Potentilla erecta" |
                  Species == "Ranunculus repens" |
                  Species == "Achillea millefolium")
  
# Log transformation
grass <- grass |> 
  mutate(PlantSp_logcover = log1p(PlantSp_cover))
forb <- forb |> 
  mutate(PlantSp_logcover = log1p(PlantSp_cover))

# Removal IG3 outlier
grass <- subset(grass, !SiteID == "IG3")
forb <- subset(forb, !SiteID == "IG3")

# Contingency tables
contin_grass <- xtabs(formula = PlantSp_logcover ~ SiteID + Species, data = grass)
contin_forb <- xtabs(formula = PlantSp_logcover ~ SiteID + Species, data = forb)

# Wide table grass
grass <- subset(grass, select = -c(PlantSp_cover))
grass <- grass |> 
  pivot_wider(names_from = Species, values_from = PlantSp_logcover)
grass <- as.data.frame(grass)

# Wide table forbs
forb <- subset(forb, select = -c(PlantSp_cover))
forb <- forb |> 
  pivot_wider(names_from = Species, values_from = PlantSp_logcover)
forb <- as.data.frame(forb)

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

# Removal IG3 outlier
arthro_grass <- subset(arthro_grass, !SiteID == "IG3")

# Contingency table and wide table
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

# Removal IG3 outlier
fjordsys <- subset(fjordsys, !SiteID == "IG3")

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

# Removal IG3 outlier
landscape <- subset(landscape, !SiteID == "IG3")

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

# Removal IG3 outlier
grazing <- subset(grazing, !SiteID == "IG3")

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

# Removal IG3 outlier
locenvi_vege <- subset(locenvi_vege, !SiteID == "IG3")

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

# Removal IG3 outlier
locenvi_beetle <- subset(locenvi_beetle, !SiteID == "IG3")


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
DCA_grass <- decorana(contin_grass)
DCA_grass # Axis length of 2.3 -> keep DCA
DCA_forb <- decorana(contin_forb)
DCA_forb # Axis length of 3.4 -> keep DCA
DCA_beetle <- decorana(contin_beetle)
DCA_beetle # Axis length under 1 -> run PCA
PCA_beetle <- prcomp(contin_beetle)

# Extraction grass species scores from DCA
scores_grass <- as.data.frame(scores(DCA_grass, choices = c(1,2), display = "sites"))
names(scores_grass) <- gsub("DCA1", "GrassDCA1", names(scores_grass))
names(scores_grass) <- gsub("DCA2", "GrassDCA2", names(scores_grass))
scores_grass <- scores_grass |>  
  mutate(SiteID = row.names(scores_grass)) # PlotID as column for future binding

# Extraction forb species scores from DCA (sup)
scores_forb <- as.data.frame(scores(DCA_forb, choices = c(1,2), display = "sites"))
names(scores_forb) <- gsub("DCA1", "ForbDCA1", names(scores_forb))
names(scores_forb) <- gsub("DCA2", "ForbDCA2", names(scores_forb))
scores_forb <- scores_forb |>  
  mutate(SiteID = row.names(scores_forb)) # PlotID as column for future binding

# Extraction beetle family scores from PCA
scores_beetle <- as.data.frame(scores(PCA_beetle, choices = c(1,2), display = "sites"))
names(scores_beetle) <- gsub("PC1", "BeetlePCA1", names(scores_beetle))
names(scores_beetle) <- gsub("PC2", "BeetlePCA2", names(scores_beetle))
scores_beetle <- scores_beetle |>  
  mutate(SiteID = row.names(scores_beetle)) # PlotID as column for future binding

# Table binding
allvar_grass <- purrr::reduce(list(scores_grass, fjordsys, landscape, grazing, locenvi_vege), dplyr::left_join)
allvar_forb <- purrr::reduce(list(scores_forb, fjordsys, landscape, grazing, locenvi_vege), dplyr::left_join)
allvar_beetle <- purrr::reduce(list(scores_beetle, fjordsys, landscape, grazing, locenvi_beetle), dplyr::left_join)

# Residuals LM for grasses
#residualPlots(lm(GrassDCA1~Elevation_max, data = allvar_grass)) # Tukey=0.98 ; p-value=0.33 -> validated
#residualPlots(lm(GrassDCA1~General_slope, data = allvar_grass)) # Tukey=-0.79 ; p-value=0.43 -> validated
#residualPlots(lm(GrassDCA1~AspectDegree, data = allvar_grass)) # Tukey=0.5 ; p-value=0.62 -> validated
#residualPlots(lm(GrassDCA1~DistanceToSea_m, data = allvar_grass)) # Tukey=0.34 ; p=0.74 -> validated
#residualPlots(lm(GrassDCA1~TotCultivatedLand_percent, data = allvar_grass)) # Tukey=0.3 ; p=0.76 -> validated
#residualPlots(lm(GrassDCA1~TotForest_percent, data = allvar_grass)) # Tukey=-0.55 ; p=0.58 -> validated
#residualPlots(lm(GrassDCA1~Infield_percent, data = allvar_grass)) # Tukey=-0.45 ; p=0.65 -> validated
#residualPlots(lm(GrassDCA1~Outfield_percent, data = allvar_grass)) # Tukey=0.71 ; p=0.48 -> validated
residualPlots(lm(GrassDCA1~Wetland_percent, data = allvar_grass)) # Tukey=2.05 ; p=0.04 -> rejected
#residualPlots(lm(GrassDCA1~Sheep, data = allvar_grass)) # binary, cannot be checked
#residualPlots(lm(GrassDCA1~Cow, data = allvar_grass)) # binary, cannot be checked
#residualPlots(lm(GrassDCA1~Goat, data = allvar_grass)) # binary, cannot be checked
residualPlots(lm(GrassDCA1~FlockSize1_adults, data = allvar_grass)) # Tukey=2.66 ; p=0.0078 -> rejected (one outlier? -> no)
#residualPlots(lm(GrassDCA1~GrazingSurface_ha, data = allvar_grass)) # Tukey=0.073 ; p=0.94 -> validated
#residualPlots(lm(GrassDCA1~TotalInfieldSurface, data = allvar_grass)) # Tukey=1.58 ; p=0.11 -> validated
#residualPlots(lm(GrassDCA1~Grazingdensity_perha, data = allvar_grass)) # Tukey=0.66 ; p=0.51 -> validated
#residualPlots(lm(GrassDCA1~MeanBD, data = allvar_grass)) # Tukey=0.34 ; p=0.73 -> validated
#residualPlots(lm(GrassDCA1~MeanPT, data = allvar_grass)) # Tukey=-1.7 ; p=0.088 -> a bit tight
#residualPlots(lm(GrassDCA1~MeanpH, data = allvar_grass)) # Tukey=0.35 ; p=0.73 -> validated
#residualPlots(lm(GrassDCA1~MeanPhosphorus, data = allvar_grass)) # Tukey=0.48 ; p=0.63 -> validated

# Residuals LM for forbs (sup)
#residualPlots(lm(ForbDCA1~Elevation_max, data = allvar_forb)) # Tukey=.0.68 ; p-value=0.5 -> validated
#residualPlots(lm(ForbDCA1~General_slope, data = allvar_forb)) # Tukey=-0.63 ; p-value=0.53 -> validated
residualPlots(lm(ForbDCA1~AspectDegree, data = allvar_forb)) # Tukey=-2.46 ; p-value=0.0.014 -> rejected
#residualPlots(lm(ForbDCA1~DistanceToSea_m, data = allvar_forb)) # Tukey=1.65 ; p=0.1 -> a bit tight
#residualPlots(lm(ForbDCA1~TotCultivatedLand_percent, data = allvar_forb)) # Tukey=0.051 ; p=0.96 -> validated
#residualPlots(lm(ForbDCA1~TotForest_percent, data = allvar_forb)) # Tukey=0.51 ; p=0.61 -> validated
#residualPlots(lm(ForbDCA1~Infield_percent, data = allvar_forb)) # Tukey=1.16 ; p=0.25 -> validated
#residualPlots(lm(ForbDCA1~Outfield_percent, data = allvar_forb)) # Tukey=0.32 ; p=0.75 -> validated
#residualPlots(lm(ForbDCA1~Wetland_percent, data = allvar_forb)) # Tukey=0.4 ; p=0.69 -> validated
#residualPlots(lm(ForbDCA1~Sheep, data = allvar_forb)) # binary, cannot be checked
#residualPlots(lm(ForbDCA1~Cow, data = allvar_forb)) # binary, cannot be checked
#residualPlots(lm(ForbDCA1~Goat, data = allvar_forb)) # binary, cannot be checked
residualPlots(lm(ForbDCA1~FlockSize1_adults, data = allvar_forb)) # Tukey=-2.38 ; p=0.017 -> rejected (one outlier? -> no)
#residualPlots(lm(ForbDCA1~GrazingSurface_ha, data = allvar_forb)) # Tukey=-1.34 ; p=0.18 -> validated
#residualPlots(lm(ForbDCA1~TotalInfieldSurface, data = allvar_forb)) # Tukey=-1.16 ; p=0.25 -> validated
#residualPlots(lm(ForbDCA1~Grazingdensity_perha, data = allvar_forb)) # Tukey=-1.06 ; p=0.29 -> validated
#residualPlots(lm(ForbDCA1~MeanBD, data = allvar_forb)) # Tukey=-0.33 ; p=0.74 -> validated
#residualPlots(lm(ForbDCA1~MeanPT, data = allvar_forb)) # Tukey=1.23 ; p=0.22 -> validated
#residualPlots(lm(ForbDCA1~MeanpH, data = allvar_forb)) # Tukey=0.12 ; p=0.9 -> validated
#residualPlots(lm(ForbDCA1~MeanPhosphorus, data = allvar_forb)) # Tukey=1.18 ; p=0.07 -> a bit tight

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
#residualPlots(lm(BeetlePCA1~Goat, data = allvar_beetle)) # binary, cannot be checked
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
fjordsys_forb <- subset(fjordsys, select = -c(AspectDegree))
landscape_grass <- subset(landscape, select = -c(Wetland_percent)) # linear relationship rejected
landscape_beetle <- subset(landscape, select = -c(Outfield_percent)) # very tight on linear relationship assumption
grazing_vege <- subset(grazing, select = -c(FlockSize1_adults)) # linear relationship rejected for both grass & forb
#locenvi vege validated for both grass & forb
locenvi_beetle <- subset(locenvi_beetle, select = -c(MeanExposedGround, MeanLitter)) # linear relationship rejected


#### Canonical correlation preparation ####

#
## Data preparation - Contingency tables

# Fjord system - grass & beetles
fjordsys_long <- fjordsys |> 
  pivot_longer(
    cols = c(Elevation_max, General_slope, AspectDegree, DistanceToSea_m),
    names_to = "Factors",
    values_to = "Values")
contin_fjordsys <- xtabs(formula = Values ~ SiteID + Factors, data = fjordsys_long)

# Fjord system - forbs
fjordsys_long_forb <- fjordsys_forb |> 
  pivot_longer(
    cols = c(Elevation_max, General_slope, DistanceToSea_m),
    names_to = "Factors",
    values_to = "Values")
contin_fjordsys_forb <- xtabs(formula = Values ~ SiteID + Factors, data = fjordsys_long_forb)

# Landscape matrix - grass
landscape_long_grass <- landscape_grass |> 
  pivot_longer(
    cols = c(TotCultivatedLand_percent, TotForest_percent, Outfield_percent, Infield_percent),
    names_to = "Factors",
    values_to = "Values")
contin_landscape_grass <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long_grass)

# Landscape matrix - forbs
landscape_long_forb <- landscape |> 
  pivot_longer(
    cols = c(TotCultivatedLand_percent, TotForest_percent, Outfield_percent, Infield_percent, Wetland_percent),
    names_to = "Factors",
    values_to = "Values")
contin_landscape_forb <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long_forb)

# Landscape matrix - beetle
landscape_long_beetle <- landscape_beetle |> 
  pivot_longer(
    cols = c(TotCultivatedLand_percent, TotForest_percent, Infield_percent, Wetland_percent),
    names_to = "Factors",
    values_to = "Values")
contin_landscape_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long_beetle)

# Grazing - vege (grass + forbs)
grazing_long_vege <- grazing_vege |> 
  pivot_longer(
    cols = c(Sheep, Cow, GrazingSurface_ha, TotalInfieldSurface, Grazingdensity_perha),
    names_to = "Factors",
    values_to = "Values")
contin_grazing_vege <- xtabs(formula = Values ~ SiteID + Factors, data = grazing_long_vege)

# Grazing - beetle
grazing_long_beetle <- grazing |> 
  pivot_longer(
    cols = c(Sheep, Cow, FlockSize1_adults, GrazingSurface_ha, TotalInfieldSurface, Grazingdensity_perha),
    names_to = "Factors",
    values_to = "Values")
contin_grazing_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = grazing_long_beetle)

# Local environment - vege (grass + forbs)
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
#correl_fjordxlandscape_vege <- matcor(contin_fjordsys, contin_landscape_vege)
#img.matcor(correl_fjordxlandscape_vege, type = 2)
#correl_fjordxgrazing_vege <- matcor(contin_fjordsys, contin_grazing_vege)
#img.matcor(correl_fjordxgrazing_vege, type = 2)
#correl_fjordxlocenvi_vege <- matcor(contin_fjordsys, contin_locenvi_vege)
#img.matcor(correl_fjordxlocenvi_vege, type = 2)
#correl_fjordxvege <- matcor(contin_fjordsys, contin_vege)
#img.matcor(correl_fjordxvege, type = 2)

# Landscape effect
#correl_landscapexgrazing_vege <- matcor(contin_landscape_vege, contin_grazing_vege)
#img.matcor(correl_landscapexgrazing_vege, type = 2)
#correl_landscapexlocalenvi_vege <- matcor(contin_landscape_vege, contin_locenvi_vege)
#img.matcor(correl_landscapexlocalenvi_vege, type = 2)
#correl_landscapexvege <- matcor(contin_landscape_vege, contin_vege)
#img.matcor(correl_landscapexvege, type = 2)

# Grazing effect
#correl_grazingxlocenvi_vege <- matcor(contin_grazing_vege, contin_locenvi_vege)
#img.matcor(correl_grazingxlocenvi_vege, type = 2)
#correl_grazingxvege <- matcor(contin_grazing_vege, contin_vege)
#img.matcor(correl_grazingxlocenvi_vege, type = 2)

# Local environment effect
#correl_locenvixvege <- matcor(contin_locenvi_vege, contin_vege)
#img.matcor(correl_grazingxlocenvi_vege, type = 2)

#
## Graphical representation of correlations beetle community

# Fjord effect
#correl_fjordxlandscape_beetle <- matcor(contin_fjordsys, contin_landscape_beetle)
#img.matcor(correl_fjordxlandscape_beetle, type = 2)
#correl_fjordxgrazing_beetle <- matcor(contin_fjordsys, contin_grazing_beetle)
#img.matcor(correl_fjordxgrazing_beetle, type = 2)
#correl_fjordxlocenvi_beetle <- matcor(contin_fjordsys, contin_locenvi_beetle)
#img.matcor(correl_fjordxlocenvi_beetle, type = 2)
#correl_fjordxbeetle <- matcor(contin_fjordsys, contin_beetle)
#img.matcor(correl_fjordxbeetle, type = 2)

# Landscape effect
#correl_landscapexgrazing_beetle <- matcor(contin_landscape_beetle, contin_grazing_beetle)
#img.matcor(correl_landscapexgrazing_beetle, type = 2)
#correl_landscapexlocalenvi_beetle <- matcor(contin_landscape_beetle, contin_locenvi_beetle)
#img.matcor(correl_landscapexlocalenvi_beetle, type = 2)
#correl_landscapexbeetle <- matcor(contin_landscape_beetle, contin_beetle)
#img.matcor(correl_landscapexbeetle, type = 2)

# Grazing effect
#correl_grazingxlocenvi_beetle <- matcor(contin_grazing_beetle, contin_locenvi_beetle)
#img.matcor(correl_grazingxlocenvi_beetle, type = 2)
#correl_grazingxbeetle <- matcor(contin_grazing_vege, contin_beetle)
#img.matcor(correl_grazingxbeetle, type = 2)

# Local environment effect
#correl_locenvixbeetle <- matcor(contin_locenvi_beetle, contin_beetle)
#img.matcor(correl_locenvixbeetle, type = 2)


#### Canonical correlation analysis grass community ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape_grass <- length(select_if(landscape_grass, is.numeric))
nvar_grazing_vege <- length(select_if(grazing_vege, is.numeric))
nvar_locenvi_vege <- length(select_if(locenvi_vege, is.numeric))
nvar_grass <- length(select_if(grass, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor_fjordxlandscape_grass <- cc(contin_fjordsys, contin_landscape_grass)
rho_fjordxlandscape_grass <- cancor_fjordxlandscape_grass$cor
rho_fjordxlandscape_grass # 1st axis correlation 0.83 ; 2nd axis correlation 0.7
p.asym(rho_fjordxlandscape_grass, nobs, nvar_fjordsys, nvar_landscape_grass, tstat = "Hotelling") # 1st dim significant pval = 2.26.10-5 ; 2nd dim significant pval = 1.67.10-2
plt.cc(cancor_fjordxlandscape_grass, var.label = TRUE)

# Fjord x grazing
cancor_fjordxgrazing_grass <- cc(contin_fjordsys, contin_grazing_vege)
rho_fjordxgrazing_grass <- cancor_fjordxgrazing_grass$cor
rho_fjordxgrazing_grass # 1st axis correlation 0.53
p.asym(rho_fjordxgrazing_grass, nobs, nvar_fjordsys, nvar_grazing_vege, tstat = "Hotelling") # NS pval = 0.95
plt.cc(cancor_fjordxgrazing_grass, var.label = TRUE)

# Fjord x local environment
cancor_fjordxlocenvi_grass <- cc(contin_fjordsys, contin_locenvi_vege)
rho_fjordxlocenvi_grass <- cancor_fjordxlocenvi_grass$cor
rho_fjordxlocenvi_grass # 1st axis correlation 0.52
p.asym(rho_fjordxlocenvi_grass, nobs, nvar_fjordsys, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.72
plt.cc(cancor_fjordxlocenvi_grass, var.label = TRUE)

# Fjord x grass community
cancor_fjordxgrass <- cc(contin_fjordsys, contin_grass)
rho_fjordxgrass <- cancor_fjordxgrass$cor
rho_fjordxgrass # 1st axis correlation 0.76
p.asym(rho_fjordxgrass, nobs, nvar_fjordsys, nvar_grass, tstat = "Hotelling") # NS with pval=0.54
plt.cc(cancor_fjordxgrass, var.label = TRUE)

#
## Landscape effect

# Landscape x grazing
cancor_landscapexgrazing_grass <- cc(contin_landscape_grass, contin_grazing_vege)
rho_landscapexgrazing_grass <- cancor_landscapexgrazing_grass$cor
rho_landscapexgrazing_grass # 1st axis correlation 0.6
p.asym(rho_landscapexgrazing_grass, nobs, nvar_landscape_grass, nvar_grazing_vege, tstat = "Hotelling") # NS pval = 0.44
plt.cc(cancor_landscapexgrazing_grass, var.label = TRUE)

# Landscape x local environment
cancor_landscapexlocalenvi_grass <- cc(contin_landscape_grass, contin_locenvi_vege)
rho_landscapexlocalenvi_grass <- cancor_landscapexlocalenvi_grass$cor
rho_landscapexlocalenvi_grass # 1st axis correlation 0.45
p.asym(rho_landscapexlocalenvi_grass, nobs, nvar_landscape_grass, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.82
plt.cc(cancor_landscapexlocalenvi_grass, var.label = TRUE)

# Landscape x grass community
cancor_landscapexgrass <- cc(contin_landscape_grass, contin_grass)
rho_landscapexgrass <- cancor_landscapexgrass$cor
rho_landscapexgrass # 1st axis correlation 0.71
p.asym(rho_landscapexgrass, nobs, nvar_landscape_grass, nvar_grass, tstat = "Hotelling") # NS with pval=0.77
plt.cc(cancor_landscapexgrass, var.label = TRUE)

#
## Grazing effect

# Grazing x local environment
cancor_grazingxlocenvi_grass <- cc(contin_grazing_vege, contin_locenvi_vege)
rho_grazingxlocenvi_grass <- cancor_grazingxlocenvi_grass$cor
rho_grazingxlocenvi_grass # 1st axis correlation 0.59
p.asym(rho_grazingxlocenvi_grass, nobs, nvar_grazing_vege, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.71
plt.cc(cancor_grazingxlocenvi_grass, var.label = TRUE)

# Grazing x grass community
cancor_grazingxgrass <- cc(contin_grazing_vege, contin_grass)
rho_grazingxgrass <- cancor_grazingxgrass$cor
rho_grazingxgrass # 1st axis correlation 0.89 : 2nd axis correlation 0.87
p.asym(rho_grazingxgrass, nobs, nvar_grazing_vege, nvar_grass, tstat = "Hotelling") # 1st axis significant with pval = 1.44.10-4 ; 2nd axis significant with pval = 5.34.10-3
plt.cc(cancor_grazingxgrass, var.label = TRUE)

#
## Local environment effect

cancor_locenvixgrass <- cc(contin_locenvi_vege, contin_grass)
rho_locenvixgrass <- cancor_locenvixgrass$cor
rho_locenvixgrass # 1st axis correlation 0.83 ; 2nd axis correlation 0.77
p.asym(rho_locenvixgrass, nobs, nvar_locenvi_vege, nvar_grass, tstat = "Hotelling") # 1st axis significant with pval=6.52.10-3 ; 2nd axis almost significant with pval = 0.06
plt.cc(cancor_locenvixgrass, var.label = TRUE)


#### Canonical correlation analysis forb community ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

# Number of variables in each set
nvar_fjordsys_forb <- length(select_if(fjordsys_forb, is.numeric))
nvar_landscape_forb <- length(select_if(landscape, is.numeric))
nvar_grazing_vege <- length(select_if(grazing_vege, is.numeric))
nvar_locenvi_vege <- length(select_if(locenvi_vege, is.numeric))
nvar_forb <- length(select_if(forb, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor_fjordxlandscape_forb <- cc(contin_fjordsys_forb, contin_landscape_forb)
rho_fjordxlandscape_forb <- cancor_fjordxlandscape_forb$cor
rho_fjordxlandscape_forb # 1st axis correlation 0.8 ; 2nd axis correlation 0.71 ; 3rd axis correlation 0.59
p.asym(rho_fjordxlandscape_forb, nobs, nvar_fjordsys_forb, nvar_landscape_forb, tstat = "Hotelling") # 1st dim significant pval = 2.09.10-5 ; 2nd dim significant pval = 3.8.10-4 ; 3rd dim significant pval = 9.17.10-3
plt.cc(cancor_fjordxlandscape_forb, var.label = TRUE) # dim 1 & 2
plt.cc(cancor_fjordxlandscape_forb, d1 = 1, d2 = 3, var.label = TRUE) # dim 1 & 3

# Fjord x grazing
cancor_fjordxgrazing_forb <- cc(contin_fjordsys_forb, contin_grazing_vege)
rho_fjordxgrazing_forb <- cancor_fjordxgrazing_forb$cor
rho_fjordxgrazing_forb # 1st axis correlation 0.53
p.asym(rho_fjordxgrazing_forb, nobs, nvar_fjordsys_forb, nvar_grazing_vege, tstat = "Hotelling") # NS pval = 0.84
plt.cc(cancor_fjordxgrazing_forb, var.label = TRUE)

# Fjord x local environment
cancor_fjordxlocenvi_forb <- cc(contin_fjordsys_forb, contin_locenvi_vege)
rho_fjordxlocenvi_forb <- cancor_fjordxlocenvi_forb$cor
rho_fjordxlocenvi_forb # 1st axis correlation 0.45
p.asym(rho_fjordxlocenvi_forb, nobs, nvar_fjordsys_forb, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.86
plt.cc(cancor_fjordxlocenvi_forb, var.label = TRUE)

# Fjord x forb community
cancor_fjordxforb <- cc(contin_fjordsys_forb, contin_forb)
rho_fjordxforb <- cancor_fjordxforb$cor
rho_fjordxforb # 1st axis correlation 0.67
p.asym(rho_fjordxforb, nobs, nvar_fjordsys_forb, nvar_forb, tstat = "Hotelling") # NS with pval = 0.4
plt.cc(cancor_fjordxforb, var.label = TRUE)

#
## Landscape effect

# Landscape x grazing
cancor_landscapexgrazing_forb <- cc(contin_landscape_forb, contin_grazing_vege)
rho_landscapexgrazing_forb <- cancor_landscapexgrazing_forb$cor
rho_landscapexgrazing_forb # 1st axis correlation 0.67
p.asym(rho_landscapexgrazing_forb, nobs, nvar_landscape_forb, nvar_grazing_vege, tstat = "Hotelling") # NS pval = 0.44
plt.cc(cancor_landscapexgrazing_forb, var.label = TRUE)

# Landscape x local environment
cancor_landscapexlocalenvi_forb <- cc(contin_landscape_forb, contin_locenvi_vege)
rho_landscapexlocalenvi_forb <- cancor_landscapexlocalenvi_forb$cor
rho_landscapexlocalenvi_forb # 1st axis correlation 0.63
p.asym(rho_landscapexlocalenvi_forb, nobs, nvar_landscape_forb, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.49
plt.cc(cancor_landscapexlocalenvi_forb, var.label = TRUE)

# Landscape x forb community
cancor_landscapexforb <- cc(contin_landscape_forb, contin_forb)
rho_landscapexforb <- cancor_landscapexforb$cor
rho_landscapexforb # 1st axis correlation 0.7
p.asym(rho_landscapexforb, nobs, nvar_landscape_forb, nvar_forb, tstat = "Hotelling") # NS with pval=0.37
plt.cc(cancor_landscapexforb, var.label = TRUE)

#
## Grazing effect

# Grazing x local environment
cancor_grazingxlocenvi_forb <- cc(contin_grazing_vege, contin_locenvi_vege)
rho_grazingxlocenvi_forb <- cancor_grazingxlocenvi_forb$cor
rho_grazingxlocenvi_forb # 1st axis correlation 0.59
p.asym(rho_grazingxlocenvi_forb, nobs, nvar_grazing_vege, nvar_locenvi_vege, tstat = "Hotelling") # NS pval = 0.71
plt.cc(cancor_grazingxlocenvi_forb, var.label = TRUE)

# Grazing x forb community
cancor_grazingxforb <- cc(contin_grazing_vege, contin_forb)
rho_grazingxforb <- cancor_grazingxforb$cor
rho_grazingxforb # 1st axis correlation 0.76
p.asym(rho_grazingxforb, nobs, nvar_grazing_vege, nvar_forb, tstat = "Hotelling") # NS with pval=0.18
plt.cc(cancor_grazingxforb, var.label = TRUE)

#
## Local environment effect

cancor_locenvixforb <- cc(contin_locenvi_vege, contin_forb)
rho_locenvixforb <- cancor_locenvixforb$cor
rho_locenvixforb # 1st axis correlation 0.69
p.asym(rho_locenvixforb, nobs, nvar_locenvi_vege, nvar_forb, tstat = "Hotelling") # 1st axis significant with pval=0.028
plt.cc(cancor_locenvixforb, var.label = TRUE)


#### Canonical correlation analysis beetle assemblage ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape_beetle <- length(select_if(landscape_beetle, is.numeric))
nvar_grazing_beetle <- length(select_if(grazing, is.numeric))
nvar_locenvi_beetle <- length(select_if(locenvi_beetle, is.numeric))
nvar_beetle <- length(select_if(beetle, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor_fjordxlandscape_beetle <- cc(contin_fjordsys, contin_landscape_beetle)
rho_fjordxlandscape_beetle <- cancor_fjordxlandscape_beetle$cor
rho_fjordxlandscape_beetle # 1st axis correlation 0.76, 2nd axis 0.72, 3rd axis 0.58
p.asym(rho_fjordxlandscape_beetle, nobs, nvar_fjordsys, nvar_landscape_beetle, tstat = "Hotelling") # 1st axis significant pval 5.63.10-5 ; 2nd axis significant pval 4.10-4 ; 3rd axis significant pval 2.25.10-2
plt.cc(cancor_fjordxlandscape_beetle, var.label = TRUE) # dim 1 & 2
plt.cc(cancor_fjordxlandscape_beetle, d1 = 1, d2 = 3, var.label = TRUE) # dim 1 & 3

# Fjord x grazing
cancor_fjordxgrazing_beetle <- cc(contin_fjordsys, contin_grazing_beetle)
rho_fjordxgrazing_beetle <- cancor_fjordxgrazing_beetle$cor
rho_fjordxgrazing_beetle # 1st axis correlation 0.54
p.asym(rho_fjordxgrazing_beetle, nobs, nvar_fjordsys, nvar_grazing_beetle, tstat = "Hotelling") # NS pval=0.94
plt.cc(cancor_fjordxgrazing_beetle, var.label = TRUE)

# Fjord x local environment
cancor_fjordxlocenvi_beetle <- cc(contin_fjordsys, contin_locenvi_beetle)
rho_fjordxlocenvi_beetle <- cancor_fjordxlocenvi_beetle$cor
rho_fjordxlocenvi_beetle # 1st axis correlation 0.48
p.asym(rho_fjordxlocenvi_beetle, nobs, nvar_fjordsys, nvar_locenvi_beetle, tstat = "Hotelling") # NS pval=0.87
plt.cc(cancor_fjordxlocenvi_beetle, var.label = TRUE)

# Fjord x beetle assemblage
cancor_fjordxbeetle <- cc(contin_fjordsys, contin_beetle)
rho_fjordxbeetle <- cancor_fjordxbeetle$cor
rho_fjordxbeetle # 1st axis correlation 0.62
p.asym(rho_fjordxbeetle, nobs, nvar_fjordsys, nvar_beetle, tstat = "Hotelling") # NS pval=0.99
plt.cc(cancor_fjordxbeetle, var.label = TRUE)

#
## Landscape effect

# Landscape x grazing
cancor_landscapexgrazing_beetle <- cc(contin_landscape_beetle, contin_grazing_beetle)
rho_landscapexgrazing_beetle <- cancor_landscapexgrazing_beetle$cor
rho_landscapexgrazing_beetle # 1st axis correlation 0.69
p.asym(rho_landscapexgrazing_beetle, nobs, nvar_landscape_beetle, nvar_grazing_beetle, tstat = "Hotelling") # NS pval=0.44
plt.cc(cancor_landscapexgrazing_beetle, var.label = TRUE)

# Landscape x local environment
cancor_landscapexlocalenvi_beetle <- cc(contin_landscape_beetle, contin_locenvi_beetle)
rho_landscapexlocalenvi_beetle <- cancor_landscapexlocalenvi_beetle$cor
rho_landscapexlocalenvi_beetle # 1st axis correlation 0.49
p.asym(rho_landscapexlocalenvi_beetle, nobs, nvar_landscape_beetle, nvar_locenvi_beetle, tstat = "Hotelling") # NS pval=0.79
plt.cc(cancor_landscapexlocalenvi_beetle, var.label = TRUE)

# Landscape x beetle assemblage
cancor_landscapexbeetle <- cc(contin_landscape_beetle, contin_beetle)
rho_landscapexbeetle <- cancor_landscapexbeetle$cor
rho_landscapexbeetle # 1st axis correlation 0.65
p.asym(rho_landscapexbeetle, nobs, nvar_landscape_beetle, nvar_beetle, tstat = "Hotelling") # NS pval=0.87
plt.cc(cancor_landscapexbeetle, var.label = TRUE)

#
## Grazing effect

# Grazing x local environment
cancor_grazingxlocenvi_beetle <- cc(contin_grazing_beetle, contin_locenvi_beetle)
rho_grazingxlocenvi_beetle <- cancor_grazingxlocenvi_beetle$cor
rho_grazingxlocenvi_beetle # 1st axis correlation 0.77
p.asym(rho_grazingxlocenvi_beetle, nobs, nvar_grazing_beetle, nvar_locenvi_beetle, tstat = "Hotelling") # NS with pval=0.094
plt.cc(cancor_grazingxlocenvi_beetle, var.label = TRUE)

# Grazing x beetle community
cancor_grazingxbeetle <- cc(contin_grazing_beetle, contin_beetle)
rho_grazingxbeetle <- cancor_grazingxbeetle$cor
rho_grazingxbeetle # 1st axis correlation 0.76
p.asym(rho_grazingxbeetle, nobs, nvar_grazing_beetle, nvar_beetle, tstat = "Hotelling") # NS pval=0.51
plt.cc(cancor_grazingxbeetle, var.label = TRUE)

#
## Local environment effect

cancor_locenvixbeetle <- cc(contin_locenvi_beetle, contin_beetle)
VIScancor_locenvixbeetle <- cancor(contin_locenvi_beetle, contin_beetle)
rho_locenvixbeetle <- cancor_locenvixbeetle$cor
rho_locenvixbeetle # 1st axis correlation 0.76
p.asym(rho_locenvixbeetle, nobs, nvar_locenvi_beetle, nvar_beetle, tstat = "Hotelling") # NS pval=0.18
plt.cc(cancor_locenvixbeetle, var.label = TRUE)
#VIScancor_locenvixbeetle <- cancor(contin_locenvi_beetle, contin_beetle)
#plot(cancor_locenvixbeetle, which = 1)