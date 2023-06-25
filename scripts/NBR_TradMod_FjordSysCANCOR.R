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
library(ggord) # Ordination plot with ggplot2
library(ggpubr) # Function ggarrange for several plots on same file
library(GGally) # Extension ggplot
library(car) # Visualisation lm residuals plots
library(CCA) # Canonical correlation
library(CCP) # Statistical test canonical correlation
#library(candisc) # visual representation canonical correlation


#### DATA LOADING ####

siteinfo_full <- read.csv("data/cleandata/NBR_FullSiteInfo.csv", sep=",") # Clean site info data
climate_full <- read.csv("data/cleandata/NBR_FullClimate.csv", sep=",") # Clean climate data
landscape_full <- read.csv("data/cleandata/NBR_FullLandscape.csv", sep=",") # Clean landscape matrix data
landuse_full <- read.csv("data/cleandata/NBR_FullLanduse.csv", sep=",") # Clean field management data
area20x20_full <- read.csv("data/cleandata/NBR_FullArea20x20.csv", sep=",") # Clean sampling area 20mx20m data
groundcover_full <- read.csv("data/cleandata/NBR_FullGroundCover.csv", sep=",") # Clean aboveground cover data
soilbulk_full <- read.csv("data/cleandata/NBR_FullSoilBulk.csv", sep=",") # Clean soil bulk density data
soilchem_full <- read.csv("data/cleandata/NBR_FullSoilChem.csv", sep=",") # Clean soil chemistry data
soilpene_full <- read.csv("data/cleandata/NBR_FullSoilPene.csv", sep=",") # Clean soil bulk density data
vege_full <- read.csv("data/cleandata/NBR_FullPlantComm.csv", sep=",") # Clean plant community data
beetle_full <- read.csv("data/cleandata/NBR_FullBeetleComm.csv", sep=",") # Clean arthropod community data


#### DATA PREPARATION ####

#
## Filter grassland sites

# Site selection & removal geographical outlier IG3
siteinfo_grass <- siteinfo_full |>  
  filter(Habitat == "permanent grassland" & SiteID != "IG3") |> 
  dplyr::select(SiteID, Type_livestock)

# Extraction in other datasets
climate_grass <- filter(climate_full, SiteID %in% siteinfo_grass$SiteID)
landscape_grass <- filter(landscape_full, SiteID %in% siteinfo_grass$SiteID)
landuse_grass <- filter(landuse_full, SiteID %in% siteinfo_grass$SiteID)
area20x20_grass <- filter(area20x20_full, SiteID %in% siteinfo_grass$SiteID)
groundcover_grass <- filter(groundcover_full, SiteID %in% siteinfo_grass$SiteID)
soilbulk_grass <- filter(soilbulk_full, SiteID %in% siteinfo_grass$SiteID)
soilchem_grass <- filter(soilchem_full, SiteID %in% siteinfo_grass$SiteID)
soilpene_grass <- filter(soilpene_full, SiteID %in% siteinfo_grass$SiteID)
vege_grass <- filter(vege_full, SiteID %in% siteinfo_grass$SiteID)
beetle_grass <- filter(beetle_full, SiteID %in% siteinfo_grass$SiteID)

#
## Summarise data at site level -> should be 30 observations for non community data

# Site info - validated

# Climate - validated

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

# Beetle community - current at pitfall level -> summary by average
beetle_grass <- beetle_grass |> 
  group_by(SiteID, BeetleFamilies) |> 
  summarise(BeetleFam_abundance = mean(BeetleFam_abundance))

#
## Selection & transformation plant community data

# Data distribution
#hist(vege_grass$PlantSp_cover) # Poisson, highly skewed

# Species frequencies across sites
vege_freq <- filter(vege_grass, PlantSp_cover>0) |>
  group_by(Species) |>
  count() |>
  dplyr::arrange(desc(n)) # filtering by frequency
filter(vege_freq, n > 19)
# 6 grass species in at least 20 sites A. capillaris; F. rubra; H. lanatus; P. pratensis; A. odoratum; D. cespitosa
# 4 forb species in at least 20 sites T. repens; C. fontanum; C. palustre; R. acetosa
filter(vege_freq, n > 9 & n < 20)
# 2 grass species in at least 10 sites D. flexuosa; P. trivialis
# 10 forb species in at least 10 sites R. acris; P. erecta; G. saxatile; r. repens; A. millefolium; C. majus; H. radicata; C. pratensis; V. palustris; C. rotundifolia

# Species average percent cover across site
vege_average <- vege_grass |>
  group_by(Species) |>
  summarise_if(is.numeric, mean, na.rm = TRUE) |>
  dplyr::arrange(desc(PlantSp_cover)) # filtering by average value across sites
filter(vege_average, PlantSp_cover >= 3)
# 7 grass species over 3 % A. capillaris; F. rubra; H. lanatus; P. pratensis; D. cespitosa; A. odoratum; L. perenne
# 3 forb species over 3 % T. repens; R. acetosa; G. saxatile
filter(vege_average, PlantSp_cover < 3 & PlantSp_cover > 1)
# 2 grass species between 1 & 3% D. flexuosa; P. trivialis
# 4 forb species over 1%  P. erecta; R. repens; A. millefolium; R acris
# filter(vege_freq, Species == "Lolium perenne") only present in 6 sites

# Selection main grass species - at least present in 10 sites AND min average cover 1% -> 8 species
grass <- subset(vege_grass,
                Species == "Agrostis capillaris" |
                Species == "Festuca rubra" | 
                Species == "Holcus lanatus" | 
                Species == "Poa pratensis" | 
                Species == "Deschampsia cespitosa" | 
                Species == "Anthoxantum odoratum" | 
                Species == "Deschampsia flexuosa" | 
                Species == "Poa trivialis")

# Selection main forb species - at least present in 10 sites AND min average cover 1% -> 7 species
forb <- subset(vege_grass,
                Species == "Trifolium repens" |
                  Species == "Rumex acetosa" |
                  Species == "Galium saxatile" |
                  Species == "Potentilla erecta" |
                  Species == "Ranunculus acris" |
                  Species == "Ranunculus repens" |
                  Species == "Achillea millefolium")

# Log transformation
grass <- grass |> 
  mutate(PlantSp_logcover = log1p(PlantSp_cover))
forb <- forb |> 
  mutate(PlantSp_logcover = log1p(PlantSp_cover))

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

# Family frequencies across sites
beetle_freq <- filter(beetle_grass, BeetleFam_abundance>0) |>
  group_by(BeetleFamilies) |>
  count() |>
  dplyr::arrange(desc(n)) # filtering by frequency
filter(beetle_freq, n > 10)
# 8 beetle families in at least 10 sites Carab; Hydro; Scara; Staph; Ptili; Curcu; Elat; Silph

# Species average percent cover across site
beetle_average <- beetle_grass |>
  group_by(BeetleFamilies) |>
  summarise_if(is.numeric, mean, na.rm = TRUE) |>
  dplyr::arrange(desc(BeetleFam_abundance)) # filtering by average value across sites
filter(beetle_average, BeetleFam_abundance > 3)
# 6 beetles families with at least 3 individuals on average Staph; Ptili; Hyrdo; Scara; Carab; Silph

# Selection main dung beetle families -> at least present in 10 sites + min 3 individuals on average -> 6 families
beetle <- subset(beetle_grass,
                 BeetleFamilies == "Carabidae" |
                   BeetleFamilies == "Staphylinidae" |
                   BeetleFamilies == "Hydrophilidae" |
                   BeetleFamilies == "Ptiliidae" |
                   BeetleFamilies == "Scarabaeidae")#|
                   #BeetleFamilies == "Siphidae")

# Log transformation
beetle <- beetle |> 
  mutate(BeetleFam_logabundance = log1p(BeetleFam_abundance))

# Contingency table and wide table
contin_beetle <- xtabs(formula = BeetleFam_logabundance ~ SiteID + BeetleFamilies, data = beetle)
beetle <- subset(beetle, select = -c(BeetleFam_abundance))
beetle <- beetle |> 
  pivot_wider(names_from = BeetleFamilies, values_from = BeetleFam_logabundance)
beetle <- as.data.frame(beetle)

#
## Creation fjord system matrix

# Desired variables
## Elevation (num), from area 20x20
## Annual precipitation (num), from climate
## Average July temp (num), from climate
## Average Jan temp (num), from climate
## Distance to sea (num), from area 20x20

# Selection variables & scaling
fjordsys <- full_join(area20x20_grass, climate_grass)
fjordsys <- subset(fjordsys, select = c(SiteID, Elevation_max, AnnualPrecipitation, MaxTempJuly, MinTempJan, DistanceToSea_m))
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
## Average stocking density, from landuse

# Dummy numeric for original character variables
landuse_grass <- landuse_grass |> 
  mutate(Sheep = ifelse(Livestock1 == "sheep",1,0)) |> 
  mutate(Cow = ifelse(Livestock1 == "cow", 1,0))

# Selection variables & scaling
grazing <- subset(landuse_grass, select = c(SiteID, Sheep, Cow, FlockSize1_adults, GrazingSurface_ha, TotalInfieldSurface, AvgStockingDensity_perha))
grazing <- grazing |> 
  mutate(across(where(is.numeric), scale))

#
## Creation plant local environment matrix

# Desired variables
## Slope angle (num), from area 20x20
## Aspect angle (num), from area 20x20
## Bulk density (num), from soilbulk
## Moisture content (num), from soilbulk
## Penetration rate (num), from soilpene
## LOI (num), from soilchem
## Nitrogen content (num), from soilchem
## Phosphorous content (num), from soilchem
## pH (num), from soilchem
## Humus content (num), from soilchem

# Selection variables & scaling
locenvi_vege <- purrr::reduce(list(soilbulk_grass, soilpene_grass, soilchem_grass, area20x20_grass), dplyr::left_join)
locenvi_vege <- subset(locenvi_vege, select = c(SiteID, MeanBD, MeanMoisture, MeanLOI, MeanNitrogen, MeanPhosphorus, MeanpH, MeanHumus, MeanPT, AspectDegree, General_slope))
locenvi_vege <- locenvi_vege |> 
  mutate(across(where(is.numeric), scale))

#
## Creation beetle local environment matrix

# Desired variables
## Slope angle (num), from area 20x20
## Aspect angle (num), from area 20x20
## Litter cover (num), from groundcover
## Mean vegetation height (num), from groundcover
## Bulk density (num), from soilbulk
## Moisture content (num), from soilbulk
## Humus content (num), from soilchem

# Selection variables & scaling
locenvi_beetle <- purrr::reduce(list(groundcover_grass, soilbulk_grass, soilpene_grass, soilchem_grass, area20x20_grass), dplyr::left_join)
locenvi_beetle <- locenvi_beetle |> 
  mutate(MeanExposedGround = MeanBareSoil + MeanRocks)
locenvi_beetle <- subset(locenvi_beetle, select = c(SiteID, MeanLitter, MeanBryo, MeanHeight, MeanBD, MeanPT, MeanMoisture, MeanHumus, AspectDegree, General_slope))
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
      lower.panel = panel.smooth) # high colinearity between max temp and min temp

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
fjordsys <- subset(fjordsys, select = -c(MinTempJan))
locenvi_vege <- subset(locenvi_vege, select = -c(MeanMoisture, MeanLOI, MeanHumus, MeanNitrogen))
locenvi_beetle <- subset(locenvi_beetle, select = -c(MeanMoisture, MeanHumus))

#
## Linear relationships of residuals

# Ordination community data
DCA_grass <- decorana(contin_grass)
DCA_grass # Axis length of 2.1 -> keep DCA
plot(DCA_grass)
DCA_forb <- decorana(contin_forb)
DCA_forb # Axis length of 2.7 -> keep DCA
plot(DCA_forb)
DCA_beetle <- decorana(contin_beetle)
DCA_beetle # Axis length at 1.3 -> run PCA
PCA_beetle <- prcomp(contin_beetle)
PCA_beetle
biplot(PCA_beetle)

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
# residualPlots(lm(GrassDCA1~Elevation_max, data = allvar_grass))
# residualPlots(lm(GrassDCA1~DistanceToSea_m, data = allvar_grass))
# residualPlots(lm(GrassDCA1~AnnualPrecipitation, data = allvar_grass))
# residualPlots(lm(GrassDCA1~MaxTempJuly, data = allvar_grass))
# residualPlots(lm(GrassDCA1~TotCultivatedLand_percent, data = allvar_grass))
# residualPlots(lm(GrassDCA1~TotForest_percent, data = allvar_grass))
# residualPlots(lm(GrassDCA1~Infield_percent, data = allvar_grass))
# residualPlots(lm(GrassDCA1~Outfield_percent, data = allvar_grass))
residualPlots(lm(GrassDCA1~Wetland_percent, data = allvar_grass)) #rejected
#residualPlots(lm(GrassDCA1~Sheep, data = allvar_grass)) # binary, cannot be checked
#residualPlots(lm(GrassDCA1~Cow, data = allvar_grass)) # binary, cannot be checked
#residualPlots(lm(GrassDCA1~Goat, data = allvar_grass)) # binary, cannot be checked
residualPlots(lm(GrassDCA1~FlockSize1_adults, data = allvar_grass)) #rejected
# residualPlots(lm(GrassDCA1~GrazingSurface_ha, data = allvar_grass))
# residualPlots(lm(GrassDCA1~TotalInfieldSurface, data = allvar_grass))
# residualPlots(lm(GrassDCA1~AvgStockingDensity_perha, data = allvar_grass))
# residualPlots(lm(GrassDCA1~General_slope, data = allvar_grass))
# residualPlots(lm(GrassDCA1~AspectDegree, data = allvar_grass))
# residualPlots(lm(GrassDCA1~MeanBD, data = allvar_grass))
# residualPlots(lm(GrassDCA1~MeanPT, data = allvar_grass))
# residualPlots(lm(GrassDCA1~MeanpH, data = allvar_grass))
# residualPlots(lm(GrassDCA1~MeanPhosphorus, data = allvar_grass))

# Residuals LM for forbs
# residualPlots(lm(ForbDCA1~Elevation_max, data = allvar_forb))
# residualPlots(lm(ForbDCA1~DistanceToSea_m, data = allvar_forb))
residualPlots(lm(ForbDCA1~AnnualPrecipitation, data = allvar_forb)) #rejected
# residualPlots(lm(ForbDCA1~MaxTempJuly, data = allvar_forb))
# residualPlots(lm(ForbDCA1~TotCultivatedLand_percent, data = allvar_forb))
# residualPlots(lm(ForbDCA1~TotForest_percent, data = allvar_forb))
# residualPlots(lm(ForbDCA1~Infield_percent, data = allvar_forb))
# residualPlots(lm(ForbDCA1~Outfield_percent, data = allvar_forb))
# residualPlots(lm(ForbDCA1~Wetland_percent, data = allvar_forb))
# residualPlots(lm(ForbDCA1~Sheep, data = allvar_forb)) # binary, cannot be checked
# residualPlots(lm(ForbDCA1~Cow, data = allvar_forb)) # binary, cannot be checked
# residualPlots(lm(ForbDCA1~Goat, data = allvar_forb)) # binary, cannot be checked
residualPlots(lm(ForbDCA1~FlockSize1_adults, data = allvar_forb)) #rejected
# residualPlots(lm(ForbDCA1~GrazingSurface_ha, data = allvar_forb))
# residualPlots(lm(ForbDCA1~TotalInfieldSurface, data = allvar_forb))
# residualPlots(lm(ForbDCA1~AvgStockingDensity_perha, data = allvar_forb))
# residualPlots(lm(ForbDCA1~General_slope, data = allvar_forb))
residualPlots(lm(ForbDCA1~AspectDegree, data = allvar_forb)) #rejected
# residualPlots(lm(ForbDCA1~MeanBD, data = allvar_forb))
# residualPlots(lm(ForbDCA1~MeanPT, data = allvar_forb))
# residualPlots(lm(ForbDCA1~MeanpH, data = allvar_forb))
# residualPlots(lm(ForbDCA1~MeanPhosphorus, data = allvar_forb))

# Residuals LM for beetles
# residualPlots(lm(BeetlePCA1~Elevation_max, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~DistanceToSea_m, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~AnnualPrecipitation, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~MaxTempJuly, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~TotCultivatedLand_percent, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~TotForest_percent, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~Infield_percent, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~Outfield_percent, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~Wetland_percent, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~Sheep, data = allvar_beetle)) # binary, cannot be checked
# residualPlots(lm(BeetlePCA1~Cow, data = allvar_beetle)) # binary, cannot be checked
# residualPlots(lm(BeetlePCA1~Goat, data = allvar_beetle)) # binary, cannot be checked
# residualPlots(lm(BeetlePCA1~FlockSize1_adults, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~GrazingSurface_ha, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~TotalInfieldSurface, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~AvgStockingDensity_perha, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~General_slope, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~AspectDegree, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~MeanBD, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~MeanPT, data = allvar_beetle))
residualPlots(lm(BeetlePCA1~MeanLitter, data = allvar_beetle)) # rejected
# residualPlots(lm(BeetlePCA1~MeanBryo, data = allvar_beetle))
# residualPlots(lm(BeetlePCA1~MeanHeight, data = allvar_beetle))



#### Canonical correlation preparation ####

#
## Data preparation - Contingency tables

# Fjord system
fjordsys_long <- fjordsys |> 
  pivot_longer(
    cols = c(Elevation_max, AnnualPrecipitation, MaxTempJuly, DistanceToSea_m),
    names_to = "Factors",
    values_to = "Values")
contin_fjordsys_grass <- xtabs(formula = Values ~ SiteID + Factors, data = fjordsys_long)
contin_fjordsys_forb <- xtabs(formula = Values ~ SiteID + Factors, data = filter(fjordsys_long, Factors != "AnnualPrecipitation")) # non-linear relationship
contin_fjordsys_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = fjordsys_long)

# Landscape
landscape_long <- landscape |> 
  pivot_longer(
    cols = c(TotCultivatedLand_percent, TotForest_percent, Infield_percent, Outfield_percent, Wetland_percent),
    names_to = "Factors",
    values_to = "Values"
  )
contin_landscape_grass <- xtabs(formula = Values ~ SiteID + Factors, data = filter(landscape_long, Factors != "Wetland_percent"))
contin_landscape_forb <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long)
contin_landscape_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long)

# Grazing
grazing_long <- grazing |> 
  pivot_longer(
    cols = c(Sheep, Cow, FlockSize1_adults, GrazingSurface_ha, TotalInfieldSurface, AvgStockingDensity_perha),
    names_to = "Factors",
    values_to = "Values")
contin_grazing_grass <- xtabs(formula = Values ~ SiteID + Factors, data = filter(grazing_long, Factors != "FlockSize1_adults")) # non-linear relationship
contin_grazing_forb <- xtabs(formula = Values ~ SiteID + Factors, data = filter(grazing_long, Factors != "FlockSize1_adults")) # non-linear relationship
contin_grazing_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = grazing_long)

# Local environment - plants
locenvi_long_vege <- locenvi_vege |> 
  pivot_longer(
    cols = c(MeanBD, MeanPT, MeanpH, MeanPhosphorus, General_slope, AspectDegree),
    names_to = "Factors",
    values_to = "Values")
contin_locenvi_grass <- xtabs(formula = Values ~ SiteID + Factors, data = locenvi_long_vege)
contin_locenvi_forb <- xtabs(formula = Values ~ SiteID + Factors, data = filter(locenvi_long_vege, Factors != "AspectDegree"))

# Local environment - beetle
locenvi_long_beetle <- locenvi_beetle |> 
  pivot_longer(
    cols = c(MeanHeight, MeanPT, MeanBD, MeanBryo, MeanLitter, AspectDegree, General_slope),
    names_to = "Factors",
    values_to = "Values")
contin_locenvi_beetle <- xtabs(formula = Values ~ SiteID + Factors, data = filter(locenvi_long_beetle, Factors != "MeanLitter"))

#
## Graphical representation of correlations plant community

# Fjord effect
# correl_fjordxlandscape_vege <- matcor(contin_fjordsys, contin_landscape_vege)
# img.matcor(correl_fjordxlandscape_vege, type = 2)
# correl_fjordxgrazing_vege <- matcor(contin_fjordsys, contin_grazing_vege)
# img.matcor(correl_fjordxgrazing_vege, type = 2)
# correl_fjordxlocenvi_vege <- matcor(contin_fjordsys, contin_locenvi_vege)
# img.matcor(correl_fjordxlocenvi_vege, type = 2)
# correl_fjordxvege <- matcor(contin_fjordsys, contin_vege)
# img.matcor(correl_fjordxvege, type = 2)

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


#### Canonical correlation analysis between driver sets ####

# Number of observations (same for all sets)
# nobs <- dim(contin_fjordsys)[1]
# 
# nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
# nvar_landscape <- length(select_if(landscape, is.numeric))
# nvar_grazing <- length(select_if(grazing, is.numeric))
# nvar_locenvi_vege <- length(select_if(locenvi_vege, is.numeric))
# nvar_locenvi_beetle <- length(select_if(locenvi_vege, is.numeric))
# 
# #
# ## Fjord effect
# 
# # Fjord x landscape
# cancor_fjordxlandscape <- cc(contin_fjordsys, contin_landscape)
# rho_fjordxlandscape <- cancor_fjordxlandscape$cor
# rho_fjordxlandscape
# # 1st axis correlation 0.79
# # 2nd axis correlation 0.74
# # 3rd axis correlation 0.55
# p.asym(rho_fjordxlandscape, nobs, nvar_fjordsys, nvar_landscape, tstat = "Hotelling") 
# # 1st dim significant - stat 3.45 - df1 20 - df2 74 - pval 1.5.10-4
# # 2nd dim significant - stat 1.77 - df1 12 - df2 82 - pval 1.5.10-3
# # 3rd dim marginally significant - stat 0.5 - df1 6 - df2 90 - pval 0.09
# plt.cc(cancor_fjordxlandscape, var.label = TRUE)
# 
# # Fjord x grazing
# cancor_fjordxgrazing <- cc(contin_fjordsys, contin_grazing)
# rho_fjordxgrazing <- cancor_fjordxgrazing$cor
# rho_fjordxgrazing 
# # 1st axis correlation 0.61
# p.asym(rho_fjordxgrazing, nobs, nvar_fjordsys, nvar_grazing, tstat = "Hotelling") 
# # 1st dim NS - stat 1.11 - df1 24 - df2 70 - pval 0.72
# plt.cc(cancor_fjordxgrazing, var.label = TRUE)
# 
# # Fjord x plant local environment
# cancor_fjordxlocenvi_vege <- cc(contin_fjordsys, contin_locenvi_vege)
# rho_fjordxlocenvi_vege <- cancor_fjordxlocenvi_vege$cor
# rho_fjordxlocenvi_vege 
# # 1st axis correlation 0.58
# p.asym(rho_fjordxlocenvi_vege, nobs, nvar_fjordsys, nvar_locenvi_vege, tstat = "Hotelling") 
# # 1st dim NS - stat 0.95 - df1 24 - df2 70 - pval 0.84
# plt.cc(cancor_fjordxlocenvi_vege, var.label = TRUE)
# 
# # Fjord x beetle local environment
# cancor_fjordxlocenvi_beetle <- cc(contin_fjordsys, contin_locenvi_beetle)
# rho_fjordxlocenvi_beetle <- cancor_fjordxlocenvi_beetle$cor
# rho_fjordxlocenvi_beetle 
# # 1st axis correlation 0.61
# p.asym(rho_fjordxlocenvi_beetle, nobs, nvar_fjordsys, nvar_locenvi_beetle, tstat = "Hotelling") 
# # 1st dim NS - stat 0.91 - df1 24 - df2 70 - pval 0.87
# plt.cc(cancor_fjordxlocenvi_beetle, var.label = TRUE)
# 
# #
# ## Landscape effect
# 
# cancor_landscapexgrazing <- cc(contin_landscape, contin_grazing)
# rho_landscapexgrazing <- cancor_landscapexgrazing$cor
# rho_landscapexgrazing 
# # 1st axis correlation 0.69
# p.asym(rho_landscapexgrazing, nobs, nvar_landscape, nvar_grazing, tstat = "Hotelling") 
# # 1st dim marginally significant - stat 1.96 - df1 30 - df2 82 - pval 0.39
# plt.cc(cancor_landscapexgrazing, var.label = TRUE)
# 
# # Landscape x plant local environment
# cancor_landscapexlocalenvi_vege <- cc(contin_landscape, contin_locenvi_vege)
# rho_landscapexlocalenvi_vege <- cancor_landscapexlocalenvi_vege$cor
# rho_landscapexlocalenvi_vege 
# # 1st axis correlation 0.77
# p.asym(rho_landscapexlocalenvi_vege, nobs, nvar_landscape, nvar_locenvi_vege, tstat = "Hotelling") 
# # 1st dim marginally significant - stat 2.68 - df1 30 - df2 82 - pval 0.09
# plt.cc(cancor_landscapexlocalenvi_vege, var.label = TRUE)
# 
# # Landscape x beetle local environment
# cancor_landscapexlocalenvi_beetle <- cc(contin_landscape, contin_locenvi_beetle)
# rho_landscapexlocalenvi_beetle <- cancor_landscapexlocalenvi_beetle$cor
# rho_landscapexlocalenvi_beetle 
# # 1st axis correlation 0.79
# p.asym(rho_landscapexlocalenvi_beetle, nobs, nvar_landscape, nvar_locenvi_beetle, tstat = "Hotelling") # 1st dim marginally significant - stat 2.67 - df1 30 - df2 82 - pval 0.092
# plt.cc(cancor_landscapexlocalenvi_beetle, var.label = TRUE)
# 
# #
# ## Grazing effect
# 
# # Grazing x plant local environment
# cancor_grazingxlocenvi_vege <- cc(contin_grazing, contin_locenvi_vege)
# rho_grazingxlocenvi_vege <- cancor_grazingxlocenvi_vege$cor
# rho_grazingxlocenvi_vege 
# # 1st axis correlation 0.67
# p.asym(rho_grazingxlocenvi_vege, nobs, nvar_grazing_vege, nvar_locenvi_vege, tstat = "Hotelling") 
# # 1st dim NS - stat 1.52 - df1 36 - df2 92 - pval 0.93
# plt.cc(cancor_grazingxlocenvi_vege, var.label = TRUE)
# 
# # Grazing x beetle local environment
# cancor_grazingxlocenvi_beetle <- cc(contin_grazing, contin_locenvi_beetle)
# rho_grazingxlocenvi_beetle <- cancor_grazingxlocenvi_beetle$cor
# rho_grazingxlocenvi_beetle 
# # 1st axis correlation 0.77
# p.asym(rho_grazingxlocenvi_beetle, nobs, nvar_grazing, nvar_locenvi_beetle, tstat = "Hotelling") 
# # 1st dim NS - stat 2.42 - df1 36 - df2 92 - pval 0.44
# plt.cc(cancor_grazingxlocenvi_beetle, var.label = TRUE)


#### Canonical correlation analysis grass community ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys_grass)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape <- length(select_if(subset(landscape, select = -c(Wetland_percent)), is.numeric))
nvar_grazing <- length(select_if(subset(grazing, select = -c(FlockSize1_adults)), is.numeric))
nvar_locenvi <- length(select_if(locenvi_vege, is.numeric))
nvar_grass <- length(select_if(grass, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor_fjordxlandscape_grass <- cc(contin_fjordsys_grass, contin_landscape_grass)
rho_fjordxlandscape_grass <- cancor_fjordxlandscape_grass$cor
rho_fjordxlandscape_grass
# 1st axis correlation 0.79
# 2nd axis correlation 0.63
p.asym(rho_fjordxlandscape_grass, nobs, nvar_fjordsys, nvar_landscape, tstat = "Hotelling")
# 1st dim significant - stat 2.44 - df1 16 - df2 78 - pval 7.1.10-4
# 2nd dim marginally significant - stat 0.77 - df1 9 - df2 86 - pval 0.07
plt.cc(cancor_fjordxlandscape_grass, var.label = TRUE)

# Fjord x grazing
cancor_fjordxgrazing_grass <- cc(contin_fjordsys_grass, contin_grazing_grass)
rho_fjordxgrazing_grass <- cancor_fjordxgrazing_grass$cor
rho_fjordxgrazing_grass
# 1st axis correlation 0.6
p.asym(rho_fjordxgrazing_grass, nobs, nvar_fjordsys, nvar_grazing, tstat = "Hotelling")
# 1st dim NS - stat 0.88 - df1 20 - df2 74 - pval 0.69
plt.cc(cancor_fjordxgrazing_grass, var.label = TRUE)

# Fjord x local environment
cancor_fjordxlocenvi_grass <- cc(contin_fjordsys_grass, contin_locenvi_grass)
rho_fjordxlocenvi_grass <- cancor_fjordxlocenvi_grass$cor
rho_fjordxlocenvi_grass
# 1st axis correlation 0.58
p.asym(rho_fjordxlocenvi_grass, nobs, nvar_fjordsys, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 0.95 - df1 24 - df2 70 - pval 0.84
plt.cc(cancor_fjordxlocenvi_grass, var.label = TRUE)

# Fjord x grass community
cancor_fjordxgrass <- cc(contin_fjordsys_grass, contin_grass)
rho_fjordxgrass <- cancor_fjordxgrass$cor
rho_fjordxgrass 
# 1st axis correlation 0.8
p.asym(rho_fjordxgrass, nobs, nvar_fjordsys, nvar_grass, tstat = "Hotelling") 
# 1st dim marginally significant - stat 2.59 - df1 28 - df2 66 - pval 0.081
plt.cc(cancor_fjordxgrass, var.label = TRUE)

#
## Landscape effect

# Landscape x grazing
cancor_landscapexgrazing_grass <- cc(contin_landscape_grass, contin_grazing_grass)
rho_landscapexgrazing_grass <- cancor_landscapexgrazing_grass$cor
rho_landscapexgrazing_grass
# 1st axis correlation 0.6
p.asym(rho_landscapexgrazing_grass, nobs, nvar_landscape, nvar_grazing, tstat = "Hotelling")
# 1st dim marginally significant - stat 1.12 - df1 20 - df2 74 - pval 0.44
plt.cc(cancor_landscapexgrazing_grass, var.label = TRUE)

# Landscape x local environment
cancor_landscapexlocalenvi_grass <- cc(contin_landscape_grass, contin_locenvi_grass)
rho_landscapexlocalenvi_grass <- cancor_landscapexlocalenvi_grass$cor
rho_landscapexlocalenvi_grass
# 1st axis correlation 0.77
p.asym(rho_landscapexlocalenvi_grass, nobs, nvar_landscape, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 2 - df1 24 - df2 70 - pval 0.11
plt.cc(cancor_landscapexlocalenvi_grass, var.label = TRUE)

# Landscape x grass community
cancor_landscapexgrass <- cc(contin_landscape_grass, contin_grass)
rho_landscapexgrass <- cancor_landscapexgrass$cor
rho_landscapexgrass 
# 1st axis correlation 0.71
p.asym(rho_landscapexgrass, nobs, nvar_landscape, nvar_grass, tstat = "Hotelling") 
# 1st dim NS - stat 1.53 - df1 28 - df2 66 - pval 0.61
plt.cc(cancor_landscapexgrass, var.label = TRUE)

#
## Grazing effect

# Grazing x local environment
cancor_grazingxlocenvi_grass <- cc(contin_grazing_grass, contin_locenvi_grass)
rho_grazingxlocenvi_grass <- cancor_grazingxlocenvi_grass$cor
rho_grazingxlocenvi_grass
# 1st axis correlation 0.6
p.asym(rho_grazingxlocenvi_grass, nobs, nvar_grazing, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 1.12 - df1 30 - df2 82 - pval 0.94
plt.cc(cancor_grazingxlocenvi_grass, var.label = TRUE)

# Grazing x grass community
cancor_grazingxgrass <- cc(contin_grazing_grass, contin_grass)
rho_grazingxgrass <- cancor_grazingxgrass$cor
rho_grazingxgrass 
# 1st axis correlation 0.87 
p.asym(rho_grazingxgrass, nobs, nvar_grazing, nvar_grass, tstat = "Hotelling") 
# 1st axis significant - stat 4.6 - df1 35 - df2 77 - pval 5.2.10-3
plt.cc(cancor_grazingxgrass, var.label = TRUE)
plt.cc(cancor_grazingxgrass, d1 = 2, d2 = 3, var.label = TRUE) # dim 2 & 3

#
## Local environment effect

cancor_locenvixgrass <- cc(contin_locenvi_grass, contin_grass)
rho_locenvixgrass <- cancor_locenvixgrass$cor
rho_locenvixgrass 
# 1st axis correlation 0.87
sigsum <- p.asym(rho_locenvixgrass, nobs, nvar_locenvi, nvar_grass, tstat = "Hotelling")
sigsum
# 1st axis significant - stat 4.63 - df1 42 - df2 86 - pval 0.038
summary_locenvixgrass <- cbind(rho_locenvixgrass, sigsum$approx, sigsum$df1, sigsum$df2, sigsum$p.value)
summary_locenvixgrass
# plt.cc(cancor_locenvixgrass, var.label = TRUE)
# plt.cc(cancor_locenvixgrass, d1 = 2, d2 = 3, var.label = TRUE) # dim 2 & 3




#### Canonical correlation analysis forb community ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys_forb)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(subset(fjordsys, select = -c(AnnualPrecipitation)), is.numeric))
nvar_landscape <- length(select_if(landscape, is.numeric))
nvar_grazing <- length(select_if(subset(grazing, select = -c(FlockSize1_adults)), is.numeric))
nvar_locenvi <- length(select_if(subset(locenvi_vege, select = -c(AspectDegree)), is.numeric))
nvar_forb <- length(select_if(forb, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor_fjordxlandscape_forb <- cc(contin_fjordsys_forb, contin_landscape_forb)
rho_fjordxlandscape_forb <- cancor_fjordxlandscape_forb$cor
rho_fjordxlandscape_forb
# 1st axis correlation 0.79
# 2nd axis correlation 0.75
# 3rd axis correlation 0.45
p.asym(rho_fjordxlandscape_forb, nobs, nvar_fjordsys, nvar_landscape, tstat = "Hotelling")
# 1st dim significant - stat 3.22 - df1 15 - df2 59 - pval 3.18.10-5
# 2nd dim significant - stat 1.55 - df1 8 - df2 65 - pval 4.23.10-4
# 3rd dim marginally significant - stat 0.28 - df1 3 - df2 71 - pval 0.092
plt.cc(cancor_fjordxlandscape_forb, var.label = TRUE) # dim 1 & 2
plt.cc(cancor_fjordxlandscape_forb, d1 = 1, d2 = 3, var.label = TRUE) # dim 1 & 3

# Fjord x grazing
cancor_fjordxgrazing_forb <- cc(contin_fjordsys_forb, contin_grazing_forb)
rho_fjordxgrazing_forb <- cancor_fjordxgrazing_forb$cor
rho_fjordxgrazing_forb
# 1st axis correlation 0.52
p.asym(rho_fjordxgrazing_forb, nobs, nvar_fjordsys, nvar_grazing, tstat = "Hotelling")
# 1st dim NS - stat 0.58 - df1 15 - df2 59 - pval 0.72
plt.cc(cancor_fjordxgrazing_forb, var.label = TRUE)

# Fjord x local environment
cancor_fjordxlocenvi_forb <- cc(contin_fjordsys_forb, contin_locenvi_forb)
rho_fjordxlocenvi_forb <- cancor_fjordxlocenvi_forb$cor
rho_fjordxlocenvi_forb
# 1st axis correlation 0.58
p.asym(rho_fjordxlocenvi_forb, nobs, nvar_fjordsys, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 0.71 - df1 15 - df2 59 - pval 0.55
plt.cc(cancor_fjordxlocenvi_forb, var.label = TRUE)

# Fjord x forb community
cancor_fjordxforb <- cc(contin_fjordsys_forb, contin_forb)
rho_fjordxforb <- cancor_fjordxforb$cor
rho_fjordxforb 
# 1st axis correlation 0.63
p.asym(rho_fjordxforb, nobs, nvar_fjordsys, nvar_forb, tstat = "Hotelling") 
# 1st dim NS - stat 1.57 - df1 21 - df2 53 - pval 0.78
plt.cc(cancor_fjordxforb, var.label = TRUE)

#
## Landscape effect

# Landscape x grazing
cancor_landscapexgrazing_forb <- cc(contin_landscape_forb, contin_grazing_forb)
rho_landscapexgrazing_forb <- cancor_landscapexgrazing_forb$cor
rho_landscapexgrazing_forb
# 1st axis correlation 0.67
p.asym(rho_landscapexgrazing_forb, nobs, nvar_landscape, nvar_grazing, tstat = "Hotelling")
# 1st dim NS - stat 1.48 - df1 25 - df2 87 - pval 0.44
plt.cc(cancor_landscapexgrazing_forb, var.label = TRUE)

# Landscape x local environment
cancor_landscapexlocalenvi_forb <- cc(contin_landscape_forb, contin_locenvi_forb)
rho_landscapexlocalenvi_forb <- cancor_landscapexlocalenvi_forb$cor
rho_landscapexlocalenvi_forb
# 1st axis correlation 0.77
p.asym(rho_landscapexlocalenvi_forb, nobs, nvar_landscape, nvar_locenvi, tstat = "Hotelling")
# 1st dim significant - stat 2.37 - df1 25 - df2 87 - pval 0.046
plt.cc(cancor_landscapexlocalenvi_forb, var.label = TRUE)

# Landscape x forb community
cancor_landscapexforb <- cc(contin_landscape_forb, contin_forb)
rho_landscapexforb <- cancor_landscapexforb$cor
rho_landscapexforb 
# 1st axis correlation 0.75
p.asym(rho_landscapexforb, nobs, nvar_landscape, nvar_forb, tstat = "Hotelling") 
# 1st dim NS - stat 3 - df1 35 - df2 77 - pval 0.15
plt.cc(cancor_landscapexforb, var.label = TRUE)

#
## Grazing effect

# Grazing x local environment
cancor_grazingxlocenvi_forb <- cc(contin_grazing_forb, contin_locenvi_forb)
rho_grazingxlocenvi_forb <- cancor_grazingxlocenvi_forb$cor
rho_grazingxlocenvi_forb
# 1st axis correlation 0.59
p.asym(rho_grazingxlocenvi_forb, nobs, nvar_grazing, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 0.96 - df1 25 - df2 87 - pval 0.87
plt.cc(cancor_grazingxlocenvi_forb, var.label = TRUE)

# Grazing x forb community
cancor_grazingxforb <- cc(contin_grazing_forb, contin_forb)
rho_grazingxforb <- cancor_grazingxforb$cor
rho_grazingxforb 
# 1st axis correlation 0.79
p.asym(rho_grazingxforb, nobs, nvar_grazing, nvar_forb, tstat = "Hotelling") 
# 1st dim NS - stat 2.89 - df1 35 - df2 77 - pval 0.18
plt.cc(cancor_grazingxforb, var.label = TRUE)

#
## Local environment effect

cancor_locenvixforb <- cc(contin_locenvi_forb, contin_forb)
rho_locenvixforb <- cancor_locenvixforb$cor
rho_locenvixforb 
# 1st axis correlation 0.77
p.asym(rho_locenvixforb, nobs, nvar_locenvi, nvar_forb, tstat = "Hotelling") 
# 1st dim marginally significant - stat 3.4 - df1 35 - df2 77 - pval 0.07
plt.cc(cancor_locenvixforb, var.label = TRUE)
plt.cc(cancor_locenvixforb, d1 = 1, d2 = 3, var.label = TRUE) # dim 1 & 3


#### Canonical correlation analysis beetle assemblage ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys_beetle)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape <- length(select_if(landscape, is.numeric))
nvar_grazing <- length(select_if(grazing, is.numeric))
nvar_locenvi <- length(select_if(subset(locenvi_beetle, select = -c(MeanLitter)), is.numeric))
nvar_beetle <- length(select_if(beetle, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor_fjordxlandscape_beetle <- cc(contin_fjordsys_beetle, contin_landscape_beetle)
rho_fjordxlandscape_beetle <- cancor_fjordxlandscape_beetle$cor
rho_fjordxlandscape_beetle
# 1st axis correlation 0.79
# 2nd axis correlation 0.75
# 3rd axis correlation 0.55
p.asym(rho_fjordxlandscape_beetle, nobs, nvar_fjordsys, nvar_landscape, tstat = "Hotelling")
# 1st dim significant - stat 3.45 - df1 20 - df2 74 - pval 1.48.10-4
# 2nd dim significant - stat 1.77 - df1 12 - df2 82 - pval 1.48-3
# 3rd dim marginally significant - stat 0.5 - df1 6 - df2 90 - pval 0.093
plt.cc(cancor_fjordxlandscape_beetle, var.label = TRUE) # dim 1 & 2
plt.cc(cancor_fjordxlandscape_beetle, d1 = 1, d2 = 3, var.label = TRUE) # dim 1 & 3

# Fjord x grazing
cancor_fjordxgrazing_beetle <- cc(contin_fjordsys_beetle, contin_grazing_beetle)
rho_fjordxgrazing_beetle <- cancor_fjordxgrazing_beetle$cor
rho_fjordxgrazing_beetle
# 1st axis correlation 0.61
p.asym(rho_fjordxgrazing_beetle, nobs, nvar_fjordsys, nvar_grazing, tstat = "Hotelling")
# 1st dim NS - stat 1.11 - df1 24 - df2 70 - pval 0.72
plt.cc(cancor_fjordxgrazing_beetle, var.label = TRUE)

# Fjord x local environment
cancor_fjordxlocenvi_beetle <- cc(contin_fjordsys_beetle, contin_locenvi_beetle)
rho_fjordxlocenvi_beetle <- cancor_fjordxlocenvi_beetle$cor
rho_fjordxlocenvi_beetle
# 1st axis correlation 0.61
p.asym(rho_fjordxlocenvi_beetle, nobs, nvar_fjordsys, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 0.91 - df1 24 - df2 70 - pval 0.87
plt.cc(cancor_fjordxlocenvi_beetle, var.label = TRUE)

# Fjord x beetle assemblage
cancor_fjordxbeetle <- cc(contin_fjordsys_beetle, contin_beetle)
rho_fjordxbeetle <- cancor_fjordxbeetle$cor
rho_fjordxbeetle
# 1st axis correlation 0.69
p.asym(rho_fjordxbeetle, nobs, nvar_fjordsys, nvar_beetle, tstat = "Hotelling")
# 1st dim NS - stat 1.38 - df1 20 - df2 74 - pval 0.22
plt.cc(cancor_fjordxbeetle, var.label = TRUE)

#
## Landscape effect

# Landscape x grazing
cancor_landscapexgrazing_beetle <- cc(contin_landscape_beetle, contin_grazing_beetle)
rho_landscapexgrazing_beetle <- cancor_landscapexgrazing_beetle$cor
rho_landscapexgrazing_beetle
# 1st axis correlation 0.69
p.asym(rho_landscapexgrazing_beetle, nobs, nvar_landscape, nvar_grazing, tstat = "Hotelling")
# 1st dim NS - stat 1.96 - df1 30 - df2 82 - pval 0.39
plt.cc(cancor_landscapexgrazing_beetle, var.label = TRUE)

# Landscape x local environment
cancor_landscapexlocalenvi_beetle <- cc(contin_landscape_beetle, contin_locenvi_beetle)
rho_landscapexlocalenvi_beetle <- cancor_landscapexlocalenvi_beetle$cor
rho_landscapexlocalenvi_beetle
# 1st axis correlation 0.79
p.asym(rho_landscapexlocalenvi_beetle, nobs, nvar_landscape, nvar_locenvi, tstat = "Hotelling")
# 1st dim marginally significant - stat 2.67 - df1 30 - df2 82 - pval 0.092
plt.cc(cancor_landscapexlocalenvi_beetle, var.label = TRUE)

# Landscape x beetle assemblage
cancor_landscapexbeetle <- cc(contin_landscape_beetle, contin_beetle)
rho_landscapexbeetle <- cancor_landscapexbeetle$cor
rho_landscapexbeetle
# 1st axis correlation 0.71
p.asym(rho_landscapexbeetle, nobs, nvar_landscape, nvar_beetle, tstat = "Hotelling")
# 1st dim NS - stat 1.76 - df1 25 - df2 87 - pval 0.24
plt.cc(cancor_landscapexbeetle, var.label = TRUE)

#
## Grazing effect

# Grazing x local environment
cancor_grazingxlocenvi_beetle <- cc(contin_grazing_beetle, contin_locenvi_beetle)
rho_grazingxlocenvi_beetle <- cancor_grazingxlocenvi_beetle$cor
rho_grazingxlocenvi_beetle
# 1st axis correlation 0.77
p.asym(rho_grazingxlocenvi_beetle, nobs, nvar_grazing, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 2.42 - df1 36 - df2 92 - pval 0.44
plt.cc(cancor_grazingxlocenvi_beetle, var.label = TRUE)

# Grazing x beetle community
cancor_grazingxbeetle <- cc(contin_grazing_beetle, contin_beetle)
rho_grazingxbeetle <- cancor_grazingxbeetle$cor
rho_grazingxbeetle
# 1st axis correlation 0.72
p.asym(rho_grazingxbeetle, nobs, nvar_grazing, nvar_beetle, tstat = "Hotelling")
# 1st dim NS - stat 1.73 - df1 30 - df2 82 - pval 0.55
plt.cc(cancor_grazingxbeetle, var.label = TRUE)

#
## Local environment effect

cancor_locenvixbeetle <- cc(contin_locenvi_beetle, contin_beetle)
VIScancor_locenvixbeetle <- cancor(contin_locenvi_beetle, contin_beetle)
rho_locenvixbeetle <- cancor_locenvixbeetle$cor
rho_locenvixbeetle
# 1st axis correlation 0.73
p.asym(rho_locenvixbeetle, nobs, nvar_locenvi, nvar_beetle, tstat = "Hotelling")
# 1st dim NS - stat 2.5 - df1 30 - df2 82 - pval 0.14
plt.cc(cancor_locenvixbeetle, var.label = TRUE)
#VIScancor_locenvixbeetle <- cancor(contin_locenvi_beetle, contin_beetle)
#plot(cancor_locenvixbeetle, which = 1)


#### Redundancy analysis on significant results ####

## Fjord x landscape highly correlated & significant
## Landscape x fine-scale correlation significant to marginally significant
## Fjord x grass highly correlated but marginally significant
## Grazing x grass significantly correlated
## Fine-scale x grass significantly correlated
## Fine-scale x forbs marginally significant correlation
## Fine-scale x beetles NS

#
## Fjord x landscape

# Redundancy analysis and plot
rda_fjordxlandscape <- rda(subset(landscape, select = -c(SiteID)) ~ Elevation_max + AnnualPrecipitation + MaxTempJuly + DistanceToSea_m, data = fjordsys)
plotrda_fjordxlandscape <- ggord(rda_fjordxlandscape,
      txt = 4,
      ptslab = TRUE,
      addsize = 3,
      size = 1,
      arrow = 0.3,
      #repel = TRUE,
      vec_lab = list(MaxTempJuly = "JulTemp", Elevation_max = "Elev", DistanceToSea_m = "SeaDist", AnnualPrecipitation = "AnnPreci"),
      xlims = c(-1.3, 1),
      ylims = c(-0.8, 1.4))
plotrda_fjordxlandscape

# Summary rda
#summary(rda_fjordxlandscape)
# 1st axis driven by maximum July temperature (-0.87)
# 2nd axis driven by distance to sea (0.9)
# 1st axis drives forest cover (-0.71) and outfield cover (0.75)
# 2nd axis drives wetland cover (0.87)

# Total variance explained by RDA
RsquareAdj(rda_fjordxlandscape)
# R2 0.26
# Adj R2 0.14

# Global RDA significance by permutation
anova.cca(rda_fjordxlandscape)
# RDA model significant - df 4 - var 1.32 - F 2.14 - pval 0.011
# residual - df 23 - variance 3.68

# Individual axis significance
anova.cca(rda_fjordxlandscape, by = "axis")
# 1st axis marginally significant - df 1 - var 0.65 - F 4.26 - pval 0.062

# Individual term significance
anova.cca(rda_fjordxlandscape, by = "term")
# Max temp July significant - df 1 - var 0.44 - F 2.89 - pval 0.035
# Distance to sea significant - df 1 - var 0.4 - F 2.63 - pval 0.041

# Residual variation using Kaiser-Guttman criterion
rda_fjordxlandscape$CA$eig[rda_fjordxlandscape$CA$eig > mean(rda_fjordxlandscape$CA$eig)]
# unexplained variance PC1 1.89

#
## Landscape x local environment

locenvi <- full_join(locenvi_vege, locenvi_beetle)

# Redundancy analysis and plot
rda_landscapexlocenvi <- rda(subset(locenvi, select = -c(SiteID)) ~ ., data = subset(landscape, select = -c(SiteID)))
plotrda_landscapexlocenvi <- ggord(rda_landscapexlocenvi,
                                 text = 4,
                                 ptslab = TRUE,
                                 addsize = 3,
                                 size = 1,
                                 arrow = 0.3,
                                 repel = TRUE,
                                 vec_lab = list(Outfield_percent = "OutfArea", Infield_percent = "InfArea", TotForest_percent = "ForArea", Wetland_percent = "WetArea", TotCultivatedLand_percent = "CultArea"),
                                 xlims = c(-1.1, 1.1),
                                 ylims = c(-0.9, 1.2))
plotrda_landscapexlocenvi

# Summary rda
#summary(rda_landscapexlocenvi)
# 1st axis driven by forest cover (-0.72) and outfield cover (0.61)
# 2nd axis driven by wetland cover (0.79)
# 1st axis drives aspect (-0.62)
# 2nd axis drives bulk density (-0.48)

# Total variance explained by RDA
RsquareAdj(rda_landscapexlocenvi)
# R2 0.16
# Adj R2 -0.02

# Global RDA significance by permutation
anova.cca(rda_landscapexlocenvi)
# RDA model NS - df 5 - var 1.13 - F 0.89 - pval 0.66
# residual - df 23 - variance 5.87

# Individual axis significance
anova.cca(rda_landscapexlocenvi, by = "axis")
# 1st axis NS - df 1 - var 0.5 - F 1.94 - pval 0.72

# Individual term significance
anova.cca(rda_landscapexlocenvi, by = "term")
# NS

# Residual variation using Kaiser-Guttman criterion
rda_landscapexlocenvi$CA$eig[rda_landscapexlocenvi$CA$eig > mean(rda_landscapexlocenvi$CA$eig)]
# unexplained variance PC1 1.76 - PC2 1.19 - PC3 1.05 - PC4 0.85

#
## Fjord x grass

# Hellinger transformation of response matrix, suited for abundance data (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
hellinger_grass <- decostand(subset(grass, select = -c(SiteID)), method = "hellinger")

# Remove space within plant names (otherwise parse error in plot)
names(hellinger_grass) <- gsub(" ", "_", names(hellinger_grass))

# Redundancy analysis
rda_fjordxgrass <- rda(hellinger_grass ~ ., data = subset(fjordsys, select = -c(SiteID)))
plotrda_fjordxgrass <- ggord(rda_fjordxgrass,
                               text = 4,
                               #parse = FALSE,
                               ptslab = TRUE,
                               addsize = 3,
                               size = 1,
                             arrow = 0.3,
                             vec_lab = list(MaxTempJuly = "JulTemp", Elevation_max = "Elev", DistanceToSea_m = "SeaDist", AnnualPrecipitation = "AnnPreci"),
                               #repel = TRUE,
                               xlims = c(-1, 1.3),
                               ylims = c(-1.1, 1.3))
plotrda_fjordxgrass

# Summary rda
#summary(rda_fjordxgrass)
# 1st axis driven by July temperature (-0.9) and maximum elevation (0.31)
# 2nd axis driven by maximum elevation (0.9) and annual precipitation (0.58)
# 1st axis drives Poa pratensis (0.26)
# 2nd axis drives Holcus lanatus (0.23)

# Total variance explained by RDA
RsquareAdj(rda_fjordxgrass)
# R2 0.11
# Adj R2 -0.04

# Global RDA significance by permutation
anova.cca(rda_fjordxgrass)
# RDA model NS - df 4 - var 0.026 - F 0.76 - pval 0.74
# residual - df 24 - var 0.2

# Individual axis significance
anova.cca(rda_fjordxgrass, by = "axis") 
# 1st axis NS  - df 1 - var 0.014 - F 1.7 - pval 0.76

# Individual term significance
anova.cca(rda_fjordxgrass, by = "term") # NS

# Residual variation using Kaiser-Guttman criterion
rda_fjordxgrass$CA$eig[rda_fjordxgrass$CA$eig > mean(rda_fjordxgrass$CA$eig)]
# unexplained variance PC1 0.07 - PC2 0.05

#
## Grazing management x grass

# Hellinger transformation of response matrix, suited for abundance data (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
hellinger_grass <- decostand(subset(grass, select = -c(SiteID)), method = "hellinger")

# Remove space within plant names (otherwise parse error in plot)
names(hellinger_grass) <- gsub(" ", "_", names(hellinger_grass))

# Redundancy analysis
rda_grazingxgrass <- rda(hellinger_grass ~ ., data = subset(grazing, select = -c(SiteID, FlockSize1_adults)))
plotrda_grazingxgrass <- ggord(rda_grazingxgrass,
                               text = 4,
                               ptslab = TRUE,
                               addsize = 3,
                               size = 1,
                               arrow = 0.3,
                               vec_lab = list(TotalInfieldSurface = "FarmInfArea", GrazingSurface_ha = "GrazArea", AvgStockingDensity_perha = "StockDens", Sheep = "Sheep", Cow = "Cow"),
                               #repel = TRUE,
                               xlims = c(-1.1, 1.1),
                               ylims = c(-1.2, 1))
plotrda_grazingxgrass

# Summary rda
#summary(rda_grazingxgrass)
# 1st axis driven by sheep (-0.73), cow (0.76) and average stocking density (0.6)
# 2nd axis driven by grazing surface (-0.65), average stocking density (-0.61) and goats
# 1st axis drives Poa pratensis (-0.29) and Poa trivialis (0.29)
# 2nd axis drives Holcus lanatus (-0.25)

# Total variance explained by RDA
RsquareAdj(rda_grazingxgrass)
# R2 0.15
# Adj R2 -0.04

# Global RDA significance by permutation
anova.cca(rda_grazingxgrass)
# RDA model NS - df 5 - var 0.034 - F 0.79 - pval 0.77
# residual - df 23 - var 0.2

# Individual axis significance
anova.cca(rda_grazingxgrass, by = "axis") 
# 1st axis NS  - df 1 - var 0.02 - F 2.24 - pval 0.67

# Individual term significance
anova.cca(rda_grazingxgrass, by = "term") # NS

# Residual variation using Kaiser-Guttman criterion
rda_grazingxgrass$CA$eig[rda_grazingxgrass$CA$eig > mean(rda_grazingxgrass$CA$eig)]
# unexplained variance PC1 0.07 - PC2 0.04 - PC3 0.03

#
## Fine-scale environment x grass

# Hellinger transformation of response matrix, suited for abundance data (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
hellinger_grass <- decostand(subset(grass, select = -c(SiteID)), method = "hellinger")

# Remove space within plant names (otherwise parse error in plot)
names(hellinger_grass) <- gsub(" ", "_", names(hellinger_grass))

# Redundancy analysis
rda_locenvixgrass <- rda(hellinger_grass ~ ., data = subset(locenvi_vege, select = -c(SiteID)))
plotrda_locenvixgrass <- ggord(rda_locenvixgrass,
                               text = 4,
                               ptslab = TRUE,
                               addsize = 3,
                               size = 1,
                               arrow = 0.3,
                               vec_lab = list(GeneralSlope = "Slope", MeanBD = "BulkDens", MeanPT = "PenRate", MeanPhosphorus = "Phosp", AspectDegree = "Aspect", MeanpH = "pH"),
                               #repel = TRUE,
                               xlims = c(-1, 1),
                               ylims = c(-1, 0.8))
plotrda_locenvixgrass

# Summary RDA
#summary(rda_locenvixgrass)
# 1st axis driven by penetration rate (-0.58) and phosphorus (0.62)
# 2nd axis driven by phosphorus (-0.54) and penetration rate (-0.5)
# 1st axis drives Holcus lanatus (0.36), Deschampsia cespitosa (0.23), Deschampsia flexuosa (-0.26) and Poa trivialis (0.23)
# 2nd axis drives Deschampsia flexuosa (-0.26), Festuca rubra (0.29) and Poa pratensis (0.28)

# Total variance explained by RDA
RsquareAdj(rda_locenvixgrass)
# R2 0.26
# Adj R2 0.06

# Global RDA significance by permutation
anova.cca(rda_locenvixgrass)
# RDA model significant - df 6 - variance 0.061 - F 1.32 - pval 0.14
# residual df 22 - variance 0.17

# Individual axis significance
anova.cca(rda_locenvixgrass, by = "axis")
# 1st axis NS  - df 1 - var 0.03 - F 3.88 - pval 0.23

# Individual term significance
anova.cca(rda_locenvixgrass, by = "term")
# Mean Phosphorus significant - df 1 - variance 0.022 - F 2.92 - pval 0.027
# Mean PT marginally significant - df 1 - variance 0.017 - F 2.16 - pval 0.061

# Residual variation using Kaiser-Guttman criterion
rda_locenvixgrass$CA$eig[rda_locenvixgrass$CA$eig > mean(rda_locenvixgrass$CA$eig)]
# unexplained variance PC1 0.06 - PC2 0.05

#
## Fine-scale environment x forb

# Hellinger transformation of response matrix, suited for abundance data (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
hellinger_forb <- decostand(subset(forb, select = -c(SiteID)), method = "hellinger")

# Remove space within plant names (otherwise parse error in plot)
names(hellinger_forb) <- gsub(" ", "_", names(hellinger_forb))

# Redundancy analysis
rda_locenvixforb <- rda(hellinger_forb ~ ., data = subset(locenvi_vege, select = -c(SiteID, AspectDegree)))
plotrda_locenvixforb <- ggord(rda_locenvixforb,
                               text = 4,
                               ptslab = TRUE,
                               addsize = 3,
                               size = 1,
                              arrow = 0.3,
                              vec_lab = list(GeneralSlope = "Slope", MeanBD = "BulkDens", MeanPT = "PenRate", MeanPhosphorus = "Phosp", MeanpH = "pH"),
                               repel = TRUE,
                              xlims = c(-1, 1),
                              ylims = c(-1, 0.9))
plotrda_locenvixforb

# Summary RDA
#summary(rda_locenvixforb)
# 1st axis driven by penetration rate (0.73), phosphorus (-0.63), pH (-0.58) and bulk density (-0.47)
# 2nd axis driven by slope (-0.56), phosphorus (0.46), pH (0.44) and bulk density (-0.39)
# 1st axis drives Potentilla erecta (0.51), Trifolium repens (-0.38), Galium saxatile (0.34) and Ranunculus acris (0.3)
# 2nd axis drives Achillea millefolium (-0.38) and Galium saxatile (-0.24)

# Total variance explained by RDA
RsquareAdj(rda_locenvixforb)
# R2 0.3
# Adj R2 0.15

# Global RDA significance by permutation
anova.cca(rda_locenvixforb)
# RDA model significant - df 5 - var 0.14 - F 1.96 - pval 0.005
# residual df 23 - variance 0.33

# Individual axis significance
anova.cca(rda_locenvixforb, by = "axis")
# 1st axis significant - df 1 - var 0.08 - F 5.84 - pval 0.007

# Individual term significance
anova.cca(rda_locenvixforb, by = "term")
# Mean PT significant - df 1 - var 0.035 - F 2.43 - pval 0.028
# Mean Phosphorus significant - df 1 - variance 0.035 - F 2.45 - pval 0.038
# Mean BD marginally significant - df 1 - variance 0.03 - F 2.08 - pval 0.07

# Residual variation using Kaiser-Guttman criterion
rda_locenvixforb$CA$eig[rda_locenvixforb$CA$eig > mean(rda_locenvixforb$CA$eig)]
# unexplained variance PC1 0.11 - PC2 0.08 - PC3 0.05

#
## Fine-scale environment x beetle

# Hellinger transformation of response matrix, suited for abundance data (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
hellinger_beetle <- decostand(subset(beetle, select = -c(SiteID)), method = "hellinger")

# Redundancy analysis
rda_locenvixbeetle <- rda(hellinger_beetle ~ ., data = subset(locenvi_beetle, select = -c(SiteID)))
plotrda_locenvixbeetle <- ggord(rda_locenvixbeetle,
                              text = 4,
                              #parse = FALSE,
                              ptslab = TRUE,
                              addsize = 3,
                              size = 1,
                              arrow = 0.3,
                              vec_lab = list(GeneralSlope = "Slope", AspectDegree = "Aspect", MeanBD = "BulkDens", MeanPT = "PenRate", MeanLitter = "LitCov", MeanBryo = "BryoCov", MeanHeight = "AvgHeight"),
                              #repel = TRUE,
                              xlims = c(-1.2, 1),
                              ylims = c(-1.1, 1.1))
plotrda_locenvixbeetle

# Summary rda
#summary(rda_locenvixbeetle)
# 1st axis driven by bulk density (-0.69), aspect (0.43) and vegetation height (0.4)
# 2nd axis driven by bryophyte layer (0.71), vegetation height (-0.76) and slope (0.4)
# 1st axis drives Carabidae (-0.34) and Ptilidae (0.35)
# 2nd axis drives Scarabaeidae (0.21)

# Total variance explained by RDA
RsquareAdj(rda_locenvixbeetle)
# R2 0.31
# Adj R2 0.08

# Global RDA significance by permutation
anova.cca(rda_locenvixbeetle)
# RDA model NS - df 7 - var 0.013 - F 1.34 - pval 0.20
# residual df 21 - variance 0.03

# Individual axis significance
anova.cca(rda_locenvixbeetle, by = "axis")
# 1st axis NS  - df 1 - var 0.01 - F 7.93 - pval 0.16

# Individual term significance
anova.cca(rda_locenvixbeetle, by = "term")
# Mean BD significant - df 1 - variance 0.005 - F 3.63 - pval 0.019
# Mean Height marginally significant - df 1 - variance 0.004 - F 2.47 - pval 0.094

# Residual variation using Kaiser-Guttman criterion
rda_locenvixbeetle$CA$eig[rda_locenvixbeetle$CA$eig > mean(rda_locenvixbeetle$CA$eig)]
# unexplained variance PC1 0.018

#
## Arrangement all plots

rdall <- ggarrange(plotrda_fjordxlandscape, plotrda_landscapexlocenvi, plotrda_grazingxgrass, plotrda_locenvixgrass, plotrda_locenvixforb, plotrda_locenvixbeetle, 
                   labels = c("A", "B", "C", "D", "E", "F"),
                   ncol = 2,
                   nrow = 3)
rdall
ggsave("outputs/RDAresults.png", plot = rdall, width = 19, height = 26, units = "cm", bg = "white")