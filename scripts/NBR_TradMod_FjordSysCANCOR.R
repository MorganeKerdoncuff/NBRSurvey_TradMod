# TradMod WP2 - NBR survey grazer loading script
#Description of the data
#Year - 2023
#Who - Morgane Kerdoncuff
#Project - TradMod
#Funding - NFR
#Place - University of Bergen, Norway

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
siteinfo_infield <- siteinfo_full |>  
  filter(Habitat == "permanent grassland" & SiteID != "IG3") |> 
  dplyr::select(SiteID, Livestock)

# Extraction in other datasets
climate_infield <- filter(climate_full, SiteID %in% siteinfo_infield$SiteID)
landscape_infield <- filter(landscape_full, SiteID %in% siteinfo_infield$SiteID)
landuse_infield <- filter(landuse_full, SiteID %in% siteinfo_infield$SiteID)
area20x20_infield <- filter(area20x20_full, SiteID %in% siteinfo_infield$SiteID)
groundcover_infield <- filter(groundcover_full, SiteID %in% siteinfo_infield$SiteID)
soilbulk_infield <- filter(soilbulk_full, SiteID %in% siteinfo_infield$SiteID)
soilchem_infield <- filter(soilchem_full, SiteID %in% siteinfo_infield$SiteID)
soilpene_infield <- filter(soilpene_full, SiteID %in% siteinfo_infield$SiteID)
vege_infield <- filter(vege_full, SiteID %in% siteinfo_infield$SiteID)
beetle_infield <- filter(beetle_full, SiteID %in% siteinfo_infield$SiteID)

#
## Summarise data at site level -> should be 29 observations for non community data

# Site info - validated

# Climate - validated

# Landscape matrix - validated

# Field management - validated

# Area 20 x 20 - validated

# Ground cover - current at sample level -> summary by average
groundcover_infield <- groundcover_infield |>  
  group_by(SiteID) |>  
  summarise(BareGrd = mean(Bare_soil, na.rm=TRUE),
            Rocks = mean(Rocks, na.rm=TRUE),
            Litter = mean(Litter, na.rm=TRUE),
            DeadWood = mean(Dead_wood, na.rm=TRUE),
            Bryo = mean(Bryophytes, na.rm=TRUE),
            Lichen = mean(Lichen, na.rm=TRUE),
            Vasc = mean(Vascular, na.rm=TRUE),
            Blossom = mean(Blossom_cover, na.rm=TRUE),
            MeanHeight = mean(VG_mean_height_cm, na.rm=TRUE),
            MaxHeight = mean(VG_max_height_cm, na.rm=TRUE),
            PlantRich = mean(Plant_species_richness, na.rm=TRUE))

# Bulk density - current at sample level -> summary by average
soilbulk_infield <- soilbulk_infield |> 
  group_by(SiteID) |> 
  summarise(BD = mean(BD),
            GWC = mean(GWC_48))

# Soil chemistry - current at plot level -> summary by average
soilchem_infield <- soilchem_infield |> 
  group_by(SiteID) |> 
  summarise(LOI = mean(LOI),
            Humus = mean(Humus_percentDM),
            pH = mean(pH),
            Phosph = mean(P.Al_mg.100g),
            Potass = mean(K.Al_mg.100g),
            Nitrog = mean(TotalN_percentDM),
            Sodium = mean(Na.Al_mg.100g))

# Soil penetration - current at sample level -> summary by average
soilpene_infield <- soilpene_infield |> 
  group_by(SiteID) |> 
  summarise(SoilPene = mean(AveragePT_cm))

# Plant community - current at quadrat level -> summary by average
vege_infield <- vege_infield |> 
  group_by(SiteID, Species) |> 
  summarise(PlantSp_cover = mean(Abundance))

# Beetle community - current at pitfall level -> summary by average
beetle_infield <- beetle_infield |> 
  group_by(SiteID, BeetleFamilies) |> 
  summarise(BeetleFam_abundance = mean(BeetleFam_abundance))

#
## Selection & transformation plant community data

# Data distribution
#hist(vege_grass$PlantSp_cover) # Poisson, highly skewed

# Species frequencies across sites
vege_freq <- filter(vege_infield, PlantSp_cover>0) |>
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
vege_average <- vege_infield |>
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

#
## Selection and transformation grass assemblage

# Selection main grass species - at least present in 10 sites AND min average cover 1% -> 8 species
grass <- subset(vege_infield,
                Species == "Agrostis capillaris" |
                Species == "Festuca rubra" | 
                Species == "Holcus lanatus" | 
                Species == "Poa pratensis" | 
                Species == "Deschampsia cespitosa" | 
                Species == "Anthoxantum odoratum" | 
                Species == "Deschampsia flexuosa" | 
                Species == "Poa trivialis")

# Hellinger transformation on contingency table (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
contin_grass <- xtabs(formula = PlantSp_cover ~ SiteID + Species, data = grass)
contin_grass <- decostand(contin_grass, method = "hellinger")

# Wide table
grass <- as.data.frame(contin_grass)
grass <- grass |> 
  pivot_wider(names_from = Species, values_from = Freq)

#
## Selection and transformation forb assemblage

# Selection main forb species - at least present in 10 sites AND min average cover 1% -> 7 species
forb <- subset(vege_infield,
                Species == "Trifolium repens" |
                  Species == "Rumex acetosa" |
                  Species == "Galium saxatile" |
                  Species == "Potentilla erecta" |
                  Species == "Ranunculus acris" |
                  Species == "Ranunculus repens" |
                  Species == "Achillea millefolium")

# Hellinger transformation on contingency table (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
contin_forb <- xtabs(formula = PlantSp_cover ~ SiteID + Species, data = forb)
contin_forb <- decostand(contin_forb, method = "hellinger")

# Wide table
forb <- as.data.frame(contin_forb)
forb <- forb |> 
  pivot_wider(names_from = Species, values_from = Freq)

#
## Transformation beetle assemblage data

#hist(arthro_grass$BeetleFam_abundance) # Poisson, highly skewed

# Family frequencies across sites
beetle_freq <- filter(beetle_infield, BeetleFam_abundance>0) |>
  group_by(BeetleFamilies) |>
  count() |>
  dplyr::arrange(desc(n)) # filtering by frequency
filter(beetle_freq, n > 10)
# 8 beetle families in at least 10 sites Carab; Hydro; Scara; Staph; Ptili; Curcu; Elat; Silph

# Species average percent cover across site
beetle_average <- beetle_infield |>
  group_by(BeetleFamilies) |>
  summarise_if(is.numeric, mean, na.rm = TRUE) |>
  dplyr::arrange(desc(BeetleFam_abundance)) # filtering by average value across sites
filter(beetle_average, BeetleFam_abundance > 3)
# 5 beetles families with at least 3 individuals on average Staph; Ptili; Hydro; Scara; Carab

# Selection main dung beetle families -> at least present in 10 sites + min 3 individuals on average -> 6 families
beetle <- subset(beetle_infield,
                 BeetleFamilies == "Carabidae" |
                   BeetleFamilies == "Staphylinidae" |
                   BeetleFamilies == "Hydrophilidae" |
                   BeetleFamilies == "Ptiliidae" |
                   BeetleFamilies == "Scarabaeidae")

# Log transformation
#beetle <- beetle |> 
#  mutate(BeetleFam_logabundance = log1p(BeetleFam_abundance))

# Hellinger transformation on contingency table (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
contin_beetle <- xtabs(formula = BeetleFam_abundance ~ SiteID + BeetleFamilies, data = beetle)
contin_beetle <- decostand(contin_beetle, method = "hellinger")

# Wide table
beetle <- as.data.frame(contin_beetle)
beetle <- beetle |> 
  pivot_wider(names_from = BeetleFamilies, values_from = Freq)

#
## Creation fjord system matrix

# Desired variables
## Elevation (num), from area 20x20
## Annual precipitation (num), from climate
## Average July temp (num), from climate
## Average Jan temp (num), from climate
## Distance to sea (num), from area 20x20

# Selection variables
fjordsys <- full_join(area20x20_infield, climate_infield)
fjordsys <- subset(fjordsys, select = c(SiteID, Elevation, annualprecipitation, maxtempJuly, mintempJan, DistanceToSea_m))

# Variable renaming
names(fjordsys) <- gsub("Elevation", "Elev", names(fjordsys))
names(fjordsys) <- gsub("annualprecipitation", "AnnPreci", names(fjordsys))
names(fjordsys) <- gsub("maxtempJuly", "JulTemp", names(fjordsys))
names(fjordsys) <- gsub("mintempJan", "JanTemp", names(fjordsys))
names(fjordsys) <- gsub("DistanceToSea_m", "SeaDist", names(fjordsys))

# Check distribution
#hist(fjordsys$Elev) # Poisson distribution, not skewed -> validated
#hist(fjordsys$AnnPreci) # Normal-like distribution, with some gaps -> validated
#hist(fjordsys$JulTemp) # Normal-like distribution -> validated
#hist(fjordsys$JanTemp) # Reverse Poisson distribution -> validated
#hist(fjordsys$SeaDist) # Poisson-like distribution, with one gap -> validated

# Scaling
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
landscape_infield <- landscape_infield |> 
  mutate(Cultivated = FullyCultivatedLand_percent + SuperficiallyCultivatedLand_percent) |> 
  mutate(Forest = ProductiveForest_percent + NonProductiveForest_percent)

# Selection variables
landscape <- subset(landscape_infield, select = c(SiteID, Cultivated, Forest, Infield_percent, Outfield_percent, Wetland_percent))

# Variable renaming
names(landscape) <- gsub("Infield_percent", "Infield", names(landscape))
names(landscape) <- gsub("Outfield_percent", "Outfield", names(landscape))
names(landscape) <- gsub("Wetland_percent", "Wetland", names(landscape))

# Check distribution
#hist(landscape$Cultivated) # Poisson-like distribution -> validated
#hist(landscape$Forest) # Normal-like distribution, one gap in the middle -> validated
#hist(landscape$Infield) # unbalanced normal-like, but no outlier -> validated
#hist(landscape$Outfield) # Poisson, not skewed -> validated
#hist(landscape$Wetland) # Poisson, not skewed -> validated

# Scaling
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
landuse_infield <- landuse_infield |> 
  mutate(Sheep = ifelse(Livestock1 == "sheep",1,0)) |> 
  mutate(Cow = ifelse(Livestock1 == "cow", 1,0))

# Selection variables
grazing <- subset(landuse_infield, select = c(SiteID, Sheep, Cow, FlockSize1_adults, SelectedFieldArea_ha, FarmInfieldArea_ha, AvgStockDensity_perha))

# Variable renaming
names(grazing) <- gsub("FlockSize1_adults", "NbAdults", names(grazing))
names(grazing) <- gsub("SelectedFieldArea_ha", "FieldSize", names(grazing))
names(grazing) <- gsub("FarmInfieldArea_ha", "FarmInfield", names(grazing))
names(grazing) <- gsub("AvgStockDensity_perha", "StockDens", names(grazing))

# Check distribution
#hist(grazing$NbAdults) # one outlier above 150 -> rejected
#hist(grazing$FieldSize) # Poisson distribution, with few gaps -> validated
#hist(grazing$FarmInfield) # Poisson-like distribution, a bit unbalanced but no outlier
#hist(grazing$StockDens) # Poisson distribution -> validated

# Scaling
grazing <- subset(grazing, select = -c(NbAdults))
grazing <- grazing |> 
  mutate(across(where(is.numeric), scale))

#
## Creation general fine-scale environment matrix

# Desired variables
## 

# Selection variables
locenvi <- purrr::reduce(list(groundcover_infield, soilbulk_infield, soilpene_infield, soilchem_infield, area20x20_infield), dplyr::left_join)
locenvi <- subset(locenvi, select = c(SiteID, Litter, Bryo, MeanHeight, BD, GWC, LOI, Nitrog, Phosph, pH, Humus, SoilPene, Aspect_degree, Slope_degree))

# Variables renaming
names(locenvi) <- gsub("Aspect_degree", "Aspect", names(locenvi))
names(locenvi) <- gsub("Slope_degree", "Slope", names(locenvi))

# Check distribution
#hist(locenvi$Litter) # Highly skewed Poisson, very small range with 2 outliers -> rejected
#hist(locenvi$Bryo) # Unbalanced Normal-like but no outlier -> validated
#hist(locenvi$MeanHeight) # Poisson like, no outlier -> validated
#hist(locenvi$BD) # Normal-like, a bit unbalanced but no outlier -> validated
#hist(locenvi$GWC) # poisson like, one gap between 50 & 55 -> validated
#hist(locenvi$LOI) # Poisson like, a bit skewed but no outlier -> validated
#hist(locenvi$Nitrog) # Poisson like, with a gap but no outlier -> validated
#hist(locenvi$Phosph) # poisson like, a bit skewed but no outlier -> validated
#hist(locenvi$pH) # good Normal distribution -> validated
#hist(locenvi$Humus) # Poisson distribution, a bit skewed -> validated
#hist(locenvi$SoilPene) # Normal distribution a bit unbalanced, no outlier -> validated
#hist(locenvi$Aspect) # inverse Normal distribution, but no outlier -> validated
#hist(locenvi$Slope) # Poisson distribution, one gap but no outlier -> validated

# Data scaling
locenvi <- subset(locenvi, select = -c(Litter))
locenvi <- locenvi |> 
  mutate(across(where(is.numeric), scale))

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
pairs(select_if(fjordsys, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # high colinearity between max temp and min temp

# Landscape variables
pairs(select_if(landscape, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # one correlation a bit strong (0.74) between total cultivated land and infield percent

# Grazing variables
pairs(select_if(grazing, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # no correlation stronger than 0.64 -> validated

# Local environmental variables
pairs(select_if(locenvi, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # 5 variables highly correlated (>0.8): BD, moisture, LOI, humus, Nitrogen

# Removal strongly correlated variables - BD kept as with higher number of replication for each site
fjordsys <- subset(fjordsys, select = -c(JanTemp))
locenvi <- subset(locenvi, select = -c(GWC, LOI, Humus, Nitrog))

#
## Data preparation - Contingency tables

# Fjord system
fjordsys_long <- fjordsys |> 
  pivot_longer(
    cols = c(Elev, AnnPreci, JulTemp, SeaDist),
    names_to = "Factors",
    values_to = "Values")
contin_fjordsys <- xtabs(formula = Values ~ SiteID + Factors, data = fjordsys_long)

# Landscape
landscape_long <- landscape |> 
  pivot_longer(
    cols = c(Cultivated, Forest, Infield, Outfield, Wetland),
    names_to = "Factors",
    values_to = "Values"
  )
contin_landscape <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long)

# Grazing
grazing_long <- grazing |> 
  pivot_longer(
    cols = c(Sheep, Cow, FieldSize, FarmInfield, StockDens),
    names_to = "Factors",
    values_to = "Values")
contin_grazing <- xtabs(formula = Values ~ SiteID + Factors, data = grazing_long)

# Local environment
locenvi_long <- locenvi |> 
  pivot_longer(
    cols = c(Bryo, MeanHeight, BD, Phosph, pH, SoilPene, Aspect, Slope),
    names_to = "Factors",
    values_to = "Values")
contin_locenvi <- xtabs(formula = Values ~ SiteID + Factors, data = locenvi_long)



#### Canonical correlation analysis between driver sets ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape <- length(select_if(landscape, is.numeric))
nvar_grazing <- length(select_if(grazing, is.numeric))
nvar_locenvi <- length(select_if(locenvi, is.numeric))

#
## Fjord effect

# Fjord x landscape
cancor <- cc(contin_fjordsys, contin_landscape)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.79
# 2nd axis correlation 0.75
# 3rd axis correlation 0.55
test <- p.asym(rhocancor, nobs, nvar_fjordsys, nvar_landscape, tstat = "Hotelling")
# 1st dim significant - stat 3.45 - df1 20 - df2 74 - pval 1.5.10-4
# 2nd dim significant - stat 1.77 - df1 12 - df2 82 - pval 1.5.10-3
# 3rd dim marginally significant - stat 0.5 - df1 6 - df2 90 - pval 0.09
fjordlandscape <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fjord", response = "Landscape", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
fjordlandscape

# Fjord x grazing
cancor <- cc(contin_fjordsys, contin_grazing)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.61
test <- p.asym(rhocancor, nobs, nvar_fjordsys, nvar_grazing, tstat = "Hotelling")
# 1st dim NS - stat 1.11 - df1 24 - df2 70 - pval 0.72
fjordgrazing <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fjord", response = "Grazing", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
fjordgrazing

# Fjord x local environment
cancor <- cc(contin_fjordsys, contin_locenvi)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.63
test <- p.asym(rhocancor, nobs, nvar_fjordsys, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 1.43 - df1 36 - df2 58 - pval 0.96
fjordlocenvi <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fjord", response = "Finescale", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
fjordlocenvi

#
## Landscape effect

# Landscape x grazing
cancor <- cc(contin_landscape, contin_grazing)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.69
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_grazing, tstat = "Hotelling")
# 1st dim NS - stat 1.96 - df1 30 - df2 82 - pval 0.39
landscapegrazing <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Grazing", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
landscapegrazing

# Landscape x local environment
cancor <- cc(contin_landscape, contin_locenvi)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.80
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 3.64 - df1 45 - df2 67 - pval 0.38
landscapelocenvi <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Finescale", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
landscapelocenvi

#
## Grazing effect

# Grazing x local environment
cancor <- cc(contin_grazing, contin_locenvi)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.79
test <- p.asym(rhocancor, nobs, nvar_grazing, nvar_locenvi, tstat = "Hotelling")
# 1st dim NS - stat 4.26 - df1 54 - df2 74 - pval 0.54
grazinglocenvi <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Grazing", response = "Finescale", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
grazinglocenvi

#
## Merging output

envar <- purrr::reduce(list(fjordlandscape, fjordgrazing, fjordlocenvi, landscapegrazing, landscapelocenvi, grazinglocenvi), dplyr::full_join)
envar <- envar |> 
  filter(dim == "dim1" | p.value < 0.1)


#### Canonical correlation analysis grass community ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape <- length(select_if(landscape, is.numeric))
nvar_grazing <- length(select_if(grazing, is.numeric))
nvar_locenvi <- length(select_if(locenvi, is.numeric))
nvar_grass <- length(select_if(grass, is.numeric))

# Fjord x grass community
cancor <- cc(contin_fjordsys, contin_grass)
rhocancor <- cancor$cor
rhocancor 
# 1st axis correlation 0.8
test <- p.asym(rhocancor, nobs, nvar_fjordsys, nvar_grass, tstat = "Hotelling") 
# 1st dim marginally significant - stat 2.59 - df1 28 - df2 66 - pval 0.081
fjordgrass <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fjord", response = "Grass", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
fjordgrass

# Landscape x grass community
cancor <- cc(contin_landscape, contin_grass)
rhocancor <- cancor$cor
rhocancor 
# 1st axis correlation 0.82
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_grass, tstat = "Hotelling") 
# 1st dim NS - stat 3.15 - df1 35 - df2 77 - pval 0.12
landscapegrass <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Grass", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
landscapegrass

# Grazing x grass community
cancor <- cc(contin_grazing, contin_grass)
rhocancor <- cancor$cor
rhocancor 
# 1st axis correlation 0.87 
test <- p.asym(rhocancor, nobs, nvar_grazing, nvar_grass, tstat = "Hotelling") 
# 1st axis significant - stat 4.6 - df1 35 - df2 77 - pval 5.2.10-3
grazinggrass <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Grazing", response = "Grass", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
grazinggrass

# Local environment effect x grass community
cancor <- cc(contin_locenvi, contin_grass)
rhocancor <- cancor$cor
rhocancor 
# 1st axis correlation 0.87
test <- p.asym(rhocancor, nobs, nvar_locenvi, nvar_grass, tstat = "Hotelling")
# 1st axis significant - stat 4.63 - df1 42 - df2 86 - pval 0.038
locenvigrass <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Finescale", response = "Grass", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
locenvigrass

# Merging output
grassvar <- purrr::reduce(list(fjordgrass, landscapegrass, grazinggrass, locenvigrass), dplyr::full_join)
grassvar <- grassvar |> 
  filter(dim == "dim1" | p.value < 0.1)


#### Canonical correlation analysis forb community ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape <- length(select_if(landscape, is.numeric))
nvar_grazing <- length(select_if(grazing, is.numeric))
nvar_locenvi <- length(select_if(locenvi, is.numeric))
nvar_forb <- length(select_if(forb, is.numeric))


# Fjord x forb community
cancor <- cc(contin_fjordsys, contin_forb)
rhocancor <- cancor$cor
rhocancor 
# 1st axis correlation 0.69
test <- p.asym(rhocancor, nobs, nvar_fjordsys, nvar_forb, tstat = "Hotelling") 
# 1st dim NS - stat 1.57 - df1 28 - df2 66 - pval 0.58
fjordforb <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fjord", response = "Forb", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
fjordforb

# Landscape x forb community
cancor <- cc(contin_landscape, contin_forb)
rhocancor <- cancor$cor
rhocancor 
# 1st axis correlation 0.75
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_forb, tstat = "Hotelling") 
# 1st dim NS - stat 3 - df1 35 - df2 77 - pval 0.15
landscapeforb <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Forb", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
landscapeforb

# Grazing x forb community
cancor <- cc(contin_grazing, contin_forb)
rhocancor <- cancor$cor
rhocancor 
# 1st axis correlation 0.79
test <- p.asym(rhocancor, nobs, nvar_grazing, nvar_forb, tstat = "Hotelling") 
# 1st dim NS - stat 2.89 - df1 35 - df2 77 - pval 0.19
grazingforb <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Grazing", response = "Forb", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
grazingforb

# Local environment effectx forb
cancor <- cc(contin_locenvi, contin_forb)
rhocancor <- cancor$cor
rhocancor 
# 1st axis correlation 0.77
test <- p.asym(rhocancor, nobs, nvar_locenvi, nvar_forb, tstat = "Hotelling") 
# 1st dim NS - stat 4.3 - df1 56 - df2 86 - pval 0.6
locenviforb <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Finescale", response = "Forb", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
locenviforb

# Merging output
forbvar <- purrr::reduce(list(fjordforb, landscapeforb, grazingforb, locenviforb), dplyr::full_join)
forbvar <- forbvar |> 
  filter(dim == "dim1" | p.value < 0.1)


#### Canonical correlation analysis beetle assemblage ####

# Number of observations (same for all sets)
nobs <- dim(contin_fjordsys)[1]

# Number of variables in each set
nvar_fjordsys <- length(select_if(fjordsys, is.numeric))
nvar_landscape <- length(select_if(landscape, is.numeric))
nvar_grazing <- length(select_if(grazing, is.numeric))
nvar_locenvi <- length(select_if(locenvi, is.numeric))
nvar_beetle <- length(select_if(beetle, is.numeric))

# Fjord x beetle assemblage
cancor <- cc(contin_fjordsys, contin_beetle)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.69
test <- p.asym(rhocancor, nobs, nvar_fjordsys, nvar_beetle, tstat = "Hotelling")
# 1st dim NS - stat 1.38 - df1 20 - df2 74 - pval 0.22
fjordbeetle <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fjord", response = "Beetle", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
fjordbeetle

# Landscape x beetle assemblage
cancor <- cc(contin_landscape, contin_beetle)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.71
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_beetle, tstat = "Hotelling")
# 1st dim NS - stat 1.76 - df1 25 - df2 87 - pval 0.24
landscapebeetle <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Beetle", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
landscapebeetle

# Grazing x beetle community
cancor <- cc(contin_grazing, contin_beetle)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.72
test <- p.asym(rhocancor, nobs, nvar_grazing, nvar_beetle, tstat = "Hotelling")
# 1st dim NS - stat 1.64 - df1 25 - df2 87 - pval 0.32
grazingbeetle <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Grazing", response = "Beetle", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
grazingbeetle

# Local environment effect
cancor <- cc(contin_locenvi, contin_beetle)
rhocancor <- cancor$cor
rhocancor
# 1st axis correlation 0.73
test <- p.asym(rhocancor, nobs, nvar_locenvi, nvar_beetle, tstat = "Hotelling")
# 1st dim NS - stat 2.5 - df1 30 - df2 82 - pval 0.14
locenvibeetle <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Finescale", response = "Beetle", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
locenvibeetle

# Merging output
beetlevar <- purrr::reduce(list(fjordbeetle, landscapebeetle, grazingbeetle, locenvibeetle), dplyr::full_join)
beetlevar <- beetlevar |> 
  filter(dim == "dim1" | p.value < 0.1)

# Merge all community outputs
comvar <- purrr::reduce(list(grassvar, forbvar, beetlevar), dplyr::full_join)


#### RDA PLOTTING FUNCTION ####

# Modification ggord function so that axes show non-cumulated explained variance proportion
ggrda <- function (ord_in, grp_in = NULL, axes = c("1", "2"), ...) 
{
  axes <- paste0("RDA", axes)
  obs <- data.frame(ord_in$CCA$wa[, axes])
  obs$Groups <- grp_in
  addpts <- data.frame(ord_in$CCA$v[, axes])
  constr <- data.frame(ord_in$CCA$biplot[, axes])
  exp_var <- summary(ord_in)$cont$importance[2, axes]
  axes <- paste0(axes, " (", round(100 * exp_var, 2), "%)")
  names(obs)[1:2] <- axes
  ggord:::ggord.default(obs, vecs = constr, axes, addpts = addpts, 
                ...)
}

#
## Example regular plotting

# Extraction axis significance
#perc <- round(100*(summary(rda)$cont$importance[2, 1:2]), 2)

# # Extraction scores
# sc_si <- scores(rda, display = "sites", choices = c(1,2), scaling=2)
# sc_sp <- scores(rda, display = "species", choices = c(1,2), scaling=2)
# sc_bp <- scores(rda, display = "bp", choices = c(1, 2), scaling=2)
# 
# 
# plot(rda,
#      scaling = 2, # set scaling type 
#      type = "none", # this excludes the plotting of any points from the results
#      frame = FALSE,
#      # set axis limits
#      #xlim = c(-1,1), 
#      #ylim = c(-1,1),
#      # label the plot (title, and axes)
#      main = "",
#      xlab = paste0("RDA1 (", perc[1], "%)"), 
#      ylab = paste0("RDA2 (", perc[2], "%)") 
# )
# # add points for site scores
# points(sc_si, 
#        #pch = 21, # set shape (here, circle with a fill colour)
#        #col = "black", # outline colour
#        #bg = "steelblue", # fill colour
#        cex = 0.8) # size
# # add points for species scores
# points(sc_sp, 
#        pch = 22, # set shape (here, square with a fill colour)
#        col = "black",
#        bg = "chartreuse4", 
#        cex = 1.2)
# # add text labels for species abbreviations
# text(sc_sp + c(0.1, 0.1), # adjust text coordinates to avoid overlap with points 
#      labels = rownames(sc_sp), 
#      col = "chartreuse4", 
#      font = 1, # bold
#      cex = 0.6)
# # add arrows for effects of the explanatory variables
# arrows(0,0, # start them from (0,0)
#        sc_bp[,1]*2, sc_bp[,2]*2, # end them at the score value
#        col = "cyan4", 
#        lwd = 2)
# # add text labels for arrows
# text(x = sc_bp[,1]*2 -0.1, # adjust text coordinate to avoid overlap with arrow tip
#      y = sc_bp[,2]*2 - 0.1, 
#      labels = rownames(sc_bp), 
#      col = "cyan4", 
#      cex = 0.8, 
#      font = 2)



#### Redundancy analysis on environmental drivers ####

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
rda <- rda(subset(landscape, select = -c(SiteID)) ~ ., data = subset(fjordsys, select = -c(SiteID)))
plotrda_fjordxlandscape <- ggrda(rda,
      txt = 4,
      ptslab = TRUE,
      addcol = "chartreuse4",
      addsize = 3,
      size = 1,
      arrow = 0.3,
      #repel = TRUE,
      veccol = "cyan4",
      labcol = "cyan4",
      xlims = c(-1.3, 1),
      ylims = c(-0.8, 1.4))
plotrda_fjordxlandscape

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |>
  # "variance" output from cca is actually eigenval
  mutate(eigenval = Variance) |>
  # calculate proportion explained variance from eigvenval and inertia
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance in comparison with all other terms at same level (marginal)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafjordlandscape <- rdatable |> 
  mutate(model = "FjordxLandscape", explanatory = "Fjord", response = "Landscape") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Fjord x grazing (CANCOR NS)

# Redundancy analysis and plot
rda <- rda(subset(grazing, select = -c(SiteID)) ~ ., data = subset(fjordsys, select = -c(SiteID)))
plotrda_fjordxgrazing <- ggrda(rda,
                               txt = 4,
                               ptslab = TRUE,
                               addcol = "darkgoldenrod3",
                               addsize = 3,
                               size = 1,
                               arrow = 0.3,
                               #repel = TRUE,
                               veccol = "cyan4",
                               labcol = "cyan4",
                                xlims = c(-1.5, 1.3),
                                ylims = c(-1.3, 1.6))
plotrda_fjordxgrazing

# Summary rda
#summary(rda_fjordxgrazing)
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafjordgrazing <- rdatable |>
  mutate(model = "FjordxGrazing", explanatory = "Fjord", response = "Grazing") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Fjord x locenvi (CANCOR NS)

# Redundancy analysis and plot
rda <- rda(subset(locenvi, select = -c(SiteID)) ~ ., data = subset(fjordsys, select = -c(SiteID)))
plotrda_fjordxlocenvi <- ggrda(rda,
                               txt = 4,
                               ptslab = TRUE,
                               addcol = "darkred",
                               addsize = 3,
                               size = 1,
                               arrow = 0.3,
                               #repel = TRUE,
                               veccol = "cyan4",
                               labcol = "cyan4",
                               xlims = c(-1.1, 0.9),
                               ylims = c(-0.9, 1.1))
plotrda_fjordxlocenvi

# Summary rda
#summary(rda_fjordxlocenvi)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafjordlocenvi <- rdatable |> 
  mutate(model = "FjordxFinescale", explanatory = "Fjord", response = "Finescale") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Landscape x grazing

# Redundancy analysis and plot
rda <- rda(subset(grazing, select = -c(SiteID)) ~ ., data = subset(landscape, select = -c(SiteID)))
plotrda_landscapexgrazing <- ggrda(rda,
                                   txt = 4,
                                   ptslab = TRUE,
                                   addcol = "darkgoldenrod3",
                                   addsize = 3,
                                   size = 1,
                                   arrow = 0.3,
                                   repel = TRUE,
                                   veccol = "chartreuse4",
                                   labcol = "chartreuse4",
                                   xlims = c(-1.2, 1.1),
                                   ylims = c(-1.2, 1.1))
plotrda_landscapexgrazing

# Summary rda
#summary(rda_landscapexgrazing)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapegrazing <- rdatable |> 
  mutate(model = "LandscapexGrazing", explanatory = "Landscape", response = "Grazing") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Landscape x local environment

# Redundancy analysis and plot
rda <- rda(subset(locenvi, select = -c(SiteID)) ~ ., data = subset(landscape, select = -c(SiteID)))
plotrda_landscapexlocenvi <- ggrda(rda,
                                  txt = 4,
                                   ptslab = TRUE,
                                   addcol = "darkred",
                                   addsize = 3,
                                   size = 1,
                                   arrow = 0.3,
                                   repel = TRUE,
                                   veccol = "chartreuse4",
                                   labcol = "chartreuse4",
                                 xlims = c(-1.1, 1.1),
                                 ylims = c(-0.9, 1.2))
plotrda_landscapexlocenvi

# Summary rda
#summary(rda_landscapexlocenvi)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapelocenvi <- rdatable |> 
  mutate(model = "LandscapexFinescale", explanatory = "Landscape", response = "Finescale") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Grazing x local environment (CANCOR NS)

# Redundancy analysis and plot
rda <- rda(subset(locenvi, select = -c(SiteID)) ~ ., data = subset(grazing, select = -c(SiteID)))
plotrda_grazingxlocenvi <- ggrda(rda,
                                 txt = 4,
                                 ptslab = TRUE,
                                 addcol = "darkred",
                                 addsize = 3,
                                 size = 1,
                                 arrow = 0.3,
                                 repel = TRUE,
                                 veccol = "darkgoldenrod3",
                                 labcol = "darkgoldenrod3",
                                   xlims = c(-1, 1),
                                   ylims = c(-1.1, 0.9))
plotrda_grazingxlocenvi

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdagrazinglocenvi <- rdatable |>
  mutate(model = "GrazingxFinescale", explanatory = "Grazing", response = "Finescale") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Environmental data summary

# Plot Grouping
rda_envi <- ggarrange(plotrda_fjordxlandscape, plotrda_fjordxgrazing, plotrda_fjordxlocenvi, plotrda_landscapexgrazing, plotrda_landscapexlocenvi, plotrda_grazingxlocenvi, 
                   labels = c("A", "B", "C", "D", "E", "F"),
                   ncol = 2,
                   nrow = 3)
rda_envi
ggsave("outputs/RDA_envi.png", plot = rda_envi, width = 19, height = 26, units = "cm", bg = "white")


# Stat summary
rdastat_envi <- purrr::reduce(list(rdafjordlandscape, rdafjordgrazing, rdafjordlocenvi, rdalandscapegrazing, rdalandscapelocenvi, rdagrazinglocenvi), dplyr::full_join)



#### RDA grass community data ####

#
## Fjord x grass

# Remove space within plant names (otherwise parse error in plot)
names(grass) <- gsub(" ", ".", names(grass))

# Redundancy analysis
rda <- rda(subset(grass, select = -c(SiteID)) ~ ., data = subset(fjordsys, select = -c(SiteID)))
plotrda_fjordxgrass <- ggrda(rda,
                             txt = 4,
                             ptslab = TRUE,
                             addcol = "black",
                             addsize = 3,
                             size = 1,
                             arrow = 0.3,
                             #repel = TRUE,
                             veccol = "cyan3",
                             labcol = "cyan3",
                               xlims = c(-1, 1.3),
                               ylims = c(-1.1, 1.2))
plotrda_fjordxgrass

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafjordgrass <- rdatable |>
  mutate(model = "FjordxGrass", explanatory = "Fjord", response = "Grass") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Landscape x grass

# Remove space within plant names (otherwise parse error in plot)
names(grass) <- gsub(" ", ".", names(grass))

# Redundancy analysis
rda <- rda(subset(grass, select = -c(SiteID)) ~ ., data = subset(landscape, select = -c(SiteID)))
plotrda_landscapexgrass <- ggrda(rda,
                                 txt = 4,
                                 ptslab = TRUE,
                                 addcol = "black",
                                 addsize = 3,
                                 size = 1,
                                 arrow = 0.3,
                                 #repel = TRUE,
                                 veccol = "chartreuse4",
                                 labcol = "chartreuse4",
                             xlims = c(-1, 1),
                             ylims = c(-1, 1))
plotrda_landscapexgrass

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapegrass <- rdatable |>
  mutate(model = "LandscapexGrass", explanatory = "Landscape", response = "Grass") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Grazing management x grass

# Remove space within plant names (otherwise parse error in plot)
names(grass) <- gsub(" ", ".", names(grass))

# Redundancy analysis
rda <- rda(subset(grass, select = -c(SiteID)) ~ ., data = subset(grazing, select = -c(SiteID)))
plotrda_grazingxgrass <- ggrda(rda,
                               txt = 4,
                               ptslab = TRUE,
                               addcol = "black",
                               addsize = 3,
                               size = 1,
                               arrow = 0.3,
                               repel = TRUE,
                               veccol = "darkgoldenrod3",
                               labcol = "darkgoldenrod3",
                               xlims = c(-1.1, 1),
                               ylims = c(-1.2, 0.9))
plotrda_grazingxgrass

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdagrazinggrass <- rdatable |>
  mutate(model = "GrazingxGrass", explanatory = "Grazing", response = "Grass") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Fine-scale environment x grass

# Remove space within plant names (otherwise parse error in plot)
names(grass) <- gsub(" ", ".", names(grass))

# Redundancy analysis
rda <- rda(subset(grass, select = -c(SiteID)) ~ ., data = subset(locenvi, select = -c(SiteID)))
plotrda_locenvixgrass <- ggrda(rda,
                               txt = 4,
                               ptslab = TRUE,
                               addcol = "black",
                               addsize = 3,
                               size = 1,
                               arrow = 0.3,
                               repel = TRUE,
                               veccol = "darkred",
                               labcol = "darkred",
                               xlims = c(-1.2, 0.9),
                               ylims = c(-1.1, 1))
plotrda_locenvixgrass

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "term")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "term") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalocenvigrass <- rdatable |>
  mutate(model = "FinescalexGrass", explanatory = "Finescale", response = "Grass") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Grass model summary

# Plot Grouping
rda_grass <- ggarrange(plotrda_fjordxgrass, plotrda_landscapexgrass, plotrda_grazingxgrass, plotrda_locenvixgrass, 
                      labels = c("A", "B", "C", "D"),
                      ncol = 2,
                      nrow = 2)
rda_grass
ggsave("outputs/RDA_grass.png", plot = rda_grass, width = 15, height = 13, units = "cm", bg = "white")

# Stat summary
rdastat_grass <- purrr::reduce(list(rdafjordgrass, rdalandscapegrass, rdagrazinggrass, rdalocenvigrass), dplyr::full_join)



#### RDA forb community data ####

#
## Fjord x forb (CANCOR NS)

# Remove space within plant names (otherwise parse error in plot)
names(forb) <- gsub(" ", ".", names(forb))

# Redundancy analysis
rda <- rda(subset(forb, select = -c(SiteID)) ~ ., data = subset(fjordsys, select = -c(SiteID)))
plotrda_fjordxforb <- ggrda(rda,
                            txt = 4,
                            ptslab = TRUE,
                            addcol = "black",
                            addsize = 3,
                            size = 1,
                            arrow = 0.3,
                            #repel = TRUE,
                            veccol = "cyan3",
                            labcol = "cyan3",
                             xlims = c(-1.1, 1),
                             ylims = c(-1.2, 0.9))
plotrda_fjordxforb

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |>
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |>
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafjordforb <- rdatable |>
  mutate(model = "FjordxForb", explanatory = "Fjord", response = "Forb") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Landscape x forb (NS)

# Remove space within plant names (otherwise parse error in plot)
names(forb) <- gsub(" ", ".", names(forb))

# Redundancy analysis
rda <- rda(subset(forb, select = -c(SiteID)) ~ ., data = subset(landscape, select = -c(SiteID)))
plotrda_landscapexforb <- ggrda(rda,
                                txt = 4,
                                ptslab = TRUE,
                                addcol = "black",
                                addsize = 3,
                                size = 1,
                                arrow = 0.3,
                                #repel = TRUE,
                                veccol = "chartreuse4",
                                labcol = "chartreuse4",
                                 xlims = c(-0.9, 1.2),
                                 ylims = c(-0.9, 1.2))
plotrda_landscapexforb

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |>
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapeforb <- rdatable |>
  mutate(model = "LandscapexForb", explanatory = "Landscape", response = "Forb") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Grazing management x forb

# Remove space within plant names (otherwise parse error in plot)
names(forb) <- gsub(" ", ".", names(forb))

# Redundancy analysis
rda <- rda(subset(forb, select = -c(SiteID)) ~ ., data = subset(grazing, select = -c(SiteID)))
plotrda_grazingxforb <- ggrda(rda,
                              txt = 4,
                              ptslab = TRUE,
                              addcol = "black",
                              addsize = 3,
                              size = 1,
                              arrow = 0.3,
                              #repel = TRUE,
                              veccol = "darkgoldenrod3",
                              labcol = "darkgoldenrod3",
                               xlims = c(-1.2, 0.8),
                               ylims = c(-1.1, 0.9))
plotrda_grazingxforb

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |>
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdagrazingforb <- rdatable |>
  mutate(model = "GrazingxForb", explanatory = "Grazing", response = "Forb") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Fine-scale environment x forb

# Remove space within plant names (otherwise parse error in plot)
names(forb) <- gsub(" ", ".", names(forb))

# Redundancy analysis
rda <- rda(subset(forb, select = -c(SiteID)) ~ ., data = subset(locenvi, select = -c(SiteID)))
plotrda_locenvixforb <- ggrda(rda,
                              txt = 4,
                              ptslab = TRUE,
                              addcol = "black",
                              addsize = 3,
                              size = 1,
                              arrow = 0.3,
                              #repel = TRUE,
                              veccol = "darkred",
                              labcol = "darkred",
                              xlims = c(-1, 1.1),
                              ylims = c(-0.9, 1.2))
plotrda_locenvixforb

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalocenviforb <- rdatable |>
  mutate(model = "FinescalexForb", explanatory = "Finescale", response = "Forb") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Forb model summary

# Plot Grouping
rda_forb <- ggarrange(plotrda_fjordxforb, plotrda_landscapexforb, plotrda_grazingxforb, plotrda_locenvixforb, 
                       labels = c("A", "B", "C", "D"),
                       ncol = 2,
                       nrow = 2)
rda_forb
ggsave("outputs/RDA_forb.png", plot = rda_forb, width = 15, height = 13, units = "cm", bg = "white")

# Stat summary
rdastat_forb <- purrr::reduce(list(rdafjordforb, rdalandscapeforb, rdagrazingforb, rdalocenviforb), dplyr::full_join)



#### RDA beetle community data ####

#
## Fjord x forb (CANCOR NS)

# Redundancy analysis
rda <- rda(subset(beetle, select = -c(SiteID)) ~ ., data = subset(fjordsys, select = -c(SiteID)))
plotrda_fjordxbeetle <- ggrda(rda,
                              txt = 4,
                              ptslab = TRUE,
                              addcol = "black",
                              addsize = 3,
                              size = 1,
                              arrow = 0.3,
                              #repel = TRUE,
                              veccol = "cyan3",
                              labcol = "cyan3",
                            xlims = c(-1.65, 1.65),
                            ylims = c(-1.5, 1.8))
plotrda_fjordxbeetle

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafjordbeetle <- rdatable |>
  mutate(model = "FjordxBeetle", explanatory = "Fjord", response = "Beetle") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Landscape x beetle (NS)

# Redundancy analysis
rda <- rda(subset(beetle, select = -c(SiteID)) ~ ., data = subset(landscape, select = -c(SiteID)))
plotrda_landscapexbeetle <- ggrda(rda,
                                  txt = 4,
                                  ptslab = TRUE,
                                  addcol = "black",
                                  addsize = 3,
                                  size = 1,
                                  arrow = 0.3,
                                  #repel = TRUE,
                                  veccol = "chartreuse4",
                                  labcol = "chartreuse4",
                                xlims = c(-1.7, 1.7),
                                ylims = c(-2.3, 1.1))
plotrda_landscapexbeetle

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "term")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "term") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapebeetle <- rdatable |>
  mutate(model = "LandscapexBeetle", explanatory = "Landscape", response = "Beetle") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Grazing management x beetle (NS)

# Redundancy analysis
rda <- rda(subset(beetle, select = -c(SiteID)) ~ ., data = subset(grazing, select = -c(SiteID)))
plotrda_grazingxbeetle <- ggrda(rda,
                                txt = 4,
                                ptslab = TRUE,
                                addcol = "black",
                                addsize = 3,
                                size = 1,
                                arrow = 0.3,
                                #repel = TRUE,
                                veccol = "darkgoldenrod3",
                                labcol = "darkgoldenrod3",
                              xlims = c(-1, 1),
                              ylims = c(-1, 1))
plotrda_grazingxbeetle

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdagrazingbeetle <- rdatable |>
  mutate(model = "GrazingxBeetle", explanatory = "Grazing", response = "Beetle") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Fine-scale environment x beetle

# Redundancy analysis
rda <- rda(subset(beetle, select = -c(SiteID)) ~ ., data = subset(locenvi, select = -c(SiteID)))
plotrda_locenvixbeetle <- ggrda(rda,
                                txt = 4,
                                ptslab = TRUE,
                                addcol = "black",
                                addsize = 3,
                                size = 1,
                                arrow = 0.3,
                                repel = TRUE,
                                veccol = "darkred",
                                labcol = "darkred",
                              xlims = c(-1.1, 0.9),
                              ylims = c(-0.9, 1.1))
plotrda_locenvixbeetle

# Summary rda
#summary(rda)

# Total variance explained by RDA
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

# Global RDA significance by permutation
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

# Individual axis significance
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Individual term significance
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

# Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

# Summary table
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalocenvibeetle <- rdatable |>
  mutate(model = "FinescalexBeetle", explanatory = "Finescale", response = "Beetle") |> 
  filter(type == "model" | PrF < 0.1 | dim == "RDA1")

#
## Beetle model summary

# Plot Grouping
rda_beetle <- ggarrange(plotrda_fjordxbeetle, plotrda_landscapexbeetle, plotrda_grazingxbeetle, plotrda_locenvixbeetle, 
                        labels = c("A", "B", "C", "D"),
                        ncol = 2,
                        nrow = 2)
rda_beetle
ggsave("outputs/RDA_beetle.png", plot = rda_beetle, width = 15, height = 13, units = "cm", bg = "white")

# Stat summary
rdastat_beetle <- purrr::reduce(list(rdafjordbeetle, rdalandscapebeetle, rdagrazingbeetle, rdalocenvibeetle), dplyr::full_join)



#### Significant RDA results ####

#
## RDA summary for paper

# Plot
rdall <- ggarrange(plotrda_fjordxlandscape, plotrda_grazingxlocenvi, plotrda_locenvixgrass, plotrda_locenvixforb, 
                   labels = c("A", "B", "C", "D"),
                   ncol = 2,
                   nrow = 2)
rdall
ggsave("outputs/RDAresults.png", plot = rdall, width = 15, height = 15, units = "cm", bg = "white")

# RDA stat
rdastat_all <- purrr::reduce(list(rdafjordlandscape, rdagrazinglocenvi, rdafjordgrass, rdafjordbeetle, rdagrazinggrass, rdagrazingbeetle, rdalocenvigrass, rdalocenviforb), dplyr::full_join)

#
## Distribution residuals for significant relationships with terms

# Ordination community data
DCA_grass <- decorana(contin_grass)
DCA_grass # Axis length of 1.9 -> run PCA
PCA_grass <- prcomp(contin_grass)
PCA_grass
biplot(PCA_grass)
DCA_forb <- decorana(contin_forb)
DCA_forb # Axis length of 3.1 -> keep DCA
plot(DCA_forb)

# Extraction grass species scores from PCA
scores_grass <- as.data.frame(scores(PCA_grass, choices = c(1,2), display = "sites"))
names(scores_grass) <- gsub("PC1", "GrassPCA1", names(scores_grass))
names(scores_grass) <- gsub("PC2", "GrassPCA2", names(scores_grass))
scores_grass <- scores_grass |>
  mutate(SiteID = row.names(scores_grass)) # PlotID as column for future binding

# Extraction forb species scores from DCA
scores_forb <- as.data.frame(scores(DCA_forb, choices = c(1,2), display = "sites"))
names(scores_forb) <- gsub("DCA1", "ForbDCA1", names(scores_forb))
names(scores_forb) <- gsub("DCA2", "ForbDCA2", names(scores_forb))
scores_forb <- scores_forb |>
  mutate(SiteID = row.names(scores_forb)) # PlotID as column for future binding

# Table binding
allvar <- purrr::reduce(list(scores_grass, scores_forb, fjordsys, landscape, grazing, locenvi), dplyr::left_join)

# Distribution residuals fjord x landscape
# residualPlots(lm(Forest~JulTemp, data = allvar)) # random -> validated
# residualPlots(lm(Forest~SeaDist, data = allvar)) # random -> validated
# residualPlots(lm(Wetland~JulTemp, data = allvar)) # random -> validated
# residualPlots(lm(Wetland~SeaDist, data = allvar)) # random -> validated
# residualPlots(lm(Outfield~JulTemp, data = allvar)) # random -> validated
# residualPlots(lm(Outfield~SeaDist, data = allvar)) # random -> validated

# Distribution residuals finescale x grass
# residualPlots(lm(GrassPCA1~pH, data = allvar)) # random -> validated
# residualPlots(lm(GrassPCA1~Bryo, data = allvar)) # random -> validated
# residualPlot(lm(GrassPCA1~MeanHeight, data = allvar)) # random -> validated
# residualPlots(lm(GrassPCA1~SoilPene, data = allvar)) # random -> validated

# Distribution residuals finescale x forb
# residualPlots(lm(ForbDCA1~SoilPene, data = allvar)) # random -> validated
# residualPlot(lm(ForbDCA1~MeanHeight, data = allvar)) # random -> validated