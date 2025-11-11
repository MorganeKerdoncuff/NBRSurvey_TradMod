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
library(ggvegan) # Function autoplot.rda for ordination plot with ggplot2
library(ggpubr) # Function ggarrange for several plots on same file
library(GGally) # Extension ggplot
library(ggimage) # Include image in figures
library(rsvg) # SVG format for figures
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

# Filter grassland sites

## Site selection & removal geographical outlier IG3
siteinfo_infield <- siteinfo_full |>  
  filter(Habitat == "permanent grassland" & SiteID != "IG3") |> 
  dplyr::select(SiteID, Livestock)

## Extraction in other datasets
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

# Summarise data at site level -> should be 29 observations for non community data

## Site info - validated
## Climate - validated
## Landscape matrix - validated
## Field management - validated
## Area 20 x 20 - validated

## Ground cover - current at sample level -> summary by average
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

## Bulk density - current at sample level -> summary by average
soilbulk_infield <- soilbulk_infield |> 
  group_by(SiteID) |> 
  summarise(BD = mean(BD),
            GWC = mean(GWC_48))

## Soil chemistry - current at plot level -> summary by average
soilchem_infield <- soilchem_infield |> 
  # Post-analysis verification to determine if soil penetration is related to soil texture -> answer is no
  # mutate(SoilType = ifelse(
  #   SoilType == "Medium_sand", 2, ifelse(
  #     SoilType == "Fine_sand", 3, ifelse(
  #       SoilType == "Silty_medium_sand", 5, ifelse(
  #         SoilType == "Silty_fine_sand", 6, ifelse(
  #           SoilType == "Mineral_mixed_humus_soil", 13, 14
  #         )))))) |>
  group_by(SiteID) |> 
  summarise(LOI = mean(LOI),
            #SoilType = mean(SoilType),
            Humus = mean(Humus_percentDM),
            pH = mean(pH),
            Phosph = mean(P.Al_mg.100g),
            Potass = mean(K.Al_mg.100g),
            Nitrog = mean(TotalN_percentDM),
            Sodium = mean(Na.Al_mg.100g))

## Soil penetration - current at sample level -> summary by average
soilpene_infield <- soilpene_infield |> 
  group_by(SiteID) |> 
  summarise(SoilPene = mean(AveragePT_cm))

## Plant community - current at quadrat level -> summary by average
vege_infield <- vege_infield |> 
  group_by(SiteID, Species) |> 
  summarise(PlantSp_cover = mean(Abundance))

## Beetle community - current at pitfall level -> summary by average
beetle_infield <- beetle_infield |> 
  group_by(SiteID, BeetleFamilies) |> 
  summarise(BeetleFam_abundance = mean(BeetleFam_abundance))

# Selection & transformation plant community data

## Data distribution
#hist(vege_grass$PlantSp_cover) # Poisson, highly skewed

## Species frequencies across sites
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

## Species average percent cover across site
vege_average <- vege_infield |>
  group_by(Species) |>
  summarise_if(is.numeric, mean, na.rm = TRUE) |>
  dplyr::arrange(desc(PlantSp_cover)) |> 
  # filtering by average value across sites
  mutate(FunctionalGroup = case_when(
    grepl("Achillea|Alchemilla|Anagalis|Anemone|Angelica|Anthriscus|Armeria|Bartsia|Campanula|Cardamine|Cerastium|Cirsium|Conopodium|Dactylorhiza|Digitalis|Epilobium|Euphrasia|Fraxinus|Galeopsis|Galium|Geranium|Gnaphalium|Hieracium|Hypericum|Hypochaeris|Lathyrus|Leontodon|Lotus|Melampyrum|Moneses|Myosotis|Narthecium|Oxalis|Pedicularis|Pinguicula|Plantago|Potentilla|Prunella|Ranunculus|Rumex|Sagina|Sedum|Senecio|Silene|Solidago|Stellaria|Succisa|Taraxacum|Trientalis|Trifolium|Urtica|Valeriana|Veronica|Vicia|Viola", Species) ~ "forbs",
    grepl("Agrostis|Aira|Alopecurus|Anthoxanthum|Bromopsis|Calamagrostis|Dactylis|Danthonia|Deschampsia|Festuca|Holcus|Lolium|Molinia|Nardus|Phleum|Poa", Species) ~ "grasses",
    grepl("Carex|Eriophorum|Juncus|Luzula|Trichophorum", Species) ~ "monocotyledons",
    grepl("Andromeda|Arctostaphylos|Betula|Calluna|Chamaepericlymenum|Empetrum|Erica|Juniperus|Loiseleuria|Picea|Polygala|Polygonum|Populus|Prunus|Rubus|Salix|Sorbus|Ulmus|Vaccinium", Species) ~ "woody",
    grepl("Athyrium|Blechnum|Dryopteris|Gymnocarpium|Phegopteris|Polypodium|Pteridium", Species) ~ "ferns",
    .default = "cryptogams"
  ))

filter(vege_average, PlantSp_cover >= 3)
# 7 grass species over 3 % A. capillaris; F. rubra; H. lanatus; P. pratensis; D. cespitosa; A. odoratum; L. perenne
# 3 forb species over 3 % T. repens; R. acetosa; G. saxatile
filter(vege_average, PlantSp_cover < 3 & PlantSp_cover > 1)
# 2 grass species between 1 & 3% D. flexuosa; P. trivialis
# 4 forb species over 1%  P. erecta; R. repens; A. millefolium; R acris
# filter(vege_freq, Species == "Lolium perenne") only present in 6 sites

# Selection and transformation grass assemblage

## Selection main grass species - at least present in 10 sites AND min average cover 1% -> 8 species
grass <- subset(vege_infield,
                Species == "Agrostis capillaris" |
                Species == "Festuca rubra" | 
                Species == "Holcus lanatus" | 
                Species == "Poa pratensis" | 
                Species == "Deschampsia cespitosa" | 
                Species == "Anthoxantum odoratum" | 
                Species == "Deschampsia flexuosa" | 
                Species == "Poa trivialis")

## Hellinger transformation on contingency table (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
contin_grass <- xtabs(formula = PlantSp_cover ~ SiteID + Species, data = grass)
contin_grass <- decostand(contin_grass, method = "hellinger")

## Wide table
grass <- as.data.frame(contin_grass)
grass <- grass |>
  pivot_wider(names_from = Species, values_from = Freq)

## Suitable variable names
names(grass) <- gsub("Agrostis capillaris", "A.capillaris", names(grass))
names(grass) <- gsub("Festuca rubra", "F.rubra", names(grass))
names(grass) <- gsub("Holcus lanatus", "H.lanatus", names(grass))
names(grass) <- gsub("Poa pratensis", "P.pratensis", names(grass))
names(grass) <- gsub("Deschampsia cespitosa", "D.cespitosa", names(grass))
names(grass) <- gsub("Anthoxantum odoratum", "A.odoratum", names(grass))
names(grass) <- gsub("Deschampsia flexuosa", "D.flexuosa", names(grass))
names(grass) <- gsub("Poa trivialis", "P.trivialis", names(grass))

# Selection and transformation forb assemblage

## Selection main forb species - at least present in 10 sites AND min average cover 1% -> 7 species
forb <- subset(vege_infield,
                Species == "Trifolium repens" |
                  Species == "Rumex acetosa" |
                  Species == "Galium saxatile" |
                  Species == "Potentilla erecta" |
                  Species == "Ranunculus acris" |
                  Species == "Ranunculus repens" |
                  Species == "Achillea millefolium")

## Hellinger transformation on contingency table (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
contin_forb <- xtabs(formula = PlantSp_cover ~ SiteID + Species, data = forb)
contin_forb <- decostand(contin_forb, method = "hellinger")

## Wide table
forb <- as.data.frame(contin_forb)
forb <- forb |>
  pivot_wider(names_from = Species, values_from = Freq)

## Suitable variable names
names(forb) <- gsub("Trifolium repens", "T.repens", names(forb))
names(forb) <- gsub("Rumex acetosa", "R.acetosa", names(forb))
names(forb) <- gsub("Galium saxatile", "G.saxatile", names(forb))
names(forb) <- gsub("Potentilla erecta", "P.erecta", names(forb))
names(forb) <- gsub("Ranunculus acris", "R.acris", names(forb))
names(forb) <- gsub("Ranunculus repens", "R.repens", names(forb))
names(forb) <- gsub("Achillea millefolium", "A.millefolium", names(forb))

## Transformation beetle assemblage data

#hist(arthro_grass$BeetleFam_abundance) # Poisson, highly skewed

## Family frequencies across sites
beetle_freq <- filter(beetle_infield, BeetleFam_abundance>0) |>
  group_by(BeetleFamilies) |>
  count() |>
  dplyr::arrange(desc(n)) # filtering by frequency
filter(beetle_freq, n > 10)
# 8 beetle families in at least 10 sites Carab; Hydro; Scara; Staph; Ptili; Curcu; Elat; Silph

## Species average percent cover across site
beetle_average <- beetle_infield |>
  group_by(BeetleFamilies) |>
  summarise_if(is.numeric, mean, na.rm = TRUE) |>
  dplyr::arrange(desc(BeetleFam_abundance)) # filtering by average value across sites
filter(beetle_average, BeetleFam_abundance > 3)
# 5 beetles families with at least 3 individuals on average Staph; Ptili; Hydro; Scara; Carab

## Selection main dung beetle families -> at least present in 10 sites + min 3 individuals on average -> 6 families
beetle <- subset(beetle_infield,
                 BeetleFamilies == "Carabidae" |
                   BeetleFamilies == "Staphylinidae" |
                   BeetleFamilies == "Hydrophilidae" |
                   BeetleFamilies == "Ptiliidae" |
                   BeetleFamilies == "Scarabaeidae")

## Hellinger transformation on contingency table (Borcard, Gillet and Legendre 2011; Legendre and Gallagher 2001)
contin_beetle <- xtabs(formula = BeetleFam_abundance ~ SiteID + BeetleFamilies, data = beetle)
contin_beetle <- decostand(contin_beetle, method = "hellinger")

## Wide table
beetle <- as.data.frame(contin_beetle)
beetle <- beetle |>
  pivot_wider(names_from = BeetleFamilies, values_from = Freq)

# Explanatory set regional scale

## Desired variables
### Elevation (num), from area 20x20
### Annual precipitation (num), from climate
### Average July temp (num), from climate
### Average Jan temp (num), from climate
### Distance to marine water body (num), from area 20x20

## Selection variables
regional <- full_join(area20x20_infield, climate_infield)
regional <- regional |>
  mutate(temprange = maxtempJuly - mintempJan)
regional <- subset(regional, select = c(SiteID, Elevation, annualprecipitation, avgtempJuly, avgtempJan, temprange, DistanceToSea_m))
# regional <- subset(regional, select = c(SiteID, Elevation, annualprecipitation, maxtempJuly, mintempJan, DistanceToSea_m))

## Suitable variable names
names(regional) <- gsub("Elevation", "Elev", names(regional))
names(regional) <- gsub("annualprecipitation", "AnnPreci", names(regional))
# names(regional) <- gsub("maxtempJuly", "JulTemp", names(regional))
# names(regional) <- gsub("mintempJan", "JanTemp", names(regional))
names(regional) <- gsub("avgtempJuly", "JulTemp", names(regional))
names(regional) <- gsub("avgtempJan", "JanTemp", names(regional))
names(regional) <- gsub("temprange", "TempRange", names(regional))
names(regional) <- gsub("DistanceToSea_m", "FjordDist", names(regional))

## Variable distribution
#hist(regional$Elev) # Poisson distribution, not skewed -> validated
#hist(regional$AnnPreci) # Normal-like distribution, with some gaps -> validated
#hist(regional$JulTemp) # Normal-like distribution -> validated
#hist(regional$JanTemp) # Reverse Poisson distribution -> validated
#hist(regional$TempRange) # Normal-like distribution -> validated
#hist(regional$SeaDist) # Poisson-like distribution, with one gap -> validated

## Scaling numerical variables
regional_sc <- regional |> 
  mutate(across(where(is.numeric), scale))

# Explanatory set landscape scale

## Desired variables
### Sum forest area (num), from landscape
### Sum cultivated area (num), from landscape
### Sum infield area (num), from landscape
### Sum outfield area (num), from landscape
### Sum wetland area (num), from landscape

## Percent of land cover
landscape_infield <- landscape_infield |> 
  mutate(TotLand_ha = TotalArea_ha - Sea_ha) |> 
  mutate(Forest = (ProductiveForest_ha + NonProductiveForest_ha)/TotLand_ha*100) |> 
  mutate(Cultivated = (FullyCultivatedLand_ha + SuperficiallyCultivatedLand_ha)/TotLand_ha*100) |> 
  mutate(Infield = Infield_ha/TotLand_ha*100) |> 
  mutate(Outfield = Outfield_ha/TotLand_ha*100) |> 
  mutate(Wetland = Wetland_ha/TotLand_ha*100) |> 
  mutate(Artificial = Infrastructure_ha/TotLand_ha*100)

## Selection variables
landscape <- subset(landscape_infield, select = c(SiteID, Cultivated, Forest, Infield, Outfield, Wetland, Artificial))

## Variable distribution
# hist(landscape$Cultivated) # Normal-like distribution, one gap in the middle -> validated
# hist(landscape$Forest) # uneven distribution but no outlier -> validated
# hist(landscape$Infield) # unbalanced normal-like, but no outlier -> validated
# hist(landscape$Outfield) # Poisson, not skewed -> validated
# hist(landscape$Wetland) # Poisson, not skewed -> validated
# hist(landscape$Artificial) # Poisson, no zero -> valiated

## Scaling numerical variables
landscape_sc <- landscape |> 
  mutate(across(where(is.numeric), scale))

# Explanatory set field scale

## Desired variables
### Type of livestock (char), from landuse
### Number of adult animals (num), from landuse
### Grazing surface of the collected field (num), from landuse
### Average stocking density, from landuse

## Dummy numeric for livestock variable
landuse_infield <- landuse_infield |> 
  mutate(Sheep = ifelse(Livestock1 == "sheep",1,0)) |> 
  mutate(Cattle = ifelse(Livestock1 == "cattle", 1,0))

## Selection variables
field <- subset(landuse_infield, select = c(SiteID, Sheep, Cattle, FlockSize1_adults, SelectedFieldArea_ha, AvgStockDensity_perha))

## Suitable variable names
names(field) <- gsub("FlockSize1_adults", "NbAdults", names(field))
names(field) <- gsub("SelectedFieldArea_ha", "FieldArea", names(field))
# names(field) <- gsub("FarmInfieldArea_ha", "FarmInfield", names(field))
names(field) <- gsub("AvgStockDensity_perha", "StockDens", names(field))

## Variable distribution
#hist(field$NbAdults) # one outlier above 150 -> rejected
#hist(field$FieldArea) # Poisson distribution, with few gaps -> validated
#hist(field$FarmInfield) # Poisson-like distribution, a bit unbalanced but no outlier
#hist(field$StockDens) # Poisson distribution -> validated

## Scaling numeric variables
field <- subset(field, select = -c(NbAdults))
field_sc <- field |> 
  mutate(across(where(is.numeric), scale))

# Explanatory set fine scale

## Desired variables
### litter and bryophyte percent cover from groundcover
### sward average height from groundcover
### bulk density and gravimetric water content from soilbulk
### LOI, nitrogen, phosphorus, pH & humus content from soilchem
### soil resistance to penetration from soilpene
### aspect & slope from area20x20

## Selection variables
fine <- purrr::reduce(list(groundcover_infield, soilbulk_infield, soilpene_infield, soilchem_infield, area20x20_infield), dplyr::left_join)
fine <- subset(fine, select = c(SiteID, Litter, Bryo, MeanHeight, BD, GWC, LOI, Nitrog, Phosph, pH, Humus, SoilPene, Aspect_degree, Slope_degree))

## Suitable variable names
names(fine) <- gsub("Aspect_degree", "Aspect", names(fine))
names(fine) <- gsub("Slope_degree", "Slope", names(fine))

# Variable distribution
#hist(fine$Litter) # Highly skewed Poisson, very small range with 2 outliers -> rejected
#hist(fine$Bryo) # Unbalanced Normal-like but no outlier -> validated
#hist(fine$MeanHeight) # Poisson like, no outlier -> validated
#hist(fine$BD) # Normal-like, a bit unbalanced but no outlier -> validated
#hist(fine$GWC) # poisson like, one gap between 50 & 55 -> validated
#hist(fine$LOI) # Poisson like, a bit skewed but no outlier -> validated
#hist(fine$Nitrog) # Poisson like, with a gap but no outlier -> validated
#hist(fine$Phosph) # poisson like, a bit skewed but no outlier -> validated
#hist(fine$pH) # good Normal distribution -> validated
#hist(fine$Humus) # Poisson distribution, a bit skewed -> validated
#hist(fine$SoilPene) # Normal distribution a bit unbalanced, no outlier -> validated
#hist(fine$Aspect) # inverse Normal distribution, but no outlier -> validated
#hist(fine$Slope) # Poisson distribution, one gap but no outlier -> validated

## Scaling numeric variables
fine <- subset(fine, select = -c(Litter))
fine_sc <- fine |> 
  mutate(across(where(is.numeric), scale))

# Variable filtering - colinearity

## Function for paired panels
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y, method = "spearman")) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

## Regional variables
pairs(select_if(regional_sc, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # high colinearity between Jan temp and temp range

# Landscape variables
pairs(select_if(landscape_sc, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # one correlation a bit strong (0.74) between total cultivated land and infield percent

# Field variables
pairs(select_if(field_sc, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # no correlation stronger than 0.64 -> validated

# Fine variables
pairs(select_if(fine_sc, is.numeric),
      upper.panel = panel.cor,
      lower.panel = panel.smooth) # 5 variables highly correlated (>0.8): BD, moisture, LOI, humus, Nitrogen

# Removal strongly correlated variables - BD kept as with higher number of replication for each site
regional_sc <- subset(regional_sc, select = -c(JanTemp))
fine_sc <- subset(fine_sc, select = -c(GWC, LOI, Humus, Nitrog))

# Contingency tables for correlation analyses

# Regional set
regional_long <- regional_sc |> 
  pivot_longer(
    cols = c(Elev, AnnPreci, JulTemp, TempRange, FjordDist),
    names_to = "Factors",
    values_to = "Values")
contin_regional <- xtabs(formula = Values ~ SiteID + Factors, data = regional_long)

# Landscape set
landscape_long <- landscape_sc |> 
  pivot_longer(
    cols = c(Cultivated, Forest, Infield, Outfield, Wetland, Artificial),
    names_to = "Factors",
    values_to = "Values"
  )
contin_landscape <- xtabs(formula = Values ~ SiteID + Factors, data = landscape_long)

# Field set
field_long <- field_sc |> 
  pivot_longer(
    cols = c(Sheep, Cattle, FieldArea, StockDens),
    names_to = "Factors",
    values_to = "Values")
contin_field <- xtabs(formula = Values ~ SiteID + Factors, data = field_long)

# Fine set
fine_long <- fine_sc |> 
  pivot_longer(
    cols = c(Bryo, MeanHeight, BD, Phosph, pH, SoilPene, Aspect, Slope),
    names_to = "Factors",
    values_to = "Values")
contin_fine <- xtabs(formula = Values ~ SiteID + Factors, data = fine_long)


#### Canonical correlation analyses between explanatory sets ####

## Number of observations (same for all sets)
nobs <- dim(contin_regional)[1]

## Number of variables per set
nvar_regional <- length(select_if(regional_sc, is.numeric))
nvar_landscape <- length(select_if(landscape_sc, is.numeric))
nvar_field <- length(select_if(field_sc, is.numeric))
nvar_fine <- length(select_if(fine_sc, is.numeric))

# Correlations with the regional explanatory set

## Regional x landscape
cancor <- cc(contin_regional, contin_landscape) # Canonical correlation
rhocancor <- cancor$cor
rhocancor # Canonical correlation coefficient
test <- p.asym(rhocancor, nobs, nvar_regional, nvar_landscape, tstat = "Hotelling") # Hotelling test
ccaregiolandscape <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Regional", response = "Landscape", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccaregiolandscape # Summary statistics

## Regional x field
cancor <- cc(contin_regional, contin_field)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_regional, nvar_field, tstat = "Hotelling")
ccaregiofield <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Regional", response = "Field", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccaregiofield

## Regional x fine
cancor <- cc(contin_regional, contin_fine)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_regional, nvar_fine, tstat = "Hotelling")
ccaregiofine <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Regional", response = "Fine", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccaregiofine

# Correlations with the landscape explanatory set

## Landscape x field
cancor <- cc(contin_landscape, contin_field)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_field, tstat = "Hotelling")
ccalandscapefield <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Field", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccalandscapefield

## Landscape x fine
cancor <- cc(contin_landscape, contin_fine)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_fine, tstat = "Hotelling")
ccalandscapefine <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Fine", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccalandscapefine

# Correlations with the field explanatory set

## Field x fine
cancor <- cc(contin_field, contin_fine)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_field, nvar_fine, tstat = "Hotelling")
ccafieldfine <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Field", response = "Fine", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccafieldfine

# Summary statistics

## All dimensions included
envar_SI <- purrr::reduce(list(ccaregiolandscape, ccaregiofield, ccaregiofine, ccalandscapefield, ccalandscapefine, ccafieldfine), dplyr::full_join)

## First dimension + other dimensions if marginally significant only
envar <- envar_SI |> 
  filter(dim == "dim1" | p.value < 0.1)


#### Canonical correlation analyses between explanatory sets and dominant grass assemblage ####

## Number of observations (same for all sets)
nobs <- dim(contin_regional)[1]

## Number of variables in each set
nvar_regional <- length(select_if(regional_sc, is.numeric))
nvar_landscape <- length(select_if(landscape_sc, is.numeric))
nvar_field <- length(select_if(field_sc, is.numeric))
nvar_fine <- length(select_if(fine_sc, is.numeric))
nvar_grass <- length(select_if(grass, is.numeric))

## Regional x grass
cancor <- cc(contin_regional, contin_grass)
rhocancor <- cancor$cor
rhocancor 
test <- p.asym(rhocancor, nobs, nvar_regional, nvar_grass, tstat = "Hotelling") 
ccaregiograss <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Regional", response = "Grass", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccaregiograss

## Landscape x grass
cancor <- cc(contin_landscape, contin_grass)
rhocancor <- cancor$cor
rhocancor 
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_grass, tstat = "Hotelling") 
ccalandscapegrass <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Grass", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccalandscapegrass

## Field x grass
cancor <- cc(contin_field, contin_grass)
rhocancor <- cancor$cor
rhocancor 
test <- p.asym(rhocancor, nobs, nvar_field, nvar_grass, tstat = "Hotelling") 
ccafieldgrass <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Field", response = "Grass", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccafieldgrass

## Fine x grass
cancor <- cc(contin_fine, contin_grass)
rhocancor <- cancor$cor
rhocancor 
test <- p.asym(rhocancor, nobs, nvar_fine, nvar_grass, tstat = "Hotelling")
ccafinegrass <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fine", response = "Grass", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccafinegrass

# Summary statistics

## All dimensions included
grassvar_SI <- purrr::reduce(list(ccaregiograss, ccalandscapegrass, ccafieldgrass, ccafinegrass), dplyr::full_join)

## First dimension + other dimensions if marginally significant only
grassvar <- grassvar_SI |> 
  filter(dim == "dim1" | p.value < 0.1)


#### Canonical correlation analyses between explanatory sets and dominant forb assemblage ####

# Number of observations (same for all sets)
nobs <- dim(contin_regional)[1]

# Number of variables in each set
nvar_regional <- length(select_if(regional_sc, is.numeric))
nvar_landscape <- length(select_if(landscape_sc, is.numeric))
nvar_field <- length(select_if(field_sc, is.numeric))
nvar_fine <- length(select_if(fine_sc, is.numeric))
nvar_forb <- length(select_if(forb, is.numeric))

## Regional x forb
cancor <- cc(contin_regional, contin_forb)
rhocancor <- cancor$cor
rhocancor 
test <- p.asym(rhocancor, nobs, nvar_regional, nvar_forb, tstat = "Hotelling") 
ccaregioforb <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Regional", response = "Forb", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccaregioforb

## Landscape x forb
cancor <- cc(contin_landscape, contin_forb)
rhocancor <- cancor$cor
rhocancor 
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_forb, tstat = "Hotelling") 
ccalandscapeforb <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Forb", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccalandscapeforb

## Field x forb
cancor <- cc(contin_field, contin_forb)
rhocancor <- cancor$cor
rhocancor 
test <- p.asym(rhocancor, nobs, nvar_field, nvar_forb, tstat = "Hotelling") 
ccafieldforb <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Field", response = "Forb", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccafieldforb

## Fine x forb
cancor <- cc(contin_fine, contin_forb)
rhocancor <- cancor$cor
rhocancor 
test <- p.asym(rhocancor, nobs, nvar_fine, nvar_forb, tstat = "Hotelling") 
ccafineforb <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fine", response = "Forb", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccafineforb

# Summary statistics

## All dimensions included
forbvar_SI <- purrr::reduce(list(ccaregioforb, ccalandscapeforb, ccafieldforb, ccafineforb), dplyr::full_join)

## First dimension + other dimensions if marginally significant only
forbvar <- forbvar_SI |> 
  filter(dim == "dim1" | p.value < 0.1)


#### Canonical correlation analyses between explanatory sets and dominant beetle assemblage ####

# Number of observations (same for all sets)
nobs <- dim(contin_regional)[1]

# Number of variables in each set
nvar_regional <- length(select_if(regional_sc, is.numeric))
nvar_landscape <- length(select_if(landscape_sc, is.numeric))
nvar_field <- length(select_if(field_sc, is.numeric))
nvar_fine <- length(select_if(fine_sc, is.numeric))
nvar_beetle <- length(select_if(beetle, is.numeric))

## Regional x beetle
cancor <- cc(contin_regional, contin_beetle)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_regional, nvar_beetle, tstat = "Hotelling")
ccaregiobeetle <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Regional", response = "Beetle", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccaregiobeetle

## Landscape x beetle
cancor <- cc(contin_landscape, contin_beetle)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_landscape, nvar_beetle, tstat = "Hotelling")
ccalandscapebeetle <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Landscape", response = "Beetle", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccalandscapebeetle

## Field x beetle
cancor <- cc(contin_field, contin_beetle)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_field, nvar_beetle, tstat = "Hotelling")
ccafieldbeetle <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Field", response = "Beetle", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccafieldbeetle

## Fine x beetle
cancor <- cc(contin_fine, contin_beetle)
rhocancor <- cancor$cor
rhocancor
test <- p.asym(rhocancor, nobs, nvar_fine, nvar_beetle, tstat = "Hotelling")
ccafinebeetle <- as.data.frame(test) |> 
  mutate(rho = rhocancor, explanatory = "Fine", response = "Beetle", dim = "dim") |>  
  mutate(dim = paste(dim, row_number(), sep = ""))
ccafinebeetle

# Summary statistics

## All dimensions included
beetlevar_SI <- purrr::reduce(list(ccaregiobeetle, ccalandscapebeetle, ccafieldbeetle, ccafinebeetle), dplyr::full_join)

## First dimension + other dimensions if marginally significant only
beetlevar <- beetlevar_SI |> 
  filter(dim == "dim1" | p.value < 0.1)

## All significant community outputs
comvar <- purrr::reduce(list(grassvar, forbvar, beetlevar), dplyr::full_join)


#### RDA plotting function ####

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


#### Redundancy analyses between explanatory sets ####

# Region x landscape

## Redundancy analysis
rda <- rda(subset(landscape_sc, select = -c(SiteID)) ~ ., data = subset(regional_sc, select = -c(SiteID)))

## RDA plot
plotrda_regiolandscape <- ggrda(rda,
      txt = 4,
      ptslab = TRUE,
      addcol = "chartreuse4",
      addsize = 3.5,
      size = 2,
      arrow = 0.2,
      repel = TRUE,
      veccol = "cyan4",
      labcol = "cyan4",
      xlims = c(-1.1, 1.1),
      ylims = c(-1.1, 1.1)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = 0.9, image = "illustrations/Icons/icon_regional.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = -0.9, image = "illustrations/Icons/icon_landscape.png"),
    size = 0.2
  )
plotrda_regiolandscape

# plot <- autoplot(
#   rda,
#   axes = c(1, 2),
#   geom = c("point", "text"),
#   layers = c("species", "sites", "biplot", "centroids"),
#   arrows = TRUE,
#   const = 1
# ) +
#   geom_image(
#     # data = tibble(x = 1, y = 1),
#     aes(x = 1, y = 1, image = "illustrations/Icons/icon_regional.png")
#   )
# plot

ggsave("outputs/singleRDA/plotrda_regiolandscape.png", plot = plotrda_regiolandscape, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |>
  # "variance" output from cca is actually eigenval
  mutate(eigenval = Variance) |>
  # calculate proportion explained variance from eigvenval and inertia
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Residual variation using Kaiser-Guttman criterion
#resid <- as.data.frame(rda$CA$eig[rda$CA$eig > mean(rda$CA$eig)])

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdaregiolandscape <- rdatable |> 
  mutate(model = "RegionalxLandscape", explanatory = "Regional", response = "Landscape") |>
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Regional x field

## Redundancy analysis
rda <- rda(subset(field_sc, select = -c(SiteID)) ~ ., data = subset(regional_sc, select = -c(SiteID)))

## RDA plot
plotrda_regiofield <- ggrda(rda,
                               txt = 4,
                               ptslab = TRUE,
                               addcol = "darkgoldenrod3",
                               addsize = 3.5,
                               size = 1,
                               arrow = 0.2,
                               #repel = TRUE,
                               veccol = "cyan4",
                               labcol = "cyan4",
                                xlims = c(-1.4, 1.4),
                                ylims = c(-1.4, 1.4)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1.05, y = 1.15, image = "illustrations/Icons/icon_regional.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1.05, y = -1.15, image = "illustrations/Icons/icon_field.png"),
    size = 0.2
  )
plotrda_regiofield
ggsave("outputs/singleRDA/plotrda_regiofield.png", plot = plotrda_regiofield, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |>
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdaregiofield <- rdatable |>
  mutate(model = "RegionalxField", explanatory = "Regional", response = "Field") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Regional x fine

## Redundancy analysis
rda <- rda(subset(fine_sc, select = -c(SiteID)) ~ ., data = subset(regional_sc, select = -c(SiteID)))

## RDA plot
plotrda_regiofine <- ggrda(rda,
                               txt = 4,
                               ptslab = TRUE,
                               addcol = "darkred",
                               addsize = 3.5,
                               size = 1,
                               arrow = 0.2,
                               #repel = TRUE,
                               veccol = "cyan4",
                               labcol = "cyan4",
                               xlims = c(-1.1, 1.1),
                               ylims = c(-1.1, 1.1)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = 0.9, image = "illustrations/Icons/icon_regional.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = -0.9, image = "illustrations/Icons/icon_fine.png"),
    size = 0.2
  )
plotrda_regiofine
ggsave("outputs/singleRDA/plotrda_regiofine.png", plot = plotrda_regiofine, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdaregiofine <- rdatable |> 
  mutate(model = "RegionalxFine", explanatory = "Regional", response = "Fine") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Landscape x field

## Redundancy analysis
rda <- rda(subset(field_sc, select = -c(SiteID)) ~ ., data = subset(landscape_sc, select = -c(SiteID)))

## RDA plot
plotrda_landscapefield <- ggrda(rda,
                                   txt = 4,
                                   ptslab = TRUE,
                                   addcol = "darkgoldenrod3",
                                   addsize = 3.5,
                                   size = 1,
                                   arrow = 0.2,
                                   repel = TRUE,
                                   veccol = "chartreuse4",
                                   labcol = "chartreuse4",
                                   xlims = c(-1.1, 1.1),
                                   ylims = c(-1.1, 1.1)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = 0.9, image = "illustrations/Icons/icon_landscape.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = -0.9, image = "illustrations/Icons/icon_field.png"),
    size = 0.2
  )
plotrda_landscapefield
ggsave("outputs/singleRDA/plotrda_landscapefield.png", plot = plotrda_landscapefield, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapefield <- rdatable |> 
  mutate(model = "LandscapexField", explanatory = "Landscape", response = "Field") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Landscape x fine

## Redundancy analysis
rda <- rda(subset(fine_sc, select = -c(SiteID)) ~ ., data = subset(landscape_sc, select = -c(SiteID)))

## RDA plot
plotrda_landscapefine <- ggrda(rda,
                                  txt = 4,
                                   ptslab = TRUE,
                                   addcol = "darkred",
                                   addsize = 3.5,
                                   size = 1,
                                   arrow = 0.2,
                                   repel = TRUE,
                                   veccol = "chartreuse4",
                                   labcol = "chartreuse4",
                                 xlims = c(-1.2, 1.2),
                                 ylims = c(-1.2, 1.2)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = 1, image = "illustrations/Icons/icon_landscape.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = -1, image = "illustrations/Icons/icon_fine.png"),
    size = 0.2
  )
plotrda_landscapefine
ggsave("outputs/singleRDA/plotrda_landscapefine.png", plot = plotrda_landscapefine, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapefine <- rdatable |> 
  mutate(model = "LandscapexFine", explanatory = "Landscape", response = "Fine") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Field x fine

## Redundancy analysis
rda <- rda(subset(fine_sc, select = -c(SiteID)) ~ ., data = subset(field_sc, select = -c(SiteID)))

## Plot RDA
plotrda_fieldfine <- ggrda(rda,
                                 txt = 4,
                                 ptslab = TRUE,
                                 addcol = "darkred",
                                 addsize = 3.5,
                                 size = 1,
                                 arrow = 0.2,
                                 repel = TRUE,
                                 veccol = "darkgoldenrod3",
                                 labcol = "darkgoldenrod3",
                                   xlims = c(-1.1, 1.1),
                                   ylims = c(-1.1, 1.1)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = 0.9, image = "illustrations/Icons/icon_field.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = -0.9, image = "illustrations/Icons/icon_fine.png"),
    size = 0.2
  )
plotrda_fieldfine
ggsave("outputs/singleRDA/plotrda_fieldfine.png", plot = plotrda_fieldfine, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafieldfine <- rdatable |>
  mutate(model = "FieldxFine", explanatory = "Field", response = "Fine") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Summary statistics for RDA analyses between explanatory sets

## Stat summary
rdastat_envi <- purrr::reduce(list(rdaregiolandscape, rdaregiofield, rdaregiofine, rdalandscapefield, rdalandscapefine, rdafieldfine), dplyr::full_join)

## Plot summary
# rda_envi <- ggarrange(plotrda_regiolandscape, plotrda_regiofield, plotrda_regiofine, plotrda_landscapefield, plotrda_landscapefine, plotrda_fieldfine,
#                    labels = c("A", "B", "C", "D", "E", "F"),
#                    ncol = 2,
#                    nrow = 3)
# rda_envi
# ggsave("outputs/RDA_envi.png", plot = rda_envi, width = 15, height = 22, units = "cm", bg = "white")

#### RDA grass community data ####

# Regional x grass

# ## Remove space within plant names to avoid parse error in plot
# names(grass) <- gsub(" ", ".", names(grass))

## Redundancy analysis
rda <- rda(subset(grass, select = -c(SiteID)) ~ ., data = subset(regional_sc, select = -c(SiteID)))

## RDA plot
plotrda_regiograss <- ggrda(rda,
                             txt = 4,
                             ptslab = TRUE,
                             addcol = "black",
                             addsize = 3.5,
                             size = 1,
                             arrow = 0.2,
                             repel = TRUE,
                             veccol = "cyan4",
                             labcol = "cyan4",
                               xlims = c(-1.2, 1.2),
                               ylims = c(-1.2, 1.2)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = 1, image = "illustrations/Icons/icon_regional.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = -1, image = "illustrations/Icons/icon_grass.png"),
    size = 0.2
  )
plotrda_regiograss
ggsave("outputs/singleRDA/plotrda_regiograss.png", plot = plotrda_regiograss, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdaregiograss <- rdatable |>
  mutate(model = "RegionalxGrass", explanatory = "Regional", response = "Grass") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Landscape x grass

## Redundancy analysis
rda <- rda(subset(grass, select = -c(SiteID)) ~ ., data = subset(landscape_sc, select = -c(SiteID)))

## RDA plot
plotrda_landscapegrass <- ggrda(rda,
                                 txt = 4,
                                 ptslab = TRUE,
                                 addcol = "black",
                                 addsize = 3.5,
                                 size = 1,
                                 arrow = 0.2,
                                 repel = TRUE,
                                 veccol = "chartreuse4",
                                 labcol = "chartreuse4",
                             xlims = c(-1.1, 1.1),
                             ylims = c(-1.1, 1.1)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = 0.9, image = "illustrations/Icons/icon_landscape.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = -0.9, image = "illustrations/Icons/icon_grass.png"),
    size = 0.2
  )
plotrda_landscapegrass
ggsave("outputs/singleRDA/plotrda_landscapegrass.png", plot = plotrda_landscapegrass, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapegrass <- rdatable |>
  mutate(model = "LandscapexGrass", explanatory = "Landscape", response = "Grass") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Field x grass

## Redundancy analysis
rda <- rda(subset(grass, select = -c(SiteID)) ~ ., data = subset(field_sc, select = -c(SiteID)))

## RDA plot
plotrda_fieldgrass <- ggrda(rda,
                               txt = 4,
                               ptslab = TRUE,
                               addcol = "black",
                               addsize = 3.5,
                               size = 1,
                               arrow = 0.2,
                               repel = TRUE,
                               veccol = "darkgoldenrod3",
                               labcol = "darkgoldenrod3",
                               xlims = c(-1.3, 1.3),
                               ylims = c(-1.3, 1.3)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1, y = 1.05, image = "illustrations/Icons/icon_field.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1, y = -1.05, image = "illustrations/Icons/icon_grass.png"),
    size = 0.2
  )
plotrda_fieldgrass
ggsave("outputs/singleRDA/plotrda_fieldgrass.png", plot = plotrda_fieldgrass, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafieldgrass <- rdatable |>
  mutate(model = "FieldxGrass", explanatory = "Field", response = "Grass") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Fine x grass

## Redundancy analysis
rda <- rda(subset(grass, select = -c(SiteID)) ~ ., data = subset(fine_sc, select = -c(SiteID)))

## Plot RDA
plotrda_finegrass <- ggrda(rda,
                           txt = 4,
                           ptslab = TRUE,
                           addcol = "black",
                           addsize = 3.5,
                           size = 1,
                           arrow = 0.2,
                           repel = TRUE,
                           force = 2,
                           veccol = "darkred",
                           labcol = "darkred",
                           xlims = c(-1.2, 1.2),
                           ylims = c(-1.2, 1.2)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = 1, image = "illustrations/Icons/icon_fine.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = -1, image = "illustrations/Icons/icon_grass.png"),
    size = 0.2
  )
plotrda_finegrass
ggsave("outputs/singleRDA/plotrda_finegrass.png", plot = plotrda_finegrass, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafinegrass <- rdatable |>
  mutate(model = "FinexGrass", explanatory = "Fine", response = "Grass") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Summary statistics for RDA analyses between explanatory sets and dominant grass assemblage

## Stat summary
rdastat_grass <- purrr::reduce(list(rdaregiograss, rdalandscapegrass, rdafieldgrass, rdafinegrass), dplyr::full_join)

## Plot summary
# rda_grass <- ggarrange(plotrda_regiograss, plotrda_landscapegrass, plotrda_fieldgrass, plotrda_finegrass,
#                       labels = c("A", "B", "C", "D"),
#                       ncol = 2,
#                       nrow = 2)
# rda_grass
# ggsave("outputs/RDA_grass.png", plot = rda_grass, width = 15, height = 15, units = "cm", bg = "white")


#### RDA forb community data ####

# Regional x forb

## Redundancy analysis
rda <- rda(subset(forb, select = -c(SiteID)) ~ ., data = subset(regional_sc, select = -c(SiteID)))

## Plot RDA
plotrda_regioforb <- ggrda(rda,
                            txt = 4,
                            ptslab = TRUE,
                            addcol = "black",
                            addsize = 3.5,
                            size = 1,
                            arrow = 0.2,
                            # repel = TRUE,
                            veccol = "cyan4",
                            labcol = "cyan4",
                             xlims = c(-1.2, 1.2),
                             ylims = c(-1.2, 1.2)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = 1, image = "illustrations/Icons/icon_regional.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = -1, image = "illustrations/Icons/icon_forb.png"),
    size = 0.2
  )
plotrda_regioforb
ggsave("outputs/singleRDA/plotrda_regioforb.png", plot = plotrda_regioforb, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdaregioforb <- rdatable |>
  mutate(model = "RegionalxForb", explanatory = "Regional", response = "Forb") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Landscape x forb

## Redundancy analysis
rda <- rda(subset(forb, select = -c(SiteID)) ~ ., data = subset(landscape_sc, select = -c(SiteID)))

## Plot RDA
plotrda_landscapeforb <- ggrda(rda,
                                txt = 4,
                                ptslab = TRUE,
                                addcol = "black",
                                addsize = 3.5,
                                size = 1,
                                arrow = 0.2,
                                #repel = TRUE,
                                veccol = "chartreuse4",
                                labcol = "chartreuse4",
                                 xlims = c(-1.1, 1.1),
                                 ylims = c(-1.1, 1.1)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = 0.9, image = "illustrations/Icons/icon_landscape.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = -0.9, image = "illustrations/Icons/icon_forb.png"),
    size = 0.2
  )
plotrda_landscapeforb
ggsave("outputs/singleRDA/plotrda_landscapeforb.png", plot = plotrda_landscapeforb, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapeforb <- rdatable |>
  mutate(model = "LandscapexForb", explanatory = "Landscape", response = "Forb") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Field x forb

## Redundancy analysis
rda <- rda(subset(forb, select = -c(SiteID)) ~ ., data = subset(field_sc, select = -c(SiteID)))

## RDA plot
plotrda_fieldforb <- ggrda(rda,
                              txt = 4,
                              ptslab = TRUE,
                              addcol = "black",
                              addsize = 3.5,
                              size = 1,
                              arrow = 0.2,
                              #repel = TRUE,
                              veccol = "darkgoldenrod3",
                              labcol = "darkgoldenrod3",
                               xlims = c(-1.2, 1.2),
                               ylims = c(-1.2, 1.2)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = 1, image = "illustrations/Icons/icon_field.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = -1, image = "illustrations/Icons/icon_forb.png"),
    size = 0.2
  )
plotrda_fieldforb
ggsave("outputs/singleRDA/plotrda_fieldforb.png", plot = plotrda_fieldforb, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafieldforb <- rdatable |>
  mutate(model = "FieldxForb", explanatory = "Field", response = "Forb") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Fine x forb

## Redundancy analysis
rda <- rda(subset(forb, select = -c(SiteID)) ~ ., data = subset(fine_sc, select = -c(SiteID)))

## RDA plot
plotrda_fineforb <- ggrda(rda,
                              txt = 4,
                              ptslab = TRUE,
                              addcol = "black",
                              addsize = 3.5,
                              size = 1,
                              arrow = 0.2,
                              repel = TRUE,
                              veccol = "darkred",
                              labcol = "darkred",
                              xlims = c(-1.2, 1.2),
                              ylims = c(-1.2, 1.2)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = 1, image = "illustrations/Icons/icon_fine.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.9, y = -1, image = "illustrations/Icons/icon_forb.png"),
    size = 0.2
  )
plotrda_fineforb
ggsave("outputs/singleRDA/plotrda_fineforb.png", plot = plotrda_fineforb, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafineforb <- rdatable |>
  mutate(model = "FinexForb", explanatory = "Fine", response = "Forb") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Summary statistics for RDA analyses between explanatory sets and dominant forb assemblage

## Stat summary
rdastat_forb <- purrr::reduce(list(rdaregioforb, rdalandscapeforb, rdafieldforb, rdafineforb), dplyr::full_join)

## Plot summary
# rda_forb <- ggarrange(plotrda_regioforb, plotrda_landscapeforb, plotrda_fieldforb, plotrda_fineforb,
#                        labels = c("A", "B", "C", "D"),
#                        ncol = 2,
#                        nrow = 2)
# rda_forb
# ggsave("outputs/RDA_forb.png", plot = rda_forb, width = 15, height = 15, units = "cm", bg = "white")

#### RDA beetle community data ####

# Regional x beetle

## Redundancy analysis
rda <- rda(subset(beetle, select = -c(SiteID)) ~ ., data = subset(regional_sc, select = -c(SiteID)))

## RDA plot
plotrda_regiobeetle <- ggrda(rda,
                              txt = 4,
                              ptslab = TRUE,
                              addcol = "black",
                              addsize = 3.5,
                              size = 1,
                              arrow = 0.2,
                              #repel = TRUE,
                              veccol = "cyan4",
                              labcol = "cyan4",
                            xlims = c(-1.8, 1.8),
                            ylims = c(-1.8, 1.8)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1.4, y = 1.5, image = "illustrations/Icons/icon_regional.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1.4, y = -1.5, image = "illustrations/Icons/icon_beetle.png"),
    size = 0.2
  )
plotrda_regiobeetle
ggsave("outputs/singleRDA/plotrda_regiobeetle.png", plot = plotrda_regiobeetle, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdaregiobeetle <- rdatable |>
  mutate(model = "RegionalxBeetle", explanatory = "Regional", response = "Beetle") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Landscape x beetle

## Redundancy analysis
rda <- rda(subset(beetle, select = -c(SiteID)) ~ ., data = subset(landscape_sc, select = -c(SiteID)))

## RDA plot
plotrda_landscapebeetle <- ggrda(rda,
                                  txt = 4,
                                  ptslab = TRUE,
                                  addcol = "black",
                                  addsize = 3.5,
                                  size = 1,
                                  arrow = 0.2,
                                  repel = TRUE,
                                 force = 1,
                                  veccol = "chartreuse4",
                                  labcol = "chartreuse4",
                                xlims = c(-1.4, 1.4),
                                ylims = c(-1.4, 1.4)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1.05, y = 1.15, image = "illustrations/Icons/icon_landscape.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1.05, y = -1.15, image = "illustrations/Icons/icon_beetle.png"),
    size = 0.2
  )
plotrda_landscapebeetle
ggsave("outputs/singleRDA/plotrda_landscapebeetle.png", plot = plotrda_landscapebeetle, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdalandscapebeetle <- rdatable |>
  mutate(model = "LandscapexBeetle", explanatory = "Landscape", response = "Beetle") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Field x beetle

## Redundancy analysis
rda <- rda(subset(beetle, select = -c(SiteID)) ~ ., data = subset(field_sc, select = -c(SiteID)))

## RDA plot
plotrda_fieldbeetle <- ggrda(rda,
                                txt = 4,
                                ptslab = TRUE,
                                addcol = "black",
                                addsize = 3.5,
                                size = 1,
                                arrow = 0.2,
                                repel = TRUE,
                                veccol = "darkgoldenrod3",
                                labcol = "darkgoldenrod3",
                              xlims = c(-1.3, 1.3),
                              ylims = c(-1.3, 1.3)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1, y = 1.05, image = "illustrations/Icons/icon_field.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 1, y = -1.05, image = "illustrations/Icons/icon_beetle.png"),
    size = 0.2
  )
plotrda_fieldbeetle
ggsave("outputs/singleRDA/plotrda_fieldbeetle.png", plot = plotrda_fieldbeetle, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafieldbeetle <- rdatable |>
  mutate(model = "FieldxBeetle", explanatory = "Field", response = "Beetle") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Fine x beetle

## Redundancy analysis
rda <- rda(subset(beetle, select = -c(SiteID)) ~ ., data = subset(fine_sc, select = -c(SiteID)))

## RDA plot
plotrda_finebeetle <- ggrda(rda,
                                txt = 4,
                                ptslab = TRUE,
                                addcol = "black",
                                addsize = 3.5,
                                size = 1,
                                arrow = 0.2,
                                repel = TRUE,
                            force = 2,
                                veccol = "darkred",
                                labcol = "darkred",
                              xlims = c(-1.1, 1.1),
                              ylims = c(-1.1, 1.1)) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = 0.9, image = "illustrations/Icons/icon_fine.png"),
    size = 0.2
  ) +
  geom_image(
    data = tibble(x = 1, y = 1),
    aes(x = 0.85, y = -0.9, image = "illustrations/Icons/icon_beetle.png"),
    size = 0.2
  )
plotrda_finebeetle
ggsave("outputs/singleRDA/plotrda_finebeetle.png", plot = plotrda_finebeetle, width = 6, height = 6, units = "cm", bg = "white")

## Adjusted variance explained by the RDA model
rsq <- as.data.frame(RsquareAdj(rda)) |> 
  mutate(type = "model", dim = "Model")

## ANOVA permutation test on RDA model
testrda <- as.data.frame(anova.cca(rda)) 
testrda <- testrda |> 
  mutate(dim = rownames(testrda), type = "model") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi)

## ANOVA permutation test on individual axes
testaxis <- as.data.frame(anova.cca(rda, by = "axis"))
testaxis <- testaxis |> 
  mutate(dim = row.names(testaxis), type = "axis") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## ANOVA permutation test on margins (relative importance terms)
testerm <- as.data.frame(anova.cca(rda, by = "margin")) 
testerm <- testerm |> 
  mutate(dim = row.names(testerm), type = "margin") |> 
  mutate(eigenval = Variance) |> 
  mutate(Variance = eigenval/summary(rda)$tot.chi) |> 
  filter(dim != "Residual")

## Summary statistics
rdatable <- purrr::reduce(list(testrda, rsq, testaxis, testerm), dplyr::full_join)
names(rdatable) <- gsub("\\(", "", names(rdatable))
names(rdatable) <- gsub("\\)", "", names(rdatable))
names(rdatable) <- gsub(">", "", names(rdatable))
rdafinebeetle <- rdatable |>
  mutate(model = "FinexBeetle", explanatory = "Fine", response = "Beetle") |> 
  filter(type == "model" | dim == "RDA1" | dim == "RDA2" | type == "margin")

# Summary statistics for RDA analyses between explanatory sets and dominant forb assemblage

# Stat summary
rdastat_beetle <- purrr::reduce(list(rdaregiobeetle, rdalandscapebeetle, rdafieldbeetle, rdafinebeetle), dplyr::full_join)
# 
# ## Plot summary
# rda_beetle <- ggarrange(plotrda_regiobeetle, plotrda_landscapebeetle, plotrda_fieldbeetle, plotrda_finebeetle,
#                         labels = c("A", "B", "C", "D"),
#                         font.label = list(size = 12),
#                         ncol = 2,
#                         nrow = 2)
# rda_beetle
# ggsave("outputs/RDA_beetle.png", plot = rda_beetle, width = 15, height = 15, units = "cm", bg = "white")

#### Significant RDA results ####

# Significant models
# filter(rdastat_envi, PrF < 0.1 & type == "model")
# filter(rdastat_grass, PrF < 0.1 & type == "model")
# filter(rdastat_forb, PrF < 0.1 & type == "model")
# filter(rdastat_beetle, PrF < 0.1 & type == "model")
# 
# # Summary significant RDA stat
# rdastat_all <- purrr::reduce(list(rdafjordlandscape, rdagrazinglocenvi, rdalocenvigrass, rdalocenviforb), dplyr::full_join)

# Plot
# rdall <- ggarrange(plotrda_regiolandscape, plotrda_fieldfine, plotrda_finegrass, plotrda_fineforb, plotrda_landscapebeetle,
#                    labels = c("A", "B", "C", "D", "E"),
#                    ncol = 2,
#                    nrow = 3,
#                    font.label = list(size = 12),
#                    hjust = -2.7)
# rdall
# ggsave("outputs/RDAresults.png", plot = rdall, width = 16, height = 20, units = "cm", bg = "white")

#
## Distribution residuals for significant relationships with terms

# Ordination community data
# DCA_grass <- decorana(contin_grass)
# DCA_grass # Axis length of 1.9 -> run PCA
# PCA_grass <- prcomp(contin_grass)
# PCA_grass
# biplot(PCA_grass)
# DCA_forb <- decorana(contin_forb)
# DCA_forb # Axis length of 3.1 -> keep DCA
# plot(DCA_forb)
# 
# # Extraction grass species scores from PCA
# scores_grass <- as.data.frame(scores(PCA_grass, choices = c(1,2), display = "sites"))
# names(scores_grass) <- gsub("PC1", "GrassPCA1", names(scores_grass))
# names(scores_grass) <- gsub("PC2", "GrassPCA2", names(scores_grass))
# scores_grass <- scores_grass |>
#   mutate(SiteID = row.names(scores_grass)) # PlotID as column for future binding
# 
# # Extraction forb species scores from DCA
# scores_forb <- as.data.frame(scores(DCA_forb, choices = c(1,2), display = "sites"))
# names(scores_forb) <- gsub("DCA1", "ForbDCA1", names(scores_forb))
# names(scores_forb) <- gsub("DCA2", "ForbDCA2", names(scores_forb))
# scores_forb <- scores_forb |>
#   mutate(SiteID = row.names(scores_forb)) # PlotID as column for future binding
# 
# # Table binding
# allvar <- purrr::reduce(list(scores_grass, scores_forb, fjordsys_sc, landscape_sc, grazing_sc, locenvi_sc), dplyr::left_join)

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