# TradMod WP2 - NBR survey grazer farmer outreach
#Description of the data
#Date
#Who
#Project
#Funding
#Place

#### PACKAGE LOADING ####

library(tidyverse) # R language
library(purrr) # Data manipulation: function "reduce" to bind several tables at the same time
library(ggplot2) # Visual representation
library(GGally) # Extension ggplot
library(xlsx) #read & turn into xl
library(vegan) # Species richness
library(tibble) # row names to columns


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
mesobio_full <- read.csv("data/cleandata/NBR_FullMesobio.csv", sep=",") # Clean mesofauna community data


#### DATA PREPARATION ####

#
## Variables extraction

siteinfo_outreach <- subset(siteinfo_full, select = c(SiteID, Type_livestock, Habitat, EPSG.25832_X, EPSG.25832_Y))
climate_outreach <- climate_full
area20x20_outreach <- subset(area20x20_full, select = c(SiteID, Elevation_max))
VGrichness_outreach <- subset(groundcover_full, select = c(SiteID, Plant_species_richness))
soilbulk_outreach <- subset(soilbulk_full, select = c(SiteID, BD))
soilchem_outreach <- subset(soilchem_full, select = c(SiteID, LOI, pH, P.Al_mg.100g, K.Al_mg.100g, Mg.Al_mg.100g, Ca.Al_mg.100g, Na.Al_mg.100g, TotalN_percentDM, SoilDensity_kg.L, Humus_percentDM))
vege_outreach <- subset(vege_full, select = c(SiteID, Species, Abundance))
beetle_outreach <- subset(beetle_full, select = c(SiteID, BeetleFamilies, BeetleFam_abundance))

#
## Summarise data at site level

# Plant richness
VGrichness_outreach <- VGrichness_outreach %>% 
  group_by(SiteID) %>% 
  summarise(MeanPlantRichness = mean(Plant_species_richness, na.rm=TRUE))
write.xlsx(VGrichness_outreach, "outreach/OutreachPlantRichness.xlsx")

# Bulk density
soilbulk_outreach <- soilbulk_outreach %>% 
  group_by(SiteID) %>% 
  summarise(MeanBulkDensity = mean(BD, na.rm=TRUE))

# Soil chemistry
soilchem_outreach <- soilchem_outreach %>%
  group_by(SiteID) %>% 
  summarise(LOI = mean(LOI), pH = mean(pH), P = mean(P.Al_mg.100g), K = mean(K.Al_mg.100g), Mg = mean(Mg.Al_mg.100g), Ca = mean(Ca.Al_mg.100g), Na = mean(Na.Al_mg.100g), TotN = mean(TotalN_percentDM), JORDTETTLEIK = mean(SoilDensity_kg.L), HUMUS = mean(Humus_percentDM)) #missing OV1 (Oygarden) and UC1

# Plant abundance cover
vege_outreach <- vege_outreach |> 
  group_by(SiteID, Species) |> 
  summarise(PlantSp_cover = mean(Abundance))

# Beetle abundance
beetle_outreach <- beetle_outreach %>% 
  group_by(SiteID, BeetleFamilies) %>% 
  summarise(BeetleFam_abundance = sum(BeetleFam_abundance, na.rm = TRUE))

# Mesofauna
mesobio_outreach <- mesobio_full |> 
  group_by(SiteID) |> 
  summarise(MeanAcari = mean(Acari), MeanCollembola = mean(Collembola))

#
## Plant species richness

vege_outreach <- vege_full |> 
  group_by(SiteID, Species) |> 
  summarise(PlantSp_cover = mean(Abundance))
  
vege_contin <- xtabs(formula = PlantSp_cover ~ SiteID + Species, data = vege_outreach)
speciesrichness <- specnumber(vege_contin)
speciesrichness <- as.data.frame(speciesrichness)
speciesrichness <- tibble::rownames_to_column(speciesrichness, "SiteID")

#  
## Summary dataset for GIS mapping

Map_All <- purrr::reduce(list(siteinfo_outreach, climate_outreach, area20x20_outreach, landscape_full, VGrichness_outreach, speciesrichness, soilbulk_outreach, soilchem_outreach, beetle_outreach, mesobio_outreach), dplyr::left_join)


### Species richness ####

#
## Boxplot plant species richness against animal for all habitats

Map_All <- Map_All |> 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "cow" = "ku")) |> 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "sheep" = "sau")) |> 
  mutate(Type_livestock = dplyr::recode(Type_livestock, "goat" = "geit")) 
  

# Total plant richness on site
totsrxanimal <- ggplot(Map_All, aes(x = as.factor(Type_livestock), y = speciesrichness, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Antall artar\n") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
totsrxanimal

# Mean plant richness per square meter
meansrxanimal <- ggplot(Map_All, aes(x = as.factor(Type_livestock), y = MeanPlantRichness, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Mean number of species per square meter") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")
meansrxanimal

#
## Boxplot plant species richness against habitat for sheep only

# Total plant richness on heathland site
sheep_totsrxhabitat <- ggplot(filter(Map_All, Type_livestock == "sau"), aes(x = as.factor(Habitat), y = speciesrichness, fill = Habitat)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Antall artar\n") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.position = "none")
sheep_totsrxhabitat

#
## Boxplot plant species richness against animal for grassland only

# Total plant richness on grassland site
grass_totsrxanimal <- ggplot(filter(Map_All, Habitat == "permanent grassland"), aes(x = as.factor(Type_livestock), y = speciesrichness, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Antall artar\n") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "none")
grass_totsrxanimal
ggsave("outreach/grass_totsrxanimal.png", plot = grass_totsrxanimal, width = 15, height = 17, units = "cm")

# Mean plant richness per square meter in grasslands
grass_meansrxanimal <- ggplot(filter(Map_All, Habitat == "permanent grassland"), aes(x = as.factor(Type_livestock), y = MeanPlantRichness, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Mean number of species per square meter\n") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "none")
grass_meansrxanimal
ggsave("outreach/grass_meansrxanimal.png", plot = grass_meansrxanimal, width = 15, height = 17, units = "cm")



#
## Boxplot beetle abundance against animal for grassland only

beetle_outreach <- beetle_full %>% 
  group_by(SiteID, BeetleFamilies) %>% 
  summarise(BeetleFam_abundance = sum(BeetleFam_abundance, na.rm = TRUE)) |> 
  pivot_wider(names_from = BeetleFamilies, values_from = BeetleFam_abundance)
beetle_outreach <- purrr::reduce(list(siteinfo_outreach, area20x20_outreach, beetle_outreach, landuse_full, landscape_full), dplyr::left_join)

# Total beetle x animals
beetlexanimal <- ggplot(filter(beetle_outreach, Habitat == "permanent grassland"), aes(x = as.factor(Type_livestock), y = Ptiliidae, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Total number of species\n") +
  ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.position = "none")
beetlexanimal
#ggsave("outreach/beetle_totsrxanimal.png", plot = beetle_totsrxanimal, width = 15, height = 17, units = "cm")

#
## Boxplot mesofauna abundance in sheep fields against habitat

# Springtails
springtailxhabitat <- ggplot(filter(Map_All, Type_livestock == "sau"), aes(x = as.factor(Habitat), y = MeanCollembola, fill = Habitat)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Total number of species\n") +
  #ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 20),
        #axis.text = element_text(size = 20),
        legend.position = "none")
springtailxhabitat

# Acari
acarixhabitat <- ggplot(filter(Map_All, Type_livestock == "sau"), aes(x = as.factor(Habitat), y = MeanAcari, fill = Habitat)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Total number of species\n") +
  #ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 20),
        #axis.text = element_text(size = 20),
        legend.position = "none")
acarixhabitat

#
## Boxplot mesofauna abundance in grasslands against animal

# Springtails
springtailxlivestock <- ggplot(filter(Map_All, SiteID == "IC1"| SiteID == "IC2"| SiteID == "IG1"| SiteID == "IG2"| SiteID == "IS1"| SiteID == "IS2"), aes(x = as.factor(Type_livestock), y = MeanCollembola, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Total number of species\n") +
  #ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 20),
        #axis.text = element_text(size = 20),
        legend.position = "none")
springtailxlivestock

# Acari
acarixlivestock <- ggplot(filter(Map_All, SiteID == "IC1"| SiteID == "IC2"| SiteID == "IG1"| SiteID == "IG2"| SiteID == "IS1"| SiteID == "IS2"), aes(x = as.factor(Type_livestock), y = MeanAcari, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Total number of species\n") +
  #ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 20),
        #axis.text = element_text(size = 20),
        legend.position = "none")
acarixlivestock

#
## Boxplot nutrients in grasslands against animal

# Nitrogen
nitrogenxlivestock <- ggplot(filter(Map_All, Habitat == "permanent grassland"), aes(x = as.factor(Type_livestock), y = TotN, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Total number of species\n") +
  #ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 20),
        #axis.text = element_text(size = 20),
        legend.position = "none")
nitrogenxlivestock

# Phosphorus
phosphorusxlivestock <- ggplot(filter(Map_All, Habitat == "permanent grassland"), aes(x = as.factor(Type_livestock), y = P, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  ylab("Total number of species\n") +
  #ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 20),
        #axis.text = element_text(size = 20),
        legend.position = "none")
phosphorusxlivestock

# Mean bulk density
BDxlivestock <- ggplot(filter(Map_All, Habitat == "permanent grassland"), aes(x = as.factor(Type_livestock), y = MeanBulkDensity, fill = Type_livestock)) +
  geom_boxplot(alpha = 0.3) +
  scale_fill_grey(start = 0.1, end = 0.9) +
  #ylab("Total number of species\n") +
  #ylim(0,1000) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_blank(),
        #axis.title.y = element_text(size = 20),
        #axis.text = element_text(size = 20),
        legend.position = "none")
BDxlivestock

# Replace NA by zeros -> needed for QGIS, otherwise it treats the layer as character and does not want to apply a graduated aesthetics
Map_All[is.na(Map_All)] <- 0

# csv for QGIS map -> should write a new script !
write.csv(Map_All, file = "outreach/NBRenv_AttributeTable.csv")


#### Fjord x landscape ####

# Make table
FjoLand <- purrr::reduce(list(siteinfo_full, area20x20_full, landscape_full), dplyr::left_join)

# Trend Forest x elevation
ElevForest <- ggplot(filter(FjoLand, Habitat == "permanent grassland"), aes(x = Elevation_max, y = ProductiveForest_percent)) +
  geom_point() +
  stat_smooth(method = "lm")
ElevForest

# Trend Infield x elevation
ElevInfield <- ggplot(filter(FjoLand, Habitat == "permanent grassland"), aes(x = Elevation_max, y = Infield_percent)) +
  geom_point() +
  stat_smooth(method = "lm")
ElevInfield

# Trend Outfield x elevation
ElevOutfield <- ggplot(filter(FjoLand, Habitat == "permanent grassland"), aes(x = Elevation_max, y = Outfield_percent)) +
  geom_point() +
  stat_smooth(method = "lm")
ElevOutfield

# Trend freshwater x elevation
ElevWetland <- ggplot(filter(FjoLand, Habitat == "permanent grassland"), aes(x = Elevation_max, y = Freshwater_percent)) +
  geom_point() +
  stat_smooth(method = "lm")
ElevWetland

#### Plant composition ####

vege_outreach <- vege_full |>
  group_by(Species) |>
  summarise(PlantSp_cover = mean(Abundance)) |>
  #filter(SiteID == "US6") |>
  filter(PlantSp_cover>0) |>
  arrange(desc(PlantSp_cover))
summary(vege_outreach)
head(vege_outreach)

#
## Species trends against gradients

# Plant vegetation table
mainvgsp <- subset(vege_full,
                Species == "Agrostis capillaris" |
                  Species == "Festuca rubra" | 
                  Species == "Holcus lanatus" |
                  Species == "Trifolium repens" |
                  Species == "Rhytidiadelphus squarrosus"
              )
mainvgsp <- mainvgsp |> 
  group_by(SiteID, Species) |> 
  summarise(Abundance = mean(Abundance)) #|> 
  #pivot_wider(names_from = Species, values_from = Abundance)
mainvgsp <- left_join(mainvgsp, Map_All)
mainvgsp <- filter(mainvgsp, Habitat == "permanent grassland")

# Agrostis capillaris
Agcap <- ggplot(filter(mainvgsp, Species == "Agrostis capillaris"), aes(x = ProductiveForest_percent, y = Abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
Agcap

# Festuca rubra
FesRu <- ggplot(filter(mainvgsp, Species == "Festuca rubra"), aes(x = Infield_percent, y = Abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
FesRu

# Holcus lanatus
HolLa <- ggplot(filter(mainvgsp, Species == "Holcus lanatus"), aes(x = Outfield_percent, y = Abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
HolLa

# Trifolium repens
trirep <- ggplot(filter(mainvgsp, Species == "Trifolium repens"), aes(x = ProductiveForest_percent, y = Abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
trirep


# Rhytidiadelphus sq
rhysq <- ggplot(filter(mainvgsp, Species == "Rhytidiadelphus squarrosus"), aes(x = ProductiveForest_percent, y = Abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
rhysq


#### Beetle composition ####

# Dataset
beetle_outreach <- beetle_full |> 
  group_by(SiteID, BeetleFamilies) |> 
  summarise(BeetleFam_abundance = sum(BeetleFam_abundance, na.rm = TRUE)) |> 
  #filter(SiteID == "US6") |> 
  arrange(desc(BeetleFam_abundance))
#head(beetle_outreach)
#sum(beetle_outreach$BeetleFam_abundance)
beetle_outreach <- left_join(beetle_outreach, soilchem_outreach)

# Staph
staph <- ggplot(filter(beetle_outreach, BeetleFamilies == "Staphylinidae"), aes(x = TotN, y = BeetleFam_abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
staph

# Hydro
hydro <- ggplot(filter(beetle_outreach, BeetleFamilies == "Hydrophilidae"), aes(x = TotN, y = BeetleFam_abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
hydro

# Ptili
ptili <- ggplot(filter(beetle_outreach, BeetleFamilies == "Ptiliidae"), aes(x = TotN, y = BeetleFam_abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
ptili

# Scara
scara <- ggplot(filter(beetle_outreach, BeetleFamilies == "Scarabaeidae"), aes(x = TotN, y = BeetleFam_abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
scara

# Carab
carab <- ggplot(filter(beetle_outreach, BeetleFamilies == "Carabidae"), aes(x = pH, y = BeetleFam_abundance)) +
  geom_point() +
  stat_smooth(method = "lm")
carab


#### Soil mesofauna ####

mesobio_outreach <- left_join(mesobio_outreach, soilchem_outreach)

# Acari
acari <- ggplot(mesobio_outreach, aes(x = P, y = MeanAcari)) +
  geom_point() +
  stat_smooth(method = "lm")
acari

# Collembola
collemb <- ggplot(mesobio_outreach, aes(x = P, y = MeanCollembola)) +
  geom_point() +
  stat_smooth(method = "lm")
collemb



#### Soil analysis violins ####

#
## Filtering infields and outfields

#Rename factors
names(soilchem_outreach)<-gsub("P.Al_mg.100g", "FOSFOR", names(soilchem_outreach))
names(soilchem_outreach)<-gsub("Ca.Al_mg.100g", "KALSIUM", names(soilchem_outreach))
names(soilchem_outreach)<-gsub("Mg.Al_mg.100g", "MAGNESIUM", names(soilchem_outreach))
names(soilchem_outreach)<-gsub("TotalN_percentDM", "NITROGEN", names(soilchem_outreach))
names(soilchem_outreach)<-gsub("LOI", "GLØDETAP", names(soilchem_outreach))
write.xlsx(soilchem_outreach, "outreach/OutreachChem.xlsx")

# Site selection
grassland <- siteinfo_outreach |>  
  filter(Habitat == "permanent grassland") |> 
  dplyr::select(SiteID, Type_livestock)
coastalheath <- siteinfo_outreach |>  
  filter(Habitat == "coastal heathland") |> 
  dplyr::select(SiteID, Type_livestock)
alpineheath <- siteinfo_outreach |>  
  filter(Habitat == "subalpine heathland") |> 
  dplyr::select(SiteID, Type_livestock)

# Extraction in other datasets
soilchem_grassland <- filter(soilchem_outreach, SiteID %in% grassland$SiteID)
soilchem_coastalheath <- filter(soilchem_outreach, SiteID %in% coastalheath$SiteID)
soilchem_alpineheath <- filter(soilchem_outreach, SiteID %in% alpineheath$SiteID)

#
## Bulk density

# All sites
PlotBD_all <- subset(Map_All, select = c(MeanBulkDensity)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinBD_all <- ggplot(PlotBD_all, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotBD_all$Rates), color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinBD_all
ggsave("outreach/NBRFarms_BDall.png", ViolinBD_all, bg = "transparent", width = 18, height = 6, units = "cm")

# Grasslands only
PlotBD_grassland <- subset(filter(Map_All, Habitat == "permanent grassland"), select = c(MeanBulkDensity)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinBD_grassland <- ggplot(PlotBD_grassland, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotBD_grassland$Rates), color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinBD_grassland
ggsave("outreach/NBRFarms_BDgrassland.png", ViolinBD_grassland, bg = "transparent", width = 18, height = 6, units = "cm")

#
## Soil density

# All sites
PlotSD_all <- subset(soilchem_outreach, select = c(JORDTETTLEIK)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinSD_all <- ggplot(PlotSD_all, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotSD_all$Rates), color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinSD_all
ggsave("outreach/NBRFarms_BDall.png", ViolinBD_all, bg = "transparent", width = 18, height = 6, units = "cm")

# Grasslands only
PlotSD_grassland <- subset(filter(Map_All, Habitat == "permanent grassland"), select = c(JORDTETTLEIK)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinSD_grassland <- ggplot(PlotSD_grassland, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotSD_grassland$Rates), color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinSD_grassland
ggsave("outreach/NBRFarms_BDgrassland.png", ViolinBD_grassland, bg = "transparent", width = 18, height = 6, units = "cm")

#
## Humus content

# All sites
PlotHumus_all <- subset(soilchem_outreach, select = c(HUMUS)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinHumus_all <- ggplot(PlotHumus_all, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotHumus_all$Rates), color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinHumus_all
ggsave("outreach/NBRFarms_Humusall.png", ViolinHumus_all, bg = "transparent", width = 18.5, height = 6, units = "cm")

# Grasslands only
PlotSD_grassland <- subset(filter(Map_All, Habitat == "permanent grassland"), select = c(JORDTETTLEIK)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinSD_grassland <- ggplot(PlotSD_grassland, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotSD_grassland$Rates), color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinSD_grassland
ggsave("outreach/NBRFarms_BDgrassland.png", ViolinBD_grassland, bg = "transparent", width = 18, height = 6, units = "cm")


#
## Infields violins

# Total N
PlotN_grassland <- subset(soilchem_grassland, select = c(NITROGEN)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinN_grassland <- ggplot(PlotN_grassland, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=0.82, color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinN_grassland
ggsave("outreach/NBRFarms_NGrassland.png", ViolinN_grassland, bg = "transparent", width = 18, height = 6, units = "cm")

# pH
PlotpH_grassland <- subset(soilchem_grassland, select = c(pH)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinpH_grassland <- ggplot(PlotpH_grassland, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  geom_hline(yintercept=4.83, color="chartreuse3", size=5) +
  scale_fill_grey() +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinpH_grassland
ggsave("outreach/NBRFarms_pHGrassland.png", ViolinpH_grassland, bg = "transparent", width = 18.5, height = 6, units = "cm")

# Phosphorus
PlotP_grassland <- subset(soilchem_grassland, select = c(FOSFOR)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinP_grassland <- ggplot(PlotP_grassland, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=4.00, color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinP_grassland
ggsave("outreach/NBRFarms_PhosphorusGrassland.png", ViolinP_grassland, bg = "transparent", width = 18.5, height = 6, units = "cm")

# Magnesium
PlotMg_grassland <- subset(soilchem_grassland, select = c(MAGNESIUM)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinMg_grassland <- ggplot(PlotMg_grassland, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=12.33, color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinMg_grassland
ggsave("outreach/NBRFarms_MagGrassland.png", ViolinMg_grassland, bg = "transparent", width = 18.5, height = 6, units = "cm")

# LOI
PlotLOI_grassland <- subset(soilchem_grassland, select = c(GLØDETAP)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinLOI_grassland <- ggplot(PlotLOI_grassland, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=28.13, color="chartreuse3", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinLOI_grassland
ggsave("outreach/NBRFarms_LOIGrassland.png", ViolinLOI_grassland, bg = "transparent", width = 18.5, height = 6, units = "cm")

#
## All violins

# Total N
PlotN_all <- subset(soilchem_outreach, select = c(NITROGEN)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinN_all <- ggplot(PlotN_all, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotN_all$Rates), color="chartreuse3", size=5) +
  #geom_hline(yintercept=median(PlotN_all$Rates), color="cornflowerblue", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinN_all
ggsave("outreach/NBRFarms_NAll.png", ViolinN_all, bg = "transparent", width = 18, height = 6, units = "cm")

# pH
PlotpH_all <- subset(soilchem_outreach, select = c(pH)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinpH_all <- ggplot(PlotpH_all, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  geom_hline(yintercept=mean(PlotpH_all$Rates), color="chartreuse3", size=5) +
  #geom_hline(yintercept=median(PlotpH_all$Rates), color="cornflowerblue", size=5) +
  scale_fill_grey() +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinpH_all
ggsave("outreach/NBRFarms_pHAll.png", ViolinpH_all, bg = "transparent", width = 18.5, height = 6, units = "cm")

# Phosphorus
PlotP_all <- subset(soilchem_outreach, select = c(FOSFOR)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinP_all <- ggplot(PlotP_all, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotP_all$Rates), color="chartreuse3", size=5) +
  #geom_hline(yintercept=median(PlotP_all$Rates), color="cornflowerblue", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinP_all
ggsave("outreach/NBRFarms_PhosphorusAll.png", ViolinP_all, bg = "transparent", width = 18.5, height = 6, units = "cm")

# Magnesium
PlotMg_all <- subset(soilchem_outreach, select = c(MAGNESIUM)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinMg_all <- ggplot(PlotMg_all, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotMg_all$Rates), color="chartreuse3", size=5) +
  #geom_hline(yintercept=median(PlotMg_all$Rates), color="cornflowerblue", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinMg_all
ggsave("outreach/NBRFarms_MagAll.png", ViolinMg_all, bg = "transparent", width = 18.5, height = 6, units = "cm")

# LOI
PlotLOI_all <- subset(soilchem_outreach, select = c(GLØDETAP)) %>% 
  gather(key = "Variables", value = "Rates")
ViolinLOI_all <- ggplot(PlotLOI_all, aes(x=Variables, y=Rates, fill=Variables)) +
  geom_violin(width=1, size=0.2) +
  scale_fill_grey() +
  geom_hline(yintercept=mean(PlotLOI_all$Rates), color="chartreuse3", size=5) +
  #geom_hline(yintercept=median(PlotLOI_all$Rates), color="cornflowerblue", size=5) +
  #scale_color_viridis(discrete=TRUE) +
  theme(
    legend.position="none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  coord_flip()
ViolinLOI_all
ggsave("outreach/NBRFarms_LOIAll.png", ViolinLOI_all, bg = "transparent", width = 18.5, height = 6, units = "cm")


# Outreach beetles
Map_Art <- Map_Art |> 
  mutate(PercentCarab = Carabidae/Beetle*100) |> 
  mutate(PercentStaph = Staphylinidae/Beetle*100) |> 
  mutate(PercentHydro = Hydrophilidae/Beetle*100) |> 
  mutate(PercentPtili = Ptiliidae/Beetle*100) |> 
  mutate(PercentScara = Scarabaeidae/Beetle*100) |> 
  mutate(PercentCurcu = Curculionidae/Beetle*100) |> 
  mutate(PercentElat = Elateridae/Beetle*100) |> 
  mutate(PercentSilph = Silphidae/Beetle*100) |> 
  mutate(PercentHist = Histeridae/Beetle*100) |> 
  mutate(percentGeot = Geotrupidae/Beetle*100)
write.xlsx(Map_Art, "OutreachBeetle.xlsx")

