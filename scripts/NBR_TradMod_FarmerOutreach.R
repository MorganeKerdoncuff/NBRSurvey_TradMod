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
beetle_full <- read.csv("data/cleandata/NBR_FullBeetleComm.csv", sep=",") # Clean arthropod community data


#### DATA PREPARATION ####

#
## Variables extraction

siteinfo_outreach <- subset(siteinfo_full, select = c(SiteID, Type_livestock, Habitat, EPSG.25832_X, EPSG.25832_Y))
area20x20_outreach <- subset(area20x20_full, select = c(SiteID, Elevation_max))
VGrichness_outreach <- subset(groundcover_full, select = c(SiteID, Plant_species_richness))
soilbulk_outreach <- subset(soilbulk_full, select = c(SiteID, BD))
soilchem_outreach <- subset(soilchem_full, select = c(SiteID, LOI, pH, P.Al_mg.100g, K.Al_mg.100g, Mg.Al_mg.100g, Ca.Al_mg.100g, Na.Al_mg.100g, TotalN_percentDM))
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
  summarise(LOI = mean(LOI), pH = mean(pH), P.Al_mg.100g = mean(P.Al_mg.100g), K.Al_mg.100g = mean(K.Al_mg.100g), Mg.Al_mg.100g = mean(Mg.Al_mg.100g), Ca.Al_mg.100g = mean(Ca.Al_mg.100g), Na.Al_mg.100g = mean(Na.Al_mg.100g), TotalN_percentDM = mean(TotalN_percentDM)) #missing OV1 (Oygarden) and UC1

# Plant abundance cover
vege_outreach <- vege_outreach |> 
  group_by(SiteID, Species) |> 
  summarise(PlantSp_cover = mean(Abundance))

# Beetle abundance
arthro_outreach <- arthro_outreach %>% 
  group_by(SiteID, BeetleFamilies) %>% 
  summarise(BeetleFam_abundance = sum(BeetleFam_abundance, na.rm = TRUE))

#  
## Summary dataset for GIS mapping

Map_All <- purrr::reduce(list(siteinfo_outreach, area20x20_outreach, VGrichness_outreach, soilbulk_outreach, soilchem_outreach), dplyr::left_join)

#Replace NA by zeros -> needed for QGIS, otherwise it treats the layer as character and does not want to apply a graduated aesthetics
Map_All[is.na(Map_All)] <- 0

#csv for QGIS map -> should write a new script !
write.csv(Map_All, file = "NBRenv_AttributeTable.csv")


#### Plant composition ####

vege_outreach <- vege_full |> 
  group_by(SiteID, Species) |> 
  summarise(PlantSp_cover = mean(Abundance)) |> 
  filter(SiteID == "IC3") |>
  filter(PlantSp_cover>0) |> 
  arrange(desc(PlantSp_cover))
summary(vege_outreach)

#### Beetle composition ####

beetle_outreach <- beetle_full |> 
  group_by(SiteID, BeetleFamilies) |> 
  summarise(BeetleFam_abundance = sum(BeetleFam_abundance, na.rm = TRUE)) |> 
  filter(SiteID == "IC1") |> 
  arrange(desc(BeetleFam_abundance))
#head(beetle_outreach)
sum(beetle_outreach$BeetleFam_abundance)

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
  geom_hline(yintercept=0.55, color="chartreuse3", size=5) +
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
  geom_hline(yintercept=5.1, color="chartreuse3", size=5) +
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
  geom_hline(yintercept=2.75, color="chartreuse3", size=5) +
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
  geom_hline(yintercept=15.5, color="chartreuse3", size=5) +
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
  geom_hline(yintercept=32.05, color="chartreuse3", size=5) +
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

