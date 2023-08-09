# TradMod WP2 - NBR survey grazer loading script
#Description of the data
## Biodiversity - Plant species richness per functional group
## Productivity - Plant aboveground biomass
## Carbon storage - LOI
## Nutrient cycling - abundance Acari & Collembola
#Year - 2023
#Who - Morgane Kerdoncuff
#Project - TradMod
#Funding - NFR
#Place - University of Bergen, Norway


#### PACKAGE LOADING ####

library(tidyverse) # R language
library(purrr) # Data manipulation: function "reduce" to bind several tables at the same time
library(vegan) # Community ecology analysis
library(dendextend) # Visual dendrogram
library(cluster) # Fuzzy clustering
library(ggplot2) # Visual representation
# library(ggord) # Ordination plot with ggplot2
# library(ggpubr) # Function ggarrange for several plots on same file
# library(GGally) # Extension ggplot


#### DATA LOADING ####

siteinfo_full <- read.csv("data/cleandata/NBR_FullSiteInfo.csv", sep=",") # Clean site info data
plantrichness_full <- read.csv("data/cleandata/NBR_FullPlantComm.csv", sep=",") # Clean aboveground cover data
soilchem_full <- read.csv("data/cleandata/NBR_FullSoilChem.csv", sep=",") # Clean soil chemistry data
mesobio_full <- read.csv("data/cleandata/NBR_FullMesobio.csv", sep=",") # Clean soil mesofauna data
biomass_full <- read.csv("data/cleandata/NBR_FullBiomass.csv", sep=",") # Clean plant biomass data



#### DATA PREPARATION ####

#
## Filter grassland sites

# Site selection
siteinfo_sheep <- siteinfo_full |>  
  filter(Type_livestock == "sheep") |> 
  dplyr::select(SiteID)

# Extraction in other datasets
plantrichness_sheep <- filter(plantrichness_full, SiteID %in% siteinfo_sheep$SiteID)
soilchem_sheep <- filter(soilchem_full, SiteID %in% siteinfo_sheep$SiteID)
mesobio_sheep <- filter(mesobio_full, SiteID %in% siteinfo_sheep$SiteID)
biomass_sheep <- filter(biomass_full, SiteID %in% siteinfo_sheep$SiteID)


#
## Plant species richness according to functional group

# Summarise at site level, by average
plantrichness_sheep <- plantrichness_sheep |> 
  group_by(SiteID, Species) |> 
  summarise(Abundance = mean(Abundance)) |> 
  filter(Abundance > 0) |> 
  arrange(Species)

# List species
#unique(plantrichness_sheep$Species)

# New functional group variable
plantrichness_sheep <- plantrichness_sheep |>
  # assign functional group to species
  mutate(FunctionalGroup = case_when(
    grepl("Achillea|Alchemilla|Anemone|Bartsia|Campanula|Cardamine|Cerastium|Cirsium|Conopodium|Dactylorhiza|Digitalis|Epilobium|Fraxinus|Galium|Gnaphalium|Hieracium|Hypericum|Hypochaeris|Lathyrus|Leontodon|Lotus|Melampyrum|Moneses|Narthecium|Oxalis|Pedicularis|Pinguicula|Plantago|Potentilla|Ranunculus|Rumex|Sagina|Sedum|Solidago|Stellaria|Succisa|Taraxacum|Trientalis|Trifolium|Veronica|Viola", Species) ~ "SPforbs",
    grepl("Agrostis|Aira|Alopecurus|Anthoxanthum|Carex|Dactylis|Danthonia|Deschampsia|Eriophorum|Festuca|Holcus|Juncus|Lolium|Luzula|Molinia|Nardus|Poa|Trichophorum", Species) ~ "SPmonocotyledons",
    grepl("Andromeda|Arctostaphylos|Betula|Calluna|Chamaepericlymenum|Empetrum|Erica|Juniperus|Loiseleuria|Picea|Polygala|Polygonum|Populus|Prunus|Rubus|Salix|Sorbus|Ulmus|Vaccinium", Species) ~ "SPwoody",
    grepl("Athyrium|Blechnum|Phegopteris|Polypodium|Pteridium", Species) ~ "SPferns",
    .default = "SPbryophytes"
  )) |> 
  # number of species per functional group
  group_by(SiteID, FunctionalGroup) |> 
  summarise(S = n()) |> 
  # make functional groups as variables
  pivot_wider(names_from = FunctionalGroup, values_from = S) |> 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
#
## Summarise data at site level -> should be 23 observations for non community data

# Site info - validated

# Plant species richness
## From groundcover data
## Same number of replicates per site
## Currently at sample level -> summary by average
# groundcover_sheep <- groundcover_sheep |>  
#   group_by(SiteID) |>  
#   summarise(MeanRichness = mean(Plant_species_richness, na.rm=TRUE))

## !! need total site species richness, from community data

# Loss of ignition
## From chemistry data
## Same number of replicates per site
## Currently at plot level -> summary by average
soilchem_sheep <- soilchem_sheep |> 
  group_by(SiteID) |> 
  summarise(MeanLOI = mean(LOI))

# Mite and springtail abundance
## From mesofauna data
## Different numbers of replicates per site -> need 3 replicates, 1 per plot
## Currently at sample level -> summary by average
mesobio_sheep <- mesobio_sheep |> 
  group_by(SiteID, PlotID) |> 
  slice(1) |>
  # Select randomly one row which match unique combination of site & plot IDs
  ungroup() |> 
  group_by(SiteID) |> 
  summarise(Acari = mean(Acari),
            Collembola = mean(Collembola))

# Woody, grass and forb biomass
## From biomass data
## Different numbers of replicates per site -> need 3 replicates, 1 per plot
## Currently at sample level -> summary by average
biomass_sheep <- biomass_sheep |>
  group_by(SiteID, PlotID) |>
  pivot_wider(names_from = FunctionalType, values_from = DWbiomass_g) |> 
  # Select randomly one row which match unique combination of site & plot IDs
  slice(1) |>
  ungroup() |>
  # make functional groups as variables
  group_by(SiteID) |>
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
  summarise(woody = mean(woody),
            monocotyledons = mean(monocotyledons),
            forbs = mean(forbs))

#
## Initiate base dataset for hierachical clustering

data_sheep <- purrr::reduce(list(siteinfo_sheep, subset(plantrichness_sheep, select = -c(SPferns)), biomass_sheep, mesobio_sheep, soilchem_sheep), dplyr::left_join)


#### HIERARCHICAL CLUSTERING ####

#
## Clustering & dendrogram

# Extraction SiteID
siteID <- data_sheep$SiteID
# data_sheep <- subset(data_sheep, select = -c(SiteID))

# Site ID as rowname
rownames(data_sheep) <- data_sheep$SiteID

# Data scaling
data_sheep_sc <- as.data.frame(scale(subset(data_sheep, select = -c(SiteID))))

# Distance matrix - environmental continuous numeric -> Euclidean
dist_sheep <- dist(data_sheep_sc, method = "euclidean")

# Clustering using Ward's (1963) clustering criterion
hclustwd_sheep <- hclust(dist_sheep, method = "ward.D2")

# Plot HC tree
plot(hclustwd_sheep)

# Show clusters (k=5) with base
rect.hclust(hclustwd_sheep, k = 5, border = 2:6)

# Show clusters (k=5) with color_branches
dendro_sheep <- color_branches(as.dendrogram(hclustwd_sheep), 
                               k = 5, 
                               #groupLabels = TRUE,
                               col = c(2, 3, 4, 5, 6))
plot(dendro_sheep)

# Cluster group
clusterwd <- stats::cutree(hclustwd_sheep, k = 5)

#
## Ordination clusters

# NMDS -> one-dimension
nmds_sheep <- metaMDS(dist_sheep)
# ordiplot(nmds_sheep, type = 'n')
# points(nmds_sheep, pch = clusterwd+20, bg = clusterwd)
# legend('topright', legend = 1:5, pch = (1:5)+20, pt.bg = 1:5)

# Colour grouping
#clusterwd # see order
col_clusters <- c("green3", "purple", "cyan2", "blue1", "red2")
col_clusters[clusterwd]

# Plot ordination with base
ordiplot(nmds_sheep, type = 'n')
points(nmds_sheep, 
       col = col_clusters[clusterwd], 
       pch = clusterwd+20)
text(nmds_sheep, 
     col = col_clusters[clusterwd], 
     labels=rownames(data_sheep_sc),
     cex = 0.7,
     pos = 1)
legend("bottomright", 
       legend = paste("Cluster", 1:5),
       col = col_clusters, 
       pt.bg = col_clusters,
       bty = "n", 
       pch = (1:5)+20)
ordihull(nmds_sheep, 
         groups = clusterwd, 
         display = "sites")

#### CLUSTER CHARTS ####

#
## Data preparation

# Extraction cluster group
data_sheep <- data_sheep |>
  mutate(cluster = clusterwd)

# Make variable group - Productivity, Diversity, Carbon, Soil
data_sheep <- data_sheep |> 
  mutate(group = ifelse())