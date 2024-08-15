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
library(FSA) # Dunn's Test after KS
library(dunn.test) # Dunn's Test after KS -> "classic" but less handy to use, but to refer to
library(dendextend) # Visual dendrogram
library(pvclust) # Multiscale bootstrap resampling
library(ggplot2) # Visual representation
library(concaveman) # better hull
library(ggforce) # better ggplot
library(forcats) # reorder groups for plots
library(ggpubr) # Function ggarrange for several plots on same file
# library(GGally) # Extension ggplot


#### DATA LOADING ####

siteinfo_full <- read.csv("data/cleandata/NBR_FullSiteInfo.csv", sep=",") # Clean site info data
landuse_full <- read.csv("data/cleandata/NBR_FullLandUse.csv", sep=",") # Clean site info data
plantrichness_full <- read.csv("data/cleandata/NBR_FullPlantComm.csv", sep=",") # Clean aboveground cover data
soilchem_full <- read.csv("data/cleandata/NBR_FullSoilChem.csv", sep=",") # Clean soil chemistry data
mesobio_full <- read.csv("data/cleandata/NBR_FullMesobio.csv", sep=",") # Clean soil mesofauna data
biomass_full <- read.csv("data/cleandata/NBR_FullBiomass.csv", sep=",") # Clean plant biomass data
climate_full <- read.csv("data/cleandata/NBR_FullClimate.csv", sep=",") # Clean climate data
soilchem_full <- read.csv("data/cleandata/NBR_FullSoilChem.csv", sep=",") # Clean soil chemistry data


#### DATA PREPARATION ####

#
## Filter grassland sites

# Site selection
siteinfo_sheep <- siteinfo_full |>  
  filter(Livestock == "sheep") |> 
  dplyr::select(SiteID, EcoZone) 
  
# Extraction in other datasets
landuse_sheep <- filter(landuse_full, SiteID %in% siteinfo_sheep$SiteID)
plantrichness_sheep <- filter(plantrichness_full, SiteID %in% siteinfo_sheep$SiteID)
soilchem_sheep <- filter(soilchem_full, SiteID %in% siteinfo_sheep$SiteID)
mesobio_sheep <- filter(mesobio_full, SiteID %in% siteinfo_sheep$SiteID)
biomass_sheep <- filter(biomass_full, SiteID %in% siteinfo_sheep$SiteID)
climate_sheep <- filter(climate_full, SiteID %in% siteinfo_sheep$SiteID)
soilchem_sheep <- filter(soilchem_full, SiteID %in% siteinfo_sheep$SiteID)

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
totplantsp <- length(unique(plantrichness_sheep$Species))

# New functional group variable
plantrichness_sheep <- plantrichness_sheep |>
  # assign functional group to species
  mutate(FunctionalGroup = case_when(
    grepl("Achillea|Alchemilla|Anemone|Bartsia|Campanula|Cardamine|Cerastium|Cirsium|Conopodium|Dactylorhiza|Digitalis|Epilobium|Euphrasia|Fraxinus|Galium|Gnaphalium|Hieracium|Hypericum|Hypochaeris|Lathyrus|Leontodon|Lotus|Melampyrum|Moneses|Narthecium|Oxalis|Pedicularis|Pinguicula|Plantago|Potentilla|Ranunculus|Rumex|Sagina|Sedum|Solidago|Stellaria|Succisa|Taraxacum|Trientalis|Trifolium|Veronica|Viola", Species) ~ "div_forbs",
    grepl("Agrostis|Aira|Alopecurus|Anthoxanthum|Carex|Dactylis|Danthonia|Deschampsia|Eriophorum|Festuca|Holcus|Juncus|Lolium|Luzula|Molinia|Nardus|Poa|Trichophorum", Species) ~ "div_monocotyledons",
    grepl("Andromeda|Arctostaphylos|Betula|Calluna|Chamaepericlymenum|Empetrum|Erica|Juniperus|Loiseleuria|Picea|Polygala|Polygonum|Populus|Prunus|Rubus|Salix|Sorbus|Ulmus|Vaccinium", Species) ~ "div_woody",
    grepl("Athyrium|Blechnum|Phegopteris|Polypodium|Pteridium", Species) ~ "div_ferns",
    .default = "div_cryptogams"
  ))

# Total number species per functional group (writing)
totforbsp <- plantrichness_sheep |> 
  group_by(Species, FunctionalGroup) |>
  summarise(Abundance = mean(Abundance)) |> 
  filter(FunctionalGroup == "div_forbs")
totmonocsp <- plantrichness_sheep |> 
  group_by(Species, FunctionalGroup) |>
  summarise(Abundance = mean(Abundance)) |> 
  filter(FunctionalGroup == "div_monocotyledons")
totwoodysp <- plantrichness_sheep |> 
  group_by(Species, FunctionalGroup) |>
  summarise(Abundance = mean(Abundance)) |> 
  filter(FunctionalGroup == "div_woody")
totbryosp <- plantrichness_sheep |> 
  group_by(Species, FunctionalGroup) |>
  summarise(Abundance = mean(Abundance)) |> 
  filter(FunctionalGroup == "div_cryptogams")

# Number of species per functional group per site
plantrichness_sheep <- plantrichness_sheep |>   
  # number of species per functional group
  group_by(SiteID, FunctionalGroup) |> 
  summarise(S = n()) |> 
  # make functional groups as variables
  pivot_wider(names_from = FunctionalGroup, values_from = S) |> 
  mutate_if(is.numeric, ~replace(., is.na(.), 0))
  
#
## Summarise data at site level -> should be 23 observations for non community data

# Site info - ok

# Loss of ignition
## From chemistry data
## Same number of replicates per site
## Currently at plot level -> summary by average
soilchem_sheep <- soilchem_sheep |> 
  group_by(SiteID) |> 
  summarise(LOI = mean(LOI))

# Mite and springtail abundance
## From mesofauna data
## Different numbers of replicates per site -> need 3 replicates, 1 per plot
## Currently at sample level -> summary by average
mesobio_sheep <- mesobio_sheep |> 
  group_by(SiteID, PlotID) |>
  # Select randomly one row which match unique combination of site & plot IDs
  slice(1) |>
  ungroup()

# Total acari and collembola collected
totacari <- sum(mesobio_sheep$Acari)
totcollem <- sum(mesobio_sheep$Collembola)

# Summary by average
mesobio_sheep <- mesobio_sheep |> 
  group_by(SiteID) |> 
  summarise(ab_acari = mean(Acari),
            ab_collembola = mean(Collembola))

# Woody, grass, forb and bryo biomass
## From biomass data
## Currently at sample level -> summary by average
biomass_sheep <- biomass_sheep |>
  pivot_wider(names_from = FunctionalType, values_from = DWbiomass_g) |> 
  # make functional groups as variables
  group_by(SiteID) |>
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
  summarise(biom_woody = mean(woody),
            biom_monocotyledons = mean(monocotyledons),
            biom_forbs = mean(forbs),
            biom_cryptogams = mean(cryptogams))

#
## Initial base dataset

# Merge all data
data_sheep <- purrr::reduce(list(
  siteinfo_sheep,
  subset(landuse_sheep, select = c(SiteID, FieldType)),
  subset(plantrichness_sheep, select = -c(div_ferns)), 
  biomass_sheep, 
  mesobio_sheep, 
  soilchem_sheep), dplyr::left_join)

# Intuitive ID according to following code - COx (costal outfield), CIx (coastal infield), FIx (fjord infield), MOx (mountain oufield)
data_sheep <- data_sheep |> 
  mutate(fieldclass = ifelse(
    EcoZone == "mountain" & FieldType == "outfield", "MO", ifelse(
      EcoZone == "coastal" & FieldType == "infield", "CI", ifelse(
        EcoZone == "fjord" & FieldType == "infield", "FI", ifelse(
          EcoZone == "fjord" & FieldType == "outfield", "FO", "CO"
        )
      )
    )
  )) |> 
  group_by(fieldclass) |> 
  mutate(FieldID = paste(fieldclass, row_number(), sep = ""))

# # Intuitive ID according to following code - COx (costal outfield), Ix (infield), MOx (mountain oufield)
# data_sheep <- data_sheep |>
#   mutate(fieldclass_b = case_when(
#     grepl("U", SiteID) ~ "MO",
#     grepl("IS2|OV|OS5|OS7|OS9", SiteID) ~ "CO",
#     .default = "IN"
#   )) |>
#   group_by(fieldclass_b) |>
#   mutate(FieldID = paste(fieldclass_b, row_number(), sep = ""))

# Site ID as rowname
# rownames(data_sheep) <- data_sheep$FieldID



#### HIERARCHICAL CLUSTERING ####

# Data scaling
data_sheep_sc <- as.data.frame(scale(subset(data_sheep, select = -c(SiteID, EcoZone, FieldType, fieldclass, FieldID))))
rownames(data_sheep_sc) <- data_sheep$FieldID
#rownames(data_sheep_sc) <- data_sheep$SiteID

#
## Clustering & dendrogram

# Distance matrix - environmental continuous numeric -> Euclidean
dist_sheep <- dist(data_sheep_sc, method = "euclidean")

# Clustering using Ward's (1963) clustering criterion
hclustwd_sheep <- hclust(dist_sheep, method = "ward.D2")

# Plot HC tree
plot(hclustwd_sheep)

# Show primary (k=2) & secondary (k=4) clusters with base
rect.hclust(hclustwd_sheep, k = 4, border = 2:6)

# Show primary (k=2) & secondary (k=4) clusters with color_branches
dendro_sheep <- color_branches(as.dendrogram(hclustwd_sheep), 
                               k = 4, 
                               #groupLabels = TRUE,
                               col = c("orange", "orange3", "pink", "pink4"))
                               # col = c("pink", "pink4", "orange", "orange3"))
plot(dendro_sheep)

# Primary and secondary cluster groups
clusterwd_prim <- stats::cutree(hclustwd_sheep, k = 2)
clusterwd_sec <- stats::cutree(hclustwd_sheep, k = 4)

# Save dendrogram
png(filename = "outputs/cluster_tree.png", width = 14, height = 8, units = "cm", res = 600, pointsize = 9, bg = "transparent")
plot(dendro_sheep)
dev.off()

#
## Bootstrap method with pvclust

# Contingency table
# data_sheep_long <- pivot_longer(data_sheep, cols = -c(SiteID), names_to = "variables", values_to = "val")
# data_sheep_contin <- xtabs(formula = val ~ variables + SiteID, data = data_sheep_long)
# 
# # Data scaling
# data_sc <- scale(data_sheep_contin)
# 
# # Calculate pval of clusters
# cluster_pval <- pvclust(data = data_sc, method.hclust = "ward.D2", method.dist = 'euclidean')
# 
# # Plot dendrogram
# plot(cluster_pval)

#
## Comparison theory distribution (classification) vs observed (primary cluster)


#
## Ordination clusters

# NMDS -> one-dimension
nmds_sheep <- metaMDS(dist_sheep)
# ordiplot(nmds_sheep, type = 'n')
# points(nmds_sheep, pch = clusterwd+20, bg = clusterwd)
# legend('topright', legend = 1:5, pch = (1:5)+20, pt.bg = 1:5)

# Colour grouping
#clusterwd # see order
# col_clusters <- c("pink", "orange4", "orange", "pink4")
# col_clusters[clusterwd_sec]
# 
# # Plot ordination with base
# cluster_ord <- ordiplot(nmds_sheep, type = 'n')
# cluster_ord <- points(nmds_sheep, 
#        col = col_clusters[clusterwd_sec], 
#        pch = clusterwd_sec+20)
# cluster_ord <- text(nmds_sheep, 
#      col = col_clusters[clusterwd_sec], 
#      labels=rownames(data_sheep_sc),
#      cex = 0.7,
#      pos = 1)
# cluster_ord <- legend( 
#        legend = paste("Cluster", 1:4),
#        x = 1.2,
#        y = 3,
#        #horiz = TRUE,
#        text.width = 1,
#        x.intersp = 1,
#        y.intersp = 0.3,
#        col = col_clusters, 
#        pt.bg = col_clusters,
#        bty = "n", 
#        pch = (1:4)+20)
# cluster_ord <- ordihull(nmds_sheep, 
#          groups = clusterwd_sec, 
#          display = "sites")
# cluster_ord

#
## Ordination cluster ggplot

# Extraction scores
nmds_scores <- as.data.frame(scores(nmds_sheep))
nmds_scores$FieldID <- rownames(nmds_scores)

# Add groups
nmds_scores <- left_join(nmds_scores, subset(data_sheep, select = c(FieldType, EcoZone, FieldID)))
nmds_scores <- nmds_scores |> 
  mutate(clusterID = paste("CL", clusterwd_sec, sep = ""))

# Plotting
ggcluster_ord <- ggplot(data = nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(
    shape = EcoZone,
    colour = FieldType
  ), size = 3) +
  scale_colour_manual(values = c("infield" = "gray", "outfield" = "black")) +
  geom_mark_ellipse(aes(
    fill = clusterID
  )) +
  scale_fill_manual(values = c("CL1" = "pink", "CL2" = "orange4", "CL3" = "orange", "CL4" = "pink4")) +
  # geom_mark_hull(aes(
  #   fill = clustercol
  # )) +
  # geom_text(data = nmds_scores, aes(
  #   x = NMDS1,
  #   y = NMDS2,
  #   label = FieldID
  # ))
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        panel.background = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
ggcluster_ord
ggsave(filename = "outputs/cluster_NMDS.png", plot = ggcluster_ord, width = 18, height = 10, units = "cm")

#### ANALYSIS OF CLUSTERS ####

# Extraction primary & secondary cluster group
data_sheep_sc <- data_sheep_sc |>
  mutate(cluster_prim = as.factor(clusterwd_prim)) |> 
  mutate(cluster_sec = as.factor(clusterwd_sec))

#
## Comparison primary clusters - unpaired non-parametric Mann-Whitney test

# Empty table for summary statistics
wcxprim_stat <- data.frame(var = as.character(),
                       W = as.numeric(),
                       pval = as.numeric(),
                       stringsAsFactors = FALSE)

# LOI
wcx <- wilcox.test(LOI ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "LOI", W = wcx$statistic, pval = wcx$p.value)

# Abundance collembola
wcx <- wilcox.test(ab_collembola ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "ab_collembola", W = wcx$statistic, pval = wcx$p.value)

# Abundance acari
wcx <- wilcox.test(ab_acari ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "ab_acari", W = wcx$statistic, pval = wcx$p.value)

# Biomass cryptogams
wcx <- wilcox.test(biom_cryptogams ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "biom_cryptogams", W = wcx$statistic, pval = wcx$p.value)

# Biomass monocotyledons
wcx <- wilcox.test(biom_monocotyledons ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "biom_monocotyledons", W = wcx$statistic, pval = wcx$p.value)

# Biomass forbs
wcx <- wilcox.test(biom_forbs ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "biom_forbs", W = wcx$statistic, pval = wcx$p.value)

# Biomass woody
wcx <- wilcox.test(biom_woody ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "biom_woody", W = wcx$statistic, pval = wcx$p.value)

# Diversity cryptogams
wcx <- wilcox.test(div_cryptogams ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "div_cryptogams", W = wcx$statistic, pval = wcx$p.value)

# Diversity monocotyledons
wcx <- wilcox.test(div_monocotyledons ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "div_monocotyledons", W = wcx$statistic, pval = wcx$p.value)

# Diversity forbs
wcx <- wilcox.test(div_forbs ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "div_forbs", W = wcx$statistic, pval = wcx$p.value)

# Diversity woody
wcx <- wilcox.test(div_woody ~ cluster_prim, data = data_sheep_sc)
wcxprim_stat <- wcxprim_stat |> 
  add_row(var = "div_woody", W = wcx$statistic, pval = wcx$p.value)

#
## Comparison secondary clusters infield group (n=2) - unpaired non-parametric Mann-Whitney test

# Empty table for summary statistics
wcxsecB_stat <- data.frame(var = as.character(),
                       W = as.numeric(),
                       pval = as.numeric(),
                       stringsAsFactors = FALSE)

# LOI
wcx <- wilcox.test(LOI ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "LOI", W = wcx$statistic, pval = wcx$p.value)

# Abundance collembola
wcx <- wilcox.test(ab_collembola ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "ab_collembola", W = wcx$statistic, pval = wcx$p.value)

# Abundance acari
wcx <- wilcox.test(ab_acari ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "ab_acari", W = wcx$statistic, pval = wcx$p.value)

# Biomass cryptogams
wcx <- wilcox.test(biom_cryptogams ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "biom_cryptogams", W = wcx$statistic, pval = wcx$p.value)

# Biomass monocotyledons
wcx <- wilcox.test(biom_monocotyledons ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "biom_monocotyledons", W = wcx$statistic, pval = wcx$p.value)

# Biomass forbs
wcx <- wilcox.test(biom_forbs ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "biom_forbs", W = wcx$statistic, pval = wcx$p.value)

# Biomass woody
wcx <- wilcox.test(biom_woody ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "biom_woody", W = wcx$statistic, pval = wcx$p.value)

# Diversity cryptogams
wcx <- wilcox.test(div_cryptogams ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "div_cryptogams", W = wcx$statistic, pval = wcx$p.value)

# Diversity monocotyledons
wcx <- wilcox.test(div_monocotyledons ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "div_monocotyledons", W = wcx$statistic, pval = wcx$p.value)

# Diversity forbs
wcx <- wilcox.test(div_forbs ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "div_forbs", W = wcx$statistic, pval = wcx$p.value)

# Diversity woody
wcx <- wilcox.test(div_woody ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3))
wcxsecB_stat <- wcxsecB_stat |> 
  add_row(var = "div_woody", W = wcx$statistic, pval = wcx$p.value)

#
## Comparison secondary clusters outfield group (n=2) - unpaired non-parametric Mann-Whitney test

# Empty table for summary statistics
wcxsec_statA <- data.frame(var = as.character(),
                          W = as.numeric(),
                          pval = as.numeric(),
                          stringsAsFactors = FALSE)

# LOI
wcx <- wilcox.test(LOI ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "LOI", W = wcx$statistic, pval = wcx$p.value)

# Abundance collembola
wcx <- wilcox.test(ab_collembola ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "ab_collembola", W = wcx$statistic, pval = wcx$p.value)

# Abundance acari
wcx <- wilcox.test(ab_acari ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "ab_acari", W = wcx$statistic, pval = wcx$p.value)

# Biomass cryptogams
wcx <- wilcox.test(biom_cryptogams ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "biom_cryptogams", W = wcx$statistic, pval = wcx$p.value)

# Biomass monocotyledons
wcx <- wilcox.test(biom_monocotyledons ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "biom_monocotyledons", W = wcx$statistic, pval = wcx$p.value)

# Biomass forbs
wcx <- wilcox.test(biom_forbs ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "biom_forbs", W = wcx$statistic, pval = wcx$p.value)

# Biomass woody
wcx <- wilcox.test(biom_woody ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "biom_woody", W = wcx$statistic, pval = wcx$p.value)

# Diversity cryptogams
wcx <- wilcox.test(div_cryptogams ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "div_cryptogams", W = wcx$statistic, pval = wcx$p.value)

# Diversity monocotyledons
wcx <- wilcox.test(div_monocotyledons ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "div_monocotyledons", W = wcx$statistic, pval = wcx$p.value)

# Diversity forbs
wcx <- wilcox.test(div_forbs ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "div_forbs", W = wcx$statistic, pval = wcx$p.value)

# Diversity woody
wcx <- wilcox.test(div_woody ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4))
wcxsec_statA <- wcxsec_statA |> 
  add_row(var = "div_woody", W = wcx$statistic, pval = wcx$p.value)


#
## Comparison secondary clusters outfield group (n=3) - Kruskal-Wallis followed by Dunn's test

# Empty table for summary stat
# kssec_stat <- data.frame(var = as.character(),
#                         KWchi = as.numeric(),
#                         df = as.numeric(),
#                         pval = as.numeric(),
#                         stringsAsFactors = FALSE)
# 
# # LOI
# ks <- kruskal.test(LOI ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "LOI", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(LOI ~ cluster_sec, 
#           data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#           method = "bonferroni")
# dunnsig_LOI <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "LOI")
# 
# # Abundance collembola
# ks <- kruskal.test(ab_collembola ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "ab_collembola", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(ab_collembola ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_collembola <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "ab_collembola")
# 
# # Abundance acari
# ks <- kruskal.test(ab_acari ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "ab_acari", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(ab_acari ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_acari <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "ab_acari")
# 
# # Biomass cryptogams
# ks <- kruskal.test(biom_cryptogams ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "biom_cryptogams", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(biom_cryptogams ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_biom.cryptogams <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "biom_cryptogams")
# 
# # Biomass monocotyledons
# ks <- kruskal.test(biom_monocotyledons ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "biom_monocotyledons", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(biom_monocotyledons ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_biom.monocotyledons <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "biom_monocotyledons")
# 
# # Biomass forbs
# ks <- kruskal.test(biom_forbs ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "biom_forbs", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(biom_forbs ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_biom.forbs <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "biom_forbs")
# 
# # Biomass woody
# ks <- kruskal.test(biom_woody ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "biom_woody", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(biom_woody ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_biom.woody <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "biom_woody")
# 
# # Diversity cryptogams
# ks <- kruskal.test(div_cryptogams ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "div_cryptogams", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(div_cryptogams ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_div.cryptogams <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "div_cryptogams")
# 
# # Diversity monocotyledons
# ks <- kruskal.test(div_monocotyledons ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "div_monocotyledons", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(div_monocotyledons ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_div.monocotyledons <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "div_monocotyledons")
# 
# # Diversity forbs
# ks <- kruskal.test(div_forbs ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "div_forbs", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(div_forbs ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_div.forbs <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "div_forbs")
# 
# # Diversity woody
# ks <- kruskal.test(div_woody ~ cluster_sec, data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5))
# kssec_stat <- kssec_stat |> 
#   add_row(var = "div_woody", KWchi = ks$statistic, df = ks$parameter, pval = ks$p.value)
# dunn <- dunnTest(div_woody ~ cluster_sec, 
#                  data = filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4 | cluster_sec == 5),
#                  method = "bonferroni")
# dunnsig_div.woody <- filter(dunn$res, dunn$res$P.adj < 0.1) |> 
#   mutate(Var = "div_woody")

# Stat summary
# dunsig_summary <- purrr::reduce(list(dunnsig_LOI, dunnsig_acari, dunnsig_collembola, dunnsig_biom.cryptogams, dunnsig_biom.forbs, dunnsig_biom.monocotyledons, dunnsig_biom.woody, dunnsig_div.cryptogams, dunnsig_div.forbs, dunnsig_div.monocotyledons, dunnsig_div.woody), dplyr::full_join)

#
## Visual representation

# Make variable group - Productivity, Diversity, Carbon, Soil nutrients
data_sheep_plot <- data_sheep_sc |> 
  pivot_longer(cols = -c(cluster_prim, cluster_sec),
               names_to = "variables",
               values_to = "sc_val") |> 
  mutate(var_group = ifelse(
    grepl("div_", variables), "diversity", ifelse(
      variables == "ab_acari" | variables == "ab_collembola" | variables == "LOI", "decomposition", "production"
    )
  )) |> 
  # rearrange by variable group for ggplot
  arrange(sc_val) |>
  mutate(variables = factor(variables, levels = c("div_woody", "div_forbs", "div_monocotyledons", "div_cryptogams", "biom_woody", "biom_forbs", "biom_monocotyledons", "biom_cryptogams", "ab_acari", "ab_collembola", "LOI")))

# Make cluster variable as factor
data_sheep_plot <- data_sheep_plot |> 
  mutate(cluster_prim = as.factor(cluster_prim)) |> 
  mutate(cluster_sec = as.factor(cluster_sec))

#
## Primary clusters

# x-axis plot arrange for primary clusters
clusterplot_axis <- filter(data_sheep_plot, cluster_prim == 1) |> 
  ggplot(aes(x = variables)) +
  xlab("") +
  ylab("") +
  labs("") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.title = element_blank()
  ) +
  scale_x_discrete(limits = c("div_woody", "div_forbs", "div_monocotyledons", "div_cryptogams", "biom_woody", "biom_forbs", "biom_monocotyledons", "biom_cryptogams", "ab_acari", "ab_collembola", "LOI")) +
  annotate("text", x = "LOI", y = 0, label = "***", size = 2.5) +
  annotate("text", x = "ab_collembola", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "ab_acari", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_cryptogams", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_monocotyledons", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_forbs", y = 0, label = "**", size = 2.5) +
  annotate("text", x = "biom_woody", y = 0, label = "***", size = 2.5) +
  annotate("text", x = "div_cryptogams", y = 0, label = "***", size = 2.5) +
  annotate("text", x = "div_monocotyledons", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "div_forbs", y = 0, label = "**", size = 2.5) +
  annotate("text", x = "div_woody", y = 0, label = "***", size = 2.5)
clusterplot_axis

# Primary cluster A
clusterplot_A <- filter(data_sheep_plot, cluster_prim == 1) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster A") +
  labs("") +
  ylim(-4, 4) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(colour = "pink3", size = 8),
    legend.title = element_blank()
  ) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_A

# Primary cluster B
clusterplot_B <- filter(data_sheep_plot, cluster_prim == 2) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster B") +
  labs("") +
  ylim(-4, 4) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(colour = "orange2", size = 8),
    legend.title = element_blank()
  ) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_B

# Arrange primary plots
clusterplot_prim <- ggarrange(clusterplot_axis, clusterplot_A, clusterplot_B, nrow = 1, common.legend = TRUE)
clusterplot_prim
ggsave(filename = "outputs/clusterplot_prim.png", plot = clusterplot_prim, width = 16, height = 9, units = "cm")

#
## Secondary clusters - Group A

# Customised x-axis
clusterplot_axis <- filter(data_sheep_plot, cluster_prim == 1) |> 
  ggplot(aes(x = variables)) +
  xlab("") +
  ylab("") +
  labs("") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    legend.title = element_blank()
  ) +
  scale_x_discrete(limits = c("div_woody", "div_forbs", "div_monocotyledons", "div_cryptogams", "biom_woody", "biom_forbs", "biom_monocotyledons", "biom_cryptogams", "ab_acari", "ab_collembola", "LOI")) +
  annotate("text", x = "LOI", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "ab_collembola", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "ab_acari", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_cryptogams", y = 0, label = "*", size = 2.5) +
  annotate("text", x = "biom_monocotyledons", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_forbs", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_woody", y = 0, label = "*", size = 2.5) +
  annotate("text", x = "div_cryptogams", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "div_monocotyledons", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "div_forbs", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "div_woody", y = 0, label = "°", size = 2.5)
clusterplot_axis

# Secondary cluster 1
clusterplot_1 <- filter(data_sheep_plot, cluster_sec == 1) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster 1") +
  labs("") +
  ylim(-4, 4) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(colour = "pink", size = 8),
    legend.title = element_blank(),
    legend.position = "none"
    ) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_1

# Secondary Cluster 4
clusterplot_4 <- filter(data_sheep_plot, cluster_sec == 4) |>
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster 4") +
  ylim(-4,4) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(colour = "pink4", size = 8),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_4

# # Secondary Cluster 5
# clusterplot_5 <- filter(data_sheep_plot, cluster_sec == 5) |> 
#   ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_boxplot() +
#   xlab("") +
#   ylab("Cluster 5") +
#   labs("") +
#   ylim(-4, 4) +
#   coord_flip() +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text.y = element_blank(),
#     axis.text.x = element_text(size = 7),
#     axis.title.x = element_text(colour = "pink2", size = 8),
#     legend.title = element_blank(),
#     legend.position = "none"
#   ) +
#   scale_fill_grey(start = 0.5, end = 1)
# clusterplot_5

# Arrange plots
clusterplot_secA <- ggarrange(clusterplot_axis, clusterplot_1, clusterplot_4, nrow = 1)
clusterplot_secA
ggsave(filename = "outputs/clusterplot_secA.png", plot = clusterplot_secA, width = 16, height = 8, units = "cm")

#
## Secondary clusters - Group B

# Customised x-axis
clusterplot_axis <- filter(data_sheep_plot, cluster_prim == 1) |> 
  ggplot(aes(x = variables)) +
  xlab("") +
  ylab("") +
  labs("") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.title = element_blank()
  ) +
  scale_x_discrete(limits = c("div_woody", "div_forbs", "div_monocotyledons", "div_cryptogams", "biom_woody", "biom_forbs", "biom_monocotyledons", "biom_cryptogams", "ab_acari", "ab_collembola", "LOI")) +
  annotate("text", x = "LOI", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "ab_collembola", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "ab_acari", y = 0, label = "*", size = 2.5) +
  annotate("text", x = "biom_cryptogams", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_monocotyledons", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_forbs", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "biom_woody", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "div_cryptogams", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "div_monocotyledons", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "div_forbs", y = 0, label = "NS", size = 2.5) +
  annotate("text", x = "div_woody", y = 0, label = "NS", size = 2.5)
clusterplot_axis


# Secondary Cluster 2
clusterplot_2 <- filter(data_sheep_plot, cluster_sec == 2) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster 2") +
  ylim(-4, 4) +
  labs("") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(colour = "orange4", size = 8),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_2

# Secondary Cluster 3
clusterplot_3 <- filter(data_sheep_plot, cluster_sec == 3) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster 3") +
  labs("") +
  ylim(-4, 4) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(colour = "orange", size = 8),
    legend.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_3

# Arrange plots
clusterplot_secB <- ggarrange(clusterplot_axis, clusterplot_2, clusterplot_3, nrow = 1)
clusterplot_secB
ggsave(filename = "outputs/clusterplot_secB.png", plot = clusterplot_secB, width = 16, height = 8, units = "cm")

# All plots
clusterplot_all <- ggarrange(clusterplot_prim, clusterplot_secA, clusterplot_secB, labels = c("A", "B", "C"), ncol = 1)
ggsave(filename = "outputs/clusterplot_all.png", 
       plot = clusterplot_all,
       bg = "white",
       width = 16, height = 24, units = "cm")


#### DATA SUMMARY ####

# Extraction column with SiteID
data_sheep_sc <- data_sheep_sc |>
  mutate(FieldID = row.names(data_sheep_sc))

# Join data
sheep_summary <- full_join(data_sheep, subset(data_sheep_sc, select = c(FieldID, cluster_prim, cluster_sec)))

# Summary primary comparison -> need to fix the "across col" to exclude non-numeric
sheep_summary_prim <- sheep_summary |>
  dplyr::group_by(cluster_prim) |>
  select_if(is.numeric) |>
  summarise(across(.cols = everything(), list(mean = mean, sd = sd)))

# Summary secondary comparison -> need to fix the "across col" to exclude non-numeric
sheep_summary_sec <- sheep_summary |>
  dplyr::group_by(cluster_sec) |>
  select_if(is.numeric) |>
  summarise(across(.cols = everything(), list(mean = mean, sd = sd)))


# #### LM AGAINST TEMP, PRECI & PHOSPHORUS ####
# 
# # New table for LM
# data_sheep_lm <- subset(data_sheep, select = -c(fieldtype, SiteID))
# 
# # Distribution explanatory variables
# hist(data_sheep_lm$annualprecipitation)
# hist(data_sheep_lm$avgtempJuly)
# hist(data_sheep_lm$phosphorus) # one outlier -> to be removed
# 
# #
# ## Correlation explanatory variables
# 
# # Precipitation x temperature
# lm <- lm(annualprecipitation ~ avgtempJuly, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$avgtempJuly, data_sheep_lm$annualprecipitation, pch = 16, col = "red")
# abline(lm)
# 
# # Precipitation x phosphorus
# lm <- lm(phosphorus ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$phosphorus, pch = 16, col = "red")
# abline(lm)
# 
# # Remove phosphorus outlier & temperature (highly correlated with precipitation)
# data_sheep_lm <- filter(data_sheep_lm, phosphorus < 20)
# data_sheep_lm <- subset(data_sheep_lm, select = -c(avgtempJuly))
# 
# #
# ## LOI
# 
# # Against precipitation simple
# lm <- lm(LOI ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$LOI, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(LOI ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$LOI, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Abundance collembola
# 
# # Against precipitation simple
# lm <- lm(ab_collembola ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$ab_collembola, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(ab_collembola ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$ab_collembola, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Abundance acari
# 
# # Against precipitation simple
# lm <- lm(ab_acari ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$ab_acari, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(ab_acari ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$ab_acari, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# # Summary plot lm soil
# 
# 
# #
# ## Biomass cryptogams
# 
# # Against precipitation simple
# lm <- lm(biom_cryptogams ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$biom_cryptogams, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(biom_cryptogams ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$biom_cryptogams, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Biomass monocotyledons
# 
# # Against precipitation simple
# lm <- lm(biom_monocotyledons ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$biom_monocotyledons, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(biom_monocotyledons ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$biom_monocotyledons, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Biomass forbs
# 
# # Against precipitation simple
# lm <- lm(biom_forbs ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$biom_forbs, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(biom_forbs ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$biom_forbs, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Biomass woody
# 
# # Against precipitation simple
# lm <- lm(biom_woody ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$biom_woody, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(biom_woody ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$biom_woody, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Diversity bryophytes
# 
# # Against precipitation simple
# lm <- lm(div_cryptogams ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$div_cryptogams, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(div_cryptogams ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$div_cryptogams, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Diversity monocotyledons
# 
# # Against precipitation simple
# lm <- lm(div_monocotyledons ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$div_monocotyledons, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(div_monocotyledons ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$div_monocotyledons, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Diversity forbs
# 
# # Against precipitation simple
# lm <- lm(div_forbs ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$div_forbs, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(div_forbs ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$div_forbs, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Diversity woody
# 
# # Against precipitation simple
# lm <- lm(div_woody ~ annualprecipitation, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$annualprecipitation, data_sheep_lm$div_woody, pch = 16, col = "blue")
# abline(lm)
# plot(lm$residuals)
# 
# # Against phosphorus simple
# lm <- lm(div_woody ~ phosphorus, data = data_sheep_lm)
# summary(lm)
# plot(data_sheep_lm$phosphorus, data_sheep_lm$div_woody, pch = 16, col = "green")
# abline(lm)
# plot(lm$residuals)
# 
# #
# ## Plot summary
# 
