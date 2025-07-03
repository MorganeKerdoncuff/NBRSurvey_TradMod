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
# library(FSA) # Dunn's Test after KS
# library(dunn.test) # Dunn's Test after KS -> "classic" but less handy to use, but to refer to
library(fpc) # clustering validation statistics
library(mri) # Adjusted Wallace index, non sensitive to cluster size (contrary to adjusted Rand index)
library(rcompanion) # effect size statistics
library(factoextra) # clustering visualisation
library(dendextend) # Visual dendrogram
library(ggplot2) # Visual representation
library(concaveman) # better hull
library(ggforce) # better ggplot
library(forcats) # reorder groups for plots
library(ggpubr) # Function ggarrange for several plots on same file
library(GGally) # Extension ggplot


#### DATA LOADING ####

siteinfo_full <- read.csv("data/cleandata/NBR_FullSiteInfo.csv", sep=",") # Clean site info data
landuse_full <- read.csv("data/cleandata/NBR_FullLandUse.csv", sep=",") # Clean site info data
plantspecies_full <- read.csv("data/cleandata/NBR_FullPlantComm.csv", sep=",") # Clean aboveground cover data
soilchem_full <- read.csv("data/cleandata/NBR_FullSoilChem.csv", sep=",") # Clean soil chemistry data
mesobio_full <- read.csv("data/cleandata/NBR_FullMesobio.csv", sep=",") # Clean soil mesofauna data
biomass_full <- read.csv("data/cleandata/NBR_FullBiomass.csv", sep=",") # Clean plant biomass data
climate_full <- read.csv("data/cleandata/NBR_FullClimate.csv", sep=",") # Clean climate data
soilchem_full <- read.csv("data/cleandata/NBR_FullSoilChem.csv", sep=",") # Clean soil chemistry data


#### DATA PREPARATION ####

## Filter sheep sites

# Site selection
siteinfo_sheep <- siteinfo_full |>  
  filter(Livestock == "sheep") |> 
  dplyr::select(SiteID, EcoZone) 
  
# Extraction in other datasets
landuse_sheep <- filter(landuse_full, SiteID %in% siteinfo_sheep$SiteID)
plantspecies_sheep <- filter(plantspecies_full, SiteID %in% siteinfo_sheep$SiteID)
soilchem_sheep <- filter(soilchem_full, SiteID %in% siteinfo_sheep$SiteID)
mesobio_sheep <- filter(mesobio_full, SiteID %in% siteinfo_sheep$SiteID)
biomass_sheep <- filter(biomass_full, SiteID %in% siteinfo_sheep$SiteID)
climate_sheep <- filter(climate_full, SiteID %in% siteinfo_sheep$SiteID)
soilchem_sheep <- filter(soilchem_full, SiteID %in% siteinfo_sheep$SiteID)

## Diversity - plant species richness according to functional group

# Summarise at site level, by average
plantspecies_sheep <- plantspecies_sheep |> 
  group_by(SiteID, Species) |> 
  summarise(Abundance = mean(Abundance)) |> 
  filter(Abundance > 0) |> 
  arrange(Species)

# New functional group variable
plantspecies_sheep <- plantspecies_sheep |>
  # assign functional group to species
  mutate(FunctionalGroup = case_when(
    grepl("Achillea|Alchemilla|Anemone|Bartsia|Campanula|Cardamine|Cerastium|Cirsium|Conopodium|Dactylorhiza|Digitalis|Epilobium|Euphrasia|Fraxinus|Galium|Gnaphalium|Hieracium|Hypericum|Hypochaeris|Lathyrus|Leontodon|Lotus|Melampyrum|Moneses|Narthecium|Oxalis|Pedicularis|Pinguicula|Plantago|Potentilla|Ranunculus|Rumex|Sagina|Sedum|Solidago|Stellaria|Succisa|Taraxacum|Trientalis|Trifolium|Veronica|Viola", Species) ~ "DIV_forbs",
    grepl("Agrostis|Aira|Alopecurus|Anthoxanthum|Carex|Dactylis|Danthonia|Deschampsia|Eriophorum|Festuca|Holcus|Juncus|Lolium|Luzula|Molinia|Nardus|Poa|Trichophorum", Species) ~ "DIV_monocotyledons",
    grepl("Andromeda|Arctostaphylos|Betula|Calluna|Chamaepericlymenum|Empetrum|Erica|Juniperus|Loiseleuria|Picea|Polygala|Polygonum|Populus|Prunus|Rubus|Salix|Sorbus|Ulmus|Vaccinium", Species) ~ "DIV_woody",
    grepl("Athyrium|Blechnum|Phegopteris|Polypodium|Pteridium", Species) ~ "DIV_ferns",
    .default = "DIV_cryptogams"
  ))

# Number of species per functional group per site
plantrichness_sheep <- plantspecies_sheep |>   
  # number of species per functional group
  group_by(SiteID, FunctionalGroup) |> 
  summarise(S = n()) |> 
  # make functional groups as variables
  pivot_wider(names_from = FunctionalGroup, values_from = S) |> 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
  # Total species richness
  mutate(DIV_tot = rowSums(across(where(is.numeric)), na.rm=TRUE))

## Aboveground production - plant biomass according to functional groups

# From biomass data
# Currently at sample level -> summary by average
biomass_sheep <- biomass_sheep |>
  pivot_wider(names_from = FunctionalType, values_from = Biomass.m2) |> 
  # make functional groups as variables
  group_by(SiteID) |>
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
  summarise(AGP_woody = mean(woody),
            AGP_monocotyledons = mean(monocotyledons),
            AGP_forbs = mean(forbs),
            AGP_cryptogams = mean(cryptogams))

  
## Carbon sequestration - LOI

# From chemistry data
# Same number of replicates per site
# Currently at plot level -> summary by average
soilchem_sheep <- soilchem_sheep |> 
  group_by(SiteID) |> 
  summarise(C_LOI = mean(LOI))

## Belowground production - mite & springtail abundance

# From mesofauna data
# Different numbers of replicates per site -> need 3 replicates, 1 per plot
# Currently at sample level -> summary by average
mesobio_sheep <- mesobio_sheep |> 
  group_by(SiteID, PlotID) |>
  # Select randomly one row which match unique combination of site & plot IDs
  slice(1) |>
  ungroup()

# average per site
meso_sheep <- mesobio_sheep |> 
  group_by(SiteID) |> 
  summarise(BGP_acari = mean(Acari),
            BGP_collembola = mean(Collembola))

## Initial base dataset

# Merge all data
data_sheep <- purrr::reduce(list(
  siteinfo_sheep,
  subset(landuse_sheep, select = c(SiteID, FieldType)),
  subset(plantrichness_sheep, select = -c(DIV_ferns)), 
  biomass_sheep, 
  meso_sheep, 
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

# Data scaling
data_sheep_sc <- as.data.frame(scale(subset(data_sheep, select = -c(SiteID, EcoZone, FieldType, fieldclass, FieldID, DIV_tot))))
rownames(data_sheep_sc) <- data_sheep$FieldID
#rownames(data_sheep_sc) <- data_sheep$SiteID


#### HIERARCHICAL CLUSTERING ####

## Distance matrix & clustering dendrogram

# Distance matrix - environmental continuous numeric -> Euclidean
dist_sheep <- dist(data_sheep_sc, method = "euclidean")

# Clustering using Ward's (1963) clustering criterion
hclustwd_sheep <- hclust(dist_sheep, method = "ward.D2")

# Show primary (k=2) & secondary (k=4) clusters with color_branches
dendro_sheep <- color_branches(as.dendrogram(hclustwd_sheep), 
                               k = 4, 
                               #groupLabels = TRUE,
                               col = c("pink", "pink4", "orange", "orange3"))
plot(dendro_sheep)

# Save dendrogram
png(filename = "outputs/cluster_tree.png", width = 14, height = 8, units = "cm", res = 600, pointsize = 9, bg = "transparent")
plot(dendro_sheep)
dev.off()

## Primary clusters quality check

# Observed & expected clusters
clusterwd_prim <- stats::cutree(hclustwd_sheep, k = 2)
clusterexp_prim <- ifelse(data_sheep$FieldType == "outfield", 1, 2)

# Clustering statistics
clusterstat_prim <- cluster.stats(dist_sheep, clusterwd_prim, clusterexp_prim)

## Secondary clusters quality check

# Observed & expected cluster
clusterwd_sec <- stats::cutree(hclustwd_sheep, k = 4)
clusterexp_sec <- ifelse(
  data_sheep$fieldclass == "CO" | data_sheep$fieldclass == "FO", 1,
  ifelse(
    data_sheep$fieldclass == "CI", 2,
    ifelse(
      data_sheep$fieldclass == "FI", 3, 4
    )
  )
)

# Clustering statistics
clusterstat_sec <- cluster.stats(dist_sheep, clusterwd_sec, clusterexp_sec)

## Summary table

# Empty table defining variables
clusterstat <- data.frame(
  Clustering = as.character(),
  NbCluster = as.numeric(),
  AvgBtwCluster = as.numeric(), # aim for highest
  AvgWthCluster = as.numeric(), # aim for lowest
  WthClusterSs = as.numeric(), # aim for lowest
  DunnIndex = as.numeric(), # aim for highest
  AdjRandIndex = as.numeric(), # aim for closest to 1 - sensible to cluster size
  AdjWallaceIndex = as.numeric() # aim for closest to 1 - non sensible to cluster size
)

# Adding rows
clusterstat <- clusterstat |>
  add_row(
    Clustering = "primary",
    NbCluster = clusterstat_prim$cluster.number,
    AvgBtwCluster = clusterstat_prim$average.between,
    AvgWthCluster = clusterstat_prim$average.within,
    WthClusterSs = clusterstat_prim$within.cluster.ss,
    DunnIndex = clusterstat_prim$dunn,
    AdjRandIndex = clusterstat_prim$corrected.rand,
    AdjWallaceIndex = mri::MAW1(clusterwd_prim, clusterexp_prim)
  ) |>
  add_row(
    Clustering = "secondary",
    NbCluster = clusterstat_sec$cluster.number,
    AvgBtwCluster = clusterstat_sec$average.between,
    AvgWthCluster = clusterstat_sec$average.within,
    WthClusterSs = clusterstat_sec$within.cluster.ss,
    DunnIndex = clusterstat_sec$dunn,
    AdjRandIndex = clusterstat_sec$corrected.rand,
    AdjWallaceIndex = mri::MAW1(clusterwd_sec, clusterexp_sec)
  )

## Ordination clusters

# NMDS
nmds_sheep <- metaMDS(dist_sheep)

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
  # ), 
  # nudge_y = 0.15,
  # size = 0.3,
  # size.unit = "cm") +
  theme(
    # axis.text.x = element_blank(),  # remove x-axis text
    #     axis.text.y = element_blank(), # remove y-axis text
    #     axis.ticks = element_blank(),  # remove axis ticks
        panel.background = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank()
    )
ggcluster_ord
ggsave(filename = "outputs/cluster_NMDS.png", plot = ggcluster_ord, width = 18, height = 10, units = "cm")


#### VARIABLE CORRELATION ####

# Data table preparation
data_sheep_cor <- data_sheep_sc |> 
  mutate(FieldID = rownames(data_sheep_sc))
data_sheep_cor <- left_join(data_sheep_cor, subset(data_sheep, select = c(FieldID, EcoZone, FieldType, DIV_tot)))

## Correlation panel

# Base lm plot function
plot_cor <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(
    aes(
      shape = EcoZone,
      colour = FieldType
      ), 
    size = 3) + 
    scale_colour_manual(values = c("infield" = "gray", "outfield" = "black")) +
    geom_smooth(method = "lm", color="navy", fill = "navy", se = TRUE) +
    theme(panel.background = element_blank(),
          legend.background = element_blank(),
          panel.grid.major = element_blank(),  #remove major-grid labels
          panel.grid.minor = element_blank(),  #remove minor-grid labels
          plot.background = element_blank())
  p
}

# Correlation panel
cor_panel <- ggpairs(
  data = data_sheep_cor,
  columns = c("AGP_woody", "AGP_forbs", "AGP_monocotyledons", "AGP_cryptogams", "DIV_woody", "DIV_forbs", "DIV_monocotyledons", "DIV_cryptogams", "DIV_tot", "BGP_acari", "BGP_collembola", "C_LOI"),
  columnLabels = c("AGPw", "AGPf", "AGPm", "AGPc", "DIVw", "DIVf", "Divm", "DIVc", "DIVt", "BGPa", "BGPc", "CLOI"),
  switch = "both",
  upper = list(continuous = wrap("cor", method = "spearman")),
  lower = list(continuous = plot_cor)
)
cor_panel

# Only significant correlation panel
corsig_panel <- ggpairs(
  data = data_sheep_cor,
  columns = c("AGP_woody", "AGP_forbs", "DIV_woody", "DIV_forbs", "DIV_cryptogams", "DIV_tot", "C_LOI"),
  columnLabels = c("AGPw", "AGPf", "DIVw", "DIVf", "DIVc", "DIVt", "CLOI"),
  switch = "both",
  upper = list(continuous = wrap("cor", method = "spearman")),
  lower = list(continuous = plot_cor)
)
corsig_panel
ggsave(filename = "outputs/corpanel_sigonly.png", plot = corsig_panel, width = 18, height = 12, units = "cm")

## Correlation stat for significant results

# Empty table defining variables
corstat <- data.frame(
  x = as.character(),
  y = as.character(),
  S.stat = as.numeric(),
  p.value = as.numeric(),
  rho = as.numeric()
)

# AGPw x DIVw
cor <- cor.test(data_sheep_cor$AGP_woody, data_sheep_cor$DIV_woody, method = "spearman")
corstat <- corstat |> 
  add_row(
    x = "AGP_woody",
    y = "DIV_woody",
    S.stat = cor$statistic,
    p.value = cor$p.value,
    rho = cor$estimate
  )

# AGPw x DIVf
cor <- cor.test(data_sheep_cor$AGP_woody, data_sheep_cor$DIV_forbs, method = "spearman")
corstat <- corstat |> 
  add_row(
    x = "AGP_woody",
    y = "DIV_forbs",
    S.stat = cor$statistic,
    p.value = cor$p.value,
    rho = cor$estimate
  )

# AGPw x DIVc
cor <- cor.test(data_sheep_cor$AGP_woody, data_sheep_cor$DIV_cryptogams, method = "spearman")
corstat <- corstat |> 
  add_row(
    x = "AGP_woody",
    y = "DIV_cryptogams",
    S.stat = cor$statistic,
    p.value = cor$p.value,
    rho = cor$estimate
  )

# AGPw x DIVtot
cor <- cor.test(data_sheep_cor$AGP_woody, data_sheep_cor$DIV_tot, method = "spearman")
corstat <- corstat |> 
  add_row(
    x = "AGP_woody",
    y = "DIV_tot",
    S.stat = cor$statistic,
    p.value = cor$p.value,
    rho = cor$estimate
  )

# AGPw x LOI
cor <- cor.test(data_sheep_cor$AGP_woody, data_sheep_cor$C_LOI, method = "spearman")
corstat <- corstat |> 
  add_row(
    x = "AGP_woody",
    y = "C_LOI",
    S.stat = cor$statistic,
    p.value = cor$p.value,
    rho = cor$estimate
  )

# Posthoc
# LOI outlier CI5 - very high LOI for infield category (area type verified) -> Risøyna ~ wet grassland
# verification relationships LOI with other factors excluding CI5 -> (+)* correlated with woody & crypto div & (-)* correlated with forb div
# filter(data_sheep_lm, FieldID != "CI5")
#ggpairs
# data_sheep_panel <- subset(filter(data_sheep_lm, FieldID != "CI5"), select = c(biom_woody, biom_forbs, div_woody, div_forbs, div_cryptogams, LOI, EcoZone, FieldType))


#### ANALYSIS OF CLUSTERS ####

# Extraction primary & secondary cluster group
data_sheep_sc <- data_sheep_sc |>
  mutate(cluster_prim = as.factor(clusterwd_prim)) |> 
  mutate(cluster_sec = as.factor(clusterwd_sec))

#
## Comparison primary clusters - unpaired non-parametric Mann-Whitney test

# Empty table for summary statistics
wcx_stat <- data.frame(cluster = as.character(), 
                       var = as.character(),
                       W = as.numeric(),
                       pval = as.numeric(),
                       rg = as.numeric(),
                       stringsAsFactors = FALSE)

# LOI
wcx <- wilcox.test(C_LOI ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "C_LOI", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$C_LOI, g = data_sheep_sc$cluster_prim))

# Abundance collembola
wcx <- wilcox.test(BGP_collembola ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "BGP_collembola", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$BGP_collembola, g = data_sheep_sc$cluster_prim))

# Abundance acari
wcx <- wilcox.test(BGP_acari ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "BGP_acari", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$BGP_acari, g = data_sheep_sc$cluster_prim))

# Biomass cryptogams
wcx <- wilcox.test(AGP_cryptogams ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "AGP_cryptogams", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$AGP_cryptogams, g = data_sheep_sc$cluster_prim))

# Biomass monocotyledons
wcx <- wilcox.test(AGP_monocotyledons ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "AGP_monocotyledons", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$AGP_monocotyledons, g = data_sheep_sc$cluster_prim))

# Biomass forbs
wcx <- wilcox.test(AGP_forbs ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "AGP_forbs", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$AGP_forbs, g = data_sheep_sc$cluster_prim))

# Biomass woody
wcx <- wilcox.test(AGP_woody ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "AGP_woody", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$AGP_woody, g = data_sheep_sc$cluster_prim))

# Diversity cryptogams
wcx <- wilcox.test(DIV_cryptogams ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "DIV_cryptogams", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$DIV_cryptogams, g = data_sheep_sc$cluster_prim))

# Diversity monocotyledons
wcx <- wilcox.test(DIV_monocotyledons ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "DIV_monocotyledons", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$DIV_monocotyledons, g = data_sheep_sc$cluster_prim))

# Diversity forbs
wcx <- wilcox.test(DIV_forbs ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "DIV_forbs", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$DIV_forbs, g = data_sheep_sc$cluster_prim))

# Diversity woody
wcx <- wilcox.test(DIV_woody ~ cluster_prim, data = data_sheep_sc)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "primary", var = "DIV_woody", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_sc$DIV_woody, g = data_sheep_sc$cluster_prim))

#
## Comparison secondary clusters outfield group (n=2) - unpaired non-parametric Mann-Whitney test

# # Empty table for summary statistics
# wcxsec_statA <- data.frame(var = as.character(),
#                            W = as.numeric(),
#                            pval = as.numeric(),
#                            rg = as.numeric(),
#                            stringsAsFactors = FALSE)

# Format data table adapted to rg test (group# 1-2)
data_sheep_scA <- filter(data_sheep_sc, cluster_sec == 1 | cluster_sec == 4)
data_sheep_scA <- data_sheep_scA |> 
  mutate(cluster = ifelse(data_sheep_scA$cluster_sec == 1, 1, 2))

# LOI
wcx <- wilcox.test(C_LOI ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "C_LOI", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$C_LOI, g = data_sheep_scA$cluster))

# Abundance collembola
wcx <- wilcox.test(BGP_collembola ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "BGP_collembola", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$BGP_collembola, g = data_sheep_scA$cluster))

# Abundance acari
wcx <- wilcox.test(BGP_acari ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "BGP_acari", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$BGP_acari, g = data_sheep_scA$cluster))

# Biomass cryptogams
wcx <- wilcox.test(AGP_cryptogams ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "AGP_cryptogams", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$AGP_cryptogams, g = data_sheep_scA$cluster))

# Biomass monocotyledons
wcx <- wilcox.test(AGP_monocotyledons ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "AGP_monocotyledons", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$AGP_monocotyledons, g = data_sheep_scA$cluster))

# Biomass forbs
wcx <- wilcox.test(AGP_forbs ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "AGP_forbs", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$AGP_forbs, g = data_sheep_scA$cluster))

# Biomass woody
wcx <- wilcox.test(AGP_woody ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "AGP_woody", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$AGP_woody, g = data_sheep_scA$cluster))

# Diversity cryptogams
wcx <- wilcox.test(DIV_cryptogams ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "DIV_cryptogams", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$DIV_cryptogams, g = data_sheep_scA$cluster))

# Diversity monocotyledons
wcx <- wilcox.test(DIV_monocotyledons ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "DIV_monocotyledons", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$DIV_monocotyledons, g = data_sheep_scA$cluster))

# Diversity forbs
wcx <- wilcox.test(DIV_forbs ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "DIV_forbs", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$DIV_forbs, g = data_sheep_scA$cluster))

# Diversity woody
wcx <- wilcox.test(DIV_woody ~ cluster_sec, data = data_sheep_scA)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryA", var = "DIV_woody", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scA$DIV_woody, g = data_sheep_scA$cluster))

#
## Comparison secondary clusters infield group (n=2) - unpaired non-parametric Mann-Whitney test

# Empty table for summary statistics
# wcxsec_statB <- data.frame(var = as.character(),
#                        W = as.numeric(),
#                        pval = as.numeric(),
#                        rg = as.numeric(),
#                        stringsAsFactors = FALSE)

# Format data table adapted to rg test (group# 1-2)
data_sheep_scB <- filter(data_sheep_sc, cluster_sec == 2 | cluster_sec == 3)
data_sheep_scB <- data_sheep_scB |> 
  mutate(cluster = ifelse(data_sheep_scB$cluster_sec == 2, 1, 2))

# LOI
wcx <- wilcox.test(C_LOI ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "C_LOI", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$C_LOI, g = data_sheep_scB$cluster))

# Abundance collembola
wcx <- wilcox.test(BGP_collembola ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "BGP_collembola", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$BGP_collembola, g = data_sheep_scB$cluster))

# Abundance acari
wcx <- wilcox.test(BGP_acari ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "BGP_acari", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$BGP_acari, g = data_sheep_scB$cluster))

# Biomass cryptogams
wcx <- wilcox.test(AGP_cryptogams ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "AGP_cryptogams", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$AGP_cryptogams, g = data_sheep_scB$cluster))

# Biomass monocotyledons
wcx <- wilcox.test(AGP_monocotyledons ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "AGP_monocotyledons", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$AGP_monocotyledons, g = data_sheep_scB$cluster))

# Biomass forbs
wcx <- wilcox.test(AGP_forbs ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "AGP_forbs", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$AGP_forbs, g = data_sheep_scB$cluster))

# Biomass woody
wcx <- wilcox.test(AGP_woody ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "AGP_woody", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$AGP_woody, g = data_sheep_scB$cluster))

# Diversity cryptogams
wcx <- wilcox.test(DIV_cryptogams ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "DIV_cryptogams", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$DIV_cryptogams, g = data_sheep_scB$cluster))

# Diversity monocotyledons
wcx <- wilcox.test(DIV_monocotyledons ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "DIV_monocotyledons", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$DIV_monocotyledons, g = data_sheep_scB$cluster))

# Diversity forbs
wcx <- wilcox.test(DIV_forbs ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "DIV_forbs", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$DIV_forbs, g = data_sheep_scB$cluster))

# Diversity woody
wcx <- wilcox.test(DIV_woody ~ cluster_sec, data = data_sheep_scB)
wcx_stat <- wcx_stat |> 
  add_row(cluster = "secondaryB", var = "DIV_woody", W = wcx$statistic, pval = wcx$p.value, rg = wilcoxonRG(data_sheep_scB$DIV_woody, g = data_sheep_scB$cluster))


#
## Comparison secondary clusters outfield group (n=3) - Kruskal-Wallis followed by Dunn's test

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

## Visual representation

# Make variable group - Aboveground productivity, Diversity, Belowground productivity, Carbon storage
data_sheep_plot <- data_sheep_sc |> 
  pivot_longer(cols = -c(cluster_prim, cluster_sec),
               names_to = "variables",
               values_to = "sc_val") |> 
  mutate(var_group = ifelse(
    grepl("DIV_", variables), "diversity", ifelse(
      grepl("BGP_", variables), "belowground production", ifelse(
        variables == "C_LOI", "C sequestration", "aboveground production"
      )
    )
  )) |> 
  # rearrange by variable group for ggplot
  arrange(sc_val) |>
  mutate(variables = factor(variables, levels = c("DIV_woody", "DIV_forbs", "DIV_monocotyledons", "DIV_cryptogams", "AGP_woody", "AGP_forbs", "AGP_monocotyledons", "AGP_cryptogams", "BGP_acari", "BGP_collembola", "C_LOI")))

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
  scale_x_discrete(limits = c("DIV_woody", "DIV_forbs", "DIV_monocotyledons", "DIV_cryptogams", "AGP_woody", "AGP_forbs", "AGP_monocotyledons", "AGP_cryptogams", "BGP_acari", "BGP_collembola", "C_LOI")) +
  annotate("text", x = "C_LOI", y = 0, label = "***", size = 2.5) +
  annotate("text", x = "BGP_collembola", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "BGP_acari", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_cryptogams", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_monocotyledons", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_forbs", y = 0, label = "**", size = 2.5) +
  annotate("text", x = "AGP_woody", y = 0, label = "***", size = 2.5) +
  annotate("text", x = "DIV_cryptogams", y = 0, label = "***", size = 2.5) +
  annotate("text", x = "DIV_monocotyledons", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "DIV_forbs", y = 0, label = "**", size = 2.5) +
  annotate("text", x = "DIV_woody", y = 0, label = "***", size = 2.5)
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
  scale_x_discrete(limits = c("DIV_woody", "DIV_forbs", "DIV_monocotyledons", "DIV_cryptogams", "AGP_woody", "AGP_forbs", "AGP_monocotyledons", "AGP_cryptogams", "BGP_acari", "BGP_collembola", "C_LOI")) +
  annotate("text", x = "C_LOI", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "BGP_collembola", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "BGP_acari", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_cryptogams", y = 0, label = "*", size = 2.5) +
  annotate("text", x = "AGP_monocotyledons", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_forbs", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_woody", y = 0, label = "*", size = 2.5) +
  annotate("text", x = "DIV_cryptogams", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "DIV_monocotyledons", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "DIV_forbs", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "DIV_woody", y = 0, label = "°", size = 2.5)
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
  scale_x_discrete(limits = c("DIV_woody", "DIV_forbs", "DIV_monocotyledons", "DIV_cryptogams", "AGP_woody", "AGP_forbs", "AGP_monocotyledons", "AGP_cryptogams", "BGP_acari", "BGP_collembola", "C_LOI")) +
  annotate("text", x = "C_LOI", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "BGP_collembola", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "BGP_acari", y = 0, label = "*", size = 2.5) +
  annotate("text", x = "AGP_cryptogams", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_monocotyledons", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_forbs", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "AGP_woody", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "DIV_cryptogams", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "DIV_monocotyledons", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "DIV_forbs", y = 0, label = "ns", size = 2.5) +
  annotate("text", x = "DIV_woody", y = 0, label = "ns", size = 2.5)
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
