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
library(forcats) # reorder groups for plots
library(ggpubr) # Function ggarrange for several plots on same file
# library(GGally) # Extension ggplot


#### DATA LOADING ####

siteinfo_full <- read.csv("data/cleandata/NBR_FullSiteInfo.csv", sep=",") # Clean site info data
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
  dplyr::select(SiteID)

# Extraction in other datasets
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
    grepl("Achillea|Alchemilla|Anemone|Bartsia|Campanula|Cardamine|Cerastium|Cirsium|Conopodium|Dactylorhiza|Digitalis|Epilobium|Fraxinus|Galium|Gnaphalium|Hieracium|Hypericum|Hypochaeris|Lathyrus|Leontodon|Lotus|Melampyrum|Moneses|Narthecium|Oxalis|Pedicularis|Pinguicula|Plantago|Potentilla|Ranunculus|Rumex|Sagina|Sedum|Solidago|Stellaria|Succisa|Taraxacum|Trientalis|Trifolium|Veronica|Viola", Species) ~ "div_forbs",
    grepl("Agrostis|Aira|Alopecurus|Anthoxanthum|Carex|Dactylis|Danthonia|Deschampsia|Eriophorum|Festuca|Holcus|Juncus|Lolium|Luzula|Molinia|Nardus|Poa|Trichophorum", Species) ~ "div_monocotyledons",
    grepl("Andromeda|Arctostaphylos|Betula|Calluna|Chamaepericlymenum|Empetrum|Erica|Juniperus|Loiseleuria|Picea|Polygala|Polygonum|Populus|Prunus|Rubus|Salix|Sorbus|Ulmus|Vaccinium", Species) ~ "div_woody",
    grepl("Athyrium|Blechnum|Phegopteris|Polypodium|Pteridium", Species) ~ "div_ferns",
    .default = "div_bryophytes"
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
  filter(FunctionalGroup == "div_bryophytes")

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

# Site info - validated
# Climate - validated

# Plant species richness
## From groundcover data
## Same number of replicates per site
## Currently at sample level -> summary by average
# groundcover_sheep <- groundcover_sheep |>  
#   group_by(SiteID) |>  
#   summarise(MeanRichness = mean(Plant_species_richness, na.rm=TRUE))

# Loss of ignition
## From chemistry data
## Same number of replicates per site
## Currently at plot level -> summary by average
soilchem_sheep <- soilchem_sheep |> 
  group_by(SiteID) |> 
  summarise(LOI = mean(LOI), phophorus = mean(P.Al_mg.100g))

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
  summarise(biom_woody = mean(woody),
            biom_monocotyledons = mean(monocotyledons),
            biom_forbs = mean(forbs),
            biom_cryptogams = mean(cryptogams))

#
## Initial base dataset

# Merge all data
data_sheep <- purrr::reduce(list(
  siteinfo_sheep, 
  subset(plantrichness_sheep, select = -c(div_ferns)), 
  biomass_sheep, 
  mesobio_sheep, 
  soilchem_sheep, 
  subset(climate_sheep, select = c(SiteID, annualprecipitation, avgtempJuly))
  ), dplyr::left_join)

# Intuitive ID according to following code - COx (costal outfield), Ix (infield), MOx (mountain oufield)
data_sheep <- data_sheep |> 
  mutate(fieldtype = case_when(
    grepl("U", SiteID) ~ "MO",
    grepl("IS2|OV|OS5|OS7|OS9", SiteID) ~ "CO",
    .default = "IN"
  )) |> 
  group_by(fieldtype) |>
  mutate(FieldID = paste(fieldtype, row_number(), sep = ""))

# Site ID as rowname
rownames(data_sheep) <- data_sheep$FieldID



#### HIERARCHICAL CLUSTERING ####

# Data scaling
data_sheep_sc <- as.data.frame(scale(subset(data_sheep, select = -c(SiteID, fieldtype, FieldID, annualprecipication, avgtempJuly, phosphorus))))
rownames(data_sheep_sc) <- data_sheep$FieldID

#
## Clustering & dendrogram

# Distance matrix - environmental continuous numeric -> Euclidean
dist_sheep <- dist(data_sheep_sc, method = "euclidean")

# Clustering using Ward's (1963) clustering criterion
hclustwd_sheep <- hclust(dist_sheep, method = "ward.D2")

# Plot HC tree
plot(hclustwd_sheep)

# Show clusters (k=5) with base
rect.hclust(hclustwd_sheep, k = 4, border = 2:6)

# Show clusters (k=5) with color_branches
dendro_sheep <- color_branches(as.dendrogram(hclustwd_sheep), 
                               k = 4, 
                               #groupLabels = TRUE,
                               col = c("pink", "pink3", "orange", "orange3"))
plot(dendro_sheep)

# Cluster group
clusterwd <- stats::cutree(hclustwd_sheep, k = 4)

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
## Ordination clusters

# NMDS -> one-dimension
nmds_sheep <- metaMDS(dist_sheep)
# ordiplot(nmds_sheep, type = 'n')
# points(nmds_sheep, pch = clusterwd+20, bg = clusterwd)
# legend('topright', legend = 1:5, pch = (1:5)+20, pt.bg = 1:5)

# Colour grouping
#clusterwd # see order
col_clusters <- c("pink3", "orange4", "orange", "pink")
col_clusters[clusterwd]

# Plot ordination with base
cluster_ord <- ordiplot(nmds_sheep, type = 'n')
cluster_ord <- points(nmds_sheep, 
       col = col_clusters[clusterwd], 
       pch = clusterwd+20)
cluster_ord <- text(nmds_sheep, 
     col = col_clusters[clusterwd], 
     labels=rownames(data_sheep_sc),
     cex = 0.7,
     pos = 1)
cluster_ord <- legend( 
       legend = paste("Cluster", 1:4),
       x = 2.3,
       y = 3.5,
       #horiz = TRUE,
       text.width = 1,
       x.intersp = 1,
       y.intersp = 1,
       col = col_clusters, 
       pt.bg = col_clusters,
       bty = "n", 
       pch = (1:4)+20)
cluster_ord <- ordihull(nmds_sheep, 
         groups = clusterwd, 
         display = "sites")
cluster_ord



#### ANALYSIS OF CLUSTERS ####

#
## Data preparation

# Extraction cluster group
data_sheep_sc <- data_sheep_sc |>
  mutate(cluster = as.numeric(clusterwd))

# Make variable group - Productivity, Diversity, Carbon, Soil nutrients
data_sheep_plot <- data_sheep_sc |> 
  pivot_longer(cols = -c(cluster),
               names_to = "variables",
               values_to = "sc_val") |> 
  mutate(var_group = ifelse(
    grepl("div_", variables), "diversity", ifelse(
      variables == "ab_acari" | variables == "ab_collembola" | variables == "LOI", "decomposition", "production"
    )
  )) |> 
  # rearrange by variable group for ggplot
  arrange(sc_val) |>
  mutate(variables = factor(variables, levels = c("div_woody", "div_forbs", "div_monocotyledons", "div_bryophytes", "biom_woody", "biom_forbs", "biom_monocotyledons", "biom_cryptogams", "ab_acari", "ab_collembola", "LOI")))

#
## Comparison with factor cluster - Kruskal-Wallis followed by Dunn's test

# Make cluster variable as factor
data_sheep_plot <- data_sheep_plot |> 
  mutate(cluster = as.factor(cluster))

# LOI
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "LOI")) # significant **
dunnTest(sc_val ~ cluster, 
          data = filter(data_sheep_plot, variables == "LOI"),
          method = "bonferroni")

# Abundance collembola
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "ab_collembola")) #NS
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "ab_collembola"),
         method = "bonferroni")

# Abundance acari
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "ab_acari")) # significant *
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "ab_acari"),
         method = "bonferroni")

# Biomass cryptogams
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "biom_cryptogams")) # significant *
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "biom_cryptogams"),
         method = "bonferroni")

# Biomass monocotyledons
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "biom_monocotyledons")) #NS
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "biom_monocotyledons"),
         method = "bonferroni")

# Biomass forbs
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "biom_forbs")) # significant *
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "biom_forbs"),
         method = "bonferroni")

# Biomass woody
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "biom_woody")) # significant **
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "biom_woody"),
         method = "bonferroni")

# Diversity bryophytes
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "div_bryophytes")) # significant **
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "div_bryophytes"),
         method = "bonferroni")

# Diversity monocotyledons
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "div_monocotyledons")) # NS
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "div_monocotyledons"),
         method = "bonferroni")

# Diversity forbs
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "div_forbs")) # significant *
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "div_forbs"),
         method = "bonferroni")

# Diversity woody
kruskal.test(sc_val ~ cluster, data = filter(data_sheep_plot, variables == "div_woody")) # significant ***
dunnTest(sc_val ~ cluster, 
         data = filter(data_sheep_plot, variables == "div_woody"),
         method = "bonferroni")

#
## Boxplots

# Cluster 1
clusterplot_1 <- filter(data_sheep_plot, cluster == 1) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster 1") +
  labs("") +
  ylim(-3.5, 3.5) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "pink3"),
    legend.title = element_blank()
    ) +
  annotate("text", x = "LOI", y = 3.5, label = "a*", size = 3) +
  annotate("text", x = "ab_acari", y = 3.5, label = "ab", size = 3) +
  annotate("text", x = "biom_cryptogams", y = 3.5, label = "a", size = 3) +
  annotate("text", x = "biom_forbs", y = 3.5, label = "ab", size = 3) +
  annotate("text", x = "biom_woody", y = 3.5, label = "a", size = 3) +
  annotate("text", x = "div_bryophytes", y = 3.5, label = "a", size = 3) +
  annotate("text", x = "div_forbs", y = 3.5, label = "a", size = 3) +
  annotate("text", x = "div_woody", y = 3.5, label = "a", size = 3) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_1

# Cluster 2
clusterplot_2 <- filter(data_sheep_plot, cluster == 2) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster 2") +
  ylim(-4,4) +
  labs("") +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "orange4"),
    legend.title = element_blank()
  ) +
  annotate("text", x = "LOI", y = 4, label = "ab", size = 3) +
  annotate("text", x = "ab_acari", y = 4, label = "b*", size = 3) +
  annotate("text", x = "biom_cryptogams", y = 4, label = "ab", size = 3) +
  annotate("text", x = "biom_forbs", y = 4, label = "b¤", size = 3) +
  annotate("text", x = "biom_woody", y = 4, label = "b*", size = 3) +
  annotate("text", x = "div_bryophytes", y = 4, label = "bc", size = 3) +
  annotate("text", x = "div_forbs", y = 4, label = "ab", size = 3) +
  annotate("text", x = "div_woody", y = 4, label = "b*c", size = 3) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_2

# Cluster 3
clusterplot_3 <- filter(data_sheep_plot, cluster == 3) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster 3") +
  labs("") +
  ylim(-3,3) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "orange"),
    legend.title = element_blank()
  ) +
  annotate("text", x = "LOI", y = 2.5, label = "b", size = 3) +
  annotate("text", x = "ab_acari", y = 2.5, label = "a", size = 3) +
  annotate("text", x = "biom_cryptogams", y = 2.5, label = "ab", size = 3) +
  annotate("text", x = "biom_forbs", y = 2.5, label = "ab", size = 3) +
  annotate("text", x = "biom_woody", y = 2.5, label = "b**", size = 3) +
  annotate("text", x = "div_bryophytes", y = 2.5, label = "b**", size = 3) +
  annotate("text", x = "div_forbs", y = 2.5, label = "b¤", size = 3) +
  annotate("text", x = "div_woody", y = 2.5, label = "b***", size = 3) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_3

# Former Cluster 4
# clusterplot_4 <- filter(data_sheep_plot, cluster == 4) |> 
#   ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   geom_boxplot() +
#   xlab("") +
#   ylab("Cluster 4") +
#   ylim(-3,3) +
#   coord_flip() +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.x = element_text(colour = "pink4"),
#     legend.title = element_blank()
#   ) +
#   scale_fill_grey(start = 0.5, end = 1)
# clusterplot_4

# Cluster 4
clusterplot_4 <- filter(data_sheep_plot, cluster == 4) |> 
  ggplot(aes(x = variables, y = sc_val, fill = var_group)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot() +
  xlab("") +
  ylab("Cluster 4") +
  labs("") +
  ylim(-3,3) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(colour = "pink"),
    legend.title = element_blank()
  ) +
  annotate("text", x = "LOI", y = 2.5, label = "a**", size = 3) +
  annotate("text", x = "ab_acari", y = 2.5, label = "ab", size = 3) +
  annotate("text", x = "biom_cryptogams", y = 2.5, label = "b*", size = 3) +
  annotate("text", x = "biom_forbs", y = 2.5, label = "a", size = 3) +
  annotate("text", x = "biom_woody", y = 2.5, label = "ab", size = 3) +
  annotate("text", x = "div_bryophytes", y = 2.5, label = "ac**", size = 3) +
  annotate("text", x = "div_forbs", y = 2.5, label = "ab", size = 3) +
  annotate("text", x = "div_woody", y = 2.5, label = "ab*c¤", size = 3) +
  scale_fill_grey(start = 0.5, end = 1)
clusterplot_4

# x-axis for plot arrange
clusterplot_axis <- filter(data_sheep_plot, cluster == 1) |> 
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
    legend.title = element_blank()
  )
# scale_fill_grey()
clusterplot_axis

# Arrange all plots
clusterplot_all <- ggarrange(clusterplot_axis, clusterplot_1, clusterplot_4, clusterplot_3, clusterplot_2, nrow = 1, common.legend = TRUE)
clusterplot_all
ggsave(filename = "outputs/clusterplots2.png", plot = clusterplot_all, height = 12, width = 20, units = "cm")


#### LM AGAINST TEMP, PRECI & PHOSPHORUS ####

# Significant differences between groups

# New table for LM
data_sheep_lm <- subset(data_sheep, select = -c(fieldtype, SiteID))

# LOI
lm <- lm(LOI ~ annualprecipitation, data = data_sheep_lm)
summary(lm)

# Abundance collembola

# Abundance acari

# Biomass cryptogams

# Biomass monocotyledons

# Biomass forbs

# Biomass woody

# Diversity bryophytes

# Diversity monocotyledons

# Diversity forbs

# Diversity woody


#### DATA SUMMARY ####

#
## Data preparation

# Extraction column with SiteID
data_sheep_sc <- data_sheep_sc |> 
  mutate(FieldID = row.names(data_sheep_sc))

# Join data
sheep_summary <- full_join(data_sheep, subset(data_sheep_sc, select = c(FieldID, cluster)))

#
## Summary -> need to fix the "across col" to exclude non-numeric

sheep_summary <- sheep_summary |> 
  dplyr::group_by(cluster) |> 
  select_if(is.numeric) |> 
  summarise(across(.cols = everything(), list(mean = mean, sd = sd)))