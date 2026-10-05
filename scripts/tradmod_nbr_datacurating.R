##
#### DESCRIPTION ####
##
## Purpose: Data curating script for the ecological survey led in summers 2019 and 2020 in the Nordhordland UNESCO Biosphere Reserve
## Author: Morgane KERDONCUFF
## ORCID: 0000-0003-2223-1857
## github
## Date created: 06/2026
## Last time modified: 06/2026
## Project: TradMod
## Funding: Norges forskningsråd (NFR)
## Institution: University of Bergen, Norway
##
##


#### PACKAGES ####

library(tidyverse) # R language
library(purrr) # Merge tables

#### SITE DESCRIPTION ####

## Numeric variables - min/max, distribution, potential outliers

### Min/max
# test <- site_description |>  
#   summarise(
#     tibble(
#       across(
#         where(is.numeric),
#         ~min(.x, na.rm = TRUE),
#         .names = "min_{.col}"
#         ),
#       across(
#         where(is.numeric),
#         ~max(.x, na.rm = TRUE),
#         .names = "max_{.col}")
#       )
#     ) |>  
#   transpose() #validated

### Variable distribution & outliers
# hist(site_description$numberAnimalsAdult) # Poisson distribution, one outlier over 150 animals
# site_description[site_description$numberAnimalsAdult>150,] # is4, no observed inconsistency with farm characteristics
# hist(site_description$numberAnimalsYoung) # Poisson distribution, no outlier
hist(site_description$fieldAreaHa) # One outlier over 1000 ha crushing the distribution
site_description[site_description$fieldAreaHa>1000,] # ug2 outfield site in upland area, value non applicable for analysis
# hist(filter(site_description, siteID != "ug2")$fieldAreaHa) # Poisson distribution
# hist(site_description$farmGrazingAreaHa) # dominance small farms, no visible outliers

#### 20x20 SAMPLING AREA ####

## Numeric variables - min/max, distribution, potential outliers

## Min/max
# test <- sampling_area |>
#   summarise(
#     tibble(
#       across(
#         where(is.numeric),
#         ~min(.x, na.rm = TRUE),
#         .names = "min_{.col}"
#         ),
#       across(
#         where(is.numeric),
#         ~max(.x, na.rm = TRUE),
#         .names = "max_{.col}")
#       )
#     ) |>
#   transpose() #validated - need further check for max herbs 97% & max lichens 20%

### NA check
# colnames(sampling_area)[apply(sampling_area, 2, anyNA)] # col: numberLivestockPaths, lengthLivestockPaths & all percent cover
sampling_area[!complete.cases(sampling_area),] # missing data for three sites (us1, ug1, oc4)

### Variable distribution & outliers
table(sampling_area$numberLivestockPaths) # dominance 0, variable to be taken out
hist(sampling_area$lengthLivestockPathM) # dominance 0, variable to be taken out
# hist(sampling_area$elevationMasl) # Poisson distribution, no outlier
# hist(sampling_area$slopeAngleDegree) # Normal distribution, no outlier
hist(sampling_area$slopeAspectDegree) # Uneven distribution, no outlier
# hist(sampling_area$percentRock) # Poisson distribution, further check for sites over 7%
# sampling_area[sampling_area$percentRock>7,] # subalpine heathlands (ug2, us3, us5), coastal heatland (ov1) & fjord at higher elevation (ig3)
# hist(sampling_area$percentMud) # very skewed Poisson distribution, no outlier
# hist(sampling_area$percentTreesTallShrubs) # Skewed Poisson distribution, no outlier
# hist(sampling_area$percentLowShrubs) # Poisson distribution, further check all grassland sites should be under 10%
# sampling_area[sampling_area$percentLowShrubs>10,] # 11 heathland sites
hist(sampling_area$percentForbs) # Poisson distribution, further check for sites >50%
sampling_area[sampling_area$percentForbs>50,] # os1, oc1, ig1, ig2, is2, iv1, ic1, og2, ic4 -> all first year/starting sites, check on vegetation quadrats + site & plot pictures
# OS1 80% - average 20% & no cover over 55% in quadrats, estimation from pictures 35%-40%
# OC1 80% - average 25% & no cover over 50% in quadrats, estimation from pictures 10%-15% 
# IG1 97% - average 20% & no cover over 30% in quadrats, estimation from pictures 15%-20%
# IS2 80% - average 70% & no cover over 90% in quadrats, estimation from pictures 45%-50%
# IC1 70% - average 45% & no cover over 70% in quadrats, estimation from pictures 60%-65%
hist(sampling_area$percentMonocotyledons) # further check for sites with odd forb distribution
hist(sampling_area$percentBryophytes) # uneven distribution, further check for sites with odd forb distribution
# hist(sampling_area$percentLichens) # highly skewed Poisson distribution, one site over 10%
# sampling_area[sampling_area$percentLichens>10,] # US4 in subalpine area, average of 12% & max 24% in quadrats - validated

# Add/remove variables

## Removal numberLivestockPaths and lengthLivestockPathM due to 
# sampling_area <- subset(sampling_area, select = -c(numberLivestockPaths, lengthLivestockPathM))

## Heat Load Index
# sampling_area <- sampling_area |> 
#   mutate(heatLoadIndex = cos(slopeAspectDegree-225)*tan(slopeAngleDegree))
# hist(sampling_area$heatLoadIndex) # 3 outliers: one under 200, two over 100
# sampling_area[sampling_area$heatLoadIndex>100,] #OG4 & IS3 -> both 11 degree slope with SW & SE exposition
# sampling_area[sampling_area$heatLoadIndex<0,] #OS6 -> 11 degree slope with NE exposition


#### NON-DESTRUCTIVE SUBPLOTS - GROUND COVER ####

## Numeric variables - min/max, distribution, potential outliers

## Min/max
# test <- ground_cover |>
#   summarise(
#     tibble(
#       across(
#         where(is.numeric),
#         ~min(.x, na.rm = TRUE),
#         .names = "min_{.col}"
#         ),
#       across(
#         where(is.numeric),
#         ~max(.x, na.rm = TRUE),
#         .names = "max_{.col}")
#       )
#     ) |>
#   transpose() #validated - need further check for max lichens 80%

### NA check
# colnames(ground_cover)[apply(ground_cover, 2, anyNA)] #validated - only blossom species ID

### Variable distribution & outliers
# hist(ground_cover$percentBareGround) # skewed Poisson distribution, check subplots > 20%
# filter(ground_cover, percentBareGround>20) #validated - subplots from recently burnt heathland ov1
# hist(ground_cover$percentRock) # skewed Poisson distribution
# hist(ground_cover$percentLitter) # skewed Poisson distribution, check subplots > 30%
# filter(ground_cover, percentLitter>30) #validated - subplots from burnt (ov1) & mountain heathlands (us2, ug1)
# hist(ground_cover$percentDeadWood) # skewed Poisson distribution, check subplots > 2%
# filter(ground_cover, percentDeadWood>2) #validated - subplots from is5, in the middle of a wood clearing
# hist(ground_cover$percentBryophytes) # Poisson distribution, no outliers
# hist(ground_cover$percentLichens) # Poisson distribution, one outlier above 40%
# filter(ground_cover, percentLichens>40) #validated - mountain site (us1-p1-n5) with high lichen cover
# hist(ground_cover$percentVascular) # Exponential distribution, no outlier
# hist(ground_cover$percentBlossom) # Skewed Poisson distribution, no outlier
# hist(ground_cover$percentDung) # Skewed Poisson distribution, check subplots > 10%
# filter(ground_cover, percentDung>10) #validated - cattle site (ic1)
# hist(ground_cover$avgVegetationHeightCm) # Poisson distribution, no outliers
# hist(ground_cover$maxVegetationHeightCm) # Normal distribution, no outliers


#### DESTRUCTIVE SUBPLOTS - SOIL PENETRATION TESTS ####

### Categories & distribution
# table(soil_pene$siteID) #validated - uc1 site (bog) to be removed
soil_pene <- filter(soil_pene, siteID != "uc1")
# table(soil_pene$plotID) #validated
# table(soil_pene$bedrockHit) # 39 failed tests over 1017 due to bedrock hit
table(filter(soil_pene, bedrockHit == "y")$siteID) # 8 sites with up to 12 failures

## Numeric variables - min/max, distribution, potential outliers

### Min/max
# test <- soil_pene |>
#   summarise(
#     tibble(
#       across(
#         where(is.numeric),
#         ~min(.x, na.rm = TRUE),
#         .names = "min_{.col}"
#         ),
#       across(
#         where(is.numeric),
#         ~max(.x, na.rm = TRUE),
#         .names = "max_{.col}")
#       )
#     ) |>
#   transpose() #validated - no visible height above maximum stick length

### NA check
# colnames(soil_pene)[apply(soil_pene, 2, anyNA)] #validated

### Variable distribution & outliers
# hist(soil_pene$visibleHeightCm) # Normal distribution, no outliers

# Add/remove variables
## New variable soilPenetrationDepth
# soil_pene <- soil_pene %>% 
#   mutate(soilPeneDepthCm = stickLengthCm - visibleHeightCm)
## Removal stickLengthCm and visibleHeightCm
# soil_pene <- subset(soil_pene, select = -c(stickLengthCm, visibleHeightCm))


#### DESTRUCTIVE SUBPLOTS - BULK DENSITY & GRAVIMETRIC WATER CONTENT ####

## Numeric variables - min/max, distribution, potential outliers

### Min/max
# test <- soil_bulk |>
#   summarise(
#     tibble(
#       across(
#         where(is.numeric),
#         ~min(.x, na.rm = TRUE),
#         .names = "min_{.col}"
#         ),
#       across(
#         where(is.numeric),
#         ~max(.x, na.rm = TRUE),
#         .names = "max_{.col}")
#       )
#     ) |>
#   transpose()

### NA check
# colnames(soil_bulk)[apply(soil_bulk, 2, anyNA)] # all variable, check row identification
soil_bulk[!complete.cases(soil_bulk),] # two missing records (is1-p3-d4-r2 & og1-p3-d2-r1) -> discarded due to lab incident
# soil_bulk <- filter(soil_bulk, recordID != "is1-p3-d4-r2" & recordID != "og1-p3-d2-r1")

### Quality check

#### Soil core volume - low core volume affect weight measurements
# hist(soil_bulk$coreVolCm3) # Visible threshold around 50 cm3
# filter(soil_bulk, coreVolCm3<40) # 37 or 2% samples unfit
# filter(soil_bulk, coreVolCm3<45) # 105 or 7% samples unfit
# filter(soil_bulk, coreVolCm3<50) # 330 or 21% samples unfit -> cores should be minimum vol of 50 cm3
# soil_bulk <- filter(soil_bulk, coreVolCm3 >= 50)

### Soil core weight - weight loss should be consistent with drying processes
# qualitycheck <- filter(soil_bulk,
#                         weightSatG - weight0hG < 0 |
#                         weight24hG - weightSatG > 0 |
#                         weight48hG - weight24hG > 0 |
#                         weightDryG - weight0hG > 0)

### Variable distribution & outliers

# W0 - Weight of fresh soil before saturation
#soilbulk_full[is.na(soilbulk_full$W0g),] # same NAs -> validated
#hist(soilbulk_full$W0g) # Weights range from 50 to 180 g in a normal distribution -> very low weights likely to be linked to low volumes

# WSAT - Soil weight after water saturation
#soilbulk_full[is.na(soilbulk_full$WSAT),] # same NAs -> validated
#hist(soilbulk_full$WSAT) # Weights range from 60 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes

# W24H - Soil weight after 24h of drying
#soilbulk_full[is.na(soilbulk_full$W24H),] # same NAs -> validated
#hist(soilbulk_full$W24H) # Weights range from 50 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes, not so much difference compared to WSAT

# W48H - Soil weight after 48h of drying
#soilbulk_full[is.na(soilbulk_full$W48H),] # same NAs -> validated
#hist(soilbulk_full$W48H) # Weights range from 50 to 190 g in a normal distribution -> very low weights likely to be linked to low volumes

# WDRY - Soil weight after over at 105C
#soilbulk_full[is.na(soilbulk_full$WDRY),] # same NAs -> validated
#hist(soilbulk_full$WDRY) # Weights range from 20 to 140 g in a normal distribution -> very low weights likely to be linked to low volumes

# Percent water loss in 24h
#soilbulk_full[is.na(soilbulk_full$percent_Waterloss24h),] # same NAs -> validated
#hist(soilbulk_full$percent_Waterloss24h) # % range from -40% to 40%, main between 0 and 10% -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Percent water loss in 48h
#soilbulk_full[is.na(soilbulk_full$percent_Waterloss48h),] # same NAs -> validated
#hist(soilbulk_full$percent_Waterloss48h) # % range from -70% to 70%, main between 0 and 20% -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Bulk density
#soilbulk_full[is.na(soilbulk_full$BD),] # same NAs -> validated
#hist(soilbulk_full$BD) # % range from -0.2% to 3, normal distribution -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Soil moisture in percentage weight = gravimetric water content
#soilbulk_full[is.na(soilbulk_full$Weightpercent_Soilmoisture),] # same NAs -> validated
#hist(soilbulk_full$Weightpercent_Soilmoisture) # % range from -50% to 100%, normal distribution -> negative and extreme values might be linked to processing issue (scale) or low soil volume

# Soil moisture in percentage volume
#soilbulk_full[is.na(soilbulk_full$Volpercent_Soilmoisture),] # same NAs -> validated
#hist(soilbulk_full$Volpercent_Soilmoisture) # % range from -50% to 100%, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

# Soil porosity
#soilbulk_full[is.na(soilbulk_full$percent_Soilporosity),] # same NAs -> validated
#hist(soilbulk_full$percent_Soilporosity) # % range from -70% to 80%, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

# WFPS
#soilbulk_full[is.na(soilbulk_full$percent_WFPS),] # same two NAs -> validated
#hist(soilbulk_full$percent_WFPS) # % range from -50% to one outlier over 10000, normal distribution -> negative values might be linked to processing issue (scale) or low soil volume

#
## New variables - gravimetric & volumetric water content from standardised W+48h dried soil

# New variables
# soilbulk_full <- soilbulk_full |> 
#   mutate(GWC_48 = (W48H - WDRY)/W48H*100) |> 
#   mutate(VWC_48 = GWC_48*BD)

# Distribution
# hist(soilbulk_full$GWC_48) # Normal distribution, from 20% to 80% -> some very high values
# hist(soilbulk_full$VWC_48) # Normal distribution, from 5% to 60%

#
## Data filtering

# Min soil core volume
#filter(soilbulk_full, CoreVol<40 & !is.na(CoreVol)) # 40 or 2% samples unfit
#filter(soilbulk_full, CoreVol<45 & !is.na(CoreVol)) # 111 or 7% samples unfit
#filter(soilbulk_full, CoreVol<50 & !is.na(CoreVol)) # 344 or 21% samples unfit -> cores should be minimum vol of 50 cm3

# Negative water loss values
#filter(soilbulk_full, percent_Waterloss24h<0 & !is.na(percent_Waterloss24h)) #132 samples with water loss 24h negative
#filter(soilbulk_full, percent_Waterloss48h<0 & !is.na(percent_Waterloss48h)) #82 samples with water loss 48h negative

# Selection data with min 50 cm3 soil volume and positive water loss
# soilbulk_full <- subset(soilbulk_full, CoreVol>50)
# soilbulk_full <- subset(soilbulk_full, percent_Waterloss24h>0)
# soilbulk_full <- subset(soilbulk_full, percent_Waterloss48h>0)

# Check new variable distribution - water loss 24h
#hist(soilbulk_full$percent_Waterloss24h) # still some extreme values over 20%
#filter(soilbulk_full, percent_Waterloss24h>20) # 6 cores with more than 20% over 24h
# 2 cores from UC1, which is excluded from the analysis -> should be removed
# 1 cores from OC3, concerned with scale issue (lots of negative values which are already removed). Water loss between 0-24 and 24-48 not coherent -> should be removed
# 2 cores from OC2, concerned with scale issue. Water loss between 0-24 and 24-48 not coherent with other samples from same plot (W48h>W24h for P1-D1_2) -> should be removed
# 1 cores from OC5, concerned with scale issue. Water loss between 0-24 and 24-48 not coherent with other samples from same site -> should be removed
# soilbulk_full <- subset(soilbulk_full, percent_Waterloss24h<20)

# Check new variable distribution - water loss 48h
#hist(soilbulk_full$percent_Waterloss48h) # still some extreme values over 25%
#filter(soilbulk_full, percent_Waterloss48h>25) # 2 cores with more than 25% over 48h
# OG6-P1-D3_1, not concerned by the scale issue and with values from other cores coherent -> to be kept
# OC2-P2-D1_3, concerned with scale issue - value not coherent with water loss 24h and with other cores -> to be removed
# soilbulk_full <- subset(soilbulk_full, BDcoreID != "OC2-P2-D1_3")

# Check new variable distribution - bulk density
#hist(soilbulk_full$BD) # no negative values anymore, quite nice normal distribution -> validated

# Check new variable distribution - soil moisture in percent weight
#hist(soilbulk_full$Weightpercent_Soilmoisture) # still some negative and extreme values (100%)
#filter(soilbulk_full, Weightpercent_Soilmoisture<20) # 5 cores with less than 20% soil moisture
# 3 cores from OC2, concerned with scale issue. OC2-P2-D2_2 negative value, OC2-P1-D1_3 very low not coherent with other samples from the plot -> to be removed - OC2-P3-D3_3 just under 20, not extreme compared with the other samples -> to be kept
# OG4-P3-D3_1, concerned with scale issue -> values are coherent within the plot and relatively close to what is find in other plots (10%-30%) -> to be kept
# soilbulk_full <- subset(soilbulk_full, BDcoreID != "OC2-P2-D2_2")
# soilbulk_full <- subset(soilbulk_full, BDcoreID != "OC2-P1-D1_3")
#filter(soilbulk_full, Weightpercent_Soilmoisture>90) # IS3-P3-D4_1, with P3 concerned with scale issue. Incoherent with other samples from same plot -> to be removed
# soilbulk_full <- subset(soilbulk_full, Weightpercent_Soilmoisture<90)

# Check new variable distribution - soil moisture in percent volume
#hist(soilbulk_full$Volpercent_Soilmoisture) # no extreme. nice normal distribution

# Check new variable distribution - standardised gravimetric water content
#hist(soilbulk_full$GWC_48) # no extreme. nice normal distribution

# Check new variable distribution - standardised volumetric water content
#hist(soilbulk_full$VWC_48) # no extreme. nice normal distribution

# Check new variable distribution - WFPS
#hist(soilbulk_full$percent_WFPS) # no extremes, nice normal distribution

# Check new number of replicates per site
#sort(table(soilbulk_full$SiteID)) 
# 9 sites with less than 20 replicates and lowest IC3 with 9 replicates (due to missing values) -> validated

#### DESTRUCTIVE SUBPLOTS - SOIL CHEMISTRY ####

### NA check
# colnames(soil_chem)[apply(soil_chem, 2, anyNA)] #validated

### Variable distribution & outliers
# hist(soil_chem$lossOnIgnitionPercentDM) # 

# Soil density
#soilchem_full[is.na(soilchem_full$SoilDensity_kg.L),] # no NA
#hist(soilchem_full$SoilDensity_kg.L) # range from 0 to 1.4 -> quite wide, but include both heathland and grassland. No visible outlier. Distribution a bit hectic

# Percent of humus in dry matter
#soilchem_full[is.na(soilchem_full$Humus_percentDM),] # no NA
#hist(soilchem_full$Humus_percentDM) # range from 0 to 90 -> matching with LOI

# pH
#soilchem_full[is.na(soilchem_full$pH),] # no NA
#hist(soilchem_full$pH) # range from 4 to 7, Normal distribution -> one outlier over 6.5
#filter(soilchem_full, pH>6.5) # OC2-P1 with the calcium outlier -> should be removed

# Phosphorus
#soilchem_full[is.na(soilchem_full$P.Al_mg.100g),] # no NA
#hist(soilchem_full$P.Al_mg.100g) # range from 0 to 40, Poisson distribution -> check high values
#filter(soilchem_full, P.Al_mg.100g>20) # 10 plots among 5 sites over 20 mg/100g
# 2 sites with all values over 20 (IS4, OC2)
# other sites (IC2, OC5, OG4), values not to far from other plots

# Potassium
#soilchem_full[is.na(soilchem_full$K.Al_mg.100g),] # no NA
#hist(soilchem_full$K.Al_mg.100g) # range from 0 to 30, Normal distribution -> check high values
#filter(soilchem_full, K.Al_mg.100g>20) # 4 plots among 2 sites (IC2, OC1) over 20 mg/100g -> coherent with other values

# Magnesium
#soilchem_full[is.na(soilchem_full$Mg.Al_mg.100g),] # no NA
#hist(soilchem_full$Mg.Al_mg.100g) # range from 0 to 35, Normal distribution -> check high values
#filter(soilchem_full, Mg.Al_mg.100g>20) # 4 plots among 2 sites (IC2, OC1) over 20 mg/100g, same as for Potassium

# Calcium
#soilchem_full[is.na(soilchem_full$Ca.Al_mg.100g),] # no NA
#hist(soilchem_full$Ca.Al_mg.100g) # range from 0 to 1000, one clear outlier
#filter(soilchem_full, Ca.Al_mg.100g>1000) # OC2-P1, not coherent with other samples -> to be removed
#filter(soilchem_full, Ca.Al_mg.100g>200) # 2 plots from same site (IC5) over 200 mg/100g

# Sodium
#soilchem_full[is.na(soilchem_full$Na.Al_mg.100g),] # no NA
#hist(soilchem_full$Na.Al_mg.100g) # range from 0 to 21, one clear outlier over 20
#filter(soilchem_full, Na.Al_mg.100g>20) # OC2-P1, same as Calcium -> to be removed
#filter(soilchem_full, Na.Al_mg.100g>12) # IS4-P3 & IS5-P2 -> coherent with rest of the samples

# Percent Dry Matter
#soilchem_full[is.na(soilchem_full$DryMatter_percent),] # no NA
#hist(soilchem_full$DryMatter_percent) # range from 10 to 100, distribution a bit hectic

# Total N in percent dry matter
#soilchem_full[is.na(soilchem_full$TotalN_percentDM),] # no NA
#hist(soilchem_full$TotalN_percentDM) # range from 0.2 to 2.2, distribution a bit hectic

#
## Data filtering/removal

# OC2-P1 outlier in several parameter -> farmer fertilizes in spring and summmer, maybe samples taken on a chunk
# soilchem_full <- subset(soilchem_full, !PlotID == "OC2-P1")

#### DESTRUCTIVE SUBPLOTS - SOIL MESOFAUNA ####

# test <- soil_meso |>  
#   summarise(
#     tibble(
#       across(
#         where(is.numeric),
#         ~min(.x, na.rm = TRUE),
#         .names = "min_{.col}"
#       ),
#       across(
#         where(is.numeric),
#         ~max(.x, na.rm = TRUE),
#         .names = "max_{.col}")
#     )
#   ) |>  
#   transpose() # All good

# Check distribution of quantitative variable
#hist(mesobio_raw$Acari) # Poisson distribution
#hist(mesobio_raw$Collembola) # Poisson distribution
#hist(soilmeso_raw$CoreDepth_cm) # Most around 14 cm -> validated

#
## New variable - abundance per soil area with correction soil volume

# Extraction survey data only
mesobio_full <- mesobio_raw |> 
  filter(SiteID != "ØY-")

# Merging datasets according to sorted fauna
mesobio_full <- left_join(mesobio_full, soilmeso_raw)
mesobio_full <- subset(mesobio_full, select = c(SampleID, Acari, Collembola, PlotID, SiteID, CoreDepth_cm))

# New variable with correction for soil volume
mesobio_full <- mesobio_full |> 
  # corrected abundance = (measured_abundance*standard_coreheight)/measured_coreheight
  mutate(Acari.m2 = ((Acari*14)/CoreDepth_cm)/(3.14*(0.105/2)^2)) |> 
  mutate(Collembola.m2 = ((Collembola*14)/CoreDepth_cm)/(3.14*(0.105/2)^2))
