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
