#### PACKAGES ####

library(tidyverse) #R language
library(readxl) #read xl files
library(lubridate) #standard date

#### FIELD MANAGEMENT ####

# Site management, from farmer interviews
landuse_raw <- read_excel(path = "data/rawdata/NBR_RawFarmerSurvey.xlsx", sheet="Farms Information_R")

## List of variables

# [1] Field identification code for data collection
# [2] Geographical location of the field, hamlet
# [3] Geographical location of the field, postcode
# [4] Geographical location of the field, municipality (both former and current classification)
# [5] Date of the interview with the farmer
# [6] Site productivity - infield high productivity, outfield low productivity
# [7] Main livestock, currently grazing in the field - ! rotational grazing management
# [8] Other livestock, grazing in other fields - ! rotational grazing management
# [9] Name of cow breed
# [10] Name of goat breed
# [11] Name of sheep breed
# [12] Number of adult animals in the main livestock
# [13] Number of young animals in the main livestock
# [14] Number of adult animals in other livestock
# [15] Number of young animals in other livestock
# [16] Size of the field containing the sampling area (ha) ! rotational grazing management
# [17] FROM GÅRDSKSART - Total grazing area of the farm (ha) ! rotational grazing management
# [18] Grazing density on the field - not collected, empty column
# [19] Period since the farmer has had the current livestock
# [20] Period without grazing management on the field
# [21] Number of months the main livestock grazes in the field in a year
# [22] Number of months other livestock graze in the field in a year ! rotational grazing
# [23] If the animals are kept inside or outside at night
# [24] Farmer impression of grazing pressure on the site
# [25] Former livestock which used to graze in the field (1)
# [26] Former livestock which used to graze in the field (2)
# [27] Former livestock which used to graze in the field (3)
# [28] Former livestock grazing period (1)
# [29] Former livestock grazing period (2)
# [30] Former livestock grazing period (3)
# [31] Type of farm management during the past 10 years (1)
# [32] Type of farm management during the past 10 years (2)
# [33] Type of farm management during the past 10 years (3)
# [34] If applicable, frequency of grass cutting
# [35] If applicable, frequency of mulching
# [36] If applicable, frequency of fertilization with manure
# [37] If applicable, type of manure used
# [38] If applicable, volume of manure used (m3)
# [39] If applicable, in which season the manure is used
# [40] If applicable, frequency of fertilization with artificial fertilizer
# [41] If applicable, weight of artificial fertilizer used (kg)
# [42] If applicable, in which season the artificial fertilizer is used
# [43] If applicable, frequency of fertilization with shell-sand/lime
# [44] If applicable, weight of shell-sand/limer used (kg)
# [45] If applicable, in which season the shell-sand/lime is used
# [46] If applicable, last time the field was sowed
# [47] If applicable, last time the field was plowed
# [48] If applicable, last time the field was drained
# [49] Type of land use prior to grazing
# [50] If applicable, type of production prior to grazing (1)
# [51] If applicable, type of production prior to grazing (2)
# [52] Time period of land use type prior to grazing
# [53] Comments

#
## Summary land use - Check table size, list of variables, variable types (num/chr)

#str(landuse_raw) # All good

#
## Name & character cleaning land use

# R friendly variable names
names(landuse_raw) <-  gsub(" ", "", names(landuse_raw))
names(landuse_raw) <-  gsub("\\)", "", names(landuse_raw))
names(landuse_raw) <-  gsub("\\(", "_", names(landuse_raw))
names(landuse_raw) <-  gsub("/", "_", names(landuse_raw))
names(landuse_raw) <-  gsub(";", "_", names(landuse_raw))
names(landuse_raw) <-  gsub("-", "_", names(landuse_raw))
names(landuse_raw) <-  gsub("\\?", "", names(landuse_raw))
names(landuse_raw) <-  gsub(":", "", names(landuse_raw))
names(landuse_raw) <- gsub("Sitecode", "SiteID", names(landuse_raw)) #rename siteID so it matches with other sheets
names(landuse_raw) <- gsub("Interviewdate", "InterviewDate", names(landuse_raw))
names(landuse_raw) <- gsub("Infield_outfield", "FieldType", names(landuse_raw))
names(landuse_raw) <- gsub("Typeoflivestock", "Livestock", names(landuse_raw))
names(landuse_raw) <- gsub("_ifapplicable", "", names(landuse_raw))
names(landuse_raw) <- gsub("_villsaubreed", "breed", names(landuse_raw))
names(landuse_raw) <- gsub("Numberofanimals", "FlockSize", names(landuse_raw))
names(landuse_raw) <- gsub("Surveysitegrazingsurface_ha", "SelectedFieldArea_ha", names(landuse_raw))
names(landuse_raw) <- gsub("Inmarksbeite_Gardskart", "FarmInfieldArea_ha", names(landuse_raw))
names(landuse_raw) <- gsub("Currentlivestock_fromdate_year", "LivestockSinceYear", names(landuse_raw))
names(landuse_raw) <- gsub("Gap_swithoutgrazing_timeperiod", "NoGrazingPeriod", names(landuse_raw))
names(landuse_raw) <- gsub("Yearlygrazingtimesurvey_year1_currentlivestock_monthsperyear", "YearlyGrazing1_month", names(landuse_raw))
names(landuse_raw) <- gsub("Yearlygrazingtimesurvey_year2_currentlivestock_monthsperyear", "YearlyGrazing2_month", names(landuse_raw))
names(landuse_raw) <- gsub("Impressionofgrazingpressureonsurveysite", "GrazingPressureAccordingToFarmer", names(landuse_raw))
names(landuse_raw) <- gsub("_type", "", names(landuse_raw))
names(landuse_raw) <- gsub("_Ifseveral", "", names(landuse_raw))
names(landuse_raw) <- gsub("Type_soffarmmanagement_soiltreatmentusedthelast10years", "FieldManagement1", names(landuse_raw))
names(landuse_raw) <- gsub("1...32", "2", names(landuse_raw))
names(landuse_raw) <- gsub("1...33", "3", names(landuse_raw))
names(landuse_raw) <- gsub("frequencyoffertilizing", "_freq", names(landuse_raw))
names(landuse_raw) <- gsub("frequency", "_freq", names(landuse_raw))
names(landuse_raw) <- gsub("type", "_type", names(landuse_raw))
names(landuse_raw) <- gsub("amount_vol_m3", "_volm3", names(landuse_raw))
names(landuse_raw) <- gsub("Artificialfertilizer", "ArtificialFert", names(landuse_raw))
names(landuse_raw) <- gsub("Artificialfertiliser", "ArtificialFert", names(landuse_raw))
names(landuse_raw) <- gsub("Art.fertilizer", "ArtificialFert", names(landuse_raw))
names(landuse_raw) <- gsub("amount_mass_kg", "_masskg", names(landuse_raw))
names(landuse_raw) <- gsub("season", "_season", names(landuse_raw))
names(landuse_raw) <- gsub("Shell_sand_lime", "ShellSandLime", names(landuse_raw))
names(landuse_raw) <- gsub("...53", "", names(landuse_raw))
names(landuse_raw) <- gsub("...54", "2", names(landuse_raw))
names(landuse_raw) <- gsub("Comments", "Comments_landuse", names(landuse_raw))

# Removal empty columns
#landuse_raw <- subset(landuse_raw, select = -c(Grazingdensity_perha)) # will be calculated later in the script

#
## Data cleaning - new R object

landuse_full <- subset(landuse_raw, select = -c(Locality, Postcode, Grazingdensity_perha, GrazingPressureAccordingToFarmer, Inside_outsideatnight))

#
## Char var land use - Check if all sites/samples are present, categories, doubletons, NAs, misprints...

# Site ID
#table(landuse_full$SiteID) # Unique ID for each site - validated

# Livestock1
#unique(landuse_full$Livestock1) # villsau should be in sheep category
landuse_full <- landuse_full |> 
  mutate(Livestock1 = dplyr::recode(Livestock1, "villsau" = "sheep")) |> 
  mutate(Livestock1 = dplyr::recode(Livestock1, "cow" = "cattle"))
#table(landuse_full$Livestock1) # correct number of sites per livestock category - validated

# Livestock2
#unique(landuse_full$Livestock2) # villsau should be in sheep category
landuse_full <- landuse_full |> 
  mutate(Livestock2 = dplyr::recode(Livestock2, "villsau" = "sheep")) |> 
  mutate(Livestock1 = dplyr::recode(Livestock1, "cow" = "cattle"))
#table(landuse_full$Livestock2) # Only 5 farmers with another livestock -> not to be included in the analysis

# Period without grazing
#unique(landuse_full$NoGrazing_period) # NAs ?
#landuse_full[is.na(landuse_full$NoGrazing_period),] # 8 NAs, including the 6 farmers who did not respond to the survey
#table(landuse_full$NoGrazing_period) # Only 7 sites with potential grazing interruption, in a quite hectic way (some 2-3 years, other 10 years), 30 sites with no interruption -> not to be included into the analysis, but interesting for the material and methods
#landuse_full[landuse_full$NoGrazing_period != "no",] # Interruption are on OS1, IG1, IS2, OS5, OS7, OC2 and OS3

# If the animals are kept inside or outside at night
#unique(landuse_full$Inside_outsideatnight) # some NAs
#landuse_full[is.na(landuse_full$Inside_outsideatnight),] #18 sites with missing values -> not to be included in the analysis

# Farmer impression of grazing pressure on the site
#unique(landuse_full$Impressionofgrazingpressureonsurveysite) # NAs
#landuse_full[is.na(landuse_full$Impressionofgrazingpressureonsurveysite),] # 31 missing values -> not to be included in the analysis

# Former livestock (1), before the current livestock
#unique(landuse_full$Formerlivestock1) # NAs + villsau should be sheep category
landuse_full <- landuse_full |> 
  mutate(Formerlivestock1 = dplyr::recode(Formerlivestock1, "villsau" = "sheep"))
#table(landuse_full$Formerlivestock1) # 11 former cows, 13 former sheep, 2 former horse and 1 "other"
#landuse_full[is.na(landuse_full$Formerlivestock1),] # 18 missing values, including the 6 farmers which did not reply the survey -> 12 with no former livestock known -> could be used in the analysis?

# Former livestock (2), before livestock (1)
#unique(landuse_full$Formerlivestock2) # NAs + villsau should be sheep category
landuse_full <- landuse_full |> 
  mutate(Formerlivestock2 = dplyr::recode(Formerlivestock2, "villsau" = "sheep"))
#table(landuse_full$Formerlivestock2) # 5 former cows, 2 former sheep, 2 former horse
#landuse_full[is.na(landuse_full$Formerlivestock2),] # 36 missing values -> not to be included in the analysis

# Former livestock (3), before livestock (2)
#unique(landuse_full$Formerlivestock3) # NAs
#table(landuse_full$Formerlivestock3) # 1 former cows, 1 other
#landuse_full[is.na(landuse_full$Formerlivestock3),] # 43 missing values -> not to be included in the analysis

# Former livestock (1) grazing period
#unique(landuse_full$Formerlivestock1_timeperiod) # quite hectic, need to be transformed into how many years the former livestock has grazed -> not to be included in the analysis

# Management or soil treatment used in the past 10 years
unique(landuse_full$FieldManagement1) # NAs
unique(landuse_full$FieldManagement2)
unique(landuse_full$FieldManagement3)
#landuse_full[is.na(landuse_full$FieldManagement1),] # 6 missing values -> farmers who did not reply to the survey
#landuse_full[is.na(landuse_full$FieldManagement2),] # 36 missing values
#landuse_full[is.na(landuse_full$FieldManagement3),] # 41 missing values
#table(landuse_full$FieldManagement1) # 9 sites with no treatments at all, 21 with at least fertilization
#table(landuse_full$FieldManagement2) # 8 sites with at least 2 treatments
#table(landuse_full$FieldManagement3) # 3 sites with 3 treatments

# Sites with no treatment
#subset(landuse_full, FieldManagement1 == "none") # 4 mountain sites (US2, US3, US5, US6), 2 coastal heathlands (IS2, OS5), 3 grasslands (IG1, IG2, IS3)

# Sites with burning
#subset(landuse_full, FieldManagement1 == "burning") # 3 coastal heathlands only treatment (OV1, OV2, OS7) -> validated
#subset(landuse_full, FieldManagement3 == "burning") # 1 grassland (IS1), with also mulching and fertilisation -> validated

# Sites with tree cutting
#subset(landuse_full, FieldManagement1 == "sitka spruce removal") # 1 coastal heathland only treatment (OS9) -> validated
#subset(landuse_full, FieldManagement1 == "tree cutting") # 1 mountain site only treatment (UG1), 1 mountain site with also grass cutting (UG2) -> validated

# Sites with herbicide
#subset(landuse_full, FieldManagement1 == "herbicides") # 1 grassland (OS5) with also tree cutting as second treatment
#subset(landuse_full, FieldManagement3 == "herbicides") # 1 grassland (OC2)

# Sites with mulching
#subset(landuse_full, FieldManagement1 == "mulching") # 2 grassland (OS3, OS1) only treatment -> validated
#subset(landuse_full, FieldManagement2 == "mulching") # 6 grasslands (OC1, IS1, IV1, IC5, OS4, OC2) with mulching + fertilization. OC2 also has herbicides.

# Frequency of mulching
#table(landuse_full$Mulching_freq) # 8 sites with mulching
#subset(landuse_full, Mulching_freq == "every year") # 3 sites with annual mulching (IS1, OS4, OS3) -> validated
#subset(landuse_full, Mulching_freq == "sometimes") # 5 sites with occasional mulching (OC1, OS1, IV1, IC5, OC2) -> validated

# Types of fertilization
#unique(landuse_full$Fertilizing_type1)
#unique(landuse_full$Fertilizing_type2)
#unique(landuse_full$Fertilizing_type3) # 3 types of fertilization: artificial fertilizer, organic (manure) and inorganic (shell-sand/lime)
#table(landuse_full$Fertilizing_type3) # 5 sites with 3 types of fertilization
#subset(landuse_full, Fertilizing_type3 != is.na(Fertilizing_type3)) # 5 grassland sites (IC1, OC3, IC4, IC5, OC2)
#table(landuse_full$Fertilizing_type2) # 11-5 = 6 sites with 2 types of fertilization
#filter(landuse_full, !is.na(Fertilizing_type2) & is.na(Fertilizing_type3)) # 6 grasslands (OC1, IS1, OC4, OG6, OC5, OS8)
#table(landuse_full$Fertilizing_type1) # 21-11 = 10 sites with only 1 type of fertilization
#filter(landuse_full, !is.na(Fertilizing_type1) & is.na(Fertilizing_type2)) # 10 grasslands (OG1, OS2, IV1, OG3, IC3, OG2, OG4, OS4, OS6, IS4)

# Frequency of manure
#unique(landuse_full$Manure_freq) # possibilities are "yearly" or "occasionally"
#subset(landuse_full, Manure_freq == "every year") # 9 grasslands (OC1, IC1, IS1, OC3, OC4, OG6, IC4, IC5, OC5) with annual manure
#subset(landuse_full, Manure_freq == "sometimes") # 2 grasslands (OC2, OS8) with occasional manure

# Frequency of artificial fertilization
#unique(landuse_full$ArtificialFert_freq) # possibilities are "yearly" or "occasionally"
#subset(landuse_full, ArtificialFert_freq == "every year") # 15 grasslands (OC1, OS2, IC1, OG3, IC3, OC3, OC4, OG2, OG4, IC4, IC5, OC5, OS6, OC2, IS4) with annual fertilization with artificial fertiliser
#subset(landuse_full, ArtificialFert_freq == "sometimes") # 4 grasslands (OC2, OG6, IV1, OS8) with occasional fertilization with artificial fertiliser

# Frequency of inorganic fertilization with shell-sand/lime
#unique(landuse_full$ShellSandLime_freq) # possibilities are "occasionally" only
#subset(landuse_full, ShellSandLime_freq == "sometimes") # 7 grassland (IC1, OG3, IC3, OC3, IC4, IC5, OC2)

#
## New variables for field management/treatments for cleaner and standardised output

# Herbicide
landuse_full <- landuse_full |> 
  mutate(FM_herbicide = ifelse(
    FieldManagement1 == "herbicides" |
      FieldManagement3 == "herbicides",
    "yes",
    "no"))

# Tree cutting
landuse_full <- landuse_full |> 
  mutate(FM_treecutting = ifelse(
    FieldManagement1 == "tree cutting" |
      FieldManagement1 == "sitka spruce removal" |
      FieldManagement2 == "tree cutting" |
      FieldManagement3 == "tree cutting",
    "yes",
    "no"))

# Burning
landuse_full <- landuse_full |> 
  mutate(FM_burning = ifelse(
    FieldManagement1 == "burning" |
      FieldManagement3 == "burning",
    "yes",
    "no"))

# Mowing
landuse_full <- landuse_full |> 
  mutate(FM_mowing = ifelse(FieldManagement2 == "grass cutting", "yes", "no"))

# Mulching
landuse_full <- landuse_full |> 
  mutate(FM_mulching = ifelse(
    Mulching_freq == "every year", "yearly", ifelse(
      Mulching_freq == "sometimes", "occasionally", "no"
    )))

# Manure
landuse_full <- landuse_full |> 
  mutate(FM_manure = ifelse(
    Manure_freq == "every year", "yearly", ifelse(
      Manure_freq == "sometimes", "occasionally", "no"
    )))

# Artificial fertilizer
landuse_full <- landuse_full |> 
  mutate(FM_artificialfert = ifelse(
    ArtificialFert_freq == "every year", "yearly", ifelse(
      ArtificialFert_freq == "sometimes", "occasionally", "no"
    )))

# Shell-sand-lime
landuse_full <- landuse_full |> 
  mutate(FM_shellsandlime = ifelse(
    ShellSandLime_freq == "every year", "yearly", ifelse(
      ShellSandLime_freq == "sometimes", "occasionally", "no"
    )))

# Replace NAs by "no" in summary management columns only
landuse_full <- landuse_full |> 
  mutate_at(vars(starts_with("FM_")), ~replace(., is.na(.), "no"))

#
## Numeric var land use - Check min/max, distribution and potential outliers

# Check min/max
test <- landuse_full |>  
  summarise(
    tibble(
      across(
        where(is.numeric),
        ~min(.x, na.rm = TRUE),
        .names = "min_{.col}"
      ),
      across(
        where(is.numeric),
        ~max(.x, na.rm = TRUE),
        .names = "max_{.col}")
    )
  ) |>  
  transpose() # All good

# Flock size adults - missing values
#landuse_full[is.na(landuse_full$FlockSize1_adults),] # some farmers would not or could not respond to our interview. The missing value were therefore taken from another series of interviews made in 2018
landuse_full <- landuse_full |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "IC2", 37, FlockSize1_adults)) |> 
  mutate(Cowbreed = ifelse(SiteID == "IC2", "norsk rødt fe", Cowbreed)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "UC1", 13, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "US1", 4, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "US4", 66, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "IS5", 90, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "OS7", 39, FlockSize1_adults)) |> 
  mutate(FlockSize1_adults = ifelse(SiteID == "IG3", 109, FlockSize1_adults)) # no information for US6

# Flock size adults - distribution
#hist(landuse_full$FlockSize1_adults) # wide range but most farms under 50 animals - one farm over 150 animals
#landuse_full[landuse_full$FlockSize1_adults>150,] # IS4 over 150 animals + US6 as NA

# Flock size young - missing values
landuse_full <- landuse_full |> 
  mutate(FlockSize1_young = ifelse(SiteID == "IC2", 17, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "UC1", 13, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "US1", 10, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "US4", 63, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "IS5", 181, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "OS7", 14, FlockSize1_young)) |> 
  mutate(FlockSize1_young = ifelse(SiteID == "IG3", 135, FlockSize1_young)) # no information for US6

# Flock size young - distribution
#hist(landuse_full$FlockSize1_young)
#landuse_full[landuse_full$FlockSize1_young>150,] # US3 & IS5 + 4 farms with no young animals + US6 as NA

# Grazing surface
#landuse_full[is.na(landuse_full$GrazingSurface_ha),] # No missing value - validated
#hist(landuse_full$GrazingSurface_ha) # big range, one outlier over 1000 ha
#landuse_full[landuse_full$GrazingSurface_ha>1000,] # UG2 in the mountain -> extended area which serves for several farmers

# Total infield surface - missing values & distribution
#landuse_full[is.na(landuse_full$TotalInfieldSurface),] # All missing values should be outfields/heathland sites - validated
#hist(landuse_full$FarmInfieldArea_ha) # big range but no visible outlier

# Date since the current livestock has been grazing
#landuse_full[is.na(landuse_full$LivestockFrom_year),] # 7 missing values from farmers who did not reply the survey
#hist(landuse_full$LivestockFrom_year) # ranges from 1920 to 2020, but ! not equivalent to no grazing, neither to grazing with another animal -> too hectic, not to be included in the analysis

# How many months main livestock graze on site during the year
#landuse_full[is.na(landuse_full$YearlyGrazing1_month),] # 6 missing values from farmers who did not reply the survey
#hist(landuse_full$YearlyGrazing1_month) # coherent values (between 2 and 12), no visible outliers -> validated

#
## New variable - average stocking density

# Livestock unit for adults - depends on type of livestock and breed (for cows)
#unique(landuse_full$Cowbreed) # Milking cows (norsk rødt fe) are 1 SSU - Beef cattle (Limousin, Aberdeen angus, Highland) are 0.8 LSU - Sheep and goats are 0.1 LSu
landuse_full <- landuse_full |> 
  mutate(LSU_adult = ifelse(Livestock1 == "sheep" | Livestock1 == "goat", 0.1, ifelse(Cowbreed == "norsk rødt fe", 1, 0.8)))
#hist(landuse_full$LSU_adultind) # validated

# Livestock unit for young - half for goats and sheep (0.05), young cows on field (0.7)
landuse_full <- landuse_full |> 
  mutate(LSU_young = ifelse(Livestock1 == "sheep" | Livestock1 == "goat", 0.05, 0.7))
#hist(landuse_full$LSU_youngind) # validated

# Grazing intensity over the year - two goat flocks in Osterøy were known to be moved to summer farm (outfield)
landuse_full <- landuse_full |>  
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) |> 
  mutate(AvgStockDensity_perha = ifelse(
    Livestock1 == "goat" & Municipality_old == "Osterøy", 
    ((FlockSize1_adults*LSU_adult+FlockSize1_young*LSU_young)/FarmInfieldArea_ha)*(YearlyGrazing1_month/12), 
    (FlockSize1_adults*LSU_adult+FlockSize1_young*LSU_young)/FarmInfieldArea_ha)
  )
#hist(landuse_full$Grazingdensity_perha) # Range from 0 to 0.6, no visible outlier
#landuse_full[is.na(landuse_full$Grazingdensity_perha),] # US6 as NA

#
## Export clean data in new excel file

write_csv(landuse_full, "data/cleandata/NBR_FullLanduse.csv")



