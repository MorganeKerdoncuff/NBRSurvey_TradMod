# TradMod WP2 - NBR survey cleaning script
#Description of the data
#Date
#Who
#Project
#Funding
#Place

####Package loading####

library(tidyverse) #R language
library(readxl) #read xl files
library(lubridate) #standard date data

####Raw data loading####

siteinfo_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="SiteInfo") #Farm information, field location and habitat type

landuse_raw <- read_excel(path = "data/Farmer_Survey.xlsx", sheet="Farms Information_R") #Field management data from farmer interview

landscape_raw <- read_excel(path = "data/NBR_LandscapeMatrix.xlsx", sheet="MatrixProportion") #Land cover data around the fields from Geonorge

site20x20_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="20mX20m") #Sampling area description, vegetation cover at the site level

vgcover_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="SoilCover") #Vegetation cover at the quadrat (subplot) level

soilpene_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="SoilPenetration") #Penetration rate in the soil

soilbulk_raw <- read_excel(path = "data/NBR_RawBD2019-2020.xlsx", na="NA") #Soil bulk density

soilmeso_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="Mesofauna") #Height of mesofauna soil core

chem2019 <- read.csv(path = "data/SoilChemistry_NBR2019.txt", sep=";") #2019 Soil chemistry data from Eurofins

chem2020 <- read.csv(path = "data/SoilChemistry_NBR2020.txt", sep=";") #2020 soil chemistry data from EUrofins

poo_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="Poo", na="NA") #poo data at subplot (quadrat) level

vege_raw <- read_excel(path = "data/NBR_RawAll.xlsx", sheet="PlantRichness") #plant community data, at species level

arthro_raw <- read_excel(path = "data/NBR_RawArthro2019-2020.xlsx", na="NA") #arthropod community data, at family level

