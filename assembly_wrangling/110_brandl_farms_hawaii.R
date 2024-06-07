

# Setup
require(tidyverse)
require(lubridate)
require(readxl)

rm(list=ls())
dataid <- '110_brandl_farms_hawaii'

# make sure your working directory is set before running these lines

dt <- read_excel('data_originals/110_brandl_farms/farms_alldata.xlsx', sheet = 2)
colnames(dt)
head(dt)

dt<- dt %>% select(!c(id, sl, mass:notes))

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 191 records
str(dt) # check structure
summary(dt)

dt <- dt %>% rename(Length = tl) %>%
  mutate(Abundance = 1, .before = Length, Latitude = 21.45, 
  Longitude = -157.81) # add approximate coordinates

dt$Length <- as.numeric(dt$Length)
summary(dt$Length)
# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. Y
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES/LENGTH
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Length))
dt  %>% filter(is.na(Length))
dt <- dt %>% filter(is.na(Length) == F) # omit these records
min(dt$Length, na.rm=T) # no zeroes? Y
dt$Length <- dt$Length / 10 # convert to cm

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
sum(is.na(dt$genspe)) == 0
sum(str_detect(dt$genspe, '^\\s+$')) == 0

sort(unique(dt$genspe))
dt$Genus <- word(dt$genspe, 1)
dt$Species <- word(dt$genspe, 2, -1)

# everything is ok!

# Prepare raw data --------------------------------------------------------
dt %>% group_by(site) %>% summarise(unitN = n_distinct(unit))
# mostly triplicate, but some losses according to paper

dt_final <- dt %>% group_by_at(vars(-Abundance)) %>% 
  summarise(Abundance = sum(Abundance)) %>% ungroup %>% #aggregate by length
  mutate(DataID = dataid,
         SamplingEventID = paste(site, unit, sep = "_"),
         Year = 2019, Month = NA, Day = NA, Family = NA, Depth = NA) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)

View(dt_final)

# Export final ------------------------------------------------------------

write_csv(dt_final, paste0('data_outputs/', dataid, '_obs.csv'))

# Convex Hull for centroid ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt %>% filter(is.na(Latitude) == F) %>% distinct(Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
