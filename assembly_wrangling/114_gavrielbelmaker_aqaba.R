

# Curation Script ---------------------------------------------------------

# Dataset: Reef fish at Gulf of Aqaba
# Location: Aqaba, Red Sea
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(readxl)
require(sf)
require(tidyverse)
require(lubridate)

rm(list=ls())

# already did formatting so this is just a cursory check for formatting etc.

# make sure your working directory is set before running these lines
dt <- read_excel('data_originals/114_gavrielbelmaker_aqaba/114_GavrielBelmaker_Aqaba_obs.xlsx', sheet = 2)

# Structure organisation + check --------------------------------------------------

summary(dt)
str(dt)
nrow(dt) # check dimensions

dt$StudyID <- '114_GavrielBelmaker_Aqaba'
dt %>% filter(is.na(Longitude)) %>% View
dt$SampleDescription <- paste(dt$Plot, dt$SampleDescription, sep = "_")
dt$Plot <- NULL
dt$Longitude[which(is.na(dt$Longitude))] <- mean(dt$Longitude, na.rm = T)

summary(dt)
str(dt)

# Secondary fields relating to sample must be factors or integers? Y
# Taxonomic fields must be characters or factors? Y


# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields. Y
dt %>% filter(Abundance == 0 | Length == 0) %>% View
dt <- dt %>% filter(Abundance > 0, Length > 0) # omit these zero observation records

summary(dt)

# YEAR MONTH DAY, all sensical

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid()
world <- st_break_antimeridian(world, lon_0 = 160)
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) + # confine to trpical view
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)

# Prepare raw data --------------------------------------------------------
dt_final <- dt %>% 
  rename(SamplingEventID = SampleDescription, DataID = StudyID, Depth = DepthElevation) %>% 
  arrange(Year, Month, Day, SamplingEventID) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID) 

# Export final ------------------------------------------------------------

readr::write_csv(dt_final, 'data_outputs/114_gavrielbelmaker_aqaba_obs.csv')

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
