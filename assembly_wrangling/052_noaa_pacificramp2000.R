

# Curation Script ---------------------------------------------------------

# Dataset: Pacific Reef Assessment and Monitoring Program: Rapid Ecological Assessments of Fish Large-Area Stationary Point Count Surveys (SPC) at Coral Reef Sites across the Pacific Ocean from 2000 to 2007
# Location: Pacific Islands
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(tidyverse)
require(sf)

rm(list=ls())
dataID <- '052_noaa_pacificramp2000'

# make sure your working directory is set before running these lines
dt <- read_delim('data_originals/052_noaa_pacificramp2000/occurrence.txt')
og_dt <- dt
events <- read_delim('data_originals/052_noaa_pacificramp2000/event.txt')
measurements <- read_delim('data_originals/052_noaa_pacificramp2000/extendedmeasurementorfact.txt')

unique(measurements$measurementType)
length <- measurements %>% filter(measurementType == "Fish Total Length, TL")
length <- length %>% select(id, occurrenceID, measurementValue) %>% rename(Length = measurementValue) # TL in cm

length %>% filter(is.na(Length))
length %>% filter(Length == 0)

# check for method inconsistencies
unique(events$samplingProtocol)
unique(events$habitat)
unique(events$dataGeneralizations)
unique(events$sampleSizeValue) # uneven?
events %>% group_by(sampleSizeValue) %>% summarise(nsurvey = n())
events %>% filter(sampleSizeValue < 10) %>% distinct(sampleSizeValue, eventDate) # someone messed up in 2005
# keep only 10 m radius

# merge in abundances and lengths
dt <- dt %>% select(occurrenceID, eventID, individualCount, family, genus, specificEpithet, 
                    decimalLatitude, decimalLongitude, eventDate) %>% 
  full_join(., length, by = "occurrenceID")
summary(dt) # check no NAs
# add in sampling event information
dt <- events %>% select(eventID, sampleSizeValue,islandGroup, island, stateProvince, 
                  minimumDepthInMeters, maximumDepthInMeters) %>% 
  left_join(dt, ., by = "eventID") %>% 
  filter(sampleSizeValue == 10) # and filter out the stationary counts that were outliers

View(dt)
summary(dt) # a successful join would have given you no NAs
str(dt)
rm(length, measurements, events)
nrow(og_dt) - nrow(dt) # dropped 738 records
og_dt <- dt # make this the original backup now

# modify before doing checks
dt$eventID <- str_remove(dt$eventID, '^PIFSC\\:ESD\\:') # remove this, uninformative
dt$stateProvince <- str_remove_all(dt$stateProvince, '\\s')
dt$island <- str_remove_all(dt$island, '\\s')
dt$SamplingEventID <- with(dt, paste(stateProvince, island, str_replace(eventID, '\\:', '\\_'), sep = '_'))
dt <- dt %>% mutate(Depth = (minimumDepthInMeters + maximumDepthInMeters) / 2) %>% 
  select(!c(id, occurrenceID, islandGroup, island, stateProvince, eventID, minimumDepthInMeters, maximumDepthInMeters)) # the region/site codes are already baked into the event ID

colnames(dt)[1:6] <- c('Abundance', 'Family', 'Genus', 'Species', 'Latitude', 'Longitude')

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 28384 records
str(dt)
# Abundance and/or length, latitude and longitude numeric? N
unique(dt$Length)
dt$Length <- as.numeric(dt$Length)
dt$sampleSizeValue <- NULL

# Dates
require(lubridate)
dt$Day <- day(dt$eventDate)
dt$Month <- month(dt$eventDate)
dt$Year <- year(dt$eventDate)
summary(dt)
dt$eventDate <- NULL

# Secondary fields relating to sample must be factors or integers? Y

# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)+0) # any NAs
dt %>% filter(is.na(Abundance)) %>% View # no records at all. omit
dt <- dt %>% filter(is.na(Abundance) == F & is.na(Length) == F)
min(dt$Abundance, na.rm=T) # no zeroes? Y

sum(is.na(dt$Length)+0) # no NAs
min(dt$Length, na.rm=T) # no zeroes? Y
summary(dt$Length) # range seems reasonable
dt %>% filter(Length > 100) %>% View # yep, all sharks and megafauna species

# YEAR MONTH DAY
# already checked

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, Year), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160)
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, aes(fill = Year), shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
dt %>% distinct(Family, Genus, Species) %>% filter(is.na(Family), is.na(Genus), is.na(Species))
dt %>% filter(is.na(Family)) %>% nrow
dt <- dt %>% filter(!is.na(Family))
dt %>% distinct(Family, Genus, Species) %>% filter(is.na(Genus)) %>% View
dt %>% distinct(Family, Genus, Species) %>% filter(is.na(Species)) %>% View
dt$Species[which(is.na(dt$Family) == F & is.na(dt$Genus) == F & is.na(dt$Species) == T)] <- 'sp'
sum(str_detect(na.omit(dt$Species), '^\\s+$') + 0) == 0

dt %>% distinct(Family, Genus, Species) %>% arrange(Genus, Species) %>% View
dt$Species[which(is.na(dt$Species) & is.na(dt$Genus) == F)] <- 'sp'

# Prepare raw data --------------------------------------------------------

dt %>% mutate(site = str_extract(SamplingEventID, '^[A-Za-z]+\\_[A-Za-z]+')) %>% 
  group_by(site) %>% summarise(years = n_distinct(Year)) %>% View

dt_final <- dt %>% mutate(DataID = dataID) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)

# Export final ------------------------------------------------------------

write_csv(dt_final, paste0('data_outputs/', dataID, '_obs.csv'))

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt_final %>% distinct(Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
