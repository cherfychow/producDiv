

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

events <- events %>% filter(sampleSizeValue == 10)
dt <- dt %>% filter(eventID %in% events$eventID) # filter
# only dropped ~600 records

# merge in abundances and lengths
dt <- dt %>% select(occurrenceID, individualCount, eventID, family, genus, specificEpithet, decimalLatitude, decimalLongitude, eventDate) %>% 
  left_join(., length, by = "occurrenceID") %>% 
  left_join(., events[c(13,25, 26, 29, 31, 32)], by = "eventID")

View(dt)
summary(dt)
str(dt)

rm(length, measurements)

# modify before doing checks
dt <- dt %>% mutate(Depth = (minimumDepthInMeters +maximumDepthInMeters) / 2,
                    SamplingEventID = str_replace_all(eventID, '\\s\\:\\s', '_')) %>% 
  select(!c(id, occurrenceID, islandGroup, island, stateProvince, eventID, minimumDepthInMeters, maximumDepthInMeters)) # the region/site codes are already baked into the event ID

dt$SamplingEventID <- str_remove(dt$SamplingEventID, '(?<=(replicate))ID\\s')
colnames(dt)[3:7] <- c('Family', 'Genus', 'Species', 'Latitude', 'Longitude')

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 698933 records
str(dt)
# Abundance and/or biomass, latitude and longitude numeric? Y

# Dates
require(lubridate)
dt$Day <- day(dt$eventDate)
dt$Month <- month(dt$eventDate)
dt$Year <- year(dt$eventDate)
summary(dt)
dt %>% filter(is.na(eventDate)) %>% View # event date is in the event ID. all the same survey
dt$Day[which(is.na(dt$Day))] <- 24
dt$Month[which(is.na(dt$Month))] <- 10
dt$Year[which(is.na(dt$Year))] <- 2014
summary(dt)
dt$eventDate <- NULL

# Secondary fields relating to sample must be factors or integers?

# Taxonomic fields must be characters or factors? Y

# fix NA coordinates

dt %>% filter(is.na(Latitude)) %>% View
dt %>% filter(str_detect(SamplingEventID, 'GUA\\-1310')) %>% distinct(Latitude, Longitude)
events %>% filter(str_detect(eventID, 'GUA\\-1310')) %>% distinct(decimalLatitude, decimalLongitude)
# just zero coordinates!
# assign it the mean of proximate Guam sites..
dt %>% filter(str_detect(SamplingEventID, 'GUA\\-13')) %>% distinct(Latitude, Longitude, SamplingEventID) %>% View
dt$Latitude[which(is.na(dt$Latitude))] <- mean(c(13.24146, 13.24117))
dt$Longitude[which(is.na(dt$Longitude))] <- mean(c(144.7070, 144.7083))

summary(dt)

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)+0)
min(dt$Abundance, na.rm=T) # no zeroes? Y
dt %>% filter(Abundance == 0) %>% View
dt <- dt %>% filter(Abundance > 0)
sum(is.na(dt$Length)+0) # 237 NAs for length
min(dt$Length, na.rm=T) # no zeroes? Y
dt %>% filter(is.na(Length)) %>% View
dt <- dt %>% filter(is.na(Length) == F) # have to omit these records
dt %>% filter(Length > 100) %>% arrange(desc(Length)) %>% View # max 999 length??
dt$Length[which(dt$Length == 999)] <- 350 # replace this with common length in FishBase

# YEAR MONTH DAY
# already checked

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, Year), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid()
world <- st_break_antimeridian(world, lon_0 = 160)
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
dt <- dt %>% filter(is.na(Family) == F)# remove
dt %>% distinct(Family, Genus, Species) %>% filter(is.na(Genus)) %>% View

sum(str_detect(na.omit(dt$Species), '^\\s+$') + 0) == 0

dt %>% distinct(Family, Genus, Species) %>% arrange(Genus, Species) %>%  View
dt$Species[which(is.na(dt$Species) & is.na(dt$Genus) == F)] <- 'sp'

# Prepare raw data --------------------------------------------------------

dt %>% mutate(site = str_extract(SamplingEventID, '^[A-Z]+(?=\\-)')) %>% 
  group_by(site) %>% summarise(years = n_distinct(Year)) %>% View

dt_final <- dt %>% mutate(DataID = dataID) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)

# Export final ------------------------------------------------------------

write_csv(dt_final, paste0('data_outputs/',dataID, '/', dataID, '_obs.csv'))

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
