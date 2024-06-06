

# Curation Script ---------------------------------------------------------

# Dataset: NOAA coral reef monitoring fish data from sites across Micronesia from 2009-09-29 to 2015-09-25
# Location: Micronesia
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(sf)
require(world)
require(ggplot2)

rm(list=ls())
dataid = '127_philippines_liloan'

# make sure your working directory is set before running these lines
dt <- read.delim('data_originals/127_philippines_liloan/occurrence.txt', header=T, sep = "\t")
dt <- dt %>% filter(class == 'Actinopterygii') # fish only

# Structure organisation + check --------------------------------------------------

# extract lengths
str_extract(dt$organismQuantityType, '^[0-9]+(?=cm)')
dt$Length <- str_extract(dt$organismQuantityType, '^[0-9]+(?=cm)') %>% as.numeric
summary(dt$Length)

require(lubridate)
unique(dt$eventDate)
dt$eventDate <- parse_date_time(dt$eventDate, "Ym")
dt$Year <- year(dt$eventDate)
dt$Month <- month(dt$eventDate)
dt$eventDate <- NULL

dt <- dt %>% select(individualCount, organismQuantity, Length, family, genus, scientificName, eventID, decimalLatitude, decimalLongitude, Year, Month)

colnames(dt)[c(1,2,7:9)] <- c('Abundance', 'Abundance', 'SamplingEventID', 'Latitude', 'Longitude')

nrow(dt) # check dimensions

# Primary field check -----------------------------------------------------

summary(dt)
dt <- dt[-2]

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Length)+0)
min(dt$Length, na.rm=T) # no zeroes? Y

sum(is.na(dt$Abundance)+0)
min(dt$Abundance, na.rm=T) # no zeroes? Y

# YEAR MONTH DAY
# again, no negative values, 0s or NAs, and all are logical
# Year < 2020, month < 12, day < 31

# LAT LONG
sum(dt[c('Latitude', 'Longitude')] == "") == 0 # no blanks?
sum(is.na(dt[c('Latitude', 'Longitude')])) == 0 # no NAs?

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good


# Sampling effort + coord checks --------------------------------------------------

require(ggplot2)
world_map <- map_data('world') # check whether the GPS coordinates match expectations

ggplot(world_map) + 
  geom_polygon(aes(x=long, y=lat, group=group), color='grey', fill='transparent') +
  coord_fixed() + 
  labs(x='Longitude', y='Latitude') +
  theme_bw() + theme(panel.grid=element_blank()) + 
  geom_point(data=dt %>% distinct(Longitude, Latitude), 
                   aes(x=Longitude, y=Latitude), shape=21, size = 2) + 
  coord_fixed(xlim=c(min(dt$Longitude) - 10, max(dt$Longitude) + 10), ylim=c(min(dt$Latitude) - 10, max(dt$Latitude) + 10))

# check sampling space x time
# just so there is not a lone color by itself
ggplot(world_map) + 
  geom_polygon(aes(x=long, y=lat, group=group), color='grey', fill='transparent') +
  coord_fixed() + 
  labs(x='Longitude', y='Latitude') +
  theme_bw() + theme(panel.grid=element_blank()) + 
  geom_point(data=dt, aes(x=Longitude, y=Latitude, color = Year), size = 2, alpha = 0.5) + 
  coord_fixed(xlim=c(min(dt$Longitude) - 10, max(dt$Longitude) + 10), ylim=c(min(dt$Latitude) - 10, max(dt$Latitude) + 10))
rm('world_map')

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
length(which(dt$Family == '' & dt$Genus == '' & dt$Species == ''))
dt %>% distinct(family, genus, scientificName) %>% View

dt$Species <- str_replace(dt$scientificName, 'sp\\.$', 'sp') %>% word(., 2)
dt$scientificName <- NULL

# Final restructuring ---------------------------------------------------

colnames(dt)[-5] <- str_to_title(colnames(dt)[-5])
dt_final <- dt %>% mutate(Day = NA, Depth = NA) %>% arrange(Year, Month, SamplingEventID, Family, Genus, Species) %>% 
  select(Abundance, Length, Family,	Genus, Species, SamplingEventID, Latitude, Longitude, Depth, Day, Month, Year)


# Export final ------------------------------------------------------------

readr::write_csv(dt_final, paste0('data_outputs/', dataid, '_obs.csv'))

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# MPA
dt_coord <- dt_final %>% distinct(Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long

