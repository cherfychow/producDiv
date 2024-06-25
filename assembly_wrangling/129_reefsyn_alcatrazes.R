

# Curation Script ---------------------------------------------------------

# Dataset: ReefSYN
# Location: Alcatrazes
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(readr)
require(sf)
require(lubridate)

rm(list=ls())
dataID <- '129_reefsyn_alcatrazes'

# make sure your working directory is set before running these lines
dt <- read_delim('data_originals/129_ReefSYN/alcatrazes/occurrence.txt')
events <- read_delim('data_originals/129_ReefSYN/alcatrazes/event.txt')
measurements <- read_delim('data_originals/129_ReefSYN/alcatrazes/measurementorfacts.txt')

unique(measurements$measurementType)
length <- measurements %>% filter(measurementType == "total length") %>%
  select(measurementValue) %>% rename(Length = measurementValue) # TL in cm

length %>% filter(is.na(Length))
length %>% filter(Length == 0)

# check for method inconsistencies
colnames(events)
unique(events$samplingProtocol)
unique(events$sampleSizeUnit)
unique(events$sampleSizeValue) # uneven?
unique(events$samplingEffort)

abund  <- measurements %>% filter(measurementType == "abundance") %>% 
  select(measurementValue) %>% rename(Abundance = measurementValue)

# merge in abundances and lengths
colnames(dt)
dt <- dt %>% select(eventID, family, genus, scientificNameID, taxonRank) %>% 
  bind_cols(., length, abund)

# add in sampling event information
colnames(events)

dt <- events %>% select(eventID, locality, decimalLatitude, decimalLongitude,
  minimumDepthInMeters, maximumDepthInMeters, samplingEffort, eventDate) %>%
  left_join(dt, ., by = "eventID")

View(dt)
summary(dt) # a successful join would have given you no NAs
str(dt)
rm(length, measurements, events)

og_dt <- dt # backup copy

## Species name match ------------------------------------------------------------------
require(worrms)
ids <- unique(dt$scientificNameID) %>% str_extract(., '(?<=\\:)\\d+$') %>% as.numeric
ids <-  tibble(AphiaID = ids, scientificNameID = unique(dt$scientificNameID))
taxa <- wm_id2name_(ids$AphiaID)
ids$Species <- unlist(taxa)
dt <- ids %>% select(!AphiaID) %>% right_join(., dt, by = "scientificNameID")

# modify before doing checks
dt$eventID <- str_remove(dt$eventID, 'BR\\:ReefSYN\\:ALCATRAZES\\-TIME\\-SERIES\\:BrazilianCoast\\:alcatrazes\\_') # remove this, uninformative
dt <- dt %>% mutate(Depth = (minimumDepthInMeters + maximumDepthInMeters) / 2, # generate depth and date fields
  Year = year(eventDate), Month = month(eventDate), Day = day(eventDate)) %>% 
  select(!c(scientificNameID, minimumDepthInMeters, maximumDepthInMeters)) %>% # remove unnecessary
  rename(Latitude = decimalLatitude, Longitude = decimalLongitude) # renames

View(dt) #check
dt$eventDate <- NULL

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 7223 records
str(dt)
# Abundance and/or length, latitude and longitude numeric?
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)) == 0 # any NAs
sum(dt$Abundance == 0) == 0 # zeroes
sum(dt$Abundance == '') == 0 # blanks

sum(is.na(dt$Length)) == 0 # any NAs
sum(dt$Length == 0) == 0 # zeroes
sum(dt$Length == '') == 0 # blanks

summary(dt$Abundance)
summary(dt$Length) # range seems reasonable
dt %>% filter(Length > 100) %>% View # yep, all sharks and megafauna species

# YEAR MONTH DAY
dt %>% select(Year:Day) %>% summary

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

library(rnaturalearth)
library(ggplot2)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, Year), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160) # run twice and the error should be gone
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
dt$Species[which(dt$taxonRank == 'Genus')] <- 'sp'
colnames(dt) <- str_to_title(colnames(dt))

# NAs
dt %>% distinct(Family, Genus, Species) %>% filter(is.na(Family), is.na(Genus), is.na(Species))
unique(dt$Taxonrank)
unique(dt$Species) %>% sort
dt$Species[which(dt$Taxonrank == 'Species')] <- word(dt$Species[which(dt$Taxonrank == 'Species')], 2)
dt %>% distinct(Family, Genus, Species)

# Prepare raw data --------------------------------------------------------

View(dt)

dt_final <- dt %>% mutate(DataID = dataID) %>% 
  rename(SamplingEventID = Eventid) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
         Latitude, Longitude, Depth, Day, Month, Year, DataID) %>% 
  arrange(Year, Month, Day, SamplingEventID)

View(dt_final)

# Export final ------------------------------------------------------------

write_csv(dt_final, paste0('data_outputs/129_reefSyn/', dataID, '_obs.csv'))

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
