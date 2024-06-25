

# Curation Script ---------------------------------------------------------

# Dataset: ReefSYN
# Location: Guarapari
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(readr)
require(sf)
require(lubridate)
require(ggplot2)
require(rnaturalearth)

rm(list=ls())
dataID <- '129_reefsyn_trindade'

# make sure your working directory is set before running these lines
dt <- read_delim('data_originals/129_reefsyn/trindade/occurrence.txt')
events <- read_delim('data_originals/129_reefsyn/trindade/event.txt')
measurements <- read_delim('data_originals/129_reefsyn/trindade/measurementorfacts.txt')

unique(measurements$measurementType)
length <- measurements %>% filter(measurementType == "total length") %>%
  select(measurementValue) %>% rename(Length = measurementValue) # TL in cm

abund  <- measurements %>% filter(measurementType == "abundance") %>% 
  select(measurementValue) %>% rename(Abundance = measurementValue)

nrow(dt)*2 == nrow(measurements) # if there are no IDs, only proceed if this is true!!!
# merge in abundances and lengths
colnames(dt)

dt <- dt %>% select(eventID, family, genus, verbatimIdentification, taxonRank) %>% 
  bind_cols(., length, abund)

# add in sampling event information
colnames(events)
events %>% distinct(sampleSizeValue, sampleSizeUnit, samplingEffort, 
  samplingProtocol) %>% View

dt <- events %>% select(eventID, decimalLatitude, decimalLongitude,
  ends_with('Meters'), year) %>% left_join(dt, ., by = "eventID")

View(dt)
summary(dt) # a successful join would have given you no NAs
str(dt)
rm(length, measurements, events)

og_dt <- dt # backup copy

# modify before doing checks
dt$eventID <- str_remove(dt$eventID, 'BR\\:ReefSYN\\:TrindadeMartinVaz\\:BrazilianOceanicIslands\\:trindade\\_') # remove this, uninformative
dt <- dt %>% rename(Latitude = decimalLatitude, Longitude = decimalLongitude) %>% 
  mutate(Depth = (minimumDepthInMeters + maximumDepthInMeters) / 2) %>% 
  select(!ends_with('Meters')) # renames
View(dt) #check

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 5708 records
str(dt)

# Abundance and/or length, latitude and longitude numeric? Y
# Taxonomic fields must be characters or factors?

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

# YEAR MONTH DAY, already checked

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, year), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160) # run twice and the error should be gone
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, aes(fill = year), shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
table(dt$taxonRank) # chceck certainty
dt %>% filter(taxonRank %in% c('Family','Hybrid','Subfamily')) %>%
  select(family, genus, verbatimIdentification) %>% distinct %>% View

dt %>% distinct(family, genus, verbatimIdentification) %>%
  arrange(family, genus, verbatimIdentification) %>% View

dt$Species <- word(dt$verbatimIdentification, 2)
dt$Species[which(dt$taxonRank == 'Family')] <- NA
dt$Species[which(dt$taxonRank == 'Subfamily')] <- NA
dt$Species[which(dt$taxonRank == 'Genus')] <- 'sp'
dt$genus[dt$taxonRank == 'Hybrid'] <- 'Menephorus'
dt$Species[dt$taxonRank == 'Hybrid'] <- 'sp1'
dt %>% distinct(family, genus, Species) %>%
  arrange(family, genus, Species) %>% View

colnames(dt)[c(2,3,10)] <- str_to_title(colnames(dt)[c(2,3,10)])

# Prepare raw data --------------------------------------------------------

View(dt)

dt_final <- dt %>% mutate(DataID = dataID, Month = NA, Day = NA) %>% 
  rename(SamplingEventID = eventID) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID,
         Latitude, Longitude, Depth, Day, Month, Year, DataID) %>%
  arrange(Year, Month, Day, SamplingEventID, Family, Genus, Species)

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
