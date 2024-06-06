

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

# make sure your working directory is set before running these lines
dt <- read.delim('./data_originals/020_noaa_micronesia/occurrence.txt', header=T, sep = "\t")
og_dt <- read.delim('./data_originals/020_noaa_micronesia/occurrence.txt', header=T, sep = "\t") # original copy for later comparisons
lengths <- read.delim('./data_originals/020_noaa_micronesia/extendedmeasurementorfact.txt', header=T, sep = "\t")
events <- read.delim('./data_originals/020_noaa_micronesia/event.txt', header=T, sep = "\t")
sites <- read.csv('./data_originals/020_noaa_micronesia/sitemeta.csv', header = T)

sort(unique(dt$Site))
sort(unique(sites$Site))

# Structure organisation + check --------------------------------------------------

# merge in lengths
summary(dt)
lengths <- lengths %>% filter(measurementType == 'Size') # only use the total length size records
unique(lengths$measurementUnit)
dt <- dt %>% select(occurrenceID, individualCount, eventID, eventDate, decimalLatitude, decimalLongitude, family, genus, specificEpithet) %>% 
  left_join(., lengths[c(2,4)], by = "occurrenceID")
View(dt)
dt$occurrenceID <- NULL
rm(lengths)

require(lubridate)
dt$Year <- year(dt$eventDate)
dt$Month <- month(dt$eventDate)
dt$Day <- day(dt$eventDate)
dt$eventDate <- NULL
colnames(dt)[1:8] <- c('Abundance', 'SampleEventID', 'Latitude', 'Longitude', 'Family', 'Genus', 'Species', 'Length')

nrow(dt) # check dimensions

sort(unique(dt$Length)) # all numbers
dt$Length <- as.numeric(dt$Length)

dt <- events %>% mutate(Depth = minimumDepthInMeters + maximumDepthInMeters / 2) %>% 
  select(eventID, Depth) %>% left_join(dt, ., by = c('SampleEventID' = 'eventID')) # add depth

dt$SampleEventID <- as.factor(dt$SampleEventID)
sort(unique(dt$SampleEventID))
dt$SampleEventID <- str_remove_all(dt$SampleEventID, 'UoGML\\:MC\\:Fish\\:')
#extract site and transect replicate to check sampling consistency
dt$Site <- str_extract(dt$SampleEventID, '(?<=site\\_)[:upper:]{3}\\-[:digit:]+(?=\\:)') %>% as.factor
dt$Transect <- str_extract(dt$SampleEventID, '(?<=\\:SPC\\_rep\\_)[:digit:]+$') %>% as.factor
sort(unique(dt$Site))
sort(unique(dt$Transect))
dt$SampleEventID <- ''

# Primary field check -----------------------------------------------------

summary(dt)

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Length)+0)
min(dt$Length, na.rm=T) # no zeroes? Y

sum(is.na(dt$Abundance)+0)
min(dt$Abundance, na.rm=T) # no zeroes? Y

# YEAR MONTH DAY
# again, no negative values, 0s or NAs, and all are logical
# Year < 2020, month < 12, day < 31
summary(dt) # all good

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
  geom_point(data=dt %>% distinct(Longitude, Latitude, Year), 
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

dt %>% distinct(Site, Transect, Year) %>% group_by(Site) %>% summarise(samples = n()) %>% View
dt %>% distinct(Site, Transect, Year) %>% group_by(Site, Transect) %>% summarise(sampledyears = n()) %>% View
events %>% distinct(id, locationID, year) %>% group_by(locationID) %>% summarise(samples = n()) %>% View
# year samples

# same number of replicates?
dt %>% distinct(Site, Transect) %>% group_by(Site) %>% summarise(ntransect = n()) %>% View
# through time
dt %>% distinct(Site, Transect, Year) %>% group_by(Site, Year) %>% summarise(ntransect = n()) %>% View


# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
summary(dt)
length(which(dt$Family == '' & dt$Genus == '' & dt$Species == ''))
length(which(dt$Genus == '' & dt$Species == '')) # how many family level records 
# blanks or spaces in all taxonomic columns

sort(unique(paste(dt$Genus, dt$Species, sep = ' ')))
# GBIF dataset so no misspellings here
# fill in sp where specific epithet is blank
dt$Species[which(dt$Genus != '' & dt$Species == '')] <- 'sp'
sort(unique(paste(dt$Genus, dt$Species, sep = ' ')))

# Final restructuring ---------------------------------------------------

# split by site MPA status
mpa <- sites$Site[sites$mpa == "Yes"] %>% unique %>% sort %>% str_to_upper
dt_mpa <- dt %>% filter(Site %in% mpa)

dt_np <- dt %>% filter(!(Site %in% mpa)) # all others

# remove any fields that don't pertain to sample event level or final format
dt_mpa_merged <- dt_mpa %>% group_by_at(vars(-Abundance)) %>% summarise(Abundance = sum(Abundance)) %>% # aggregate records within sampling events
  mutate(SampleEventID = paste(Year, Site, Transect, sep = "_")) 

dt_np_merged <- dt_np %>% group_by_at(vars(-Abundance)) %>% summarise(Abundance = sum(Abundance)) %>% # aggregate records within sampling events
  mutate(SampleEventID = paste(Year, Site, Transect, sep = "_")) 

nrow(dt_mpa) - nrow(dt_mpa_merged)
nrow(dt_np) - nrow(dt_np_merged)

# Prepare raw data --------------------------------------------------------

# fill in blank columns, reorder, arrange rows
dt_mpa_merged <-  dt_mpa_merged %>% ungroup %>% arrange(Year, Month, Day, Site, Transect) %>% 
  select(-c(Site, Transect)) %>% 
  relocate(Abundance, Length, Family,	Genus,	Species,	SampleEventID,	Latitude,	Longitude,	Depth,	Day,	Month,	Year)

dt_np_merged <-  dt_np_merged %>% ungroup %>% arrange(Year, Month, Day, Site, Transect) %>% 
  select(-c(Site, Transect)) %>% 
  relocate(Abundance, Length, Family,	Genus,	Species,	SampleEventID,	Latitude,	Longitude,	Depth,	Day,	Month,	Year)


# Export final ------------------------------------------------------------

write.csv(dt_mpa_merged, 'data_outputs/020_noaa_micronesia_MPA_rawdata.csv', row.names = F)
write.csv(dt_np_merged, 'data_outputs/020_noaa_micronesia_rawdata.csv', row.names = F)

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# MPA
dt_coord <- dt_mpa_merged %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long

# transform to mercator to get area in sq km
area <- st_transform(dt_coord, st_crs("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")) %>%
  st_convex_hull() %>% st_area()  ##get area in sq km
write_clip(area)

# non MPA
dt_coord <- dt_np_merged %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long

# transform to mercator to get area in sq km
area <- st_transform(dt_coord, st_crs("+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")) %>%
  st_convex_hull() %>% st_area()  ##get area in sq km
write_clip(area)

# Plot the geometries -----------------------------------------------------

require(ggplot2)
world_map <- map_data('world') # check whether the GPS coordinates match expectations
ggplot(world_map) +
  geom_polygon(aes(x=long, y=lat, group=group), color='grey', fill='transparent') +
  geom_point(data = dt_merged, aes(x = Longitude, y = Latitude), size = 2, shape = 21, colour = 'blue') +
  geom_point(aes(x = centroid[1], y = centroid[2]), size = 3, colour = 'green') +
  geom_sf(data = dt_coord %>% st_convex_hull(), colour = 'blue', fill = 'transparent') +
  coord_sf(xlim = c(135, 180), ylim = c(2,10)) + # make sure these fit your points
  labs(x = NULL, y = NULL) +
  theme_minimal()
