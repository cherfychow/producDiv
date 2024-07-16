

# Curation Script ---------------------------------------------------------

# Dataset: Blue Alliance reef fish monitoring
# Location: Philippines
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
dataID <- '131_bluealliance_philippines'

# make sure your working directory is set before running these lines
dt <- readxl::read_excel('data_originals/131_bluealliance_philippines/BA_FVCraw.xlsx', sheet = 1)[1:18]
site <- read_csv('data_originals/131_bluealliance_philippines/sites.csv')

colnames(site)
dt <- site %>% select(Date, Location, Site, Long, Lat) %>%
  full_join(., dt, by = c("Date", "Location", "Site"))

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 10239 records
str(dt)

# coerce dates before anything
dt$Date <- as.Date(dt$Date, origin = "1899-12-30")
summary(dt)
# Abundance and/or length, latitude and longitude numeric? Y
dt %>% distinct %>% filter(is.na(Long)) %>% View
# Taxonomic fields must be characters or factors?

# remove unnecessary fields and split date
dt <- dt %>% mutate(Day = day(Date), Month = month(Date), Year = year(Date)) %>%
  select(!c(starts_with('Sample'), Time, `Visibility (m)`, 
  Observer, `Scientific name`, Date, Side))

# fix depth
sort(unique(dt$`Depth (m)`))
dt$Depth <- str_remove(dt$`Depth (m)`, 'm')
sum(is.na(dt$`Depth (m)`))
depth <- str_split(unique(dt$Depth), '\\-')
depth2 <- lapply(depth, FUN = function(x) {
  if (length(x) == 2) {mean(as.numeric(x))} else return(x)
})
depth2 <- unlist(depth2)
names(depth2) <- unique(dt$Depth)
depth2 <- depth2[-20]
depth2
dt$Depth <- str_replace_all(dt$Depth, depth2) %>% as.numeric
summary(dt$Depth)
dt$`Depth (m)` <- NULL

# renames
colnames(dt)
dt <- rename(dt, Longitude = Long, Latitude = Lat, Length = Size, Abundance = Count)

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
dt %>% filter(Length > 50) %>% View # yep, ok

summary(dt) # check dates

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- dt %>% filter(is.na(Latitude) == F) %>%
  distinct(Longitude, Latitude, Year) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326)
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

dt %>% distinct(Family, Genus, Species) %>%
  arrange(Family, Genus, Species) %>% View

# non reals
sum(is.na(dt$Species)) == 0
sum(dt$Species == '') == 0
# uncertains
sum(str_detect(dt$Species, '\\.')) # genus level
sum(str_detect(dt$Species, '\\d')) # known uncertains

dt %>% distinct(Genus, Species) %>%
  filter(str_detect(Species, '\\.|\\d|[:punct:]|[:symbol:]')) %>% View
dt$Species <- str_remove_all(dt$Species, '\\.\\s*|\\s\\([a-zA-Z\\s]+\\)') # clean up


# Sampling effort + sampling event IDs  -----------------------------------

str(dt)
dt %>% group_by(Location, Site) %>% 
  summarise(replicate_n = n_distinct(Replicate)) %>% View
dt %>% group_by(Location, Site, Replicate) %>% 
  summarise(interval_n = n_distinct(Interval)) %>% View
# interval is just segment of belt transect
dt %>% group_by(Location, Site) %>% 
  summarise(year_n = n_distinct(Year)) %>% View
# not actually a time series once you isolate year

# Prepare raw data --------------------------------------------------------

View(dt)

dt_final <- dt %>% mutate(DataID = dataID, 
  SamplingEventID = paste(Year, Location, Site, Replicate, sep = '_')) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID,
         Latitude, Longitude, Depth, Day, Month, Year, DataID) %>%
  arrange(Year, SamplingEventID, Family, Genus, Species)

View(dt_final)

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
