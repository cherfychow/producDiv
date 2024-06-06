

# Curation Script ---------------------------------------------------------

# Dataset: Nearshore marine biodiversity of Osa Peninsula, Costa Rica: where the ocean meets the rainforest
# Location: Costa Rica
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(sf)
require(world)
require(ggplot2)
require(readxl)

rm(list=ls())

# make sure your working directory is set before running these lines
dt <- read_excel('./data_originals/025_friedlander_costarica/Osa_fish_data_Datadryad.xlsx', sheet = 2)
og_dt <- dt

# Structure organisation + check --------------------------------------------------

summary(dt)
str(dt)
nrow(dt) # check dimensions

dt <- dt %>% select(Date, Location, STATION, Geomorph, Lat, Long, DEPTH_M, TRANSECT, 
                    FishBaseName, Family, NUMBER, Min, Max, Avg, Area)

# make this traceable
dt$index <- 1:nrow(dt)

# try to expand lengths
dt$Length <- 0
dt$Length[which(dt$NUMBER == 1)] <- dt$Min[which(dt$NUMBER == 1)]
sum(is.na(dt$Max)+0)
dt$Length[which(is.na(dt$Max))] <- dt$Min[which(is.na(dt$Max))]

dt_expanded <- dt %>% filter(Length > 0) # move over records that don't need length expansion
dt <- dt %>% filter(Length == 0)
dt$mean <- (dt$Min + dt$Max) / 2
dt %>% filter(Avg == mean, Length == 0) %>% View # records where uniform distribution could be used


for (i in which(dt$Avg == dt$mean)) {
  n = dt$NUMBER[i]
  min = dt$Min[i]
  max = dt$Max[i]
  if (exists("temp") == FALSE) {
    temp <- tibble(Length = runif(n = n, min = min, max = max) %>% round(., 0)) %>% 
      bind_cols(dt[rep(i, n), c(1:10, 15:16)], .)
  } else {
    temp <- tibble(Length = runif(n = n, min = min, max = max) %>% round(., 0)) %>% 
      bind_cols(dt[rep(i, n), c(1:10, 15:16)], .) %>% 
      bind_rows(temp, .)
  }
}

summary(temp)
temp %>% filter(is.na(Length)) %>% View
problems <- sort(unique(temp$index[which(is.na(temp$Length))]))

n_distinct(temp$index) == n_distinct(dt$index[dt$Avg == dt$mean])

# aggregate records of the same length
temp <- temp %>% filter(is.na(Length) == F) %>% 
  group_by_all() %>% summarise(NUMBER = n())
dt_expanded <- dt_expanded %>% select(!c(Min, Max, Avg)) %>% bind_rows(., temp)
summary(dt_expanded)
dt <- dt %>% filter(Avg != mean) %>% bind_rows(dt[which(dt$index %in% problems), ], .) # problem records

# rows 1-5 inverted min max
tempmax <- dt$Max[1:5]
dt$Max[1:5] <- dt$Min[1:5]
dt$Min[1:5] <- tempmax
rm(temp)
# the others are wildly wrong averages that I'm just going to ignore as entry error
for (i in 1:nrow(dt)) {
  n = dt$NUMBER[i]
  min = dt$Min[i]
  max = dt$Max[i]
  if (exists("temp") == FALSE) {
    temp <- tibble(Length = runif(n = n, min = min, max = max) %>% round(., 0)) %>% 
      bind_cols(dt[rep(i, n), c(1:10, 15:16)], .)
  } else {
    temp <- tibble(Length = runif(n = n, min = min, max = max) %>% round(., 0)) %>% 
      bind_cols(dt[rep(i, n), c(1:10, 15:16)], .) %>% 
      bind_rows(temp, .)
  }
}

temp <- temp %>% filter(is.na(Length) == F) %>% group_by_all() %>% summarise(NUMBER = n())

dt_expanded <- bind_rows(dt_expanded, temp) %>% arrange(Date, Location, STATION, DEPTH_M, TRANSECT, FishBaseName)

rm(temp)
dt <- dt_expanded %>% select(!index)
rm(dt_expanded)

require(lubridate)
dt$Year <- year(dt$Date)
dt$Month <- month(dt$Date)
dt$Day <- day(dt$Date)

str(dt)
dt %>% group_by(Area) %>% summarise(tally = n()) # differing sampling effort
# split these into a different study later
n_distinct(dt$STATION)
n_distinct(dt$Location)
n_distinct(dt$TRANSECT)

# 3 transects at 26 sites at 6 regions

# Secondary fields relating to sample must be factors or integers? Y
# Taxonomic fields must be characters or factors? Y

dt <- dt %>% rename(Abundance = NUMBER, Taxa = FishBaseName, 
         Latitude = Lat, Longitude = Long, Depth = DEPTH_M)

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)+0)
min(dt$Abundance, na.rm=T) # no zeroes? Y
sum(is.na(dt$Length)+0) # 44604 NAs for length but might be for size classes
min(dt$Length, na.rm=T) # no zeroes? Y

# YEAR MONTH DAY, all sensical

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, Location), coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(., crs = proj)
world <- sf::st_as_sf(maps::map('world', plot = FALSE, fill = TRUE))
ggplot() +
  geom_sf(data = world, color='grey40', fill='grey95') +
  geom_sf(data = points, aes(fill = Location), shape = 21, size = 2) +
  guides(fill = 'none') + theme_minimal() + 
  coord_sf(crs = proj)

# zoom in
ggplot() +
  geom_sf(data = world, color='grey40', fill='grey95') +
  geom_sf(data = points, aes(fill = Location), shape = 21, size = 2) +
  guides(fill = 'none') + theme_minimal() + 
  coord_sf(crs = 4326, xlim = c(-85,-80), ylim = c(5,10))


# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
sum(is.na(dt$Taxa) + 0) == 0
sum(str_detect(dt$Taxa, '^\\s+$') + 0) == 0

sort(unique(dt$Taxa))
dt$Genus <- word(dt$Taxa, 1)
dt$Species <- word(dt$Taxa, 2, -1)
dt$Taxa <- NULL
sort(unique(dt$Family))

# Prepare raw data --------------------------------------------------------
dt %>% distinct(Location, STATION, Geomorph, TRANSECT) %>% group_by(Location, Geomorph) %>% 
  summarise(nstation = n_distinct(STATION))
dt %>% distinct(Location, STATION, Geomorph, TRANSECT) %>% group_by(Location, Geomorph, STATION) %>% 
  summarise(ntransect = n_distinct(TRANSECT)) %>% View
# Station = distinct Location * geomorph

dt_final <- dt %>% arrange(Date, STATION, TRANSECT) %>% 
  mutate(SamplingEventID = paste(STATION, TRANSECT, sep = "_"), 
                          DataID = "025_friedlander_costarica") %>% 
  select(!c(STATION, Location, TRANSECT, Geomorph, Date)) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID) 

# Export final ------------------------------------------------------------

# Area is due to the size dependent belt transect. convert to density

dt_final$Abundance[which(dt_final$Area == 100)] <- dt_final$Abundance[which(dt_final$Area == 100)]/2

readr::write_csv(dt_final %>% select(!Area), 'data_outputs/025_friedlander_costarica/025_friedlander_costarica_obs.csv')

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
