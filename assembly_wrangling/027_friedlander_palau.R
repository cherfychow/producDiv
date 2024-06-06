

# Curation Script ---------------------------------------------------------

# Dataset: Multiple datasets from Body size, reef area and temperature predict global reef-fish species richness across spatial scales
# Location: Global
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(tidyverse)
require(sf)
require(world)
require(readxl)

rm(list=ls())

# make sure your working directory is set before running these lines
dt <- read_excel('data_originals/027_friedlander_palau.xlsx')
dt_og <- dt
View(dt)

# mostly scuba methods, but one survey is with snorkel --> split
summary(dt)
str(dt)
dt$Max <- as.numeric(dt$Max)
dt$mean <- (dt$Min + dt$Max) / 2
dt$Length <- 0
dt$Length[which(dt$NUMBER == 1)] <- dt$Min[which(dt$NUMBER == 1)]

dt_expanded <- dt %>% filter(Length > 0) # move over records that don't need length expansion
dt %>% filter(Avg == mean, Length == 0) %>% View # records where uniform distribution could be used

for (i in which(dt$Avg == dt$mean & dt$Length == 0)) {
  n = dt$NUMBER[i]
  min = dt$Min[i]
  max = dt$Max[i]
  if (exists("temp") == FALSE) {
    temp <- tibble(Length = runif(n = n, min = min, max = max) %>% round(., 0), index = i) %>% 
      bind_cols(dt[rep(i, n), 1:15], .)
  } else {
    temp <- tibble(Length = runif(n = n, min = min, max = max) %>% round(., 0), index = i) %>% 
      bind_cols(dt[rep(i, n), 1:15], .) %>% 
      bind_rows(temp, .)
  }
}

# did everyone get included
sum((sort(unique(temp$index)) == sort(unique(which(dt$Avg == dt$mean & dt$Length == 0)))) + 0) == n_distinct(which(dt$Avg == dt$mean & dt$Length == 0))
# aggregate records of the same length
temp <- temp %>% select(!index) %>% group_by_all() %>% summarise(NUMBER = n())
dt_expanded <- dt_expanded %>% select(!c(Min, Max, Avg, mean)) %>% bind_rows(., temp)
dt <- dt %>% filter(Length == 0)
dt_expanded <- dt %>% filter(is.na(Max)) %>% mutate(Length = Min) %>% select(!c(Min, Max, Avg, mean)) %>% 
  bind_rows(dt_expanded, .)
dt <- dt %>% filter(!is.na(Max))
dt <- dt %>% filter(Avg != mean) # nothing left!
rm(dt, temp)

unique(dt_expanded$DATATYP)
summary(dt_expanded)
dt <- dt_expanded
rm(dt_expanded)

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 8435 records
str(dt)
# Abundance and/or biomass, latitude and longitude numeric? Y

# Dates
require(lubridate)
dt <- dt %>% mutate(Day = day(Date), Month = month(Date), Year = year(Date)) %>% 
  select(!c(Date, Time))

# Secondary fields relating to sample must be factors or integers? Y

# Taxonomic fields must be characters or factors? Y

meta <- dt %>% distinct(STATION, Reef, Lat, Long, MPA, Habitat, Mangt)
dt <- dt %>% select(NUMBER, Length, `NAME-FISHBASE`, STATION, Reef, TRANSECT, Lat, Long, Depth_m, Day, Month, Year) %>% 
  rename(Abundance = NUMBER, Taxa = `NAME-FISHBASE`, Latitude = Lat, Longitude = Long, Depth = Depth_m)

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

proj = '+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, Reef), coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(., crs = proj)
world <- sf::st_as_sf(maps::map('world', plot = FALSE, fill = TRUE))
ggplot() +
  geom_sf(data = world, color='grey40', fill='grey95') +
  geom_sf(data = points, aes(fill = Reef), shape = 21, size = 2) +
  guides(fill = 'none') + theme_minimal() + 
  coord_sf(crs = proj)

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
sum(is.na(dt$Taxa) + 0) == 0
sum(str_detect(dt$Taxa, '^\\s+$') + 0) == 0

sort(unique(dt$Taxa))
dt$Taxa <- str_remove(dt$Taxa, '\\.$')
sort(unique(dt$Taxa))
dt$Genus <- word(dt$Taxa, 1)
dt$Species <- word(dt$Taxa, 2, -1)
dt$Taxa <- NULL

# Prepare raw data --------------------------------------------------------

dt_final <- dt %>% mutate(SamplingEventID = paste(Reef, STATION, TRANSECT, sep = "_"), 
                          Family = NA, DataID = "027_friedlander_palau") %>% 
  select(!c(STATION, Reef, TRANSECT)) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year)

# Export final ------------------------------------------------------------

write_csv(dt_final, 'data_outputs/027_friedlander_palau_obs.csv')

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt_final %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
