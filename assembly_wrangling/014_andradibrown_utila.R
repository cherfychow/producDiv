

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(ggplot2)
require(readr)
require(stringr)
require(lubridate)
require(readxl)
require(sf)

rm(list=ls())

dataid <- '014_andradibrown_utila'
# load data in
dt <- read_csv('data_originals/014_andradibrown_honduras/pone.0156641.s002.csv')
sites <- read_csv('data_originals/014_andradibrown_honduras/coords.csv')

View(dt)

# Structure organisation 
colnames(dt) <- str_to_title(colnames(dt))
dt <- left_join(dt, sites, by = 'Site')
og_dt <- dt
colnames(dt)[c(2,8)] <- c('Depth', 'Length')
dt$Weight.g <- NULL
dt$Name <- NULL

# overall checks
str(dt)
summary(dt)
nrow(dt)

# COORDINATE CHECK
library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160)
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, shape = 21, size = 2, fill = 'green') +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)


# ABUNDANCE/LENGTH CHECKS

# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Length)) # NAs
sum(dt$Length == 0, na.rm = T) # zeroes
sum(dt$Length == '', na.rm = T) # blanks

dt <- dt %>% filter(is.na(Length) == F) # omit the one NA


# TAXONOMIC CHECKS

# No NAs in taxonomic fields, remove all non-organism records
paste(dt$Genus, dt$Species, sep = ' ') %>% is.na %>% sum
dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View
dt$Species[is.na(dt$Family) == F & is.na(dt$Genus) == T] <- 'sp'
dt$Species[is.na(dt$Genus) == F & is.na(dt$Species) == T] <- 'sp'

dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View

# SAMPLING STRUCTURE, CONSISTENCY CHECK
dt %>% distinct(Site, Depth, Transect) %>% View
dt %>% distinct(Site, Depth, Transect) %>% 
  group_by(Site, Depth) %>% 
  summarise(transectN = n_distinct(Transect)) %>% View
# some replicates are inconsistent for deeper depths, but it's ok for non-time-series

# Prepare raw data --------------------------------------------------------

dt_final <- dt %>% mutate(DataID = dataid,
                          SamplingEventID = paste(Site, Depth, Transect, sep = "_")) %>% 
  # aggregate records if they're at the same sampling event, species, and size
  group_by_all() %>% summarise(Abundance = n()) %>% ungroup %>% 
  mutate(Length = Length / 10, Year = 2014, Month = NA, Day = NA) %>% # fill in unavailable fields, convert length to cm
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
         Latitude, Longitude, Depth, Day, Month, Year, DataID)


# Export final ------------------------------------------------------------

write_csv(dt_final, paste0('data_outputs/', dataid, '_obs.csv'))

# Convex Hull for centroid ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt_final %>% distinct(Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
