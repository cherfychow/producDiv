

# Setup
require(tidyverse)
require(lubridate)
require(readxl)

rm(list=ls())
dataid <- '110_brandl_farms_saudi'

# make sure your working directory is set before running these lines

dt <- read_excel('data_originals/110_brandl_farms/farms_alldata.xlsx', sheet = 4)
colnames(dt)
dt %>% distinct(site, reef, shelf, exposure, habitat, unit) %>% View

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 114 records
str(dt) # check structure
summary(dt)

dt <- dt %>% select(!c(shelf, habitat, starts_with('date'))) %>% 
  rename(Length = tl) %>%
  mutate(Abundance = 1, .before = Length) # add approximate coordinates

# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. Y
# Taxonomic fields must be characters or factors? Y

colnames(dt) <- str_to_title(colnames(dt))

# Coordinate check

library(sf)
library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Long, Lat), coords = c("Long", "Lat"), crs = 4326)
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

# Primary field check -----------------------------------------------------

# ABUNDANCES/LENGTH
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Length))
summary(dt$Length)
dt$Length <- dt$Length / 10 # convert to cm

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View

dt$Species <- word(dt$Species, 1)

# everything is ok!

# Prepare raw data --------------------------------------------------------
dt %>% group_by(Reef, Site) %>% summarise(unitN = n_distinct(Unit))

dt_final <- dt %>% group_by_at(vars(-Abundance)) %>% 
  summarise(Abundance = sum(Abundance)) %>% ungroup %>% #aggregate by length
  rename(Latitude = Lat, Longitude = Long)%>% 
  mutate(DataID = dataid,
         SamplingEventID = paste(Reef, Site, Unit, sep = "_"),
         Year = 2019, Month = NA, Day = NA) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)

View(dt_final)

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
