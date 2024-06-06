

# Curation Script ---------------------------------------------------------

# Dataset: 104_brandl_crypticflux_moorea 
# Location: Moorea
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(tidyverse)
require(lubridate)
require(readxl)

rm(list=ls())
dataid <- '104_brandl_crypticflux_belize'

# make sure your working directory is set before running these lines

dt <- read_excel('data_originals/104_brandl_crypticflux.xlsx', sheet = 2)
dt %>% distinct(Habitat, Location, Site, Lat, Long) %>% View
dt <- dt %>% select(TL, Genus, Species, Site, Lat, Long)
og_dt <- dt
View(dt)
summary(dt)


# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 65334 records. removed 21410 rows in our filtering above.
str(dt) # check structure

dt <- dt %>% rename(Length = TL, Latitude = Lat, Longitude = Long) %>%
  mutate(Abundance = 1, .before = Length)

dt %>% select(Length) %>% View
dt$Length <- as.numeric(dt$Length)
# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. Y
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
dt %>% filter(is.na(Length)) #remove
dt <- dt %>% filter(!is.na(Length))
min(dt$Length, na.rm=T) # no zeroes? Y

summary(dt$Length)
min(dt$Length)
dt$Length <- dt$Length / 10 # convert to cm

# LAT LONG
library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160)
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
sum(is.na(dt$Genus) + 0) == 0
sum(str_detect(dt$Genus, '^\\s+$') + 0) == 0
dt %>% filter(Genus == 'NA') %>% View
dt <- dt %>% filter(!Genus == 'NA') # remove completely unknown

dt %>% distinct(Genus, Species) %>% arrange(Genus, Species) %>% View
dt$Species <- str_replace(dt$Species, 'sp.$', 'sp')
dt$Family <- NA
dt$Family[str_which(dt$Genus, 'idae$')] <- dt$Genus[str_which(dt$Genus, 'idae$')]
dt$Genus[str_which(dt$Genus, 'idae$')] <- NA
dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View

# everything is ok!

# Prepare raw data --------------------------------------------------------

dt_final <- dt %>% mutate(DataID = dataid,
                          SamplingEventID = paste(Site, sep = "_"),
                          Depth = NA, Day = NA, Month = NA, Year = NA) %>% # not including swath because the true sample is both swathes together
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
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
