

# Curation Script ---------------------------------------------------------

# Dataset: NOAA SEFSC South Florida National Coral Reef Monitoring Program
# Location: Florida and Caribbean regions
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(tidyverse)
require(sf)
require(rvc)

rm(list=ls())
dataID <- '022_noaa_puertorico'

# get the SEFSC reef visual census data
dt <- read_csv('data_originals/CRCP_Reef_Fish_Surveys_Puerto_Rico_9da3_c6c5_1e35.csv')
dt <- dt[-1, ] # remove first row of units

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 20k+ records
str(dt)
dt$latitude <- as.numeric(dt$latitude)
dt$longitude <- as.numeric(dt$longitude)
dt$LEN <- as.numeric(dt$LEN)
summary(dt)

# Abundance and/or biomass, latitude and longitude numeric? Y

# Dates, all good

# Secondary fields relating to sample must be factors or integers? Y

# Taxonomic fields must be characters or factors? Y

dt <- dt %>% rename(Abundance = NUM, Length = LEN, Latitude = latitude, 
                    Longitude = longitude, Taxa = SCIENTIFIC_NAME)

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)+0)
min(dt$Abundance, na.rm=T) # no zeroes? Y
sum(is.na(dt$Length)+0) # no NAs
min(dt$Length, na.rm=T) # negative?

dt %>% filter(Length > 100) %>% View # check large fish records

# YEAR MONTH DAY
# already checked in summary
dt$Date_UTC <- NULL

colnames(dt) <- str_to_sentence(colnames(dt))

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, Year), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid()
world <- st_break_antimeridian(world, lon_0 = 160)
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, aes(fill = Year), shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(0, 4366274.76)) +
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

sort(unique(dt$Taxa))
dt[is.na(dt$Taxa), ] %>% View
dt <- dt %>% filter(is.na(dt$Taxa) == F) # omit NAs in taxonomy
dt %>% filter(str_detect(Taxa, 'species')) %>% View
dt %>% filter(str_detect(Taxa, 'species')) %>% distinct(Taxa) %>% pull(Taxa)
dt$Taxa <- str_replace(dt$Taxa, 'Halichoeres species', 'Halichoeres sp')
replace = c('sp.' = 'sp',
            'blenny' = 'Blenniidae',
            'goby' = 'Gobiidae',
            'razorfish' = 'Centriscidae',
            'porgy' = 'Sparidae',
            'snapper' = 'Lutjanidae',
            'species' = 'sp')
dt$Taxa <- str_replace_all(dt$Taxa, replace)
dt$Family <- NA
dt$Family[str_which(dt$Taxa, 'idae\\s')] <- word(dt$Taxa[str_which(dt$Taxa, 'idae\\s')], 1)

dt$Genus <- word(dt$Taxa, 1)
sort(unique(dt$Genus))
dt$Species <- word(dt$Taxa, 2)
sort(unique(dt$Species))
dt$Taxon <- NULL
dt$Species[dt$Species == 'spgleri'] <- 'spengleri'
dt$Genus[str_which(dt$Genus, 'idae$')] <- NA

dt %>% distinct(Family, Genus, Species) %>% View

# Prepare raw data --------------------------------------------------------

n_distinct(dt$Mapgrid_nr)
n_distinct(dt$Station_nr)
n_distinct(dt$Strat)
unique(dt$Strat)
n_distinct(dt$Primary_sample_unit)

dt_final <- dt %>% mutate(DataID = dataID, 
                          SamplingEventID = paste(Strat, Primary_sample_unit, sep = "_")) %>% 
  rename(Depth = Sample_depth) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)

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
