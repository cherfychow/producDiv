

# Curation Script ---------------------------------------------------------

# Dataset: 017_brandl_uaecryptic
# Location: UAE
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(tidyverse)
require(lubridate)
require(readxl)

rm(list=ls())
dataid <- '017_brandl_uaecryptic'

# make sure your working directory is set before running these lines

dt <- read_csv('data_originals/017_brandl_UAEcryptic/UAE18_RawData_Comms.csv')
dt <- dt %>% select(TL, Family, Genus, Species, Site, Sitename, Sitenumber, Lat, Long, Depth)
og_dt <- dt
View(dt)

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 1592 records
str(dt) # check structure
summary(dt)

dt <- dt %>% rename(Length = TL, Latitude = Lat, Longitude = Long) %>%
  mutate(Abundance = 1, .before = Length)

# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. Y
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
dt %>% filter(is.na(Length))
min(dt$Length, na.rm=T) # no zeroes? Y
dt$Length <- dt$Length / 10 # convert to cm

# LAT LONG
library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude) %>% filter(is.na(Latitude) == F), coords = c("Longitude", "Latitude"), crs = 4326)
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

dt %>% distinct(Family, Genus, Species) %>% View

dt$Species <- str_replace(dt$Species, 'speA', 'sp1')
dt %>% distinct(Family, Genus, Species) %>% filter(str_detect(Genus, 'idae$')) %>% View

# everything is ok!

# Prepare raw data --------------------------------------------------------
n_distinct(dt$Site)
n_distinct(dt$Sitenumber)
n_distinct(dt$Sitename)

dt_final <- dt %>% mutate(DataID = dataid,
                          SamplingEventID = paste(Site, Sitenumber, sep = "_"),
                          Year = NA, Month = NA, Day = NA) %>% # not including swath because the true sample is both swathes together
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
dt_coord <- dt %>% filter(is.na(Latitude) == F) %>% distinct(Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
