

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
dataid <- '104_brandl_crypticflux_lizard'

# make sure your working directory is set before running these lines

dt <- read_excel('data_originals/104_brandl_crypticflux.xlsx', sheet = 4)
dt <- dt %>% select(TL, Family, Genus, Species, Site, Sitenumber, Habitat, Lat, Long, Depth)
og_dt <- dt
View(dt)

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 1502 records
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

dt %>% distinct(Latitude, Longitude) %>% filter(Longitude > 0)
dt$Longitude[which(dt$Longitude < 0)] <- dt$Longitude[which(dt$Longitude < 0)] * -1

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
sum(is.na(dt$Genus) + 0) == 0
sum(str_detect(dt$Genus, '^\\s+$') + 0) == 0

# lots of distinct unknowns, need to clean with scientific name merged
dt$Taxa <- paste(dt$Genus, dt$Species, sep = ' ')
sort(unique(dt$Taxa))
replace <- unique(dt$Taxa[str_which(dt$Taxa, 'spe[A-Z]$')]) %>% sort %>% as_tibble
replace$correct <- c(rep(1, 6), 2:7, 1, 1, 1:3, 1)
replace$correct <- str_extract(replace$value, '^[A-Za-z]+\\ssp') %>% paste0(., replace$correct)
temp <- replace$correct
names(temp) <- replace$value
dt$Taxa <- str_replace_all(dt$Taxa, temp)
sort(unique(dt$Taxa))
dt$Species <- word(dt$Taxa, 2)
dt %>% distinct(Family, Genus, Species) %>% filter(str_detect(Genus, 'idae')) %>% View
dt$Genus[which(dt$Genus == 'Apogonidae' & dt$Species == 'seminigracaudus')] <- 'Apogon'
dt$Genus[which(dt$Genus == 'Apogonidae' & dt$Species == 'trimaculatus')] <- 'Pristicon'
dt$Genus[which(dt$Genus == 'Pomacentridae' & dt$Species == 'amboinensis')] <- 'Pomacentrus'
dt %>% distinct(Family, Genus, Species, Taxa) %>% arrange(Family, Genus, Species) %>% View
dt$Genus[str_which(dt$Genus, 'idae$')] <- NA

# everything is ok!

# Prepare raw data --------------------------------------------------------
n_distinct(dt$Site) == n_distinct(dt$Sitenumber)

dt_final <- dt %>% rename(SamplingEventID = Site) %>% mutate(DataID = dataid,
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
dt_coord <- dt_final %>% filter(is.na(Latitude) == F) %>% distinct(Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
