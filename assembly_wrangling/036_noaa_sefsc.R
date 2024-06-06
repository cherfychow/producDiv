

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
dataID <- '036_noaa_sefsc'

# get the SEFSC reef visual census data
dt <- getSampleData(years = 1999:2023, regions = c('DRY TORT', 'FLA KEYS', 'SEFCRI', 'FBGMNS', 'STTSTJ', 'STX'))
sp <- getTaxonomicData()
dt <- dt %>% filter(NUM > 0) %>% select(!c(X, UNDERWATER_VISIBILITY, HABITAT_CD, STATION_NR, 
                                           SPECIES_NR, TIME_SEEN, PROT, SUBREGION_NR, DEPTH_CLASS:DEPTH_STRAT, ADMIN))
summary(dt)
dt <- dt %>% group_by_at(vars(-NUM)) %>% summarise(Abundance = sum(NUM))

# check sampling effort
dt %>% group_by(ZONE_NR, REGION, SUBREGION, STRAT) %>% 
  summarise(samples = n_distinct(PRIMARY_SAMPLE_UNIT)) %>% View
dt <- left_join(dt, sp[1:3], by = 'SPECIES_CD') %>% select(!SPECIES_CD)


# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 2,171,926 records
str(dt)
# Abundance and/or biomass, latitude and longitude numeric? Y

# Dates, all good

# Secondary fields relating to sample must be factors or integers? Y

# Taxonomic fields must be characters or factors? Y


# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)+0)
min(dt$Abundance, na.rm=T) # no zeroes? Y
sum(is.na(dt$LEN)+0) # no NAs
min(dt$LEN, na.rm=T) # negative?
dt %>% filter(LEN < 0) %>% View
dt$LEN[which(dt$LEN < 0)] <- NA # assume these negative values are the fillers that should be NAs
dt %>% filter(LEN > 200) %>% View # some fish are > 2 m long??
# length fixes
dt$LEN[which(dt$LEN == 1500)] <- 150 # nurse shark max length 430 cm
dt$LEN[which(dt$LEN == 450)] <- 430 # nurse shark max length 430 cm
dt$LEN[which(dt$LEN == 1120)] <- 112 # tarpon max length 250 cm
dt <- dt[-which(dt$LEN == 536 | dt$LEN == 769),] # remove these records. no idea

# YEAR MONTH DAY
# already checked

dt$SPECIES_CD <- NULL
colnames(dt) <- str_to_sentence(colnames(dt))
colnames(dt)
colnames(dt)[c(5,6,11,17)] <- c('Latitude', 'Longitude', 'Length', 'Taxon')

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

sort(unique(dt$Taxon))
dt <- ungroup(dt)
dt %>% distinct(Family, Taxon) %>% View
# No NAs in taxonomic fields, remove all non-organism records

dt <- dt %>% filter(is.na(Family) == F & Family != '') # remove these uncertain records
replace = c('sp.$' = 'sp',
            '^[a-z]+\\sspecies$' = '',
            'grouper-sea bass species' = '', 
            'species$' = 'sp',
            'unknown sp' = '')
dt$Species <- str_replace_all(dt$Taxon, replace)
sort(unique(dt$Species))
dt$Species[str_which(dt$Species, 'idae\\s')] <- NA
dt$Species[str_count(dt$Species, '\\s') > 1 & is.na(dt$Species) == F] <- word(dt$Species[str_count(dt$Species, '\\s') > 1 & is.na(dt$Species) == F], 1, 2)

dt$Genus <- word(dt$Species, 1)
sort(unique(dt$Genus))
dt$Species <- word(dt$Species, 2)
sort(unique(dt$Family))
dt$Family <- str_remove(dt$Family, '\\\xa0')
dt$Taxon <- NULL

# Prepare raw data --------------------------------------------------------

n_distinct(dt$Mapgrid_nr)
n_distinct(dt$Strat)
unique(dt$Strat)
n_distinct(dt$Primary_sample_unit)
unique(dt$Region)
unique(dt$Subregion)
unique(dt$Zone_nr)

dt_final <- dt %>% mutate(DataID = dataID, 
                          SamplingEventID = paste(Region, Subregion, Zone_nr, Strat, Primary_sample_unit, sep = "_")) %>% 
  select(!c(Mapgrid_nr:Mpa_nr, Strat, Region, Subregion, Primary_sample_unit)) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
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
