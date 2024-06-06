

# Curation Script ---------------------------------------------------------

# Dataset: Reef Life Survey- Fish abundance and biomass
# Location: Global
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(tidyverse)
require(sf)

rm(list=ls())
dataID <- '054_RLS'

# make sure your working directory is set before running these lines
dt <- read_csv('data_originals/054_RLS_2024April/IMOS_-_National_Reef_Monitoring_Network_Sub-Facility_-_Global_reef_fish_abundance_and_biomass.csv', skip = 70)
readme <- read_csv('data_originals/054_RLS_2024April/IMOS_-_National_Reef_Monitoring_Network_Sub-Facility_-_Global_reef_fish_abundance_and_biomass.csv', n_max = 69)
View(dt)
summary(dt)
str(dt)
og_dt <- dt

# exclude anything not tropical/subtropical
dt <- dt %>% filter(latitude >= -30, latitude <= 30)
# the geom column has much more precise coordinate data

dt <- st_as_sf(dt, wkt = "geom")
dt$Latitude <- st_coordinates(dt)[,2]
dt$Longitude <- st_coordinates(dt)[,1]

dt <- st_drop_geometry(dt) %>% select(total, size_class, family, species_name, reporting_name, 
                                      survey_id, location, site_code, site_name, method, 
                                      program, Latitude, Longitude, depth, survey_date)
dt$program <- NULL

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 695793 records
str(dt)
# Abundance and/or biomass, latitude and longitude numeric? Y

# Dates
require(lubridate)
dt$Day <- day(dt$survey_date)
dt$Month <- month(dt$survey_date)
dt$Year <- year(dt$survey_date)
dt$survey_date <- NULL

# Secondary fields relating to sample must be factors or integers?

# Taxonomic fields must be characters or factors? Y

dt <- dt %>% rename(Abundance = total, Length = size_class, Family = family, Depth = depth)

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)+0)
min(dt$Abundance, na.rm=T) # no zeroes? Y
dt %>% filter(Abundance == 0) %>% View
dt$Abundance[which(dt$Abundance == 0)] <- 1 # seems to be data entry error.
sum(is.na(dt$Length)+0) # 44604 NAs for length but might be for size classes
min(dt$Length, na.rm=T) # no zeroes? Y
problemsurveys <- dt %>% filter(Length == 0) %>% distinct(survey_id) %>% pull(survey_id)
dt %>% filter(survey_id %in% problemsurveys) %>% View # are these zeroes one off cases
dt <- dt %>% filter(Length > 0) # filter these out

sum(is.na(dt$Length)+0) # no NAs in size_classes

# YEAR MONTH DAY
# non-temporal dataset. ignore.

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid()
world <- st_break_antimeridian(world, lon_0 = 160)
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
sum(is.na(dt$species_name) + 0) == 0
sum(str_detect(dt$species_name, '^\\s+$') + 0) == 0

dt %>% distinct(species_name, reporting_name) %>% arrange(species_name) %>% View
# use reporting name

# uncertains
dt %>% distinct(Family, reporting_name) %>% filter(str_detect(reporting_name, 'spp.$')) %>% View
dt <- dt %>% filter(!is.na(Family)) %>% mutate(Species = str_replace(reporting_name, 'spp\\.$', 'sp'))
dt %>% distinct(Family, reporting_name) %>% filter(str_detect(reporting_name, 'id\\sspp')) %>% View
dt$Species[str_which(dt$Species, 'id\\ssp')] <- NA # remove family level IDs from species column
dt$Genus <- word(dt$Species, 1)
dt$Species <- word(dt$Species, 2, -1)

dt$reporting_name <- NULL
dt$species_name <- NULL

# Prepare raw data --------------------------------------------------------

n_distinct(dt$site_code)
n_distinct(dt$survey_id)

# check survey effort evenness
dt %>% group_by(site_code, Year) %>% summarise(nsurvey = n_distinct(survey_id)) %>% View

dt_final <- dt %>% mutate(SamplingEventID = paste(site_code, survey_id, Year, Month, Day, sep = "_"),
                          DataID = dataID) %>% 
  select(!c(survey_id, location, site_code, site_name, method)) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)

# size class intervals
lengths <- dt_final %>% filter(Length <= 62.5) %>% distinct(Length) %>% arrange(Length)
lengths$Length_lwr <- NA
lengths$Length_upr <- NA
lengths$Length_upr[1:12] <- lengths$Length[1:12] + (lengths$Length[2:13] - lengths$Length[1:12])/2 # mid point with i + 1
lengths$Length_lwr[2:13] <- lengths$Length[2:13] - (lengths$Length[2:13] - lengths$Length[1:12])/2 # mid opint with i - 1
lengths$Length_upr[13] <- 62.5 + 6.25
lengths$Length_lwr[1] <- 2.5 - 1.25
lengths <- tibble(Length = seq(75, max(dt_final$Length), by = 12.5)) %>% 
  mutate(Length_lwr = Length - 6.25, Length_upr = Length + 6.25) %>% 
  bind_rows(lengths, .)

dt_final <- left_join(dt_final, lengths, by = "Length") %>% relocate(Length_lwr, Length_upr, .after = Length)

# Export final ------------------------------------------------------------

write_csv(dt_final, paste0('data_outputs/054_RLS/', dataID, '_obs.csv'))

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt_final %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long

