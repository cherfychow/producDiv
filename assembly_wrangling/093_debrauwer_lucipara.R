

# Curation Script ---------------------------------------------------------

# Dataset: A cause for hope: largely intact coral-reef communities with high reef-fish biomass in a remote Indonesian island group
# Location: Lucipara, Indonesia
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(readxl)
require(sf)
require(tidyverse)
require(lubridate)

rm(list=ls())

# Maarten already did formatting so this is just a cursory check for formatting etc.

# make sure your working directory is set before running these lines
dt <- read_excel('./data_originals/093_debrauwer_Lucipara data.xlsx', sheet = 2)
og_dt <- dt

# Structure organisation + check --------------------------------------------------

summary(dt)
str(dt)
nrow(dt) # check dimensions

sort(unique(dt$Month))
dt$Month <- 10
sort(unique(dt$SampleEventID)) # 16 sampling events, 4 transects at 4 sites.

# Secondary fields relating to sample must be factors or integers? Y
# Taxonomic fields must be characters or factors? Y


# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields. Y
# YEAR MONTH DAY, all sensical

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
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) + # confine to trpical view
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)

# Prepare raw data --------------------------------------------------------
dt_final <- dt %>% arrange(Year, Month, Day, SampleEventID) %>% 
  mutate(StudyID = "093_debrauwer_lucipara") %>% 
  rename(SamplingEventID = SampleEventID, DataID = StudyID) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID) 

# Export final ------------------------------------------------------------

readr::write_csv(dt_final, 'data_outputs/093_debrauwer_lucipara_obs.csv')

