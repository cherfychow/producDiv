

# Curation Script ---------------------------------------------------------

# Dataset: Karimunjawa WCS fish survey
# Location: Karimunjawa, Indonesia
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(readr)

require(sf)
require(lubridate)
require(ggplot2)
require(rnaturalearth)

rm(list=ls())
dataID <- 'BT436_hoey_karimunjawa'

# make sure your working directory is set before running these lines
dt <- read_csv('data_originals/BT436_hoey_karimunjawa/Karimunjawa_fish surveys.csv')
summary(dt)
View(dt) #check

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
str(dt)

# Abundance and/or length, latitude and longitude numeric? Y
# Taxonomic fields must be characters or factors? Y

# size class split
dt <- dt %>% mutate(Length_lwr = case_when(str_detect(Size_Class, '\\-') ~ str_extract(Size_Class, '^\\d+'),
    str_detect(Size_Class, '\\>') ~ '40') |> as.numeric(),
  Length_upr = case_when(str_detect(Size_Class, '\\-') ~ str_extract(Size_Class, '(?<=\\-)\\d+'),
    str_detect(Size_Class, '\\>') ~ NA) |> as.numeric()) %>%
  mutate(Length = Length_lwr + (Length_upr - Length_lwr)/2) %>%
  select(!c(Obs, Region, Size_Category, `Numb. of fish`, a, b, `biomass (kg per250m2)`, Size_Class)) %>%
  rename(Abundance = `abundance (per250m2)`, Year = Survey_Period, Depth = `Depth(m)`)

summary(dt)
dt$Length[is.na(dt$Length)] <- 40 # fill it in for the ones with no upper limit

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)) == 0 # any NAs
sum(dt$Abundance == 0) == 0 # zeroes
dt %>% filter(Abundance == 0) %>% View # there were definitely at least 1 if there's length
dt$Abundance[dt$Abundance == 0] <- 1
sum(dt$Abundance == '') == 0 # blanks

sum(is.na(dt$Length)) == 0 # any NAs
sum(dt$Length == 0) == 0 # zeroes
sum(dt$Length == '') == 0 # blanks

dt$Date <- as.Date(dt$Date, '%e-%m-%Y')
dt$Month <- month(dt$Date)
dt$Day <- day(dt$Date)
summary(dt)
# YEAR MONTH DAY, all good

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.

proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- dt %>% filter(is.na(Latitude) == F) %>% 
  distinct(Longitude, Latitude, Year) %>% 
  st_as_sf(., coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160) # run twice and the error should be gone
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, aes(fill = Year), shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

dt %>% distinct(Family, Genus, Species) %>%
  arrange(Family, Genus, Species) %>% View

dt %>% distinct(Family, Genus, Species) %>%
  filter(str_detect(Species, '[:punct:]|[:symbol:]')) %>% View

dt$Species[str_which(dt$Species, '(melas\\/|cf\\.)')] <- 'sp'
dt$Species <- str_remove(dt$Species, '\\.')
dt$Family <- str_to_title(dt$Family)
sort(unique(dt$Genus))
dt$Genus <- str_to_title(dt$Genus)
dt %>% distinct(Family, Genus, Species) %>% 
  filter(Genus == 'Goby')
dt$Species[dt$Genus == 'Goby'] <- NA
dt$Genus[dt$Genus == 'Goby'] <- NA

# Sampling effort + sampling event structure check ------------------------

dt %>% group_by(`Site Name`) %>% 
  summarise(transectN = n_distinct(Transect), yearN = n_distinct(Year)) %>% View
dt$`Site Name` <- str_remove_all(dt$`Site Name`, '\\s')

# Prepare raw data --------------------------------------------------------

View(dt)

dt_final <- dt %>% 
  mutate(DataID = dataID, SamplingEventID = paste(`Site Name`, Transect, Year, sep = "_")) %>% 
  select(Abundance, Length, Length_lwr, Length_upr, Family, Genus, Species, SamplingEventID,
         Latitude, Longitude, Depth, Day, Month, Year, DataID) %>%
  arrange(Year, Month, Day, SamplingEventID, Family, Genus, Species)

View(dt_final)

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
