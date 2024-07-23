

# Curation Script ---------------------------------------------------------

# Dataset: ReefSYN
# Location: Guarapari
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
dataID <- '133_hadjhammou_SERF'

# make sure your working directory is set before running these lines
dt <- read_csv('data_originals/133_hadjhammou_serf/SERF_FishSizeData_TL.csv')

# Structure check ---------------------------------------------------------

summary(dt)
head(dt, 150) |> View()
head(dt$UniqueTransect)
head(dt$UniqueSite)

n_distinct(dt$UniqueSite) == (paste0(dt$Site_Lat, dt$Site_Long) |> n_distinct())
paste0(dt$Site_Lat, dt$Site_Long) |> n_distinct()
n_distinct(dt$UniqueSite)

dt_time <- dt %>% distinct(Site_Lat, Site_Long, UniqueSite, Year) %>% group_by(Site_Lat, Site_Long) %>% 
  mutate(duration = max(Year)-min(Year) + 1) %>%
  filter(duration >= 2) %>% pull(UniqueSite)
dt <- dt %>% filter(UniqueSite %in% dt_time)

colnames(dt)
dt <- dt %>% filter(LengthTypeSurvey == 'TL') %>%
  select(Number, Size, Family, Genus, Species, UniqueSite, 
  UniqueTransect, SampArea, Site_Lat, Site_Long, CensusMethod,
  Day, Month, Year)

str(dt)
# Abundance and/or length, latitude and longitude numeric? Y
# Taxonomic fields must be characters or factors? Y
unique(dt$CensusMethod)
summary(dt)

# rename
dt <- rename(dt, Abundance = Number, Length = Size, 
  Latitude = Site_Lat, Longitude = Site_Long)

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)) == 0 # any NAs
sum(dt$Abundance == 0) == 0 # zeroes
sum(dt$Abundance == '') == 0 # blanks
dt %>% filter(Abundance == 0) %>% View
dt$Abundance[dt$Abundance == 0] <- 1 # there's length, so should be at least 1

sum(is.na(dt$Length)) == 0 # any NAs
sum(dt$Length == 0) == 0 # zeroes
sum(dt$Length == '') == 0 # blanks

summary(dt$Abundance)
summary(dt$Length) # range seems reasonable
dt %>% filter(Length > 70) %>% View
# remove these improbable records
dt <- dt[-which(with(dt, Genus == 'Monotaxis' & Species == 'grandoculis' & Length > 60)),]
dt <- dt[-which(with(dt, Genus == 'Gomphosus' & Species == 'varius' & Length > 40)),]
dt <- dt[-which(with(dt, Genus == 'Ctenochaetus' & Species == 'striatus' & Length > 30)),]

# YEAR MONTH DAY
dt$Day[dt$Day < 1] <- NA
unique(dt$Month)
month <- as.character(1:12)
names(month) <- c('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEPT*', 'OCT', 'NOV', 'DEC')
dt$Month <- str_replace_all(dt$Month, month) |> as.numeric()
dt$Month[dt$Month < 1] <- NA

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, Year), coords = c("Longitude", "Latitude"), crs = 4326)
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

is.na(with(dt, paste(Genus, Species, sep = ' '))) |> sum() == 0
(with(dt, paste(Genus, Species, sep = ' ')) == '') |> sum() == 0
is.na(dt$Genus) |> sum() == 0
is.na(dt$Family) |> sum() == 0

dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View

dt$Genus[str_which(dt$Genus, 'idae$')] <- NA
dt$Species[str_which(dt$Species, 'idae$')] <- NA
dt$Family[dt$Family == 'Scaridae'] <- 'Labridae'
dt$Species <- str_remove(dt$Species, '\\.')

with(dt, paste(Genus, Species, sep = ' ')) |> unique() |> sort()

# Check sampling effort + structure --------------------------------------------------------

dt %>% distinct(UniqueSite, UniqueTransect, SampArea) %>%
  group_by(UniqueSite) %>% summarise(areaN = n_distinct(SampArea)) %>% View
# some transects have different transect areas
# standardise abundance to density/250 sq m
dt <- mutate(dt, Abundance = 250 * Abundance/SampArea)
dt$CensusMethod[dt$CensusMethod == 'Standard belt transect'] <- 'belt'
dt$CensusMethod[dt$CensusMethod == 'Point intercept'] <- 'point'

# Prepare raw data --------------------------------------------------------

View(dt)

dt_final <- dt %>% mutate(DataID = paste(dataID, CensusMethod, sep = '_'), 
  SamplingEventID = paste(Year, UniqueTransect, sep = '_')) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID,
         Latitude, Longitude, Depth, Day, Month, Year, DataID) %>%
  arrange(Year, Month, Day, SamplingEventID, Family, Genus, Species)

View(dt_final)

# Export final ------------------------------------------------------------

for (i in 1:n_distinct(dt_final$DataID)) {
  dt_final %>% filter(DataID == unique(dt_final$DataID)[i]) %>%
    write_csv(., paste0('data_outputs/', dataID, '/', unique(dt_final$DataID)[i], '_obs.csv'))
}

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt_final %>% distinct(DataID, Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% 
  group_by(DataID) %>%
  summarise(Latitude = st_union(geometry) %>% st_centroid %>% unlist %>% .[2],
    Longitude = st_union(geometry) %>% st_centroid %>% unlist %>% .[1]) %>%
  st_drop_geometry()
