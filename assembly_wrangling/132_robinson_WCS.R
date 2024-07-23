

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
dataID <- '132_robinson_WCS'

# make sure your working directory is set before running these lines
dt <- read_csv('./data_originals/132_robinson_WCS/fish_individuals.csv')

# Structure check ---------------------------------------------------------

summary(dt)
dt %>% group_by(project) %>% summarise(duration = max(year)-min(year)+1)

colnames(dt)
dt <- dt %>% select(count, size, fish_family, fish_genus, fish_taxon, year, project, site, 
  transect_number, depth, transect_width, transect_length, latitude, longitude, sample_date)

str(dt)
# Abundance and/or length, latitude and longitude numeric? Y
# Taxonomic fields must be characters or factors?
summary(dt)

dt <- dt %>% mutate(Month = month(sample_date), Day = day(sample_date)) %>% select(!sample_date)
colnames(dt)[1:4] <- c('Abundance', 'Length', 'Family', 'Genus')
colnames(dt) <- str_to_title(colnames(dt))

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
dt %>% filter(Length > 70) %>% View # yep, all sharks and megafauna species

# YEAR MONTH DAY, already checked

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

is.na(dt$Fish_taxon) |> sum() == 0
(dt$Fish_taxon == '') |> sum() == 0
is.na(dt$Genus) |> sum() == 0
is.na(dt$Family) |> sum() == 0

dt$Species <- word(dt$Fish_taxon, 2)
dt %>% distinct(Family, Genus, Species, Fish_taxon) %>% 
  arrange(Family, Genus, Species) %>% View

dt %>% filter(is.na(Species)) %>% distinct(Fish_taxon) %>% pull %>% sort
dt$Species[is.na(dt$Species) & is.na(dt$Genus) == F] <- 'sp'
dt$Fish_taxon <- NULL
dt %>% distinct(Family, Genus, Species) %>% 
  arrange(Family, Genus, Species) %>% View

# Check sampling effort + structure --------------------------------------------------------

dt %>% group_by(Project, Site) %>% summarise(ntransect = n_distinct(Transect_number)) %>% View
summary(dt)
dt %>% group_by(Project, Site) %>% summarise(nyear = n_distinct(Year)) %>% View

dt %>% distinct(Project, Transect_width, Transect_length) %>% View
# inconsistent within project, standardise it all to density
dt <- mutate(dt, Abundance = 100 * Abundance / (Transect_width * Transect_length))

unique(dt$Project)
dt$Project <- str_remove_all(dt$Project, '\\s')
dt$Project <- str_remove_all(dt$Project, "\\d{4}[:punct:]*|\\'|(June)|(July)|(December)")

# Prepare raw data --------------------------------------------------------

View(dt)

dt_final <- dt %>% mutate(DataID = paste(dataID, Project, sep = '_'), 
  SamplingEventID = paste(Year, Site, Transect_number, sep = '_')) %>% 
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
