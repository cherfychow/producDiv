

# Curation Script ---------------------------------------------------------

# Dataset: ReefSYN
# Location: Guarapari
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(sf)
require(lubridate)
require(ggplot2)
require(rnaturalearth)

rm(list=ls())
dataID <- '130_protectWPS_philippines'

# make sure your working directory is set before running these lines
dt <- as.list(rep('', 5))
dt[[1]] <- readxl::read_excel('data_originals/130_protectWPS_Philippines/ProducDiv_datatemplate_ph_jpc.xlsx', sheet = 2)
dt[[2]] <- readxl::read_excel('data_originals/130_protectWPS_Philippines/ProducDiv_datatemplate_ph_CORVA3_jpc.xlsx', sheet = 2)
dt[[3]] <- readxl::read_excel('data_originals/130_protectWPS_Philippines/ProducDiv_datatemplate_ph_DENRRIC_jpc.xlsx', sheet = 2)
dt[[4]] <- readxl::read_excel('data_originals/130_protectWPS_Philippines/ProducDiv_datatemplate_ph_Protect-WPS_jpc.xlsx', sheet = 2)
dt[[5]] <- readxl::read_excel('data_originals/130_protectWPS_Philippines/ProducDiv_datatemplate_ph_WWF_jpc.xlsx', sheet = 2)

merged <- bind_rows(dt)
dt <- merged
rm(merged)

summary(dt) # a successful join would have given you no NAs
str(dt)
og_dt <- dt # backup copy

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 50043 records
str(dt)

# Abundance and/or length, latitude and longitude numeric? N

# convert to decimal degree coordinates
angle2dec <- function(angle) {
  angle <- as.character(angle)
  
  if (sum(str_detect(angle, '\\s*[:upper:]$')) > 0) {
    angle <- str_remove(angle, '\\s*[:upper:]$')
    warning("You passed coordinates containing N, E, S, or W. Stringr may have accounted for it but double check your outputs")
  }
  
  angle <- str_remove(angle, '[^\\d\\.]+$') # remove the last symbol
  # use punctuation marks and symbols that aren't decimal points to split degrees minutes seconds apart
  x <- do.call(rbind, str_split(angle, '[^\\d\\.]+'))
  if (ncol(x) == 3) {
    x <- apply(x, 1L, function(y) {
    y <- as.numeric(y)
    y[1] + y[2]/60 + y[3]/3600
  })
  } else {
    x <- apply(x, 1L, function(y) {
      y <- as.numeric(y)
      y[1] + y[2]/60
    })
  }
  return(x)
}

dt <- mutate(dt, Latitude2 = angle2dec(Latitude), Longitude2 = angle2dec(Longitude))
summary(dt[15:16])
dt <- dt %>% select(!c(Latitude, Longitude)) %>%
  rename(Latitude = Latitude2, Longitude = Longitude2)
# Taxonomic fields must be characters or factors?
str(dt)

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)) == 0 # any NAs
sum(dt$Abundance == 0) == 0 # zeroes
sum(dt$Abundance == '') == 0 # blanks

sum(is.na(dt$Length)) == 0 # any NAs
sum(dt$Length == 0) == 0 # zeroes
sum(dt$Length == '') == 0 # blanks

summary(dt$Abundance)
summary(dt$Length) # range seems reasonable

# YEAR MONTH DAY, already checked

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.

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

# No NAs in taxonomic fields, remove all non-organism records
dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View
dt$Species[dt$Species == 'sp.'] <- 'sp'

# Check sampling event IDs ---------------------------------------------------

n_distinct(dt$SamplingEventID)
sort(unique(dt$SamplingEventID))

# Prepare raw data --------------------------------------------------------

View(dt)

dt_final <- dt %>% mutate(DataID = dataID) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID,
         Latitude, Longitude, Depth, Day, Month, Year, DataID) %>%
  arrange(Year, Month, Day, SamplingEventID, Family, Genus, Species)

View(dt_final)

# Export final ------------------------------------------------------------

readr::write_csv(dt_final, file.path('data_outputs', paste0(dataID, '_obs.csv')))

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt_final %>% distinct(Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
