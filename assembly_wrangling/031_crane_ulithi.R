

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
dataid <- '031_crane_ulithi'

# make sure your working directory is set before running these lines

dt <- read_excel('data_originals/031_crane_ulithi/fish_transects_2012_13_14.xlsx')
sites <- read_csv('data_originals/031_crane_ulithi/table1_sites.csv') %>% 
  select(Jurisdiction, Code, `Site name`, Coordinates)
sitename <- read_excel('data_originals/031_crane_ulithi/fish_transects_2012_13_14.xlsx')[2,7:108] %>% 
  pivot_longer(., cols = everything(), names_to = 'transect', values_to = 'Site')

dt <- dt[-c(1,2), -c(2:5)]
dt <- dt[-c(105:110)]
colnames(dt)[1:2] <- c('Taxa', 'Length')
dt <- pivot_longer(dt, cols = contains('Transect'), names_to = "transect", values_to = "Abundance")
dt <- dt %>% filter(is.na(Abundance) == F)

# for site fields
temp <- read.csv('data_originals/031_crane_ulithi/manualsitematch.csv')
colnames(temp) <- c('old', 'Location')
sitename <- left_join(sitename, temp, by=c("Site" = "old"))
sitename$Site <- NULL
dt <- left_join(dt, sitename, by="transect") # get locations field

 # import function for cooridnate conversioon
source('https://gist.githubusercontent.com/cherfychow/5473d5ba1d4ceb2fab68722373ad5ed1/raw/90aa9e438044d0c8689d18f2cf07f83dafb860a7/dms%2520to%2520decimal%2520degree')

# coord fields
sites$Latitude <- str_split(sites$Coordinates, '(?<=N)\\,\\s') %>% unlist %>%
  .[seq(1, 65, by = 2)] %>%
  str_replace_all(., '[:punct:]\\s*|[:symbol:]\\s*', ' ') %>% 
  angle2dec
sites$Longitude <- str_split(sites$Coordinates, '(?<=N)\\,\\s') %>% unlist %>% 
  .[seq(2, 66, by = 2)] %>% 
  str_replace_all(., '[:punct:]\\s*|[:symbol:]\\s*', ' ') %>% 
  angle2dec

sites <- sites %>% select(!c(Coordinates, `Site name`))
dt <- left_join(dt, sites, by=c("Location" = "Code"))

# date and transect sample fields
dt$Date <- str_extract(dt$transect, '\\d{1,2}\\/\\d{1,2}\\/\\d{2,4}') # get dates
dt$transectID <- str_extract(dt$transect, '(?<=Transect\\s)\\d+') # get transect number
sort(unique(dt$transect[is.na(dt$transectID)]))
dt$transectID[is.na(dt$transectID)] <- str_remove(dt$transect[is.na(dt$transectID)], 'Number\\:\\s{2,}') %>% str_extract(., '(?<=Transect\\s)\\d+')
dt$transectID <- as.numeric(dt$transectID)
dt$transect <- NULL

dt$Date <- as.Date(dt$Date, tryFormats = c('%m/%d/%y','%m/%d/%Y'))
summary(dt$Date)
dt$Day <- day(dt$Date)
dt$Month <- month(dt$Date)
dt$Year <- year(dt$Date)
dt$Date <- NULL
View(dt)

rm(temp, angle2dec, sitename, sites)

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 5588 records
str(dt) # check structure
summary(dt)

dt$Length = as.numeric(dt$Length)
dt$Abundance = as.numeric(dt$Abundance)

# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. Y
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
dt %>% filter(is.na(Length))
min(dt$Length, na.rm=T) # no zeroes? Y
dt$Length <- dt$Length / 10 # convert to cm

dt %>% filter(Length > 100) %>% View

summary(dt)

# LAT LONG
library(rnaturalearth)
require(sf)
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

sort(unique(dt$Taxa))

replace <- c('punctofasciatus' = 'punctatofasciatus',
             'trifasciatus' = 'trifascialis', 
             'sp\\.$' = 'sp',
             'sp\\. Blue eye' = 'sp1',
             '(?<=sp)\\.\\s(?=\\d)' = '',
             '\\*' = '',
             'scopas x flav' = 'sp1')
dt$Species <- str_replace_all(dt$Taxa, replace)
sort(unique(dt$Species))

dt$Genus <- word(dt$Species, 1)
dt$Species <- word(dt$Species, 2)
# everything is ok!

# Prepare raw data --------------------------------------------------------
n_distinct(dt$transectID)

dt %>% group_by(Jurisdiction, Location) %>% summarise(transectN = n_distinct(transectID), yearN = n_distinct(Year))
# just want to keep the time series sites
key <- dt %>% group_by(Jurisdiction, Location) %>% 
  summarise(yearN = n_distinct(Year)) %>% 
  filter(yearN > 1) %>% pull(Location)

dt_final <- dt %>% filter(Location %in% key) %>% 
  mutate(DataID = dataid, SamplingEventID = paste(Jurisdiction, Location, transectID, Year, sep = "_"),
         Family = NA, Depth = NA) %>% 
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
centroid <- dt_coord %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
