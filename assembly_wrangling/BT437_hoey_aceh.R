

# Curation Script ---------------------------------------------------------

# Dataset: Aceh WCS fish survey
# Location: Aceh, Indonesia
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(readxl)
require(sf)
require(lubridate)
require(ggplot2)
require(rnaturalearth)

rm(list=ls())
dataID <- 'BT436_hoey_aceh'

# make sure your working directory is set before running these lines
dt <- read_excel('data_originals/BT437_hoey_aceh/Aceh fish census_2010-19_data for Cher_headerfix.xlsx')
coords <- read.csv('data_originals/BT437_hoey_aceh/BT437_coords.csv', header = T)
colnames(coords) <- c('Latitude', 'Longitude', 'Site')
summary(dt[1:3])
View(dt) #check

dt <- filter(dt, is.na(Family) == F, is.na(Species) == F)
dt <- tidyr::pivot_longer(dt, cols = !c(Family, Species, size_class), names_to = "transect", values_to = "Abundance")
dt <- dt %>% filter(Abundance != 0)
View(dt)

fields <- str_split_fixed(dt$transect, '\\_', n = 3)
dt <- as_tibble(fields) %>% rename(Site = V1, Survey = V2, Date = V3) %>%
  bind_cols(., dt)
dt$transect <- NULL
dt$Date <- as.Date(as.numeric(dt$Date), origin = "1899-12-30")
summary(dt$Date)

dt <- dt %>% mutate(Day = day(Date), Month = month(Date), Year = year(Date))

# deal with size classes
unique(dt$size_class)
dt <- str_split_fixed(dt$size_class, '\\-', 2) %>% as_tibble %>%
  rename(Length_lwr = V1, Length_upr = V2) %>% 
  mutate(across(everything(), as.numeric)) %>%
  mutate(Length = (Length_lwr + Length_upr) / 2) %>% bind_cols(., dt)
View(dt)
dt$Date <- NULL
dt$size_class <- NULL

coords$Site <- str_remove_all(coords$Site, '\\_')
sort(coords$Site)
sort(unique(dt$Site))
# fix to match
replace <- c('N$' = 'north',
             'S$' = 'south',
             'Leunbalee' = 'LeunBalee',
             'LeunBallee' = 'LeunBalee',
             'Lucas$' = 'LucasReef',
             'Lhokgna' = 'Lhoknga',
             'Graden' = 'Garden',
             'Hitam' = 'Itam')
dt$Site <- str_replace_all(dt$Site, replace)
sort(unique(dt$Site)) == sort(coords$Site)

dt <- left_join(dt, coords, by = "Site")

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
str(dt)
summary(dt)
# Abundance and/or length, latitude and longitude numeric? Y
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)) == 0 # any NAs
sum(dt$Abundance == 0) == 0 # zeroes
sum(dt$Abundance == '') == 0 # blanks

sum(is.na(dt$Length)) == 0 # any NAs
sum(dt$Length == 0) == 0 # zeroes
sum(dt$Length == '') == 0 # blanks

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
