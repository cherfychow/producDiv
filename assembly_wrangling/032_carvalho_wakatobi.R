

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(ggplot2)
require(readr)
require(stringr)
require(lubridate)
require(readxl)
require(sf)

rm(list=ls())
dataid <- '032_carvalho_indonesia'

# make sure your working directory is set before running these lines
dt <- read_excel('data_originals/032_carvalho_indonesia/fish_size_spectra_data.xlsx')
sites <- read_csv('data_originals/032_carvalho_indonesia/drivers_size_spectra_data.csv')
og_dt <- dt
View(dt)
str(dt)
summary(dt)



# Structure organisation --------------------------------------------------

dt <- dt %>% select(abundance, size_cm, family, genus, species, region, site_name, transect)
dt <- sites %>% select(region, site_name, lat_dd, lon_dd) %>% 
  left_join(dt, ., by = c('region', 'site_name'))
colnames(dt) <- str_to_title(colnames(dt))

summary(dt)

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 18776 rows
str(dt) # check structure and data types
# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. NA
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
min(dt$Abundance) # all 1

# Lengths
sum(dt$Size_cm == '')
# looked fine in summary too

dt <- dt %>% rename(Latitude = Lat_dd, Longitude = Lon_dd)

# LAT LONG
library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160)
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, shape = 21, size = 2, fill = 'green') +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal()

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
sum(is.na(dt$Species)) == 0
sum(str_detect(dt$Species, '^\\s+$')) == 0

dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View
dt %>% filter(Genus == 'all') %>% View
dt$Genus[dt$Genus == 'all'] <- NA

dt$Species <- str_replace_all(dt$Species, 'spp', 'sp')
dt$Species[str_which(dt$Species, 'idae$')] <- 'sp'
dt$Species[dt$Species == 'scarinae'] <- 'sp'
dt$Family <- str_to_title(dt$Family)
sort(unique(dt$Genus))

# final check
dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View

# Prepare raw data --------------------------------------------------------

# check sampling structure
dt %>% distinct(Region, Site_name, Transect) %>% View
dt %>% distinct(Region, Site_name, Transect) %>% group_by(Region, Site_name) %>% 
  summarise(TransectN = n_distinct(Transect)) %>% View

dt_final <- dt %>% mutate(DataID = dataid,
                          SamplingEventID = paste(Region, str_remove(Site_name, '\\s+'), Transect, sep = "_")) %>% 
  # aggregate records if they're at the same sampling event, species, and size
  group_by_at(vars(-Abundance)) %>% summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
  rename(Length = Size_cm) %>% 
  mutate(Year = 2018, Month = NA, Day = NA, Depth = NA) %>% # fill in unavailable fields
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
dt_coord <- dt_final %>% distinct(Longitude, Latitude) %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
