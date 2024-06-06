
## Curation script for Wilson and Graham: reef fish
## curator VB

library(ggplot2)
library(stringr)
library(lubridate)
library(dplyr)
library(sf)
library(rnaturalearth)
library(readr)

ds <- read_csv("data_originals/BT692_grahamwilson_seychelles/Graham_Wilson Seychelles FISH data_BioTIME.csv")
dataid <- '692_grahamwilson_seychelles'

# transform to long format
str(ds)

names(ds)
ds <- ds %>% select(abundance.500m2, length, Family, Species, count, Location, Habitat, Island, lat, lon, Year)
str(ds)

# get to know the data 
table(ds$count, ds$Year)# some years less replicates
table(ds$Location, ds$Year) # some years less sites
unique(ds$lat)
unique(ds$lon)
names(ds)
dim(ds)

summary(ds)

colnames(ds)[c(1,2,9,10)] <- c('Abundance', 'Length', 'Latitude', 'Longitude')

# Check abundances and lengths

sum(ds$Abundance == 0)
sum(is.na(ds$Abundance))
sum(ds$Abundance == '')

sum(ds$Length == 0)
sum(is.na(ds$Length))
sum(ds$Length == '')

# Check coordinates
library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(ds %>% distinct(Longitude, Latitude, Year), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160)
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal()

ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points, aes(color = Year), shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(crs = 4326, 
           ylim = c(min(ds$Latitude) - 2, max(ds$Latitude) + 2), 
           xlim = c(min(ds$Longitude) - 2, max(ds$Longitude) + 2)) +
  theme_minimal()

rm(proj, points, world)
detach('package:rnaturalearth')
detach('package:ggplot2')

# check taxonomy
sort(unique(ds$Species))
table(ds$Species)
replace <- c('Ambyglyphidodon' = 'Amblyglyphidodon',
             'Cantherines' = 'Cantherhines',
             'zanzibariensis' = 'zanzibarensis',
             'prasiognathos' = 'prasiognathus',
             'siganus' = 'Siganus',
             'chrysopterus' = 'chrysopterum')

ds$Species <- str_replace_all(ds$Species, replace)
ds$Genus <- word(ds$Species, 1)
sort(unique(ds$Genus))
ds$Species <- word(ds$Species, 2)

ds %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View

# Sampling event check

ds %>% distinct(Year, Island, Location, count) %>% 
  arrange(Island, Location, Year, count) %>% View

ds %>% group_by(Island, Location) %>% 
  summarise(YearN = n_distinct(Year), transectN = n_distinct(count)) %>% View

# Cousin only visited for 5 years


# prepare for export

dt_final <- ds %>% mutate(DataID = dataid,
                          SamplingEventID = paste(str_remove(Location, '\\s'), count, Year, sep = "_"),
                          Depth = NA, Day = NA, Month = NA) %>% # not including swath because the true sample is both swathes together
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
         Latitude, Longitude, Depth, Day, Month, Year, DataID)

n_distinct(dt_final$SamplingEventID) # 1647 sampling events

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
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
