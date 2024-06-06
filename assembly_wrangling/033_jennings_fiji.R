
## Curation script for Wilson and Graham: reef fish
## curator VB

library(stringr)
library(lubridate)
library(dplyr)
library(readr)

dt <- read_csv("data_originals/033_jennings_fiji/Fiji_reef_fish_Jennings_Polunin_1997.csv")
coords <- read_csv('data_originals/033_jennings_fiji/sitescoords.csv')
colnames(coords) <- c('Longitude', 'Latitude', 'ground')
dataid <- '033_jennings_fiji'

str(dt)
sites <- read_csv('data_originals/033_jennings_fiji/Fiji_habitat_Jennings_Polunin_1997.csv', skip = 1)
dt <- sites %>% select(area, site, depth) %>% left_join(dt, ., by = c('ground' = 'area', 'area' = 'site'))
dt <- left_join(dt, coords[1:3], by = "ground")

colnames(dt) <- str_to_title(colnames(dt))
dt$Aphia_id <- NULL
dt$Abundance <- 1

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

# Check abundances and lengths
sum(dt$Length == 0)
sum(is.na(dt$Length))
sum(dt$Length == '')

# check taxonomy
sort(unique(dt$Species))

dt$Genus <- word(dt$Species, 1)
sort(unique(dt$Genus))
dt$Species <- word(dt$Species, 2)

# Sampling event check
unique(dt$Replicate) # 12 replicates at 10 sites

# prepare for export

dt_final <- dt %>% group_by_at(vars(-Abundance)) %>% summarise(Abundance = sum(Abundance)) %>% 
  mutate(DataID = dataid,
         SamplingEventID = paste(Ground, Area, Replicate, sep = "_"),
         Year = NA, Day = NA, Month = NA, Family = NA) %>% ungroup %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
         Latitude, Longitude, Depth, Day, Month, Year, DataID) %>% 
  arrange(SamplingEventID, Genus, Species)

n_distinct(dt_final$SamplingEventID) # 840 sampling events

# Export final ------------------------------------------------------------

write_csv(dt_final, paste0('data_outputs/', dataid, '_obs.csv'))

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

