

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
dataid <- '72_alvarezfilip_cozumel'

# make sure your working directory is set before running these lines

dt <- read_excel('data_originals/072_alvarezfilip_caribbean/072_AlvarezFilip_Caribbean.xlsx')
dt <- dt %>% select(n, `MEAN size (cm)`, `Scientific Name`, transect, reef1)
sites <- read_excel('data_originals/072_alvarezfilip_caribbean/sites coordinates_Alvarez-Filip et al 2011 ecosphere.xlsx')[1:4]
dt <- left_join(dt, sites[-2], by=c("transect" = 'Transect_ID'))
og_dt <- dt
View(dt)
summary(dt)

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 12696 records. removed 21410 rows in our filtering above.
str(dt) # check structure

dt %>% distinct(reef1, transect) %>% group_by(reef1) %>% summarise(transectn = n())
dt <- dt %>% rename(Length = `MEAN size (cm)`, Abundance = n, Taxa = `Scientific Name`)

# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. Y
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
dt %>% filter(is.na(Length))
min(dt$Length, na.rm=T) # no zeroes? Y
summary(dt$Length)

# LAT LONG
library(rnaturalearth)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(sites %>% distinct(Longitude, Latitude), coords = c("Longitude", "Latitude"), crs = 4326)
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
dt$Taxa <- str_replace(dt$Taxa, 'sectatrix\\s\\/\\sincisor', 'sp')
dt$Genus <- word(dt$Taxa, 1)
dt$Species <- word(dt$Taxa, 2)
# everything is ok!

# Prepare raw data --------------------------------------------------------

dt_final <- dt %>% mutate(DataID = dataid,
                          SamplingEventID = paste(reef1, transect, sep = "_"),
                          Year = NA, Month = NA, Day = NA, Depth = NA, Family = NA) %>%
  group_by_at(vars(-Abundance)) %>% summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)

lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length
lengths$Length_lwr = c(2.5, 12.5, 22.5, 32.5, 47.5)
lengths$Length_upr = lengths$Length_lwr + 5
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)


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
