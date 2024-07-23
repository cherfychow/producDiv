

# Curation Script ---------------------------------------------------------

# Dataset: Blue Alliance reef fish monitoring
# Location: Philippines
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
dataID <- '131_bluealliance_pemba'

# make sure your working directory is set before running these lines
dt <- readxl::read_excel('./data_originals/131_bluealliance_pemba/BA_FVC raw_TZN_for CC.xlsx', sheet = 1)[1:18]
site <- read_csv('./data_originals/131_bluealliance_pemba/sites.csv')

colnames(site)
unique(dt$Site)
dt <- site %>% select(Site, Long, Lat) %>% left_join(dt, ., by = "Site")

# Structure check ---------------------------------------------------------

summary(dt)
nrow(dt) # check dimensions
# 2543 records
str(dt)

# Abundance and/or length, latitude and longitude numeric? Y
# Taxonomic fields must be characters or factors? Y

# remove unnecessary fields and split date
dt <- dt %>% mutate(Day = day(Date), Month = month(Date), Year = year(Date)) %>%
  select(!c(starts_with('Sample'), Time, `Visibility (m)`, 
  Observer, `Scientific name`, Date, Side))

# fix depth
sort(unique(dt$`Depth (m)`))
colnames(dt)[3] <- 'Depth'
dt$Depth <- str_remove(dt$`Depth`, 'm')
sum(is.na(dt$`Depth`))
depth <- str_split(unique(dt$Depth), '\\-')
depth2 <- lapply(depth, FUN = function(x) {
  if (length(x) == 2) {mean(as.numeric(x))} else return(x)
})
depth2 <- unlist(depth2)
names(depth2) <- unique(dt$Depth)
depth2
dt$Depth <- str_replace_all(dt$Depth, depth2) %>% as.numeric
summary(dt$Depth)

# renames
colnames(dt)
dt <- rename(dt, Longitude = Long, Latitude = Lat, Length = Size, Abundance = Count)

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Abundance)) == 0 # any NAs
dt %>% filter(is.na(Abundance)) %>% View
dt$Abundance[is.na(dt$Abundance)] <- 1 # there should be at least 1 if there's length
sum(dt$Abundance == 0) == 0 # zeroes
sum(dt$Abundance == '') == 0 # blanks

sum(is.na(dt$Length)) == 0 # any NAs
sum(dt$Length == 0) == 0 # zeroes
sum(dt$Length == '') == 0 # blanks

summary(dt$Abundance)
summary(dt$Length) # range seems reasonable
dt %>% filter(Length > 40) %>% View # yep, ok
dt <- dt %>% filter(Length < 180) # remove the huge hemigymnus melapterus
dt$Genus[dt$Genus == 'Thalassoma' & dt$Length > 40] <- NA

summary(dt) # check dates

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('Latitude', 'Longitude')]) # good

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

dt %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>%
  filter(str_detect(Species, '^sp')) %>% View

dt$Species <- str_remove(dt$Species, '\\.\\s*') |> 
  str_replace('\\(pale\\)', '1')
dt %>% distinct(Family, Genus, Species) %>% filter(is.na(Genus))
dt$Species[is.na(dt$Genus)] <- NA

# non reals
sum(dt$Species == '') == 0

# Sampling effort + sampling event IDs  -----------------------------------

str(dt)
dt %>% group_by(Site) %>% 
  summarise(replicate_n = n_distinct(Replicate)) %>% View
dt %>% group_by(Site, Replicate) %>% 
  summarise(interval_n = n_distinct(Interval)) %>% View
# interval is just segment of belt transect
dt %>% group_by(Site) %>% 
  summarise(year_n = n_distinct(Year)) %>% View
# not actually a time series once you isolate year

# Prepare raw data --------------------------------------------------------

View(dt)

dt_final <- dt %>% mutate(DataID = dataID, 
  SamplingEventID = paste(Year, Site, Replicate, sep = '_')) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID,
         Latitude, Longitude, Depth, Day, Month, Year, DataID) %>%
  arrange(Year, SamplingEventID, Family, Genus, Species)

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