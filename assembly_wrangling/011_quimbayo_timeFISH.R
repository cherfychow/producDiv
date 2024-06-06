

# Curation Script ---------------------------------------------------------

# Dataset: TimeFISH: Long-term assessment of reef fish assemblages in a transition zone in the Southwestern Atlantic
# Location: SW Atlantic  
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(sf)
require(world)
require(ggplot2)

rm(list=ls())

# make sure your working directory is set before running these lines
dt <- read.csv('./data_originals/011_quimbayo_timeFISH/ECY22-0768_TimeFISH_census_data.csv', header=T)
sites <- read.csv('./data_originals/011_quimbayo_timeFISH/ECY22-0768_TimeFISH_location_information.csv', header=T)
taxa <- read.csv('./data_originals/011_quimbayo_timeFISH/ECY22-0768_TimeFISH_taxonomic_information.csv', header=T)
View(dt)

# mostly scuba methods, but one survey is with snorkel --> split
summary(dt)
summary(sites)
sites$sampling_method <- as.factor(sites$sampling_method)

# Structure organisation --------------------------------------------------

# add the metadata from the sites sheet
# first check if transect ids line up
n_distinct(dt$transect_id)
n_distinct(sites$transect_id)
sum((sort(unique(dt$transect_id)) == sort(unique(sites$transect_id))) + 0)
dt <- sites %>% select(transect_id, province, state, location, site, 
                       no_take_zone, longitude, latitude, depth, sampling_method, day, month, year) %>% 
  distinct(transect_id, .keep_all = T) %>% left_join(dt, ., by="transect_id")

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 27377 records
str(dt) # check structure
summary(dt)
# Abundance and/or biomass, latitude and longitude numeric? Y

# dates need to be coerced
require(lubridate)
dt$month <- with(dt, paste(day, month, year, sep = " ")) %>% dmy() %>% month
summary(dt) # logical? yes
# Year, month and day must be integers or factors?
dt$month <- as.integer(dt$month)

# Secondary fields relating to sample must be factors or integers?
for (i in c(6:10)) {
  dt[,i] <- as.factor(dt[,i])
} # all the location + transect fields

# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. Y
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$total_length)+0)
min(dt$total_length, na.rm=T) # no zeroes? Y
sum(is.na(dt$abundance)+0)
min(dt$abundance, na.rm=T) # no zeroes? Y

unique(dt$total_length) %>% sort

# YEAR MONTH DAY
# again, no negative values, 0s or NAs, and all are logical
# Year < 2020, month < 12, day < 31
# already checked
# LAT LONG
sum(dt[c('latitude', 'longitude')] == "") == 0 # no blanks?
sum(is.na(dt[c('latitude', 'longitude')])) == 0 # no NAs?

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('latitude', 'longitude')]) # good

world_map <- map_data('world') # check whether the GPS coordinates match expectations
world <- ggplot(world_map) + 
  geom_polygon(aes(x=long, y=lat, group=group), color='grey', fill='transparent') +
  coord_fixed() + 
  labs(x='Longitude', y='Latitude') +
  theme_bw() + theme(panel.grid=element_blank())
points <- world + geom_point(data=sites %>% distinct(longitude, latitude), aes(x=longitude, y=latitude), shape=21)
points
points_zoom <- points + coord_fixed(xlim=c(-55,-40), ylim=c(-40,-15))
points_zoom # all looks good

# check sampling space x time
# just so there is not a lone color by itself
world + geom_point(data=sites, aes(x=longitude, y=latitude, color = year), size = 2, alpha = 0.5) + 
  coord_fixed(xlim=c(-55,-40), ylim=c(-30,-20))
rm(list=c('world_map', 'world', 'points', 'points_zoom'))
sites %>% group_by(province, state, location, site) %>% summarise(duration = max(year) - min(year) + 1) %>% View 
# duration ranges 1 to 16 years by site

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
sum(is.na(dt$species_name) + 0) == 0
sum(str_detect(dt$species_name, '^\\s+$') + 0) == 0

dt$species_name <- str_replace_all(dt$species_name, '\\_', ' ')
sort(unique(dt$species_name))
sort(unique(dt$species_name)) %>% word(., 2)
sort(unique(word(dt$species_name, 1)))

require(worrms)
taxonclean <- wm_records_names(unique(dt$species_name), fuzzy = F, marine_only = T)
taxonclean <- bind_rows(taxonclean)
taxonclean %>% filter(rank != 'Species') # any species that aren't species level
replace <- taxonclean$valid_name[which(taxonclean$status == 'unaccepted')] # replace unaccepted names
names(replace) <- taxonclean$scientificname[which(taxonclean$status == 'unaccepted')]
dt$species_name <- str_to_sentence(dt$species_name)
taxa$genus <- str_to_sentence(taxa$genus)
dt$species_code <- NULL
# replace
dt$species_name <- str_replace_all(dt$species_name, replace)

sort(unique(dt$species_name))

# check that genera are genera, not family names (-idae/eae)
dt$Genus <- word(dt$species_name, 1)
dt$Genus %>% str_detect(., 'idae$|ini$') %>% sum(. + 0)

dt$Species <- word(dt$species_name, start=2, end=-1) # keep only specific epithet
dt <- taxonclean %>% distinct(family, genus) %>% left_join(dt, ., by=c("Genus" = "genus"))
dt$species_name <- NULL

# Prepare raw data --------------------------------------------------------
colnames(dt) <- str_to_title(colnames(dt))
colnames(dt)[2] <- 'Length'

dt_sc <- dt %>% filter(Sampling_method == 'scuba') %>% 
  mutate(DataID = '011_quimbayo_timefish',
         SamplingEventID = paste(Province, State, Location, Site, Transect_id, sep = "_")) %>% # not including swath because the true sample is both swathes together
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
         Latitude, Longitude, Depth, Day, Month, Year, DataID)

dt_sn <- dt %>% filter(Sampling_method == 'snorkel') %>% 
  mutate(DataID = '011_quimbayo_timefish',
         SamplingEventID = paste(Province, State, Location, Site, Transect_id, sep = "_")) %>% # not including swath because the true sample is both swathes together
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
         Latitude, Longitude, Depth, Day, Month, Year, DataID)

# Export final ------------------------------------------------------------

write.csv(dt_sc, 'data_outputs/011_quimbayo_timeFISH_scuba_obs.csv', row.names=F)
write.csv(dt_sn, 'data_outputs/011_quimbayo_timeFISH_snorkel_obs.csv', row.names=F)

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt_sc_p %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long

dt_coord <- dt_sc %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long

dt_coord <- dt_sn_p %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long