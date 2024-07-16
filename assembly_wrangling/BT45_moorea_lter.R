

# Curation Script ---------------------------------------------------------

# Dataset: STUDY 45 MCR LTER: Coral Reef: Long-term Population and Community Dynamics: Fishes, ongoing since 2005
# Location: Moorea
# Curator: Cher Chow
# Date: 25-Oct-2021

# Set up ------------------------------------------------------------------
# load the necessary packages
require(tidyverse)
require(lubridate)
require(readxl)

rm(list=ls())
dataid <- 'BT45_moorea_lter'

# make sure your working directory is set before running these lines
dt <- read_csv('data_originals/BT45_moorea_lter/MCR_LTER_Annual_Fish_Survey_20230615.csv')
og_dt <- dt
View(dt)
summary(dt)

# Structure organisation --------------------------------------------------

unique(dt$Site) # see an NA
n_distinct(dt$Location)
dt[6:9] <- dt[6:9] %>% mutate(across(.cols=everything(), .fns=as.factor))
summary(dt[6:9]) # habitat has an error?

dt[which(is.na(dt$Site)),] # just the row total at the end. remove.
dt <- dt[-which(is.na(dt$Site)),]

dt %>% filter(Habitat == 'error') %>% select(Location, Habitat) # all forereef.

# add approximate depth information from methods
depth <- data.frame(Habitat=sort(unique(dt$Habitat)), Depth=NA)
depth[2] <- c(1.5, 12, 10)

# I have to take the bounding coordinates for each site and get a centroid. This is only applicable for the site level
sites <- data.frame(Site=sort(unique(dt$Site)), Latitude=NA, Longitude=NA)
sites[1,2:3] <- c((-17.47185366 + -17.48641792)/2, (-149.8455917 + -149.829821) /2) #LTER 0
sites[2,2:3] <- c((-17.46576169 + -17.48131958)/2, (-149.8116849 + -149.7961685)/2)
sites[3,2:3] <- c((-17.50382025 + -17.52087158)/2, (-149.7708619 + -149.7519968)/2)
sites[4,2:3] <- c((-17.53305021 + -17.55064263)/2, (-149.7772857 + -149.7566866)/2)
sites[5,2:3] <- c((-17.56818162 + -17.59182383)/2, (-149.8869755 + -149.8561009)/2)
sites[6,2:3] <- c((-17.50735955 + -17.52839766)/2, (-149.934537 + -149.9115336)/2)


dt <- left_join(dt, sites, by="Site")
dt <- left_join(dt, depth, by="Habitat")
rm(sites, depth)

# just grab the fields we need
dt <- dt %>% select(Count, Total_Length, Family, Taxonomy, Location, Site, Habitat, Transect, Swath, Latitude, Longitude, Depth, Date, Year)
View(dt) # check the joining and filtering

# rename fields
colnames(dt)[1:2] <- c('Abundance', 'Length')
str(dt)

# dates need to be coerced
dt$Day <- day(dt$Date) %>% as.integer()
dt$Month <- month(dt$Date) %>% as.integer()
dt$Date <- NULL

# before we run any checks, the data contributor noted sampling inconsistencies in 2005 as well as variability in data recorder/diver in transect 3
# I'd omit these before running any checks
dt %>% group_by(Site, Habitat, Year) %>% summarise(transects = length(unique(Transect))) %>% View
dt <- dt %>% filter(!Year == 2015) # omit 2015 from the time series

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 65334 records. removed 21410 rows in our filtering above.
str(dt) # check structure


# logical?
summary(dt[,c('Year', 'Month', 'Day')]) # yes

summary(dt)
# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. NA
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# transform to density first
dt <- mutate(dt, Abundance = case_when(Swath == 1 ~ Abundance * 5, .default = Abundance * 1))
# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
dt %>% filter(is.na(Length))
dt %>% filter(is.na(Abundance)) # blanks in biomass? biomass more complete
min(dt$Abundance, na.rm=T) # no zeroes? Y
min(dt$Length, na.rm=T) # no zeroes? Y
dt %>% filter(Length <= 0) %>% View
# a lot of unmeasurable lengths
dt <- dt %>% filter(Length > 0)
summary(dt$Length)
min(dt$Length)
dt$Length <- dt$Length / 10 # convert to cm
dt %>% filter(Length > 120) %>% View
dt <- dt %>% filter(Length < 500) # remove those crazy length records

# LAT LONG
library(rnaturalearth)
require(sf)
proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(Longitude, Latitude, Year), coords = c("Longitude", "Latitude"), crs = 4326)
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% st_make_valid %>% 
  st_break_antimeridian(., lon_0 = 160)
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
sum(is.na(dt$Taxonomy) + 0) == 0
sum(str_detect(dt$Taxonomy, '^\\s+$') + 0) == 0

# genus level add sp
dt[str_which(dt$Taxonomy, 'sp.$'),]
dt$Taxonomy <- str_replace_all(dt$Taxonomy, 'sp.$', 'sp')

# misspellings check, but not taxonomic cleaning
# check that genera are genera, not family names (-idae/eae)
dt$Taxonomy %>% str_detect(., 'idae$|ini$') %>% sum(. + 0)
dt[str_which(dt$Taxonomy, 'idae$|ini$'),] %>% View # not a family name. just a species

# Check the species records
sort(unique(dt$Taxonomy))

dt[str_which(dt$Taxonomy, 'unidentified *'),] %>% unique() %>% View
dt$Species <- dt$Taxonomy # make a copy
dt$Species[str_which(dt$Species, 'idae\\sunidentified')]
dt$Family[str_which(dt$Species, 'idae\\sunidentified')] <- word(dt$Species[str_which(dt$Species, 'idae\\sunidentified')], 1)
dt$Family[str_which(dt$Species, 'Parrotfish|Scaridae|parrotfish')] <- 'Labridae'
dt$Species <- str_replace(dt$Species, 'Parrotfish|Scaridae|parrotfish', 'Scarini')

# put the families in their place!
dt$Species[str_which(dt$Species, 'idae\\sunidentified')] <- str_replace(dt$Species[str_which(dt$Species, 'idae\\sunidentified')], 'unidentified\\s*', 'sp')
dt <- dt %>% filter(!str_detect(dt$Species, '^Osteichthyes')) # remove the record that is too uncertain
dt <- dt %>% filter(!str_detect(dt$Species, 'No fish observed')) # remove zeroes
dt[str_which(dt$Family, 'ormes'),] # check that family column doesn't have any more non-family records
dt$Species <- str_replace(dt$Species, 'sp.\\s', 'sp')

# check again
dt  %>% distinct(Family, Species) %>% 
  arrange(Family, Species) %>% View
dt$Species <- str_remove(dt$Species, '(U|u)nidentified\\s')
dt %>% filter(str_detect(Species, '\\d')) %>% distinct(Species) %>% View
dt$Species <- str_replace(dt$Species, '\\s8$', '\\ssp8')
dt %>% filter(str_detect(Species, '\\(cf')) %>% distinct(Species) %>% View
# I'll just remove these cfs... 
dt$Species <- str_remove(dt$Species, '\\s\\(cf\\)')

dt$Genus <- word(dt$Species, 1) # migrate the first word in the species column to genus
dt$Species <- word(dt$Species, start=2, end=-1)
sort(unique(dt$Genus))
dt$Genus[str_which(dt$Genus, 'idae$')] <- NA
dt$Taxonomy <- NULL

# Prepare raw data --------------------------------------------------------
with(dt, paste(Site, Habitat, Transect, Year, Month, Day, sep = "_")) %>% n_distinct
with(dt, paste(Site, Habitat, Transect, Year, Month, sep = "_")) %>% n_distinct
with(dt, paste(Site, Habitat, Transect, Year, sep = "_")) %>% n_distinct # month and day matter

dt_final <- dt %>% mutate(DataID = dataid,
                          SamplingEventID = paste(Site, Habitat, Transect, Year, Month, Day, sep = "_")) %>% # not including swath because the true sample is both swathes together
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
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long
