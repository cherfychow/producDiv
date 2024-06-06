

# Curation Script ---------------------------------------------------------

# Dataset: Benthic and reef fish surveys from the Kingdom of Tonga in 2016-2018 (375 sites)
# Location: Tonga

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(ggplot2)
require(readr)
require(stringr)
require(lubridate)
require(readxl)

rm(list=ls())
dataid <- '060_smallhornwest_tonga'

# make sure your working directory is set before running these lines
dt <- read_delim('data_originals/060_smallhornwest_tonga/Smallhorn-West-etal_2020_rawdata.txt')
og_dt <- dt
View(dt)
str(dt)
summary(dt)

# Structure organisation --------------------------------------------------

dt <- dt %>% select(`Ind No [#]`, `l [cm]`, Family, Species, Site, `Transect no`, `Length [m]`, `Width [m]`, `Date/Time (Start)`, Latitude, Longitude, `Depth water [m]`)

colnames(dt)[c(1,2,12)] <- c('Abundance', 'Length', 'Depth')

# dates need to be coerced
dt$Day <- day(dt$`Date/Time (Start)`) %>% as.integer()
dt$Month <- month(dt$`Date/Time (Start)`) %>% as.integer()
dt$Year <- year(dt$`Date/Time (Start)`) %>% as.integer
dt$`Date/Time (Start)` <- NULL

summary(dt) # check transect width consistency...
dt %>% filter(is.na(`Width [m]`)) %>% View # at the same site in 2017
dt %>% filter(Site == 'Matamaka MS1', `Transect no` == 365) %>% View
dt %>% filter(`Transect no` == 365) %>% View
dt <- dt %>% filter(`Transect no` != 365) # filter out this transect
dt$Abundance <- dt$Abundance / dt$`Width [m]` * 5

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 80985 rows
str(dt) # check structure


# logical?
summary(dt[,c('Year', 'Month', 'Day')]) # yes

summary(dt)
# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. NA
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
dt %>% filter(is.na(Length))
dt %>% filter(is.na(Abundance)) %>% View # blanks in biomass? biomass more complete
dt$Abundance[which(is.na(dt$Abundance))] <- 1 / dt$`Width [m]`[which(is.na(dt$Abundance))] * 5 # there is length data, so put 1 for the occurrence at least
min(dt$Abundance, na.rm=T) # no zeroes? Y
min(dt$Length, na.rm=T) # no zeroes? Y
dt %>% filter(Length <= 0) %>% View
# a lot of unmeasurable lengths
dt %>% filter(Length > 100) %>% View
# only one unreasonable record
dt <- dt[-which(dt$Species == 'Zebrasoma veliferum' & dt$Length == 211), ]

# LAT LONG
library(rnaturalearth)
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
sum(is.na(dt$Species) + 0) == 0
sum(str_detect(dt$Species, '^\\s+$') + 0) == 0

sort(unique(dt$Species))

replace = c('sp\\.' = 'sp',
            '\\?' = ' ',
            'arueus' = 'aureus',
            'Red with tailspot' = 'sp1',
            'Cheiliodipterus' = 'Cheilodipterus',
            'Chormis' = 'Chromis',
            'Juvenile scarid' = 'Scaridae sp',
            'dicki$' = 'dickii', 
            'Remora' = 'Echeneidae',
            'Thalossoma' = 'Thalassoma')

dt$Species <- str_replace_all(dt$Species, replace)
dt$Genus <- word(dt$Species, 1) %>% str_to_title
sort(unique(dt$Genus))
dt$Genus[str_which(dt$Genus, 'idae$')] <- NA
dt$Species <- word(dt$Species, 2)

# Prepare raw data --------------------------------------------------------

# check sampling structure
dt %>% group_by(Site) %>% summarise(transectN = n_distinct(`Transect no`), yearN = n_distinct(Year)) %>% View
# not actually a time series

dt_final <- dt %>% mutate(DataID = dataid,
                          SamplingEventID = paste(Site, `Transect no`, sep = "_")) %>% 
  mutate(Year = NA, Month = NA, Day = NA) %>% # remove since it's not a time-series
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
