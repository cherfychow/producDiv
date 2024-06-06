

# Curation Script ---------------------------------------------------------

# Dataset: Pacific Reef Assessment and Monitoring Program: Rapid Ecological Assessments of Fish Large-Area Stationary Point Count Surveys (SPC) at Coral Reef Sites across the Pacific Ocean from 2000 to 2007
# Location: Pacific Islands
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(dplyr)
require(stringr)
require(lubridate)
require(readxl)
require(sf)
require(rnaturalearth)
require(ggplot2)

rm(list=ls())
dataID <- '122_SIFkoester_aldabra'

dt <- read_excel('data_originals/122_SIFkoester_aldabra/Data_extended_surveys.xlsx')
og_dt <- dt
names(dt)
# get more site information from OBIS
sites <- readr::read_delim('data_originals/122_SIFkoester_aldabra/occurrence.txt')
sites <- sites %>% distinct(fieldNumber, decimalLatitude, decimalLongitude) %>% 
  filter(is.na(fieldNumber) == F)

dt <- dt %>% select(Amount, Actual_size, Family, Species, Transect, Site, Survey_Date, Year)
dt <- sites %>% rename(Latitude = decimalLatitude, Longitude = decimalLongitude) %>% 
  left_join(dt, ., by = c('Site' = 'fieldNumber'))
# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 3402 records
str(dt)
# Abundance and/or biomass, latitude and longitude numeric? Y

# Dates
dt$Day <- day(dt$Survey_Date)
dt$Month <- month(dt$Survey_Date)
dt$Survey_Date <- NULL
summary(dt)

# Secondary fields relating to sample must be factors or integers?

# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Amount)+0)
min(dt$Amount, na.rm=T) # no zeroes? Y

min(dt$Actual_size, na.rm=T) # no zeroes? Y
dt %>% filter(is.na(Actual_size))

dt %>% filter(Actual_size > 100) %>% arrange(desc(Actual_size)) %>% View 

# YEAR MONTH DAY
# already checked

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
detach('package:sf')
detach('package:ggplot2')
detach('package:rnaturalearth')

# Taxonomic field check ---------------------------------------------------

dt$Species <- str_replace(dt$Species, '\\.', ' ')
sum(is.na(dt$Species)) # number of NAs 
sum(dt$Species == '') # or blanks

sort(unique(dt$Species))
dt %>% filter(str_detect(Species, 'other')) %>% distinct(Family, Species)
dt$Species[str_which(dt$Species, '^other')] <- NA

dt$Genus <- word(dt$Species, 1)
dt$Species <- word(dt$Species, 2)
sort(unique(dt$Genus))

# Prepare raw data --------------------------------------------------------

# check sampling structure

dt %>% group_by(Site) %>% summarise(transectN = n_distinct(Transect))
dt %>% distinct(Site, Year) %>% arrange(Site)

dt <- dt %>% filter(Site != 'ARM07')
dt$Depth <- str_extract(dt$Transect, '(?<=\\_)\\d+(?=m)') %>% as.numeric

dt_final <- dt %>% mutate(DataID = dataID,
                          SamplingEventID = paste(Transect, Year, sep = "_")) %>% 
  rename(Abundance = Amount, Length = Actual_size) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID) %>% 
  arrange(Year, SamplingEventID, Genus, Species)

# Export final ------------------------------------------------------------

readr::write_csv(dt_final, paste0('data_outputs/', dataID, '_obs.csv'))


