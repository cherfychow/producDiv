

# Setup
require(tidyverse)
require(lubridate)
require(readxl)

rm(list=ls())
dataid <- '110_brandl_farms_panama'

# make sure your working directory is set before running these lines

dt <- read_excel('data_originals/110_brandl_farms/farms_alldata.xlsx', sheet = 3)
colnames(dt)
View(dt)

# two methods here, need to split

dt<- dt %>% select(-c(surfacearea, sl, genspe))

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 646 records
str(dt) # check structure
summary(dt)

dt <- dt %>% rename(Length = tl) %>%
  mutate(Abundance = 1, .before = Length, Latitude = 9.35143, 
  Longitude = (-82.2571)) # add coordinates from Brandl et al. 2017

# Secondary fields such as elevation or depth will normally be numeric unless in a treatment format. Y
# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES/LENGTH
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$Length)) == 0
summary(dt$Length)
dt$Length <- dt$Length / 10 # convert to cm

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
dt %>% distinct(genus, species) %>% View
dt$species <- str_remove(dt$species, '\\.')

# everything is ok!

# Prepare raw data --------------------------------------------------------
colnames(dt) <- str_to_title(colnames(dt))
dt_final <- dt %>% group_by_at(vars(-Abundance)) %>% 
  summarise(Abundance = sum(Abundance)) %>% ungroup %>% #aggregate by length
  mutate(DataID = dataid, SamplingEventID = str_remove(Site, '\\s'),
         Year = 2016, Month = NA, Day = NA, Family = NA, Depth = NA) %>% 
  select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)

View(dt_final)

# Export final ------------------------------------------------------------

write_csv(dt_final, paste0('data_outputs/', dataid, '_obs.csv'))

