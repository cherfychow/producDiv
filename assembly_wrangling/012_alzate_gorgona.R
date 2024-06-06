

# Curation Script ---------------------------------------------------------

# Dataset: A cause for hope: largely intact coral-reef communities with high reef-fish biomass in a remote Indonesian island group
# Location: Lucipara, Indonesia
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(readxl)
require(sf)
require(tidyverse)
require(lubridate)

rm(list=ls())

# Maarten already did formatting so this is just a cursory check for formatting etc.

# make sure your working directory is set before running these lines
dt <- read_excel('data_originals/012_template_Alzate_Gorgona_AA.xlsx', sheet = 2)
og_dt <- dt
dataid <- '012_alzate_gorgona'
dt$StudyID <- dataid
dt <- dt[1:579,]
dt$Plot <- NULL

# Structure organisation + check --------------------------------------------------

summary(dt)
str(dt)
nrow(dt) # check dimensions

sort(unique(dt$Month))
replace <- c('December' = '12',
             'November' = '11',
             'October' = '10')
dt$Month <- str_replace_all(dt$Month, replace) %>% as.numeric

# Secondary fields relating to sample must be factors or integers? Y
# Taxonomic fields must be characters or factors? Y


# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields. Y
sort(unique(dt$Length))
summary(dt)
# YEAR MONTH DAY, all sensical

sort(unique(dt$Species))
dt$Genus <- word(dt$Species, 1)
dt$Species <- word(dt$Species, 2, -1)

dt$DepthElevation[str_which(dt$DepthElevation, 'meter')] <- str_remove(dt$DepthElevation[str_which(dt$DepthElevation, 'meter')], '\\smeters')
feet <- dt$DepthElevation[str_which(dt$DepthElevation, 'feet')] %>% str_remove(., '\\sfeet') %>% as.numeric
dt$DepthElevation[str_which(dt$DepthElevation, 'feet')] <- feet * 0.3048
dt$DepthElevation <- as.numeric(dt$DepthElevation)
dt$Latitude <- 2.97
dt$Longitude <- -78.18

# Prepare raw data --------------------------------------------------------
dt_final <- dt %>% mutate(SampleDescription = paste(SampleDescription, Year, Month, Day, sep = "_")) %>% 
  arrange(Year, Month, Day, Genus, Species) %>% 
  rename(Depth = DepthElevation, SamplingEventID = SampleDescription, DataID = StudyID) %>% 
  relocate(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID) 

lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length_upr = c(5, 10, seq(15, 55, by = 5))
lengths$Length_lwr = c(1, 5, seq(10, 50, by = 5))
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)

n_distinct(dt_final$SamplingEventID)

# Export final ------------------------------------------------------------

readr::write_csv(dt_final, 'data_outputs/012_alzate_gorgona.csv')

