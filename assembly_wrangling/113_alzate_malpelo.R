

# Set up ------------------------------------------------------------------
# load the necessary packages
require(readxl)
require(sf)
require(tidyverse)
require(lubridate)

rm(list=ls())

# make sure your working directory is set before running these lines
dt <- read_excel('data_originals/Fish surveys_MMPalacios.xlsx', sheet = 3)
dataid <- '113_alzate_malpelo'



# Structure organisation + check --------------------------------------------------

dt <- pivot_longer(dt, cols = ends_with('cm'), names_to = "sizeclass", values_to = "Abundance")
dt <- dt %>% filter(is.na(Abundance) == F)
dt$sizeclass <- str_remove(dt$sizeclass, '\\s(?=\\-)')
dt$Length_lwr <- str_extract(dt$sizeclass, '^\\d+') %>% as.numeric
dt$Length_upr <- str_extract(dt$sizeclass, '\\d+(?=cm)') %>% as.numeric
dt$Length_lwr[which(dt$sizeclass == '>60cm')] <- 60
dt$Length <- with(dt, Length_lwr + (Length_upr - Length_lwr)/2)
dt$Length_upr[which(dt$Length_upr == 60)] <- NA

summary(dt)
str(dt)
nrow(dt) # check dimensions

# quick tidy
dt <- dt %>% select(!c(Location, SurveyDuration_min, Visibility_m, Diver, ReefZone,
                       SurveyArea_m2, FishDens_100m2, sizeclass, Abundance_total))
dt$Abundance <- (dt$Abundance / dt$...27 ) * 100 # density / 100 m2
dt$...27 <- NULL

# check if dates are replaceable with site and replicates
dt %>% mutate(Day = day(Date), Month = month(Date), Year = year(Date)) %>% 
  distinct(Year, Month, Day, Reef, Site, Rep) %>% 
  arrange(Site, Rep) %>% View

# ok dates do matter for sampling events

dt <- dt %>% mutate(Day = day(Date), Month = month(Date), Year = year(Date))
dt$Reef <- NULL # all at the same
dt$Date <- NULL

# Secondary fields relating to sample must be factors or integers? Y
# Taxonomic fields must be characters or factors? Y


# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields. Y
sum(is.na(dt$Abundance))
sum(dt$Abundance == '')
summary(dt$Abundance)

sort(unique(dt$Length))
summary(dt$Length)
sum(is.na(dt$Length))
sum(dt$Length == '')

sort(unique(dt$Species))
dt$Species <- str_replace(dt$Species, 'leweini', 'lewini')
dt$Species <- str_replace(dt$Species, 'sp\\.', 'sp')
dt$Genus <- word(dt$Species, 1)
dt$Genus <- str_replace(dt$Genus, 'Cantherines', 'Cantherhines')
dt$Genus <- str_replace(dt$Genus, 'Johnrandalia', 'Johnrandallia')
sort(unique(dt$Genus))
dt$Species <- word(dt$Species, 2)

dt %>% distinct(Genus, Species) %>% arrange(Genus, Species) %>% View


# Prepare raw data --------------------------------------------------------
dt_final <- dt %>% mutate(SamplingEventID = paste(Site, Rep, Day, sep = "_"),
                          DataID = dataid, Depth = NA, Family = NA, 
                          Latitude = NA, Longitude = NA) %>% 
  arrange(Year, Month, Day, Site, Rep, Genus, Species) %>% 
  select(Abundance, Length, Length_lwr, Length_upr, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID, SurveyMethod) 

n_distinct(dt_final$SamplingEventID)
unique(dt_final$SurveyMethod)

# Export final ------------------------------------------------------------

dt_final %>% filter(SurveyMethod == 'EDGAR') %>% select(!SurveyMethod) %>% 
  readr::write_csv(., paste0('data_outputs/', dataid, '_RLS_obs.csv'))

dt_final %>% filter(SurveyMethod != 'EDGAR') %>% select(!SurveyMethod) %>% 
  readr::write_csv(., paste0('data_outputs/', dataid, '_SIMAC_obs.csv'))
