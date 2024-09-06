
# ProducDiv
# Data checks, gridding, sample-based rarefaction


# SET UP ===========================================================================================
require(dplyr)
require(stringr)
require(readr)
require(ggplot2)
require(sf)
require(rnaturalearth)
require(dggridR)
require(beepr)

# global parameters
G = 0.04 # grid cell size in sq km 

# load in data
dt_assemblage <- read_csv('./src/macro_reef_fish_clean.csv')
dt_meta <- read_csv('./src/macro_reef_fish_meta.csv')

# GRID DATA ========================================================================================
# not all datasets use permanent locations
# there is no way to determine time-series based on coordinates with multiple sampling regimes
# so we need to grid by a reasonable resolution
# sampling events per year per grid cell per dataset = time series

proj = '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
# construct grid
grid <- dgconstruct(area = G, metric = TRUE, resround = 'nearest')

summary(dt_assemblage)

# isolate out coordinates
dt_coord <- dt_assemblage %>% distinct(DataID, Latitude, Longitude)
grid_match <- dgGEO_to_SEQNUM(grid, dt_coord$Longitude, dt_coord$Latitude)
length(grid_match$seqnum) == nrow(dt_coord)

dt_coord$gridcell <- grid_match$seqnum # match grid with coordinates
# merge this in
# View(dt_coord)
dt_assemblage <- left_join(dt_assemblage, dt_coord, by = c("DataID", "Latitude", "Longitude"))
summary(dt_assemblage)

# IDENTIFY TIME SERIES =============================================================================

# pull out rareIDs that are time-series
timeseries <- dt_assemblage %>% group_by(DataID, gridcell) %>%
  summarise(yearN = max(Year) - min(Year) + 1) %>% 
  filter(yearN >= 2) %>% distinct(DataID, gridcell) %>% 
  mutate(rareID = paste(DataID, gridcell, sep = '_'))

# length(timeseries) # 5108 study*gridcells

# split the dataset between single and multiple year studies

dt_single <- dt_assemblage[-which(dt_assemblage$DataID %in% timeseries$DataID & dt_assemblage$gridcell %in% timeseries$gridcell),]
dt_time <- inner_join(dt_assemblage, timeseries, by = c('DataID', 'gridcell'))
dt_time$rareID <- with(dt_time, paste(DataID, gridcell, sep = "_"))

# full unrarefied raw data
# dt_single
# dt_time
# dt_assemblage

# LOCATION FIELD ===================================================================================

dt_meta %>% filter(DataID %in% dt_time$DataID) %>% 
  select(DataID, SamplingEventStructure) %>% View
View(dt_time %>% group_by(DataID, rareID, Year) %>% 
  mutate(samplingN = n_distinct(SamplingEventID)) %>% ungroup %>% 
  group_by(DataID, rareID) %>% 
  summarise(coordN = n_distinct(paste0(Latitude, Longitude)), mean_sampling_yr = mean(samplingN)))
dt_time$Location <- NA

dt_time <- mutate(dt_time, 
                     Location = case_when(.default = NA,
                                          DataID == '001_dyiigurra_macro' ~ str_extract(SamplingEventID, '(?<=\\_)[A-Z0-9a-z\\_]+(?=\\_)'),
                                          DataID == '011_quimbayo_timefish_scuba' ~ str_remove(SamplingEventID, '\\_\\d{4}\\_\\d+$'),
                                          DataID == '020_noaa_micronesia' ~ str_extract(SamplingEventID, '(?<=^\\d{4}\\_)[A-Z0-9a-z\\-]+(?=\\_)'),
                                          DataID == '031_crane_ulithi' ~ str_extract(SamplingEventID, '(?<=\\_)[A-Za-z]+(?=\\_)'),
                                          DataID == '054_RLS' ~ str_extract(SamplingEventID, '^[A-Z0-9]+(?=\\_)'),
                                          DataID == '122_SIFkoester_aldabra' ~ str_extract(SamplingEventID, '^[A-Z0-9]+(?=\\_)'),
                                          str_detect(DataID, 'SERF') ~ str_extract(SamplingEventID, '(?<=\\_)\\d{3,4}(?=\\_)'),
                                          str_detect(DataID, 'richardson') ~ str_extract(SamplingEventID, '(?<=\\_)[A-Z0-9]{3,}(?=\\_)'),
                                          str_detect(SamplingEventID, '^anequim') ~ str_extract(SamplingEventID, '(?<=\\d{4}\\_)[a-z]+'),
                                          str_detect(DataID, 'Aqaba') ~ str_extract(SamplingEventID, '^[A-Za-z\\_](?=\\_\\d)'),
                                          DataID == 'BT45_moorea_lter' ~ str_extract(SamplingEventID, '^\\d\\_[A-Z]+')
                     ))

dt_time %>% filter(is.na(Location) == F) %>% distinct(DataID, SamplingEventID, Location) %>% 
  group_by(DataID) %>% sample_n(20) %>% View

dt_time %>% filter(is.na(Location) == F) %>% group_by(rareID) %>% 
  summarise(locationN = n_distinct(Location)) %>% 
  filter(locationN >1) %>% View

# RAREFY TIME SERIES ===============================================================================

# data structure now: dataset > grid cell/rareID > year > sampling event
# if there's any discrepancy in sampling effort within the grid cell
# find minimum number of sampling events in a year across the time series
# downsample the time-series years to match that number
# consistent number of sampling events

# this is only run on time series studies

set.seed(24)

# rarefy each gridcell for each study by lowest number of sampling events
for (D in 1:nrow(timeseries)) { # per rareID

    dataid <- timeseries$rareID[D]
    loc <- dt_time %>% filter(rareID == dataid) %>% distinct(Location) %>% pull
    
    if (length(loc) > 1) { # per location if that rareID needs it
      for (L in 1:length(loc)) {
        if (is.na(loc[L])) next
        # get the lowest common denominator across years rN
        rN <- dt_time %>% filter(rareID == dataid, Location == loc[L]) %>%
          group_by(Year) %>% 
          summarise(eventN = n_distinct(SamplingEventID))  %>%
          pull(eventN) %>% min() 
        
        # randomly select the sampling events
        events <- dt_time %>% filter(rareID == dataid, Location == loc[L]) %>%
          distinct(rareID, Location, SamplingEventID, Year) %>%
          group_by(rareID, Location, Year) %>% sample_n(rN, replace = FALSE) %>%
          pull(SamplingEventID)
        
        # check if the new object exists yet
        if (D == 1) {
          dt_rareass <- dt_time %>%
            filter(rareID == dataid, Location == loc[L], SamplingEventID %in% events)
        } else {
          dt_rareass <- dt_time %>% 
            filter(rareID == dataid, Location == loc[L], SamplingEventID %in% events) %>%
            bind_rows(dt_rareass, .)
        }
        if (exists('complete') == F) { # progress object
          complete = 1
          } else complete = complete + 1
        print(paste0("Iteration ", complete, " complete"))
      }
    } else {
      
      # get the lowest common denominator across years rN
      rN <- dt_time %>% filter(rareID == dataid) %>%
        group_by(Year) %>% 
        summarise(eventN = n_distinct(SamplingEventID))  %>%
        pull(eventN) %>% min() 
      
      # randomly select the sampling events
      events <- dt_time %>% filter(rareID == dataid) %>%
        distinct(rareID, SamplingEventID, Year) %>%
        group_by(rareID, Year) %>% sample_n(rN, replace = FALSE) %>%
        pull(SamplingEventID)
      
      # check if the new object exists yet
      if (D == 1) {
        dt_rareass <- dt_time %>%
          filter(rareID == dataid, SamplingEventID %in% events)
      } else {
        dt_rareass <- dt_time %>% 
          filter(rareID == dataid, SamplingEventID %in% events) %>%
          bind_rows(dt_rareass, .)
      }
      if (exists('complete') == F) { # progress object
        complete = 1
      } else complete = complete + 1
      print(paste0("Iteration ", complete, " complete"))
    }
    if (D == nrow(timeseries)) {print("Loop complete")}
}

# ITERATE RICHNESS =================================================================================

set.seed(24)

# number of resampling iterations
# creates sub datasets to calculate richness
iter = 25
dt_div <- dt_rareass %>% group_by(rareID, Location, SamplingEventID) %>% 
  summarise(S = n_distinct(worms_valid), iter = 1) %>% ungroup

for (D in 1:nrow(timeseries)) { # per rareID
  dataid <- timeseries$rareID[D]
  L = dt_time %>% filter(rareID == dataid) %>% pull(Location) %>% unique
  
  # if this rareID has the additional location field
  if ((is.na(L) == F) |> sum() > 0 & length(L) > 1) {
    
    for (l in 1:length(L)) { # each location
      loc = L[l]
      if (is.na(loc)) next
      # get the lowest common denominator across years rN
      rN <- dt_time %>% filter(rareID == dataid, Location == loc) %>%
        group_by(Year) %>% 
        summarise(eventN = n_distinct(SamplingEventID)) %>%
        pull(eventN) %>% min()
      
      for (i in 2:iter) { # rinse and repeat
        # randomly select the sampling events
        events <- dt_time %>% filter(rareID == dataid, Location == loc) %>%
          distinct(SamplingEventID, Year) %>%
          group_by(Year) %>% sample_n(rN, replace = FALSE) %>%
          pull(SamplingEventID)
        
        dt_div <- dt_time %>% 
          filter(rareID == dataid, Location == loc, SamplingEventID %in% events) %>%
          group_by(rareID, Location, SamplingEventID) %>% 
          summarise(S = n_distinct(worms_valid), iter = i) %>% ungroup %>%
          bind_rows(dt_div, .) # name the iteration
        
        if (exists('complete') == F) {
          complete = 1
        } else complete = complete + 1
        
      } 
    }
    
  } else {
    
    rN <- dt_time %>% filter(rareID == dataid) %>%
      group_by(Year) %>% 
      summarise(eventN = n_distinct(SamplingEventID)) %>%
      pull(eventN) %>% min()
    
    for (i in 2:iter) { # rinse and repeat
      # randomly select the sampling events
      events <- dt_time %>% filter(rareID == dataid) %>%
        distinct(SamplingEventID, Year) %>%
        group_by(Year) %>% sample_n(rN, replace = FALSE) %>%
        pull(SamplingEventID)
      
      dt_div <- dt_time %>% 
        filter(rareID == dataid, SamplingEventID %in% events) %>%
        group_by(rareID, SamplingEventID) %>% 
        summarise(S = n_distinct(worms_valid), iter = i) %>% ungroup %>%
        bind_rows(dt_div, .) # name the iteration
      
      if (exists('complete') == F) {
        complete = 1
      } else complete = complete + 1
      
    }
  } 
  print(paste0("Iteration ", complete, " complete"))
  print(summary(dt_div))
  if (D == n_distinct(dt_time$rareID)) {print('Loop complete')}
}

beep()
rm(i, rN, complete, l, L, loc, D, events, dataid, G)

# calculate richness confidence intervals
# across the iterations for each sampling event
dt_div <- dt_div %>% group_by(rareID, Location, SamplingEventID) %>% 
  summarise(S = median(S), S_lwr = quantile(S, 0.025), S_upr = quantile(S, 0.975)) %>% ungroup

summary(dt_rareass)
summary(dt_div)
dt_div %>% filter(S == 1) %>% View


# DATA CHECKS ======================================================================================
# check the rarefied, gridded time series
sum(is.na(dt_rareass$gridcell))
sum(str_detect(dt_rareass$rareID, 'NA$'))
dt_rareass <- mutate(dt_rareass, 
                     seriesID = case_when(.default = rareID,
                                          is.na(Location) == F ~ paste(rareID, Location, sep = '_')))
dt_div <- mutate(dt_div, 
                     seriesID = case_when(.default = rareID,
                                          is.na(Location) == F ~ paste(rareID, Location, sep = '_')))
# seriesID = rareID + Location

# time series durations
dt_rareass %>% group_by(seriesID) %>%
    summarise(yearN = max(Year) - min(Year) + 1) %>% View

# some after gridding turned out to be single-year
single_series <- dt_rareass %>% group_by(seriesID) %>%
  summarise(yearN = max(Year) - min(Year) + 1) %>% filter(yearN == 1) %>%
  pull(seriesID)

# migrate and remove from time series data
dt_single <- dt_rareass %>% filter(seriesID %in% single_series) %>% select(!c(rareID, seriesID, Location)) %>% 
  bind_rows(dt_single, .)
dt_rareass <- dt_rareass %>% filter(!seriesID %in% single_series)
dt_div <- dt_div %>% filter(!seriesID %in% single_series)
summary(dt_rareass)
dt_rareass %>% group_by(seriesID) %>%
  summarise(yearN = max(Year) - min(Year) + 1) %>% summary

# rarefaction check
dt_rareass %>% group_by(seriesID, Year) %>%
    summarise(samplingN = n_distinct(SamplingEventID)) %>%
    ungroup() %>% group_by(seriesID) %>%
    summarise(min_samp = min(samplingN), max_samp = max(samplingN),
    consistent = case_when(min_samp == max_samp ~ 'Y',
                           min_samp != max_samp ~ 'N')) %>%
    filter(consistent == 'N') # any rarefaction failures

# removed 036_noaa_sefsc_5916040223 because the sampling events come from multiple years
dt_rareass <- dt_rareass %>% 
  filter(!rareID %in% c('036_noaa_sefsc_5916040223', '036_noaa_sefsc_5925251511'))
dt_single <- dt_time %>% 
  filter(rareID %in% c('036_noaa_sefsc_5916040223', '036_noaa_sefsc_5925251511')) %>% 
  select(!c(Location, rareID)) %>% bind_rows(dt_single, .)
