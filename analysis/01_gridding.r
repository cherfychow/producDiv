
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
  filter(yearN >= 2) %>% distinct(DataID, gridcell)

# length(timeseries) # 5108 study*gridcells

# split the dataset between single and multiple year studies

dt_single <- dt_assemblage[-which(dt_assemblage$DataID %in% timeseries$DataID & dt_assemblage$gridcell %in% timeseries$gridcell),]
dt_time <- inner_join(dt_assemblage, timeseries, by = c('DataID', 'gridcell'))

# RAREFY TIME SERIES ===============================================================================

# data structure now: dataset > grid cell/rareID > year > sampling event
# if there's any discrepancy in sampling effort within the grid cell
# find minimum number of sampling events in a year across the time series
# downsample the time-series years to match that number
# consistent number of sampling events

# this is only run on time series studies

dt_time$rareID <- with(dt_time, paste(DataID, gridcell, sep = "_"))

set.seed(24)

# rarefy each gridcell for each study by lowest number of sampling events
for (D in 1:nrow(timeseries)) {

    dataid  <-  unique(dt_time$rareID)[D]
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
}

# number of resampling iterations
# creates sub datasets to calculate richness
iter = 25
dt_div <- dt_rareass %>% group_by(rareID, SamplingEventID) %>% 
  summarise(S = n_distinct(worms_valid), iter = 1) %>% ungroup
for (D in 1:n_distinct(dt_time$rareID)) {
    dataid  <-  unique(dt_time$rareID)[D]
    # get the lowest common denominator across years rN
    rN <- dt_time %>% filter(rareID == dataid) %>%
        group_by(Year) %>% 
        summarise(eventN = n_distinct(SamplingEventID)) %>%
        pull(eventN) %>% min() 
    for (i in 2:iter) { # rinse and repeat
        # randomly select the sampling events
        events <- dt_time %>% filter(rareID == dataid) %>%
        distinct(rareID, SamplingEventID, Year) %>%
        group_by(rareID, Year) %>% sample_n(rN, replace = FALSE) %>%
        pull(SamplingEventID)

        dt_div <- dt_time %>% 
            filter(rareID == dataid, SamplingEventID %in% events) %>%
          group_by(rareID, SamplingEventID) %>% 
          summarise(S = n_distinct(worms_valid), iter = i) %>% ungroup %>%
          bind_rows(dt_div, .) # name the iteration
    }
}  
beep()
rm(i, rN, D, events, dataid, G, dt_time)

# calculate intervals for richness
# should not vary that much unless sampling effort varies a lot
dt_div <- dt_div %>% group_by(DataID, gridcell, SamplingEventID) %>% 
  summarise(S = median(S), S_lwr = quantile(S, 0.025), S_upr = quantile(S, 0.975)) %>% ungroup

# generate rarefy ID (= DataID x grid cell)
summary(dt_rareass)



# DATA CHECKS -----------------------------------------------------------
# check the rarefied, gridded time series
# sum(is.na(dt_rareass$gridcell))
# sum(str_detect(dt_rareass$rareID, 'NA$'))
# 
# summary(dt_rareass)
# # time series durations
# dt_rareass %>% group_by(rareID) %>% 
#     summarise(yearN = max(Year) - min(Year) + 1) %>% summary
# 
# dt_rareass %>% group_by(rareID) %>% 
#     summarise(yearN = max(Year) - min(Year) + 1) %>%
#     pull(yearN) %>% hist
# 
# # rarefaction check
# dt_rareass %>% group_by(rareID, Year) %>%
#     summarise(samplingN = n_distinct(SamplingEventID)) %>%
#     ungroup() %>% group_by(rareID) %>%
#     summarise(min_samp = min(samplingN), max_samp = max(samplingN),
#     consistent = case_when(min_samp == max_samp ~ 'Y', 
#                            min_samp != max_samp ~ 'N')) %>%
#     filter(consistent == 'N') # any rarefaction failures
