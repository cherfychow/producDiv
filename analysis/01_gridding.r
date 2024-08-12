
# ProducDiv
# Data checks, gridding, sample-based rarefaction


# SET UP ---------------------------------------------------------
require(dplyr)
require(stringr)
require(readr)
require(rfishbase)
require(taxize)
require(ggplot2)
require(sf)
require(rnaturalearth)
require(dggridR)
require(beepr)


# global parameters
G = 0.04 # grid cell size in sq km

# GRID DATA --------------------------------------------------------------------
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
dt_coord$gridcell <- grid_match$seqnum # match with coordinates
# merge this in
View(dt_coord)
dt_assemblage <- left_join(dt_assemblage, dt_coord, by = c("DataID", "Latitude", "Longitude"))
summary(dt_assemblage)


# RAREFY -----------------------------------------------------------

# data structure now: dataset > grid cell > year > sampling event
# if there's any discrepancy in sampling effort within the grid cell
# find minimum number of sampling events in a year across the time series
# downsample the time-series years to match that number
# consistent number of sampling events

set.seed(24)

# rarefy based on sampling events per year
for (D in 1:n_distinct(dt_assemblage$DataID)) {

    dataid  <-  unique(dt_assemblage$DataID)[D]
    # get the lowest common denominator across years rN
    rN <- dt_assemblage %>% filter(DataID == dataid)  %>%
        group_by(DataID, gridcell, Year) %>% 
        summarise(eventN = n_distinct(SamplingEventID))  %>%
        pull(eventN) %>% min() 

    # randomly select the sampling events
    events <- dt_assemblage %>% filter(DataID == dataid) %>%
    distinct(gridcell, Year, SamplingEventID) %>%
    group_by(gridcell, Year) %>% sample_n(rN, replace = FALSE) %>%
    pull(SamplingEventID)

    # check if the new object exists yet
    if (exists('dt_rareass') == F) {
        dt_rareass <- dt_assemblage %>% 
            filter(DataID == dataid, SamplingEventID %in% events)
    } else {
        dt_rareass <- dt_assemblage %>% 
            filter(DataID == dataid, SamplingEventID %in% events) %>%
            bind_rows(dt_rareass, .)
    }
}

# number of resampling iterations
# creates sub datasets
dt_rareass$iter <- 1 # if you choose to do multiple resampling iterations
iter = 10

for (D in 1:n_distinct(dt_assemblage$DataID)) {
    dataid  <-  unique(dt_assemblage$DataID)[D]
    # get the lowest common denominator across years rN
    rN <- dt_assemblage %>% filter(DataID == dataid)  %>%
        group_by(DataID, gridcell, Year) %>% 
        summarise(eventN = n_distinct(SamplingEventID))  %>%
        pull(eventN) %>% min() 
    for (i in 1:(iter - 1)) { # rinse and repeat
        # randomly select the sampling events
        events <- dt_assemblage %>% filter(DataID == dataid) %>%
        distinct(gridcell, Year, SamplingEventID) %>%
        group_by(gridcell, Year) %>% sample_n(rN, replace = FALSE) %>%
        pull(SamplingEventID)

        dt_rareass <- dt_assemblage %>% 
            filter(DataID == dataid, SamplingEventID %in% events) %>%
            mutate(iter = i+1) %>% bind_rows(dt_rareass, .) # name the iteration
    }
}  
beep()
rm(i, rN, D, events, dataid, G)
# generate rarefy ID (= DataID x grid cell)
summary(dt_rareass)
dt_rareass$rareID <- with(dt_rareass, paste(DataID, gridcell, sep = "_"))


# DATA CHECKS -----------------------------------------------------------
# check the rarefied, gridded time series
sum(is.na(dt_rareass$gridcell))
sum(str_detect(dt_rareass$rareID, 'NA$'))

summary(dt_rareass)
# time series durations
dt_rareass %>% group_by(rareID) %>% 
    summarise(yearN = max(Year) - min(Year) + 1) %>% summary

dt_rareass %>% group_by(rareID) %>% 
    summarise(yearN = max(Year) - min(Year) + 1) %>%
    pull(yearN) %>% hist

# rarefaction check
dt_rareass %>% group_by(rareID, Year) %>%
    summarise(samplingN = n_distinct(SamplingEventID)) %>%
    ungroup() %>% group_by(sxrareID) %>%
    summarise(min_samp = min(samplingN), max_samp = max(samplingN),
    consistent = case_when(min_samp == max_samp ~ 'Y', 
                           min_samp != max_samp ~ 'N')) %>%
    filter(consistent == 'N') # any rarefaction failures

detach(package:taxize)