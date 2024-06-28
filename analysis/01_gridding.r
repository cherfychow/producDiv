
# ProducDiv
# Data checks, gridding, sample-based rarefaction


# SET UP ---------------------------------------------------------
require(dplyr)
require(stringr)
require(sf)
require(readr)
require(ggplot2)
require(rnaturalearth)
require(dggridR)
require(rfishbase)

# global analysis parameters
G = 0.04 # grid cell size in sq km

# load in data
tempdt <- read_csv('data_outputs/001_dyiigurra_macro_obs.csv')
tempdt2 <- read_csv('data_outputs/45_moorea_lter_obs.csv')
dt_assemblage <- bind_rows(tempdt, tempdt2)
rm(tempdt, tempdt2)

View(dt_assemblage)

# CHECK DATA ------------------------------------------------------
# global taxonomy check

dt_assemblage %>% distinct(Family, Genus, Species) %>%
    arrange(Family, Genus, Species)%>% View
dt_assemblage %>% distinct(Genus, Species) %>% 
    filter(str_detect(Genus, 'idae')) %>% View

# NAs for non specific uncertain IDs
dt_assemblage$Genus[str_detect(dt_assemblage$Genus, 'idae$') & dt_assemblage$Species == 'sp'] <- NA
dt_assemblage$Genus[is.na(dt_assemblage$Genus) & dt_assemblage$Species == 'sp'] <- NA

# GRID DATA -------------------------------------------------------

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
dt_assemblage <- left_join(dt_assemblage, dt_coord, by = c("Latitude", "Longitude", 'DataID'))

# visual checks
# mostly years and sampling events per cell

cellcount <- dt_assemblage %>% group_by(gridcell, DataID) %>% 
    summarise(yearN = n_distinct(Year), eventN = n_distinct(SamplingEventID))
grid_count <- dgcellstogrid(grid, cellcount$gridcell) %>% 
    left_join(., cellcount, by = c('seqnum' = 'gridcell'))

world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% 
st_make_valid() %>% st_break_antimeridian(., lon_0 = 160)

# Dyiigurra
lon = c(min(dt_assemblage$Longitude[str_detect(dt_assemblage$DataID, 'macro')]),
        max(dt_assemblage$Longitude[str_detect(dt_assemblage$DataID, 'macro')]))
lat = c(min(dt_assemblage$Latitude[str_detect(dt_assemblage$DataID, 'macro')]),
        max(dt_assemblage$Latitude[str_detect(dt_assemblage$DataID, 'macro')]))
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = grid_count, aes(fill = yearN), shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(xlim = lon, ylim = lat) +
  theme_minimal()

# Mo'orea
lon = c(min(dt_assemblage$Longitude[str_detect(dt_assemblage$DataID, 'moorea')]),
        max(dt_assemblage$Longitude[str_detect(dt_assemblage$DataID, 'moorea')]))
lat = c(min(dt_assemblage$Latitude[str_detect(dt_assemblage$DataID, 'moorea')]),
        max(dt_assemblage$Latitude[str_detect(dt_assemblage$DataID, 'moorea')]))
ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = grid_count, aes(fill = yearN), shape = 21, size = 2) +
  guides(fill = 'none') + 
  coord_sf(xlim = lon, ylim = lat) +
  theme_minimal()

rm(world, lon, lat, cellcount, grid_match)
View(grid_count)
rm(grid_count)

# RAREFY -----------------------------------------------------------

# dataset > grid cell > year > sampling event
# grid cell becomes the proxy for location
# sampling events per year has to be equal

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
rm(rN, D, events, dataid, G)
