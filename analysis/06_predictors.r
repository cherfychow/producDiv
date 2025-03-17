
# ProducDiv
# matching up with human gravity data

require(tidytable)
require(sf)



# read in the global grid layer from PNAS Cinner et al
grav <- st_read('src/cinneretal_gravity.shp')
summary(grav)
str(grav)
grav <- st_make_valid(grav)

# goes through a similar process as SST matching to find an overlapping or nearest neighbour feature
# 10 x 10 km resolution

sf_use_s2(FALSE)
grav_match <- st_nearest_feature(grid_sf, grav) # gives an index of grav for each grid cell of ours


# Compile time series variables --------------------------------------------------------------

# compile a data table of predictor variables for each time series sampling event
vars_rare <- dt_rareass %>% 
  distinct(DataID, seriesID, SamplingEventID, gridcell, Year, Latitude, Longitude) %>% 
  group_by(pick(DataID:Year)) %>% 
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude)) %>% 
  # observed richness
  left_join(., dt_div, by = colnames(dt_div)[-5])

# productivity
vars_rare <- dt_prod %>% select(ends_with('ID'), starts_with('prod_kg')) %>% 
  left_join(vars_rare, ., by = c('DataID', 'seriesID', 'SamplingEventID'))

# human gravity
grid_sf$hgravity <- st_drop_geometry(grav)$Grav_tot[grav_match]
# grid_sf <- grid_sf %>% rename(gridcell = seqnum)

vars_rare <- left_join(vars_rare, st_drop_geometry(grid_sf), by = 'gridcell')

# Compile single year variables --------------------------------------------------------------

# observed richness
vars_single <- dt_sp %>% select(taxaID, worms_valid) %>% 
  left_join(dt_single[,1:17], ., by = 'taxaID')
vars_single <- vars_single %>% 
  distinct(DataID, seriesID, SamplingEventID, gridcell, Year, Latitude, Longitude, worms_valid) %>% 
  group_by(pick(DataID:Year)) %>% 
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude), S = n_distinct(worms_valid)) %>% ungroup

# productivity
vars_single <- dt_prod_single %>% select(DataID, seriesID, SamplingEventID, starts_with('prod_kg')) %>% 
  left_join(vars_single, ., by = c('DataID', 'seriesID', 'SamplingEventID'), relationship = 'one-to-one')

# human gravity
vars_single <- left_join(vars_single, st_drop_geometry(grid_sf), by = 'gridcell')
