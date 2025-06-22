
# ProducDiv
# matching up with human gravity data

require(tidytable)
require(sf)



# read in the global grid layer from PNAS Cinner et al
grav <- st_read('src/market_gravity/GravityFullPolygons.shp')
summary(grav)
str(grav)
grav <- st_make_valid(grav)

# goes through a similar process as SST matching to find an overlapping or nearest neighbour feature
# 10 x 10 km resolution

sf_use_s2(FALSE)
grav_match <- st_nearest_feature(grid_sf, grav) # gives an index of grav for each grid cell of ours
grid_sf$hgravity <- st_drop_geometry(grav)$Grav_citie[grav_match]
grid_sf$grav_nc <- st_drop_geometry(grav)$Grav_NC[grav_match]
grid_sf$grav_ID <- st_drop_geometry(grav)$reef_ID[grav_match]
grid_sf$hgravity[is.na(grid_sf$hgravity)] <- grid_sf$grav_nc[is.na(grid_sf$hgravity)] # use grav NC if cumulative market gravity isn't available
grid_sf <- grid_sf %>% rename(gridcell = seqnum)

# biogeography regions
ecoR <- st_read('src/biogeo_regions.shp')
st_crs(ecoR) <- 4326
ecoR <- ecoR %>% st_make_valid(.) %>% st_break_antimeridian(., lon_0 = 160)
biogeo_match <- st_nearest_feature(grid_sf, ecoR)
grid_sf$biogeo <- st_drop_geometry(ecoR)$region[biogeo_match]


# Compile time series variables --------------------------------------------------------------

# compile a data table of predictor variables for each time series sampling event
vars_rare <- dt_rareass %>% 
  distinct(DataID, seriesID, SamplingEventID, gridcell, Year) %>% 
  left_join(., dt_div, by = colnames(dt_div)[-5])

# productivity
vars_rare <- dt_prod_rare %>% select(ends_with('ID'), starts_with('prod_kg')) %>% 
  left_join(vars_rare, ., by = c('DataID', 'seriesID', 'SamplingEventID'))

vars_rare <- left_join(vars_rare, st_drop_geometry(grid_sf), by = 'gridcell')

# Compile single year variables --------------------------------------------------------------

# observed richness
vars_single <- dt_single %>% 
  distinct(DataID, seriesID, SamplingEventID, gridcell, Year) %>% 
  left_join(., dt_div_single, by = c('SamplingEventID', 'DataID', 'seriesID'))

# productivity
vars_single <- dt_prod_single %>% select(DataID, seriesID, SamplingEventID, starts_with('prod_kg')) %>% 
  left_join(vars_single, ., by = c('DataID', 'seriesID', 'SamplingEventID'), relationship = 'one-to-one')

# human gravity
vars_single <- left_join(vars_single, st_drop_geometry(grid_sf), by = 'gridcell')

vars_rare <- vars_rare %>% ungroup
vars_single <- vars_single %>% ungroup


# Match MPA status =================================================================================

