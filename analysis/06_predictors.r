
# ProducDiv
# matching up with human gravity data

require(tidytable)
require(sf)
require(terra)
require(stringr)

# read in the global grid layer from PNAS Cinner et al
grav <- st_read('../src/market_gravity/GravityFullPolygons.shp')
summary(grav)
str(grav)
grav <- st_make_valid(grav)

# read in MPA implementation + protection
dt_prot <- st_read('src/protectedseas_areabased.geojson')

dt_prot %>% select(level_fishing_protection, restrictions, effective_from, year_established) %>%
  head(20) %>% View

dt_prot <- dt_prot %>% dplyr::select(ogc_fid, level_fishing_protection, year_established)

# ==================================================================================================
#  Match up MPA status to surveys
# ==================================================================================================
# trim and make valid
sf_use_s2(FALSE)
dt_prot <- st_crop(dt_prot, xmin = -180, xmax = 180, ymin = -40, ymax = 40)

temp <- dt_assemblage %>% distinct(Latitude, Longitude)
mpa_match <- temp %>% 
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>%
  st_intersects(., dt_prot, sparse = T)
mpa_match[which(lengths(mpa_match) == 0)] <- NA # assign NA to gridcells with no MPA match
which(lengths(mpa_match) == 0) # check that it did it right
which(lengths(mpa_match) > 1) # multiple match?

for (i in which(lengths(mpa_match) > 1)) {
  mpa_match[[i]] <- mpa_match[[i]][runif(1, min = 1, max = length(mpa_match[[i]]))]
} # randomly keep one

nrow(temp) == length(mpa_match)

# match up MPAs spatially first
mpa_match <- unlist(mpa_match)
temp$level_fishing_protection <- NA
temp$year_established <- NA
temp$level_fishing_protection[which(is.na(mpa_match)==F)] <- st_drop_geometry(dt_prot$level_fishing_protection[na.omit(mpa_match)])
temp$year_established[which(is.na(mpa_match)==F)] <- st_drop_geometry(dt_prot$year_established[na.omit(mpa_match)])

summary(temp)
dt_prot <- temp

# ==================================================================================================
#  Gravity matching
# ==================================================================================================

# goes through a similar process as SST matching to find an overlapping or nearest neighbour feature
# 10 x 10 km resolution
grav_match <- st_nearest_feature(grid_sf, grav) # gives an index of grav for each grid cell of ours
grid_sf$hgravity <- st_drop_geometry(grav)$Grav_citie[grav_match]
grid_sf$grav_tot <- st_drop_geometry(grav)$Grav_tot[grav_match]
grid_sf$grav_ID <- st_drop_geometry(grav)$reef_ID[grav_match]
grid_sf$hgravity[is.na(grid_sf$hgravity)] <- grid_sf$grav_nc[is.na(grid_sf$hgravity)] # use grav NC if cumulative market gravity isn't available
grid_sf <- grid_sf %>% rename(gridcell = seqnum)

# biogeography regions
ecoR <- st_read('../src/biogeographical_regions/biogeo_regions.shp')
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
vars_rare <- dt_prod_all %>% select(ends_with('ID'), starts_with('prod_kg')) %>% 
  left_join(vars_rare, ., by = c('DataID', 'seriesID', 'SamplingEventID'))

vars_rare <- left_join(vars_rare, st_drop_geometry(grid_sf), by = 'gridcell')

# Compile single year variables --------------------------------------------------------------

# observed richness
vars_single <- dt_single %>% 
  distinct(DataID, seriesID, SamplingEventID, gridcell, Year) %>% 
  left_join(., dt_div_single, by = c('SamplingEventID', 'DataID', 'seriesID'))

# productivity
vars_single <- dt_prod_all %>% select(DataID, seriesID, SamplingEventID, starts_with('prod_kg')) %>% 
  left_join(vars_single, ., by = c('DataID', 'seriesID', 'SamplingEventID'))

# human gravity
vars_single <- left_join(vars_single, st_drop_geometry(grid_sf), by = 'gridcell')

vars_rare <- vars_rare %>% ungroup
vars_single <- vars_single %>% ungroup

# add N
vars_single <- dt_single %>% group_by(DataID, seriesID, SamplingEventID) %>% summarise(N = sum(Abundance)) %>% 
  left_join(vars_single, ., by = c('DataID', 'seriesID', 'SamplingEventID'))

# add sample completeness
vars_single <- left_join(vars_single, dt_complete, by = c('DataID', 'seriesID', 'SamplingEventID'))

# MPA subset
temp <- dt_assemblage %>% distinct(seriesID, SamplingEventID, Latitude, Longitude)
temp <- temp %>% group_by(seriesID, SamplingEventID) %>% slice_sample(n = 1) %>% ungroup
vars_single_mpa <- left_join(vars_single, temp, by = c('seriesID', 'SamplingEventID'))
vars_single_mpa <- left_join(vars_single_mpa, dt_prot, by = c('Latitude', 'Longitude'))

vars_single_mpa <- dt_assemblage %>% select(seriesID, SamplingEventID, Year, Month) %>%
  distinct(seriesID, SamplingEventID, Year, .keep_all = T) %>%
  left_join(ungroup(vars_single_mpa), ., by = c('seriesID', 'SamplingEventID', 'Year'))
vars_single_mpa <- vars_single_mpa %>% 
  mutate(level_fishing_protection = case_when(.default = level_fishing_protection,
    is.na(level_fishing_protection) ~ 1,
    year_established > Year ~ 1))
summary(vars_single_mpa)
# write_csv(vars_single_mpa, 'src/vars_single_mpa.csv')
