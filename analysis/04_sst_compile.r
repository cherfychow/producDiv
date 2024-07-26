
# ProducDiv
# Compile historic OISST data that matches fish data

require(dplyr)
require(sf)
require(lubridate)
require(tidync) # for netCDF

# https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html

summary(dt_assemblage)
# monthly resolution best for data

# EXTRACT OISST DATA ================================================================

OISST <- tidync::tidync('./src/noaa_oisst_v2_monthly.nc')
# filter by convex hull and time

convhull <- dt_coord %>% st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326) %>% 
  ungroup %>% st_union %>% st_convex_hull %>% st_coordinates %>% .[,1:2] %>%
  as.data.frame

# Generate set of grid cell latitudes (center of cell) from south to north
lats = seq(-89.875,89.875,0.25)
# Generate set of grid cell longitudes (center of cell)
lons = seq(0.125,359.875,0.25)

# generate coordinates that match the closest OISST grid centre
dt_coord <- dt_coord %>%
  mutate(sst_latgrid = lats[max.col(-abs(outer(Latitude,lats,"-")))],
         sst_lon = case_when(Longitude < 0 ~ Longitude + 360, .default = Longitude),
         sst_longrid = lons[max.col(-abs(outer(sst_lon,lons,"-")))])
dt_coord$sst_lon <- NULL # remove it since it was just an intermediate
# dt_coord will link which gridcell links to which OISST grid

# use dataset gridcells as the linking ID
dt_sst <- dt_coord %>% distinct(gridcell, sst_latgrid, sst_longrid)

# time ranges for each gridcell
dt_sst <- dt_assemblage %>% group_by(gridcell) %>% 
  distinct(Year, Month, gridcell) %>%
  mutate(date = paste(Year, Month, '01', sep = '-') |> as.Date()) %>%
  summarise(date_min = min(date), date_max = max(date)) %>% 
  left_join(dt_sst, ., by = "gridcell")

# spatial filter
# this is not computationally efficient for the full dataset without the server
OISST <- OISST %>% 
  hyper_filter(lon = lon %in% dt_sst$sst_longrid,
               lat = lat %in% dt_sst$sst_latgrid, # targeted spatial filter
               time >= min(dt_sst$date_min), # coarse time filter
               time <= max(dt_sst$date_max)) %>% 
  tidync::hyper_tibble() %>% 
  dplyr::select(lon, lat, time, sst) %>% 
  dplyr::mutate(date = as.Date(time, origin = "1800-1-1"))

# because we used this per coordinate spatial filter
# we shouldn't need to use the land mask data

# time filter
dt_sst <- dt_sst %>% select(!gridcell) %>% distinct %>% 
  group_by(sst_latgrid, sst_longrid) %>% 
  summarise(date_min = min(date_min), date_max = max(date_max)) %>% 
# get the ranges for the courser sst grid
  left_join(OISST, ., by = c("lat" = 'sst_latgrid', 'lon' = 'sst_longrid')) %>% 
  select(!time) %>% group_by(lat, lon) %>%
  filter(date >= date_min, date <= date_max)

dt_sst <- dt_sst %>% select(!c(date_min, date_max)) %>% 
  mutate(Year = year(date), Month = month(date))

# MERGE SITE WITH SST ================================================================

dt_coord

dt_assemblage %>% distinct(DataID, gridcell, Year, Month) %>%
  left_join(., dt_coord, by = c('DataID', 'gridcell')) %>% 
  rename(lat = sst_latgrid, lon = sst_longrid) %>% head
  left_join(., dt_sst[-4], by = c(''))
