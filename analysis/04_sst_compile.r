
# ProducDiv
# Compile historic OISST data that matches fish data

require(dplyr)
require(sf)
require(lubridate)
require(tidync) # for netCDF
require(rnaturalearth)
require(ggplot2)

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
# only coarse filters or else something goes wrong
OISST <- OISST %>% 
  hyper_filter(lon = lon >= min(dt_sst$sst_longrid) & lon <= max(dt_sst$sst_longrid), # targeted spatial filter
               lat = lat >= min(dt_sst$sst_latgrid) & lat <= max(dt_sst$sst_latgrid)) %>% 
  tidync::hyper_tibble() %>% # output as tibble
  dplyr::select(lon, lat, time, sst) %>% 
  arrange(lat, lon, time)
OISST$date <- as.Date(OISST$time, origin = "1800/01/01")
# double check that it starts with 1981-09-01!

# fine filters
OISST <- OISST %>% 
  filter(lon %in% dt_sst$sst_longrid, lat %in% dt_sst$sst_latgrid,
         date >= min(dt_sst$date_min), date <= max(dt_sst$date_max)) %>% 
  # some rearranging
  arrange(lat, lon, date) %>% select(!time) %>%
  # extract year and month to help match with our data
  mutate(Year = year(date), Month = month(date))

# because we used this per coordinate spatial filter
# we shouldn't need to use the land mask data

# now set up dt_sst to contain our dataset site information * SST layer

dt_sst <- dt_sst %>% select(!c(date_min, date_max))

# MERGE SITE WITH SST ================================================================

dt_sst <- dt_assemblage %>% distinct(DataID, gridcell, Year, Month) %>%
  # make sure there won't be many-to-many matches
  left_join(., dt_sst, by = 'gridcell') %>%
  # matches up each gridcell to the SST grid centroid
  # essentially magnifies dt_coord by time
  left_join(., OISST[-4], by = c('sst_latgrid' = 'lat', 'sst_longrid' = 'lon',
                                 'Year', 'Month')) %>%
  arrange(gridcell, Year, Month)

# visualise to check

ggplot() +
  geom_tile(data = dt_sst, aes(x = sst_longrid, y = sst_latgrid, fill = sst)) +
  theme_minimal() + guides(fill = 'none') +
  scale_fill_viridis_c(option = 'C')

rm(OISST, lats, lons, convhull)
