
# ProducDiv
# Compile historic OISST data that matches fish data

# switch to data table for speed
require(tidytable)
require(sf)
require(lubridate)
require(tidync) # for netCDF
require(rnaturalearth)
require(dggridR)
require(ggplot2)


# https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.highres.html

# make sure there is monthly resolution for data
summary(dt_rareass)
summary(dt_single)

# make time resolution column
dt_rareass <- dt_rareass %>% mutate(date_res = case_when(is.na(Month) == F & is.na(Year) == F ~ 'monthly', 
                                                         is.na(Month) & is.na(Year) == F ~ 'yearly',
                                                         .default = NA))
dt_single <- dt_single %>% mutate(date_res = case_when(is.na(Month) == F & is.na(Year) == F ~ 'monthly', 
                                                         is.na(Month) & is.na(Year) == F ~ 'yearly',
                                                         .default = NA))
# tallies
dt_single %>% group_by(date_res) %>% summarise(n = n_distinct(DataID))
dt_rareass %>% group_by(date_res) %>% summarise(n = n_distinct(DataID))

# give data records IDs to keep track
dt_single$ID <- paste0('s', 1:nrow(dt_single))
dt_rareass$ID <- paste0('r', 1:nrow(dt_rareass))

# EXTRACT OISST DATA ================================================================

OISST <- tidync::tidync('./src/noaa_oisst_monthly_2025.nc')
OISST <- OISST %>% 
  hyper_filter(lat = between(lat, -30, 30)) %>% 
  tidync::hyper_tibble() %>% # output as tibble
  dplyr::select(lon, lat, time, sst) %>% 
  arrange(lat, lon, time)
OISST$date <- as.Date(OISST$time, origin = "1800/01/01")
OISST <- OISST %>% mutate(Year = year(date), Month = month(date)) %>% 
  select(!c(date, time))

# MATCH OISST TO RAREFY GRID ========================================================

# the grids are 0.25º and only cover sea area
# make an sf grid that matches the SST data

# make an sf point object of OISST coverage
oisst_points <- OISST %>% distinct(lat, lon) %>% 
  mutate(sst_lat = lat, sst_lon = lon) %>% # make a column that keeps the 0-360° coordinates
  mutate(lon = case_when(.default = lon, 
                         lon <=180 ~ lon, # east
                         lon > 180 ~ lon-360)) %>% # and then convert these to -180-180 for matching
  st_as_sf(., coords = c('lon', 'lat'), crs = 4326) # converted ones to sf

# set the boundary for the grid
boundary = st_polygon(list(rbind(c(-180,30), c(180,30), c(180,-30), c(-180,-30), c(-180,30))), dim = 'XY')

# convert the points to a grid by covering the whole thing and then subtracting out land grids
grid_sst <- st_make_grid(boundary, cellsize = 0.25, crs = 4326, what = 'polygons', square = T)
noland <- st_intersects(grid_sst, oisst_points) # which grid cells match the sst points
grid_sst <- grid_sst[which(lengths(noland) == 1),] # use noland to index the grid and remove land areas

# sanity check
noland[[1]]
grid_sst[1]
oisst_points[noland[[1]],]

# subset oisst points to match up with grid_sst
grid_sst <- oisst_points[noland[which(lengths(noland) == 1)] %>% unlist,] %>% 
   bind_cols(grid_sst, .) # match in
colnames(grid_sst)[1] <- 'geometry'
colnames(grid_sst)[4] <- 'oisst_cent'
grid_sst <- st_as_sf(grid_sst, crs = 4326) # coerce back to sf dataframe

allcells <- unique(dt_coord$gridcell)
# get our grid cell centroids and then match that to the closest centroid of the OISST grids
grid_sf <- dgcellstogrid(grid, cells = allcells, return_sf = T)
grid_cent <- grid_sf %>% mutate(geometry = st_centroid(geometry))

# use intersect to match up rarefy grids with sst grids
# some coastal coordinates might have trouble intersecting with OISST, so we'll use nearest neighbour
raresst_match <- st_nearest_feature(grid_sf, grid_sst) # possible to intersect with multiple
# did every grid get a match
(which(lengths(raresst_match) > 0) |> length()) == nrow(grid_sf)

# SPATIAL MATCH!
dt_sst <- grid_sst[raresst_match[which(lengths(raresst_match) > 0)] |> unlist(),] %>% 
  bind_cols(st_drop_geometry(grid_sf),.)
dt_sst <- dt_sst %>% rename(gridcell = seqnum)
rm(oisst_points, boundary, noland, raresst_match, grid_sst)
gc()

# # time ranges for each gridcell
# 
# dt_sst_yr <- bind_rows(dt_rareass, dt_single) %>% filter(date_res == 'yearly') %>% 
#   distinct(Year, gridcell) %>% 
#   group_by(gridcell) %>%
#   summarise(year_min = min(Year), year_max = max(Year)) %>% distinct %>% # calculate the year ranges needed for that gridcell
#   left_join(., dt_sst, by = 'gridcell') %>% 
#   mutate(date_min = paste(year_min, 01, 01, sep = '-') %>% as.Date(., tryFormats= '%Y-%m-%d'), 
#          date_max = paste(year_max, 12, 31, sep = '-') %>% as.Date(., tryFormats= '%Y-%m-%d'))
# 
# dt_sst_mo <- bind_rows(dt_rareass, dt_single) %>% filter(date_res == 'monthly') %>% 
#   distinct(Year, Month, gridcell) %>% 
#   mutate(date = paste(Year, Month, 15, sep = '-') %>% as.Date(., tryFormats= '%Y-%m-%d')) %>% 
#   group_by(gridcell) %>%
#   summarise(date_min = min(date), date_max = max(date)) %>% 
#   left_join(., dt_sst, by = 'gridcell')
# 
# gc() # help free up memory
# 
# # check for any NAs
# dt_sst_yr %>% filter(is.na(date_min) | is.na(date_max))
# dt_sst_mo %>% filter(is.na(date_min) | is.na(date_max))
# 
# # spatial filter
# # filter for grids within our lat long range
# # this is not computationally efficient for the full dataset without the server
# 
# dt_sst <- bind_rows(dt_sst_mo, dt_sst_yr) %>% select(!c(year_min, year_max)) %>% 
#   group_by(gridcell, sst_latgrid, sst_longrid) %>% 
#   summarise(date_min = min(date_min), date_max = max(date_max))
# summary(dt_sst)
# 
# # convert the tidync object to a tibble
# # coarse spatial filter
# 
# 
# gc() # free up memory after dealing with the giant data layer

# SPACE * TIME OISST MATCH ================================================================

# monthly

# make the data join base first
# grab all the required months and years from each gridcell
dt_sst_mo <- bind_rows(dt_rareass, dt_single) %>% filter(date_res == 'monthly') %>% 
  distinct(gridcell, Year, Month)
# make sure there won't be many-to-many matches

# match in SST grid coordinates
dt_sst_mo <- dt_sst %>% st_drop_geometry() %>% select(gridcell, sst_lat, sst_lon) %>% 
  distinct %>% left_join(dt_sst_mo, ., by = 'gridcell')
summary(dt_sst_mo)

# matches up each gridcell to the SST grid centroid
# essentially magnifies dt_coord by time
# use data table for a faster join
dt_sst_mo <- left_join(dt_sst_mo, OISST, by = c('sst_lat' = 'lat', 'sst_lon' = 'lon', 'Year', 'Month'))
summary(dt_sst_mo)

# yearly

# make the data join base first
# grab all the required months and years from each gridcell
dt_sst_yr <- bind_rows(dt_rareass, dt_single) %>% filter(date_res == 'yearly') %>%
  distinct(gridcell, Year)
# make sure there won't be many-to-many matches

# match in SST grid coordinates
dt_sst_yr <- dt_sst %>% st_drop_geometry() %>% select(gridcell, sst_lat, sst_lon) %>% 
  distinct %>% left_join(dt_sst_yr, ., by = 'gridcell')
summary(dt_sst_yr)

# aggregate OISST to yearly mean and match with SST
dt_sst_yr <- OISST %>% group_by(lat, lon, Year) %>% summarise(sst = mean(sst)) %>% 
  ungroup %>% 
  # then join in
  left_join(dt_sst_yr, OISST, by = c('sst_lat' = 'lat', 'sst_lon' = 'lon', 'Year'))

dt_sst_mo$date_res <- 'monthly'
dt_sst_yr$date_res <- 'yearly'

dt_sst <- bind_rows(dt_sst_mo, dt_sst_yr) 

rm(OISST)
gc()