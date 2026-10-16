
# load working data

require(readr)
require(tidytable)
require(sf)

# raw data
dt_single <- read_csv('src/dt_single.csv')
dt_rareass <- read_csv('src/dt_rareass.csv')
dt_assemblage <- bind_rows(dt_single, dt_rareass)
dt_coord <- dt_assemblage %>% distinct(Latitude, Longitude, gridcell)

# data for fitting models
vars_rare <- read_csv('src/vars_rare.csv')
vars_time <- read_csv('src/vars_time.csv')
vars_single <- read_csv('src/vars_single.csv')
vars_single_mpa <- read_csv('src/vars_single_mpa.csv')
vars_time_flux <- read_csv('src/vars_time_flux.csv')

ecoR <- st_read('src/biogeo.shp')

# models
load('outputs/m3_timeN.rdata')
load('outputs/m4_timeflux.rdata')
load('outputs/m1_spatial.rda')
load('outputs/m2_spatialN.rda')
