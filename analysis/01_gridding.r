
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

# load in data
tempdt <- read_csv('data_outputs/001_dyiigurra_macro_obs.csv')
tempdt2 <- read_csv('data_outputs/BT45_moorea_lter_obs.csv')
dt_assemblage <- bind_rows(tempdt, tempdt2)
rm(tempdt, tempdt2)

dt_meta <- read_csv('./src/sandbox_meta.csv') # dataset level metadata

View(dt_assemblage)

# CHECK + VALIDATE DATA ------------------------------------------------------
# global taxonomy check

# make a column to help with identifying uncertain records
dt_assemblage <- dt_assemblage %>% 
    mutate(taxon_resolution = case_when(str_detect(Species, 'sp\\d*$') & is.na(Genus) == F ~ 'genus', 
                                        is.na(Family) == F & is.na(Genus) ~ 'family',
                                        str_detect(Species, 'sp\\d*$') == F ~ 'species',
                                        str_detect(Genus, 'Scari(ni|nae)') ~ 'family'))
dt_assemblage %>% distinct(Genus, Species)

sp <- dt_assemblage %>% filter(taxon_resolution == 'species') %>% 
    distinct(Genus, Species) %>% arrange(Genus, Species) %>%
    mutate(temp = paste(Genus, Species, sep = ' ')) %>%
    pull(temp)

dt_spcheck <- validate_names(sp) # check against FishBase
dt_spcheck <- tibble(input = sp, validated = dt_spcheck)
View(dt_spcheck)
dt_spcheck %>% filter(is.na(validated)) %>% View
dt_spcheck <- dt_spcheck %>% filter(validated != input)
replace = dt_spcheck$validated
names(replace) <- dt_spcheck$input

# manual uncertain misspelling replace
dt_assemblage$Species[str_which(dt_assemblage$Genus, 'Acanthurus')] <- str_replace(dt_assemblage$Species[str_which(dt_assemblage$Genus, 'Acanthurus')], 'nigros', 'nigrofuscus')
dt_assemblage %>% filter(str_detect(Genus, 'Scarini')) %>%
    distinct(Family, Genus, Species) %>% View
dt_assemblage$Species[dt_assemblage$Genus == 'Scarinissp8'] <- 'sp8'
dt_assemblage$Genus[dt_assemblage$Genus == 'Scarinissp8'] <- 'Scarini'
dt_assemblage$Species[str_which(dt_assemblage$Species, '^\\d$')] <- paste0('sp', dt_assemblage$Species[str_which(dt_assemblage$Species, '^\\d$')])
dt_assemblage %>% distinct(Family, Genus, Species) %>% filter(is.na(Species))
dt_assemblage$Species[is.na(dt_assemblage$Genus) == F & is.na(dt_assemblage$Species)] <- 'sp'
dt_assemblage$Family[which(dt_assemblage$Family == 'Scaridae')] <- 'Labridae'
dt_assemblage$Genus[which(dt_assemblage$Genus == 'Scaridae')] <- 'Scarini'

# from FishBase validation
dt_assemblage <- dt_assemblage %>% mutate(genus_sp = case_when(
    is.na(Genus) & is.na(Species) ~ Family,
    is.na(Genus) & is.na(Species) == F ~ paste(Family, Species, sep = ' '),
    is.na(Genus) == F & is.na(Species) == F ~ paste(Genus, Species, sep = ' ')
))
dt_assemblage$genus_sp <- str_replace_all(dt_assemblage$genus_sp, replace)
dt_assemblage %>% filter(str_detect(genus_sp, '^NA')) %>% distinct %>% View
unique(dt_assemblage$genus_sp)

dt_assemblage %>% distinct(Family, Genus, Species, genus_sp) %>%
    arrange(Family, Genus, Species) %>% View
dt_assemblage %>% distinct(Genus, Species) %>% 
    filter(str_detect(Genus, 'idae'))

dt_assemblage %>% filter(is.na(Family)) %>% distinct(Family, genus_sp)
dt_assemblage <- filter(dt_assemblage, genus_sp != 'Teleostei unidentified') # omit
# take the validated names and resplit to fix the Genus and Species column
dt_assemblage <- mutate(dt_assemblage, Genus2 = case_when(str_detect(word(genus_sp, 1), 'idae$') ~ NA,
    .default = word(genus_sp, 1)), 
    # family level IDs should be NA in Species, unless it's Scarini or a known uncertain
    Species2 = case_when(is.na(word(genus_sp, 2)) == F ~ word(genus_sp, 2), .default = NA))

dt_assemblage %>% distinct(Family, Genus, Species, Genus2, Species2) %>%
    filter(Species2 != Species | Genus2 != Genus) %>%
    arrange(Family, Genus2, Species2) %>% View
dt_assemblage %>% distinct(Family, Genus, Species, Genus2, Species2) %>%
    filter(str_detect(Species, '\\d$|sp$')) %>%
    arrange(Family, Genus2, Species2) %>% View # check if the uncertainty information is preserved;

# looks ok. replace
dt_assemblage <- dt_assemblage %>% select(!c(Genus, Species)) %>% 
    rename(Genus = Genus2, Species = Species2)

rm(sp, dt_spcheck, replace)

## VALIDATE GENUS/FAMILY -------------------------------------------------------

# use taxize to interact with worms API better
gen <- dt_assemblage$Genus |> unique() |> sort() |> na.omit()
gen <- gen[-which(gen == 'Scarini')]
family <- classification(gen, db = 'worms', rows = 1) # only consider the top result

# tibblize
for(i in 1:length(family)) {
    # error handling
    try <- tryCatch(
        tidyr::pivot_wider(family[[i]][-3], names_from = rank, values_from = name) %>%
            select(contains("amily"), Genus),
        error = function(e) e
    )
    # move on if it doesn't work
    if(inherits(try, "error")) next 

    if (exists('family2')){
        family2 <- tidyr::pivot_wider(family[[i]][-3], names_from = rank, values_from = name) %>%
        select(contains("amily"), Genus) %>% bind_rows(family2, .)
    } else {
        family2 <- tidyr::pivot_wider(family[[i]][-3], names_from = rank, values_from = name) %>%
        select(contains("amily"), Genus)
    }
}
family2 <- family2 %>% distinct(Family, Subfamily, Genus) %>%
     arrange(Family, Subfamily, Genus)
family2 <- tibble(Family = 'Labridae', Genus = 'Scarini') %>% bind_rows(family2, .)
family2$Family[family2$Family == 'Scaridae'] <- 'Labridae' # outdated!
# technically wrong but Scarini needs to be in Genus for some

dt_assemblage <- family2 %>% rename(Family2 = Family) %>% select(!Subfamily) %>% left_join(dt_assemblage, ., by = "Genus")
dt_assemblage %>% filter(is.na(Family2), is.na(Genus) == F) %>% distinct(Family, Genus, Family2) %>% View
dt_assemblage %>% distinct(Family, Family2, Genus, genus_sp) %>% View
dt_assemblage$Genus[dt_assemblage$Genus == 'Perciformes'] <- NA
rm(gen, family2, family, try, i)
dt_assemblage$Family2[is.na(dt_assemblage$Genus) & is.na(dt_assemblage$Family2)] <- dt_assemblage$Family[is.na(dt_assemblage$Genus) & is.na(dt_assemblage$Family2)]
dt_assemblage <- dt_assemblage %>% select(!Family) %>% rename(Family = Family2)

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

## VISUAL CHECKS -------------------------------------------------------
# tread carefully... plotting everything will take forever

# cellcount <- dt_assemblage %>% group_by(gridcell, DataID) %>% 
#     summarise(yearN = n_distinct(Year), eventN = n_distinct(SamplingEventID))
# grid_count <- dgcellstogrid(grid, cellcount$gridcell) %>% 
#     left_join(., cellcount, by = c('seqnum' = 'gridcell'))

# world <- ne_countries(scale = 'medium', returnclass = 'sf') %>% 
# st_make_valid() %>% st_break_antimeridian(., lon_0 = 160)

# # Dyiigurra
# lon = c(min(dt_assemblage$Longitude[str_detect(dt_assemblage$DataID, 'macro')]),
#         max(dt_assemblage$Longitude[str_detect(dt_assemblage$DataID, 'macro')]))
# lat = c(min(dt_assemblage$Latitude[str_detect(dt_assemblage$DataID, 'macro')]),
#         max(dt_assemblage$Latitude[str_detect(dt_assemblage$DataID, 'macro')]))
# ggplot() +
#   geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
#   geom_sf(data = grid_count, aes(fill = yearN), shape = 21, size = 2) +
#   guides(fill = 'none') + 
#   coord_sf(xlim = lon, ylim = lat) +
#   theme_minimal()

# # Mo'orea
# lon = c(min(dt_assemblage$Longitude[str_detect(dt_assemblage$DataID, 'moorea')]),
#         max(dt_assemblage$Longitude[str_detect(dt_assemblage$DataID, 'moorea')]))
# lat = c(min(dt_assemblage$Latitude[str_detect(dt_assemblage$DataID, 'moorea')]),
#         max(dt_assemblage$Latitude[str_detect(dt_assemblage$DataID, 'moorea')]))
# ggplot() +
#   geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
#   geom_sf(data = grid_count, aes(fill = yearN), shape = 21, size = 2) +
#   guides(fill = 'none') + 
#   coord_sf(xlim = lon, ylim = lat) +
#   theme_minimal()

# rm(world, lon, lat, cellcount, grid_match, dt_coord)
# View(grid_count)
# rm(grid_count)

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