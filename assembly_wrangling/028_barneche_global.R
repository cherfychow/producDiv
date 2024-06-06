

# Curation Script ---------------------------------------------------------

# Dataset: Multiple datasets from Body size, reef area and temperature predict global reef-fish species richness across spatial scales
# Location: Global
# Curator: Cher Chow

# Set up ------------------------------------------------------------------
# load the necessary packages
require(tidyverse)
require(sf)
require(world)

rm(list=ls())

# make sure your working directory is set before running these lines
dt <- read_csv('data_originals/028_barneche_2018/data.csv')
View(dt)
dt <- dt %>% filter(studyName != 'rls_march_2017')

# mostly scuba methods, but one survey is with snorkel --> split
summary(dt)
str(dt)
dt %>% filter(is.na(size_cm)) %>% View
dt %>% group_by(type) %>% summarise(total = n_distinct(studyName))

# From Diego's repo readme
# studyName: Original subset of data shared via co-authors (internal use only);
# region: Province name;
# locality: Sub-province name;
# sites: Site name;
# area: Area of individual transects, in m2;
# type: Type of transects:
  # quant_class = abundance and size (cm) tallied; 
  # pres_abs = presence / absence; 
  # quant_total = total abundance, but no size tallied;
# transect_id: Individual transect ID;
# lon: Site longitude;
# lat: Site latitude;
# depth_m: Depth in meters;
# species: Species name;
# size_cm: Size of species tallied (cm);
# abun: Total number of individuals;
# size_class: Species size class (0--7 cm, 7--15 cm, 15--30 cm, 30--50 cm, 50--80 cm and > 80 cm) based on maximum size;
# max_size: Species maximum size published in the literature (cm)

dt %>% filter(type == 'pres_abs') %>% View
dt %>% filter(type == 'quant_total') %>% View
dt %>% filter(type == 'pres_abs') %>% group_by(sites) %>% summarise(transects = n_distinct(transect_id))

dt <- dt %>% filter(!type == 'pres_abs') # omit 142 transects from Bahamas, Cuba, Belize

# Structure check ---------------------------------------------------------

nrow(dt) # check dimensions
# 536158 records
# Abundance and/or biomass, latitude and longitude numeric? Y

# Dates
dt <- dt %>% mutate(Day = NA, Month = NA, Year = NA)

# Secondary fields relating to sample must be factors or integers?
# fix depth
dt %>% filter(str_detect(depth_m, '\\_')) %>% pull(depth_m) %>% unique
unique(dt$depth_m) # mostly numbers but will keep the qualitative classifications to remove later
summary(dt$area) # nested within transect?
dt %>% group_by(region, locality, sites) %>% summarise(n_distinct(area)) %>% View # yes

# Taxonomic fields must be characters or factors? Y

# Primary field check -----------------------------------------------------

# ABUNDANCES
# No negative values, zeroes, or NAs in abundance/biomass fields.
sum(is.na(dt$abun)+0)
min(dt$abun, na.rm=T) # no zeroes? Y
sum(is.na(dt$size_cm)+0) # 44604 NAs for length but might be for size classes
min(dt$size_cm, na.rm=T) # no zeroes? Y

sum(is.na(dt$size_class)+0) # no NAs in size_classes

# YEAR MONTH DAY
# non-temporal dataset. ignore.

# Latitude constrained from -90 to 90.
# Longitude constrained -180 to 180.
summary(dt[c('lat', 'lon')]) # good

proj = '+proj=eqearth +lon_0=180 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
points <- st_as_sf(dt %>% distinct(lon, lat, studyName), coords = c("lon", "lat"), crs = 4326)
world <- sf::st_as_sf(maps::map('world', plot = FALSE, fill = TRUE))
ggplot() +
  geom_sf(data = world, color='grey40', fill='grey95') +
  geom_sf(data = points, aes(fill = studyName), shape = 21, shape = 2) +
  guides(fill = 'none') + coord_sf(crs = proj)

# the locations look like they make sense
rm(points, world, proj)

# Taxonomic field check ---------------------------------------------------

# No NAs in taxonomic fields, remove all non-organism records
sum(is.na(dt$species) + 0) == 0
sum(str_detect(dt$species, '^\\s+$') + 0) == 0

dt %>% distinct(species) %>% arrange(species) %>% View
dt %>% distinct(species) %>% filter(str_count(species, '\\_') > 1) %>% View
require(worrms)
n_distinct(dt$species)

seq <- seq(1, n_distinct(dt$species), by=150)
taxonclean <- as.list(rep(NA, length(seq)-1))
for (i in 1:(length(seq)-1)) {
  taxonclean[[i]] <- wm_records_names(unique(dt$species)[seq[i]:(seq[i+1]-1)], fuzzy = F, marine_only = T)
}
taxonclean[[length(seq)]] <- wm_records_names(unique(dt$species)[seq[[length(seq)]]:n_distinct(dt$species)], fuzzy = F, marine_only = T)
for (i in 1:length(taxonclean)) {
  taxonclean[[i]] <- bind_rows(taxonclean[[i]])
}
taxonclean <- bind_rows(taxonclean)
replace <- taxonclean$valid_name[which(taxonclean$status == 'unaccepted')] # replace unaccepted names
names(replace) <- taxonclean$scientificname[which(taxonclean$status == 'unaccepted')]
dt$species <- str_to_sentence(dt$species)
# replace
dt$taxa <- str_replace_all(dt$species, '\\_', ' ') %>% str_replace_all(., replace)

dt %>% distinct(taxa) %>% arrange(taxa) %>% View
dt %>% distinct(taxa) %>% arrange(taxa) %>% filter(str_detect(taxa, '\\d')|str_count(taxa, '\\s')>1) %>% View
dt %>% distinct(taxa) %>% arrange(taxa) %>% filter(str_detect(taxa, '\\[')) %>% View
# most of these are place designations to uncertain IDs but this would be nested in dataset already
replace <- c('\\s(ao)*tep' = '', # nested in a dataset
             'Apogon zoramia' = 'Zoramia',
             '\\[2\\]' = '2',
             '\\s\\[[a-zA-Z\\s]+\\]' = '',
             '\\sn%' = '',
             '\\sn2\\sspn\\sstome' = '',
             '\\s(?=\\d$)' = '',
             '\\sn$' = '1')
dt$taxa <- str_replace_all(dt$taxa, replace)

# check that genera are genera, not family names (-idae/eae)
dt$Genus <- word(dt$taxa, 1)
# match in family
dt <- taxonclean %>% distinct(family, genus) %>% rename(Family = family) %>% 
  mutate(Family = str_replace(Family, 'Scaridae', 'Labridae')) %>% 
  left_join(dt, ., by=c("Genus" = "genus"))
dt$Genus %>% str_detect(., 'idae$|ini$') %>% sum(. + 0)
dt %>% filter(str_detect(Genus, 'idae$|ini$')) %>% distinct(Genus) %>% View

dt$Family[str_which(dt$Genus, 'idae$')] <- dt$Genus[str_which(dt$Genus, 'idae$')]
dt$Genus[str_which(dt$Genus, 'idae$')] <- ''

dt$Species <- word(dt$taxa, start=2, end=-1) # keep only specific epithet
dt$species <- NULL
dt$taxa <- NULL

# Prepare raw data --------------------------------------------------------
dt %>% group_by(studyName) %>% summarise(areas = n_distinct(area)) %>% View
# most areas are nested within study methods, but not agreen raja ampat
dt$abun[which(dt$studyName == 'agreen_raja_ampat_2012')] <- dt$abun[which(dt$studyName == 'agreen_raja_ampat_2012')] / (dt$area[which(dt$studyName == 'agreen_raja_ampat_2012')] / 500) # standardise abundances to individuals per 500 m2 (smallest sampling unit)

# which studies have size-classes only
dt %>% filter(is.na(size_cm), is.na(size_class) == F) %>% View
dt <- dt %>% filter(is.na(size_cm) == F) %>% select(!size_class)
sort(unique(dt$depth_m))
dt$depth_m[str_which(dt$depth_m, '^[a-zA-Z]+$')] <- NA
dt$depth_m[str_which(dt$depth_m, '0\\_15')] <- NA

dt_final <- dt %>% mutate(SamplingEventID = paste(region, locality, sites, depth_m, transect_id, sep = "_"),
                          DataID = paste('028_barneche_', studyName)) %>% 
  select(studyName, abun, size_cm, Family, Genus, Species, 
         SamplingEventID, lat, lon, depth_m, Day, Month, Year) %>% 
  rename(Abundance = abun, Length = size_cm, Longitude = lon, Latitude = lat, Depth = depth_m)

# Export final ------------------------------------------------------------

for (i in 1:n_distinct(dt_final$studyName)) {
  study = unique(dt_final$studyName)[i]
  write_csv(dt_final %>% filter(studyName == study) %>% select(!studyName), 
            paste0('data_outputs/028_barneche_global/028_', study, '_obs.csv'))
}

# Atlantic TEP needs more splitting

dt_tep <- dt %>% filter(studyName == 'atlantic_tep') %>% 
  mutate(SamplingEventID = paste(region, locality, sites, depth_m, transect_id, sep = "_"),
                          DataID = paste('028_barneche_', studyName)) %>% 
  select(region, locality, abun, size_cm, Family, Genus, Species, 
         SamplingEventID, lat, lon, depth_m, Day, Month, Year) %>% 
  rename(Abundance = abun, Length = size_cm, Longitude = lon, Latitude = lat, Depth = depth_m)

dt_tep$region <- paste(dt_tep$region, dt_tep$locality, sep = "_")

for (i in 1:n_distinct(dt_tep$region)) {
  study = unique(dt_tep$region)[i]
  write_csv(dt_tep %>% filter(region == study) %>% select(!c(region, locality)), 
            paste0('data_outputs/028_barneche_global/028_', study, '_obs.csv'))
}

# Size class uncertainty ------------------------------------------------

dt_final <- read_csv('data_outputs/size_class/028_cocos_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths <- lengths[lengths$Length != 22, ]
lengths$Length_lwr = c(seq(2.5, 87.5, by = 5), 95, seq(115, 175, by = 10), 205, 235, 245, 295)
lengths$Length_upr = c(seq(7.5, 92.5, by = 5), 105, seq(125, 185, by = 10), 215, 245, 255, 305)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_cocos_obs.csv')

dt_final <- read_csv('data_outputs/size_class/028_offshore_tep_malpelo_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths <- lengths[lengths$Length != 18, ]
lengths$Length_lwr = c(seq(2.5, 67.5, by = 5), 77.5, 95, 115, 145)[-4]
lengths$Length_upr = c(lengths$Length_lwr[1:14] + 5, lengths$Length_lwr[15:17] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_offshore_tep_malpelo_obs.csv')

dt_final <- read_csv('data_outputs/size_class/028_sw_atlantic_arraial_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length_lwr = c(seq(2.5, 47.5, by = 5), 57.5)[-c(2,4)]
lengths$Length_upr = c(lengths$Length_lwr + 5)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_sw_atlantic_arraial_obs.csv')

dt_final <- read_csv('data_outputs/size_class/028_sw_atlantic_ilha_gde_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length_lwr = c(2.5, 12.5, 22.5, 32.5, 47.5)
lengths$Length_upr = c(lengths$Length_lwr + 5)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_sw_atlantic_ilha_gde_obs.csv')

dt_final <- read_csv('data_outputs/size_class/028_sw_atlantic_l_santos_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length_lwr = c(2.5, 12.5, 22.5, 32.5, 47.5)
lengths$Length_upr = c(lengths$Length_lwr + 5)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_sw_atlantic_l_santos_obs.csv')

dt_final <- read_csv('data_outputs/size_class/028_sw_atlantic_rio_de_janeiro_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length_lwr = c(2.5, 12.5, 22.5, 32.5)
lengths$Length_upr = c(lengths$Length_lwr + 5)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_sw_atlantic_rio_de_janeiro_obs.csv')

dt_final <- read_csv('data_outputs/size_class/SF/028_eastern_atlantic_cverde_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length_lwr = c(2.5, 12.5, 22.5, 27.5, 32.5, 47.5)
lengths$Length_upr = c(lengths$Length_lwr + 5)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_eastern_atlantic_cverde_obs.csv')

dt_final <- read_csv('data_outputs/size_class/SF/028_eastern_atlantic_stome_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length_lwr = c(2.5, 12.5, 22.5, 27.5, 32.5, 85, 95)
lengths$Length_upr = c(lengths$Length_lwr[1:5] + 5, lengths$Length_lwr[6:7] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_eastern_atlantic_stome_obs.csv')

dt_final <- read_csv('data_outputs/size_class/SF/028_offshore_sw_atlantic_noronha_obs.csv')
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length
lengths$Length_lwr = c(2.5, 12.5, 22.5, 27.5, 32.5, 37.5, 47.5, 65, 95, )length
lengths$Length_upr = c(lengths$Length_lwr[1:5] + 5, lengths$Length_lwr[6:7] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, 'data_outputs/028_barneche_global/028_eastern_atlantic_stome_obs.csv')

fileid = '028_offshore_sw_atlantic_rocas'
dt_final <- read_csv(paste0('data_outputs/size_class/SF/', fileid, '_obs.csv'))
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length
lengths$Length_lwr = c(2.5, 12.5, seq(22.5, 62.5, by = 5), 77.5, 95, 115, 145)
lengths$Length_upr = c(lengths$Length_lwr[1:12] + 5, lengths$Length_lwr[13:15] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, paste0('data_outputs/028_barneche_global/', fileid, '_obs.csv'))

fileid = '028_offshore_sw_atlantic_trindade'
dt_final <- read_csv(paste0('data_outputs/size_class/SF/', fileid, '_obs.csv'))
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length
lengths$Length_lwr = c(2.5, 12.5, seq(22.5, 67.5, by = 5), 77.5, 82.5, 87.5, 95, 115, 145)
lengths$Length_upr = c(lengths$Length_lwr[1:15] + 5, lengths$Length_lwr[16:18] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, paste0('data_outputs/028_barneche_global/', fileid, '_obs.csv'))

fileid = '028_sw_atlantic_abrolhos'
dt_final <- read_csv(paste0('data_outputs/size_class/SF/', fileid, '_obs.csv'))
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length
lengths$Length_lwr = c(2.5, 12.5, seq(22.5, 37.5, by = 5), seq(45, 65, by = 10))
lengths$Length_upr = c(lengths$Length_lwr[1:6] + 5, lengths$Length_lwr[7:9] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, paste0('data_outputs/028_barneche_global/', fileid, '_obs.csv'))

fileid = '028_sw_atlantic_guarapari'
dt_final <- read_csv(paste0('data_outputs/size_class/SF/', fileid, '_obs.csv'))
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length
lengths$Length_lwr = c(2.5, 12.5, 22.5, 32.5, 37.5, 45)
lengths$Length_upr = c(lengths$Length_lwr[1:5] + 5, lengths$Length_lwr[6] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, paste0('data_outputs/028_barneche_global/', fileid, '_obs.csv'))

fileid = '028_sw_atlantic_salvador_bts'
dt_final <- read_csv(paste0('data_outputs/size_class/SF/', fileid, '_obs.csv'))
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length
lengths$Length_lwr = c(2.5, 12.5, seq(22.5, 47.5, by = 5), 55, 75, 145)
lengths$Length_upr = c(lengths$Length_lwr[1:8] + 5, lengths$Length_lwr[9:nrow(lengths)] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, paste0('data_outputs/028_barneche_global/', fileid, '_obs.csv'))

fileid = '028_sw_atlantic_santa_catarina'
dt_final <- read_csv(paste0('data_outputs/size_class/SF/', fileid, '_obs.csv'))
lengths <- dt_final %>% distinct(Length) %>% arrange(Length)
lengths$Length
lengths <- lengths[lengths$Length != 11,]
lengths$Length_lwr = c(2.5, 12.5, seq(22.5, 57.5, by = 5), seq(65,95, by = 10), 115)
lengths$Length_upr = c(lengths$Length_lwr[1:10] + 5, lengths$Length_lwr[11:nrow(lengths)] + 10)
dt_final <- left_join(dt_final, lengths, by="Length") %>% relocate(starts_with('Length'), .after = Abundance)
write_csv(dt_final, paste0('data_outputs/028_barneche_global/', fileid, '_obs.csv'))

# Spatial calculations ------------------------------------------------

# load libraries
library(sf)
library(clipr)

# 1. Convert data points into point spatial object
# assumes your rawdata dataframe is called rawdata
dt_coord <- dt_sc_p %>% select(Longitude, Latitude) %>% distinct() %>%
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326) %>% st_union()

# 2. Calculate convex hull, area and centroid
centroid <- dt_coord %>% st_convex_hull() %>% st_centroid() %>% unlist # get centroid
write_clip(centroid[c(2,1)]) # copy as lat-long