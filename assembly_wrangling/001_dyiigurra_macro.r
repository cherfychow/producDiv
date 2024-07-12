
# Formatting Dyiigurra data 

require(dplyr)
require(readr)
require(stringr)

dt <- read_csv('../../lizard_fish_data/output/lizard_uvc_clean.csv')
sites  <- st_read('/Users/cher/Downloads/dyiigurra_map/site_centroids.shp')
sites$site_name[19] <- 'Turtle Beach'

summary(dt)
dt %>% filter(is.na(size))
dt <- filter(dt, is.na(size) == F) # omit these 6 records
View(dt)

table(dt$size) # check if it's binned size classes
sort(unique(dt$species))
dt %>% filter(species == 'ni') %>% View # omit
dt  <- filter(dt, species != 'ni')
sort(unique(dt$site_id))

# merge in site coordinates
sites <- st_transform(sites, crs = 4326)
sites <- st_coordinates(sites) %>% as_tibble %>% rename(Longitude = X, Latitude = Y) %>% 
    bind_cols(st_drop_geometry(sites)) %>% 
    mutate(site_id = str_to_lower(site_name) %>% str_replace_all(., '\\s', '\\_'))
sites$site_id <- str_remove_all(sites$site_id, "\\'")
sites

colnames(dt) <- str_to_title(colnames(dt))
dt_final <- dt %>% left_join(., sites[-3], by = c("Site_id" = "site_id")) %>%
    mutate(SamplingEventID = paste(Year, Site_id, Survey_id, sep = '_'),
        Month = NA, Day = NA, Genus = word(Species, 1), 
        Species = word(Species, 2),  Depth = NA, DataID = '001_dyiigurra_macro') %>% 
    rename(Length = Size) %>% 
    group_by_all() %>% summarise(Abundance = n()) %>% ungroup %>% 
    select(Abundance, Length, Family, Genus, Species, SamplingEventID, 
           Latitude, Longitude, Depth, Day, Month, Year, DataID)
View(dt_final)

# species fixes
dt_final %>% filter(str_detect(Genus, 'idae')) %>% 
    distinct(Family, Genus, Species) %>% View
dt_final$Genus[dt_final$Genus == 'Scaridae'] <- 'Scarini'
dt_final$Genus[dt_final$Genus != 'Scarini' & str_detect(dt_final$Genus, 'idae$')] <- NA
distinct(dt_final, Family, Genus, Species) %>% View
dt_final$Species[which(is.na(dt_final$Genus))] <- NA

write_csv(dt_final, 'data_outputs/001_dyiigurra_macro_obs.csv')
