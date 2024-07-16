

# ProducDiv
# Compile the megadataset
# also gives a summary of the data

require(dplyr)
require(readr)
require(ggplot2)
require(sf)
require(rnaturalearth)
require(readxl)
require(stringr)

## =====================================================================
# COMPILE OBSERVATIONS

# pull the file names for all the datasets
PD_obs <- list.files(path = 'data_outputs', pattern = '_obs.csv', 
                     recursive = T, full.names = T)

PD_dt <- as.list(rep('', length(PD_obs))) # empty list to populate
for (i in 1:length(PD_obs)) {
  PD_dt[[i]] <- read_csv(PD_obs[i]) # read each file in
}

PD_data <- bind_rows(PD_dt)
rm(PD_dt, PD_obs)

## =====================================================================
# COMPILE METADATA

meta <- list.files(path = 'data_outputs', pattern = '_meta.xlsx', 
                     recursive = T, full.names = T)

PD_meta <- as.list(rep('', length(meta))) # empty list to populate

for (i in 1:length(PD_meta)) {
  PD_meta[[i]] <- read_excel(meta[i], col_names = F) # read each file in
  if (ncol(PD_meta[[i]]) == 3) {
    colnames(PD_meta[[i]]) <- c('fieldname', 'value1', 'value2') # temp colnames
  } else {
    colnames(PD_meta[[i]]) <- c('fieldname', 'value1') # temp colnames
  }
}

PD_cite <- as.list(rep("", length(meta)))
metadt <- as.list(rep("", length(meta)))

for (i in 1:length(PD_meta)) {
# temporary object of the bibtex in a vector
# metadata
  if (ncol(PD_meta[[i]]) == 3) {
    tempbib <- PD_meta[[i]][23, ] %>% select(value1, value2) %>% as.vector
    tempbib <- str_remove_all(tempbib, '^\\@article\\{\\}$')
    metadt[[i]] <- PD_meta[[i]][1:22,] %>% select(!value2) %>%
    tidyr::pivot_wider(., names_from = 'fieldname', values_from = 'value1')
  } else {
    tempbib <- PD_meta[[i]][23, ] %>% pull(value1)
    tempbib <- str_remove_all(tempbib, '^\\@article\\{\\}$')
    metadt[[i]] <- PD_meta[[i]][1:22,] %>%
    tidyr::pivot_wider(., names_from = 'fieldname', values_from = 'value1')
  }

# citation table
PD_cite[[i]] <- tibble(DataID = rep(unname(PD_meta[[i]][1,2]), length(tempbib)), 
                       bibtex = tempbib)
}

PD_cite <- bind_rows(PD_cite)
bind_rows(metadt) %>% head %>% View
PD_meta <- bind_rows(metadt)

View(PD_meta)
View(PD_cite)
## =====================================================================
# VISUAL CHECK

# make a spatial point object with the data
points <- PD_data %>% group_by(DataID) %>% 
  mutate(Duration = max(Year) - min(Year) + 1) %>%
  distinct(Longitude, Latitude, DataID, Year, Duration) %>%
  filter(is.na(Latitude) == F) %>%
  st_as_sf(., crs = 4326, coords = c("Longitude", "Latitude")) %>%
  st_make_valid(.) %>%
  st_cast(., to = "MULTIPOINT")


proj <- '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'
world <- ne_countries(scale = 'medium', returnclass = 'sf') %>%
st_make_valid(.) %>% st_break_antimeridian(., lon_0 = 160)
world <- st_cast(world, to = 'MULTIPOLYGON')

# by dataset ID
allPD <- ggplot() +
  geom_sf(data = world, aes(group = admin), color='grey40', fill='grey95') +
  geom_sf(data = points %>% filter(Duration > 1), 
          aes(fill = Year), shape = 21, size = 2) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal() + guides(fill = 'none') +
  scale_fill_viridis_c(option = 'C')

# ggsave(allPD, filename = "datamap_june.png", device = 'png', dpi = 300, width = 210, height = 150, units = "mm")
# need to examine this better visually without rarefying at this stage
plotly::ggplotly(allPD)
