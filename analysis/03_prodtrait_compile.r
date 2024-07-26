
# ProducDiv
# Compile traits for productivity modelling

require(dplyr)
require(readr)
require(stringr)
require(rfishprod)
require(rfishbase)


# COMPILE TRAITS ==========================================================
# load("workspace.RData")

dt_prodDB <- rfishprod::db

allsp <- dt_assemblage %>% filter(taxon_resolution %in% c('genus', 'species')) %>%
  distinct(genus_sp) %>%
  mutate(genus_sp = str_remove(genus_sp, '\\d+')) %>% 
  # collapse known uncertains because they don't matter for productivity
  pull(genus_sp) %>% na.omit %>% unique %>% sort # remove NAs
# View(allsp)

# the input data for productivity will eventually be dt_traits

# fetch traits from FishBase for species not covered
# L-W parameters, need to scrape
# MaxSizeTL
# Diet/Trophic group
# Position
# SSTMean

# Lmax + water column
fishbase_traits <- species(species_list = allsp[-str_which(allsp, '\\ssp$')], 
  fields = c('Species', 'Comments'))

# Bayesian L-W parameters
fishbase_traits <- estimate(species_list = fishbase_traits$Species, fields = c('Species', 'MaxLengthTL', 'a', 'b')) %>%
  left_join(fishbase_traits, ., by = "Species")

# get habitat information to help manually generate position trait
position <- ecology(species_list = fishbase_traits$Species) %>% select(Species, Macrophyte:Driftwood, BoatHulls:HardCorals) %>%
  tidyr::pivot_longer(., cols = Macrophyte:HardCorals, names_to = "Habitat", values_to = "YN") %>%
  filter(YN == -1)

fishbase_traits <- position %>% group_by(Species) %>% summarise(Habitat = paste(Habitat, collapse = ', ')) %>%
  full_join(fishbase_traits, ., by = "Species")

# Vertical position traits from GASPAR
gaspar <- read_csv('src/gaspar_traits.csv')
unique(dt_prodDB$Position)
replace <- c('Bnth', 'BtPl', 'Pelg')
names(replace) <- c('Bottom', 'Low', 'High')
gaspar$Position_v <- str_replace_all(gaspar$Level_water, replace)
fishbase_traits <- gaspar %>% select(FullSpecies, Position_v) %>% 
  rename(Species = FullSpecies) %>%
  left_join(fishbase_traits, ., by = "Species")

# use the two to get vertical and horizontal position traits
fishbase_traits <- mutate(fishbase_traits, Position_h = case_when(
  # associated if it's not exclusively coral reefs
    str_detect(Habitat, 'seagrass|Seagrass') ~ 'As', .default = "Dw"
  )) %>% mutate(Position = paste0(Position_v, Position_h)) %>%
  select(Species, MaxLengthTL, Position, a, b, Comments) %>%
  rename(MaxSizeTL = MaxLengthTL)

## MATCH UP TROPHIC GUILDS ==========================================================

# Morais and Bellwood don't use the same trophic guild groupings as Parravicini
# I'd rather predict Kmax from these groupings
# Match extrapolated guild predictions to rfishprod dt_prodDB
dt_troph <- read.csv("src/parravicini_trophic_guilds.csv") %>% 
  select(species, trophic_guild_predicted) %>% mutate(species = str_replace(species, '\\_', ' ')) %>%
  rename(Species = species)
replace <- c('InvSes', 'Herb', 'Cor', 'Pisc', 'InvMic', 'InvMac', 'Crus', 'Plank')
names(replace) <- as.character(1:8)
dt_troph$Diet_Par <- str_replace_all(dt_troph$trophic_guild_predicted |> as.character(), replace)
dt_troph <- dt_troph[-2]
fishbase_traits <- left_join(fishbase_traits, dt_troph, by = "Species")
rm(replace)

# Export compiled traits
# readr::write_csv(fishbase_traits, 'src/fish_prodtraits.csv')

# Make rfishprod training diet data match
dt_prodDB %>% distinct(Species, Diet) %>% left_join(., dt_troph, by = "Species") %>% View
dt_prodDB %>% distinct(Species, Diet) %>% left_join(., dt_troph, by = "Species") %>%
  pull(Diet_Par) %>% is.na %>% sum # how many species don't have a match with Parravicini guilds
dt_prodDB <- dt_prodDB %>% left_join(., dt_troph, by = "Species")
dt_prodDB %>% filter(is.na(Diet_Par)) %>% distinct(Diet)

# manual matches
dt_prodDB$Diet_Par[is.na(dt_prodDB$Diet_Par) & str_detect(dt_prodDB$Diet, '^Her')] <- 'Herb'
dt_prodDB$Diet_Par[is.na(dt_prodDB$Diet_Par) & dt_prodDB$Diet == 'FisCep'] <- 'Pisc'
dt_prodDB$Diet_Par[is.na(dt_prodDB$Diet_Par) & dt_prodDB$Diet == 'Plktiv'] <- 'Plank'
dt_prodDB$Diet_Par[is.na(dt_prodDB$Diet_Par) & dt_prodDB$Diet == 'FisCep'] <- 'Pisc'
# can't resolve invertivores like this
db_sp_troph <- dt_prodDB %>% filter(is.na(Diet_Par)) %>% pull(Species)
# worth manual matching? yes, about > 20% of training data
dt_prodDB %>% group_by(Diet) %>% summarise(n())
length(db_sp_troph)
db_sp_troph <- fooditems(species_list = db_sp_troph)
View(db_sp_troph)
rm(db_sp_troph)

# have to do it manually off-code
# readr::write_csv(db_sp_troph, 'src/missing_troph_sp.csv')
db_sp_troph <- readr::read_csv('src/missing_troph_sp.csv') # read back in
# merge in the manual records
dt_troph <- bind_rows(dt_troph, db_sp_troph)
# rejoin
dt_prodDB <- dt_prodDB %>% select(!c(Diet_Par, Diet)) %>% left_join(., dt_troph, by = "Species") %>%
  filter(is.na(Diet_Par) == F)
summary(dt_prodDB)
dt_prodDB <- rename(dt_prodDB, Diet = Diet_Par)

# restash the species that have complete compiled traits
dt_traits <- dt_prodDB %>% select(Species, MaxSizeTL, Diet, Position, a, b) %>%
  filter(Species %in% allsp) %>% arrange(Species, desc(MaxSizeTL), Position) %>%
  ungroup %>% group_by(Species) %>% mutate(order = 1:n()) %>%
  filter(order == 1) %>% select(!order) %>% ungroup

rm(gaspar, db_sp_troph)

## REMERGE ====================================================================

fishbase_traits <- readr::read_csv('src/fish_prodtraits.csv')

# check if there's any overlap with dt_traits from the Morais database
fishbase_traits %>% filter(Species %in% dt_traits$Species) %>% nrow
dt_traits %>% filter(Species %in% fishbase_traits$Species) %>% nrow
dt_traits <- dt_traits %>% filter(Species %in% fishbase_traits$Species == F) %>% # ensure no duplicates
  bind_rows(., fishbase_traits) %>%
  arrange(Species)


## GENUS LEVEL TRAITS =================================================
# are they worth keeping
# dt_rareass %>% filter(iter == 1, str_detect(genus_sp, '\\ssp\\d*$')) %>%
#  group_by(genus_sp) %>% summarise(sumN = sum(Abundance)) %>% View
# yes

# update species gap
genus <- allsp[allsp %in% dt_traits$Species == F] |> word(1)
# take the genus mean for these traits?
traits_genus <- dt_traits %>% mutate(Genus = word(Species, 1)) %>% 
  filter(Genus %in% genus) %>% select(Genus, MaxSizeTL, a, b) %>%
  group_by(Genus) %>% 
  summarise(MaxSizeTL = mean(MaxSizeTL), a = mean(a), b = mean(b))

# frequency based matching for categorical
diet <- dt_traits %>% mutate(Genus = word(Species, 1)) %>% 
  filter(Genus %in% genus) %>% ungroup %>% select(Genus, Diet) %>%
  group_by_all %>% summarise(freq = n()) %>% arrange(Genus, desc(freq)) %>%
  mutate(order = 1:n())

position <- dt_traits %>% mutate(Genus = word(Species, 1)) %>% 
  filter(Genus %in% genus) %>% ungroup %>% select(Genus, Position) %>%
  group_by_all %>% summarise(freq = n()) %>% arrange(Genus, desc(freq)) %>%
  mutate(order = 1:n())

diet <- diet %>% filter(order == 1) %>% ungroup %>% select(!c(freq, order))
position <- position %>% filter(order == 1) %>% ungroup %>% select(!c(freq, order))
traits_genus <- left_join(traits_genus, diet, by = 'Genus') %>%
  left_join(., position, by = "Genus") %>%
  mutate(Species = paste0(Genus, ' sp'))
rm(diet, position, genus)

dt_traits <- dt_traits %>% mutate(Genus = word(Species, 1)) %>% 
  bind_rows(., traits_genus) %>% arrange(Species)
View(dt_traits)

allsp[-which(allsp %in% dt_traits$Species)] # manual uncertain compiles
# allsp[-which(allsp %in% dt_traits$Species)]  %>% clipr::write_clip(.)
rm(traits_genus)
dt_traits <- read_csv('./src/missing_traits_genus.csv') %>% bind_rows(dt_traits, .)
rm(fishbase_traits)
