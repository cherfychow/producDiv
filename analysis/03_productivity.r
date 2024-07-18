
# ProducDiv
# Calculate productivity

require(dplyr)
require(stringr)
require(sf)
require(rfishprod)
require(rfishbase)


# Compile traits ==========================================================
# load("/Users/cherchow/Dropbox/ch4_producDiv/datarepo/workspace.RData")

dt_prodDB <- rfishprod::db

sp_gap <- dt_assemblage %>% distinct(Family, Genus, Species) %>%
  mutate(genus_sp = case_when(is.na(Genus) == F ~ paste(Genus, Species, sep = ' ') %>%
	str_remove(., '\\d+'), .default = NA)) %>% 
  # collapse known uncertains because they don't matter for productivity
  arrange(Family, genus_sp) %>% distinct %>%
  filter(!is.na(genus_sp)) # remove NAs
View(sp_gap)

# the input data for productivity will eventually be dt_traits
# put everything in right now that that already exists
dt_traits <- dt_prodDB %>% select(Species, MaxSizeTL, Diet, Position, a, b) %>%
  filter(Species %in% sp_gap$genus_sp) %>% distinct


# fetch traits from FishBase for species not covered
# L-W parameters, need to scrape
# MaxSizeTL
# Diet/Trophic group
# Position
# SSTMean

# Lmax + water column
fishbase_traits <- sp_gap %>% 
  filter(genus_sp %in% dt_traits == F, str_detect(genus_sp, '\\ssp$') == F) %>%
  pull(genus_sp) %>% species(species_list = ., 
  fields = c('Species', 'Length', 'LTypeMaxM', 'LengthFemale', 'LTypeMaxF', 'DemersPelag'))

# Bayesian L-W parameters
fishbase_traits <- sp_gap %>% 
  filter(genus_sp %in% dt_traits == F, str_detect(genus_sp, '\\ssp$') == F) %>%
  pull(genus_sp) %>% 
  estimate(species_list = ., fields = c('Species', 'MaxLengthTL', 'a', 'b')) %>%
  left_join(fishbase_traits, ., by = "Species")
sum(is.na(fishbase_traits$MaxLengthTL))
sum(is.na(fishbase_traits$Length))
fishbase_traits %>% filter(LTypeMaxM != 'TL') %>% View

# get habitat information to help manually generate position trait
position <- sp_gap %>% 
filter(genus_sp %in% fishbase_traits == F, str_detect(genus_sp, '\\ssp$') == F) %>%
  pull(genus_sp) %>% 
  ecology(species_list = .) %>% select(Species, Macrophyte:Driftwood, BoatHulls:HardCorals) %>%
  tidyr::pivot_longer(., cols = Macrophyte:HardCorals, names_to = "Habitat", values_to = "YN") %>%
  filter(YN == -1)

fishbase_traits <- position %>% group_by(Species) %>% summarise(Habitat = paste(Habitat, collapse = ', ')) %>%
  full_join(fishbase_traits, ., by = "Species")
fishbase_traits$Genus <- word(fishbase_traits$Species, 1)
rm(position)  

fishbase_traits <- select(fishbase_traits, !(Length:LTypeMaxF)) %>% relocate(Genus, .before = Species)

## Merge trophic guilds ==========================================================

# Morais and Bellwood don't use the same trophic guild groupings as Parravicini
# I'd rather predict Kmax from these groupings
# Match extrapolated guild predictions to rfishprod dt_prodDB
dt_troph <- read.csv("src/parravicini_trophic_guilds.csv") %>% 
  select(species, trophic_guild_predicted) %>% mutate(species = str_replace(species, '\\_', ' ')) %>%
  rename(Species = species)
replace <- c('InvSes', 'Herb', 'Cor', 'Pisc', 'InvMic', 'InvMac', 'Crus', 'Plank')
names(replace) <- as.character(1:8)
dt_troph$Diet_Par <- str_replace_all(dt_troph$trophic_guild_predicted |> as.character(), replace)
fishbase_traits <- left_join(fishbase_traits, dt_troph, by = "Species")
rm(replace)

# Export compiled traits
# readr::write_csv(fishbase_traits, 'src/fish_prodtraits.csv')

# Make rfishprod training diet data match
dt_prodDB %>% distinct(Species, Diet) %>% left_join(., dt_troph[-2], by = "Species") %>% View
dt_prodDB %>% distinct(Species, Diet) %>% left_join(., dt_troph[-2], by = "Species") %>%
  pull(Diet_Par) %>% is.na %>% sum # how many species don't have a match with Parravicini guilds
dt_prodDB <- dt_prodDB %>% left_join(., dt_troph[-2], by = "Species")
dt_prodDB %>% filter(is.na(Diet_Par)) %>% distinct(Diet)

# manual matches
dt_prodDB$Diet_Par[is.na(dt_prodDB$Diet_Par) & str_detect(dt_prodDB$Diet, '^Her')] <- 'Herb'
dt_prodDB$Diet_Par[is.na(dt_prodDB$Diet_Par) & dt_prodDB$Diet == 'FisCep'] <- 'Pisc'
dt_prodDB$Diet_Par[is.na(dt_prodDB$Diet_Par) & dt_prodDB$Diet == 'Plktiv'] <- 'Plank'
dt_prodDB$Diet_Par[is.na(dt_prodDB$Diet_Par) & dt_prodDB$Diet == 'FisCep'] <- 'Pisc'
# can't resolve invertivores like this
dt_prodDB_sp_troph <- dt_prodDB %>% filter(is.na(Diet_Par)) %>% pull(Species)
# worth manual matching? yes, about > 20% of training data
dt_prodDB %>% group_by(Diet) %>% summarise(n())
length(dt_prodDB_sp_troph)
dt_prodDB_sp_troph <- fooditems(species_list = dt_prodDB_sp_troph)
View(dt_prodDB_sp_troph)

# have to do it manually off-code
# readr::write_csv(dt_prodDB_sp_troph, 'src/missing_troph_sp.csv')
dt_prodDB_sp_troph <- readr::read_csv('src/missing_troph_sp.csv') # read back in
# merge in the manual records
dt_troph <- dt_troph %>% select(!trophic_guild_predicted) %>%
  bind_rows(., dt_prodDB_sp_troph)
# rejoin
dt_prodDB <- dt_prodDB %>% select(!c(Diet_Par, Diet)) %>% left_join(., dt_troph, by = "Species") %>%
  filter(is.na(Diet_Par) == F)
summary(dt_prodDB)

# restash the species that have complete compiled traits
dt_traits <- dt_prodDB %>% rename(Diet = Diet_Par) %>%
select(Species, MaxSizeTL, Diet, Position, a, b) %>%
  filter(Species %in% sp_gap$genus_sp) %>% distinct

# Compile SST ==========================================================

years <- unique(dt_assemblage$Year)

# Estimate Kmax ==========================================================
# estimate VGBF Kmax
# input trait data
# dataset = reference data to train

# Formula from Morais and Bellwood (2018)
fmod <- formula(~ sstmean + MaxSizeTL + Diet + Position + Method)
dt_kmax <- predKmax(repdata, dataset = db,
                   fmod = fmod, niter = iter,
                   return = 'pred')

