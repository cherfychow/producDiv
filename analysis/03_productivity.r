
# ProducDiv
# Calculate productivity

require(dplyr)
require(stringr)
require(sf)
require(rfishprod)
require(rfishbase)


# Compile traits ==========================================================
# load("/Users/cherchow/Dropbox/ch4_producDiv/datarepo/workspace.RData")

sp_gap <- dt_assemblage %>% distinct(Family, Genus, Species) %>%
  mutate(genus_sp = case_when(is.na(Genus) == F ~ paste(Genus, Species, sep = ' ') %>%
	str_remove(., '\\d+'), .default = NA)) %>% 
  # collapse known uncertains because they don't matter for productivity
  arrange(Family, genus_sp) %>% distinct %>%
  filter(!is.na(genus_sp)) # remove NAs
View(sp_gap)

prod_traits <- db %>% select(Species, MaxSizeTL, Diet, Position, a, b) %>%
  filter(Species %in% sp_gap$genus_sp) %>% distinct


# fetch traits from FishBase for species not covered
# L-W parameters, need to scrape
# MaxSizeTL
# Diet/Trophic group
# Position
# SSTMean

# Lmax + water column
fishbase_traits <- sp_gap %>% 
  filter(genus_sp %in% prod_traits == F, str_detect(genus_sp, '\\ssp$') == F) %>%
  pull(genus_sp) %>% species(species_list = ., 
  fields = c('Species', 'Length', 'LTypeMaxM', 'LengthFemale', 'LTypeMaxF', 'DemersPelag'))

# Bayesian L-W parameters
fishbase_traits <- sp_gap %>% 
  filter(genus_sp %in% prod_traits == F, str_detect(genus_sp, '\\ssp$') == F) %>%
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
dt_troph <- read.csv("src/parravicini_trophic_guilds.csv") %>% 
  select(species, trophic_guild_predicted) %>% mutate(species = str_replace(species, '\\_', ' ')) %>%
  rename(Species = species)
replace <- c('InvSes', 'Herb', 'Cor', 'Pisc', 'InvMic', 'InvMac', 'Crus', 'Plank')
names(replace) <- as.character(1:8)
dt_troph$Diet_Par <- str_replace_all(dt_troph$trophic_guild_predicted |> as.character(), replace)

db %>% distinct(Species, Diet) %>% left_join(., dt_troph[-2], by = "Species") %>% View
db %>% distinct(Species, Diet) %>% left_join(., dt_troph[-2], by = "Species") %>%
  pull(Diet_Par) %>% is.na %>% sum # how many species don't have a match with Parravicini guilds
db <- db %>% distinct(Species, Diet) %>% left_join(., dt_troph[-2], by = "Species")
db %>% filter(is.na(Diet_Par)) %>% distinct(Diet)

# manual matches
db$Diet_Par[is.na(db$Diet_Par) & str_detect(db$Diet, '^Her')] <- 'Herb'
db$Diet_Par[is.na(db$Diet_Par) & db$Diet == 'FisCep'] <- 'Pisc'
db$Diet_Par[is.na(db$Diet_Par) & db$Diet == 'Plktiv'] <- 'Plank'
db$Diet_Par[is.na(db$Diet_Par) & db$Diet == 'FisCep'] <- 'Pisc'
# can't resolve invertivores like this
db_sp_troph <- db %>% filter(is.na(Diet_Par)) %>% pull(Species)
# worth manual matching? yes, about > 20% of training data
db %>% group_by(Diet) %>% summarise(n())
length(db_sp_troph)
db_sp_troph <- fooditems(species_list = db_sp_troph)
View(db_sp_troph)

# have to do it manually off-code
# readr::write_csv(db_sp_troph, 'src/missing_troph_sp.csv')

# Export compiled traits
# readr::write_csv(fishbase_traits, 'src/fish_prodtraits.csv')

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

