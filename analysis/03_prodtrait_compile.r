
# ProducDiv
# Compile traits for productivity modelling

require(dplyr)
require(readr)
require(stringr)
require(rfishprod)
require(rfishbase)

# SETUP ============================================================================================

dt_prodDB <- rfishprod::db
dt_assemblage <- bind_rows(dt_rareass, dt_single) # write over with rarefied
dt_assemblage %>% filter(is.na(fishbase_valid), taxon_resolution == 'species') %>% 
  distinct(worms_valid, fishbase_valid) %>% View
dt_rareass$fishbase_valid[dt_rareass$worms_valid == 'Haemulon atlanticus'] <- 'Haemulon atlanticus'

# species list to compile traits for
allsp <- dt_assemblage %>% 
  filter(taxon_resolution %in% c('genus', 'species')) %>%
  distinct(fishbase_valid) %>%
  mutate(fishbase_valid = str_remove(fishbase_valid, '\\d+')) %>% 
  # collapse known uncertains because they don't matter for productivity
  pull(fishbase_valid) %>% na.omit %>% unique %>% sort # remove NAs
# View(allsp)

# traits compiled in earlier tests
dt_traits <- read_csv('./src/manual_traits_sandbox.csv')
allsp[-which(allsp %in% dt_traits$Species|allsp %in% dt_prodDB$Species)] %>% length
# use up existing sources first

# FILLING IN TRAITS WE ALREADY HAVE ================================================================

# Parravicini diet traits
dt_troph <- read.csv("src/parravicini_trophic_guilds.csv") %>% 
  select(species, trophic_guild_predicted) %>% mutate(species = str_replace(species, '\\_', ' ')) %>%
  rename(Species = species)
replace <- c('InvSes', 'Herb', 'Cor', 'Pisc', 'InvMic', 'InvMac', 'Crus', 'Plank')
names(replace) <- as.character(1:8)
dt_troph$Diet_Par <- str_replace_all(dt_troph$trophic_guild_predicted |> as.character(), replace)
dt_troph <- dt_troph[-2]
rm(replace)

# From Morais and Bellwood 2020
# Change db diet to Parravicini framework
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

# table to fill in with trait info
trait_compile <- tibble(Species = allsp[-which(allsp %in% dt_traits$Species)])
trait_compile <- dt_prodDB %>% select(Species, MaxSizeTL, Diet, Position, a, b) %>% distinct %>% 
  left_join(trait_compile, ., by = "Species")
trait_compile <- left_join(trait_compile, dt_troph, by = "Species")

# FishBase extraction
# fetch the Biology description

# species description
fishbase_traits <- trait_compile %>% filter(is.na(MaxSizeTL)) %>% pull(Species) %>% 
  species(species_list = ., fields = c('Species', 'Comments'))

# Vertical position traits from GASPAR (Floeter et al. )
gaspar <- read_csv('src/gaspar_traits.csv')
unique(dt_prodDB$Position)
replace <- c('Bnth', 'BtPl', 'Pelg')
names(replace) <- c('Bottom', 'Low', 'High')
gaspar$Position_v <- str_replace_all(gaspar$Level_water, replace)
fishbase_traits <- gaspar %>% select(FullSpecies, Position_v) %>% 
  rename(Species = FullSpecies) %>%
  left_join(fishbase_traits, ., by = "Species")

# export out to manually finish compiling
# write_csv(fishbase_traits, './src/fishbase_traits_extract.csv')
# write_csv(trait_compile %>% filter(!str_detect(Species, ' sp')), './src/traits_tofill.csv')

## TROPH =============================================================================

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
db_sp_troph <- readr::read_csv('src/manual_missing_troph.csv') # read back in
# merge in the manual records
dt_troph <- bind_rows(dt_troph, db_sp_troph)
# rejoin
dt_prodDB <- dt_prodDB %>% select(!c(Diet_Par, Diet)) %>% left_join(., dt_troph, by = "Species") %>%
  filter(is.na(Diet_Par) == F)
summary(dt_prodDB)
dt_prodDB <- rename(dt_prodDB, Diet = Diet_Par)

rm(gaspar, db_sp_troph)


## REMERGE ====================================================================

rm(fishbase_traits, trait_compile) # remove the objects that helped with manual filling

traits_full <- read_csv('./src/manual_traits_full.csv')
traits_full <- distinct(traits_full)
# make sure it doesn't overlap with dt_traits
traits_full <- traits_full %>% filter(!Species %in% dt_traits$Species)
traits_full$Diet <- NULL

# manual name fixes: FISHBASE
replace = c('Bodianus pulcher' = 'Semicossyphus pulcher',
            'Carangoides gymnostethus' = 'Turrum gymnostethus',
            'Compsanthias hawaiiensis' = 'Pseudanthias hawaiiensis',
            'Compsanthias ventralis' = 'Pseudanthias ventralis',
            'Neosynchiropus ocellatus' = 'Synchiropus ocellatus',
            'Ostracion cubicus' = 'Ostracion cubicum',
            'Oxycheilinus digrammus' = 'Oxycheilinus digramma',
            'Plectorhinchus multivittatum' = 'Plectorhinchus multivittatus',
            'Plectroglyphidodon lacrymatus' = 'Stegastes lacrymatus',
            'Premnas biaculeatus' = 'Amphiprion biaculeatus',
            'Sufflamen fraenatus' = 'Sufflamen fraenatum',
            'Tylosurus crocodilus crocodilus' = 'Tylosurus crocodilus',
            'Astroscopus y-graecum' = 'Astroscopus ygraecum',
            'Cantherhines dumerili' = 'Cantherhines dumerilii',
            'Carangoides orthogrammus' = 'Ferdauia orthogrammus',
            'Chromis agilis' = 'Pycnochromis agilis',
            'Chromis amboinensis' = 'Pycnochromis amboinensis',
            'Chromis atripes' = 'Pycnochromis atripes',
            'Chromis iomelas' = 'Pycnochromis iomelas',
            'Chromis margaritifer' = 'Pycnochromis margaritifer',
            'Chromis multilineata' = 'Azurina multilineata',
            'Chromis retrofasciata' = 'Pycnochromis retrofasciatus',
            'Chromis vanderbilti' = 'Pycnochromis vanderbilti',
            'Cromileptes altivelis' = 'Chromileptes altivelis',
            'Grammatorycnus' = 'Grammatorcynus',
            'Plectorhinchus unicolor' = 'Plectorhinchus schotaf',
            'Pseudanthias bartlettorum' = 'Nemanthias bartlettorum',
            'Pseudanthias dispar' = 'Nemanthias dispar')

# manual name fixes: WoRMS
replace_worms = c('Centropyge bispinosus' = 'Centropyge bispinosa',
                  'Centropyge flavissimus' = 'Centropyge flavissima',
                  'Cantherhines dumerili' = 'Cantherhines dumerilii',
                  'Archamia fucata' = 'Taeniamia fucata',
                  'Chromis agilis' = 'Pycnochromis agilis',
                  'Chromis lepidolepis' = 'Azurina lepidolepis',
                  'Chromis lineata' = 'Pycnochromis lineatus',
                  'Plectroglyphidodon lacrymatus' = 'Stegastes lacrymatus',
                  'Oxycheilinus digrammus' = 'Oxycheilinus digramma',
                  'Tylosurus crocodilus crocodilus' = 'Tylosurus crocodilus',
                  'Sufflamen fraenatus' = 'Sufflamen fraenatum',
                  'Chromis agilis' = 'Pycnochromis agilis',
                  'Chromis amboinensis' = 'Pycnochromis amboinensis',
                  'Chromis atripes' = 'Pycnochromis atripes',
                  'Chromis iomelas' = 'Pycnochromis iomelas',
                  'Chromis margaritifer' = 'Pycnochromis margaritifer',
                  'Chromis multilineata' = 'Azurina multilineata',
                  'Chromis retrofasciata' = 'Pycnohromis retrofasciatus',
                  'Chromis vanderbilti' = 'Pycnohromis vanderbilti',
                  'Grammatorycnus' = 'Grammatorcynus',
                  'Plectorhinchus unicolor' = 'Plectorhinchus schotaf')

# replace these manual catches first
dt_assemblage %>% distinct(fishbase_valid, worms_valid) %>% filter(worms_valid %in% names(replace_worms))
dt_assemblage %>% distinct(fishbase_valid, worms_valid) %>% filter(fishbase_valid %in% names(replace))
dt_assemblage %>% distinct(fishbase_valid, worms_valid) %>% filter(worms_valid %in% names(replace))

# in main data
# fishbase name function = match traits
dt_assemblage$worms_valid <- str_replace_all(dt_assemblage$worms_valid, replace_worms)
dt_assemblage$fishbase_valid <- str_replace_all(dt_assemblage$fishbase_valid, replace)
dt_rareass$worms_valid <- str_replace_all(dt_rareass$worms_valid, replace_worms)
dt_rareass$fishbase_valid <- str_replace_all(dt_rareass$fishbase_valid, replace)
dt_single$worms_valid <- str_replace_all(dt_single$worms_valid, replace_worms)
dt_single$fishbase_valid <- str_replace_all(dt_single$fishbase_valid, replace)
dt_assemblage$worms_valid[dt_assemblage$fishbase_valid == 'Oxycheilinus unifasciatus'] <- 'Oxycheilinus unifasciatus'
dt_rareass$worms_valid[dt_rareass$fishbase_valid == 'Oxycheilinus unifasciatus'] <- 'Oxycheilinus unifasciatus'
dt_single$worms_valid[dt_single$fishbase_valid == 'Oxycheilinus unifasciatus'] <- 'Oxycheilinus unifasciatus'
dt_assemblage$fishbase_valid <- str_remove(dt_assemblage$fishbase_valid, 'cf\\.*\\s')
dt_rareass$fishbase_valid <- str_remove(dt_rareass$fishbase_valid, 'cf\\.*\\s')
dt_single$fishbase_valid <- str_remove(dt_single$fishbase_valid, 'cf\\.*\\s')

# in traits
dt_traits$Species <- str_replace_all(dt_traits$Species, replace) # for fishbase compatibility
traits_full$Species <- str_replace_all(traits_full$Species, replace)
dt_traits$Species <- str_remove(dt_traits$Species, 'cf\\.*\\s')
traits_full$Species <- str_remove(traits_full$Species, 'cf\\.*\\s')

# check if stringr did its job
dt_assemblage %>% distinct(worms_valid, fishbase_valid) %>% filter(worms_valid %in% names(replace_worms))
dt_assemblage %>% distinct(worms_valid, fishbase_valid) %>% filter(fishbase_valid %in% names(replace))
dt_traits %>% filter(Species %in% names(replace))
traits_full %>% filter(Species %in% names(replace))

## MATCH L-W MAXLENGTH =================================================

# extract MaxLength, a, b estimates from FishBase
Lab <- estimate(traits_full$Species, fields = c('Species', 'MaxLengthTL', 'a', 'b'))
Lab <- distinct(Lab)
nrow(traits_full) - nrow(Lab) # 36 species not matching with fishbase

missing_Lab <- traits_full %>% filter(Species %in% Lab$Species == F) %>% select(Species)
missing_Lab$fb_valid <- validate_names(missing_Lab$Species)
# some species need revalidating because fishbase contradicts itself front and backend
missing_Lab$fb_valid[missing_Lab$fb_valid == 'Genus Species'] = NA
missing_Lab$fb_valid[str_detect(missing_Lab$Species, 'Pycno')] <- missing_Lab$Species[str_detect(missing_Lab$Species, 'Pycno')]
missing_Lab$fb_valid <- str_replace(missing_Lab$fb_valid, 'Pycnohromis', 'Pycnochromis')
missing_Lab$fb_valid[missing_Lab$Species == 'Holacanthus townsendi'] = 'Holacanthus bermudensis'
missing_Lab$fb_valid[missing_Lab$Species == 'Pomacentrus maafu'] = 'Pomacentrus moluccensis'

# for the purposes of trait matching, open nomenclature will either go to closest relative values or Genus level
missing_Lab$fb_valid[str_which(missing_Lab$Species, 'aff\\.|hybrid')] <- missing_Lab$Species[str_which(missing_Lab$Species, 'aff\\.|hybrid')] %>%
  str_replace(., 'aff\\.|hybrid', 'sp')

# match
missing_Lab <- estimate(missing_Lab$fb_valid, fields = c('Species', 'MaxLengthTL', 'a', 'b')) %>% 
  left_join(missing_Lab, ., by = c('fb_valid' = 'Species'))

rm(replace, replace_worms)
replace_fb <- missing_Lab$fb_valid
names(replace_fb) <- missing_Lab$Species

# fix these names in main and trait data
dt_assemblage$fishbase_valid <- str_replace_all(dt_assemblage$fishbase_valid, replace_fb)
dt_rareass$fishbase_valid <- str_replace_all(dt_rareass$fishbase_valid, replace_fb)
dt_single$fishbase_valid <- str_replace_all(dt_single$fishbase_valid, replace_fb)
dt_assemblage$worms_valid <- str_replace(dt_assemblage$worms_valid, 'Pycnohromis', 'Pycnochromis')
dt_rareass$worms_valid <- str_replace(dt_rareass$worms_valid, 'Pycnohromis', 'Pycnochromis')
dt_single$worms_valid <- str_replace(dt_single$worms_valid, 'Pycnohromis', 'Pycnochromis')
traits_full$Species <- str_replace_all(traits_full$Species, replace_fb)
traits_full <- distinct(traits_full)

# merge with the full Lab table
Lab <- missing_Lab %>% select(!Species) %>% rename(Species = fb_valid) %>% 
  filter(Species %in% Lab$Species == F) %>% 
  bind_rows(Lab, .) %>% 
  arrange(Species) %>% distinct
summary(Lab)

# merge with manually compiled diet and position traits
traits_full <- left_join(traits_full, Lab, by = 'Species')
rm(Lab, missing_Lab, replace_fb)

# check if there's any overlap with dt_traits
inner_join(traits_full, dt_traits, by = 'Species')
traits_full %>% filter(Species %in% dt_traits$Species) # yep 16 species overlapped

# match in
traits_full <- traits_full %>% filter(!Species %in% dt_traits$Species) %>% 
  rename(MaxSizeTL = MaxLengthTL, Diet = Diet_Par) %>% 
  mutate(Genus = word(Species, 1))

genfam <- dt_assemblage %>% distinct(fishbase_valid, worms_valid, genus_valid, family_valid) %>% 
  filter(fishbase_valid %in% traits_full$Species) %>% distinct(genus_valid, family_valid)
traits_full <- left_join(traits_full, genfam, by = c('Genus' = 'genus_valid'))
traits_missing <- traits_full %>% filter(is.na(a) | is.na(MaxSizeTL)) # fishbase hates anthias
traits_full <- traits_full %>% filter(!is.na(a))
summary(traits_full)

# give Chromis anadema the l-w for C. xanhura
traits_missing[7,5] <- 0.01413
traits_missing[7,6] <- 2.94
traits_missing <- traits_missing %>% filter(!str_detect(Species, '\\ssp$'))
# Handle Anthiadidae
traits_missing$family_valid[traits_missing$family_valid != 'Pomacentridae' | is.na(traits_missing$family_valid)] <- 'Anthiadidae'
# only 17 values for the whole family... give family level mean
traits_missing <- traits_missing %>% mutate(a = case_when(is.na(a) ~ 0.02515, .default = a),
                                           b = case_when(is.na(b) ~ 2.96048375, .default = b))

# merge back in
traits_full <- bind_rows(traits_full, traits_missing) %>% arrange(Species)
summary(traits_full)

# merge in with the FULL trait table
traits_full %>% rename(Family = family_valid) %>% bind_rows(dt_traits, .) %>% arrange(Species) %>% summary
dt_traits <- traits_full %>% rename(Family = family_valid) %>% bind_rows(dt_traits, .) %>% arrange(Species)

dt_traits %>% group_by(Species) %>% mutate(n = n()) %>% filter(n > 1) %>% View
# duplicates because worms and fishbase don't agree on some families
dt_traits <- dt_traits %>% select(!Family) %>% distinct
dt_traits <- dt_traits[-which(dt_traits$Species == 'Cantherhines dumerilii' & dt_traits$Diet == 'InvSes'),]

## GENUS LEVEL TRAITS =================================================

allsp <- dt_assemblage %>% distinct(worms_valid, fishbase_valid, taxon_resolution)
# allsp %>% filter(is.na(fishbase_valid)) %>% View

# dt_assemblage %>% filter(taxon_resolution == 'species', is.na(fishbase_valid)) %>%
#   distinct(worms_valid, fishbase_valid) %>% View

dt_assemblage <- mutate(dt_assemblage, 
                        fishbase_valid = case_when(worms_valid == 'Haemulon atlanticus' ~ 'Haemulon atlanticus',
                                                   worms_valid == 'Pomacentrus maafu' ~ 'Pomacentrus moluccensis',
                                                   worms_valid == 'Holacanthus townsendi' ~ 'Holacanthus bermudensis',
                                                   worms_valid == 'Hypoplectrus hybrid' ~ 'Hypoplectrus sp',
                                                   worms_valid == 'Malacoctenus aff.' ~ 'Malacoctenus sp',
                                                   worms_valid == 'Menaphorus punticulatus' ~ NA,
                                                   is.na(fishbase_valid) ~ worms_valid,
                                                   .default = fishbase_valid))
dt_rareass <- mutate(dt_rareass, 
                        fishbase_valid = case_when(worms_valid == 'Haemulon atlanticus' ~ 'Haemulon atlanticus',
                                                   worms_valid == 'Pomacentrus maafu' ~ 'Pomacentrus moluccensis',
                                                   worms_valid == 'Holacanthus townsendi' ~ 'Holacanthus bermudensis',
                                                   worms_valid == 'Hypoplectrus hybrid' ~ 'Hypoplectrus sp',
                                                   worms_valid == 'Malacoctenus aff.' ~ 'Malacoctenus sp',
                                                   worms_valid == 'Menaphorus punticulatus' ~ NA,
                                                   is.na(fishbase_valid) ~ worms_valid,
                                                   .default = fishbase_valid))

dt_single <- mutate(dt_single, 
                        fishbase_valid = case_when(worms_valid == 'Haemulon atlanticus' ~ 'Haemulon atlanticus',
                                                   worms_valid == 'Pomacentrus maafu' ~ 'Pomacentrus moluccensis',
                                                   worms_valid == 'Holacanthus townsendi' ~ 'Holacanthus bermudensis',
                                                   worms_valid == 'Hypoplectrus hybrid' ~ 'Hypoplectrus sp',
                                                   worms_valid == 'Malacoctenus aff.' ~ 'Malacoctenus sp',
                                                   worms_valid == 'Menaphorus punticulatus' ~ NA,
                                                   is.na(fishbase_valid) ~ worms_valid,
                                                   .default = fishbase_valid))
allsp <- dt_assemblage %>% distinct(worms_valid, fishbase_valid, taxon_resolution)
anti_join(allsp, dt_traits, by = c('fishbase_valid' = 'Species')) %>% View

rm(traits_missing)
# still some species missing?!

missing <- read_csv('src/traits_missing.csv')
dt_traits <- missing %>% select(!worms_valid) %>% rename(Species = fishbase_valid) %>% 
  bind_rows(dt_traits, .)
rm(missing)

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

## FAMILY LEVEL TRAITS =================================================

family <- sort(unique(dt_assemblage$genus_sp)) %>% str_remove(., '\\d') %>% 
  unique
family <- family[-which(family %in% dt_traits$Species)] %>% word(., 1) %>% unique

traits_family <- dt_assemblage %>% distinct(Family, Genus) %>% 
  left_join(dt_traits, ., by = "Genus") %>% 
  filter(Family %in% family) %>% select(Family, MaxSizeTL, a, b) %>%
  group_by(Family) %>% 
  summarise(MaxSizeTL = mean(MaxSizeTL), a = mean(a), b = mean(b))

# frequency based matching for categorical
diet <- dt_assemblage %>% distinct(Family, Genus) %>% 
  left_join(dt_traits, ., by = "Genus") %>% 
  filter(Family %in% family) %>% ungroup %>% select(Family, Diet) %>%
  group_by_all %>% summarise(freq = n()) %>% arrange(Family, desc(freq)) %>%
  mutate(order = 1:n())

position <- dt_assemblage %>% distinct(Family, Genus) %>% 
  left_join(dt_traits, ., by = "Genus") %>% 
  filter(Family %in% family) %>% ungroup %>% select(Family, Position) %>%
  group_by_all %>% summarise(freq = n()) %>% arrange(Family, desc(freq)) %>%
  mutate(order = 1:n())

diet <- diet %>% filter(order == 1) %>% ungroup %>% select(!c(freq, order))
position <- position %>% filter(order == 1) %>% ungroup %>% select(!c(freq, order))

traits_family <- left_join(traits_family, diet, by = 'Family') %>%
  left_join(., position, by = "Family")

traits_family$Species <- traits_family$Family

dt_traits <- dt_assemblage %>% distinct(Family, Genus) %>% 
  left_join(dt_traits, ., by = "Genus") %>% bind_rows(., traits_family)

rm(traits_family, family, diet, position)

# readr::write_csv(dt_traits, './src/traits_final.csv')
