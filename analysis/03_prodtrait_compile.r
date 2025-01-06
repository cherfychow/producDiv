
# ProducDiv
# Compile traits for productivity modelling
# runs independently from the analysis pipeline

require(dplyr)
require(readr)
require(stringr)
require(rfishprod)
require(rfishbase)

# SETUP ============================================================================================

dt_sp <- read_csv('src/dt_species.csv')

dt_prodDB <- rfishprod::db
dt_sp %>% filter(is.na(fishbase_valid), taxon_resolution == 'species') %>% View
dt_sp %>% filter(taxon_resolution == 'species') %>% 
  arrange(fishbase_valid) %>% distinct(fishbase_valid) %>% View
dt_sp %>% filter(taxon_resolution == 'lowres') %>% View
dt_sp$fishbase_valid[dt_sp$taxon_resolution == 'lowres'] <- NA

# make a backbone to compile traits
# fishbase_valid as the matching key to taxaID
dt_traits_sp <- dt_sp %>% 
  filter(is.na(fishbase_valid) == F, taxon_resolution == 'species') %>% 
  distinct(fishbase_valid) %>% arrange(fishbase_valid)

dt_traits_g <- dt_sp %>% distinct(fishbase_valid) %>% 
  filter(is.na(fishbase_valid) == F, str_detect(fishbase_valid, 'sp$')) %>% 
  arrange(fishbase_valid)

dt_traits_f <- dt_sp %>%
  filter(is.na(fishbase_valid) == F, taxon_resolution == 'family') %>% 
  distinct(fishbase_valid) %>% arrange(fishbase_valid)

# traits compiled in earlier tests
sandbox_traits <- read_csv('./src/manual_traits_sandbox.csv')
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


## TRY TO MATCH BACK IN ====================================================================

rm(fishbase_traits, trait_compile) # remove the objects that helped with manual filling

traits_full <- read_csv('./src/manual_traits_full.csv')
traits_full <- distinct(traits_full)
# make sure it doesn't overlap with sandbox_traits
traits_full <- traits_full %>% filter(!Species %in% sandbox_traits$Species)

traits_full$Diet <- NULL # remove, it was just helping us fill in
# fishbase updated while we were compiling

traits_full <- traits_full %>% filter(!str_detect(Species, 'Menaphorus|aff\\.'))
traits_full$Species <- str_remove(traits_full$Species, 'cf\\.\\s|[:punct:]')
traits_full %>% filter(Species %in% dt_traits_sp$fishbase_valid == F) %>% View

# known hybrids/too new species
traits_full <- traits_full %>% filter(!str_detect(Species, 'maafu'))
traits_full$Species[traits_full$Species == 'Holacanthus townsendi'] <- 'Holacanthus bermudensis'

# try have fishbase revalidate, but kind of manual
temp <- tibble(fishbase_valid = validate_names(traits_full$Species))
temp <- temp[-1404,]
traits_full <- bind_cols(traits_full, temp)
traits_full %>% filter(is.na(fishbase_valid))
traits_full <- traits_full[-which(traits_full$Species == 'Hypoplectrus hybrid'),]
traits_full$Species[traits_full$Species == 'Bodianus pulcher'] <- 'Semicossyphus pulcher'
traits_full$fishbase_valid[is.na(traits_full$fishbase_valid)] <- traits_full$Species[is.na(traits_full$fishbase_valid)]
traits_full[is.na(traits_full$fishbase_valid),] # any left? great!
rm(temp)

traits_full %>% select(!Species) %>% distinct %>% nrow
traits_full %>% filter(fishbase_valid %in% dt_traits_sp$fishbase_valid) %>% nrow
dt_traits_sp$fishbase_valid[which(dt_traits_sp$fishbase_valid %in% traits_full$fishbase_valid == F & dt_traits_sp$fishbase_valid %in% sandbox_traits$Species == F)]

sandbox_traits$fishbase_valid <- validate_names(sandbox_traits$Species)
dt_traits_sp <- dt_traits_sp %>% filter(!str_detect(fishbase_valid, 'sp$'))
# species level records with genus level trait matching. just delete, covered by genus table

## MATCH L-W MAXLENGTH =================================================

missLab <- dt_traits_sp$fishbase_valid[dt_traits_sp$fishbase_valid %in% sandbox_traits$fishbase_valid == F]

# extract MaxLength, a, b estimates from FishBase
Lab <- estimate(missLab, fields = c('Species', 'MaxLengthTL', 'a', 'b')) %>% distinct
Lab %>% filter(is.na(MaxLengthTL)) %>% View
Lab %>% filter(is.na(a)) %>% View # only issues with length-weight params

failab <- Lab %>% filter(is.na(a)) # isolate to fix
Lab <- Lab %>% filter(!is.na(a))

# give Chromis anadema the l-w for C. xanhura
failab[which(failab$Species == 'Chromis anadema'),3:4] <- list(0.01413, 2.94) 
# is everything else anthias? yes
dt_sp %>% filter(fishbase_valid %in% failab$Species) %>% select(family_valid, fishbase_valid) %>% View

# Handle Anthiadidae
# only 17 values for the whole family... give family level mean
failab <- failab %>% mutate(a = case_when(is.na(a) ~ 0.02515, .default = a),
                            b = case_when(is.na(b) ~ 2.96048375, .default = b))
Lab <- bind_rows(Lab, failab) %>% arrange(Species) # merge missing in
# merge in sandbox MaxLength, a, b
colnames(Lab)[2] <- 'MaxSizeTL'
colnames(Lab)[1] <- 'fishbase_valid'
Lab <- sandbox_traits %>% select(fishbase_valid, MaxSizeTL, a, b) %>% 
  filter(is.na(fishbase_valid) == F) %>% # only the species for now
  bind_rows(., Lab) # merge
  
# join in with skeleton
dt_traits_sp <- left_join(dt_traits_sp, Lab, by = 'fishbase_valid')
summary(dt_traits_sp)
rm(Lab)
rm(failab)

# temp <- read_csv('src/missing2.csv')
# traits_full <- temp %>% select(fishbase_valid, Diet, Position) %>% bind_rows(traits_full, .)

# categorical traits match test
colnames(traits_full)[3] <- 'Diet'
dt_traits_sp <- sandbox_traits %>% filter(!str_detect(fishbase_valid, '\\ssp$|idae$')) %>% 
  select(fishbase_valid, Diet, Position) %>% # just categorical
  bind_rows(., traits_full[-1]) %>% # join with manual compiled traits
  left_join(dt_traits_sp, ., by = 'fishbase_valid')

## GENUS LEVEL TRAITS =================================================

sandbox_traits %>% filter(str_detect(Species, '\\ssp$')) %>% View
genus_query <- dt_traits_g %>% filter(fishbase_valid %in% sandbox_traits$Species == F) %>% 
  pull(fishbase_valid) %>% word(., 1)

# update species gap
# take the genus mean for these traits?
# fishbase needs the species table as an intermediate to query
fb_sp <- species_names(fields = c('SpecCode', 'Species', 'Genus')) %>% filter(Genus %in% genus_query)
Lab_g <- estimate(fb_sp$Species, fields = c('Species', 'MaxLengthTL', 'a', 'b'))
# calculate genus means
Lab_g <- Lab_g %>% mutate(Genus = word(Species, 1)) %>% group_by(Genus) %>% 
  summarise(MaxSizeTL = mean(MaxLengthTL, na.rm = T), 
            a = mean(a, na.rm = T), 
            b = mean(b, na.rm = T)) %>% ungroup
Lab_g %>% filter(is.na(a)) # anthias again
Lab_g$a[is.na(Lab_g$a)] <- 0.02515
Lab_g$b[is.na(Lab_g$b)] <- 2.96048375
Lab_g$fishbase_valid <- paste(Lab_g$Genus, 'sp')
Lab_g$Genus <- NULL
rm(fb_sp)

# merge in length traits
dt_traits_g <- sandbox_traits %>% filter(str_detect(Species, '\\ssp$')) %>%
  select(Species, MaxSizeTL, a, b) %>% rename(fishbase_valid = Species) %>% 
  bind_rows(., Lab_g) %>%
  left_join(dt_traits_g, ., by = 'fishbase_valid')

# frequency based matching for categorical
diet <- dt_traits_sp %>% mutate(Genus = word(fishbase_valid, 1)) %>% 
  filter(Genus %in% genus_query) %>% ungroup %>% 
  select(Genus, Diet) %>%
  group_by_all %>% summarise(freq = n()) %>% arrange(Genus, desc(freq)) %>%
  mutate(order = 1:n())

position <- dt_traits_sp %>% mutate(Genus = word(fishbase_valid, 1)) %>% 
  filter(Genus %in% genus_query) %>% ungroup %>%
  select(Genus, Position) %>%
  group_by_all %>% summarise(freq = n()) %>% 
  arrange(Genus, desc(freq)) %>%
  mutate(order = 1:n())

diet <- diet %>% filter(order == 1) %>% ungroup %>% select(!c(freq, order))
position <- position %>% filter(order == 1) %>% ungroup %>% select(!c(freq, order))

temp <- left_join(diet, position, by = "Genus")
temp$fishbase_valid <- with(temp, paste(Genus, 'sp', sep = ' '))
temp <- sandbox_traits %>% filter(str_detect(Species, '\\ssp$')) %>%
  select(Species, Diet, Position) %>% rename(fishbase_valid = Species) %>% 
  bind_rows(temp, .)

dt_traits_g <- left_join(dt_traits_g, temp[-1], by = 'fishbase_valid')
dt_traits_g %>% filter(is.na(Diet)) %>% View
dt_traits_g <- dt_traits_g %>% filter(!is.na(Diet)) # isolate to manually fill in

dt_traits_g <- read_csv('src/missing_g.csv') %>% bind_rows(., dt_traits_g)

## FAMILY LEVEL TRAITS =================================================

# fishbase still uses Scaridae
dt_sp$fishbase_valid[with(dt_sp, which(genus_valid == 'Scarini' & str_detect(taxon_resolution, 'family')))] <- 'Scaridae'

fb_sp <- species_names(fields = c('SpecCode', 'Species', 'Family')) %>% filter(Family %in% dt_traits_f$fishbase_valid)
Lab_f <- estimate(fb_sp$Species, fields = c('Species', 'MaxLengthTL', 'a', 'b'))
Lab_f <- left_join(Lab_f, fb_sp[c(2,5)], by = 'Species') %>% group_by(Family) %>% 
  summarise(MaxSizeTL = mean(MaxLengthTL, na.rm = T), 
            a = mean(a, na.rm = T), 
            b = mean(b, na.rm = T)) %>% ungroup

# any missing?
dt_traits_f %>% filter(fishbase_valid %in% Lab_f$Family == F) # good
colnames(Lab_f)[1] <- 'fishbase_valid'

# frequency based matching for categorical
diet <- dt_sp %>% distinct(family_valid, fishbase_valid) %>% 
  left_join(dt_traits_sp, ., by = "fishbase_valid") %>% 
  filter(family_valid %in% dt_traits_f$fishbase_valid) %>% ungroup %>% 
  select(family_valid, Diet) %>%
  group_by_all %>% summarise(freq = n()) %>% arrange(family_valid, desc(freq)) %>%
  mutate(order = 1:n())

position <- dt_sp %>% distinct(family_valid, fishbase_valid) %>% 
  left_join(dt_traits_sp, ., by = "fishbase_valid") %>% 
  filter(family_valid %in% dt_traits_f$fishbase_valid) %>% ungroup %>% 
  select(family_valid, Position) %>%
  group_by_all %>% summarise(freq = n()) %>% arrange(family_valid, desc(freq)) %>%
  mutate(order = 1:n())

diet <- diet %>% filter(order == 1) %>% ungroup %>% select(!c(freq, order))
position <- position %>% filter(order == 1) %>% ungroup %>% select(!c(freq, order))

dt_traits_f %>% filter(fishbase_valid %in% diet$family_valid == F) # ok manual 2

diet <- left_join(diet, position, by = 'family_valid')
diet <- tibble(family_valid = c('Scaridae', 'Clupeidae'),
               Diet = c('Herb', 'Plank'),
               Position = c('BtPlDw', 'PelgAs')) %>% bind_rows(diet, .)

# merge in
dt_traits_f <- diet %>% rename(fishbase_valid = family_valid) %>% 
  left_join(dt_traits_f, ., by = 'fishbase_valid') %>% 
  left_join(., Lab_f, by = 'fishbase_valid')

rm(diet, position, Lab_f, Lab_g, fb_sp, dt_troph, sandbox_traits, temp)


# FINAL JOIN! ======================================================================================

summary(dt_traits_sp)
summary(dt_traits_g)
summary(dt_traits_f)

dt_traits <- bind_rows(dt_traits_sp, dt_traits_g, dt_traits_f)

# some final checks
dt_sp %>% filter(fishbase_valid %in% dt_traits$fishbase_valid == F) %>% View # not included?
n_distinct(dt_traits$fishbase_valid) # duplicates?
dt_traits <- distinct(dt_traits)
dt_traits %>% group_by(fishbase_valid) %>% mutate(n = n()) %>% filter(n > 1) %>% View
# conflicting sources?
dt_traits <- dt_traits %>% distinct(fishbase_valid, .keep_all = T)

# readr::write_csv(dt_traits, './src/traits_final.csv')
# write_csv(dt_sp, './src/dt_species_clean.csv')