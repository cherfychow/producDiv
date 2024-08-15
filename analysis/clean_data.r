

require(dplyr)
require(readr)
require(taxize)
require(rfishbase)
require(worrms)
require(stringr)
require(purrr)

# load in data
dt_assemblage <- read_csv('src/macro_reef_fish_obs.csv')
dt_meta <- read_csv('./src/macro_reef_fish_meta.csv') # dataset level metadata
gc() # keeping memory usage low

summary(dt_assemblage)

# TAXONOMY CLEAN ===================================================================================

# dt_assemblage %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>% View

# clean epithets
dt_assemblage$Species[with(dt_assemblage, which(Genus == 'Eviota' & str_detect(Species, '^cf')))] <- str_replace(dt_assemblage$Species[with(dt_assemblage, which(Genus == 'Eviota' & str_detect(Species, '^cf')))], 'cf', 'cf. ')
# remove genera, but allow the cf. records
dt_assemblage$Genus[which(str_detect(dt_assemblage$Species, '\\s') & str_detect(dt_assemblage$Species, '^cf') == F)] <- word(dt_assemblage$Species[which(str_detect(dt_assemblage$Species, '\\s') & str_detect(dt_assemblage$Species, '^cf') == F)], 1)
dt_assemblage$Species[str_which(dt_assemblage$Species, '^cf', negate = T)] <- str_remove(dt_assemblage$Species[str_which(dt_assemblage$Species, '^cf', negate = T)], '^[A-Z]{1}[a-z]+\\s')
dt_assemblage$Species[str_which(dt_assemblage$Species, '^cf', negate = T)] <- word(dt_assemblage$Species[str_which(dt_assemblage$Species, '^cf', negate = T)], 1) # just epithets, no subspecies

# uncertains
# dt_assemblage %>% distinct(Genus, Species) %>% filter(str_detect(Species, '^sp.{1,2}$')) %>% View
dt_assemblage$Species[dt_assemblage$Species %in% c('sp.', 'spp', 'species', 'spies')] <- 'sp'
dt_assemblage$Species[dt_assemblage$Species == 'spx'] <- 'sp1'
dt_assemblage$Species[dt_assemblage$Species == 'spz'] <- 'sp2'
dt_assemblage$Species[dt_assemblage$Species == 'spd'] <- 'sp1'
dt_assemblage$Family[dt_assemblage$Genus == 'Blenny'] <- 'Blenniidae'
dt_assemblage$Species[dt_assemblage$Genus == 'Blenny'] <- NA
dt_assemblage$Genus[dt_assemblage$Genus == 'Blenny'] <- NA
dt_assemblage$Species <- str_replace(dt_assemblage$Species, 'i{3}$', 'ii')
dt_assemblage$Species[which(dt_assemblage$Species == 'red')] <- 'sp1'

# parrotfish specific
dt_assemblage$Species[dt_assemblage$Genus == 'Scarinissp8'] <- 'sp8'
dt_assemblage$Genus[dt_assemblage$Genus == 'Scarinissp8'] <- 'Scarini'
dt_assemblage$Species[str_which(dt_assemblage$Species, '^\\d$')] <- paste0('sp', dt_assemblage$Species[str_which(dt_assemblage$Species, '^\\d$')])
dt_assemblage$Genus[dt_assemblage$Family == 'Scaridae' & is.na(dt_assemblage$Genus)] <- 'Scarini'
dt_assemblage$Family[dt_assemblage$Family == 'Scaridae'] <- 'Labridae'
dt_assemblage$Genus <- str_to_title(dt_assemblage$Genus)
dt_assemblage$Species <- str_to_lower(dt_assemblage$Species)
dt_assemblage$Genus[dt_assemblage$Genus == 'Labrid'] <- NA

# omits
dt_assemblage <- dt_assemblage[-which(str_detect(dt_assemblage$Species, 'crcp|franc') & is.na(dt_assemblage$Genus)),]
dt_assemblage <- dt_assemblage[-which(dt_assemblage$Genus == 'Teleostei'),]


# make a column to help with sorting uncertain records
dt_assemblage$taxon_resolution <- NA
# do the normals first
dt_assemblage$taxon_resolution[with(dt_assemblage, which(is.na(Family) == F & is.na(Genus)))] <- 'family'
dt_assemblage$taxon_resolution[with(dt_assemblage, which(is.na(Genus) == F & str_detect(Species, '^sp$')))] <- 'genus'
dt_assemblage$taxon_resolution[with(dt_assemblage, which(is.na(Genus) == F & is.na(Species) == F))] <- 'species'

# now the exceptions
dt_assemblage$taxon_resolution[with(dt_assemblage, which(Genus == 'Scarini'))] <- 'family'
dt_assemblage$taxon_resolution[with(dt_assemblage, which(is.na(Genus) & str_detect(Species, '^sp\\d{1,2}$')))] <- 'family specified'
dt_assemblage$taxon_resolution[with(dt_assemblage, which(is.na(Genus) == F & str_detect(Species, '^sp\\d{1,2}$')))] <- 'genus specified'
# Scarini special
dt_assemblage$taxon_resolution[with(dt_assemblage, which(Genus == 'Scarini' & str_detect(Species, '^sp\\d{1,2}$')))] <- 'family specified'

# reorg
dt_assemblage <- relocate(dt_assemblage, taxon_resolution, .after = Species)

# family level consistency, specified from unspecified
dt_assemblage$Species[which(dt_assemblage$taxon_resolution == 'family')] <- NA

dt_assemblage %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>%  View
dt_assemblage %>% filter(str_detect(Species, '\\s')) %>% View

## Automate misspelling fix========================================================================

# generate combined taxonomy
dt_assemblage <- dt_assemblage %>% mutate(genus_sp = case_when(
  taxon_resolution == 'family' ~ Family,
  taxon_resolution == 'family specified' ~ paste(Family, Species, sep = ' '),
  .default = paste(Genus, Species, sep = ' ')
)) %>% relocate(genus_sp, .before = Family)

# generate species list for checking
sp <- dt_assemblage %>% # only certain, species level IDs for this stage
  filter(taxon_resolution == 'species' & str_detect(Species, '\\s|[:punct:]') == F) %>% 
  distinct(genus_sp) %>% pull %>% sort

# fetch taxonomic IDs from fuzzy matching
dt_spcheck <- tibble(genus_sp = sp)
# GBIF is smarter about misspellings than WoRMS, but requires more user input
dt_spcheck$worms_id <- get_wormsid(dt_spcheck$genus_sp, searchtype = 'scientific', fuzzy = T,
                                   marine_only = T, messages = F)
dt_spcheck$gbif_id <- NA
dt_spcheck$gbif_id[is.na(dt_spcheck$worms_id)] <- get_gbifid(dt_spcheck$genus_sp[is.na(dt_spcheck$worms_id)], 
                                                             messages = F, method = 'backbone')

# first round of catching unaccepted names
dt_spcheck %>% filter(is.na(worms_id) & is.na(gbif_id)) # ok only three to fix manually

# get names. filter out synonyms
dt_spcheck <- dt_spcheck %>% filter(is.na(worms_id) == F) %>%
  mutate(across(ends_with('id'), as.numeric))
worms <- wm_record_(dt_spcheck$worms_id |> na.omit())  %>% bind_rows # will take a loong time
gc()
dt_spcheck$worms_id <- as.numeric(dt_spcheck$worms_id)
dt_spcheck <- worms %>% select(AphiaID, scientificname, status, valid_name, class, family, genus) %>% 
  rename(worms_valid = valid_name, worms_status = status, worms_name = scientificname) %>% distinct %>% 
  left_join(dt_spcheck, ., by = c('worms_id' = 'AphiaID'))

dt_spcheck %>% filter(is.na(worms_id)) %>% View
gbif <- id2name(dt_spcheck$gbif_id |> na.omit(), db = 'gbif')
gbif <- bind_rows(gbif[1:length(gbif)])
dt_spcheck2 <- gbif %>% select(!rank) %>% 
  rename(gbif_name = name, gbif_status = status, gbif_id = id) %>% distinct %>% 
  left_join(dt_spcheck, ., by = "gbif_id")
dt_spcheck <- dt_spcheck2

# gbif is really a middleman for better fuzzy matching misspellings
# put them back through worms

worms2 <- wm_records_names(dt_spcheck$gbif_name |> na.omit(), marine_only = T) %>% bind_rows
# isolate to then rejoin
gbif_match <- dt_spcheck %>% filter(is.na(gbif_id) == F) %>% select(genus_sp, gbif_id, gbif_name)
# remove bad gbif matches
gbif_match <- gbif_match %>% filter(str_detect(gbif_name, '\\s'))
gbif_match <- worms2 %>% select(scientificname, status, valid_name, class, family, genus) %>% 
  rename(worms_valid = valid_name, worms_status = status, worms_name = scientificname) %>% 
  distinct %>% left_join(gbif_match, ., by = c('gbif_name' = 'worms_name'))

# rejoin
dt_spcheck <- dt_spcheck %>% filter(!genus_sp %in% gbif_match$genus_sp) %>%
  bind_rows(., gbif_match) %>% arrange(genus_sp)

dt_spcheck <- dt_spcheck %>% select(!starts_with('gbif'))
dt_spcheck %>% filter(is.na(worms_valid)) %>% View
rm(dt_spcheck2)

# manual validation for these
dt_spcheck <- dt_spcheck %>% 
  mutate(worms_valid = case_when(.default = worms_valid,
                                 genus_sp == "Calostomus coralinus" ~ 'Calotomus carolinus',
                                 genus_sp == "Acanthochromis japonicus" ~ 'Acanthochromis polyacanthus',
                                 genus_sp == "Cantherhines maculatus" ~ 'Cantherhines maculata',
                                 genus_sp == "Ecsenius flavus" ~ 'Plagiotremus flavus',
                                 genus_sp == "Enneanectes tutuilae" ~ 'Enneapterygius tutuilae',
                                 genus_sp == "Epinephelides guttatus" ~ 'Epinephelus guttatus',
                                 genus_sp == "Holocanthus townsendi" ~ 'Holacanthus townsendi', # hybrid
                                 genus_sp == "Hypoplectrus tann" ~ 'Hypoplectrus randallorum',
                                 genus_sp == "Labroides manabei" ~ 'Labropsis manabei',
                                 genus_sp == "Ucla macrospilus" ~ 'Acanthemblemaria macrospilus',
                                 genus_sp == "Chromis tricincta" ~ 'Chrysiptera tricincta',
                                 genus_sp == "Plectorhinchus pictus" ~ 'Plectorhinchus picus',
                                 genus_sp == 'Stegostoma varium' ~ 'Stegostoma tigrinum'))


## IMPLEMENT FIXES =================================================================================

dt_cleanfish <- dt_spcheck %>% select(genus_sp, worms_valid) %>% 
  left_join(dt_assemblage, ., by = "genus_sp")

dt_cleanfish %>% filter(taxon_resolution != 'species') %>% 
  distinct(Family, Genus, Species, genus_sp, taxon_resolution) %>% arrange(genus_sp) %>% View
dt_cleanfish %>% 
  distinct(Family, Genus, Species, genus_sp, taxon_resolution) %>% 
  filter(str_detect(genus_sp, 'idae')) %>% arrange(genus_sp) %>% View

# eliminate non-fish
notfish <- dt_spcheck %>% 
  filter(is.na(class) == F) %>% 
  filter(!class %in% c('Teleostei', 'Elasmobranchii')) %>%
  pull(genus_sp)

dt_cleanfish <- dt_cleanfish %>% filter(!genus_sp %in% notfish)
notfish <- word(notfish, 1) # also use the genera to filter it out
dt_cleanfish <- dt_cleanfish %>% filter(!Genus %in% notfish)

# double check
rm(sp, nomatch, worms_valid, replace, dt_spcheck)
# create new clean columns
dt_assemblage <- dt_assemblage %>% rename(genus_sp_valid = genus_sp) %>% 
  # extract genus from validated genus_sp
  mutate(Genus_valid = case_when(is.na(Genus) == F ~ word(genus_sp_valid, 1), 
                                 .default = NA),
         # get epithet from validated genus_sp or original if it's a specified uncertain
         Species_valid = case_when(taxon_resolution == 'species' ~ word(genus_sp_valid, 2),
                                   str_detect(taxon_resolution, 'specified') ~ Species,
                                   taxon_resolution == 'genus' ~ Species,
                                   .default = NA))

## WoRMS and FISHBASE HARMONISE=====================================================================

fb_valid <- dt_assemblage %>% filter(taxon_resolution == 'species') %>% 
  distinct(genus_sp_valid) %>% arrange(genus_sp_valid)
fb_sp <- species(species_list = fb_valid$genus_sp_valid) # get fishbase valid species
fb_valid <- left_join(fb_valid, fb_sp[1:2], by = c('genus_sp_valid' = 'Species'))

## VALIDATE MATCH GENUS/FAMILY =====================================================================

# use taxize to interact with worms API better
gen <- dt_assemblage$Genus_valid |> unique() |> sort() |> na.omit()
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


