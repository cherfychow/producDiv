

require(sf)
require(rnaturalearth)
require(dplyr)
require(readr)
require(ggplot2)
require(rfishbase)
require(worrms)
require(taxize)
require(stringr)
require(purrr)

# load in data
dt_assemblage <- read_csv('src/macro_reef_fish_obs.csv')
dt_meta <- read_csv('./src/macro_reef_fish_meta.csv') # dataset level metadata
gc() # keeping memory usage low

summary(dt_assemblage)

# TAXONOMY CLEAN ===================================================================================

dt_assemblage %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>%  View

# clean epithets
# remove genera, but allow the cf. records
dt_assemblage$Species[-str_which(dt_assemblage$Species, '^cf')] <- str_remove(dt_assemblage$Species[-str_which(dt_assemblage$Species, '^cf')], '^[A-Z]{1}[a-z]+\\s')
dt_assemblage$Species[-str_which(dt_assemblage$Species, '^cf')] <- word(dt_assemblage$Species[-str_which(dt_assemblage$Species, '^cf')], 1) # just epithets, no subspecies

# uncertains
dt_assemblage %>% distinct(Genus, Species) %>% filter(str_detect(Species, '^sp.{1,2}$')) %>% View
dt_assemblage$Species[dt_assemblage$Species %in% c('sp.', 'spp', 'species', 'spies')] <- 'sp'
dt_assemblage$Species[dt_assemblage$Species == 'spx'] <- 'sp1'
dt_assemblage$Species[dt_assemblage$Species == 'spz'] <- 'sp2'
dt_assemblage$Species[dt_assemblage$Species == 'spd'] <- 'sp1'
dt_assemblage$Family[dt_assemblage$Genus == 'Blenny'] <- 'Blenniidae'
dt_assemblage$Species[dt_assemblage$Genus == 'Blenny'] <- NA
dt_assemblage$Genus[dt_assemblage$Genus == 'Blenny'] <- NA
dt_assemblage$Species <- str_replace(dt_assemblage$Species, 'i{3}$', 'ii')

# parrotfish specific
dt_assemblage$Species[dt_assemblage$Genus == 'Scarinissp8'] <- 'sp8'
dt_assemblage$Genus[dt_assemblage$Genus == 'Scarinissp8'] <- 'Scarini'
dt_assemblage$Species[str_which(dt_assemblage$Species, '^\\d$')] <- paste0('sp', dt_assemblage$Species[str_which(dt_assemblage$Species, '^\\d$')])
dt_assemblage$Genus[dt_assemblage$Family == 'Scaridae' & is.na(dt_assemblage$Genus)] <- 'Scarini'
dt_assemblage$Family[dt_assemblage$Family == 'Scaridae'] <- 'Labridae'
dt_assemblage$Genus <- str_to_title(dt_assemblage$Genus)
dt_assemblage$Species <- str_to_lower(dt_assemblage$Species)

# omits
dt_assemblage <- dt_assemblage %>% filter(!str_detect(Species, 'crcp|franc'))
dt_assemblage <- filter(dt_assemblage, Genus != 'Teleostei') # omit


# make a column to help with sorting uncertain records
dt_assemblage <- dt_assemblage %>% 
  # family = Gobiidae
  # genus = Scarus sp
  # specified = working name distinguishable from other uncertains, e.g. sp1 and sp2
  mutate(taxon_resolution = case_when(str_detect(Species, 'sp$') & is.na(Genus) == F & str_detect(Genus, 'Scari(ni|nae)') == F ~ 'genus', 
                                      str_detect(Species, 'sp\\d+$') & is.na(Genus) == F ~ 'genus specified', 
                                      is.na(Family) == F & is.na(Genus) ~ 'family',
                                      is.na(Family) == F & is.na(Genus) & str_detect(Species, '\\d') ~ 'family specified',
                                      str_detect(Species, 'sp\\d*$') == F & is.na(Genus) == F ~ 'species',
                                      str_detect(Genus, 'Scari(ni|nae)') ~ 'family'))

# family level consistency, specified from unspecified
dt_assemblage <- mutate(dt_assemblage, 
                        Species = case_when(taxon_resolution == 'family' & Species == 'sp' ~ NA, 
                                            .default = Species))
dt_assemblage$Species[str_detect(dt_assemblage$Genus, 'Scarini') & str_detect(dt_assemblage$Species, '\\d') == F] <- NA

dt_assemblage %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>%  View
dt_assemblage %>% filter(str_detect(Species, '\\s')) %>% View

## Automate misspelling fix========================================================================

# generate combined taxonomy
dt_assemblage <- dt_assemblage %>% mutate(genus_sp = case_when(
  taxon_resolution == 'family' ~ Family,
  taxon_resolution == 'family specified' ~ paste(Family, Species, sep = ' '),
  taxon_resolution %in% c('species', 'genus specified', 'genus') ~ paste(Genus, Species, sep = ' ')
))

# generate species list for checking
sp <- dt_assemblage %>% filter(taxon_resolution == 'species') %>% 
  filter(str_detect(Species, '\\s|[:punct:]|\\d') == F) %>% 
  distinct(genus_sp) %>% pull %>% sort

# fetch taxonomic IDs from fuzzy matching
dt_spcheck <- tibble(genus_sp = sp)
dt_spcheck$worms_id <- get_wormsid(dt_spcheck$genus_sp, searchtype = 'scientific', fuzzy = T,
                                   marine_only = T, messages = F)

# first round of catching unaccepted names
unaccepted <- dt_spcheck %>% filter(is.na(worms_id)) %>% select(genus_sp)

# get names. filter out synonyms
dt_spcheck <- dt_spcheck %>% filter(is.na(worms_id) == F) %>%
  mutate(worms_id = as.numeric(worms_id))
worms <- wm_record(dt_spcheck$worms_id) %>% bind_rows
gc()

# cephalopods and seals???
notfish <- worms %>% filter(!class %in% c('Teleostei', 'Elasmobranchii'))
dt_assemblage <- dt_assemblage %>% filter(!genus_sp %in% notfish$scientificname) # remove
# make dt_spcheck the correction object
dt_spcheck <- worms %>% filter(status != 'accepted') %>% inner_join(dt_spcheck, ., by = c('worms_id' = 'AphiaID'))

# alternative way to match these wildly off name misspellings
gbif_id <- get_gbifid(worms_unaccepted$genus_sp, messages = F, method = 'backbone')
gbif_unaccepted <- taxize::id2name(id = gbif_unaccepted$gbif_id)

# isolate names that still can't be matched
nomatch <- unaccepted$genus_sp[which(map(gbif_unaccepted, length) == 1)]
gbif_unaccepted <- gbif_unaccepted[which(map(gbif_unaccepted, length) > 1)] %>% bind_rows
unaccepted <- unaccepted %>% filter(!genus_sp %in% nomatch) %>% bind_cols(., gbif_unaccepted)
# add the poor matches to nomatch
nomatch <- unaccepted %>% filter(rank != 'species') %>% pull(genus_sp) %>% append(nomatch, .)
unaccepted <- unaccepted %>% filter(status != 'accepted', rank == 'species')
rm(worms, gbif_unaccepted, notfish)
nomatch

# get valid names for unaccepted
worms_valid <- wm_records_names(unaccepted$name, marine_only = T) %>% bind_rows
rm(unaccepted, worms_valid)

nomatch <- nomatch[-which(nomatch == 'Apogon red')]
nomatch <- nomatch[-which(nomatch == 'Menaphorus punticulatus')] # unresolved species
nomatch <- nomatch[-which(nomatch == 'Hypoplectrus hybrid')]
# manual validation for these

names(nomatch) <- nomatch
nomatch <- c(
  "Calostomus coralinus" = 'Calotomus carolinus',
  "Acanthochromis japonicus" = 'Acanthochromis polyacanthus',
  "Cantherhines maculatus"= 'Cantherhines maculata',
  "Ecsenius flavus" = 'Plagiotremus flavus',
  "Enneanectes tutuilae" = 'Enneapterygius tutuilae',
  "Epinephelides guttatus" = 'Epinephelus guttatus',
  "Holocanthus townsendi" = 'Holacanthus townsendi', # hybrid
  "Hypoplectrus tann" = 'Hypoplectrus randallorum',
  "Labroides manabei" = 'Labropsis manabei',
  "Pomacentrus margaritifer" = 'Pycnochromis margaritifer',
  "Ucla macrospilus" = 'Acanthemblemaria macrospilus',
  "Chromis tricincta" = 'Chrysiptera tricincta'
)

## IMPLEMENT FIXES =================================================================================

replace <- dt_spcheck$valid_name
names(replace) <- dt_spcheck$genus_sp
replace
nomatch

dt_assemblage$genus_sp <- str_replace_all(dt_assemblage$genus_sp, replace)
dt_assemblage$genus_sp <- str_replace_all(dt_assemblage$genus_sp, nomatch)

replace <- worms_valid$valid_name
names(replace) <- worms_valid$scientificname
dt_assemblage$genus_sp <- str_replace_all(dt_assemblage$genus_sp, replace)

dt_assemblage %>% distinct(Family, Genus, Species, genus_sp) %>% 
  arrange(Family, genus_sp) %>% View

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


