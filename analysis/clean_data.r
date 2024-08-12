

require(sf)
require(rnaturalearth)
require(dplyr)
require(readr)
require(ggplot2)
require(taxize) # helps handle large inputs to API better
require(worrms)
require(stringr)

# load in data
dt_assemblage <- read_csv('src/macro_reef_fish_obs.csv')
dt_meta <- read_csv('./src/macro_reef_fish_meta.csv') # dataset level metadata
gc() # keeping memory usage low

# CHECK + VALIDATE DATA ============================================================================

summary(dt_assemblage)

## TAXONOMY ========================================================================================

# make a column to help with sorting uncertain records
dt_assemblage <- dt_assemblage %>% 
  # family = Gobiidae
  # genus = Scarus sp
  # specified = working name distinguishable from other uncertains, e.g. sp1 and sp2
  mutate(taxon_resolution = case_when(str_detect(Species, 'sp$') & is.na(Genus) == F ~ 'genus', 
                                      str_detect(Species, 'sp\\d+$') & is.na(Genus) == F ~ 'genus specified', 
                                      is.na(Family) == F & is.na(Genus) & is.na(Species) ~ 'family',
                                      is.na(Family) == F & is.na(Genus) & is.na(Species) == F ~ 'family specified',
                                      str_detect(Species, 'sp\\d*$') == F ~ 'species',
                                      str_detect(Genus, 'Scari(ni|nae)') ~ 'family'))

# family level consistency, specified from unspecified
dt_assemblage <- mutate(dt_assemblage, 
                        Species = case_when(taxon_resolution == 'family' & Species == 'sp' ~ NA, 
                                            .default = Species))

dt_assemblage %>% distinct(Family, Genus, Species) %>% arrange(Family, Genus, Species) %>%  View
dt_assemblage %>% filter(str_detect(Species, '\\s')) %>% View
# clean epithets
# remove genera, but allow the cf. records
dt_assemblage$Species[-str_which(dt_assemblage$Species, '^cf')] <- str_remove(dt_assemblage$Species[-str_which(dt_assemblage$Species, '^cf')], '^[A-Z]{1}[a-z]+\\s')
dt_assemblage$Species[-str_which(dt_assemblage$Species, '^cf')] <- word(dt_assemblage$Species[-str_which(dt_assemblage$Species, '^cf')], 1) # just epithets, no subspecies

# uncertains
dt_assemblage %>% distinct(Genus, Species) %>% filter(str_detect(Species, '^sp.{1,2}$')) %>% View
dt_assemblage$Species[dt_assemblage$Species %in% c('sp.', 'spp')] <- 'sp'

dt_assemblage <- dt_assemblage %>% 
  mutate(taxon_resolution = case_when(str_detect(Species, 'sp$') & is.na(Genus) == F ~ 'genus', 
                                      str_detect(Species, 'sp\\d+$') & is.na(Genus) == F ~ 'genus specified', 
                                      is.na(Family) == F & is.na(Genus) & is.na(Species) ~ 'family',
                                      is.na(Family) == F & is.na(Genus) & is.na(Species) == F ~ 'family specified',
                                      str_detect(Species, 'sp\\d*$') == F ~ 'species',
                                      str_detect(Genus, 'Scari(ni|nae)') ~ 'family'))

### Clean misspellings =============================================================================

# manual uncertain misspelling replace
dt_assemblage %>% filter(str_detect(Genus, 'Scarini')) %>%
  distinct(Family, Genus, Species) %>% View
dt_assemblage$Species[dt_assemblage$Genus == 'Scarinissp8'] <- 'sp8'
dt_assemblage$Genus[dt_assemblage$Genus == 'Scarinissp8'] <- 'Scarini'
dt_assemblage$Species[str_which(dt_assemblage$Species, '^\\d$')] <- paste0('sp', dt_assemblage$Species[str_which(dt_assemblage$Species, '^\\d$')])

# parrotfish specific
dt_assemblage$Genus[dt_assemblage$Family == 'Scaridae' & is.na(dt_assemblage$Genus)] <- 'Scarini'
dt_assemblage$Family[dt_assemblage$Family == 'Scaridae'] <- 'Labridae'
dt_assemblage$Genus <- str_to_title(dt_assemblage$Genus)
dt_assemblage$Species <- str_to_lower(dt_assemblage$Species)
dt_assemblage$Species <- str_replace(dt_assemblage$Species, 'i{3}$', 'ii')

# generate combined taxonomy
dt_assemblage <- dt_assemblage %>% mutate(genus_sp = case_when(
  taxon_resolution == 'family' ~ Family,
  taxon_resolution == 'family specified' ~ paste(Family, Species, sep = ' '),
  taxon_resolution %in% c('species', 'genus specified', 'genus') ~ paste(Genus, Species, sep = ' ')
))

# omits
dt_assemblage$Family[dt_assemblage$genus_sp == 'Blenny species'] <- 'Blenniidae'
dt_assemblage$Genus[dt_assemblage$genus_sp == 'Blenny species'] <- NA
dt_assemblage$Species[dt_assemblage$genus_sp == 'Blenny species'] <- NA
dt_assemblage$genus_sp[dt_assemblage$genus_sp == 'Blenny species'] <- NA
dt_assemblage <- dt_assemblage %>% filter(!str_detect(genus_sp, '^NA'))
dt_assemblage <- filter(dt_assemblage, genus_sp != 'Teleostei unidentified') # omit
dt_assemblage$Species[dt_assemblage$Species == 'species'] <- 'sp'
dt_assemblage$Species[dt_assemblage$Species == 'spx'] <- 'sp1'
dt_assemblage$Species[dt_assemblage$Species == 'spz'] <- 'sp2'

# redo
dt_assemblage <- dt_assemblage %>% mutate(genus_sp = case_when(
  taxon_resolution == 'family' ~ Family,
  taxon_resolution == 'family specified' ~ paste(Family, Species, sep = ' '),
  taxon_resolution %in% c('species', 'genus specified', 'genus') ~ paste(Genus, Species, sep = ' ')
))

### Database validation ============================================================================
# generate species list for checking
sp <- dt_assemblage %>% filter(taxon_resolution == 'species') %>% 
  filter(str_detect(Species, '\\s|[:punct:]') == F) %>% 
  distinct(genus_sp) %>% pull %>% sort

# fetch taxonomic IDs from fuzzy matching
dt_spcheck <- tibble(genus_sp = sp)
# API takes a while
# will prompt selections for conflicts
worms_class <- classification(sp, messages = F, db = 'worms')

# filter out ones that failed to match
worms_unaccepted <- names(worms_class[which(map_dbl(worms_class, length) == 1)])

worms_class <- worms_class[which(map_dbl(worms_class, length) > 1)] %>% map(., function(x) {
  x %>% select(!id) %>% tidyr::pivot_wider(names_from = 'rank', values_from = 'name')
})

for (i in 1:length(worms_class)) {
  worms_class[[i]]$input <- names(worms_class)[i]
}

worms_class <- bind_rows(worms_class)
dt_spcheck <- rfishbase::synonyms(sp)
dt_spcheck %>% filter(Status != 'accepted name') %>% View

# direct fuzzy match through worms for unaccepteds
worms_unaccepted <- worms_accepted %>% filter(status != 'accepted') %>% pull(name) %>% 
  c(., dt_spcheck$genus_sp[is.na(dt_spcheck$worms_id)]) %>% 
  classification(., db = 'worms')

# cephalopods and seals???
notfish <- worms_accepted %>% 
  filter(!class %in% c('Teleostei', 'Elasmobranchii')) %>% 
  pull(scientificname)

dt_assemblage <- dt_assemblage %>% filter(!genus_sp %in% notfish)

# implement corrections
replace <- spcheck$valid
names(replace) <- spcheck$input

dt_assemblage$genus_sp <- str_replace_all(dt_assemblage$genus_sp, replace)
dt_assemblage %>% filter(str_detect(genus_sp, '^NA')) %>% distinct %>% View

# double check
rm(spcheck, replace, dt_spcheck)

# redo classification
dt_assemblage <- dt_assemblage %>% 
  # family = Gobiidae
  # genus = Scarus sp
  # specified = working name distinguishable from other uncertains, e.g. sp1 and sp2
  mutate(taxon_resolution = case_when(str_detect(Species, 'sp$') & is.na(Genus) == F ~ 'genus', 
                                      str_detect(Species, 'sp\\d+$') & is.na(Genus) == F ~ 'genus specified', 
                                      is.na(Family) == F & is.na(Genus) & is.na(Species) ~ 'family',
                                      is.na(Family) == F & is.na(Genus) & is.na(Species) == F ~ 'family specified',
                                      str_detect(Species, 'sp\\d*$') == F ~ 'species',
                                      str_detect(Genus, 'Scari(ni|nae)') ~ 'family'))

# take the validated names and resplit to fix the Genus and Species column
dt_assemblage <- mutate(dt_assemblage, Genus2 = case_when(str_detect(word(genus_sp, 1), 'idae$') ~ NA,
                                                          .default = word(genus_sp, 1)), 
                        # family level IDs should be NA in Species, unless it's Scarini or a known uncertain
                        Species2 = case_when(is.na(word(genus_sp, 2)) == F ~ word(genus_sp, 2), .default = NA))

dt_assemblage %>% distinct(Family, Genus, Species, Genus2, Species2) %>%
  filter(Species2 != Species | Genus2 != Genus) %>%
  arrange(Family, Genus2, Species2) %>% View
dt_assemblage %>% distinct(Family, Genus, Species, Genus2, Species2) %>%
  filter(str_detect(Species, '\\d$|sp$')) %>%
  arrange(Family, Genus2, Species2) %>% View # check if the uncertainty information is preserved;

# looks ok. replace
dt_assemblage <- dt_assemblage %>% select(!c(Genus, Species)) %>% 
  rename(Genus = Genus2, Species = Species2)

rm(sp, dt_spcheck, replace)

## VALIDATE GENUS/FAMILY -------------------------------------------------------

# use taxize to interact with worms API better
gen <- dt_assemblage$Genus |> unique() |> sort() |> na.omit()
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


