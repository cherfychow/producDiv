
# ProducDiv
# Calculate productivity

require(dplyr)
require(rfishprod)
require(rfishbase)

# Compile traits ==========================================================
# load("/Users/cherchow/Dropbox/ch4_producDiv/datarepo/workspace.RData")

sp_gap <- dt_assemblage %>% distinct(Family, Genus, Species) %>%
  mutate(sp_working = case_when(is.na(Genus) == F ~ paste(Genus, Species, sep = ' ') %>%
	str_remove(., '\\d+'), .default = NA)) %>% 
  # collapse known uncertains because they don't matter for productivity
  filter(!sp_working %in% db$Species) %>% select(!c(Genus, Species)) %>% 
	# filter out 
  arrange(Family, sp_working) %>% distinct
View(sp_gap)

# fetch traits from FishBase
# L-W parameters, growth
# MaxSizeTL
# Diet
# Position
# SSTMean (from Growth)
# Growth method
# Species size
# L0
View(rfishprod:::repdata)
View(db)

traits <- as.list(rep(NA, 7)) # empty list to store trait tables


# Compile SST ==========================================================

# Getting levels ready #
repdata <- tidytrait(repdata, db)

# Estimate Kmax ==========================================================
# estimate VGBF Kmax
# input trait data
# dataset = reference data to train

# Formula from Morais and Bellwood (2018)
fmod <- formula(~ sstmean + MaxSizeTL + Diet + Position + Method)
dt_kmax <- predKmax(repdata, dataset = db,
                   fmod = fmod, niter = iter,
                   return = 'pred')

