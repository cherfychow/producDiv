
# ProducDiv
# Calculate productivity from fish traits + SST with observed lengths

require(dplyr)
require(stringr)
require(rfishprod)
require(beepr)

# PREP =========================================================================

# give predKmax data on all fishes observed
# then calculate productivity per rareID

# need to keep some sampling info since SST are specific to location + time
# dt_assemblage, because it's full data
# traits
prodtrait <- dt_assemblage %>% 
  distinct(DataID, gridcell, SamplingEventID, Year, Month, genus_sp, Length) %>%
  # all family level = just family
  mutate(genus_sp = case_when(str_detect(genus_sp, '(?<=idae)\\ssp\\d*$') ~ str_remove(genus_sp, '\\ssp\\d*$'),
                              .default = str_remove(genus_sp, '\\d+'))) %>%
  # known uncertains irrelevant for productivity
  left_join(., dt_traits, by = c('genus_sp' = 'Species')) # traits

# SST
prodtrait <- left_join(prodtrait, dt_sst[-c(1,5,6)], 
                       by = c('gridcell', 'Year', 'Month')) %>% # SST
  rename(sstmean = sst, Size = Length, Species = genus_sp)

# our own version of tidytrait
sort(unique(prodtrait$Diet)) == sort(unique(dt_prodDB$Diet))
sort(unique(prodtrait$Position)) == sort(unique(dt_prodDB$Position))
# have to fit the db's factor ordering
prodtrait <- prodtrait %>% 
  mutate(Method = 'Otolith',
         Diet = factor(Diet, levels = sort(unique(prodtrait$Diet))),
         Position = factor(Position, levels = levels(dt_prodDB$Position)))
summary(prodtrait)
prodtrait <- prodtrait %>% filter(!is.na(b)) # only 2 very uncertain records
# make db match
dt_prodDB$Diet <- factor(dt_prodDB$Diet, levels = sort(unique(prodtrait$Diet)))


# check if any fish are larger than their max size
prodtrait %>% filter(Size > MaxSizeTL) %>% nrow

# limit observed size to be below max size
prodtrait <- mutate(prodtrait, 
                    Size = case_when(Size > MaxSizeTL ~ MaxSizeTL - 0.1, 
                                     .default = Size))

# Formula from Morais and Bellwood (2018) #
growthf <- formula(~ sstmean + MaxSizeTL + Diet + Position)

# Kmax, growth coef ============================================================

set.seed(24)
vbgf_par <- prodtrait %>% distinct %>% 
  predKmax(., dataset = dt_prodDB,
           fmod = growthf, niter = 100, return = 'pred')

vbgf_par <- vbgf_par$pred # just keep predictions
rm(prodtrait) # kept in vgbf_par
summary(vbgf_par)
# negative Kmax values
vbgf_par %>% filter(Kmax <= 0) %>% nrow() # work on this... see if 1000 runs solves it

# Mortality and gross somatic gain =============================================

vbgf_par$M <- with(vbgf_par, predM(Lmeas = Size, Lmax = MaxSizeTL, 
                                      Kmax = Kmax, method = 'Lorenzen'))
vbgf_par <- vbgf_par %>% filter(Kmax > 0) %>% # leave out the nonsensicals for now
  mutate(Size_t1 = applyVBGF(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax), # rate of change day^-1
         size_delta = Size_t1 - Size) # length change

vbgf_par <- vbgf_par %>% 
  mutate(biomass_gain = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax),
         biomass_loss = somaLoss(M, Lmeas = Size, a, b),
         biomass_delta = biomass_gain - biomass_loss)

summary(vbgf_par)

# Apply to assemblage data =====================================================

colnames(dt_rareass)
dt_rareass <- dt_rareass %>% 
  # known uncertain blind species names for productivity calculations
  mutate(prod_genus_sp = case_when(str_detect(genus_sp, '(?<=idae)\\ssp\\d*$') ~ str_remove(genus_sp, '\\ssp\\d*$'),
                                   .default = str_remove(genus_sp, '\\d+')))

dt_prod <- vbgf_par %>% select(DataID:Species, 
                    Size, size_delta, biomass_delta) %>% 
  left_join(dt_rareass, ., by = c('DataID', 'gridcell', 'SamplingEventID',
                                  'Year', 'Month', 'prod_genus_sp' = 'Species',
                                  'Length' = 'Size'))
# aggregate to species level
dt_prod <- dt_prod %>% 
  group_by_at(vars(-c(Abundance, size_delta, biomass_delta))) %>%
  # multiply by count
  summarise(size_delta = Abundance * size_delta,
            biomass_delta = Abundance * biomass_delta) %>% 
  ungroup %>% group_by_at(vars(-c(Length, size_delta, biomass_delta))) %>% 
  # aggregate total length and biomass change for the species
  summarise(size_delta = sum(size_delta, na.rm = T),
            biomass_delta = sum(biomass_delta, na.rm = T)) %>% ungroup

# aggregate to site level

dt_prod <- dt_prod %>% ungroup %>% 
  group_by(pick(SamplingEventID:DataID, gridcell, rareID)) %>%
  summarise(size_delta = sum(size_delta, na.rm = T), 
            biomass_delta = sum(biomass_delta, na.rm = T)) %>% 
  # convert to kg/yr
  # GET STUDY GRAIN SIZE
  mutate(biomass_delta_yr = biomass_delta * 365 / 1000) %>% ungroup

## Visualise ===================================================================

require(ggplot2)

ggplot() +
  