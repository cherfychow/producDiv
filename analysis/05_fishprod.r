
# ProducDiv
# Calculate productivity from fish traits + SST with observed lengths

require(dplyr)
require(stringr)
require(rfishprod)
require(beepr)

# PREP ====================================================================

# give predKmax data on all fishes observed
# then calculate productivity per rareID

# need to keep some sampling info since SST are specific to location + time
# dt_assemblage, because it's full data
# traits
prodtrait <- dt_assemblage %>% distinct(DataID, gridcell, SamplingEventID, Year, Month, genus_sp, Length) %>%
  # all family level = just family
  mutate(genus_sp = case_when(str_detect(genus_sp, '(?<=idae)\\ssp\\d*$') ~ str_remove(genus_sp, '\\ssp\\d*$'),
                              .default = str_remove(genus_sp, '\\d+'))) %>%
  # known uncertains irrelevant for productivity
  left_join(., dt_traits, by = c('genus_sp' = 'Species')) # traits

# SST
prodtrait <- left_join(prodtrait, dt_sst[-c(1,5,6)], by = c('gridcell', 'Year', 'Month')) %>% # SST
  rename(sstmean = sst, Size = Length, Species = genus_sp)

# our own version of tidytrait
sort(unique(prodtrait$Diet)) == sort(unique(dt_prodDB$Diet))
sort(unique(prodtrait$Position)) == sort(unique(dt_prodDB$Position))
# have to fit the db's factor ordering
prodtrait <- prodtrait %>% mutate(Method = 'Otolith',
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
                    Size = case_when(Size > MaxSizeTL ~ MaxSizeTL - 0.1, .default = Size))

# Formula from Morais and Bellwood (2018) #
growthf <- formula(~ sstmean + MaxSizeTL + Diet + Position)

# FIT MODELS ====================================================================

set.seed(24)
# Predicting Kmax, the standardised VBGF parameter (Recommendation: use 100s to 1000s iterations) #
kmax <- prodtrait %>% 
  predKmax(., dataset = dt_prodDB,
           fmod = growthf, niter = 100, return = 'pred')
beep()