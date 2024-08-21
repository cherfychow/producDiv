
# ProducDiv
# Calculate productivity from fish traits + SST with observed lengths

require(dplyr)
require(stringr)
require(rfishprod)
require(beepr)

# PREP =============================================================================================

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
prodtrait %>% filter(Size >= MaxSizeTL) %>% nrow

# limit observed size to be below max size
prodtrait <- mutate(prodtrait, 
                    Size = case_when(Size > MaxSizeTL ~ MaxSizeTL - 0.1, 
                                     .default = Size))

# Formula from Morais and Bellwood (2018) #
growthf <- formula(~ sstmean + MaxSizeTL + Diet + Position)

# Kmax, growth coef ================================================================================

# Kmax estimate per species-gridcell-year-dataID / species-samplingeventid-dataID
# recommended boostrapts on 100s to 1000s scale
set.seed(24)
vbgf_par <- prodtrait %>% distinct %>% 
  predKmax(., dataset = dt_prodDB,
           fmod = growthf, niter = 500, return = 'pred')

vbgf_par <- vbgf_par$pred # just keep predictions
rm(prodtrait, growthf)
beep()

summary(vbgf_par)
# negative Kmax values
vbgf_par %>% filter(Kmax <= 0) %>% nrow() # work on this later... see if 1000 runs solves it

# Mortality and gross somatic gain =================================================================

# mortality
vbgf_par <- vbgf_par %>% mutate(
  M = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax, method = 'Lorenzen'),
  M_l = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_lowq, method = 'Lorenzen'),
  M_u = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_uppq, method = 'Lorenzen'))

vbgf_par <- vbgf_par %>% filter(Kmax > 0) %>% # leave out the nonsensicals for now
  mutate(
    Size_t1 = applyVBGF(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax), # rate of change day^-1
    size_delta = Size_t1 - Size
  ) # length change

# lowers and uppers refer to the Kmax value used!
# not necessarily the lower or upper of that metric
vbgf_par <- vbgf_par %>% 
  mutate(biomass_gain = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax),
         biomass_loss = somaLoss(M, Lmeas = Size, a, b),
         delta_g_day = biomass_gain - biomass_loss,
         biomass_gain_l = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         biomass_loss_l = somaLoss(M_l, Lmeas = Size, a, b),
         delta_g_day_l = biomass_gain_l - biomass_loss_l,
         biomass_gain_u = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_uppq),
         biomass_loss_u = somaLoss(M_u, Lmeas = Size, a, b),
         delta_g_day_u = biomass_gain_u - biomass_loss)

summary(vbgf_par)
# don't need to keep traits
vbgf_par <- vbgf_par %>% mutate(rareID = paste(DataID, gridcell, sep = '_')) %>% 
  select(rareID, SamplingEventID, Species, Size, Kmax:delta_g_day_u)

# Apply to assemblage data =========================================================================

colnames(dt_rareass)
dt_rareass <- dt_rareass %>% 
  # known uncertain blind species names for productivity calculations
  mutate(prod_genus_sp = case_when(str_detect(genus_sp, '(?<=idae)\\ssp\\d*$') ~ str_remove(genus_sp, '\\ssp\\d*$'),
                                   .default = str_remove(genus_sp, '\\d+')))

# match up
dt_prod <- vbgf_par %>% select(rareID:Size, starts_with('delta')) %>% 
  left_join(dt_rareass, ., by = c('rareID', 'SamplingEventID', 
                                  'prod_genus_sp' = 'Species',
                                  'Length' = 'Size'))
# aggregate to species level
dt_prod <- dt_prod %>% 
  group_by(pick(!c(Abundance, starts_with('delta')))) %>%
  # multiply by count
  summarise(delta_g_day = Abundance * delta_g_day,
            delta_g_day_l = Abundance * delta_g_day_l,
            delta_g_day_u = Abundance * delta_g_day_u) %>% 
  ungroup %>% group_by(pick(!c(Length, delta_g_day, delta_g_day_l, delta_g_day_u))) %>% 
  # aggregate total length and biomass change for the species
  summarise(delta_g_day = sum(delta_g_day, na.rm = T),
            delta_g_day_l = sum(delta_g_day_l, na.rm = T),
            delta_g_day_u = sum(delta_g_day_u, na.rm = T)) %>% ungroup

# aggregate to site level

dt_prod <- dt_prod %>% ungroup %>% 
  group_by(rareID, Year, SamplingEventID) %>%
  summarise(delta_g_day = sum(delta_g_day, na.rm = T),
            delta_g_day_l = sum(delta_g_day_l, na.rm = T),
            delta_g_day_u = sum(delta_g_day_u, na.rm = T)) %>% 
  # convert to kg/yr
  mutate(delta_g_day = delta_g_day * 365 / 1000,
         delta_g_day_l = delta_g_day_l * 365 / 1000,
         delta_g_day_u = delta_g_day_u * 365 / 1000) %>% ungroup

colnames(dt_prod) <- str_replace(colnames(dt_prod), 'delta\\_g\\_day', 'prod_kg_yr')

# kg/ha/yr
dt_prod <- dt_rareass %>% distinct(rareID, DataID) %>% 
  left_join(dt_meta[c(1,7)], ., by = "DataID") %>% 
  mutate(Area = str_extract(`Grain size`, '\\d{2,}')) %>% 
  select(DataID, rareID, Area) %>% 
  left_join(dt_prod, ., by = "rareID")
dt_prod$Area <- as.numeric(dt_prod$Area)
# unit conversion
dt_prod <- dt_prod %>% mutate(across(starts_with('prod'), ~ .x * 10000 / Area))


  

